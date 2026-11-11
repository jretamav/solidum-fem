# solidum_fem/solidum/math/solvers/arclength.py
"""``ArcLengthSolver`` — método de longitud de arco cilíndrico (Crisfield).

Estructura del corrector (un ensamblaje por iteración, estado coherente)
------------------------------------------------------------------------
Cada iteración del corrector ensambla una sola vez en el iterado corriente
``U_k`` y evalúa la convergencia **antes** de resolver, con el par
``(‖R(U_k)‖, ‖δU_{k−1}‖)``. Al converger, el estado trial de los elementos
es exactamente el del ensamblaje en ``U_k`` y se comitea junto con
``U_current = U_k`` y ``λ_curr = λ_k``. La versión anterior evaluaba el
residuo en ``U_k`` pero almacenaba ``U_{k+1}``: el estado committed iba
un iterado por detrás del desplazamiento guardado (auditoría 2026-09-22).

El bucle es el :class:`~solidum.math.solvers.corrector.NewtonCorrector`
compartido (ADR 0015). Lo propio del arc-length —predictor tangente,
las dos resoluciones por iteración (``du_R`` y ``du_t``) y la restricción
que fija ``ddλ``— vive en :class:`_ArcProblem`; una variante que sólo
cambia la magnitud controlada redefine su método ``constraint`` (p. ej.
:class:`IndirectDisplacementSolver`). El paquete lo exporta con nombre
público, ``solidum.math.solvers.ArcProblem``, para las variantes de fuera.

Primer paso (2026-09-23)
------------------------
La longitud de arco es una longitud en unidades de desplazamiento, y el
usuario no conoce de antemano los desplazamientos de su modelo. Por eso el
primer paso se declara, por omisión, como **fracción de la carga de
referencia** (``initial_dlambda`` = Δλ₁, Crisfield 1991, cap. 9) y el solver
deriva la longitud con el predictor elástico: Δl₁ = Δλ₁·‖K⁻¹·F_ref‖. El
primer paso lleva así la misma fracción de carga en cualquier sistema de
unidades. ``initial_dl`` sigue disponible como longitud explícita,
excluyente con ``initial_dlambda``. Antes el valor por omisión era
``initial_dl = 0.1``: en el ejemplo 4 del manual de ejemplos eso pedía 86
veces la carga de colapso en el primer paso.

Llegada a ``max_lambda`` (2026-09-23)
-------------------------------------
Todos los pasos se dan con la restricción de arco. Si uno converge con
λ > ``max_lambda``, su estado **no se consolida**: el paso se repite desde el
último estado convergido con λ = ``max_lambda`` fijo (Newton puro), partiendo
de la interpolación lineal sobre el tramo recién recorrido. Como ese tramo
cruza ``max_lambda``, hay equilibrio y está cerca. Antes, el paso final se
decidía con el predictor: si la extrapolación tangente rebasaba
``max_lambda``, el solver fijaba λ = ``max_lambda`` sin recorrer la curva
(control de carga), lo que podía saltar un punto límite; y un paso de arco
que convergía por encima de ``max_lambda`` terminaba el trazado ahí (medido:
arco de von Mises con ``max_lambda = 0.39`` devolvía λ = 0.398).

Parada en el valor de una magnitud (ADR 0025)
---------------------------------------------
Con ``control`` (una combinación de grados de libertad, ``cᵀ·U``) y
``control_target``, el trazado se detiene cuando la magnitud alcanza ese
valor: el desplazamiento de un nodo, la deriva de un entrepiso. En el arco
cilíndrico ``control`` sólo **mide**; la restricción del paso sigue siendo
la norma del incremento. Si un paso convergido cruza el objetivo, su estado
no se consolida y el paso se repite desde el último estado convergido con
la restricción del control indirecto, ``cᵀ·ΔU`` = lo que falta, partiendo
de la interpolación lineal del tramo: el último punto coincide con el
objetivo a precisión de máquina. Es la misma llegada exacta que la de
``max_lambda`` (si un paso cruza los dos, se llega al que está antes en el
tramo). En el control indirecto, ``control`` es además la restricción.

Fin del trazado (ADR 0025)
--------------------------
El bucle termina al alcanzar ``max_lambda``, al alcanzar ``control_target``, al agotar
``max_steps`` o cuando un paso no converge ni con la longitud mínima. El
motivo queda en ``termination`` (:class:`~solidum.math.solvers.diagnostics.Termination`)
junto con ``lambda_final``, ``steps_done``, ``n_bisections`` y ``U_final``
(último estado convergido), también cuando ``solve`` termina con una
excepción; ``solve`` los reinicia al empezar. Agotar ``max_steps`` **no**
lanza excepción —trazar "hasta donde se pueda" es un uso legítimo del
arc-length—, pero deja ``termination = "max_steps"`` y ``completed = False``,
avisa con un ``WARNING`` y ``solidum.run`` lo refleja en
``SolveResult.converged`` y ``SolveResult.termination``. El fracaso lanza una
:class:`~solidum.math.solvers.diagnostics.SolverDivergedError` tipada con su
causa: :class:`~solidum.math.solvers.diagnostics.StepConstraintError` si la
restricción del paso no tiene solución (raíces imaginarias).

Registro (ADR 0025): una línea ``INFO`` por paso convergido; el detalle por
iteración, en ``DEBUG``; una bisección, en ``WARNING`` con su motivo; ``ERROR``
sólo cuando el análisis se detiene.

Etapas de carga (ADR 0026): ``solve(F_ref, u0=..., constant_load=...)`` parte
del estado ``u0`` con la carga constante ``F_c`` ya aplicada; el residuo es
``F_c + λ·F_ref − F_int`` con λ desde 0 y los valores impuestos crecen de los
de ``u0`` a los del modelo (:mod:`solidum.math.solvers.stages`). Es el empuje
lateral con la gravedad ya aplicada.
"""
from __future__ import annotations

import math

import numpy as np

from solidum.constants import (
    ARCLENGTH_DEFAULT_INITIAL_DLAMBDA,
    ARCLENGTH_MIN_DL_FACTOR,
    ZERO_TOL,
)
from solidum.math.convergence import ConvergenceCriterion
from solidum.math.solvers._shared import _log, domain_is_symmetric
from solidum.math.solvers.corrector import (
    CorrectionAborted,
    NewtonCorrector,
    default_calibration_scales,
)
from solidum.math.solvers.diagnostics import (
    IllPosedSystemError,
    SingularTangentError,
    SolverDivergedError,
    Termination,
)
from solidum.math.solvers.model_checks import (
    check_initial_tangent,
    ensure_statically_restrained,
)
from solidum.math.solvers.stages import Stage, check_initial_equilibrium, make_stage
from solidum.registry import SolverRegistry

# |cᵀ·du_t| por debajo de esta fracción de ‖c‖·‖du_t‖: la magnitud
# controlada no responde a la carga y ``Δλ`` no está determinado.
_INSENSITIVE_RTOL = 1.0e-12
# ``control_target`` se da por alcanzado a esta fracción de su valor: la
# llegada es exacta salvo el redondeo de acumular ``cᵀ·U`` paso a paso.
_TARGET_RTOL = 1.0e-12


def normalize_control(control, name: str) -> list[tuple]:
    """Lista de ``(nodo, dof, coef)`` validada en forma (no en existencia: los
    nodos y sus grados de libertad se comprueban al resolver). Cada término es
    una tupla ``(nodo, dof, coef)`` o un dict ``{node, dof, coef}``."""
    if isinstance(control, dict) or not control:
        raise ValueError(
            f"{name}: 'control' debe ser una lista no vacía de términos "
            f"(nodo, dof, coef), p. ej. [(b, 'ux', 1.0), (a, 'ux', -1.0)].")
    terms = []
    for k, term in enumerate(control):
        if isinstance(term, dict):
            unknown = set(term) - {'node', 'dof', 'coef'}
            if unknown or 'node' not in term or 'dof' not in term:
                raise ValueError(
                    f"{name}: el término {k} de 'control' debe tener las claves "
                    f"'node' y 'dof' (y opcionalmente 'coef'); recibido {term}.")
            node, dof, coef = term['node'], term['dof'], term.get('coef', 1.0)
        else:
            try:
                node, dof, coef = term
            except (TypeError, ValueError):
                raise ValueError(
                    f"{name}: el término {k} de 'control' debe ser (nodo, dof, "
                    f"coef) o {{node, dof, coef}}; recibido {term!r}.") from None
        coef = float(coef)
        if not isinstance(dof, str):
            raise ValueError(f"{name}: dof {dof!r} del término {k} debe ser un nombre ('ux', …).")
        if not math.isfinite(coef) or coef == 0.0:
            raise ValueError(f"{name}: coef={coef} del término {k} debe ser finito y no nulo.")
        terms.append((node, dof, coef))
    return terms


class _ArcProblem:
    """Paso del arc-length en el protocolo ``NewtonProblem``.

    El iterado ``x`` es ``(U_iter, λ_iter, dU_iter)`` con ``dU_iter`` el
    incremento acumulado del paso (``U_iter = U_current + dU_iter``); el
    estado es ``(K, F_int)``. La corrección ``dx`` es ``(dU_update, ddλ)``.

    ``mode`` fija la restricción que determina ``ddλ``:

    - ``"newton"`` — llegada a ``max_lambda``: ``λ`` fijo, corrección de
      Newton pura;
    - ``"cylindrical"`` — cuadrática de Crisfield ``‖ΔU‖² = dl²`` con
      selección de raíz por menor ángulo con el incremento previo;
    - ``"indirect"`` — restricción lineal ``cᵀ·ΔU = dl`` con el vector
      ``control_vector`` (control indirecto de desplazamiento, y llegada
      exacta a ``control_target``).

    Con ``max_lambda``, un paso con restricción que converge por encima de
    él no consolida el estado y marca ``overshoot``: el solver repite el
    paso con llegada exacta (``ArcLengthSolver._land``). Lo mismo si
    ``overshoot_check(U)`` es verdadero (el paso cruzó ``control_target``).
    """

    def __init__(self, assembler, U_current, lambda_curr, F_ext_ref, free_dofs,
                 *, mode: str, dl: float, max_lambda: float | None = None,
                 stage: Stage | None = None, control_vector: np.ndarray | None = None,
                 overshoot_check=None):
        self.assembler = assembler
        self.U_current = U_current
        self.lambda_curr = float(lambda_curr)
        self.F_ext_ref = F_ext_ref
        self.free_dofs = free_dofs
        self.mode = mode
        self.dl = float(dl)
        self.max_lambda = max_lambda
        self.c = control_vector
        self.overshoot_check = overshoot_check
        self.overshoot = False
        # Etapa de carga (ADR 0026): carga constante y valores impuestos de
        # partida. Sin etapa, nulos.
        staged = stage is not None and stage.staged
        self.F_const = stage.F_const if staged else None
        self.g_start = stage.g_start if staged else None
        self._F_calibration = stage.calibration_load(F_ext_ref) if staged else F_ext_ref
        # Equilibrio reducido ``TᵀR = 0`` (ADR 0004): con una restricción
        # multipunto ``R[libres]`` no se anula en la solución.
        self._T = assembler.constraint_set.build(assembler.domain.total_dofs)[0]

    # --- protocolo ----------------------------------------------------------

    def assemble(self, x):
        return self.assembler.assemble_non_linear_system(x[0])

    def residual(self, x, state):
        R = x[1] * self.F_ext_ref - state[1]
        return R if self.F_const is None else R + self.F_const

    def residual_norm(self, R):
        return float(np.linalg.norm(self._T.T @ R))

    def calibration_scales(self, x, state):
        # La carga de referencia es la escala natural del arc-length (λ
        # varía, F_ext_ref es fijo).
        return default_calibration_scales(self._F_calibration, state[1], state[0])

    def reference_force(self, x, state):
        F_ext = x[1] * self.F_ext_ref
        if self.F_const is not None:
            F_ext = F_ext + self.F_const
        return max(float(np.linalg.norm(F_ext)), float(np.linalg.norm(state[1])))

    def x_norm(self, x):
        return float(np.linalg.norm(x[0]))

    def log_context(self, x):
        return f" | lam={x[1]:.4f}"

    def correction(self, x, state, R, solve):
        U_iter, lambda_iter, dU_iter = x
        K_global = state[0]
        K_t_red, F_t_red, T_t, g_t = self.assembler.reduce(
            K_global, self.F_ext_ref.copy(), g_start=self.g_start)
        K_red, R_red, T_R, g_R = self.assembler.reduce(
            K_global, R, U_current=U_iter, load_factor=lambda_iter, g_start=self.g_start,
        )
        du_R = self.assembler.expand(solve(K_red, R_red), T_R, g_R)
        du_t = self.assembler.expand(solve(K_t_red, F_t_red), T_t, g_t)
        return self.constraint(x, du_R, du_t)

    def apply(self, x, dx, alpha):
        U_iter, lambda_iter, dU_iter = x
        dU_update, ddlambda = dx
        dU_new = dU_iter + alpha * dU_update
        return ((self.U_current + dU_new, lambda_iter + alpha * ddlambda, dU_new),
                float(np.linalg.norm(alpha * dU_update)))

    def on_converged(self, x, state):
        # Paso con restricción que acaba por encima de max_lambda: no se
        # consolida (el trial queda en U_iter y el siguiente ensamblaje lo
        # recalcula desde el estado consolidado); el solver aterriza.
        if self.mode != "newton" and (
                (self.max_lambda is not None and x[1] > self.max_lambda + ZERO_TOL)
                or (self.overshoot_check is not None and self.overshoot_check(x[0]))):
            self.overshoot = True
            return
        # Estado trial del ensamblaje en U_iter: coherente con lo guardado.
        self.assembler.commit_all_states()

    # --- restricción --------------------------------------------------------

    def constraint(self, x, du_R, du_t):
        """``(dU_update, ddλ)`` según ``mode``."""
        _, _, dU_iter = x
        if self.mode == "newton":
            return du_R, 0.0
        if self.mode == "cylindrical":
            dU_new = dU_iter + du_R
            a = np.dot(du_t, du_t)
            b = 2.0 * np.dot(dU_new, du_t)
            c = np.dot(dU_new, dU_new) - self.dl ** 2
            det = b ** 2 - 4.0 * a * c
            if det < 0:
                raise CorrectionAborted("Raíces imaginarias. La solución diverge del arco.")
            ddl1 = (-b + np.sqrt(det)) / (2.0 * a)
            ddl2 = (-b - np.sqrt(det)) / (2.0 * a)
            # Raíz que produce el menor ángulo con el incremento previo.
            theta1 = np.dot(dU_iter, dU_new + ddl1 * du_t)
            theta2 = np.dot(dU_iter, dU_new + ddl2 * du_t)
            ddlambda = ddl1 if theta1 > theta2 else ddl2
            return du_R + ddlambda * du_t, ddlambda
        if self.mode == "indirect":
            denom = float(self.c @ du_t)
            if abs(denom) <= _INSENSITIVE_RTOL * np.linalg.norm(self.c) * np.linalg.norm(du_t):
                raise CorrectionAborted(
                    "cᵀ·du_t ≈ 0: la magnitud controlada no responde a la carga en "
                    "este estado; el paso no determina Δλ.")
            ddlambda = (self.dl - float(self.c @ (dU_iter + du_R))) / denom
            return du_R + ddlambda * du_t, ddlambda
        raise ValueError(f"_ArcProblem: modo desconocido {self.mode!r}.")


@SolverRegistry.register
class ArcLengthSolver:
    """
    Solucionador no lineal con Método de Longitud de Arco Cilíndrico (Crisfield).
    Permite trazar curvas de equilibrio con fenómenos de snap-through y snap-back
    variando simultáneamente los desplazamientos y la carga externa.

    Primer paso
    -----------
    initial_dlambda : float, opcional
        Fracción de la carga de referencia del primer paso (adimensional).
        Por omisión ``ARCLENGTH_DEFAULT_INITIAL_DLAMBDA`` (0.1). El solver
        deriva la longitud de arco Δl₁ = Δλ₁·‖K⁻¹·F_ref‖ con el predictor
        elástico del primer paso.
    initial_dl : float, opcional
        Longitud de arco del primer paso, en unidades de desplazamiento.
        Excluyente con ``initial_dlambda``.

    Parada en el valor de una magnitud
    ----------------------------------
    control : list, opcional
        Magnitud vigilada, ``cᵀ·U = Σ coef·u(nodo, dof)``. Cada término es una
        tupla ``(nodo, dof, coef)`` o un dict ``{node, dof, coef}`` (``coef``
        vale 1 si se omite); ``nodo`` es un :class:`Node` o su id. La deriva
        entre los nodos ``a`` y ``b`` en ``x`` es ``[(b, 'ux', 1.0), (a, 'ux',
        -1.0)]``. Los grados de libertad no pueden estar prescritos. En el arco
        cilíndrico sólo se mide (la línea ``INFO`` de cada paso la muestra);
        en :class:`IndirectDisplacementSolver` es además la restricción.
    control_target : float, opcional
        Valor **absoluto** de ``cᵀ·U`` en el que se detiene el trazado
        (requiere ``control``). El paso que lo cruza se repite con llegada
        exacta (ver el docstring del módulo): deja ``termination =
        "control_target"`` y ``reached_control_target = True``. Puede
        combinarse con ``max_lambda``: se detiene en lo que llegue primero
        (con un objetivo, conviene un ``max_lambda`` holgado). En una etapa
        (ADR 0026) la magnitud parte de su valor en ``u0``.

    Atributos de estado tras ``solve`` (ADR 0025)
    ---------------------------------------------
    ``solve`` los reinicia al empezar y los deja también cuando termina con
    una excepción.

    termination : Termination or None
        Motivo de la parada: ``"max_lambda"``, ``"control_target"`` (variantes
        con objetivo), ``"max_steps"`` o ``"diverged"``. ``None`` antes de
        trazar o si ``solve`` falló antes del primer paso (modelo mal
        planteado).
    completed : bool
        ``True`` si el trazado llegó a donde se pidió (``max_lambda`` u
        objetivo); ``False`` con ``max_steps`` agotado o fracaso.
    lambda_final : float
        Factor de carga del último paso convergido.
    reached_max_lambda : bool
        ``True`` si el trazado alcanzó ``max_lambda``.
    reached_control_target : bool
        ``True`` si el trazado alcanzó ``control_target``.
    control_final : float
        Valor de la magnitud ``cᵀ·U`` en el último paso convergido (0 sin
        ``control``).
    steps_done : int
        Número de pasos convergidos.
    n_bisections : int
        Reducciones de la longitud de paso por pasos que no convergieron.
    U_final : np.ndarray or None
        Desplazamientos del último paso convergido (también tras un fracaso).
    dl_reference : float
        Longitud de arco del primer paso (Δl₁) usada en la corrida; es la
        referencia de ``dl_max_factor`` y del umbral de aborto.
    """

    PIPELINE_KIND = "static"

    def __init__(self, assembler, convergence: ConvergenceCriterion | None = None, max_iter=20, max_lambda=1.0,
                 initial_dl: float | None = None, max_steps=100,
                 dl_grow_factor=1.5, dl_max_factor=5.0, dl_shrink_factor=0.6,
                 dl_grow_iter_threshold=4, dl_shrink_iter_threshold=8,
                 linear_algebra: str = "auto", *, initial_dlambda: float | None = None,
                 control=None, control_target: float | None = None):
        self.assembler = assembler
        name = type(self).__name__
        # Magnitud vigilada y objetivo (ADR 0025). ``_c_raw`` es el vector de
        # la magnitud tal como se declaró; ``_c``, el mismo orientado hacia el
        # objetivo (``_direction`` = ±1). Se construyen al resolver, sobre la
        # numeración de ecuaciones vigente.
        self.control = None if control is None else normalize_control(control, name)
        if control_target is not None:
            if self.control is None:
                raise ValueError(
                    f"{name}: control_target necesita 'control', la combinación de "
                    f"grados de libertad cuyo valor se vigila.")
            control_target = float(control_target)
            if not math.isfinite(control_target):
                raise ValueError(
                    f"{name}: control_target={control_target} debe ser un número finito.")
        self.control_target: float | None = control_target
        self._c: np.ndarray | None = None
        self._c_raw: np.ndarray | None = None
        self._direction: float = 1.0
        self._target_scale: float = 0.0
        self.reached_control_target: bool = False
        self.control_final: float = 0.0
        self._landed_on: str | None = None
        # Política de convergencia (ADR 0007). Compartida con NonlinearSolver:
        # cambiar la política aquí llega automáticamente al arc-length.
        self.convergence = convergence if convergence is not None else ConvergenceCriterion()
        self.max_iter = max_iter
        self.max_lambda = max_lambda
        self.initial_dl, self.initial_dlambda = self._first_step_spec(initial_dl, initial_dlambda)
        self.dl_reference: float | None = None
        self.max_steps = max_steps
        # Factores de auto-ajuste de la longitud de arco:
        #   Si converge en < dl_grow_iter_threshold iter → ampliar dl × dl_grow_factor (max: initial_dl × dl_max_factor)
        #   Si converge en > dl_shrink_iter_threshold iter → reducir dl × dl_shrink_factor
        self.dl_grow_factor = dl_grow_factor
        self.dl_max_factor = dl_max_factor
        self.dl_shrink_factor = dl_shrink_factor
        self.dl_grow_iter_threshold = dl_grow_iter_threshold
        self.dl_shrink_iter_threshold = dl_shrink_iter_threshold
        self.linear_algebra = linear_algebra
        # Corrector compartido (ADR 0015). Régimen postcrítico: K_t puede
        # ser indefinida → no asumir SPD. Si el usuario fuerza
        # ``linear_algebra: cholesky`` y K_t se vuelve indefinida, el
        # corrector degrada a LU para que el override no rompa el análisis.
        self.corrector = NewtonCorrector(
            self.convergence,
            max_iter=self.max_iter,
            is_symmetric=domain_is_symmetric(assembler.domain),
            is_positive_definite=False,
            linear_algebra=self.linear_algebra,
            near_nullspace=assembler.near_nullspace,
        )
        # Estado del último trazado (ver docstring de la clase; ADR 0025).
        self.termination: Termination | None = None
        self.lambda_final: float = 0.0
        self.reached_max_lambda: bool = False
        self.steps_done: int = 0
        self.n_bisections: int = 0
        self.U_final: np.ndarray | None = None
        # Residuo relativo de equilibrio al empezar una etapa (ADR 0026) y la
        # etapa en curso (la leen el predictor, el paso y la llegada).
        self.initial_residual: float | None = None
        self._stage: Stage | None = None

    @classmethod
    def _first_step_spec(cls, initial_dl, initial_dlambda):
        """Valida la declaración del primer paso: ``(initial_dl, initial_dlambda)``
        con uno de los dos en ``None``."""
        name = cls.__name__
        if initial_dl is not None and initial_dlambda is not None:
            raise ValueError(
                f"{name}: declarar sólo uno de 'initial_dlambda' (fracción de la carga "
                f"de referencia en el primer paso, recomendado) o 'initial_dl' "
                f"(longitud de arco en unidades de desplazamiento)."
            )
        if initial_dl is not None:
            if not float(initial_dl) > 0.0:
                raise ValueError(f"{name}: initial_dl={initial_dl} debe ser > 0.")
            return float(initial_dl), None
        dlam = ARCLENGTH_DEFAULT_INITIAL_DLAMBDA if initial_dlambda is None else float(initial_dlambda)
        if not dlam > 0.0:
            raise ValueError(f"{name}: initial_dlambda={initial_dlambda} debe ser > 0.")
        return None, dlam

    # ------------------------------------------------------------------
    # Puntos de extensión de la restricción (variantes que sólo cambian la
    # magnitud controlada, p. ej. IndirectDisplacementSolver)
    # ------------------------------------------------------------------

    def _control_measure(self, du: np.ndarray) -> float:
        """Magnitud que la restricción del paso fija en ``Δl``, evaluada en
        un incremento ``du``: aquí su norma euclídea (arco cilíndrico)."""
        return float(np.linalg.norm(du))

    def _predictor_dlambda(self, du_t: np.ndarray, sign: float, dl: float) -> float:
        """Incremento de carga del predictor tangente ``ΔU = Δλ·du_t``."""
        return sign * dl / (self._control_measure(du_t) + ZERO_TOL)

    def _step_problem(self, U_current, lambda_curr, F_ext_ref, free_dofs, dl):
        """Problema de Newton del paso con la restricción del solver."""
        return _ArcProblem(
            self.assembler, U_current, lambda_curr, F_ext_ref, free_dofs,
            mode="cylindrical", dl=dl, max_lambda=self.max_lambda, stage=self._stage,
            overshoot_check=(self._crossed_target if self.control_target is not None
                             else None),
        )

    def _on_solve_start(self) -> None:
        """Preparación dependiente de la numeración de ecuaciones: el vector
        de la magnitud vigilada, si hay ``control``."""
        if self.control is not None:
            self._setup_control()

    def _setup_control(self) -> None:
        """Vector ``c`` de la magnitud vigilada sobre la numeración vigente.

        ``c`` actúa sobre el incremento completo ``ΔU = T·Δu_libre + Δλ·g``
        (ADR 0004), que incluye los esclavos reconstruidos por ``T``: el
        esclavo de una restricción multipunto con maestros (fila de ``T`` no
        vacía) se admite, porque su incremento lo fijan sus maestros. Un grado
        de libertad prescrito (fila de ``T`` vacía: apoyo declarado, o esclavo
        cuyos maestros están todos prescritos, que hereda el apoyo por el
        cierre transitivo) se rechaza: su incremento es ``Δλ·g_s``, no depende
        del equilibrio. También se rechaza una combinación cuya proyección
        ``Tᵀ·c`` sobre los grados de libertad libres se anula (p. ej.
        ``u_s − u_m`` con ``u_s = u_m + g``)."""
        name = type(self).__name__
        domain = self.assembler.domain
        ndof = domain.total_dofs
        T = self.assembler.constraint_set.build(ndof)[0].tocsr()
        c = np.zeros(ndof)
        for node, dof, coef in self.control:
            if not hasattr(node, 'dofs'):
                node_obj = domain.get_node(node)
                if node_obj is None:
                    raise ValueError(
                        f"{name}: el nodo {node!r} de 'control' no existe en el modelo.")
                node = node_obj
            if dof not in node.dofs:
                raise ValueError(
                    f"{name}: el nodo {node.id} no tiene el grado de libertad {dof!r}; "
                    f"tiene {sorted(node.dofs)}.")
            index = node.dofs[dof]
            row = T.data[T.indptr[index]:T.indptr[index + 1]]
            if not np.any(row):
                raise ValueError(
                    f"{name}: el grado de libertad {dof!r} del nodo {node.id} está "
                    f"prescrito (apoyo, o esclavo de una restricción multipunto cuyos "
                    f"maestros están todos prescritos): su incremento no depende del "
                    f"equilibrio y no puede formar parte de la magnitud controlada. "
                    f"Para controlar el desplazamiento de un nodo respecto de un apoyo "
                    f"fijo basta con el nodo solo.")
            c[index] += coef
        # Proyección sobre los grados de libertad libres: cᵀ·ΔU = (Tᵀc)ᵀ·Δu_libre
        # + Δλ·cᵀg. Sin restricciones multipunto, Tᵀc = 0 sólo si c = 0.
        c_free = T.T @ c
        scale = abs(T).T @ np.abs(c)
        if np.linalg.norm(c_free) <= _INSENSITIVE_RTOL * np.linalg.norm(scale):
            raise ValueError(
                f"{name}: los términos de 'control' se anulan entre sí (Tᵀ·c = 0: "
                f"tras aplicar las restricciones multipunto, la combinación no "
                f"depende de ningún grado de libertad libre).")
        self._c_raw = c
        self._direction = 1.0
        self._c = c
        if self.control_target is not None:
            # Sentido del trazado: del valor de partida (u0 de la etapa, o
            # cero) hacia el objetivo (ADR 0025, ampliado para etapas).
            start = float(c @ self._stage.U0)
            self._target_scale = max(abs(self.control_target), abs(start))
            gap = self.control_target - start
            if abs(gap) <= _TARGET_RTOL * self._target_scale:
                _log.warning(
                    f"{name}: control_target = {self.control_target:.6e} es el valor "
                    f"de partida de la magnitud controlada; no hay nada que trazar.")
            self._direction = 1.0 if gap >= 0.0 else -1.0
            self._c = self._direction * c

    def _first_dl(self, du_t: np.ndarray) -> float:
        """Longitud de arco del primer paso a partir del predictor elástico
        ``du_t = K⁻¹·F_ref`` (misma magnitud que la restricción del paso)."""
        if self.control_target is not None and self._c is not None:
            m = float(self._c @ du_t)
            if np.isfinite(m) and m < 0.0:
                _log.warning(
                    f"{type(self).__name__}: la carga de referencia aleja la magnitud "
                    f"vigilada de control_target; el arco avanza con λ creciente y "
                    f"puede no alcanzarlo (lo alcanzaría tras un punto límite).")
        if self.initial_dl is not None:
            dl = self.initial_dl
        else:
            norm = self._control_measure(du_t)
            if not np.isfinite(norm) or norm == 0.0:
                raise ValueError(
                    f"{type(self).__name__}: la carga de referencia no produce "
                    f"desplazamiento (‖K⁻¹·F_ref‖ = {norm}); no hay escala para el "
                    f"primer paso. Revisar F_ref o declarar 'initial_dl'."
                )
            dl = self.initial_dlambda * norm
            _log.info(
                f"  Primer paso: Δλ₁ = {self.initial_dlambda:g} de la carga de "
                f"referencia ⇒ Δl₁ = {dl:.4e}"
            )
        self.dl_reference = float(dl)
        return float(dl)

    def _land(self, U_current, lambda_curr, x_over, F_ext_ref, free_dofs):
        """Llegada exacta a ``max_lambda`` o a ``control_target`` dentro del
        tramo recién recorrido.

        ``x_over`` es el paso con restricción que convergió más allá (sin
        consolidar). Si cruzó los dos, se llega al que está antes en el tramo.
        A ``max_lambda``: Newton puro con λ fijo. A ``control_target``:
        restricción indirecta ``cᵀ·ΔU`` = lo que falta, con λ libre. En los dos
        casos desde el último estado consolidado, arrancando de la
        interpolación lineal del tramo. Devuelve el ``CorrectorResult`` y deja
        en ``_landed_on`` a cuál se llegó.
        """
        U_over, lambda_over, _ = x_over
        t_lambda = (math.inf if lambda_over <= self.max_lambda + ZERO_TOL
                    else (self.max_lambda - lambda_curr) / (lambda_over - lambda_curr))
        t_target = math.inf
        if self.control_target is not None and self._crossed_target(U_over):
            before, after = self._remaining(U_current), self._remaining(U_over)
            t_target = before / (before - after)
        if t_target < t_lambda:
            self._landed_on = "control_target"
            dU0 = t_target * (U_over - U_current)
            lambda0 = lambda_curr + t_target * (lambda_over - lambda_curr)
            _log.debug(
                f"  {self._step_summary(U_over).strip(' |')} rebasa control_target = "
                f"{self.control_target:.6e}: llegada exacta dentro del tramo recorrido.")
            problem = _ArcProblem(
                self.assembler, U_current, lambda_curr, F_ext_ref, free_dofs,
                mode="indirect", dl=self._remaining(U_current), stage=self._stage,
                control_vector=self._c,
            )
        else:
            self._landed_on = "max_lambda"
            if not math.isfinite(t_lambda):
                # Llamada sin cruce detectado (una variante que aterriza por su
                # cuenta): la interpolación de siempre hacia max_lambda.
                t_lambda = (self.max_lambda - lambda_curr) / (lambda_over - lambda_curr)
            dU0 = t_lambda * (U_over - U_current)
            lambda0 = self.max_lambda
            _log.debug(
                f"  λ = {lambda_over:.4f} rebasa max_lambda = {self.max_lambda:.4f}: "
                f"llegada exacta dentro del tramo recorrido."
            )
            problem = _ArcProblem(
                self.assembler, U_current, lambda_curr, F_ext_ref, free_dofs,
                mode="newton", dl=0.0, stage=self._stage,
            )
        return self.corrector.run(
            problem, (U_current + dU0, lambda0, dU0), check_initial=True,
            initial_delta_norm=float(np.linalg.norm(dU0)),
        )

    # ------------------------------------------------------------------
    # Estado del trazado (ADR 0025)
    # ------------------------------------------------------------------

    def _begin_trace(self) -> None:
        """Reinicia el estado del trazado al empezar ``solve``: un solver
        reutilizado no debe mostrar los valores de la corrida anterior."""
        self.termination = None
        self.lambda_final = 0.0
        self.reached_max_lambda = False
        self.reached_control_target = False
        self.control_final = 0.0
        self.steps_done = 0
        self.n_bisections = 0
        self.U_final = None
        self.dl_reference = None
        self.initial_residual = None
        self._stage = None
        self._landed_on = None

    @property
    def completed(self) -> bool:
        """``True`` si el trazado llegó a donde se pidió (``max_lambda`` o el
        objetivo de la variante); ``False`` con ``max_steps`` agotado, tras un
        fracaso o antes de trazar."""
        return self.termination in (Termination.MAX_LAMBDA, Termination.CONTROL_TARGET)

    def _control_value(self, U) -> float:
        """Magnitud vigilada tal como se declaró, ``cᵀ·U``."""
        return float(self._c_raw @ U)

    def _remaining(self, U) -> float:
        """Lo que falta hasta ``control_target``, medido en el sentido del
        trazado (negativo si se rebasó)."""
        return self._direction * (self.control_target - self._control_value(U))

    def _target_reached(self, U: np.ndarray | None) -> bool:
        """``True`` si ``U`` está en ``control_target`` (a la tolerancia de
        redondeo); siempre ``False`` sin objetivo."""
        if self.control_target is None or U is None or self._c_raw is None:
            return False
        return self._remaining(U) <= _TARGET_RTOL * self._target_scale

    def _crossed_target(self, U) -> bool:
        """``True`` si ``U`` rebasó ``control_target`` más allá del redondeo:
        el paso que llegó a ``U`` debe repetirse con llegada exacta."""
        return self._remaining(U) < -_TARGET_RTOL * self._target_scale

    def _limit_step(self, dl: float, U: np.ndarray) -> float:
        """Longitud del paso siguiente: una variante con restricción sobre la
        magnitud vigilada la recorta para llegar exactamente al objetivo
        (aquí, sin cambio: el cruce se resuelve con ``_land``)."""
        return dl

    def _step_summary(self, U: np.ndarray) -> str:
        """Texto extra de la línea ``INFO`` del paso: la magnitud vigilada."""
        if self._c_raw is None:
            return ""
        return f" | control = {self._control_value(U):.6e}"

    def _record(self, lambda_curr: float, steps: int, U: np.ndarray | None) -> None:
        """Hasta dónde llegó el trazado: último paso convergido."""
        self.lambda_final = float(lambda_curr)
        self.steps_done = int(steps)
        self.U_final = None if U is None else np.array(U, dtype=float, copy=True)
        if U is not None and self._c_raw is not None:
            self.control_final = self._control_value(U)

    def _finish(self, lambda_curr: float, step: int, U: np.ndarray | None = None) -> None:
        """Registra el estado final de un trazado que terminó sin fracaso y
        avisa si quedó incompleto."""
        self._record(lambda_curr, step, U)
        self.reached_max_lambda = lambda_curr >= self.max_lambda - ZERO_TOL
        if self._target_reached(U):
            self.termination = Termination.CONTROL_TARGET
        elif self.reached_max_lambda:
            self.termination = Termination.MAX_LAMBDA
        else:
            self.termination = Termination.MAX_STEPS
        self.reached_control_target = self.termination is Termination.CONTROL_TARGET
        name = type(self).__name__
        if self.termination is Termination.MAX_STEPS:
            _log.warning(
                f"{name}: trazado detenido en λ={lambda_curr:.4f} "
                f"< max_lambda={self.max_lambda:.4f} tras agotar max_steps="
                f"{self.max_steps} (termination = 'max_steps'). El resultado "
                f"corresponde al último paso convergido; aumente max_steps o el "
                f"primer paso (initial_dlambda) para completar."
            )
            return
        bisecciones = f", {self.n_bisections} bisección(es)" if self.n_bisections else ""
        _log.info(
            f"{name}: trazado completo (termination = '{self.termination}'): "
            f"λ = {lambda_curr:.4f} en {step} pasos{bisecciones}."
        )

    def _fail(self, error: SolverDivergedError, lambda_curr: float, steps: int,
              U: np.ndarray | None) -> SolverDivergedError:
        """Deja el estado de un trazado que fracasa (``termination =
        "diverged"`` y hasta dónde llegó) y devuelve la excepción tipada, lista
        para lanzar."""
        self._record(lambda_curr, steps, U)
        self.reached_max_lambda = False
        self.termination = Termination.DIVERGED
        _log.error(
            f"{type(self).__name__}: análisis detenido ({error.mode}). Último paso "
            f"convergido: λ = {lambda_curr:.4f} tras {steps} pasos."
        )
        return error

    # ------------------------------------------------------------------
    # Predictor tangente (común a las variantes)
    # ------------------------------------------------------------------

    def _tangent_predictor(self, U_current, F_ext_ref, delta_U_step, step):
        """Ensambla en ``U_current``, calibra el criterio si hace falta y
        devuelve ``(du_t, sign)`` con ``du_t = K⁻¹·F_ref`` y el sentido de
        avance (para no regresar por donde se vino). Devuelve ``None`` si la
        tangente es singular (el llamador biseca la longitud de paso)."""
        K_global, F_int_global = self.assembler.assemble_non_linear_system(U_current)
        stage = self._stage
        staged = stage is not None and stage.staged
        if not self.convergence.is_calibrated:
            self.convergence.calibrate(*default_calibration_scales(
                stage.calibration_load(F_ext_ref) if staged else F_ext_ref,
                F_int_global, K_global))
        K_t_red, F_t_red, T_t, g_t = self.assembler.reduce(
            K_global, F_ext_ref.copy(), g_start=stage.g_start if staged else None)
        try:
            du_t_red = self.corrector.solve(K_t_red, F_t_red)
        except IllPosedSystemError:
            # Veredicto de la red de seguridad sobre la tangente inicial
            # (ADR 0019): el modelo está mal planteado; bisecar no lo arregla.
            raise
        except RuntimeError:
            return None
        du_t = self.assembler.expand(du_t_red, T_t, g_t)
        sign = 1.0
        if step > 1 and np.dot(delta_U_step, du_t) < 0:
            sign = -1.0
        return du_t, sign

    def solve(self, F_ext_ref: np.ndarray, step_callback=None, *,
              u0: np.ndarray | None = None,
              constant_load: np.ndarray | None = None) -> np.ndarray:
        """Traza la trayectoria de equilibrio con la carga ``λ·F_ext_ref``.

        Parameters
        ----------
        F_ext_ref
            Carga de referencia, la que crece con λ.
        step_callback
            ``callback(step, U, λ)`` tras cada paso convergido.
        u0, constant_load
            Etapa de carga (ADR 0026): estado de partida (los desplazamientos
            de la etapa anterior; el estado interno de los elementos lo guarda
            el ensamblador) y carga que actúa completa durante toda la etapa,
            p. ej. la gravedad ya aplicada. El residuo es
            ``F_c + λ·F_ext_ref − F_int``, con λ desde 0.
        """
        self._begin_trace()
        domain = self.assembler.domain
        if domain.total_dofs == 0:
            domain.generate_equation_numbers()
        ndof = domain.total_dofs
        self._stage = make_stage(self.assembler, u0, constant_load, type(self).__name__)
        U_current = self._stage.U0.copy()

        lambda_curr = 0.0
        step = 0
        steps_done = 0
        dl = None          # Δl₁: se fija con el predictor del primer paso
        dl_ref = None

        delta_U_step = np.zeros(ndof)  # Historial del incremento del paso para guiar el arco

        _log.info("--- INICIANDO SOLVER NO LINEAL (MÉTODO ARC-LENGTH) ---")
        # Red de seguridad, capa 1 (ADR 0019): un mecanismo rígido no tiene
        # equilibrio estático; se rechaza antes de iterar con el movimiento
        # libre en palabras, en vez de agotar bisecciones.
        ensure_statically_restrained(self.assembler, type(self).__name__)
        # Capa 3 sobre la tangente inicial (ADR 0019, 2026-09-30): el
        # predictor del primer paso factoriza K(U = 0); sus pivotes nulos
        # delatan un mecanismo interno aunque la carga no lo active.
        check_initial_tangent(self.assembler, self.corrector, type(self).__name__)
        # Etapa (ADR 0026): debe empezar en equilibrio con la carga constante.
        self.initial_residual = check_initial_equilibrium(
            self.assembler, self._stage, self.convergence, F_ext_ref, type(self).__name__)

        cs = self.assembler.constraint_set
        free_dofs = cs.free_dofs(ndof)
        self._on_solve_start()

        while (lambda_curr < self.max_lambda and step < self.max_steps
               and not self._target_reached(U_current)):
            step += 1

            # Gancho de inicio de paso (ADR 0020, P5). Evaluado con el estado
            # convergido del paso anterior para evitar chattering.
            self.assembler.prepare_all_steps(U_current)

            # --- 1. PREDICTOR ---
            pred = self._tangent_predictor(U_current, F_ext_ref, delta_U_step, step)
            if pred is None:
                if dl is None:
                    # Sin paso previo no hay longitud que reducir, y repetir
                    # desde el mismo estado daría la misma tangente.
                    raise self._fail(SingularTangentError(
                        last_residual=float("nan"), last_delta=0.0,
                        last_load_factor=lambda_curr, n_bisections=self.n_bisections,
                        extra_message=("la tangente del estado inicial es singular: "
                                       "el predictor no tiene dirección de carga"),
                    ), lambda_curr, steps_done, U_current)
                dl /= 2.0
                self.n_bisections += 1
                _log.warning(f"[PASO {step}] tangente singular en el predictor; "
                             f"bisección: dl = {dl:.4e}")
                if dl < ARCLENGTH_MIN_DL_FACTOR * dl_ref:
                    raise self._fail(SingularTangentError(
                        last_residual=float("nan"), last_delta=0.0,
                        last_load_factor=lambda_curr, n_bisections=self.n_bisections,
                        extra_message=(
                            f"paso {step}: la tangente del predictor sigue singular "
                            f"con dl = {dl:.2e} (< {ARCLENGTH_MIN_DL_FACTOR:g}·Δl₁)"),
                    ), lambda_curr, steps_done, U_current)
                continue
            du_t, sign = pred
            if dl is None:
                dl = dl_ref = self._first_dl(du_t)
            dl_step = self._limit_step(dl, U_current)
            # Pseudo-tiempo del paso = longitud de arco (o incremento de la
            # magnitud controlada), conocida tras el predictor (ADR 0034).
            self.assembler.set_pseudo_time_increment(dl_step)
            _log.debug(f"[PASO {step}] longitud de paso dl = {dl_step:.4e}")

            dlambda = self._predictor_dlambda(du_t, sign, dl_step)
            dU_iter = dlambda * du_t
            x0 = (U_current + dU_iter, lambda_curr + dlambda, dU_iter)

            # --- 2. CORRECTOR ITERATIVO (restricción del solver) ---
            problem = self._step_problem(U_current, lambda_curr, F_ext_ref, free_dofs, dl_step)
            res = self.corrector.run(
                problem, x0, check_initial=True,
                initial_delta_norm=float(np.linalg.norm(dU_iter)),
            )
            # Paso que cruza max_lambda: llegada exacta dentro del tramo.
            landed = False
            if res.converged and problem.overshoot:
                res = self._land(U_current, lambda_curr, res.x, F_ext_ref, free_dofs)
                landed = True

            if res.converged:
                U_current, lambda_curr, delta_U_step = res.x
                steps_done += 1
                _log.info(
                    f"[PASO {step}] λ = {lambda_curr:.4f} | dl = {dl_step:.3e} | "
                    f"{res.n_solves} iter{self._step_summary(U_current)}"
                    + (f" | llegada exacta a {self._landed_on}" if landed else "")
                )
                # Auto-ajuste de longitud de arco
                if res.n_solves < self.dl_grow_iter_threshold:
                    dl = min(dl * self.dl_grow_factor, dl_ref * self.dl_max_factor)
                elif res.n_solves > self.dl_shrink_iter_threshold:
                    dl *= self.dl_shrink_factor

                if step_callback:
                    step_callback(step, U_current, lambda_curr)
                continue

            dl = 0.5 * dl_step
            self.n_bisections += 1
            reason = res.failure_reason()
            _log.warning(f"[PASO {step}] sin convergencia ({reason}); "
                         f"bisección: dl = {dl:.4e}")
            if dl < ARCLENGTH_MIN_DL_FACTOR * dl_ref:
                raise self._fail(res.divergence_error(
                    last_load_factor=float(res.x[1]),
                    n_bisections=self.n_bisections,
                    extra_message=(
                        f"paso {step}: la longitud de paso bajó a {dl:.2e} "
                        f"(< {ARCLENGTH_MIN_DL_FACTOR:g}·Δl₁ = "
                        f"{ARCLENGTH_MIN_DL_FACTOR * dl_ref:.2e}) sin converger; "
                        f"último paso convergido: λ = {lambda_curr:.4f} tras "
                        f"{steps_done} pasos"),
                ), lambda_curr, steps_done, U_current)

        self._finish(lambda_curr, steps_done, U_current)
        return U_current
