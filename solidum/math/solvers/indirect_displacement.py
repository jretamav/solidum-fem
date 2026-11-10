# solidum_fem/solidum/math/solvers/indirect_displacement.py
"""``IndirectDisplacementSolver`` — control indirecto de desplazamiento.

Variante de :class:`~solidum.math.solvers.arclength.ArcLengthSolver` (de Borst
1987) que sólo cambia la restricción del paso: en vez de fijar la norma del
incremento de **todos** los grados de libertad (arco cilíndrico), fija el
incremento de una **combinación elegida** de ellos,

    cᵀ·ΔU = Δl,

por ejemplo el alargamiento de la zona que se ablanda o la apertura de la boca
de una fisura (CMOD), como se controla un ensayo de fractura en laboratorio.
La carga ``λ`` queda libre y puede bajar: sigue ramas con retroceso
(*snap-back*) mientras la magnitud controlada crezca monótonamente.

Por qué hace falta: en un retroceso por localización el incremento de la
mayoría de los grados de libertad (la parte que se descarga) apunta casi al
revés que el del paso anterior, y la selección de raíz por menor ángulo del
arco cilíndrico elige la rama equivocada o no converge. La restricción
indirecta es lineal en ``ddλ``: no hay raíces que elegir. Ver
``docs/specs/IndirectDisplacementSolver.md``.

Parada en un valor objetivo (ADR 0025): con ``control_target`` el trazado se
detiene cuando la magnitud controlada ``cᵀ·U`` alcanza ese valor. Como la
restricción fija exactamente ``cᵀ·ΔU = Δl`` en cada paso, basta recortar el
último paso a lo que falta: la llegada es exacta, a precisión de máquina, con
cualquier tamaño de paso. Es la parada natural de un ensayo o de un análisis de
empuje lateral, que terminan en un desplazamiento y no en un factor de carga
(en una meseta plástica λ casi no cambia).
"""
from __future__ import annotations

import numpy as np

from solidum.math.convergence import ConvergenceCriterion
from solidum.math.solvers._shared import _log
from solidum.math.solvers.arclength import _INSENSITIVE_RTOL, ArcLengthSolver, _ArcProblem
from solidum.registry import SolverRegistry


@SolverRegistry.register
class IndirectDisplacementSolver(ArcLengthSolver):
    """Control indirecto de desplazamiento (de Borst 1987).

    Parameters
    ----------
    control : list
        Combinación de grados de libertad controlada, ``Σ coef·u(nodo, dof)``.
        Cada término es una tupla ``(nodo, dof, coef)`` o un dict
        ``{node, dof, coef}`` (``coef`` vale 1 si se omite); ``nodo`` es un
        :class:`Node` o su id. El alargamiento entre los nodos ``a`` y ``b``
        en ``x`` es ``[(b, 'ux', 1.0), (a, 'ux', -1.0)]``. Los grados de
        libertad no pueden estar prescritos (con apoyo, o esclavos de una
        restricción multipunto cuyos maestros están todos prescritos); el
        esclavo de una restricción multipunto con maestros libres sí se
        admite, porque su incremento lo fijan sus maestros. La combinación
        debe crecer al aplicar la carga de referencia.
    initial_dlambda, initial_dl
        Primer paso, como en :class:`ArcLengthSolver`: fracción de la carga de
        referencia (por omisión 0.1), con ``Δl₁ = Δλ₁·cᵀ·K₀⁻¹·F_ref``, o
        incremento explícito ``Δl₁`` de la magnitud controlada, en sus
        unidades.
    control_target : float, opcional
        Valor de la magnitud controlada ``cᵀ·U`` en el que se detiene el
        trazado (ADR 0025), en sus unidades. Es un valor **absoluto**: en una
        etapa (ADR 0026) la magnitud parte de su valor en ``u0``, y el
        trazado avanza hacia el objetivo, sea mayor o menor (descarga). Si la
        carga de referencia aleja la magnitud del objetivo, λ se vuelve
        negativa: la carga que crece se invierte. El paso que lo alcanzaría se
        recorta a lo que falta, así que el último punto coincide con el
        objetivo a precisión de máquina. Deja ``termination =
        "control_target"`` y ``reached_control_target = True``. Puede
        combinarse con ``max_lambda``: se detiene en lo que llegue primero
        (con un objetivo, conviene un ``max_lambda`` holgado). Un objetivo
        igual al valor de partida no traza nada y lo avisa.

    dl_max_factor : float, default 1.0
        Como en :class:`ArcLengthSolver`, pero **sin crecimiento por omisión**:
        la magnitud controlada avanza a ritmo constante, como en un ensayo
        controlado por CMOD. Con ablandamiento local el problema incremental
        admite varios equilibrios, y un paso que sobrepasa el pico más que la
        imperfección que localiza el daño puede converger a otro (medido en la
        barra con un elemento debilitado un 5 %: pasos del 10 % de la carga
        dañan todos los elementos; también con control por desplazamiento).
        El resto de parámetros y el tratamiento de ``max_lambda``/``max_steps``
        son los de :class:`ArcLengthSolver`; ``dl`` es aquí el incremento de la
        magnitud controlada por paso.

    La magnitud, el objetivo y el estado del trazado (``control_final``,
    ``reached_control_target``) los define :class:`ArcLengthSolver`; aquí
    ``control`` es obligatorio y además fija la restricción del paso, así
    que el último paso se recorta a lo que falta y no hace falta repetirlo.
    """

    def __init__(self, assembler, convergence: ConvergenceCriterion | None = None,
                 max_iter=20, max_lambda=1.0, initial_dl: float | None = None,
                 max_steps=100, dl_grow_factor=1.5, dl_max_factor=1.0,
                 dl_shrink_factor=0.6, dl_grow_iter_threshold=4,
                 dl_shrink_iter_threshold=8, linear_algebra: str = "auto", *,
                 control, initial_dlambda: float | None = None,
                 control_target: float | None = None):
        super().__init__(
            assembler, convergence=convergence, max_iter=max_iter,
            max_lambda=max_lambda, initial_dl=initial_dl, max_steps=max_steps,
            dl_grow_factor=dl_grow_factor, dl_max_factor=dl_max_factor,
            dl_shrink_factor=dl_shrink_factor,
            dl_grow_iter_threshold=dl_grow_iter_threshold,
            dl_shrink_iter_threshold=dl_shrink_iter_threshold,
            linear_algebra=linear_algebra, initial_dlambda=initial_dlambda,
            control=control, control_target=control_target,
        )

    def _limit_step(self, dl: float, U) -> float:
        # La restricción fija cᵀ·ΔU = Δl: basta recortar a lo que falta.
        if self.control_target is None:
            return dl
        return min(dl, self._remaining(U))

    def _crossed_target(self, U) -> bool:
        # El recorte de ``_limit_step`` llega exactamente: no hay cruce que
        # repetir (el redondeo lo absorbe ``_target_reached``).
        return False

    def _control_measure(self, du: np.ndarray) -> float:
        return float(self._c @ du)

    def _predictor_dlambda(self, du_t: np.ndarray, sign: float, dl: float) -> float:
        # Sin regla de signo: Δλ = Δl/(cᵀ·du_t) sale con el sentido correcto
        # (negativo en un retroceso, donde la tangente hace bajar la carga
        # mientras la magnitud controlada sigue creciendo).
        m = self._control_measure(du_t)
        if abs(m) <= _INSENSITIVE_RTOL * np.linalg.norm(self._c) * np.linalg.norm(du_t):
            raise RuntimeError(
                "IndirectDisplacementSolver: la magnitud controlada no responde a "
                "la carga en el estado actual (cᵀ·K⁻¹·F_ref ≈ 0). Elegir una "
                "combinación de grados de libertad que cambie con la carga.")
        return dl / m

    def _step_problem(self, U_current, lambda_curr, F_ext_ref, free_dofs, dl):
        return _ArcProblem(
            self.assembler, U_current, lambda_curr, F_ext_ref, free_dofs,
            mode="indirect", dl=dl, max_lambda=self.max_lambda,
            control_vector=self._c, stage=self._stage,
        )

    def _first_dl(self, du_t: np.ndarray) -> float:
        m = self._control_measure(du_t)
        if self.control_target is not None and np.isfinite(m) and m < 0.0:
            # Con objetivo, el sentido lo fija el objetivo: si la carga de
            # referencia aleja la magnitud de él, λ disminuye (Δλ = Δl/m < 0).
            _log.info(
                "  La carga de referencia aleja la magnitud controlada del "
                "objetivo: el trazado avanza hacia él con λ decreciente.")
            if self.initial_dl is not None:
                dl = self.initial_dl
            else:
                dl = self.initial_dlambda * abs(m)
                _log.info(f"  Primer paso: |Δλ₁| = {self.initial_dlambda:g} de la "
                          f"carga de referencia ⇒ Δl₁ = {dl:.4e}")
            self.dl_reference = float(dl)
            return float(dl)
        if not (np.isfinite(m) and m > 0.0):
            raise ValueError(
                f"IndirectDisplacementSolver: con la carga de referencia la "
                f"magnitud controlada no crece (cᵀ·K⁻¹·F_ref = {m:.3e}). El control "
                f"exige una magnitud que aumente al cargar: invertir el signo de "
                f"los coeficientes o elegir otros grados de libertad (o dar un "
                f"control_target, que fija el sentido).")
        return super()._first_dl(du_t)
