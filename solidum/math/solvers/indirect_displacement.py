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

import math

import numpy as np

from solidum.math.convergence import ConvergenceCriterion
from solidum.math.solvers._shared import _log
from solidum.math.solvers.arclength import ArcLengthSolver, _ArcProblem
from solidum.math.solvers.corrector import CorrectionAborted
from solidum.math.solvers.diagnostics import Termination
from solidum.registry import SolverRegistry

# |cᵀ·du_t| por debajo de esta fracción de ‖c‖·‖du_t‖: la magnitud
# controlada no responde a la carga y ``Δλ`` no está determinado.
_INSENSITIVE_RTOL = 1.0e-12
# ``control_target`` se da por alcanzado a esta fracción de su valor: la
# llegada es exacta salvo el redondeo de acumular ``cᵀ·U`` paso a paso.
_TARGET_RTOL = 1.0e-12


class _IndirectProblem(_ArcProblem):
    """Paso con la restricción lineal ``cᵀ·ΔU = Δl`` (``mode="indirect"``)."""

    def __init__(self, *args, control_vector: np.ndarray, **kwargs):
        super().__init__(*args, **kwargs)
        self.c = control_vector

    def constraint(self, x, du_R, du_t):
        if self.mode != "indirect":
            return super().constraint(x, du_R, du_t)
        _, _, dU_iter = x
        denom = float(self.c @ du_t)
        if abs(denom) <= _INSENSITIVE_RTOL * np.linalg.norm(self.c) * np.linalg.norm(du_t):
            raise CorrectionAborted(
                "cᵀ·du_t ≈ 0: la magnitud controlada no responde a la carga en "
                "este estado; el paso no determina Δλ.")
        ddlambda = (self.dl - float(self.c @ (dU_iter + du_R))) / denom
        return du_R + ddlambda * du_t, ddlambda


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
        )
        self.control = self._normalize_control(control)
        # ``_c`` es el vector de la restricción, orientado hacia el objetivo
        # (``_direction`` = ±1); ``_c_raw``, el de la magnitud tal como se
        # declaró, con el que se reporta.
        self._c: np.ndarray | None = None
        self._c_raw: np.ndarray | None = None
        self._direction: float = 1.0
        self._target_scale: float = 0.0
        if control_target is not None:
            control_target = float(control_target)
            if not math.isfinite(control_target):
                raise ValueError(
                    f"IndirectDisplacementSolver: control_target={control_target} debe "
                    f"ser un número finito.")
        self.control_target: float | None = control_target
        # Estado propio del trazado (ADR 0025), además del de ArcLengthSolver.
        self.reached_control_target: bool = False
        self.control_final: float = 0.0

    @staticmethod
    def _normalize_control(control) -> list[tuple]:
        """Lista de ``(nodo, dof, coef)`` validada en forma (no en existencia:
        los nodos y sus grados de libertad se comprueban al resolver)."""
        name = "IndirectDisplacementSolver"
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

    # --- estado del trazado y parada en el objetivo (ADR 0025) ----------------

    def _begin_trace(self) -> None:
        super()._begin_trace()
        self.reached_control_target = False
        self.control_final = 0.0

    def _control_value(self, U) -> float:
        """Magnitud controlada tal como se declaró, ``cᵀ·U``."""
        return float(self._c_raw @ U)

    def _remaining(self, U) -> float:
        """Lo que falta hasta el objetivo, medido en el sentido del trazado."""
        return self._direction * (self.control_target - self._control_value(U))

    def _target_reached(self, U) -> bool:
        if self.control_target is None or U is None or self._c is None:
            return False
        return self._remaining(U) <= _TARGET_RTOL * self._target_scale

    def _limit_step(self, dl: float, U) -> float:
        if self.control_target is None:
            return dl
        return min(dl, self._remaining(U))

    def _step_summary(self, U) -> str:
        return f" | control = {self._control_value(U):.6e}"

    def _record(self, lambda_curr, steps, U) -> None:
        super()._record(lambda_curr, steps, U)
        if U is not None and self._c_raw is not None:
            self.control_final = self._control_value(U)

    def _finish(self, lambda_curr, step, U=None) -> None:
        super()._finish(lambda_curr, step, U)
        self.reached_control_target = self.termination is Termination.CONTROL_TARGET

    # --- puntos de extensión del padre ---------------------------------------

    def _on_solve_start(self) -> None:
        """Vector ``c`` sobre la numeración de ecuaciones vigente.

        ``c`` actúa sobre el incremento completo ``ΔU = T·Δu_libre + Δλ·g``
        (ADR 0004), que incluye los esclavos reconstruidos por ``T``: el
        esclavo de una restricción multipunto con maestros (fila de ``T`` no
        vacía) se puede controlar, porque su incremento lo fijan sus maestros.
        Un grado de libertad prescrito (fila de ``T`` vacía: apoyo declarado, o
        esclavo cuyos maestros están todos prescritos, que hereda el apoyo por
        el cierre transitivo) se rechaza: su incremento es ``Δλ·g_s``, no
        depende del equilibrio. También se rechaza una combinación cuya
        proyección ``Tᵀ·c`` sobre los grados de libertad libres se anula (p.
        ej. ``u_s − u_m`` con ``u_s = u_m + g``)."""
        domain = self.assembler.domain
        ndof = domain.total_dofs
        T = self.assembler.constraint_set.build(ndof)[0].tocsr()
        c = np.zeros(ndof)
        for node, dof, coef in self.control:
            if not hasattr(node, 'dofs'):
                node_obj = domain.get_node(node)
                if node_obj is None:
                    raise ValueError(
                        f"IndirectDisplacementSolver: el nodo {node!r} de 'control' "
                        f"no existe en el modelo.")
                node = node_obj
            if dof not in node.dofs:
                raise ValueError(
                    f"IndirectDisplacementSolver: el nodo {node.id} no tiene el grado "
                    f"de libertad {dof!r}; tiene {sorted(node.dofs)}.")
            index = node.dofs[dof]
            row = T.data[T.indptr[index]:T.indptr[index + 1]]
            if not np.any(row):
                raise ValueError(
                    f"IndirectDisplacementSolver: el grado de libertad {dof!r} del "
                    f"nodo {node.id} está prescrito (apoyo, o esclavo de una "
                    f"restricción multipunto cuyos maestros están todos prescritos): "
                    f"su incremento no depende del equilibrio y no puede formar parte "
                    f"de la magnitud controlada. Para controlar el desplazamiento de "
                    f"un nodo respecto de un apoyo fijo basta con el nodo solo.")
            c[index] += coef
        # Proyección sobre los grados de libertad libres: cᵀ·ΔU = (Tᵀc)ᵀ·Δu_libre
        # + Δλ·cᵀg. Sin restricciones multipunto, Tᵀc = 0 sólo si c = 0.
        c_free = T.T @ c
        scale = abs(T).T @ np.abs(c)
        if np.linalg.norm(c_free) <= _INSENSITIVE_RTOL * np.linalg.norm(scale):
            raise ValueError(
                "IndirectDisplacementSolver: los términos de 'control' se anulan "
                "entre sí (Tᵀ·c = 0: tras aplicar las restricciones multipunto, la "
                "combinación no depende de ningún grado de libertad libre).")
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
                    f"IndirectDisplacementSolver: control_target = "
                    f"{self.control_target:.6e} es el valor de partida de la "
                    f"magnitud controlada; no hay nada que trazar.")
            self._direction = 1.0 if gap >= 0.0 else -1.0
            self._c = self._direction * c

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
        return _IndirectProblem(
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
