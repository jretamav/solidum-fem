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
"""
from __future__ import annotations

import math

import numpy as np

from solidum.math.convergence import ConvergenceCriterion
from solidum.math.solvers.arclength import ArcLengthSolver, _ArcProblem
from solidum.math.solvers.corrector import CorrectionAborted
from solidum.registry import SolverRegistry

# |cᵀ·du_t| por debajo de esta fracción de ‖c‖·‖du_t‖: la magnitud
# controlada no responde a la carga y ``Δλ`` no está determinado.
_INSENSITIVE_RTOL = 1.0e-12


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
        libertad deben ser libres (sin apoyo ni restricción multipunto), y la
        combinación debe crecer al aplicar la carga de referencia.
    initial_dlambda, initial_dl
        Primer paso, como en :class:`ArcLengthSolver`: fracción de la carga de
        referencia (por omisión 0.1), con ``Δl₁ = Δλ₁·cᵀ·K₀⁻¹·F_ref``, o
        incremento explícito ``Δl₁`` de la magnitud controlada, en sus
        unidades.

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
                 control, initial_dlambda: float | None = None):
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
        self._c: np.ndarray | None = None

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

    # --- puntos de extensión del padre ---------------------------------------

    def _on_solve_start(self) -> None:
        """Vector ``c`` sobre la numeración de ecuaciones vigente."""
        domain = self.assembler.domain
        ndof = domain.total_dofs
        free = set(int(i) for i in self.assembler.constraint_set.free_dofs(ndof))
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
            if index not in free:
                raise ValueError(
                    f"IndirectDisplacementSolver: el grado de libertad {dof!r} del "
                    f"nodo {node.id} está restringido (apoyo o restricción "
                    f"multipunto); la magnitud controlada sólo puede combinar "
                    f"grados de libertad libres. Para controlar el desplazamiento "
                    f"de un nodo respecto de un apoyo fijo basta con el nodo solo.")
            c[index] += coef
        if not np.any(c):
            raise ValueError(
                "IndirectDisplacementSolver: los términos de 'control' se anulan "
                "entre sí (c = 0).")
        self._c = c

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
            control_vector=self._c,
        )

    def _first_dl(self, du_t: np.ndarray) -> float:
        m = self._control_measure(du_t)
        if not (np.isfinite(m) and m > 0.0):
            raise ValueError(
                f"IndirectDisplacementSolver: con la carga de referencia la "
                f"magnitud controlada no crece (cᵀ·K⁻¹·F_ref = {m:.3e}). El control "
                f"exige una magnitud que aumente al cargar: invertir el signo de "
                f"los coeficientes o elegir otros grados de libertad.")
        return super()._first_dl(du_t)
