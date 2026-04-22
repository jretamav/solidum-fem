"""API pública de resultados de una solución (ADR 0002).

Contiene los tipos de dato que exponen los resultados al exterior:

- :class:`ElementForces`: esfuerzos internos en ejes locales de un elemento, en
  los nodos i, j. Convenciones de signo definidas en ``Reglas.md §5``.
- :class:`SolveResult`: agregado inmutable producido al final de ``solver.solve``
  y retenido por el ``Domain`` en ``last_result``.

Ambos son ``frozen`` para impedir que un consumidor mute resultados y sorprenda
a otro posterior.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Literal

import numpy as np

if TYPE_CHECKING:
    from fenix.core.domain import Domain
    from fenix.math.assembly import Assembler


ElementKind = Literal["truss", "cable", "frame2d", "frame3d"]

# Claves válidas por familia. Un componente ausente del dict significa "no aplica"
# a ese tipo de elemento (p. ej. no hay M en un truss) — no se rellena con ceros.
_VALID_COMPONENTS: dict[ElementKind, frozenset[str]] = {
    "truss":   frozenset({"N"}),
    "cable":   frozenset({"N"}),
    "frame2d": frozenset({"N", "V", "M"}),
    "frame3d": frozenset({"N", "Vy", "Vz", "T", "My", "Mz"}),
}


@dataclass(frozen=True)
class ElementForces:
    """Esfuerzos internos en ejes locales, evaluados en los nodos i, j.

    Parameters
    ----------
    kind
        Familia del elemento. Determina qué claves pueden aparecer en ``components``.
    components
        Diccionario ``{nombre: array(2,)}`` con el valor en el nodo i (índice 0) y
        en el nodo j (índice 1). Las claves válidas dependen de ``kind``:

        - ``"truss"`` / ``"cable"``: ``{"N"}``.
        - ``"frame2d"``: ``{"N", "V", "M"}``. Convención de viga estructural
          (Reglas.md §5): ``V`` rota el diferencial en horario si es positivo,
          ``M`` positivo es sagging.
        - ``"frame3d"``: ``{"N", "Vy", "Vz", "T", "My", "Mz"}``. Convención
          stress-resultant / RHR pura (Reglas.md §5).
    """

    kind: ElementKind
    components: dict[str, np.ndarray]

    def __post_init__(self) -> None:
        valid = _VALID_COMPONENTS.get(self.kind)
        if valid is None:
            raise ValueError(f"ElementForces.kind desconocido: {self.kind!r}")

        extra = set(self.components) - valid
        if extra:
            raise ValueError(
                f"Componentes inválidos para kind={self.kind!r}: {sorted(extra)}. "
                f"Válidos: {sorted(valid)}"
            )

        for name, arr in self.components.items():
            if not isinstance(arr, np.ndarray) or arr.shape != (2,):
                raise ValueError(
                    f"components[{name!r}] debe ser np.ndarray shape (2,), "
                    f"recibido {type(arr).__name__} shape={getattr(arr, 'shape', None)}"
                )

    def at_node_i(self) -> dict[str, float]:
        """Valores en el nodo i como escalares."""
        return {k: float(v[0]) for k, v in self.components.items()}

    def at_node_j(self) -> dict[str, float]:
        """Valores en el nodo j como escalares."""
        return {k: float(v[1]) for k, v in self.components.items()}


@dataclass(frozen=True)
class SolveResult:
    """Resultado agregado de una solución. Inmutable; calculado eager al final
    de ``solver.solve``.

    Parameters
    ----------
    U
        Desplazamientos globales, shape ``(n_dof,)``.
    F_applied
        Cargas aplicadas tal como las ve el sistema: nodales del usuario más
        equivalentes de elemento (distribuidas, térmicas) ya ensambladas a
        nodos. Shape ``(n_dof,)``. No incluye reacciones.
    R
        Reacciones globales, shape ``(n_dof,)``. Ceros en DOF libres; valor no
        nulo solo en DOF restringidos. Se expone redundantemente con
        ``reactions_by_node`` por conveniencia (vector para GUIs, dict para scripts).
    reactions_by_node
        Vista ``{node_id: {dof_name: valor}}`` filtrada a nodos con al menos un
        DOF restringido.
    element_forces
        ``{elem_id: ElementForces}`` con los esfuerzos internos de cada elemento
        que implementa ``internal_forces``. Elementos sólidos pueden omitirse.
    converged
        ``True`` si el solver alcanzó el criterio de convergencia. Para solvers
        lineales, siempre ``True`` tras una resolución exitosa.
    num_steps
        Número de pasos/iteraciones realizados. ``1`` para solver lineal; ``N``
        para Newton-Raphson / arc-length.
    """

    U: np.ndarray
    F_applied: np.ndarray
    R: np.ndarray
    reactions_by_node: dict[int, dict[str, float]] = field(default_factory=dict)
    element_forces: dict[int, ElementForces] = field(default_factory=dict)
    converged: bool = True
    num_steps: int = 1


def build_solve_result(
    domain: "Domain",
    assembler: "Assembler",
    U: np.ndarray,
    F_applied: np.ndarray,
    *,
    converged: bool = True,
    num_steps: int = 1,
) -> SolveResult:
    """Post-procesa una solución ``U`` para construir un ``SolveResult`` completo.

    Calcula:
    - Reacciones ``R`` globales: ``F_int(U) − F_applied`` en DOFs restringidos,
      cero en DOFs libres (por equilibrio al convergir).
    - ``reactions_by_node``: vista filtrada a nodos con al menos un DOF con
      condición de Dirichlet.
    - ``element_forces``: ``element.internal_forces(U)`` para cada elemento
      que implemente el contrato (los que devuelvan ``None`` se omiten).

    El llamador (típicamente ``fenix.run``) asigna el resultado a
    ``domain.last_result``.
    """
    ndof = U.shape[0]

    _, F_int = assembler.assemble_non_linear_system(U)

    if not getattr(assembler, "_bc_built", False):
        assembler._build_bc_arrays()
    bc_dofs = assembler._bc_dofs

    R = np.zeros(ndof)
    if bc_dofs.size > 0:
        R[bc_dofs] = F_int[bc_dofs] - F_applied[bc_dofs]

    reactions_by_node: dict[int, dict[str, float]] = {}
    for node_id, node in domain.nodes.items():
        if not node.boundary_conditions:
            continue
        per_dof = {
            dof_name: float(R[node.dofs[dof_name]])
            for dof_name in node.boundary_conditions
            if dof_name in node.dofs
        }
        if per_dof:
            reactions_by_node[node_id] = per_dof

    element_forces: dict[int, ElementForces] = {}
    for elem_id, elem in domain.elements.items():
        ef = elem.internal_forces(U)
        if ef is not None:
            element_forces[elem_id] = ef

    return SolveResult(
        U=U,
        F_applied=F_applied,
        R=R,
        reactions_by_node=reactions_by_node,
        element_forces=element_forces,
        converged=converged,
        num_steps=num_steps,
    )
