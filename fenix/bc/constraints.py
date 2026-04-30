"""Restricciones afines sobre los DOFs del sistema (ADR 0004).

Una restricción afín tiene la forma::

    u_s = g_s + Σ α_si · u_mi

donde ``s`` es el DOF esclavo (eliminado del sistema reducido), ``m_i`` los
DOFs maestro (sobreviven), ``α_si`` los coeficientes y ``g_s`` el término
independiente. Para Dirichlet pura no hay maestros y ``g_s`` es el valor
prescrito.

Fase 1 expone solo Dirichlet (``add_dirichlet``). MPC lineales y cierre
transitivo de cadenas master-slave entran en la fase 2.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class AffineConstraint:
    """Restricción afín ``u_s = g + Σ α_i · u_mi``.

    Para Dirichlet pura ``masters`` es vacío y ``g`` es el valor prescrito.
    """

    slave_dof: int
    masters: tuple[int, ...]
    coefficients: tuple[float, ...]
    g: float

    def is_dirichlet(self) -> bool:
        return len(self.masters) == 0


class ConstraintSet:
    """Conjunto de restricciones afines (ADR 0004 fase 1).

    Construido típicamente por el ``Assembler`` a partir de
    ``Node.boundary_conditions``. Producto principal: el par ``(slave_dofs,
    slave_values)`` que permite construir la matriz ``T`` y el vector ``g``
    del modelo de eliminación directa.
    """

    def __init__(self) -> None:
        self._constraints: dict[int, AffineConstraint] = {}

    def add_dirichlet(self, slave_dof: int, value: float = 0.0) -> None:
        """Restricción ``u[slave_dof] = value``.

        Si el DOF ya tenía una restricción Dirichlet con el mismo valor, la
        operación es idempotente. Si el valor difiere, se levanta
        ``ValueError`` (validación temprana, Reglas.md §1).
        """
        existing = self._constraints.get(slave_dof)
        if existing is not None:
            if existing.masters or existing.g != value:
                raise ValueError(
                    f"DOF {slave_dof} ya tiene una restricción incompatible: "
                    f"existente g={existing.g}, masters={existing.masters}; "
                    f"nueva g={value}."
                )
            return
        self._constraints[slave_dof] = AffineConstraint(
            slave_dof=slave_dof, masters=(), coefficients=(), g=float(value)
        )

    @property
    def slave_dofs(self) -> np.ndarray:
        """DOFs prescritos en orden ascendente."""
        return np.array(sorted(self._constraints.keys()), dtype=np.int64)

    @property
    def slave_values(self) -> np.ndarray:
        """Términos independientes ``g_s`` alineados con ``slave_dofs``."""
        return np.array(
            [self._constraints[s].g for s in sorted(self._constraints)],
            dtype=np.float64,
        )

    def free_dofs(self, ndof: int) -> np.ndarray:
        """DOFs libres (no esclavos), ordenados."""
        bc = set(self._constraints.keys())
        return np.array([i for i in range(ndof) if i not in bc], dtype=np.int64)

    def __len__(self) -> int:
        return len(self._constraints)
