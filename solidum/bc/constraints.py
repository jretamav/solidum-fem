"""Restricciones afines sobre los DOFs del sistema (ADR 0004).

Una restricción afín tiene la forma::

    u_s = g_s + Σ α_si · u_mi

donde ``s`` es el DOF esclavo (eliminado del sistema reducido), ``m_i`` los
DOFs maestro (sobreviven), ``α_si`` los coeficientes y ``g_s`` el término
independiente. Para Dirichlet pura no hay maestros y ``g_s`` es el valor
prescrito.

Fase 1 expone solo Dirichlet (``add_dirichlet``); fase 2 añade ``add_linear``
y el cierre transitivo de cadenas master-slave.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
import scipy.sparse as sp


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
    """Conjunto de restricciones afines (ADR 0004).

    Construido típicamente por el ``Assembler`` a partir de
    ``Node.boundary_conditions`` (Dirichlet) y ``Domain.linear_constraints``
    (multipoint constraints). Producto principal: el operador disperso ``T`` y
    el vector ``g_indep`` tales que ``u = T·u_libre + g_indep``.

    El cierre transitivo de cadenas master-slave (``s ← m`` con ``m ← m'``
    ⇒ ``s ← m'``) se aplica de forma diferida en :meth:`build`.
    """

    def __init__(self) -> None:
        self._constraints: dict[int, AffineConstraint] = {}
        self._closed: bool = False
        self._cache: tuple[sp.csr_matrix, np.ndarray, int] | None = None

    # ------------------------------------------------------------------
    # API de declaración
    # ------------------------------------------------------------------

    def add_dirichlet(self, slave_dof: int, value: float = 0.0) -> None:
        """Restricción ``u[slave_dof] = value``.

        Si el DOF ya tenía una restricción Dirichlet con el mismo valor, la
        operación es idempotente. Si difiere, se levanta ``ValueError``.
        """
        self._add(slave_dof, masters=(), coefficients=(), g=float(value))

    def add_linear(
        self,
        slave_dof: int,
        masters: Sequence[int],
        coefficients: Sequence[float],
        g: float = 0.0,
    ) -> None:
        """Restricción afín lineal ``u[slave_dof] = g + Σ α_i · u[m_i]``.

        Útil para apoyos en plano oblicuo, periodicidad, uniones rígidas y
        simetrías no alineadas con los ejes globales.

        Validaciones inmediatas (Reglas.md §1):
        - ``masters`` y ``coefficients`` deben tener la misma longitud.
        - ``slave_dof`` no puede aparecer entre sus propios maestros directos.
        - Si ``slave_dof`` ya estaba declarado y la nueva restricción difiere,
          error.

        Las cadenas ``s ← m`` con ``m`` también esclavo se resuelven en
        :meth:`build` mediante cierre transitivo.
        """
        masters = tuple(int(m) for m in masters)
        coefficients = tuple(float(c) for c in coefficients)
        if len(masters) != len(coefficients):
            raise ValueError(
                "masters y coefficients deben tener la misma longitud "
                f"({len(masters)} vs {len(coefficients)})."
            )
        if slave_dof in masters:
            raise ValueError(
                f"DOF {slave_dof} aparece como maestro directo de sí mismo."
            )
        self._add(slave_dof, masters=masters, coefficients=coefficients, g=float(g))

    def _add(
        self,
        slave_dof: int,
        masters: tuple[int, ...],
        coefficients: tuple[float, ...],
        g: float,
    ) -> None:
        existing = self._constraints.get(slave_dof)
        if existing is not None:
            same = (
                existing.masters == masters
                and existing.coefficients == coefficients
                and existing.g == g
            )
            if not same:
                raise ValueError(
                    f"DOF {slave_dof} ya tiene una restricción incompatible: "
                    f"existente (masters={existing.masters}, coef={existing.coefficients}, g={existing.g}); "
                    f"nueva (masters={masters}, coef={coefficients}, g={g})."
                )
            return
        self._constraints[slave_dof] = AffineConstraint(
            slave_dof=slave_dof, masters=masters, coefficients=coefficients, g=g
        )
        self._closed = False
        self._cache = None

    # ------------------------------------------------------------------
    # Cierre transitivo y construcción de T
    # ------------------------------------------------------------------

    def _close(self) -> None:
        """Resuelve cadenas master-slave hasta que ningún maestro sea esclavo.

        Detecta ciclos (``a ← b``, ``b ← a``) y referencias auto-recursivas
        sobrevenidas durante la sustitución.
        """
        if self._closed:
            return

        max_iters = max(len(self._constraints), 1) + 1
        for _ in range(max_iters):
            slave_set = set(self._constraints.keys())
            changed = False
            new_constraints: dict[int, AffineConstraint] = {}
            for s, c in self._constraints.items():
                new_g = c.g
                expanded_pairs: list[tuple[int, float]] = []
                for m, alpha in zip(c.masters, c.coefficients):
                    if m == s:
                        raise ValueError(
                            f"Restricción cíclica: DOF {s} se referencia a sí "
                            "mismo (directa o indirectamente)."
                        )
                    if m in slave_set:
                        sub = self._constraints[m]
                        new_g += alpha * sub.g
                        for mm, beta in zip(sub.masters, sub.coefficients):
                            expanded_pairs.append((mm, alpha * beta))
                        changed = True
                    else:
                        expanded_pairs.append((m, alpha))
                # Combinar coeficientes duplicados (mismo maestro varias veces).
                combined: dict[int, float] = {}
                for m, coef in expanded_pairs:
                    combined[m] = combined.get(m, 0.0) + coef
                new_masters = tuple(combined.keys())
                new_coefs = tuple(combined.values())
                new_constraints[s] = AffineConstraint(
                    slave_dof=s,
                    masters=new_masters,
                    coefficients=new_coefs,
                    g=new_g,
                )
            self._constraints = new_constraints
            if not changed:
                break
        else:
            raise ValueError(
                "No se pudo aplicar el cierre transitivo en "
                f"{max_iters} iteraciones: posible ciclo en las restricciones."
            )

        # Validación final: ningún maestro debe seguir siendo esclavo.
        slave_set = set(self._constraints.keys())
        for s, c in self._constraints.items():
            for m in c.masters:
                if m in slave_set:
                    raise ValueError(
                        f"Tras cierre transitivo, el DOF {m} sigue figurando "
                        f"como maestro de {s} y como esclavo en otra restricción."
                    )

        self._closed = True

    def build(self, ndof: int) -> tuple[sp.csr_matrix, np.ndarray]:
        """Devuelve ``(T, g_indep)`` tras cierre transitivo.

        ``T`` tiene shape ``(ndof, n_libre)`` y para cada DOF libre ``f``
        cumple ``T[f, k] = 1`` con ``k`` la posición de ``f`` en la lista
        ordenada de DOFs libres. Para cada DOF esclavo ``s`` lleva los
        coeficientes ``α_si`` en las columnas correspondientes a sus
        maestros. ``g_indep`` es no nulo solo en filas esclavas con término
        independiente.
        """
        self._close()
        if self._cache is not None and self._cache[2] == ndof:
            return self._cache[0], self._cache[1]

        slave_set = set(self._constraints.keys())
        free_dofs = [i for i in range(ndof) if i not in slave_set]
        free_idx = {dof: k for k, dof in enumerate(free_dofs)}
        n_free = len(free_dofs)

        rows: list[int] = []
        cols: list[int] = []
        vals: list[float] = []

        # Identidad rectangular para los DOFs libres.
        for f, k in free_idx.items():
            rows.append(f)
            cols.append(k)
            vals.append(1.0)

        # Coeficientes α_si en las filas esclavas.
        g_indep = np.zeros(ndof)
        for s, c in self._constraints.items():
            g_indep[s] = c.g
            for m, alpha in zip(c.masters, c.coefficients):
                if m not in free_idx:
                    raise ValueError(
                        f"DOF maestro {m} no es libre tras cierre — bug interno."
                    )
                rows.append(s)
                cols.append(free_idx[m])
                vals.append(alpha)

        T = sp.csr_matrix(
            (vals, (rows, cols)),
            shape=(ndof, n_free),
        )
        self._cache = (T, g_indep, ndof)
        return T, g_indep

    # ------------------------------------------------------------------
    # Inspección
    # ------------------------------------------------------------------

    @property
    def slave_dofs(self) -> np.ndarray:
        """DOFs prescritos en orden ascendente."""
        return np.array(sorted(self._constraints.keys()), dtype=np.int64)

    @property
    def slave_values(self) -> np.ndarray:
        """Términos independientes ``g_s`` alineados con ``slave_dofs``.

        Para restricciones lineales no Dirichlet ``g_s`` puede ser cero
        aunque la restricción no sea trivial (los coeficientes ``α_si``
        capturan el resto).
        """
        return np.array(
            [self._constraints[s].g for s in sorted(self._constraints)],
            dtype=np.float64,
        )

    def free_dofs(self, ndof: int) -> np.ndarray:
        """DOFs libres (no esclavos), ordenados."""
        slave_set = set(self._constraints.keys())
        return np.array([i for i in range(ndof) if i not in slave_set], dtype=np.int64)

    def __len__(self) -> int:
        return len(self._constraints)
