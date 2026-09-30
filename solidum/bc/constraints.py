"""Restricciones afines sobre los DOFs del sistema (ADR 0004).

Una restricción afín tiene la forma::

    u_s = g_s + Σ α_si · u_mi

donde ``s`` es el DOF esclavo (eliminado del sistema reducido), ``m_i`` los
DOFs maestro (sobreviven), ``α_si`` los coeficientes y ``g_s`` el término
independiente. Para Dirichlet pura no hay maestros y ``g_s`` es el valor
prescrito.

Fase 1 expone solo Dirichlet (``add_dirichlet``); fase 2 añade ``add_linear``
y el cierre transitivo de cadenas master-slave.

Dos inspecciones independientes del cierre sirven a los solvers que necesitan
distinguir el Dirichlet declarado de las MPC: :attr:`ConstraintSet.dirichlet_dofs`
(los DOF con valor prescrito tal como se declararon) y
:meth:`ConstraintSet.build_g` (el ``g`` completo para unos valores de Dirichlet
dados, con la composición del cierre ``g_s = g_mpc,s + Σ α·ḡ_d``).
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
        # Restricciones tal como se declararon; el cierre no las modifica.
        self._declared: dict[int, AffineConstraint] = {}
        # Restricciones vigentes: las declaradas hasta el cierre, las cerradas
        # después (mismas claves en ambos diccionarios).
        self._constraints: dict[int, AffineConstraint] = {}
        # Descomposición afín del g cerrado de cada esclavo en la parte propia
        # de las MPC y la heredada de los Dirichlet declarados:
        # g_s = const_s + Σ_d w_sd · ḡ_d. La rellena ``_close``.
        self._g_affine: dict[int, tuple[float, dict[int, float]]] = {}
        self._closed: bool = False
        self._cache: tuple[sp.csr_matrix, np.ndarray, int] | None = None
        self._g_cache: tuple[sp.csr_matrix, np.ndarray, int] | None = None

    # ------------------------------------------------------------------
    # API de declaración
    # ------------------------------------------------------------------

    def add_dirichlet(self, slave_dof: int, value: float = 0.0) -> None:
        """Restricción ``u[slave_dof] = value``.

        Si el DOF ya tenía una restricción Dirichlet con el mismo valor, la
        operación es idempotente. Si difiere, se levanta ``ValueError``. El DOF
        queda en :attr:`dirichlet_dofs`.
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

        Con ``masters`` vacío la restricción es un Dirichlet ``u[s] = g`` y
        cuenta como tal en :attr:`dirichlet_dofs`.
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
        # Se compara con lo declarado, no con lo cerrado: redeclarar la misma
        # restricción tras un ``build`` sigue siendo idempotente.
        existing = self._declared.get(slave_dof)
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
        constraint = AffineConstraint(
            slave_dof=slave_dof, masters=masters, coefficients=coefficients, g=g
        )
        self._declared[slave_dof] = constraint
        self._constraints[slave_dof] = constraint
        self._closed = False
        self._cache = None
        self._g_cache = None

    # ------------------------------------------------------------------
    # Cierre transitivo y construcción de T
    # ------------------------------------------------------------------

    def _close(self) -> None:
        """Resuelve cadenas master-slave hasta que ningún maestro sea esclavo.

        Detecta ciclos (``a ← b``, ``b ← a``) y referencias auto-recursivas
        sobrevenidas durante la sustitución. Parte siempre de las restricciones
        declaradas. En paralelo con ``g`` lleva su descomposición afín
        (``_g_affine``): la parte propia de las MPC y los pesos con que cada
        Dirichlet declarado entra en el ``g`` de cada esclavo, que usa
        :meth:`build_g`.
        """
        if self._closed:
            return

        constraints = dict(self._declared)
        g_affine: dict[int, tuple[float, dict[int, float]]] = {
            s: ((0.0, {s: 1.0}) if c.is_dirichlet() else (c.g, {}))
            for s, c in constraints.items()
        }
        max_iters = max(len(constraints), 1) + 1
        for _ in range(max_iters):
            slave_set = set(constraints.keys())
            changed = False
            new_constraints: dict[int, AffineConstraint] = {}
            new_affine: dict[int, tuple[float, dict[int, float]]] = {}
            for s, c in constraints.items():
                new_g = c.g
                new_const, weights = g_affine[s]
                new_weights = dict(weights)
                expanded_pairs: list[tuple[int, float]] = []
                for m, alpha in zip(c.masters, c.coefficients):
                    if m == s:
                        raise ValueError(
                            f"Restricción cíclica: DOF {s} se referencia a sí "
                            "mismo (directa o indirectamente)."
                        )
                    if m in slave_set:
                        sub = constraints[m]
                        new_g += alpha * sub.g
                        sub_const, sub_weights = g_affine[m]
                        new_const += alpha * sub_const
                        for d, w in sub_weights.items():
                            new_weights[d] = new_weights.get(d, 0.0) + alpha * w
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
                new_affine[s] = (new_const, new_weights)
            constraints = new_constraints
            g_affine = new_affine
            if not changed:
                break
        else:
            raise ValueError(
                "No se pudo aplicar el cierre transitivo en "
                f"{max_iters} iteraciones: posible ciclo en las restricciones."
            )

        # Validación final: ningún maestro debe seguir siendo esclavo.
        slave_set = set(constraints.keys())
        for s, c in constraints.items():
            for m in c.masters:
                if m in slave_set:
                    raise ValueError(
                        f"Tras cierre transitivo, el DOF {m} sigue figurando "
                        f"como maestro de {s} y como esclavo en otra restricción."
                    )

        self._constraints = constraints
        self._g_affine = g_affine
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

    def build_g(self, ndof: int, dirichlet_values: np.ndarray) -> np.ndarray:
        """Vector ``g`` (``ndof``) del modelo afín para unos valores de Dirichlet dados.

        Es el ``g_indep`` de :meth:`build` con los valores declarados de los
        DOF de :attr:`dirichlet_dofs` sustituidos por
        ``dirichlet_values[dirichlet_dofs]``; el resto de entradas de
        ``dirichlet_values`` no se lee. Aplica la composición del cierre
        transitivo: un DOF Dirichlet recibe su valor, y un esclavo de MPC
        recibe ``g_s = g_mpc,s + Σ_d α_sd·ḡ_d``, con ``g_mpc,s`` la parte
        propia de las MPC de su cadena y ``α_sd`` el peso con que el Dirichlet
        ``d`` le llega a través de ella (un esclavo atado a un maestro con
        Dirichlet hereda su valor). Los DOF libres reciben ``0``.

        Sirve a los solvers con Dirichlet variable en el tiempo
        (``ThetaMethodSolver``): basta el vector de valores prescritos y las
        filas esclavas se reconstruyen con su restricción. Con los valores
        declarados reproduce ``build(ndof)[1]`` salvo redondeo en las filas
        esclavas de cadenas; sin MPC, exactamente.

        Parameters
        ----------
        ndof : int
            Número total de DOF.
        dirichlet_values : array-like, shape ``(ndof,)``
            Valores prescritos, en la posición global de cada DOF.
        """
        values = np.asarray(dirichlet_values, dtype=np.float64).ravel()
        if values.size != ndof:
            raise ValueError(
                f"build_g: dirichlet_values tiene tamaño {values.size}; se "
                f"esperaban {ndof} (un valor por DOF)."
            )
        self._close()
        if self._g_cache is None or self._g_cache[2] != ndof:
            rows: list[int] = []
            cols: list[int] = []
            vals: list[float] = []
            g_mpc = np.zeros(ndof)
            for s, (const, weights) in self._g_affine.items():
                g_mpc[s] = const
                for d, w in weights.items():
                    rows.append(s)
                    cols.append(d)
                    vals.append(w)
            G = sp.csr_matrix((vals, (rows, cols)), shape=(ndof, ndof))
            self._g_cache = (G, g_mpc, ndof)
        G, g_mpc = self._g_cache[0], self._g_cache[1]
        return g_mpc + G @ values

    # ------------------------------------------------------------------
    # Inspección
    # ------------------------------------------------------------------

    @property
    def dirichlet_dofs(self) -> np.ndarray:
        """DOFs con Dirichlet declarado, en orden ascendente.

        Son los declarados con :meth:`add_dirichlet` o con :meth:`add_linear`
        sin maestros (``u_s = g``, que es un Dirichlet). No incluye los
        esclavos de MPC, tampoco el esclavo que el cierre transitivo deja sin
        maestros por estar atado a un maestro con Dirichlet: ese DOF queda
        determinado, pero su valor lo fija la MPC, no un valor prescrito. El
        resultado no depende de que se haya llamado a :meth:`build`.
        """
        return np.array(
            sorted(s for s, c in self._declared.items() if c.is_dirichlet()),
            dtype=np.int64,
        )

    @property
    def slave_dofs(self) -> np.ndarray:
        """DOFs esclavos en orden ascendente: los Dirichlet y los esclavos de las
        MPC, es decir, todos los que salen del sistema reducido. Para sólo el
        Dirichlet declarado, :attr:`dirichlet_dofs`."""
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
