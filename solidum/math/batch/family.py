"""Familias de lote: agrupación de elementos derivada de los contratos (ADR 0014).

Una *familia* es el conjunto de elementos del dominio que comparten
clase de elemento, **instancia** de material, regla de cuadratura y
número de nodos. La clave se deriva de atributos que ya existen; un
elemento nuevo cae en su propia familia sin tocar el ensamblador.

Un elemento participa en el camino por lotes si su clase declara
``BATCH_KINEMATICS`` y su material declara ``STATE_SCHEMA`` y un kernel
puntual (``batch_kernel()``). Si falta cualquiera de las dos cosas, el
elemento queda en el camino por elemento de siempre, sin penalización
de diseño: el kernel por lotes es una capacidad que cada componente
declara, no un requisito del contrato.

Por qué la instancia de material y no su clase: los parámetros del
material (``E``, ``σ_y``, …) entran al kernel como constantes de la
familia, así que dos instancias distintas de la misma clase son dos
familias.
"""
from __future__ import annotations

import numpy as np

from solidum.constants import BATCH_MEMORY_BUDGET_BYTES, BATCH_PARALLEL_DEFAULT
from solidum.logging import get_logger
from solidum.math.batch.kernels import (
    FLAG_BAD_JACOBIAN,
    solid_family_gauss_kernel,
    solid_family_gauss_kernel_parallel,
    solid_family_kernel,
    solid_family_kernel_parallel,
)
from solidum.math.batch.schema import StateSchema
from solidum.math.batch.state import BatchedElementState, FamilyState

_log = get_logger("assembly.batch")


def material_is_batchable(material) -> bool:
    """El material declara esquema de estado y kernel puntual."""
    if material is None:
        return False
    if getattr(material, "STATE_SCHEMA", None) is None:
        return False
    kernel = getattr(material, "batch_kernel", None)
    return callable(kernel) and kernel() is not None


def element_is_batchable(elem) -> bool:
    """El elemento declara cinemática compilada y su material es batchable."""
    if getattr(type(elem), "BATCH_KINEMATICS", None) is None:
        return False
    if getattr(elem, "points", None) is None:
        return False
    return material_is_batchable(getattr(elem, "material", None))


def _quadrature_signature(elem) -> tuple:
    key = getattr(elem, "quadrature_key", None)
    if key is not None:
        return ("key", key)
    pts = tuple(tuple(float(c) for c in p) for p in elem.points)
    wts = tuple(float(w) for w in elem.weights)
    return ("tuple", pts, wts)


def family_key(elem) -> tuple:
    return (type(elem), id(elem.material), _quadrature_signature(elem), len(elem.nodes))


class Family:
    """Elementos que comparten receta y sus datos por filas."""

    def __init__(self, key: tuple, elements: list, positions: list[int],
                 dof_indices: list[list[int]],
                 memory_budget: int | None = None,
                 parallel: bool | None = None):
        self.key = key
        self.elements = elements
        self.positions = np.asarray(positions, dtype=np.int64)
        self.parallel = BATCH_PARALLEL_DEFAULT if parallel is None else bool(parallel)
        self.kernel = solid_family_kernel_parallel if self.parallel else solid_family_kernel
        self.gauss_kernel = (solid_family_gauss_kernel_parallel if self.parallel
                             else solid_family_gauss_kernel)
        elem0 = elements[0]
        # Familia térmica: el "esfuerzo" del kernel es k·∇T = −q.
        self.is_thermal = hasattr(elem0, "FLUX_DIM")
        self.element_class = type(elem0)
        self.material = elem0.material
        N = len(elements)

        pts, w = elem0.batch_quadrature()
        self.gp_points = np.ascontiguousarray(pts, dtype=np.float64)
        self.gp_weights = np.ascontiguousarray(w, dtype=np.float64)
        self.n_gp = int(self.gp_points.shape[0])
        self.spatial_dim = int(self.gp_points.shape[1])
        self.n_nodes = len(elem0.nodes)
        self.n_dof = len(elem0.DOF_NAMES) * self.n_nodes
        strain_dim = getattr(elem0, "STRAIN_DIM", None)
        self.n_sigma = int(strain_dim if strain_dim is not None else elem0.FLUX_DIM)

        self.kin_fn = self.element_class.BATCH_KINEMATICS
        self.mat_fn = self.material.batch_kernel()
        self.mat_params = np.ascontiguousarray(
            np.asarray(self.material.batch_params(), dtype=np.float64).reshape(-1))
        C = np.asarray(self.material.batch_matrix(), dtype=np.float64)
        if C.shape != (self.n_sigma, self.n_sigma):
            raise ValueError(
                f"{type(self.material).__name__}.batch_matrix() devuelve forma "
                f"{C.shape}; se esperaba ({self.n_sigma}, {self.n_sigma}) para "
                f"{self.element_class.__name__}."
            )
        self.C = np.ascontiguousarray(C)

        self.X = np.ascontiguousarray(np.stack(
            [e.batch_reference_coordinates(self.spatial_dim) for e in elements]
        ), dtype=np.float64)
        scale_attr = getattr(self.element_class, "BATCH_SCALE", None)
        if scale_attr is None:
            self.scale = np.ones(N)
        else:
            self.scale = np.ascontiguousarray(
                [float(getattr(e, scale_attr)) for e in elements], dtype=np.float64)

        self.dof_indices = np.ascontiguousarray(
            np.asarray(dof_indices, dtype=np.int64).reshape(N, self.n_dof))
        self.dof_flat = self.dof_indices.reshape(-1)

        self.schema = StateSchema(self.material.STATE_SCHEMA)
        initial = np.zeros(self.schema.n_state)
        self.schema.pack(self.material.initial_state(), initial)
        self.state = FamilyState(N, self.n_gp, self.schema, self.n_sigma, initial)

        # Posición del bloque contiguo de esta familia en el vector `data`
        # de la COO cacheada; la fija el Assembler al construir la topología.
        # El kernel escribe las K_e directamente en ese bloque (vista
        # ``(n, n_dof, n_dof)`` sobre ``data``): no hay copia ni buffer
        # aparte para las matrices elementales.
        self.ptr_start: int | None = None
        self.n_entries = N * self.n_dof * self.n_dof

        # El único temporal por trozo es F_e (8·n_dof bytes por elemento).
        budget = BATCH_MEMORY_BUDGET_BYTES if memory_budget is None else int(memory_budget)
        per_elem = 8 * self.n_dof
        self.chunk = int(max(1, min(N, budget // per_elem)))
        self._F_out = np.zeros((self.chunk, self.n_dof))

        # Funciones de forma en los puntos de Gauss del elemento de
        # referencia (iguales para toda la familia); se evalúan al primer
        # post-proceso que las pida.
        self._N_gp = None
        self._N_gp_ready = False

        self._adopted = False

    # ------------------------------------------------------------------

    @property
    def n_elements(self) -> int:
        return len(self.elements)

    def adopt_states(self) -> None:
        """Sustituye el ``ElementState`` de cada elemento por una vista sobre
        los arreglos de la familia, conservando el contenido previo, y deja
        en cada elemento la referencia ``_batch_family`` que usa el
        post-proceso por lotes para encontrar la familia desde el dominio."""
        if self._adopted:
            return
        for e, elem in enumerate(self.elements):
            elem._batch_family = self
            plain = getattr(elem, "state", None)
            if plain is None:
                # Elementos sin estado (térmicos): no se les impone uno.
                continue
            view = BatchedElementState(self.state, e)
            view.adopt(plain)
            elem.state = view
        self._adopted = True

    def release(self) -> None:
        """Devuelve a cada elemento un ``ElementState`` clásico con el
        contenido actual de los arreglos (al invalidar la topología) y
        retira la referencia a la familia."""
        if not self._adopted:
            return
        for elem in self.elements:
            if getattr(elem, "_batch_family", None) is self:
                elem._batch_family = None
            st = getattr(elem, "state", None)
            if isinstance(st, BatchedElementState) and st.family_state is self.state:
                elem.state = st.to_plain()
        self._adopted = False

    @property
    def adopted(self) -> bool:
        return self._adopted

    def commit(self) -> None:
        self.state.commit()

    # ------------------------------------------------------------------

    def evaluate(self, U: np.ndarray, data: np.ndarray, F_int: np.ndarray,
                 keep_trial: bool = False) -> None:
        """Evalúa ``K_e`` y ``F_int_e`` de todos los elementos en ``U`` y los
        vuelca en el vector COO ``data`` (contiguo, ``float64``) y en
        ``F_int``. Escribe el estado y los esfuerzos **trial** de la familia,
        salvo con ``keep_trial=True`` (evaluación auxiliar, p. ej. la rigidez
        en ``u = 0`` de ``assemble_system``): entonces el trial previo se
        restaura al salir, como hace ``Element.compute_global_stiffness``."""
        if keep_trial:
            st = self.state
            saved = (st.S_trial.copy(), st.sig_trial.copy())
            try:
                self.evaluate(U, data, F_int, keep_trial=False)
            finally:
                np.copyto(st.S_trial, saved[0])
                np.copyto(st.sig_trial, saved[1])
            return
        if self.ptr_start is None:
            raise RuntimeError("Family.evaluate: la familia no tiene posición en la COO.")
        N = self.n_elements
        n_gp = self.n_gp
        n_dof = self.n_dof
        u = np.ascontiguousarray(U[self.dof_indices], dtype=np.float64)
        st = self.state
        st.flags.fill(0)
        n2 = n_dof * n_dof
        for e0 in range(0, N, self.chunk):
            e1 = min(N, e0 + self.chunk)
            n = e1 - e0
            r0, r1 = e0 * n_gp, e1 * n_gp
            p0 = self.ptr_start + e0 * n2
            K_out = data[p0:p0 + n * n2].reshape(n, n_dof, n_dof)
            F_out = self._F_out[:n]
            self.kernel(
                self.kin_fn, self.mat_fn,
                self.X[e0:e1], u[e0:e1], self.gp_points, self.gp_weights,
                self.scale[e0:e1],
                st.S_committed[r0:r1], st.S_trial[r0:r1],
                self.mat_params, self.C,
                K_out, F_out, st.sig_trial[r0:r1], st.flags[r0:r1],
            )
            self._raise_if_bad_jacobian(st.flags[r0:r1], e0)
            F_int += np.bincount(
                self.dof_flat[e0 * n_dof:e1 * n_dof],
                weights=F_out.reshape(-1), minlength=F_int.shape[0],
            )
        self._report_material_flags(st.flags)

    # ------------------------------------------------------------------

    def _raise_if_bad_jacobian(self, flags: np.ndarray, e_offset: int = 0) -> None:
        """Convierte las marcas ``FLAG_BAD_JACOBIAN`` del kernel en el
        ``ValueError`` del camino por elemento, con los ids afectados."""
        bad = np.flatnonzero(flags == FLAG_BAD_JACOBIAN)
        if bad.size == 0:
            return
        elems = sorted({int(e_offset + r // self.n_gp) for r in bad})
        ids = [self.elements[e].id for e in elems[:10]]
        more = "" if len(elems) <= 10 else f" (y {len(elems) - 10} más)"
        raise ValueError(
            f"Jacobiano negativo o cero en {self.element_class.__name__} "
            f"id={ids if len(ids) > 1 else ids[0]}{more}. Revisa la conectividad "
            f"(orden de nodos) o la distorsión del elemento."
        )

    def _report_material_flags(self, flags: np.ndarray) -> None:
        n_flagged = int(np.count_nonzero(flags > 0))
        if n_flagged:
            self.material.batch_report(n_flagged)

    # ------------------------------------------------------------------
    # Post-proceso por familia (ADR 0014 §9)
    # ------------------------------------------------------------------

    def gauss_shape_functions(self):
        """``N`` en los puntos de Gauss del elemento de referencia,
        ``(n_gp, n_nodos)``, o ``None`` si el elemento no las declara."""
        if not self._N_gp_ready:
            N = self.elements[0].batch_shape_functions(self.gp_points)
            self._N_gp = None if N is None else np.ascontiguousarray(N, dtype=np.float64)
            self._N_gp_ready = True
        return self._N_gp

    def gauss_points_global(self):
        """Coordenadas globales de los puntos de Gauss, ``(N, n_gp, d)``:
        ``x_g = N(ξ_g)·X_e`` sobre la geometría de referencia. ``None`` si
        el elemento no declara funciones de forma."""
        N = self.gauss_shape_functions()
        if N is None:
            return None
        return np.einsum("gn,end->egd", N, self.X)

    def gauss_state(self, U: np.ndarray) -> dict:
        """``compute_gauss_state`` de toda la familia en una llamada.

        Evalúa la cinemática y la constitutiva en cada punto de Gauss desde
        el estado **committed** (como hace ``compute_gauss_state`` en el
        camino por elemento) sin tocar el trial. Devuelve arreglos con el
        eje de elemento delante, en el orden de ``self.elements``:

        - mecánica: ``strain`` y ``stress`` ``(N, n_gp, n_sigma)``;
        - térmica: ``grad_T`` y ``flux`` ``(N, n_gp, d)`` con ``q = −k·∇T``;
        - ``points_natural`` ``(n_gp, d)`` y ``points_global`` ``(N, n_gp, d)``
          (``None`` si el elemento no declara funciones de forma).

        La rodaja ``[i]`` de cada arreglo es exactamente lo que devolvería
        ``self.elements[i].compute_gauss_state(U)`` (salvo el último bit
        del orden de suma de ``B·u``).
        """
        N = self.n_elements
        n_gp = self.n_gp
        n_sig = self.n_sigma
        P = N * n_gp
        u = np.ascontiguousarray(U[self.dof_indices], dtype=np.float64)
        st = self.state
        eps = np.zeros((P, n_sig))
        sig = np.zeros((P, n_sig))
        scratch = np.empty_like(st.S_trial)
        flags = np.zeros(P, dtype=np.int8)
        self.gauss_kernel(
            self.kin_fn, self.mat_fn, self.X, u, self.gp_points,
            st.S_committed, scratch, self.mat_params, self.C, eps, sig, flags,
        )
        self._raise_if_bad_jacobian(flags)
        self._report_material_flags(flags)
        eps = eps.reshape(N, n_gp, n_sig)
        sig = sig.reshape(N, n_gp, n_sig)
        out = {
            "points_natural": self.gp_points.copy(),
            "points_global": self.gauss_points_global(),
        }
        if self.is_thermal:
            out["grad_T"] = eps
            out["flux"] = -sig
        else:
            out["strain"] = eps
            out["stress"] = sig
        return out

    def state_average(self, name: str | None):
        """Promedio sobre los puntos de Gauss de la variable interna escalar
        ``name`` del estado committed, ``(N,)``. Si ``name`` no es una
        variable escalar del esquema se toma la primera escalar; si no hay
        ninguna, ceros. Es la versión por arreglos del promedio de la
        ``PRIMARY_STATE_VAR`` que exporta el VTK."""
        schema = self.schema
        col = None
        if name is not None and name in schema.names:
            k = schema.names.index(name)
            if schema.shapes[k] == ():
                col = schema.offsets[k]
        if col is None:
            for k, shape in enumerate(schema.shapes):
                if shape == ():
                    col = schema.offsets[k]
                    break
        if col is None:
            return np.zeros(self.n_elements)
        return self.state.S_committed[:, col].reshape(self.n_elements, self.n_gp).mean(axis=1)


def build_families(elements: list, dof_indices: list[list[int]],
                   memory_budget: int | None = None,
                   parallel: bool | None = None) -> tuple[list[Family], list[int]]:
    """Agrupa los elementos batchables por clave de familia.

    Returns
    -------
    families, loose
        Familias construidas (en orden de primera aparición) y posiciones
        (índices en ``elements``) de los elementos que siguen el camino
        por elemento.
    """
    groups: dict[tuple, tuple[list, list[int]]] = {}
    loose: list[int] = []
    for pos, elem in enumerate(elements):
        if not element_is_batchable(elem):
            loose.append(pos)
            continue
        key = family_key(elem)
        if key not in groups:
            groups[key] = ([], [])
        groups[key][0].append(elem)
        groups[key][1].append(pos)

    families = []
    for key, (elems, positions) in groups.items():
        fam = Family(key, elems, positions, [dof_indices[p] for p in positions],
                     memory_budget=memory_budget, parallel=parallel)
        families.append(fam)
    if families:
        _log.debug(
            "Ensamblaje por lotes: %d familia(s), %d elemento(s) por lotes, %d por elemento.",
            len(families), sum(f.n_elements for f in families), len(loose),
        )
    return families, loose
