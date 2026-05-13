# fenix_fem/fenix/math/assembly.py
from typing import Callable

import numpy as np
import scipy.sparse as sp

from fenix.bc.constraints import ConstraintSet
from fenix.core.domain import Domain
from fenix.logging import get_logger

_log = get_logger("assembly")


class Assembler:
    def __init__(self, domain: Domain):
        self.domain = domain
        self.ndof = 0
        self.K_global = None
        self.F_global = None

        # Variables para caché de topología COO
        self._topology_built = False
        self._coo_rows = None
        self._coo_cols = None
        self._total_entries = 0
        self._elem_dof_indices = []

        # Caché del ConstraintSet (ADR 0004 fase 1).
        self._constraint_set: ConstraintSet | None = None

    def _get_element_global_indices(self, element) -> list:
        """Extrae los índices globales de DOFs del elemento en el orden correcto.

        Delega en element.get_global_dof_indices(), que itera sobre DOF_NAMES × nodos
        sin asumir uniformidad de DOFs entre nodos.
        """
        return element.get_global_dof_indices()

    def _build_topology(self):
        """Precalcula y cachea la topología de ensamblaje (filas, columnas y mapeo) una sola vez."""
        if self.domain.total_dofs == 0:
            self.domain.generate_equation_numbers()

        self.ndof = self.domain.total_dofs
        self._total_entries = sum((len(e.DOF_NAMES) * len(e.nodes))**2 for e in self.domain.elements.values())

        self._coo_rows = np.zeros(self._total_entries, dtype=np.int32)
        self._coo_cols = np.zeros(self._total_entries, dtype=np.int32)
        self._elem_dof_indices = []

        ptr = 0
        for element in self.domain.elements.values():
            global_indices = self._get_element_global_indices(element)
            self._elem_dof_indices.append(global_indices)

            n_idx = len(global_indices)
            self._coo_rows[ptr:ptr + n_idx**2] = np.repeat(global_indices, n_idx)
            self._coo_cols[ptr:ptr + n_idx**2] = np.tile(global_indices, n_idx)
            ptr += n_idx**2

        self._topology_built = True

    def assemble_system(self):
        if not self._topology_built:
            self._build_topology()

        self.F_global = np.zeros(self.ndof)
        data = np.zeros(self._total_entries, dtype=np.float64)

        ptr = 0
        for i, element in enumerate(self.domain.elements.values()):
            K_e = element.compute_global_stiffness()
            n_idx = len(self._elem_dof_indices[i])

            data[ptr:ptr + n_idx**2] = K_e.ravel()
            ptr += n_idx**2

        # CSR: eficiente para solve y para modificar entradas diagonales existentes
        self.K_global = sp.coo_matrix(
            (data, (self._coo_rows, self._coo_cols)), shape=(self.ndof, self.ndof)
        ).tocsr()

    def assemble_non_linear_system(self, U_current: np.ndarray):
        """Construye la matriz tangente global y el vector de fuerzas internas a máxima velocidad."""
        if not self._topology_built:
            self._build_topology()

        F_int_global = np.zeros(self.ndof)
        data = np.zeros(self._total_entries, dtype=np.float64)

        ptr = 0
        for i, element in enumerate(self.domain.elements.values()):
            global_indices = self._elem_dof_indices[i]
            u_local = U_current[global_indices]

            K_e, F_int_e = element.compute_element_state(u_local)

            F_int_global[global_indices] += F_int_e

            n_idx = len(global_indices)
            data[ptr:ptr + n_idx**2] = K_e.ravel()
            ptr += n_idx**2

        K_global = sp.coo_matrix(
            (data, (self._coo_rows, self._coo_cols)), shape=(self.ndof, self.ndof)
        ).tocsr()
        return K_global, F_int_global

    # ------------------------------------------------------------------
    # Imposición de Dirichlet por eliminación directa (ADR 0004 fase 1).
    # ------------------------------------------------------------------

    @property
    def constraint_set(self) -> ConstraintSet:
        """ConstraintSet derivado de ``Node.boundary_conditions`` (lazy)."""
        if self._constraint_set is None:
            self._constraint_set = self._build_constraint_set()
        return self._constraint_set

    def _build_constraint_set(self) -> ConstraintSet:
        if not self._topology_built:
            self._build_topology()
        cs = ConstraintSet()
        for node in self.domain.nodes.values():
            for dof_name, value in node.boundary_conditions.items():
                cs.add_dirichlet(node.dofs[dof_name], value)
        for spec in getattr(self.domain, "linear_constraints", []):
            slave_node_id, slave_dof_name = spec["slave"]
            slave_node = self.domain.get_node(slave_node_id)
            if slave_node is None or slave_dof_name not in slave_node.dofs:
                raise ValueError(
                    f"Restricción lineal: nodo {slave_node_id} o DOF "
                    f"'{slave_dof_name}' inexistente."
                )
            slave = slave_node.dofs[slave_dof_name]
            masters: list[int] = []
            for nid, dn in spec["masters"]:
                node = self.domain.get_node(nid)
                if node is None or dn not in node.dofs:
                    raise ValueError(
                        f"Restricción lineal: nodo maestro {nid} o DOF "
                        f"'{dn}' inexistente."
                    )
                masters.append(node.dofs[dn])
            cs.add_linear(slave, masters, spec["coefficients"], spec.get("g", 0.0))
        return cs

    def reduce(
        self,
        K: sp.spmatrix,
        F: np.ndarray,
        U_current: np.ndarray | None = None,
        load_factor: float = 1.0,
    ):
        """Reduce ``K, F`` a la subred de DOFs libres (eliminación directa).

        Modelo afín ``u = T·u_libre + g`` con ``T`` sparse de shape
        ``(ndof, n_libre)`` que selecciona DOFs libres y, en filas esclavas,
        lleva los coeficientes ``α_si`` de las restricciones lineales. ``g``
        es no nulo solo en filas esclavas con término independiente.

        Para sistemas incrementales (Newton, corrector arc-length) se
        interpreta el resultado como ``δu = T·δu_libre + g_inc`` con
        ``g_inc = T·U_current[free] + load_factor·g_indep − U_current``,
        que en DOFs libres se anula y en DOFs esclavos representa la
        corrección necesaria para satisfacer la restricción.

        Returns
        -------
        K_red, F_red, T, g
            ``K_red = TᵀKT``, ``F_red = Tᵀ(F − K·g)``, el operador ``T`` y
            el vector ``g``. La reconstrucción se hace con :meth:`expand`.
        """
        cs = self.constraint_set
        T, g_indep = cs.build(self.ndof)

        if U_current is None:
            g = load_factor * g_indep
        else:
            free_dofs = cs.free_dofs(self.ndof)
            u_free = U_current[free_dofs]
            g = T @ u_free + load_factor * g_indep - U_current

        K_red = T.T @ K @ T
        F_red = T.T @ (F - K @ g)
        return K_red, F_red, T, g

    def expand(
        self,
        u_red: np.ndarray,
        T: sp.spmatrix,
        g: np.ndarray,
    ) -> np.ndarray:
        """Reconstruye ``u`` completo: ``u = T·u_red + g``."""
        return T @ u_red + g

    def commit_all_states(self):
        """Confirma las variables internas de todos los elementos tras la convergencia del paso."""
        for elem in self.domain.elements.values():
            elem.commit_state()

    def apply_point_load(self, node_id: int, dof_name: str, value: float):
        node = self.domain.get_node(node_id)
        if node and dof_name in node.dofs:
            idx = node.dofs[dof_name]
            self.F_global[idx] += value

    @staticmethod
    def _promote_to_3d(v, name: str) -> np.ndarray:
        """Acepta vector 2D ó 3D y promueve a R³ con padding cero si es 2D."""
        arr = np.asarray(v, dtype=float).ravel()
        if arr.size == 2:
            return np.array([arr[0], arr[1], 0.0])
        if arr.size == 3:
            return arr
        raise ValueError(
            f"{name}: dimensión inesperada ({arr.size}). Esperado 2 ó 3 componentes."
        )

    def _iterate_body_force(self, get_b_for_element: Callable) -> np.ndarray:
        """Recorrido interno común a `assemble_body_load` y `assemble_self_weight`.

        Itera sobre todos los elementos que declaran `compute_body_load` y
        delega en `get_b_for_element(element) -> np.ndarray` la obtención del
        vector de fuerza de cuerpo apropiado para cada uno (en sus
        componentes 2D ó 3D según `'uz' in element.DOF_NAMES`).
        """
        if not self._topology_built:
            self._build_topology()

        F_body = np.zeros(self.ndof)
        for i, element in enumerate(self.domain.elements.values()):
            method = getattr(element, "compute_body_load", None)
            if method is None:
                continue
            b_elem = get_b_for_element(element)
            f_e = method(b_elem)
            F_body[self._elem_dof_indices[i]] += f_e
        return F_body

    def assemble_body_load(self, b) -> np.ndarray:
        """Acumula la integral consistente de fuerza de cuerpo uniforme.

        Aplica el mismo vector `b` (fuerza por unidad de volumen, en ejes
        globales) a todos los elementos. Útil cuando el usuario precalcula
        `ρ·g` para una estructura monomaterial, o cuando la fuerza de cuerpo
        no proviene de gravedad (campo electromagnético, centrífuga, etc.).

        Para peso propio en estructuras multimaterial, preferir
        :meth:`assemble_self_weight` que usa `material.density` por elemento
        (ADR 0008).

        Promoción 2D↔3D: se acepta `b` con 2 ó 3 componentes; cada elemento
        recibe el slice apropiado según su dimensionalidad
        (`'uz' in DOF_NAMES` ⇒ 3D, en otro caso 2D). La omisión silenciosa
        de elementos sin `compute_body_load` es deliberada — el método es
        opcional, no parte del contrato base de `Element`.

        Parameters
        ----------
        b : np.ndarray or list, shape (2,) or (3,)
            Fuerza de cuerpo por unidad de volumen en ejes globales.

        Returns
        -------
        np.ndarray, shape (ndof,)
            Vector global de fuerzas nodales consistentes con `∫NᵀbA dx`.
        """
        b3 = self._promote_to_3d(b, "assemble_body_load")

        def _b_for_element(element):
            return b3 if 'uz' in element.DOF_NAMES else b3[:2]

        return self._iterate_body_force(_b_for_element)

    def assemble_self_weight(self, g) -> np.ndarray:
        """Acumula peso propio usando `material.density` de cada elemento (ADR 0008).

        Cada elemento ve `b_element = material.density · g`. Materiales que
        no declaran densidad (`density = 0.0`) producen aporte nulo y emiten
        un warning identificando el material concreto — probable olvido del
        usuario, no comportamiento físico legítimo (la masa no puede ser cero
        en peso propio).

        Promoción 2D↔3D: idéntica a `assemble_body_load`.

        Parameters
        ----------
        g : np.ndarray or list, shape (2,) or (3,)
            Vector aceleración gravitatoria en ejes globales (m/s² o las
            unidades correspondientes del problema). Típicamente
            `[0, -9.81]` en 2D ó `[0, 0, -9.81]` en 3D.

        Returns
        -------
        np.ndarray, shape (ndof,)
            Vector global de fuerzas nodales consistentes con peso propio.
        """
        g3 = self._promote_to_3d(g, "assemble_self_weight")
        warned_materials: set[int] = set()

        def _b_for_element(element):
            density = getattr(element.material, "density", 0.0)
            if density == 0.0 and id(element.material) not in warned_materials:
                warned_materials.add(id(element.material))
                _log.warning(
                    f"assemble_self_weight: material '{type(element.material).__name__}' "
                    f"tiene density=0.0; su aporte al peso propio será nulo. "
                    f"¿Olvidó declarar 'density' en YAML?"
                )
            g_elem = g3 if 'uz' in element.DOF_NAMES else g3[:2]
            return density * g_elem

        return self._iterate_body_force(_b_for_element)
