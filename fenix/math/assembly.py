# fenix_fem/fenix/math/assembly.py
import numpy as np
import scipy.sparse as sp

from fenix.bc.constraints import ConstraintSet
from fenix.core.domain import Domain


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
        return cs

    def reduce(
        self,
        K: sp.spmatrix,
        F: np.ndarray,
        U_current: np.ndarray | None = None,
        load_factor: float = 1.0,
    ):
        """Reduce ``K, F`` a la subred de DOFs libres (eliminación directa).

        Modelo afín ``u = T·u_libre + g`` con ``T`` rectangular que selecciona
        DOFs libres y ``g`` no nulo solo en DOFs prescritos. Para sistemas
        incrementales (Newton) se interpreta ``u`` como ``δu`` y ``g`` como
        ``δu`` impuesto en los esclavos: ``g[bc] = load_factor·valor − U_current[bc]``.

        Parameters
        ----------
        K, F
            Sistema completo ``K·u = F``.
        U_current
            Si se pasa, ``g`` se calcula como incremento desde ``U_current``
            (uso típico: Newton, arc-length corrector). ``None`` para resolver
            ``u`` absoluto (lineal, predictor).
        load_factor
            Escala los valores prescritos. Útil para incrementos en
            Newton-Raphson y arc-length.

        Returns
        -------
        K_red, F_red, free_dofs, g_full
            ``K_red = T^T K T``, ``F_red = T^T (F − K·g_full)``, los índices
            de DOFs libres y el vector ``g_full ∈ ℝ^n``. La reconstrucción se
            hace con :meth:`expand`.
        """
        cs = self.constraint_set
        bc_dofs = cs.slave_dofs
        bc_vals = cs.slave_values
        free_dofs = cs.free_dofs(self.ndof)

        g_full = np.zeros(self.ndof)
        if bc_dofs.size > 0:
            target = bc_vals * load_factor
            if U_current is not None:
                g_full[bc_dofs] = target - U_current[bc_dofs]
            else:
                g_full[bc_dofs] = target

        K_red = K[free_dofs, :][:, free_dofs]
        if bc_dofs.size > 0 and np.any(g_full):
            F_red = F[free_dofs] - (K @ g_full)[free_dofs]
        else:
            F_red = F[free_dofs].copy()

        return K_red, F_red, free_dofs, g_full

    def expand(
        self,
        u_red: np.ndarray,
        free_dofs: np.ndarray,
        g_full: np.ndarray,
    ) -> np.ndarray:
        """Reconstruye ``u`` completo: ``u = T·u_red + g_full``."""
        u_full = g_full.copy()
        u_full[free_dofs] = u_red
        return u_full

    def commit_all_states(self):
        """Confirma las variables internas de todos los elementos tras la convergencia del paso."""
        for elem in self.domain.elements.values():
            elem.commit_state()

    def apply_point_load(self, node_id: int, dof_name: str, value: float):
        node = self.domain.get_node(node_id)
        if node and dof_name in node.dofs:
            idx = node.dofs[dof_name]
            self.F_global[idx] += value
