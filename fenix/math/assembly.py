# fenix_fem/fenix/math/assembly.py
import numpy as np
import scipy.sparse as sp
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
        
        # Caché para Condiciones de Frontera Vectorizadas
        self._bc_built = False
        self._bc_dofs = np.array([], dtype=np.int32)
        self._bc_vals = np.array([], dtype=np.float64)

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

    def _build_bc_arrays(self):
        """Precalcula los DOFs sujetos a condiciones de frontera para evitar bucles."""
        dofs = []
        vals = []
        for node in self.domain.nodes.values():
            for dof_name, prescribed_value in node.boundary_conditions.items():
                dofs.append(node.dofs[dof_name])
                vals.append(prescribed_value)
        self._bc_dofs = np.array(dofs, dtype=np.int32)
        self._bc_vals = np.array(vals, dtype=np.float64)
        self._bc_built = True

    def apply_bcs_to_system(self, K, R, penalty_value: float = 1e15, load_factor: float = 1.0, U_current: np.ndarray = None):
        """Aplica condiciones de Dirichlet usando matrices dispersas vectorizadas (Cero overhead)."""
        if not getattr(self, '_bc_built', False):
            self._build_bc_arrays()
            
        if len(self._bc_dofs) > 0:
            # Sumar matriz diagonal dispersa es instantáneo en C
            penalty_array = np.zeros(self.ndof)
            penalty_array[self._bc_dofs] = penalty_value
            K = K + sp.diags(penalty_array, format='csr')
            
            target_vals = self._bc_vals * load_factor
            if U_current is not None:
                R[self._bc_dofs] = penalty_value * (target_vals - U_current[self._bc_dofs])
            else:
                R[self._bc_dofs] = penalty_value * target_vals
                
        return K, R

    def commit_all_states(self):
        """Confirma las variables internas de todos los elementos tras la convergencia del paso."""
        for elem in self.domain.elements.values():
            elem.commit_state()

    def apply_dirichlet_bcs(self):
        """Retrocompatibilidad para scripts que ensamblan el sistema estático lineal."""
        self.K_global, self.F_global = self.apply_bcs_to_system(self.K_global, self.F_global)

    def apply_point_load(self, node_id: int, dof_name: str, value: float):
        node = self.domain.get_node(node_id)
        if node and dof_name in node.dofs:
            idx = node.dofs[dof_name]
            self.F_global[idx] += value
