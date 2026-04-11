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

    def assemble_system(self):
        if self.domain.total_dofs == 0:
            self.domain.generate_equation_numbers()
        self.ndof = self.domain.total_dofs
        self.F_global = np.zeros(self.ndof)
        
        # Optimización: Preasignación exacta de memoria y tipado estricto
        total_entries = sum((len(e.get_dofs()))**2 for e in self.domain.elements.values())
        rows = np.zeros(total_entries, dtype=np.int32)
        cols = np.zeros(total_entries, dtype=np.int32)
        data = np.zeros(total_entries, dtype=np.float64)

        ptr = 0
        for elem_id, element in self.domain.elements.items():
            K_e = element.compute_global_stiffness()
            global_indices = []
            dofs_per_node = len(element.get_dofs()) // len(element.nodes)
            
            for i, node in enumerate(element.nodes):
                for dof_name in element.get_dofs()[i*dofs_per_node : (i+1)*dofs_per_node]:
                    global_indices.append(node.dofs[dof_name])
            
            # Llenado vectorizado: evita los bucles anidados lentos de Python
            n_idx = len(global_indices)
            rows[ptr:ptr + n_idx**2] = np.repeat(global_indices, n_idx)
            cols[ptr:ptr + n_idx**2] = np.tile(global_indices, n_idx)
            data[ptr:ptr + n_idx**2] = K_e.flatten()
            ptr += n_idx**2
                    
        self.K_global = sp.coo_matrix((data[:ptr], (rows[:ptr], cols[:ptr])), shape=(self.ndof, self.ndof)).tolil()

    def assemble_non_linear_system(self, U_current: np.ndarray):
        """Construye la matriz tangente global y el vector de fuerzas internas."""
        if self.domain.total_dofs == 0:
            self.domain.generate_equation_numbers()
            
        ndof = self.domain.total_dofs
        F_int_global = np.zeros(ndof)
        
        # Misma optimización de preasignación para el bucle no lineal
        total_entries = sum((len(e.get_dofs()))**2 for e in self.domain.elements.values())
        rows = np.zeros(total_entries, dtype=np.int32)
        cols = np.zeros(total_entries, dtype=np.int32)
        data = np.zeros(total_entries, dtype=np.float64)
        
        ptr = 0
        for elem_id, element in self.domain.elements.items():
            global_indices = []
            dofs_per_node = len(element.get_dofs()) // len(element.nodes)
            
            for i, node in enumerate(element.nodes):
                for dof_name in element.get_dofs()[i*dofs_per_node : (i+1)*dofs_per_node]:
                    global_indices.append(node.dofs[dof_name])
            
            u_local = U_current[global_indices]
            K_e, F_int_e = element.compute_element_state(u_local)
            
            # Ensamblaje del vector interno
            F_int_global[global_indices] += F_int_e
            
            # Ensamblaje de la matriz tangente
            n_idx = len(global_indices)
            rows[ptr:ptr + n_idx**2] = np.repeat(global_indices, n_idx)
            cols[ptr:ptr + n_idx**2] = np.tile(global_indices, n_idx)
            data[ptr:ptr + n_idx**2] = K_e.flatten()
            ptr += n_idx**2
                    
        K_global = sp.coo_matrix((data[:ptr], (rows[:ptr], cols[:ptr])), shape=(ndof, ndof)).tolil()
        return K_global, F_int_global

    def apply_bcs_to_system(self, K, R, penalty_value: float = 1e15, load_factor: float = 1.0, U_current: np.ndarray = None):
        """Función centralizada para aplicar condiciones de Dirichlet a cualquier sistema matricial."""
        for node_id, node in self.domain.nodes.items():
            for dof_name, prescribed_value in node.boundary_conditions.items():
                idx = node.dofs[dof_name]
                K[idx, idx] += penalty_value
                target_val = prescribed_value * load_factor
                if U_current is not None:
                    R[idx] = penalty_value * (target_val - U_current[idx])
                else:
                    R[idx] = penalty_value * target_val
        return K, R

    def apply_dirichlet_bcs(self):
        """Retrocompatibilidad para scripts que ensamblan el sistema estático lineal."""
        self.K_global, self.F_global = self.apply_bcs_to_system(self.K_global, self.F_global)

    def apply_point_load(self, node_id: int, dof_name: str, value: float):
        node = self.domain.get_node(node_id)
        if node and dof_name in node.dofs:
            idx = node.dofs[dof_name]
            self.F_global[idx] += value
