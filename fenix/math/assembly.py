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
        
        rows, cols, data = [], [], []

        for elem_id, element in self.domain.elements.items():
            K_e = element.compute_global_stiffness()
            global_indices = []
            dofs_per_node = len(element.get_dofs()) // len(element.nodes)
            
            for i, node in enumerate(element.nodes):
                for dof_name in element.get_dofs()[i*dofs_per_node : (i+1)*dofs_per_node]:
                    global_indices.append(node.dofs[dof_name])
            
            for i in range(len(global_indices)):
                for j in range(len(global_indices)):
                    rows.append(global_indices[i])
                    cols.append(global_indices[j])
                    data.append(K_e[i, j])
                    
        self.K_global = sp.coo_matrix((data, (rows, cols)), shape=(self.ndof, self.ndof)).tolil()

    def assemble_non_linear_system(self, U_current: np.ndarray):
        """Construye la matriz tangente global y el vector de fuerzas internas."""
        if self.domain.total_dofs == 0:
            self.domain.generate_equation_numbers()
            
        ndof = self.domain.total_dofs
        F_int_global = np.zeros(ndof)
        
        rows, cols, data = [], [], []
        
        for elem_id, element in self.domain.elements.items():
            global_indices = []
            dofs_per_node = len(element.get_dofs()) // len(element.nodes)
            
            for i, node in enumerate(element.nodes):
                for dof_name in element.get_dofs()[i*dofs_per_node : (i+1)*dofs_per_node]:
                    global_indices.append(node.dofs[dof_name])
            
            u_local = U_current[global_indices]
            K_e, F_int_e = element.compute_element_state(u_local)
            
            for i in range(len(global_indices)):
                F_int_global[global_indices[i]] += F_int_e[i]
                for j in range(len(global_indices)):
                    rows.append(global_indices[i])
                    cols.append(global_indices[j])
                    data.append(K_e[i, j])
                    
        K_global = sp.coo_matrix((data, (rows, cols)), shape=(ndof, ndof)).tolil()
        return K_global, F_int_global

    def apply_dirichlet_bcs(self):
        penalty = 1e15
        for node_id, node in self.domain.nodes.items():
            for dof_name, prescribed_value in node.boundary_conditions.items():
                idx = node.dofs[dof_name]
                self.K_global[idx, idx] += penalty
                self.F_global[idx] += penalty * prescribed_value

    def apply_point_load(self, node_id: int, dof_name: str, value: float):
        node = self.domain.get_node(node_id)
        if node and dof_name in node.dofs:
            idx = node.dofs[dof_name]
            self.F_global[idx] += value
