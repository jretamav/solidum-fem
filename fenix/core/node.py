# fenix_fem/fenix/core/node.py
class Node:
    def __init__(self, node_id: int, coordinates: list[float]):
        self.id = node_id
        self.coordinates = coordinates
        self.dofs = {}
        self.boundary_conditions = {}

    def add_dof(self, dof_name: str):
        if dof_name not in self.dofs:
            self.dofs[dof_name] = -1

    def fix_dof(self, dof_name: str, value: float = 0.0):
        if dof_name in self.dofs:
            self.boundary_conditions[dof_name] = value
        else:
            raise ValueError(f"El DoF '{dof_name}' no existe en el nodo {self.id}.")
