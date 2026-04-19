# fenix_fem/fenix/core/domain.py
from typing import Dict
from fenix.core.node import Node

class Domain:
    def __init__(self):
        self.nodes: Dict[int, Node] = {}
        self.elements = {}
        self.total_dofs: int = 0

    def add_node(self, node_id: int, coordinates: list[float]) -> Node:
        if node_id in self.nodes:
            raise ValueError(f"El nodo {node_id} ya existe.")
        new_node = Node(node_id, coordinates)
        self.nodes[node_id] = new_node
        return new_node

    def get_node(self, node_id: int) -> Node:
        return self.nodes.get(node_id)

    def add_element(self, element) -> None:
        if element.id in self.elements:
            raise ValueError(f"El elemento {element.id} ya existe.")
        self.elements[element.id] = element

    def generate_equation_numbers(self, verbose: bool = False):
        eq_number = 0
        for node_id, node in self.nodes.items():
            for dof_name in node.dofs.keys():
                node.dofs[dof_name] = eq_number
                eq_number += 1
        self.total_dofs = eq_number
        if verbose:
            print(f"Numeración completada. Grados de libertad totales: {self.total_dofs}")
