# fenix_fem/fenix/core/element.py
from abc import ABC, abstractmethod
from typing import List
from fenix.core.node import Node

class Element(ABC):
    def __init__(self, element_id: int, nodes: List[Node]):
        self.id = element_id
        self.nodes = nodes

    @abstractmethod
    def compute_global_stiffness(self):
        pass

    @abstractmethod
    def get_dofs(self) -> List[str]:
        pass
