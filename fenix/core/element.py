# fenix_fem/fenix/core/element.py
from abc import ABC, abstractmethod
from typing import List
import numpy as np
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

    def get_coordinate_matrix(self, ndim: int = 2) -> np.ndarray:
        """Extrae de forma elegante las coordenadas de los nodos en una matriz de NumPy."""
        coords = np.zeros((len(self.nodes), ndim))
        for i, node in enumerate(self.nodes):
            for j in range(min(ndim, len(node.coordinates))):
                coords[i, j] = node.coordinates[j]
        return coords
