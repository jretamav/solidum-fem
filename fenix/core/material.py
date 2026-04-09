# fenix_fem/fenix/core/material.py
from abc import ABC, abstractmethod

class Material(ABC):
    @abstractmethod
    def compute_state(self, strain, state_vars=None):
        """
        Calcula el esfuerzo, el módulo tangente/secante y actualiza las variables internas.
        Retorna: (stress, tangent_modulus, new_state_vars)
        """
        pass
