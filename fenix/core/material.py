# fenix_fem/fenix/core/material.py
from abc import ABC, abstractmethod
from typing import ClassVar, Optional

class Material(ABC):
    """
    Clase base abstracta para todos los materiales de Fenix FEM.

    Contrato obligatorio para subclases con variables internas
    ----------------------------------------------------------
    PRIMARY_STATE_VAR : ClassVar[Optional[str]]
        Nombre de la clave dentro del dict de state_vars que representa la
        variable interna principal (ej. 'alpha', 'damage').  Se usa en
        VtkExporter para exportar la variable de estado sin hardcodear nombres.
        Los materiales puramente elásticos deben dejarlo en ``None``.
    """

    PRIMARY_STATE_VAR: ClassVar[Optional[str]] = None

    @abstractmethod
    def compute_state(self, strain, state_vars=None):
        """
        Calcula el esfuerzo, el módulo tangente/secante y actualiza las variables internas.
        Retorna: (stress, tangent_modulus, new_state_vars)
        """
        pass
