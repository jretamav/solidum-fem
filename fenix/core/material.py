# fenix_fem/fenix/core/material.py
from abc import ABC, abstractmethod
from typing import ClassVar, Optional


class Material(ABC):
    """Clase base abstracta para todos los materiales de Fenix FEM.

    Contrato de subclases
    ---------------------
    STRAIN_DIM : ClassVar[int]
        Dimensión del vector de deformación que el material espera/produce:
            1  → escalar 1D (Truss, Frame)
            3  → Voigt 2D plano [εxx, εyy, γxy]
            6  → Voigt 3D       [εxx, εyy, εzz, γxy, γyz, γxz]
        Los elementos validan al construirse que `material.STRAIN_DIM`
        coincide con la dimensión que ellos requieren, atajando errores
        que de otro modo aparecerían como `matmul` críptico en runtime.

    PRIMARY_STATE_VAR : ClassVar[Optional[str]]
        Nombre de la clave principal dentro del dict de state_vars
        (ej. 'alpha', 'damage'). Se usa en VtkExporter para exportar
        la variable de estado sin hardcodear nombres. `None` para
        materiales puramente elásticos sin variables internas.
    """

    STRAIN_DIM: ClassVar[int]
    PRIMARY_STATE_VAR: ClassVar[Optional[str]] = None

    @abstractmethod
    def compute_state(self, strain, state_vars=None):
        """Calcula esfuerzo, módulo tangente/secante y nuevas variables internas.

        Returns
        -------
        (stress, tangent_modulus, new_state_vars)
        """
        pass
