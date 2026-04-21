# fenix_fem/fenix/materials/cable_1d.py
"""Material de cable 1D con elasticidad unilateral.

Ver docs/specs/CableMaterial1D.md para la especificación completa.
"""
from fenix.core.material import Material
from fenix.registry import MaterialRegistry


@MaterialRegistry.register
class CableMaterial1D(Material):
    """Cable axial 1D: elástico lineal en tracción, nulo en compresión.

    Respuesta memoryless gobernada por:

        σ(ε)   = E·ε   si ε > 0
                 0      si ε ≤ 0

        E_t(ε) = E     si ε > 0
                 0     si ε ≤ 0

    En ε = 0 la asignación al régimen compresivo evita introducir
    rigidez espuria cuando el cable está exactamente sin tensión.

    Parameters
    ----------
    E : float
        Módulo de Young en tracción. Debe ser estrictamente positivo.
    """

    STRAIN_DIM = 1
    PRIMARY_STATE_VAR = None

    def __init__(self, E):
        if E <= 0:
            raise ValueError(f"CableMaterial1D: E debe ser > 0, se recibió {E}.")
        self.E = E

    def compute_state(self, strain, state_vars=None, **kwargs):
        if strain > 0.0:
            stress = self.E * strain
            tangent = self.E
        else:
            stress = 0.0
            tangent = 0.0
        return stress, tangent, state_vars
