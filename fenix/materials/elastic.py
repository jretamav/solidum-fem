# fenix_fem/fenix/materials/elastic.py
from fenix.core.material import Material
from fenix.registry import MaterialRegistry


@MaterialRegistry.register
class Elastic1D(Material):
    STRAIN_DIM = 1

    def __init__(self, E, density: float | None = None):
        self.E = E
        self.density = density

    def compute_state(self, strain, state_vars=None, **kwargs):
        """
        Calcula el esfuerzo y la rigidez tangente.
        Se aceptan state_vars y **kwargs para mantener la compatibilidad 
        con elementos no lineales que envían variables históricas.
        """
        stress = self.E * strain
        tangent = self.E
        
        # Retorna el esfuerzo, la tangente y las variables de estado intactas
        return stress, tangent, state_vars
