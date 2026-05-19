# fenix_fem/fenix/materials/elastic.py
from fenix.core.material import Material
from fenix.registry import MaterialRegistry


@MaterialRegistry.register
class Elastic1D(Material):
    STRAIN_DIM = 1

    def __init__(self, E: float, density: float | None = None):
        if E <= 0.0:
            raise ValueError(
                f"Elastic1D: E={E} debe ser estrictamente positivo."
            )
        if density is not None and density < 0.0:
            raise ValueError(
                f"Elastic1D: density={density} no puede ser negativa."
            )
        self.E = float(E)
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
