# fenix_fem/fenix/materials/elastic_2d.py
import numpy as np
from fenix.core.material import Material
from fenix.registry import MaterialRegistry


@MaterialRegistry.register
class Elastic2D(Material):
    STRAIN_DIM = 3

    _VALID_HYPOTHESES = frozenset({"plane_stress", "plane_strain"})

    def __init__(self, E: float, nu: float, hypothesis: str = 'plane_stress',
                 density: float | None = None):
        if E <= 0.0:
            raise ValueError(
                f"Elastic2D: E={E} debe ser estrictamente positivo."
            )
        # ν ∈ (−1, 0.5) por restricción termodinámica (módulo de Poisson
        # admisible para un sólido isótropo elástico). En el límite ν → 0.5
        # plane_strain diverge (incompresible).
        if not (-1.0 < nu < 0.5):
            raise ValueError(
                f"Elastic2D: nu={nu} fuera del rango admisible (-1, 0.5). "
                "El límite ν → 0.5 (incompresible) requiere formulación mixta "
                "(no implementada)."
            )
        if hypothesis not in self._VALID_HYPOTHESES:
            raise ValueError(
                f"Elastic2D: hypothesis={hypothesis!r} no reconocida. "
                f"Valores válidos: {sorted(self._VALID_HYPOTHESES)}."
            )
        if density is not None and density < 0.0:
            raise ValueError(
                f"Elastic2D: density={density} no puede ser negativa."
            )
        self.E, self.nu, self.hypothesis = float(E), float(nu), hypothesis
        self.density = density
        coef = self.E / (1.0 - self.nu**2) if hypothesis == 'plane_stress' else self.E / ((1.0 + self.nu)*(1.0 - 2.0*self.nu))
        if hypothesis == 'plane_stress':
            self.C = coef * np.array([[1.0, self.nu, 0.0], [self.nu, 1.0, 0.0], [0.0, 0.0, (1.0 - self.nu)/2.0]])
        else:
            self.C = coef * np.array([[1.0-self.nu, self.nu, 0.0], [self.nu, 1.0-self.nu, 0.0], [0.0, 0.0, (1.0-2.0*self.nu)/2.0]])

    def compute_state(self, strain: np.ndarray, state_vars=None):
        return self.C @ strain, self.C, state_vars
