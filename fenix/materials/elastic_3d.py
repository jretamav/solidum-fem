# fenix_fem/fenix/materials/elastic_3d.py
"""Elasticidad lineal isótropa 3D (ADR 0012, spec ``docs/specs/Elastic3D.md``).

Sin variantes de hipótesis (en 3D no aplican plane_stress/plane_strain).
Convención Voigt 6D del proyecto: ``[xx, yy, zz, xy, yz, xz]`` con
``γ_ij = 2·ε_ij`` *engineering*.
"""
import numpy as np

from fenix.core.material import Material
from fenix.registry import MaterialRegistry


@MaterialRegistry.register
class Elastic3D(Material):
    """Material elástico lineal isótropo en 3D.

    Parameters
    ----------
    E : float
        Módulo de Young (>0).
    nu : float
        Coeficiente de Poisson, ν ∈ (-1, 0.5). El límite ν → 0.5
        (incompresible) requiere formulación mixta u-p, no soportada.
    density : float, optional
        Densidad (kg/m³). Opcional al construir; obligatoria solo si se
        ensambla peso propio o matriz de masa (ADR 0008).
    """

    STRAIN_DIM = 6

    def __init__(self, E: float, nu: float, density: float | None = None):
        if E <= 0.0:
            raise ValueError(
                f"Elastic3D: E={E} debe ser estrictamente positivo."
            )
        if not (-1.0 < nu < 0.5):
            raise ValueError(
                f"Elastic3D: nu={nu} fuera del rango admisible (-1, 0.5). "
                "El límite ν → 0.5 (incompresible) requiere formulación "
                "mixta u-p (no implementada)."
            )
        if density is not None and density < 0.0:
            raise ValueError(
                f"Elastic3D: density={density} no puede ser negativa."
            )

        self.E = float(E)
        self.nu = float(nu)
        self.density = density

        coef = self.E / ((1.0 + self.nu) * (1.0 - 2.0 * self.nu))
        shear = (1.0 - 2.0 * self.nu) / 2.0
        a = 1.0 - self.nu
        b = self.nu
        self.C = coef * np.array([
            [a, b, b, 0.0, 0.0, 0.0],
            [b, a, b, 0.0, 0.0, 0.0],
            [b, b, a, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, shear, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, shear, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, shear],
        ])

    def compute_state(self, strain: np.ndarray, state_vars=None):
        return self.C @ strain, self.C, state_vars
