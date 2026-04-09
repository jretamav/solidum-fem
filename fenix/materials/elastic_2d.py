# fenix_fem/fenix/materials/elastic_2d.py
import numpy as np
from fenix.core.material import Material

class Elastic2D(Material):
    def __init__(self, E: float, nu: float, hypothesis: str = 'plane_stress'):
        self.E, self.nu, self.hypothesis = E, nu, hypothesis
        coef = self.E / (1.0 - self.nu**2) if hypothesis == 'plane_stress' else self.E / ((1.0 + self.nu)*(1.0 - 2.0*self.nu))
        if hypothesis == 'plane_stress':
            self.C = coef * np.array([[1.0, self.nu, 0.0], [self.nu, 1.0, 0.0], [0.0, 0.0, (1.0 - self.nu)/2.0]])
        else:
            self.C = coef * np.array([[1.0-self.nu, self.nu, 0.0], [self.nu, 1.0-self.nu, 0.0], [0.0, 0.0, (1.0-2.0*self.nu)/2.0]])

    def compute_state(self, strain: np.ndarray, state_vars=None):
        return self.C @ strain, self.C, state_vars
