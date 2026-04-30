# fenix_fem/fenix/materials/damage_2d.py
import numpy as np
from fenix.core.material import Material
from fenix.materials.elastic_2d import Elastic2D
from fenix.constants import DAMAGE_MAX
from fenix.registry import MaterialRegistry


@MaterialRegistry.register
class IsotropicDamage2D(Material):
    """
    Modelo de daño isotrópico continuo con ley de ablandamiento exponencial.
    """
    STRAIN_DIM = 3
    PRIMARY_STATE_VAR = 'damage'  # variable de daño en [0, 1]

    def __init__(self, E: float, nu: float, kappa_0: float, alpha: float, hypothesis: str = 'plane_stress'):
        self.elastic_base = Elastic2D(E, nu, hypothesis)
        self.kappa_0 = kappa_0
        self.alpha = alpha

    def compute_state(self, strain: np.ndarray, state_vars=None):
        # 1. Recuperar variable histórica (máxima deformación alcanzada kappa)
        kappa_old = self.kappa_0 if state_vars is None else state_vars.get('kappa', self.kappa_0)
        
        # 2. Calcular deformación equivalente (norma simple de deformación)
        eps_eq = np.sqrt(strain[0]**2 + strain[1]**2 + 0.5 * strain[2]**2)
        
        # 3. Evolución del daño (Condición de carga/descarga de Kuhn-Tucker)
        kappa_new = max(kappa_old, eps_eq)
        
        if kappa_new <= self.kappa_0:
            d = 0.0
        else:
            d = 1.0 - (self.kappa_0 / kappa_new) * np.exp(-self.alpha * (kappa_new - self.kappa_0))
            d = min(d, DAMAGE_MAX) # Evitar singularidad numérica
            
        # 4. Esfuerzo y tensor constitutivo secante
        Ce = self.elastic_base.C
        C_sec = (1.0 - d) * Ce
        sigma = C_sec @ strain
        
        new_state = {'kappa': kappa_new, 'damage': d}
        return sigma, C_sec, new_state
