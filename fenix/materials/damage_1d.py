# fenix_fem/fenix/materials/damage_1d.py
import math
from fenix.core.material import Material

class IsotropicDamage1D(Material):
    """
    Modelo de daño isotrópico continuo unidimensional con ley de ablandamiento exponencial.
    Ideal para ser utilizado con elementos de armadura (Truss2D o Truss3D).
    """
    def __init__(self, E: float, kappa_0: float, alpha: float):
        self.E = E              # Módulo de Young intacto
        self.kappa_0 = kappa_0  # Umbral de deformación elástica inicial
        self.alpha = alpha      # Parámetro que controla la velocidad de degradación (ablandamiento)

    def compute_state(self, strain: float, state_vars=None):
        # 1. Recuperar variable histórica (máxima deformación alcanzada kappa)
        kappa_old = self.kappa_0 if state_vars is None else state_vars.get('kappa', self.kappa_0)
        
        # 2. Calcular deformación equivalente (valor absoluto en 1D)
        eps_eq = abs(strain)
        
        # 3. Evolución del daño (Condición de Kuhn-Tucker)
        kappa_new = max(kappa_old, eps_eq)
        
        if kappa_new <= self.kappa_0:
            d = 0.0
        else:
            d = 1.0 - (self.kappa_0 / kappa_new) * math.exp(-self.alpha * (kappa_new - self.kappa_0))
            d = min(d, 0.999)  # Evitar singularidad numérica en la matriz de rigidez
            
        # 4. Esfuerzo y módulo secante
        E_sec = (1.0 - d) * self.E
        sigma = E_sec * strain
        
        new_state = {'kappa': kappa_new, 'damage': d}
        return sigma, E_sec, new_state