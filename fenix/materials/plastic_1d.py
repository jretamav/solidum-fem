# fenix_fem/fenix/materials/plastic_1d.py
from fenix.core.material import Material
from fenix.registry import MaterialRegistry


@MaterialRegistry.register
class Elastoplastic1D(Material):
    """
    Material elastoplástico 1D con endurecimiento isótropo lineal.
    Implementa el algoritmo de retorno (Return Mapping) clásico.
    """
    STRAIN_DIM = 1
    PRIMARY_STATE_VAR = 'alpha'  # deformación plástica acumulada equivalente

    def __init__(self, E: float, sigma_y: float, H: float = 0.0):
        self.E = E              # Módulo de Young
        self.sigma_y = sigma_y  # Esfuerzo de fluencia inicial
        self.H = H              # Módulo de endurecimiento (0 = plasticidad perfecta)

    def admissibility_scale(self, state_vars=None) -> float:
        """Esfuerzo de fluencia corriente ``σ_y + H·α`` (ADR 0006)."""
        alpha = 0.0 if state_vars is None else state_vars.get('alpha', 0.0)
        return self.sigma_y + self.H * alpha

    def compute_state(self, strain: float, state_vars=None):
        # 1. Recuperar memoria histórica
        eps_p_old = 0.0 if state_vars is None else state_vars.get('eps_p', 0.0)
        alpha_old = 0.0 if state_vars is None else state_vars.get('alpha', 0.0)

        # 2. Predictor Elástico (Estado de Prueba o Trial)
        eps_e_trial = strain - eps_p_old
        sigma_trial = self.E * eps_e_trial

        # 3. Función de Fluencia
        # f = |sigma_trial| - (sigma_y + H * alpha)
        yield_stress = self.sigma_y + self.H * alpha_old
        f = abs(sigma_trial) - yield_stress

        if self.is_admissible(f, state_vars):
            # Paso Elástico (Adentro de la superficie de fluencia)
            sigma = sigma_trial
            E_t = self.E
            new_state = {'eps_p': eps_p_old, 'alpha': alpha_old}
        else:
            # Paso Plástico (Corrección de Retorno / Return Mapping)
            delta_gamma = f / (self.E + self.H)
            sign_sigma = 1.0 if sigma_trial > 0 else -1.0

            # Actualización de esfuerzo
            sigma = sigma_trial - delta_gamma * self.E * sign_sigma
            
            # Evolución de variables internas
            eps_p_new = eps_p_old + delta_gamma * sign_sigma
            alpha_new = alpha_old + delta_gamma

            # Módulo tangente elastoplástico (Consistente)
            E_t = (self.E * self.H) / (self.E + self.H)

            new_state = {'eps_p': eps_p_new, 'alpha': alpha_new}

        return sigma, E_t, new_state
