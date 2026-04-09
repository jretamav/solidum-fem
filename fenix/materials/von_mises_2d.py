# fenix_fem/fenix/materials/von_mises_2d.py
import numpy as np
import math
from fenix.core.material import Material
from numba import njit

@njit
def _compute_j2_plasticity(strain, eps_p_old, alpha_old, sigma_y, H, K, G, C_e):
    # Deformación volumétrica (traza) y deformación desviadora
    eps_v = strain[0] + strain[1] 
    
    e_dev = np.array([
        strain[0] - eps_v / 3.0,
        strain[1] - eps_v / 3.0,
        -eps_v / 3.0,
        strain[2] / 2.0
    ])
    
    # Predictor Elástico
    e_dev_trial = e_dev - eps_p_old
    s_trial = 2.0 * G * e_dev_trial
    norm_s_trial = math.sqrt(s_trial[0]**2 + s_trial[1]**2 + s_trial[2]**2 + 2.0 * s_trial[3]**2)
    
    yield_stress = sigma_y + H * alpha_old
    f_trial = norm_s_trial - math.sqrt(2.0 / 3.0) * yield_stress
    
    if f_trial <= 1e-9:
        sigma = C_e @ strain
        return sigma, C_e.copy(), eps_p_old.copy(), alpha_old
        
    # Corrector Plástico (Return Mapping)
    delta_gamma = f_trial / (2.0 * G + (2.0 / 3.0) * H)
    N = s_trial / norm_s_trial
    
    s_new = s_trial - 2.0 * G * delta_gamma * N
    eps_p_new = eps_p_old + delta_gamma * N
    alpha_new = alpha_old + math.sqrt(2.0 / 3.0) * delta_gamma
    
    p = K * eps_v
    sigma = np.array([
        s_new[0] + p,
        s_new[1] + p,
        s_new[3]
    ])
    
    # Matriz Tangente Algorítmica Consistente
    beta = (2.0 * G * delta_gamma) / norm_s_trial
    gamma_factor = 1.0 / (1.0 + H / (3.0 * G)) - beta
    
    v = np.array([1.0, 1.0, 0.0])
    I_dev = np.array([
        [ 2.0/3.0, -1.0/3.0, 0.0],
        [-1.0/3.0,  2.0/3.0, 0.0],
        [ 0.0,      0.0,     0.5]
    ])
    
    N_voigt = np.array([N[0], N[1], N[3]])
    N_otimes_N = np.array([
        [N_voigt[0]*N_voigt[0], N_voigt[0]*N_voigt[1], N_voigt[0]*N_voigt[2]],
        [N_voigt[1]*N_voigt[0], N_voigt[1]*N_voigt[1], N_voigt[1]*N_voigt[2]],
        [N_voigt[2]*N_voigt[0], N_voigt[2]*N_voigt[1], N_voigt[2]*N_voigt[2]*0.5]
    ])
    
    C_alg = K * np.outer(v, v) + 2.0 * G * (1.0 - beta) * I_dev - 2.0 * G * gamma_factor * N_otimes_N
    
    return sigma, C_alg, eps_p_new, alpha_new

class VonMises2D(Material):
    """
    Modelo de plasticidad J2 (Von Mises) con endurecimiento isotrópico.
    Formulación estricta para Deformación Plana usando Return Mapping.
    """
    def __init__(self, E: float, nu: float, sigma_y: float, H: float = 0.0):
        self.E = E
        self.nu = nu
        self.sigma_y = sigma_y
        self.H = H
        
        # Módulos elásticos de volumen (K) y corte (G)
        self.K = E / (3.0 * (1.0 - 2.0 * nu))
        self.G = E / (2.0 * (1.0 + nu))
        
        # Matriz elástica en notación Voigt [xx, yy, xy]
        coef = E / ((1.0 + nu) * (1.0 - 2.0 * nu))
        self.C_e = coef * np.array([
            [1.0 - nu, nu,       0.0],
            [nu,       1.0 - nu, 0.0],
            [0.0,      0.0,      (1.0 - 2.0 * nu) / 2.0]
        ])

    def compute_state(self, strain: np.ndarray, state_vars=None):
        # 1. Recuperar variables internas históricas
        # eps_p_old almacena las deformaciones plásticas tensoriales [xx, yy, zz, xy]
        eps_p_old = np.zeros(4) if state_vars is None else state_vars.get('eps_p', np.zeros(4))
        alpha_old = 0.0 if state_vars is None else state_vars.get('alpha', 0.0)

        sigma, C_alg, eps_p_new, alpha_new = _compute_j2_plasticity(
            strain, eps_p_old, alpha_old, self.sigma_y, self.H, self.K, self.G, self.C_e
        )
        
        new_state = {'eps_p': eps_p_new, 'alpha': alpha_new}
        return sigma, C_alg, new_state
