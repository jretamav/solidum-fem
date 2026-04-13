import numpy as np

class ElementState:
    """
    Gestiona las variables históricas y de estado interno en los puntos de integración
    de un elemento finito, resolviendo violaciones al SRP.
    """
    def __init__(self, num_ip: int, init_stress=None):
        self.num_ip = num_ip
        self.vars = [None] * num_ip
        self.vars_trial = [None] * num_ip
        
        self.stresses = [np.copy(init_stress) if init_stress is not None else None for _ in range(num_ip)]
        self.stresses_trial = [np.copy(init_stress) if init_stress is not None else None for _ in range(num_ip)]

    def commit(self):
        """Fija las variables de estado trial cuando el paso de carga converge."""
        self.vars = self.vars_trial.copy()
        self.stresses = self.stresses_trial.copy()