import copy
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
        """Fija las variables de estado trial cuando el paso de carga converge.

        Usa copia profunda para que el estado comprometido y el trial sean
        completamente independientes, evitando corrupción por modificaciones in-place.
        """
        self.vars = copy.deepcopy(self.vars_trial)
        self.stresses = [np.copy(s) if s is not None else None for s in self.stresses_trial]

    # ------------------------------------------------------------------
    # Aislamiento del trial: evaluaciones que no deben dejar huella
    # ------------------------------------------------------------------

    def snapshot_trial(self):
        """Copia independiente del estado trial (variables y esfuerzos),
        para restaurarlo con :meth:`restore_trial` tras una evaluación
        auxiliar (p. ej. ``compute_global_stiffness`` en ``u = 0``) que no
        forma parte del paso en curso."""
        return (copy.deepcopy(self.vars_trial),
                [np.copy(s) if s is not None else None for s in self.stresses_trial])

    def restore_trial(self, snapshot) -> None:
        vars_trial, stresses_trial = snapshot
        self.vars_trial = vars_trial
        self.stresses_trial = stresses_trial