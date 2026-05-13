# fenix_fem/fenix/math/solvers/newmark.py
"""``NewmarkSolver`` — integración Newmark-β para ``M·ü + C·u̇ + K·u = F(t)``
(ADR 0009 fase 3).
"""
from typing import Callable, Optional

import numpy as np
import scipy.sparse as sp

from fenix.math.damping import rayleigh_from_modes
from fenix.math.linalg import StiffnessProperties, select_solver
from fenix.math.solvers._shared import _log, domain_is_symmetric
from fenix.registry import SolverRegistry
from fenix.results import TransientResult


@SolverRegistry.register
class NewmarkSolver:
    """Análisis dinámico transitorio lineal por integración Newmark-β (ADR 0009 fase 3).

    Resuelve ``M·ü + C·u̇ + K·u = F(t)`` con amortiguamiento Rayleigh
    ``C = α·M + β·K`` y apoyos Dirichlet constantes en el tiempo. ``K`` y
    ``M`` se asumen constantes (problema lineal); la formulación se
    reduce a DOFs libres por eliminación directa (ADR 0004) y la matriz
    efectiva ``A_eff = M_red + γΔt·C_red + βΔt²·K_red`` se factoriza una
    sola vez al inicio (ADR 0003), por lo que cada paso temporal es una
    resolución triangular barata.

    Parameters
    ----------
    assembler : Assembler
        Vinculación al modelo.
    t_end : float
        Tiempo final (s).
    dt : float
        Paso temporal constante (s).
    beta, gamma : float, default 0.25, 0.5
        Parámetros de la familia Newmark. Default = average acceleration,
        incondicionalmente estable y sin amortiguamiento numérico.
    rayleigh : dict or None, default None
        Amortiguamiento Rayleigh ``C = α·M + β·K``. Tres formas aceptadas:
        - ``None``: sin amortiguamiento (α = β = 0).
        - ``{"alpha": ..., "beta": ...}``: coeficientes directos.
        - ``{"xi1": ..., "omega1": ..., "xi2": ..., "omega2": ...}``:
          calibración modal automática.
    u0, u0_dot : np.ndarray or None
        Condiciones iniciales globales (shape ``(ndof,)``). Cero por
        defecto. Si se pasan, deben ser consistentes con los apoyos
        Dirichlet en los DOFs prescritos.
    F_func : callable or None
        ``F_func(t: float) -> np.ndarray`` devuelve el vector global de
        fuerzas externas en el instante ``t``. ``None`` ⇒ vibración libre.
    linear_algebra : str, default "auto"
        Backend para factorizar ``A_eff_red`` (ADR 0003).
    lumping : str, default "consistent"
        Discretización de masa. Solo ``"consistent"`` en fase 1.
    """

    def __init__(
        self,
        assembler,
        t_end: float,
        dt: float,
        *,
        beta: float = 0.25,
        gamma: float = 0.5,
        rayleigh: Optional[dict] = None,
        u0: Optional[np.ndarray] = None,
        u0_dot: Optional[np.ndarray] = None,
        F_func: Optional[Callable[[float], np.ndarray]] = None,
        linear_algebra: str = "auto",
        lumping: str = "consistent",
    ):
        if dt <= 0.0:
            raise ValueError(f"NewmarkSolver: dt={dt} debe ser positivo.")
        if t_end <= 0.0:
            raise ValueError(f"NewmarkSolver: t_end={t_end} debe ser positivo.")
        self.assembler = assembler
        self.t_end = float(t_end)
        self.dt = float(dt)
        self.beta = float(beta)
        self.gamma = float(gamma)
        self.rayleigh_cfg = rayleigh
        self.u0 = u0
        self.u0_dot = u0_dot
        self.F_func = F_func
        self.linear_algebra = str(linear_algebra)
        self.lumping = str(lumping)

    @staticmethod
    def _resolve_rayleigh(cfg: Optional[dict]) -> tuple[float, float]:
        """Traduce el dict del usuario a ``(α, β)``. Tres formas — ver __init__."""
        if cfg is None:
            return 0.0, 0.0
        if not isinstance(cfg, dict):
            raise ValueError(
                f"NewmarkSolver.rayleigh: esperado dict o None, recibido "
                f"{type(cfg).__name__}."
            )
        if "alpha" in cfg and "beta" in cfg:
            return float(cfg["alpha"]), float(cfg["beta"])
        required = {"xi1", "omega1", "xi2", "omega2"}
        if required.issubset(cfg):
            return rayleigh_from_modes(
                float(cfg["xi1"]), float(cfg["omega1"]),
                float(cfg["xi2"]), float(cfg["omega2"]),
            )
        raise ValueError(
            "NewmarkSolver.rayleigh: dict admitido es {'alpha', 'beta'} o "
            "{'xi1','omega1','xi2','omega2'}; recibido " + repr(set(cfg.keys()))
        )

    def solve(self) -> TransientResult:
        _log.info("--- INICIANDO SOLVER NEWMARK ---")

        # Ensamblar K (lineal en u=0) y M (consistente, cacheada).
        self.assembler.assemble_system()
        K = self.assembler.K_global
        M = self.assembler.assemble_mass_matrix(lumping=self.lumping)

        # Amortiguamiento Rayleigh: C = α·M + β·K.
        alpha_r, beta_r = self._resolve_rayleigh(self.rayleigh_cfg)
        C = alpha_r * M + beta_r * K

        # Reducción por Dirichlet (ADR 0004). T selecciona DOFs libres;
        # g recoge valores prescritos constantes en el tiempo.
        cs = self.assembler.constraint_set
        T, g = cs.build(self.assembler.ndof)
        free_dofs = cs.free_dofs(self.assembler.ndof)

        K_red = (T.T @ K @ T).tocsr()
        M_red = (T.T @ M @ T).tocsr()
        C_red = (T.T @ C @ T).tocsr() if (alpha_r != 0.0 or beta_r != 0.0) \
                else sp.csr_matrix(K_red.shape)
        # Término constante por apoyos prescritos no nulos.
        F_dir = T.T @ (K @ g)

        ndof = self.assembler.ndof
        n_free = K_red.shape[0]

        # Condiciones iniciales en globales → proyección a DOFs libres.
        u0_global = (np.zeros(ndof) if self.u0 is None
                      else np.asarray(self.u0, dtype=float).reshape(ndof))
        u0_dot_global = (np.zeros(ndof) if self.u0_dot is None
                          else np.asarray(self.u0_dot, dtype=float).reshape(ndof))
        u_free = u0_global[free_dofs].copy()
        udot_free = u0_dot_global[free_dofs].copy()

        # F(0) reducido.
        F0_global = (np.zeros(ndof) if self.F_func is None
                      else np.asarray(self.F_func(0.0), dtype=float).reshape(ndof))
        F0_red = T.T @ F0_global

        # Aceleración inicial: M_red · ü₀ = F₀ − F_dir − C_red·u̇₀ − K_red·u₀.
        rhs0 = F0_red - F_dir - C_red @ udot_free - K_red @ u_free
        props_M = StiffnessProperties(
            is_symmetric=domain_is_symmetric(self.assembler.domain),
            is_positive_definite=True,
            size=n_free,
        )
        M_solver = select_solver(props_M, override=self.linear_algebra)
        uddot_free = M_solver.solve(M_red, rhs0)

        # Factorización reutilizable de A_eff = M + γΔt·C + βΔt²·K.
        dt = self.dt
        beta = self.beta
        gamma = self.gamma
        A_eff = (M_red + gamma * dt * C_red + beta * dt * dt * K_red).tocsr()
        props_A = StiffnessProperties(
            is_symmetric=domain_is_symmetric(self.assembler.domain),
            is_positive_definite=True,
            size=n_free,
        )
        A_solver = select_solver(props_A, override=self.linear_algebra)
        A_factor = A_solver.factorize(A_eff)

        # Historiales en globales. n_steps + 1 columnas (incluye t=0).
        n_steps = int(np.ceil(self.t_end / dt))
        t_history = np.linspace(0.0, n_steps * dt, n_steps + 1)
        u_history = np.zeros((ndof, n_steps + 1))
        udot_history = np.zeros((ndof, n_steps + 1))
        uddot_history = np.zeros((ndof, n_steps + 1))
        # Estado inicial en globales (con apoyos g).
        u_history[:, 0] = T @ u_free + g
        udot_history[:, 0] = T @ udot_free
        uddot_history[:, 0] = T @ uddot_free

        # Bucle temporal.
        half_dt2 = 0.5 * dt * dt
        one_minus_2beta = 1.0 - 2.0 * beta
        one_minus_gamma = 1.0 - gamma
        beta_dt2 = beta * dt * dt
        gamma_dt = gamma * dt

        for step in range(n_steps):
            t_next = t_history[step + 1]

            # Predictores.
            u_pred = u_free + dt * udot_free + half_dt2 * one_minus_2beta * uddot_free
            udot_pred = udot_free + dt * one_minus_gamma * uddot_free

            # F(t_{n+1}) reducido.
            F_global = (np.zeros(ndof) if self.F_func is None
                         else np.asarray(self.F_func(t_next), dtype=float).reshape(ndof))
            F_red = T.T @ F_global

            # Sistema efectivo: A_eff · ü_{n+1} = F − F_dir − C·u̇_pred − K·u_pred.
            rhs = F_red - F_dir - C_red @ udot_pred - K_red @ u_pred
            uddot_free = A_factor.solve(rhs)

            # Correctores.
            u_free = u_pred + beta_dt2 * uddot_free
            udot_free = udot_pred + gamma_dt * uddot_free

            # Volcado a historial global.
            u_history[:, step + 1] = T @ u_free + g
            udot_history[:, step + 1] = T @ udot_free
            uddot_history[:, step + 1] = T @ uddot_free

        _log.info(f"  -> {n_steps} pasos completados. "
                   f"Rayleigh: α={alpha_r:.4e}, β={beta_r:.4e}.")

        return TransientResult(
            t_history=t_history,
            u_history=u_history,
            udot_history=udot_history,
            uddot_history=uddot_history,
            n_steps=n_steps,
            alpha_rayleigh=alpha_r,
            beta_rayleigh=beta_r,
            converged=True,
        )
