# fenix_fem/fenix/math/solvers/newmark.py
"""``NewmarkSolver`` — integración Newmark-β para ``M·ü + C·u̇ + K·u = F(t)``
(ADR 0009 fase 3).

``NewtonNewmarkSolver`` (ADR 0009 fase 4) — variante no lineal con Newton-Raphson
dentro de cada paso temporal. Subclase de ``NewmarkSolver`` que reusa
predictores/correctores y reducción de Dirichlet, sobrescribe ``solve()`` para
ensamblar el residuo dinámico no lineal y resolver iterativamente.
"""
from typing import Callable, Optional

import numpy as np
import scipy.sparse as sp

from fenix.math.convergence import ConvergenceCriterion, stiffness_diag_scale
from fenix.math.damping import rayleigh_from_modes
from fenix.math.linalg import StiffnessProperties, select_solver
from fenix.math.solvers._shared import (
    CholeskyNotPositiveDefiniteError,
    _log,
    domain_is_symmetric,
)
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


@SolverRegistry.register
class NewtonNewmarkSolver(NewmarkSolver):
    """Análisis dinámico transitorio **no lineal** por Newmark + Newton (ADR 0009 fase 4).

    Subclase de :class:`NewmarkSolver` para problemas con materiales con historia
    (plasticidad, daño) o no linealidad geométrica. En cada paso temporal
    introduce un bucle de Newton-Raphson sobre el residuo dinámico:

    .. math::

        \\mathbf R(\\ddot{\\mathbf u}_{n+1}) = \\mathbf F_\\text{ext}(t_{n+1}) - \\mathbf F_\\text{int}(\\mathbf u_{n+1})
                                              - \\mathbf C\\,\\dot{\\mathbf u}_{n+1} - \\mathbf M\\,\\ddot{\\mathbf u}_{n+1}

    El jacobiano dinámico
    ``J = M + γΔt·C + βΔt²·K_tangente`` se re-factoriza cada iteración (o
    cada N iter con ``freeze_tangent_after_iter``, Newton modificado ADR 0003).

    El amortiguamiento Rayleigh se calibra con la rigidez **elástica de
    referencia** ``K_0`` (al inicio del análisis, ``u = 0``) y se mantiene
    constante en el tiempo. Convención estándar (Abaqus, ANSYS, OpenSees);
    evita acoplamiento ad-hoc entre disipación viscosa y plástica.

    Parameters extra a :class:`NewmarkSolver`
    -----------------------------------------
    convergence : ConvergenceCriterion, optional
        Política de convergencia (ADR 0007). Default: ``ConvergenceCriterion()``
        con tolerancias por defecto del proyecto.
    max_iter : int, default 20
        Máximo iteraciones Newton por paso temporal. Si se agota se lanza
        ``RuntimeError``; el estado trial no se commitea.
    freeze_tangent_after_iter : int or None, default None
        Si ``int``, Newton modificado (ADR 0003 fase 2): factoriza fresco las
        primeras N iter de cada paso y reusa la factorización en las siguientes.
        ``None`` ⇒ Newton estándar (re-factoriza cada iteración).

    Notes
    -----
    Si los materiales son lineales, el solver converge en una iteración por
    paso y reproduce exactamente :class:`NewmarkSolver` (validado en
    ``tests/test_newmark_nonlinear.py``).

    Ver ``docs/specs/NewtonNewmarkSolver.md``.
    """

    def __init__(
        self,
        assembler,
        t_end: float,
        dt: float,
        *,
        convergence: Optional[ConvergenceCriterion] = None,
        max_iter: int = 20,
        freeze_tangent_after_iter: Optional[int] = None,
        beta: float = 0.25,
        gamma: float = 0.5,
        rayleigh: Optional[dict] = None,
        u0: Optional[np.ndarray] = None,
        u0_dot: Optional[np.ndarray] = None,
        F_func: Optional[Callable[[float], np.ndarray]] = None,
        linear_algebra: str = "auto",
        lumping: str = "consistent",
    ):
        super().__init__(
            assembler, t_end, dt,
            beta=beta, gamma=gamma, rayleigh=rayleigh,
            u0=u0, u0_dot=u0_dot, F_func=F_func,
            linear_algebra=linear_algebra, lumping=lumping,
        )
        self.convergence = convergence if convergence is not None else ConvergenceCriterion()
        self.max_iter = int(max_iter)
        self.freeze_tangent_after_iter = freeze_tangent_after_iter

    def solve(self) -> TransientResult:
        _log.info("--- INICIANDO SOLVER NEWMARK NO LINEAL (Newton dentro de Newmark) ---")

        # Disparar el build de topología del assembler antes de leer ndof.
        # `assembler.ndof` solo se rellena tras un primer ensamblaje.
        if self.assembler.domain.total_dofs == 0:
            self.assembler.domain.generate_equation_numbers()
        ndof = self.assembler.domain.total_dofs

        # Condiciones iniciales en globales. Aseguramos compatibilidad con apoyos.
        u_total = (np.zeros(ndof) if self.u0 is None
                    else np.asarray(self.u0, dtype=float).reshape(ndof)).copy()
        udot_total = (np.zeros(ndof) if self.u0_dot is None
                       else np.asarray(self.u0_dot, dtype=float).reshape(ndof)).copy()

        # Masa (constante).
        M = self.assembler.assemble_mass_matrix(lumping=self.lumping)

        # Sistema no lineal en u_0 para obtener K_0 (Rayleigh) y F_int(u_0).
        K_0, F_int = self.assembler.assemble_non_linear_system(u_total)

        alpha_r, beta_r = self._resolve_rayleigh(self.rayleigh_cfg)
        C = alpha_r * M + beta_r * K_0

        # Reducción por Dirichlet (apoyos constantes).
        cs = self.assembler.constraint_set
        T_op, g_vec = cs.build(ndof)
        free_dofs = cs.free_dofs(ndof)
        n_free = T_op.shape[1]

        # Forzar u_total compatible con apoyos: u_total[prescribed] = g_vec[prescribed].
        # En DOFs libres se mantiene self.u0; en prescritos se impone g_vec.
        prescribed = np.ones(ndof, dtype=bool)
        prescribed[free_dofs] = False
        u_total[prescribed] = g_vec[prescribed]
        # u̇ y ü en DOFs prescritos = 0 (apoyos constantes en el tiempo).
        udot_total[prescribed] = 0.0

        # Reducción de las matrices invariantes en el tiempo.
        M_red = (T_op.T @ M @ T_op).tocsr()
        if alpha_r != 0.0 or beta_r != 0.0:
            C_red = (T_op.T @ C @ T_op).tocsr()
        else:
            C_red = sp.csr_matrix(M_red.shape)

        # Aceleración inicial consistente: M·ü₀ = F_ext(0) - F_int(u₀) - C·u̇₀.
        F0_global = (np.zeros(ndof) if self.F_func is None
                      else np.asarray(self.F_func(0.0), dtype=float).reshape(ndof))
        rhs0 = T_op.T @ (F0_global - F_int - C @ udot_total)
        props_M = StiffnessProperties(
            is_symmetric=domain_is_symmetric(self.assembler.domain),
            is_positive_definite=True,
            size=n_free,
        )
        M_solver = select_solver(props_M, override=self.linear_algebra)
        uddot_free = M_solver.solve(M_red, rhs0)
        uddot_total = np.zeros(ndof)
        uddot_total[free_dofs] = uddot_free

        # Historiales.
        dt = self.dt
        beta = self.beta
        gamma = self.gamma
        n_steps = int(np.ceil(self.t_end / dt))
        t_history = np.linspace(0.0, n_steps * dt, n_steps + 1)
        u_history = np.zeros((ndof, n_steps + 1))
        udot_history = np.zeros((ndof, n_steps + 1))
        uddot_history = np.zeros((ndof, n_steps + 1))
        u_history[:, 0] = u_total
        udot_history[:, 0] = udot_total
        uddot_history[:, 0] = uddot_total

        # Coeficientes Newmark precomputados.
        half_dt2 = 0.5 * dt * dt
        one_minus_2beta = 1.0 - 2.0 * beta
        one_minus_gamma = 1.0 - gamma
        beta_dt2 = beta * dt * dt
        gamma_dt = gamma * dt

        # Selector del solver lineal (puede degradar Cholesky→LU si emerge no-PD).
        is_sym = domain_is_symmetric(self.assembler.domain)
        is_pd = True
        props_J = StiffnessProperties(
            is_symmetric=is_sym, is_positive_definite=is_pd, size=n_free,
        )
        linalg = select_solver(props_J, override=self.linear_algebra)
        frozen_factor = None  # Newton modificado: cache de factorización por paso

        for step in range(n_steps):
            t_next = t_history[step + 1]

            # Predictores Newmark.
            u_pred = u_total + dt * udot_total + half_dt2 * one_minus_2beta * uddot_total
            udot_pred = udot_total + dt * one_minus_gamma * uddot_total

            # Inicialización del Newton: ü^(0) = ü_n (continuidad).
            uddot_iter = uddot_total.copy()
            u_iter = u_pred + beta_dt2 * uddot_iter
            udot_iter = udot_pred + gamma_dt * uddot_iter

            F_ext_next = (np.zeros(ndof) if self.F_func is None
                           else np.asarray(self.F_func(t_next), dtype=float).reshape(ndof))

            converged = False
            for it in range(self.max_iter):
                K_t, F_int_iter = self.assembler.assemble_non_linear_system(u_iter)
                R = F_ext_next - F_int_iter - C @ udot_iter - M @ uddot_iter
                R_red = T_op.T @ R

                J = (M_red + gamma_dt * C_red + beta_dt2 * (T_op.T @ K_t @ T_op)).tocsr()

                # Calibración del criterio en el primer ensamblaje.
                if not self.convergence.is_calibrated:
                    force_scale = max(
                        np.linalg.norm(F_ext_next),
                        np.linalg.norm(F_int_iter),
                        1.0,
                    )
                    K_diag = stiffness_diag_scale(K_t)
                    disp_scale = force_scale / K_diag
                    self.convergence.calibrate(force_scale, disp_scale)

                # Resolver δü_red, con Newton modificado opcional.
                threshold = self.freeze_tangent_after_iter
                try:
                    if threshold is None:
                        delta_uddot_red = linalg.solve(J, R_red)
                    elif it < threshold or frozen_factor is None:
                        frozen_factor = linalg.factorize(J)
                        delta_uddot_red = frozen_factor.solve(R_red)
                    else:
                        delta_uddot_red = frozen_factor.solve(R_red)
                except CholeskyNotPositiveDefiniteError:
                    _log.warning("Cholesky no-PD en NewtonNewmark; degradando a LU.")
                    is_pd = False
                    props_J = StiffnessProperties(
                        is_symmetric=is_sym, is_positive_definite=False, size=n_free,
                    )
                    linalg = select_solver(props_J, override=self.linear_algebra)
                    frozen_factor = None
                    delta_uddot_red = linalg.solve(J, R_red)

                # Actualizar incógnitas con correctores Newmark.
                uddot_iter[free_dofs] += delta_uddot_red
                u_iter = u_pred + beta_dt2 * uddot_iter
                udot_iter = udot_pred + gamma_dt * uddot_iter

                # δu corresponde a βΔt² · δü (cambio en desplazamiento por la iter).
                delta_u_norm = beta_dt2 * np.linalg.norm(delta_uddot_red)
                ref_force = max(np.linalg.norm(F_ext_next), np.linalg.norm(F_int_iter))
                state = self.convergence.evaluate(
                    residual_norm=np.linalg.norm(R[free_dofs]),
                    ref_force=ref_force,
                    delta_u_norm=delta_u_norm,
                    u_norm=np.linalg.norm(u_iter[free_dofs]),
                )

                if state.converged:
                    _log.info(
                        f"  [PASO {step+1}/{n_steps}] t={t_next:.4e} | "
                        f"iter={it+1} | R/tol_F={state.ratio_force:.2e}"
                    )
                    self.assembler.commit_all_states()
                    converged = True
                    break

            # Cierre de paso: la factorización congelada solo es válida dentro
            # del paso (K_t cambia al pasar a t_{n+2}).
            frozen_factor = None

            if not converged:
                raise RuntimeError(
                    f"NewtonNewmarkSolver: paso {step+1} (t={t_next:.4e}) no "
                    f"convergió en {self.max_iter} iteraciones."
                )

            u_total = u_iter
            udot_total = udot_iter
            uddot_total = uddot_iter

            u_history[:, step + 1] = u_total
            udot_history[:, step + 1] = udot_total
            uddot_history[:, step + 1] = uddot_total

        _log.info(f"  -> {n_steps} pasos completados (no lineal). "
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
