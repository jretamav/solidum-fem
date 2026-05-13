# fenix_fem/fenix/math/solvers.py
from typing import Callable, Optional

import numpy as np
import scipy.sparse as sp

from fenix.constants import ZERO_TOL
from fenix.logging import get_logger
from fenix.math.convergence import (
    ConvergenceCriterion,
    make_convergence_from_config,
    stiffness_diag_scale,
)
from fenix.registry import SolverRegistry
from fenix.math.damping import rayleigh_from_modes
from fenix.math.linalg import EigenSolver, LUSolver, StiffnessProperties, select_solver
from fenix.results import ModalResult, TransientResult

try:
    from fenix.math.linalg.cholesky import CholeskyNotPositiveDefiniteError
except ImportError:
    class CholeskyNotPositiveDefiniteError(Exception):
        """Placeholder cuando scikit-sparse no está instalado (nunca se lanza)."""


_log = get_logger("solvers")


def _domain_is_symmetric(domain) -> bool:
    """Agrega los flags declarativos IS_SYMMETRIC / PRESERVES_SYMMETRY (ADR 0003 §2)."""
    for elem in domain.elements.values():
        if not getattr(type(elem), "PRESERVES_SYMMETRY", True):
            return False
        material = getattr(elem, "material", None)
        if material is not None and not getattr(type(material), "IS_SYMMETRIC", True):
            return False
    return True


@SolverRegistry.register
class LinearSolver:
    """Solucionador de sistemas algebraicos lineales en un solo paso."""
    def __init__(self, assembler, linear_algebra: str = "auto"):
        self.assembler = assembler
        self.linear_algebra = linear_algebra

    def solve(self, F_ext_global: np.ndarray) -> np.ndarray:
        _log.info("--- INICIANDO SOLVER LINEAL ---")
        self.assembler.assemble_system()
        K_global = self.assembler.K_global.copy()
        F = F_ext_global.copy()

        # Eliminación directa de DOFs prescritos (ADR 0004).
        K_red, F_red, T, g_full = self.assembler.reduce(K_global, F)

        props = StiffnessProperties(
            is_symmetric=_domain_is_symmetric(self.assembler.domain),
            is_positive_definite=True,
            size=K_red.shape[0],
        )
        linalg = select_solver(props, override=self.linear_algebra)
        try:
            u_red = linalg.solve(K_red, F_red)
        except CholeskyNotPositiveDefiniteError:
            # Fallback automático SPD→LU (ADR 0003 §5).
            _log.warning("Cholesky reportó no-positividad. Degradando a LU.")
            u_red = LUSolver().solve(K_red, F_red)

        U = self.assembler.expand(u_red, T, g_full)
        _log.info("  -> CONVERGENCIA ALCANZADA (1 Iteración).")
        return U

@SolverRegistry.register
class ModalSolver:
    """Análisis modal: frecuencias propias y modos de vibración (ADR 0009 fase 1).

    Resuelve el problema generalizado simétrico ``K·φ = ω²·M·φ`` sobre los
    DOFs libres tras imposición de Dirichlet por eliminación directa. La
    capa algebraica delega en :class:`EigenSolver` (ARPACK Lanczos con
    shift-invert centrado en ``sigma``).

    Parameters
    ----------
    assembler : Assembler
        Vinculación al modelo. Se reusan sus cachés de topología COO y
        ConstraintSet.
    n_modes : int
        Número de modos a calcular. Restricción ARPACK: ``1 ≤ n_modes < N_libre``.
    sigma : float, default 0.0
        Shift (rad²/s²) para shift-invert. ``0`` → frecuencias más bajas
        (uso por defecto en ingeniería estructural). Cualquier otro valor
        extrae los autovalores más cercanos a él.
    which : str, default "LM"
        Estrategia ARPACK. Con shift-invert el default es el adecuado;
        override solo para diagnóstico.
    tolerance : float, default 1.0e-9
        Tolerancia ARPACK sobre el residuo del autovalor.
    lumping : str, default "consistent"
        Discretización de la masa. Solo ``"consistent"`` implementada en
        fase 1 (ADR 0009 §1).
    linear_algebra : str, default "auto"
        Reservado. En fase 1 ARPACK usa su factorización LU interna sin
        override; el parámetro vive aquí para coherencia con los demás
        solvers y para no romper la firma cuando se cablee a la capa
        algebraica (ADR 0003).
    """

    def __init__(
        self,
        assembler,
        n_modes: int,
        *,
        sigma: float = 0.0,
        which: str = "LM",
        tolerance: float = 1.0e-9,
        lumping: str = "consistent",
        linear_algebra: str = "auto",
    ):
        self.assembler = assembler
        self.n_modes = int(n_modes)
        self.sigma = float(sigma)
        self.which = str(which)
        self.tolerance = float(tolerance)
        self.lumping = str(lumping)
        self.linear_algebra = str(linear_algebra)

    def solve(self) -> ModalResult:
        _log.info("--- INICIANDO SOLVER MODAL ---")

        self.assembler.assemble_system()
        K_global = self.assembler.K_global
        M_global = self.assembler.assemble_mass_matrix(lumping=self.lumping)

        K_red, M_red, T = self.assembler.reduce_pair(K_global, M_global)

        eig = EigenSolver(
            sigma=self.sigma, which=self.which, tol=self.tolerance
        )
        lambdas, phi_red = eig.solve(K_red, M_red, self.n_modes)

        # Clip de autovalores ligeramente negativos por ruido numérico
        # (modos de cuerpo rígido o casi). El cero exacto produce T = ∞,
        # que se reporta tal cual.
        lam_max = float(np.max(np.abs(lambdas))) if lambdas.size else 0.0
        lambdas = np.where(np.abs(lambdas) < 1e-12 * lam_max, 0.0, lambdas)
        lambdas = np.clip(lambdas, 0.0, None)

        omegas = np.sqrt(lambdas)
        frequencies_hz = omegas / (2.0 * np.pi)
        with np.errstate(divide="ignore", invalid="ignore"):
            periods = np.where(
                frequencies_hz > 0.0,
                np.where(frequencies_hz > 0.0, 1.0 / frequencies_hz, np.inf),
                np.inf,
            )

        modes = T @ phi_red

        if omegas.size:
            _log.info(
                f"  -> {self.n_modes} modo(s) calculado(s). "
                f"ω₁={omegas[0]:.4e} rad/s ({frequencies_hz[0]:.4e} Hz)"
            )

        return ModalResult(
            frequencies_rad=omegas,
            frequencies_hz=frequencies_hz,
            periods=periods,
            modes=np.asarray(modes),
            n_modes=self.n_modes,
            converged=True,
        )


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
            is_symmetric=_domain_is_symmetric(self.assembler.domain),
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
            is_symmetric=_domain_is_symmetric(self.assembler.domain),
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
class NonlinearSolver:
    """Solucionador incremental-iterativo de Newton-Raphson con paso adaptativo."""
    def __init__(self, assembler, convergence: Optional[ConvergenceCriterion] = None, max_iter=20, num_steps=10, adaptive=True, min_delta_lambda=1e-5, linear_algebra: str = "auto", freeze_tangent_after_iter: int | None = None):
        self.assembler = assembler
        # Política de convergencia (ADR 0007). El criterio se calibra una vez
        # al inicio de solve() con la escala del primer ensamblaje.
        self.convergence = convergence if convergence is not None else ConvergenceCriterion()
        self.max_iter = max_iter
        self.num_steps = num_steps
        self.adaptive = adaptive
        self.min_delta_lambda = min_delta_lambda
        self.linear_algebra = linear_algebra
        # Newton modificado (ADR 0003 fase 2): si es int, se factoriza fresca
        # las primeras N iteraciones del paso y se reusa la factorización en
        # las siguientes. None ⇒ Newton estándar (factoriza cada iteración).
        self.freeze_tangent_after_iter = freeze_tangent_after_iter
        # Se inicializa al comienzo de solve(); se degrada a LU si Cholesky aborta.
        self._linalg = None
        self._is_pd = True
        self._frozen_factor = None  # FactorizedSolver | None

    def _make_linalg(self, ndof: int):
        props = StiffnessProperties(
            is_symmetric=_domain_is_symmetric(self.assembler.domain),
            is_positive_definite=self._is_pd,
            size=ndof,
        )
        return select_solver(props, override=self.linear_algebra)

    def _solve_reduced(self, K_red, R_red, iteration: int = 0):
        """Resuelve K_red·δU_red = R_red con fallback SPD→LU y Newton modificado.

        Si ``freeze_tangent_after_iter`` está activo, factoriza fresca solo
        durante las primeras N iteraciones del paso y reusa la factorización
        cacheada en las siguientes (ADR 0003 §5 + fase 2).
        """
        if self.freeze_tangent_after_iter is None:
            try:
                return self._linalg.solve(K_red, R_red)
            except CholeskyNotPositiveDefiniteError:
                return self._fallback_to_lu_and_solve(K_red, R_red)

        threshold = self.freeze_tangent_after_iter
        if iteration < threshold or self._frozen_factor is None:
            try:
                self._frozen_factor = self._linalg.factorize(K_red)
            except CholeskyNotPositiveDefiniteError:
                self._is_pd = False
                self._linalg = self._make_linalg(K_red.shape[0])
                _log.warning("Cholesky reportó no-positividad. Degradando a LU para el resto del análisis.")
                self._frozen_factor = self._linalg.factorize(K_red)
        return self._frozen_factor.solve(R_red)

    def _fallback_to_lu_and_solve(self, K_red, R_red):
        _log.warning("Cholesky reportó no-positividad. Degradando a LU para el resto del análisis.")
        self._is_pd = False
        self._linalg = self._make_linalg(K_red.shape[0])
        return self._linalg.solve(K_red, R_red)

    def solve(self, F_ext_global: np.ndarray, step_callback=None) -> np.ndarray:
        domain = self.assembler.domain
        ndof = domain.total_dofs
        U_current = np.zeros(ndof)

        _log.info("--- INICIANDO SOLVER NO LINEAL (CONTROL DE PASO ADAPTATIVO) ---")

        # Tras reducción el solver algebraico opera sobre n_libre, no sobre ndof.
        cs = self.assembler.constraint_set
        n_free = ndof - len(cs)
        free_dofs = cs.free_dofs(ndof)
        self._linalg = self._make_linalg(n_free)

        load_factor = 0.0
        target_load = 1.0
        delta_lambda = 1.0 / self.num_steps
        step = 0

        while load_factor < target_load - 1e-9:
            step += 1

            if load_factor + delta_lambda > target_load:
                delta_lambda = target_load - load_factor

            next_load_factor = load_factor + delta_lambda
            _log.info(f"[PASO {step}] Intentando Factor de Carga: {next_load_factor:.4f} (Incremento: {delta_lambda:.4f})")

            F_ext_step = F_ext_global * next_load_factor
            U_iter = U_current.copy()
            converged = False

            for iteration in range(self.max_iter):
                K_global, F_int_global = self.assembler.assemble_non_linear_system(U_iter)
                R = F_ext_step - F_int_global

                K_red, R_red, T_op, g_inc = self.assembler.reduce(
                    K_global, R, U_current=U_iter, load_factor=next_load_factor
                )

                # Calibración del criterio en el primer ensamblaje de la corrida
                # (ADR 0007). Las escalas se derivan del estado inicial real:
                # ‖F_ext_global‖ para fuerza, ‖F_ext‖/max|diag(K)| para
                # desplazamiento. Si ‖F_ext_global‖ = 0 (control en desplazamiento
                # puro), usar la reacción interna del primer paso o fallback 1.0.
                if not self.convergence.is_calibrated:
                    force_scale = max(
                        np.linalg.norm(F_ext_global),
                        np.linalg.norm(F_int_global),
                        1.0,
                    )
                    K_diag = stiffness_diag_scale(K_global)
                    disp_scale = force_scale / K_diag
                    self.convergence.calibrate(force_scale, disp_scale)

                try:
                    delta_U_red = self._solve_reduced(K_red, R_red, iteration=iteration)
                except RuntimeError:
                    _log.error("Matriz Singular detectada.")
                    break

                delta_U = self.assembler.expand(delta_U_red, T_op, g_inc)
                U_iter += delta_U

                # Criterio dual fuerza + desplazamiento (ADR 0007). El residuo se
                # evalúa en DOFs libres (en prescritos R contiene la reacción
                # física, no un fallo de equilibrio).
                ref_force = max(np.linalg.norm(F_ext_step), np.linalg.norm(F_int_global))
                state = self.convergence.evaluate(
                    residual_norm=np.linalg.norm(R[free_dofs]),
                    ref_force=ref_force,
                    delta_u_norm=np.linalg.norm(delta_U),
                    u_norm=np.linalg.norm(U_iter),
                )
                _log.info(
                    f"  Iteración {iteration+1:2d} | "
                    f"R/tol_F: {state.ratio_force:.4e} | "
                    f"dU/tol_d: {state.ratio_disp:.4e}"
                )

                if state.converged:
                    _log.info("  -> CONVERGENCIA ALCANZADA.")
                    self.assembler.commit_all_states()

                    U_current = U_iter
                    load_factor = next_load_factor
                    converged = True

                    if self.adaptive and iteration < 4 and delta_lambda < (1.0 / self.num_steps):
                        delta_lambda = min(delta_lambda * 1.5, 1.0 / self.num_steps)
                        _log.info(f"  -> Acelerando el próximo incremento a {delta_lambda:.4f}")

                    if step_callback:
                        step_callback(step, U_current, load_factor)

                    break

            # Cierre de paso (converja o no): la K_t cambia con U_current,
            # así que la factorización congelada deja de ser válida para el
            # siguiente paso (Newton modificado, ADR 0003 fase 2).
            self._frozen_factor = None

            if not converged:
                if self.adaptive:
                    delta_lambda /= 2.0
                    _log.warning(f"NO CONVERGIÓ. Bisección: reduciendo incremento a {delta_lambda:.4f}")
                    if delta_lambda < self.min_delta_lambda:
                        raise RuntimeError(f"El incremento de carga ({delta_lambda:.2e}) cayó por debajo del mínimo ({self.min_delta_lambda:.2e}). El solver ha divergido.")
                else:
                    raise RuntimeError(f"El solucionador no convergió en el paso {step}.")

        return U_current

@SolverRegistry.register
class ArcLengthSolver:
    """
    Solucionador no lineal con Método de Longitud de Arco Cilíndrico (Crisfield).
    Permite trazar curvas de equilibrio con fenómenos de snap-through y snap-back
    variando simultáneamente los desplazamientos y la carga externa.
    """
    def __init__(self, assembler, convergence: Optional[ConvergenceCriterion] = None, max_iter=20, max_lambda=1.0, initial_dl=0.1, max_steps=100,
                 dl_grow_factor=1.5, dl_max_factor=5.0, dl_shrink_factor=0.6,
                 dl_grow_iter_threshold=4, dl_shrink_iter_threshold=8,
                 linear_algebra: str = "auto"):
        self.assembler = assembler
        # Política de convergencia (ADR 0007). Compartida con NonlinearSolver:
        # cambiar la política aquí llega automáticamente al arc-length.
        self.convergence = convergence if convergence is not None else ConvergenceCriterion()
        self.max_iter = max_iter
        self.max_lambda = max_lambda
        self.dl = initial_dl
        self.max_steps = max_steps
        # Factores de auto-ajuste de la longitud de arco:
        #   Si converge en < dl_grow_iter_threshold iter → ampliar dl × dl_grow_factor (max: initial_dl × dl_max_factor)
        #   Si converge en > dl_shrink_iter_threshold iter → reducir dl × dl_shrink_factor
        self.dl_grow_factor = dl_grow_factor
        self.dl_max_factor = dl_max_factor
        self.dl_shrink_factor = dl_shrink_factor
        self.dl_grow_iter_threshold = dl_grow_iter_threshold
        self.dl_shrink_iter_threshold = dl_shrink_iter_threshold
        self.linear_algebra = linear_algebra
        self._linalg = None

    def _make_linalg(self, ndof: int):
        # Régimen postcrítico: K_t puede ser indefinida → no asumir SPD.
        props = StiffnessProperties(
            is_symmetric=_domain_is_symmetric(self.assembler.domain),
            is_positive_definite=False,
            size=ndof,
        )
        return select_solver(props, override=self.linear_algebra)

    def _solve(self, K, b):
        """Resuelve K·x = b con fallback Cholesky→LU.

        El default ``auto`` ya elige LU para régimen postcrítico, pero si el
        usuario fuerza ``linear_algebra: cholesky`` desde YAML para
        diagnóstico y ``K_t`` se vuelve indefinida en algún paso, degradamos
        limpiamente a LU para que el override no rompa el análisis.
        """
        try:
            return self._linalg.solve(K, b)
        except CholeskyNotPositiveDefiniteError:
            _log.warning("Cholesky reportó no-positividad. Degradando a LU para el resto del análisis.")
            self._linalg = LUSolver()
            return self._linalg.solve(K, b)

    def _negative_pivots(self, K) -> int | None:
        """Diagnóstico de bifurcación vía Sturm sequence (ADR 0003 fase 2).

        Cuando ``LDLTSolver`` esté implementado (backend con LDLᵀ verdadero
        de Bunch-Kaufman, p. ej. pypardiso), este método factoriza ``K_t`` y
        retorna ``factor.n_negative_pivots``. Por la ley de inercia de
        Sylvester equivale al número de autovalores negativos: 0 antes de
        bifurcación, 1 tras el primer punto crítico, etc.

        En fase 2 retorna ``None`` (LDLᵀ es placeholder; LU general no
        preserva inercia y dar un conteo aproximado sería engañoso).
        """
        return None

    def solve(self, F_ext_ref: np.ndarray, step_callback=None) -> np.ndarray:
        domain = self.assembler.domain
        ndof = domain.total_dofs
        U_current = np.zeros(ndof)

        lambda_curr = 0.0
        step = 0
        dl = self.dl

        delta_U_step = np.zeros(ndof)  # Historial del incremento del paso para guiar el arco

        _log.info("--- INICIANDO SOLVER NO LINEAL (MÉTODO ARC-LENGTH) ---")

        cs = self.assembler.constraint_set
        n_free = ndof - len(cs)
        self._linalg = self._make_linalg(n_free)

        while lambda_curr < self.max_lambda and step < self.max_steps:
            step += 1
            _log.info(f"[PASO {step}] Longitud de Arco (dl): {dl:.4e}")

            U_iter = U_current.copy()
            lambda_iter = lambda_curr
            converged = False

            # --- 1. PREDICTOR ---
            K_global, F_int_global = self.assembler.assemble_non_linear_system(U_iter)

            # Calibración del criterio en el primer ensamblaje (ADR 0007).
            # La carga de referencia F_ext_ref es la escala natural en arc-length
            # (el factor de carga lambda se irá ajustando, pero F_ext_ref es fijo).
            if not self.convergence.is_calibrated:
                force_scale = max(
                    np.linalg.norm(F_ext_ref),
                    np.linalg.norm(F_int_global),
                    1.0,
                )
                K_diag = stiffness_diag_scale(K_global)
                disp_scale = force_scale / K_diag
                self.convergence.calibrate(force_scale, disp_scale)

            K_t_red, F_t_red, T_t, g_t = self.assembler.reduce(K_global, F_ext_ref.copy())

            try:
                du_t_red = self._solve(K_t_red, F_t_red)
            except RuntimeError:
                _log.error("Matriz singular en predictor. Bisección de dl...")
                dl /= 2.0
                continue

            du_t = self.assembler.expand(du_t_red, T_t, g_t)

            # Determinar el sentido del avance (evitar regresar por donde vinimos)
            sign = 1.0
            if step > 1 and np.dot(delta_U_step, du_t) < 0:
                sign = -1.0

            dlambda = sign * dl / (np.linalg.norm(du_t) + ZERO_TOL)

            # Si el paso predictor sobrepasaría max_lambda, fijar lambda exactamente
            final_step = (sign > 0 and lambda_curr + dlambda >= self.max_lambda - ZERO_TOL)
            if final_step:
                dlambda = self.max_lambda - lambda_curr

            lambda_iter += dlambda
            dU_iter = dlambda * du_t
            U_iter += dU_iter

            # --- 2. CORRECTOR ITERATIVO ---
            for iteration in range(self.max_iter):
                K_global, F_int_global = self.assembler.assemble_non_linear_system(U_iter)
                R = lambda_iter * F_ext_ref - F_int_global

                K_t_red, F_t_red, T_t, g_t = self.assembler.reduce(K_global, F_ext_ref.copy())
                K_red, R_red, T_R, g_R = self.assembler.reduce(
                    K_global, R, U_current=U_iter, load_factor=lambda_iter
                )

                try:
                    du_R_red = self._solve(K_red, R_red)
                    du_t_red = self._solve(K_t_red, F_t_red)
                except RuntimeError:
                    _log.error("Matriz Singular en corrector.")
                    break

                du_R = self.assembler.expand(du_R_red, T_R, g_R)
                du_t = self.assembler.expand(du_t_red, T_t, g_t)

                if final_step:
                    # Último paso: lambda fijo, solo corrección de desplazamientos (Newton-Raphson puro)
                    ddlambda = 0.0
                    dU_update = du_R
                else:
                    # Ecuación cuadrática de restricción de Crisfield
                    dU_new = dU_iter + du_R
                    a = np.dot(du_t, du_t)
                    b = 2.0 * np.dot(dU_new, du_t)
                    c = np.dot(dU_new, dU_new) - dl**2

                    det = b**2 - 4.0 * a * c
                    if det < 0:
                        _log.error("Raíces imaginarias. La solución diverge del arco.")
                        break

                    ddl1 = (-b + np.sqrt(det)) / (2.0 * a)
                    ddl2 = (-b - np.sqrt(det)) / (2.0 * a)

                    # Elegir la raíz que produzca el menor ángulo con el incremento previo
                    theta1 = np.dot(dU_iter, dU_new + ddl1 * du_t)
                    theta2 = np.dot(dU_iter, dU_new + ddl2 * du_t)
                    ddlambda = ddl1 if theta1 > theta2 else ddl2
                    dU_update = du_R + ddlambda * du_t
                    dU_iter = dU_new + ddlambda * du_t

                # Actualizar iteraciones
                lambda_iter += ddlambda
                if not final_step:
                    dU_iter = dU_iter  # ya actualizado arriba
                else:
                    dU_iter = dU_iter + dU_update
                U_iter = U_current + dU_iter

                ref_force = max(
                    np.linalg.norm(F_ext_ref) * abs(lambda_iter),
                    np.linalg.norm(F_int_global),
                )
                state = self.convergence.evaluate(
                    residual_norm=np.linalg.norm(R[cs.free_dofs(ndof)]),
                    ref_force=ref_force,
                    delta_u_norm=np.linalg.norm(dU_update),
                    u_norm=np.linalg.norm(U_iter),
                )
                _log.info(
                    f"  Iter. {iteration+1:2d} | lam={lambda_iter:.4f} | "
                    f"R/tol_F: {state.ratio_force:.4e} | "
                    f"dU/tol_d: {state.ratio_disp:.4e}"
                )

                if state.converged:
                    _log.info(f"  -> CONVERGENCIA. (Lambda alcanzado: {lambda_iter:.4f})")
                    self.assembler.commit_all_states()

                    U_current = U_iter; lambda_curr = lambda_iter; delta_U_step = dU_iter
                    converged = True
                    # Auto-ajuste de longitud de arco
                    if iteration < self.dl_grow_iter_threshold:
                        dl = min(dl * self.dl_grow_factor, self.dl * self.dl_max_factor)
                    elif iteration > self.dl_shrink_iter_threshold:
                        dl *= self.dl_shrink_factor

                    if step_callback:
                        step_callback(step, U_current, lambda_curr)

                    break

            if not converged:
                dl *= 0.5
                _log.warning(f"Bisección: reduciendo longitud de arco a {dl:.4e}")
                if dl < 1e-6 * self.dl:
                    raise RuntimeError("Arc-Length fracasó irreparablemente.")

        return U_current
