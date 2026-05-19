# fenix_fem/fenix/math/solvers/newmark.py
"""``NewmarkSolver`` — integración Newmark-β para ``M·ü + C·u̇ + K·u = F(t)``
(ADR 0009 fase 3).

``NewtonNewmarkSolver`` (ADR 0009 fase 4) — variante no lineal con Newton-Raphson
dentro de cada paso temporal. Subclase de ``NewmarkSolver`` que reusa
predictores/correctores y reducción de Dirichlet, sobrescribe ``solve()`` para
ensamblar el residuo dinámico no lineal y resolver iterativamente.
"""
from __future__ import annotations

from typing import Callable

import numpy as np
import scipy.sparse as sp

from fenix.constants import (
    LINE_SEARCH_MAX_BACKTRACKS,
    LINE_SEARCH_RHO,
)
from fenix.math.convergence import ConvergenceCriterion, stiffness_diag_scale
from fenix.math.damping import resolve_rayleigh_config
from fenix.math.linalg import StiffnessProperties, select_solver
from fenix.math.solvers._shared import (
    CholeskyNotPositiveDefiniteError,
    _log,
    domain_is_symmetric,
)
from fenix.math.solvers.diagnostics import (
    classify_divergence,
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
        Discretización de masa. ``"consistent"`` (default) o ``"lumped"``
        (ADR 0009 fase 2).
    """

    PIPELINE_KIND = "transient"

    def __init__(
        self,
        assembler,
        t_end: float,
        dt: float,
        *,
        beta: float = 0.25,
        gamma: float = 0.5,
        rayleigh: dict | None = None,
        u0: np.ndarray | None = None,
        u0_dot: np.ndarray | None = None,
        F_func: Callable[[float], np.ndarray] | None = None,
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


    def solve(self) -> TransientResult:
        _log.info("--- INICIANDO SOLVER NEWMARK ---")

        # Ensamblar K (lineal en u=0) y M (consistente, cacheada).
        self.assembler.assemble_system()
        K = self.assembler.K_global
        M = self.assembler.assemble_mass_matrix(lumping=self.lumping)

        # Amortiguamiento Rayleigh: C = α·M + β·K.
        alpha_r, beta_r = resolve_rayleigh_config(
            self.rayleigh_cfg, source=type(self).__name__,
        )
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
        convergence: ConvergenceCriterion | None = None,
        max_iter: int = 20,
        freeze_tangent_after_iter: int | None = None,
        line_search: bool = False,
        beta: float = 0.25,
        gamma: float = 0.5,
        rayleigh: dict | None = None,
        u0: np.ndarray | None = None,
        u0_dot: np.ndarray | None = None,
        F_func: Callable[[float], np.ndarray] | None = None,
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
        # Line search por descenso no monótono (ADR 0011). Default ``False``
        # por la misma razón que ``NonlinearSolver``: contraproducente en
        # regímenes con tangente consistente cuasi-cuadrática. La masa del
        # jacobiano dinámico además estabiliza el Newton sin necesidad de
        # globalización en la mayoría de los casos (ver auditoría fase A).
        self.line_search = bool(line_search)

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

        alpha_r, beta_r = resolve_rayleigh_config(
            self.rayleigh_cfg, source=type(self).__name__,
        )
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

            # Historial del paso para clasificación de divergencia (ADR 0011).
            residual_history: list[float] = []
            delta_history: list[float] = []
            singular_tangent_seen = False
            last_residual = float("inf")
            last_delta = 0.0
            last_alpha = 1.0

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

                # Resolver δü_red, con Newton modificado opcional. Si la
                # tangente dinámica J = M + γΔt·C + βΔt²·K_t resulta singular
                # (RuntimeError de `splu` tras degradar a LU — cerca de
                # bifurcaciones dinámicas o tangentes patológicas), se flipea
                # el flag para que `classify_divergence` devuelva
                # `SingularTangentError` (ADR 0011).
                threshold = self.freeze_tangent_after_iter
                try:
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
                except RuntimeError:
                    _log.error("Tangente dinámica singular en NewtonNewmark.")
                    singular_tangent_seen = True
                    break

                # Line search por descenso no monótono (ADR 0011). Escala δü
                # (y por ende δu, δu̇ vía correctores Newmark) por α ∈ (0, 1].
                R_norm_before = float(np.linalg.norm(R[free_dofs]))
                alpha = self._armijo_step_dynamic(
                    uddot_iter, delta_uddot_red, u_pred, udot_pred,
                    beta_dt2, gamma_dt, free_dofs,
                    F_ext_next, C, M, R_norm_before,
                )
                last_alpha = alpha

                # Actualizar incógnitas con correctores Newmark (α·δü).
                uddot_iter[free_dofs] += alpha * delta_uddot_red
                u_iter = u_pred + beta_dt2 * uddot_iter
                udot_iter = udot_pred + gamma_dt * uddot_iter

                # Re-ensamblar para obtener R y F_int coherentes con el U avanzado.
                # (El last_alpha != 1 en su ensamblaje interno habrá dejado el
                # state trial coherente, pero re-ensamblar aquí mantiene
                # simetría con la rama sin line search.)
                _, F_int_after = self.assembler.assemble_non_linear_system(u_iter)
                R_after = F_ext_next - F_int_after - C @ udot_iter - M @ uddot_iter
                R_norm_after = float(np.linalg.norm(R_after[free_dofs]))

                # δu corresponde a βΔt² · α · δü (cambio en desplazamiento por la iter).
                delta_u_norm = beta_dt2 * float(np.linalg.norm(alpha * delta_uddot_red))
                ref_force = max(np.linalg.norm(F_ext_next), np.linalg.norm(F_int_after))
                state = self.convergence.evaluate(
                    residual_norm=R_norm_after,
                    ref_force=ref_force,
                    delta_u_norm=delta_u_norm,
                    u_norm=np.linalg.norm(u_iter[free_dofs]),
                )
                residual_history.append(R_norm_after)
                delta_history.append(delta_u_norm)
                last_residual = R_norm_after
                last_delta = delta_u_norm

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
                # Clasificar el modo y lanzar excepción tipada (ADR 0011).
                err_cls = classify_divergence(
                    residual_history, delta_history,
                    singular_tangent_detected=singular_tangent_seen,
                )
                raise err_cls(
                    last_residual=last_residual,
                    last_delta=last_delta,
                    last_load_factor=t_next,
                    n_bisections=0,
                    extra_message=(
                        f"NewtonNewmarkSolver: paso {step+1} (t={t_next:.4e}) no "
                        f"convergió en {self.max_iter} iteraciones; "
                        f"line_search={'on' if self.line_search else 'off'}; "
                        f"último α={last_alpha:.3f}."
                    ),
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

    def _armijo_step_dynamic(self, uddot_iter: np.ndarray,
                              delta_uddot_red: np.ndarray,
                              u_pred: np.ndarray, udot_pred: np.ndarray,
                              beta_dt2: float, gamma_dt: float,
                              free_dofs,
                              F_ext_next: np.ndarray,
                              C, M,
                              R_norm_before: float) -> float:
        """Line search por descenso no monótono para el residuo dinámico (ADR 0011).

        Aplica el mismo patrón que ``NonlinearSolver._armijo_step`` al
        residuo dinámico ``R = F_ext − F_int(u) − C·u̇ − M·ü``. Escala
        ``δü`` por ``α ∈ (0, 1]`` y propaga la consistencia a ``δu``,
        ``δu̇`` vía los correctores Newmark.

        Devuelve solo ``α``: el bucle exterior re-ensambla para obtener
        F_int y R coherentes con el U avanzado.

        Cuando ``self.line_search=False`` devuelve 1.0 directamente sin
        evaluar.
        """
        if not self.line_search:
            return 1.0

        rho = LINE_SEARCH_RHO
        max_bt = LINE_SEARCH_MAX_BACKTRACKS

        alpha = 1.0
        for _ in range(max_bt + 1):
            uddot_trial = uddot_iter.copy()
            uddot_trial[free_dofs] += alpha * delta_uddot_red
            u_trial = u_pred + beta_dt2 * uddot_trial
            udot_trial = udot_pred + gamma_dt * uddot_trial

            _, F_int_trial = self.assembler.assemble_non_linear_system(u_trial)
            R_trial = F_ext_next - F_int_trial - C @ udot_trial - M @ uddot_trial
            R_trial_norm = float(np.linalg.norm(R_trial[free_dofs]))

            if R_trial_norm <= R_norm_before:
                return alpha
            alpha *= rho

        return alpha


def _hht_autoderive_beta_gamma(alpha: float) -> tuple[float, float]:
    """Auto-derivación de β, γ desde α para HHT-α (Hilber 1977).

    Preserva orden 2 de precisión y estabilidad incondicional:
        β = (1 − α)² / 4
        γ = (1 − 2α) / 2
    Con α=0 recupera Newmark trapezoidal (β=1/4, γ=1/2).
    """
    beta = (1.0 - alpha) ** 2 / 4.0
    gamma = (1.0 - 2.0 * alpha) / 2.0
    return beta, gamma


@SolverRegistry.register
class HHTSolver(NewmarkSolver):
    """Análisis dinámico transitorio lineal por Hilber-Hughes-Taylor α (1977).

    Variante de :class:`NewmarkSolver` con **disipación numérica controlada**
    en altas frecuencias. La ecuación de equilibrio temporal evalúa fuerzas
    elásticas, viscosas y externas en un instante intermedio ``t_{n+1−α}``,
    mientras la fuerza inercial se mantiene en ``t_{n+1}``:

    .. math::

        \\mathbf M\\,\\ddot{\\mathbf u}_{n+1}
          + (1+\\alpha)\\mathbf C\\dot{\\mathbf u}_{n+1} - \\alpha\\mathbf C\\dot{\\mathbf u}_n
          + (1+\\alpha)\\mathbf K\\mathbf u_{n+1} - \\alpha\\mathbf K\\mathbf u_n
          = (1+\\alpha)\\mathbf F_{n+1} - \\alpha\\mathbf F_n

    Para preservar segundo orden y estabilidad incondicional, ``β`` y ``γ``
    se auto-derivan desde ``α`` (Hilber 1977): ``β=(1−α)²/4, γ=(1−2α)/2``.
    Con ``α=0`` recupera Newmark trapezoidal exactamente.

    Parameters
    ----------
    alpha : float, default −0.05
        Parámetro HHT en [−1/3, 0]. Valor canónico ``−0.05`` (Hilber 1977,
        Abaqus, OpenSees). ``α=0`` → sin disipación numérica. ``α=−1/3`` →
        disipación máxima manteniendo orden 2. Radio espectral en alta
        frecuencia: ``ρ_∞ = (1+α)/(1−α)``.
    beta, gamma : float or None
        Si ``None`` (default), se autoderivan desde ``alpha``. Override
        explícito posible pero no recomendado: combinaciones arbitrarias
        pueden perder estabilidad incondicional u orden 2.
    (resto heredado de NewmarkSolver — ver docstring del padre)

    Notes
    -----
    La matriz efectiva ``A_eff = M + (1+α)γΔt·C + (1+α)βΔt²·K`` es constante
    (depende solo de M, C, K y de los coeficientes) y se factoriza una sola
    vez al inicio. Cada paso temporal es una resolución triangular barata,
    igual que `NewmarkSolver`.

    Ver ``docs/specs/HHTSolver.md``.
    """

    def __init__(
        self,
        assembler,
        t_end: float,
        dt: float,
        *,
        alpha: float = -0.05,
        beta: float | None = None,
        gamma: float | None = None,
        rayleigh: dict | None = None,
        u0: np.ndarray | None = None,
        u0_dot: np.ndarray | None = None,
        F_func: Callable[[float], np.ndarray] | None = None,
        linear_algebra: str = "auto",
        lumping: str = "consistent",
    ):
        if not (-1.0 / 3.0 - 1.0e-12 <= alpha <= 0.0 + 1.0e-12):
            raise ValueError(
                f"HHTSolver: alpha={alpha} fuera del rango [-1/3, 0]. "
                "Valores fuera de este rango pierden estabilidad incondicional u orden 2."
            )
        beta_eff, gamma_eff = _hht_autoderive_beta_gamma(alpha)
        if beta is not None or gamma is not None:
            # Override explícito de la autoderivación Hilber 1977: la propiedad
            # "estabilidad incondicional + orden 2" sólo se preserva con
            # β=(1-α)²/4, γ=(1-2α)/2. Cualquier otra combinación queda fuera
            # del régimen probado; avisamos en lugar de fallar silenciosamente.
            import warnings
            warnings.warn(
                f"HHTSolver: override explícito de beta/gamma "
                f"(beta={beta}, gamma={gamma}) reemplaza los valores "
                f"autoderivados desde alpha={alpha} "
                f"(β={beta_eff}, γ={gamma_eff}). Combinaciones distintas "
                "de las Hilber 1977 (β=(1-α)²/4, γ=(1-2α)/2) pueden romper "
                "la estabilidad incondicional o el orden 2.",
                RuntimeWarning,
                stacklevel=2,
            )
        if beta is not None:
            beta_eff = float(beta)
        if gamma is not None:
            gamma_eff = float(gamma)

        super().__init__(
            assembler, t_end, dt,
            beta=beta_eff, gamma=gamma_eff, rayleigh=rayleigh,
            u0=u0, u0_dot=u0_dot, F_func=F_func,
            linear_algebra=linear_algebra, lumping=lumping,
        )
        self.alpha = float(alpha)

    def solve(self) -> TransientResult:
        _log.info(f"--- INICIANDO SOLVER HHT-α (alpha={self.alpha:.4f}) ---")

        self.assembler.assemble_system()
        K = self.assembler.K_global
        M = self.assembler.assemble_mass_matrix(lumping=self.lumping)

        alpha_r, beta_r = resolve_rayleigh_config(
            self.rayleigh_cfg, source=type(self).__name__,
        )
        C = alpha_r * M + beta_r * K

        cs = self.assembler.constraint_set
        T, g = cs.build(self.assembler.ndof)
        free_dofs = cs.free_dofs(self.assembler.ndof)

        K_red = (T.T @ K @ T).tocsr()
        M_red = (T.T @ M @ T).tocsr()
        C_red = (T.T @ C @ T).tocsr() if (alpha_r != 0.0 or beta_r != 0.0) \
                else sp.csr_matrix(K_red.shape)
        F_dir = T.T @ (K @ g)

        ndof = self.assembler.ndof
        n_free = K_red.shape[0]

        u0_global = (np.zeros(ndof) if self.u0 is None
                      else np.asarray(self.u0, dtype=float).reshape(ndof))
        u0_dot_global = (np.zeros(ndof) if self.u0_dot is None
                          else np.asarray(self.u0_dot, dtype=float).reshape(ndof))
        u_free = u0_global[free_dofs].copy()
        udot_free = u0_dot_global[free_dofs].copy()

        F0_global = (np.zeros(ndof) if self.F_func is None
                      else np.asarray(self.F_func(0.0), dtype=float).reshape(ndof))
        F0_red = T.T @ F0_global

        # Aceleración inicial consistente con la ecuación de movimiento en t=0
        # (idéntica a Newmark: no involucra α porque no hay paso anterior).
        rhs0 = F0_red - F_dir - C_red @ udot_free - K_red @ u_free
        props_M = StiffnessProperties(
            is_symmetric=domain_is_symmetric(self.assembler.domain),
            is_positive_definite=True,
            size=n_free,
        )
        M_solver = select_solver(props_M, override=self.linear_algebra)
        uddot_free = M_solver.solve(M_red, rhs0)

        # Sistema efectivo HHT-α: A_eff = M + (1+α)γΔt·C + (1+α)βΔt²·K.
        # Constante en el tiempo → factorización única reutilizable.
        alpha_hht = self.alpha
        one_plus_alpha = 1.0 + alpha_hht
        dt = self.dt
        beta = self.beta
        gamma = self.gamma
        A_eff = (M_red
                 + one_plus_alpha * gamma * dt * C_red
                 + one_plus_alpha * beta * dt * dt * K_red).tocsr()
        props_A = StiffnessProperties(
            is_symmetric=domain_is_symmetric(self.assembler.domain),
            is_positive_definite=True,
            size=n_free,
        )
        A_solver = select_solver(props_A, override=self.linear_algebra)
        A_factor = A_solver.factorize(A_eff)

        n_steps = int(np.ceil(self.t_end / dt))
        t_history = np.linspace(0.0, n_steps * dt, n_steps + 1)
        u_history = np.zeros((ndof, n_steps + 1))
        udot_history = np.zeros((ndof, n_steps + 1))
        uddot_history = np.zeros((ndof, n_steps + 1))
        u_history[:, 0] = T @ u_free + g
        udot_history[:, 0] = T @ udot_free
        uddot_history[:, 0] = T @ uddot_free

        half_dt2 = 0.5 * dt * dt
        one_minus_2beta = 1.0 - 2.0 * beta
        one_minus_gamma = 1.0 - gamma
        beta_dt2 = beta * dt * dt
        gamma_dt = gamma * dt

        # Cache de F_n (carga del paso anterior) para los términos α·F_n.
        F_n_red = F0_red

        for step in range(n_steps):
            t_next = t_history[step + 1]

            # Predictores Newmark (idénticos al padre).
            u_pred = u_free + dt * udot_free + half_dt2 * one_minus_2beta * uddot_free
            udot_pred = udot_free + dt * one_minus_gamma * uddot_free

            # u_n, u̇_n son el estado al inicio del paso (antes del avance).
            u_n = u_free
            udot_n = udot_free

            F_global = (np.zeros(ndof) if self.F_func is None
                         else np.asarray(self.F_func(t_next), dtype=float).reshape(ndof))
            F_next_red = T.T @ F_global

            # RHS HHT-α: (1+α)·F_{n+1} − α·F_n − F_dir
            #            − (1+α)·C·ũ̇ + α·C·u̇_n
            #            − (1+α)·K·ũ + α·K·u_n
            rhs = (one_plus_alpha * F_next_red - alpha_hht * F_n_red
                   - F_dir
                   - one_plus_alpha * (C_red @ udot_pred) + alpha_hht * (C_red @ udot_n)
                   - one_plus_alpha * (K_red @ u_pred) + alpha_hht * (K_red @ u_n))

            uddot_free = A_factor.solve(rhs)

            # Correctores Newmark.
            u_free = u_pred + beta_dt2 * uddot_free
            udot_free = udot_pred + gamma_dt * uddot_free

            u_history[:, step + 1] = T @ u_free + g
            udot_history[:, step + 1] = T @ udot_free
            uddot_history[:, step + 1] = T @ uddot_free

            F_n_red = F_next_red

        _log.info(f"  -> {n_steps} pasos completados (HHT-α). "
                   f"alpha={alpha_hht:.4f}, ρ_∞={(one_plus_alpha)/(1.0-alpha_hht):.4f}. "
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
class NewtonHHTSolver(NewtonNewmarkSolver):
    """Análisis dinámico transitorio **no lineal** por HHT-α + Newton.

    Variante de :class:`NewtonNewmarkSolver` que aplica el esquema temporal
    HHT-α al residuo dinámico no lineal. Residuo y jacobiano modificados:

    .. math::

        \\mathbf R(\\ddot{\\mathbf u}_{n+1}) =
          (1+\\alpha)\\mathbf F_{n+1} - \\alpha\\mathbf F_n
          - [(1+\\alpha)\\mathbf F_\\text{int}(\\mathbf u_{n+1}) - \\alpha\\mathbf F_\\text{int}(\\mathbf u_n)]
          - [(1+\\alpha)\\mathbf C\\dot{\\mathbf u}_{n+1} - \\alpha\\mathbf C\\dot{\\mathbf u}_n]
          - \\mathbf M\\ddot{\\mathbf u}_{n+1}

    .. math::

        \\mathbf J = \\mathbf M + (1+\\alpha)\\gamma\\Delta t\\,\\mathbf C
                                  + (1+\\alpha)\\beta\\Delta t^2\\,\\mathbf K_\\text{t}

    Todo lo demás (criterio de convergencia, line search opt-in, telemetría
    tipada, commit/rollback de estado, Rayleigh con K_0 constante) se hereda
    de ``NewtonNewmarkSolver`` sin cambios.

    Parameters
    ----------
    alpha : float, default −0.05
        Parámetro HHT (ver `HHTSolver`).
    beta, gamma : float or None
        Auto-derivados desde alpha si `None`.
    (resto heredado de NewtonNewmarkSolver)

    Notes
    -----
    Con materiales lineales, este solver reproduce exactamente `HHTSolver`
    (el residuo de Newton se anula en una iteración). Validado en tests.

    Ver ``docs/specs/HHTSolver.md``.
    """

    def __init__(
        self,
        assembler,
        t_end: float,
        dt: float,
        *,
        alpha: float = -0.05,
        convergence: ConvergenceCriterion | None = None,
        max_iter: int = 20,
        freeze_tangent_after_iter: int | None = None,
        line_search: bool = False,
        beta: float | None = None,
        gamma: float | None = None,
        rayleigh: dict | None = None,
        u0: np.ndarray | None = None,
        u0_dot: np.ndarray | None = None,
        F_func: Callable[[float], np.ndarray] | None = None,
        linear_algebra: str = "auto",
        lumping: str = "consistent",
    ):
        if not (-1.0 / 3.0 - 1.0e-12 <= alpha <= 0.0 + 1.0e-12):
            raise ValueError(
                f"NewtonHHTSolver: alpha={alpha} fuera del rango [-1/3, 0]. "
                "Valores fuera de este rango pierden estabilidad incondicional u orden 2."
            )
        beta_eff, gamma_eff = _hht_autoderive_beta_gamma(alpha)
        if beta is not None or gamma is not None:
            # Mismo aviso que HHTSolver (override fuera del régimen Hilber 1977).
            import warnings
            warnings.warn(
                f"NewtonHHTSolver: override explícito de beta/gamma "
                f"(beta={beta}, gamma={gamma}) reemplaza los valores "
                f"autoderivados desde alpha={alpha} "
                f"(β={beta_eff}, γ={gamma_eff}). Combinaciones distintas "
                "de las Hilber 1977 (β=(1-α)²/4, γ=(1-2α)/2) pueden romper "
                "la estabilidad incondicional o el orden 2.",
                RuntimeWarning,
                stacklevel=2,
            )
        if beta is not None:
            beta_eff = float(beta)
        if gamma is not None:
            gamma_eff = float(gamma)

        super().__init__(
            assembler, t_end, dt,
            convergence=convergence, max_iter=max_iter,
            freeze_tangent_after_iter=freeze_tangent_after_iter,
            line_search=line_search,
            beta=beta_eff, gamma=gamma_eff, rayleigh=rayleigh,
            u0=u0, u0_dot=u0_dot, F_func=F_func,
            linear_algebra=linear_algebra, lumping=lumping,
        )
        self.alpha = float(alpha)

    def solve(self) -> TransientResult:
        _log.info(f"--- INICIANDO SOLVER HHT-α NO LINEAL (alpha={self.alpha:.4f}) ---")

        if self.assembler.domain.total_dofs == 0:
            self.assembler.domain.generate_equation_numbers()
        ndof = self.assembler.domain.total_dofs

        u_total = (np.zeros(ndof) if self.u0 is None
                    else np.asarray(self.u0, dtype=float).reshape(ndof)).copy()
        udot_total = (np.zeros(ndof) if self.u0_dot is None
                       else np.asarray(self.u0_dot, dtype=float).reshape(ndof)).copy()

        M = self.assembler.assemble_mass_matrix(lumping=self.lumping)

        K_0, F_int = self.assembler.assemble_non_linear_system(u_total)

        alpha_r, beta_r = resolve_rayleigh_config(
            self.rayleigh_cfg, source=type(self).__name__,
        )
        C = alpha_r * M + beta_r * K_0

        cs = self.assembler.constraint_set
        T_op, g_vec = cs.build(ndof)
        free_dofs = cs.free_dofs(ndof)
        n_free = T_op.shape[1]

        prescribed = np.ones(ndof, dtype=bool)
        prescribed[free_dofs] = False
        u_total[prescribed] = g_vec[prescribed]
        udot_total[prescribed] = 0.0

        M_red = (T_op.T @ M @ T_op).tocsr()
        if alpha_r != 0.0 or beta_r != 0.0:
            C_red = (T_op.T @ C @ T_op).tocsr()
        else:
            C_red = sp.csr_matrix(M_red.shape)

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

        # Estado del paso anterior (al inicio del análisis: instante 0).
        # Usado en los términos α·X_n de HHT.
        F_int_prev = F_int.copy()
        udot_prev = udot_total.copy()
        F_prev_global = F0_global.copy()

        dt = self.dt
        beta = self.beta
        gamma = self.gamma
        alpha_hht = self.alpha
        one_plus_alpha = 1.0 + alpha_hht
        n_steps = int(np.ceil(self.t_end / dt))
        t_history = np.linspace(0.0, n_steps * dt, n_steps + 1)
        u_history = np.zeros((ndof, n_steps + 1))
        udot_history = np.zeros((ndof, n_steps + 1))
        uddot_history = np.zeros((ndof, n_steps + 1))
        u_history[:, 0] = u_total
        udot_history[:, 0] = udot_total
        uddot_history[:, 0] = uddot_total

        half_dt2 = 0.5 * dt * dt
        one_minus_2beta = 1.0 - 2.0 * beta
        one_minus_gamma = 1.0 - gamma
        beta_dt2 = beta * dt * dt
        gamma_dt = gamma * dt

        is_sym = domain_is_symmetric(self.assembler.domain)
        is_pd = True
        props_J = StiffnessProperties(
            is_symmetric=is_sym, is_positive_definite=is_pd, size=n_free,
        )
        linalg = select_solver(props_J, override=self.linear_algebra)
        frozen_factor = None

        for step in range(n_steps):
            t_next = t_history[step + 1]

            u_pred = u_total + dt * udot_total + half_dt2 * one_minus_2beta * uddot_total
            udot_pred = udot_total + dt * one_minus_gamma * uddot_total

            uddot_iter = uddot_total.copy()
            u_iter = u_pred + beta_dt2 * uddot_iter
            udot_iter = udot_pred + gamma_dt * uddot_iter

            F_ext_next = (np.zeros(ndof) if self.F_func is None
                           else np.asarray(self.F_func(t_next), dtype=float).reshape(ndof))

            residual_history: list[float] = []
            delta_history: list[float] = []
            singular_tangent_seen = False
            last_residual = float("inf")
            last_delta = 0.0

            converged = False
            for it in range(self.max_iter):
                K_t, F_int_iter = self.assembler.assemble_non_linear_system(u_iter)

                # Residuo HHT-α no lineal:
                # R = (1+α)·F_{n+1} − α·F_n
                #     − [(1+α)·F_int(u_{n+1}) − α·F_int_n]
                #     − [(1+α)·C·u̇_{n+1} − α·C·u̇_n]
                #     − M·ü_{n+1}
                R = (one_plus_alpha * F_ext_next - alpha_hht * F_prev_global
                     - one_plus_alpha * F_int_iter + alpha_hht * F_int_prev
                     - one_plus_alpha * (C @ udot_iter) + alpha_hht * (C @ udot_prev)
                     - M @ uddot_iter)
                R_red = T_op.T @ R

                # Jacobiano: J = M + (1+α)·γΔt·C + (1+α)·βΔt²·K_t
                J = (M_red
                     + one_plus_alpha * gamma_dt * C_red
                     + one_plus_alpha * beta_dt2 * (T_op.T @ K_t @ T_op)).tocsr()

                if not self.convergence.is_calibrated:
                    force_scale = max(
                        np.linalg.norm(F_ext_next),
                        np.linalg.norm(F_int_iter),
                        1.0,
                    )
                    K_diag = stiffness_diag_scale(K_t)
                    disp_scale = force_scale / K_diag
                    self.convergence.calibrate(force_scale, disp_scale)

                # Tangente dinámica singular (RuntimeError de `splu` tras
                # degradar a LU): flipea el flag para `classify_divergence`
                # → `SingularTangentError` (ADR 0011).
                threshold = self.freeze_tangent_after_iter
                try:
                    try:
                        if threshold is None:
                            delta_uddot_red = linalg.solve(J, R_red)
                        elif it < threshold or frozen_factor is None:
                            frozen_factor = linalg.factorize(J)
                            delta_uddot_red = frozen_factor.solve(R_red)
                        else:
                            delta_uddot_red = frozen_factor.solve(R_red)
                    except CholeskyNotPositiveDefiniteError:
                        _log.warning("Cholesky no-PD en NewtonHHT; degradando a LU.")
                        is_pd = False
                        props_J = StiffnessProperties(
                            is_symmetric=is_sym, is_positive_definite=False, size=n_free,
                        )
                        linalg = select_solver(props_J, override=self.linear_algebra)
                        frozen_factor = None
                        delta_uddot_red = linalg.solve(J, R_red)
                except RuntimeError:
                    _log.error("Tangente dinámica singular en NewtonHHT.")
                    singular_tangent_seen = True
                    break

                # Line search (opt-in, ADR 0011). Mismo helper que NewtonNewmark.
                R_norm_before = float(np.linalg.norm(R[free_dofs]))
                alpha_ls = self._armijo_step_dynamic(
                    uddot_iter, delta_uddot_red, u_pred, udot_pred,
                    beta_dt2, gamma_dt, free_dofs,
                    F_ext_next, C, M, R_norm_before,
                )

                uddot_iter[free_dofs] += alpha_ls * delta_uddot_red
                u_iter = u_pred + beta_dt2 * uddot_iter
                udot_iter = udot_pred + gamma_dt * uddot_iter

                # Re-ensamblar para R coherente con el U avanzado.
                _, F_int_after = self.assembler.assemble_non_linear_system(u_iter)
                R_after = (one_plus_alpha * F_ext_next - alpha_hht * F_prev_global
                           - one_plus_alpha * F_int_after + alpha_hht * F_int_prev
                           - one_plus_alpha * (C @ udot_iter) + alpha_hht * (C @ udot_prev)
                           - M @ uddot_iter)
                R_norm_after = float(np.linalg.norm(R_after[free_dofs]))

                delta_u_norm = beta_dt2 * float(np.linalg.norm(alpha_ls * delta_uddot_red))
                ref_force = max(np.linalg.norm(F_ext_next), np.linalg.norm(F_int_after))
                state = self.convergence.evaluate(
                    residual_norm=R_norm_after,
                    ref_force=ref_force,
                    delta_u_norm=delta_u_norm,
                    u_norm=np.linalg.norm(u_iter[free_dofs]),
                )
                residual_history.append(R_norm_after)
                delta_history.append(delta_u_norm)
                last_residual = R_norm_after
                last_delta = delta_u_norm

                if state.converged:
                    _log.info(
                        f"  [PASO {step+1}/{n_steps}] t={t_next:.4e} | "
                        f"iter={it+1} | R/tol_F={state.ratio_force:.2e}"
                    )
                    self.assembler.commit_all_states()
                    converged = True
                    # Actualizar F_int_iter para almacenarlo como F_int_prev
                    F_int_iter = F_int_after
                    break

            frozen_factor = None

            if not converged:
                err_cls = classify_divergence(
                    residual_history, delta_history,
                    singular_tangent_detected=singular_tangent_seen,
                )
                raise err_cls(
                    last_residual=last_residual,
                    last_delta=last_delta,
                    last_load_factor=t_next,
                    n_bisections=0,
                    extra_message=(
                        f"NewtonHHTSolver (alpha={alpha_hht:.4f}): "
                        f"paso {step+1} (t={t_next:.4e}) no convergió en "
                        f"{self.max_iter} iteraciones; "
                        f"line_search={'on' if self.line_search else 'off'}."
                    ),
                )

            u_total = u_iter
            udot_total = udot_iter
            uddot_total = uddot_iter

            u_history[:, step + 1] = u_total
            udot_history[:, step + 1] = udot_total
            uddot_history[:, step + 1] = uddot_total

            # Cache para el siguiente paso (términos α·X_n).
            F_int_prev = F_int_iter
            udot_prev = udot_total.copy()
            F_prev_global = F_ext_next

        _log.info(f"  -> {n_steps} pasos completados (HHT-α no lineal). "
                  f"alpha={alpha_hht:.4f}, ρ_∞={(one_plus_alpha)/(1.0-alpha_hht):.4f}. "
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
