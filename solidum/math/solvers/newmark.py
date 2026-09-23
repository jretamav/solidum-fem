# solidum_fem/solidum/math/solvers/newmark.py
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

from solidum.math.convergence import ConvergenceCriterion
from solidum.math.damping import resolve_rayleigh_config
from solidum.math.linalg import StiffnessProperties, select_solver
from solidum.math.solvers._shared import (
    ElementForcesRecorder,
    _log,
    domain_is_symmetric,
    number_of_steps,
    solve_mass_system,
)
from solidum.math.solvers.corrector import NewtonCorrector, default_calibration_scales
from solidum.registry import SolverRegistry
from solidum.results import TransientResult


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
        record_internal_forces: bool = False,
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
        # Registro de fuerzas internas por paso (ADR 0002 / deuda #15):
        # opt-in porque cuesta una evaluacion por elemento 1D y paso.
        self.record_internal_forces = bool(record_internal_forces)


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
            near_nullspace=self.assembler.near_nullspace,
        )
        M_solver = select_solver(props_M, override=self.linear_algebra)
        uddot_free = solve_mass_system(M_solver, M_red, rhs0, type(self).__name__)

        # Factorización reutilizable de A_eff = M + γΔt·C + βΔt²·K.
        dt = self.dt
        beta = self.beta
        gamma = self.gamma
        A_eff = (M_red + gamma * dt * C_red + beta * dt * dt * K_red).tocsr()
        props_A = StiffnessProperties(
            is_symmetric=domain_is_symmetric(self.assembler.domain),
            is_positive_definite=True,
            size=n_free,
            near_nullspace=self.assembler.near_nullspace,
        )
        A_solver = select_solver(props_A, override=self.linear_algebra)
        A_factor = A_solver.factorize(A_eff)

        # Historiales en globales. n_steps + 1 columnas (incluye t=0).
        n_steps = number_of_steps(self.t_end, dt)
        t_history = np.linspace(0.0, n_steps * dt, n_steps + 1)
        u_history = np.zeros((ndof, n_steps + 1))
        udot_history = np.zeros((ndof, n_steps + 1))
        uddot_history = np.zeros((ndof, n_steps + 1))
        # Estado inicial en globales (con apoyos g).
        u_history[:, 0] = T @ u_free + g
        forces_rec = ElementForcesRecorder(self.assembler.domain, self.record_internal_forces)
        forces_rec.record(u_history[:, 0])
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
            forces_rec.record(u_history[:, step + 1])
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
            element_forces_history=forces_rec.result(),
        )



class _DynamicNewtonProblem:
    """Paso de Newton-Newmark / Newton-HHT en el protocolo ``NewtonProblem``
    (ADR 0015).

    El iterado ``x`` es ``(u, u̇, ü, ü_libre)``; la incógnita del Newton es
    ``ü_libre`` y ``u``, ``u̇`` se reconstruyen con los correctores de Newmark
    (``u = ũ + βΔt²·ü``, ``u̇ = ũ̇ + γΔt·ü``) y con ``T`` (Dirichlet y MPC).
    El estado es ``(K_t, F_int)``. ``residual_fn(F_int, u̇, ü)`` es la forma
    del residuo dinámico del solver (Newmark o HHT-α) y ``jac_scale`` el
    factor ``(1+α)`` que HHT aplica a los términos de rigidez y
    amortiguamiento del jacobiano.
    """

    def __init__(self, assembler, T_op, free_dofs, M, C, M_red, C_red,
                 u_pred, udot_pred, beta_dt2, gamma_dt, F_ext_next, residual_fn,
                 jac_scale: float = 1.0):
        self.assembler = assembler
        self.T_op = T_op
        self.free_dofs = free_dofs
        self.M = M
        self.C = C
        self.M_red = M_red
        self.C_red = C_red
        self.u_pred = u_pred
        self.udot_pred = udot_pred
        self.beta_dt2 = float(beta_dt2)
        self.gamma_dt = float(gamma_dt)
        self.F_ext_next = F_ext_next
        self.residual_fn = residual_fn
        self.jac_scale = float(jac_scale)

    def initial_iterate(self, uddot_free):
        """``ü^(0) = ü_n`` (continuidad) y los ``u``, ``u̇`` correspondientes."""
        return self._iterate(np.array(uddot_free, dtype=float, copy=True))

    def _iterate(self, uddot_free):
        uddot = np.asarray(self.T_op @ uddot_free, dtype=float)
        u = self.u_pred + self.beta_dt2 * uddot
        udot = self.udot_pred + self.gamma_dt * uddot
        return (u, udot, uddot, uddot_free)

    # --- protocolo ----------------------------------------------------------

    def assemble(self, x):
        return self.assembler.assemble_non_linear_system(x[0])

    def residual(self, x, state):
        return self.residual_fn(state[1], x[1], x[2])

    def residual_norm(self, R):
        return float(np.linalg.norm(self.T_op.T @ R))

    def calibration_scales(self, x, state):
        return default_calibration_scales(self.F_ext_next, state[1], state[0])

    def reference_force(self, x, state):
        return max(float(np.linalg.norm(self.F_ext_next)), float(np.linalg.norm(state[1])))

    def x_norm(self, x):
        return float(np.linalg.norm(x[0][self.free_dofs]))

    def correction(self, x, state, R, solve):
        # Jacobiano dinámico: J = M + s·γΔt·C + s·βΔt²·K_t (s = 1 en
        # Newmark, 1+α en HHT). Una tangente singular la señala el backend
        # con RuntimeError, que el corrector clasifica (ADR 0011).
        T_op = self.T_op
        J = (self.M_red
             + self.jac_scale * self.gamma_dt * self.C_red
             + self.jac_scale * self.beta_dt2 * (T_op.T @ state[0] @ T_op)).tocsr()
        return solve(J, T_op.T @ R)

    def apply(self, x, delta_uddot_red, alpha):
        step = alpha * delta_uddot_red
        # δu corresponde a βΔt² · α · δü (cambio en desplazamiento por la iter).
        return self._iterate(x[3] + step), self.beta_dt2 * float(np.linalg.norm(step))

    def on_converged(self, x, state):
        self.assembler.commit_all_states()

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
        record_internal_forces: bool = False,
    ):
        super().__init__(
            assembler, t_end, dt,
            beta=beta, gamma=gamma, rayleigh=rayleigh,
            u0=u0, u0_dot=u0_dot, F_func=F_func,
            linear_algebra=linear_algebra, lumping=lumping,
            record_internal_forces=record_internal_forces,
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
        # Corrector compartido (ADR 0015): bucle de Newton, backend
        # algebraico con degradación a LU, Newton modificado y line search.
        # Sin línea por iteración: el solver registra una por paso.
        self.corrector = NewtonCorrector(
            self.convergence,
            max_iter=self.max_iter,
            is_symmetric=domain_is_symmetric(assembler.domain),
            is_positive_definite=True,
            linear_algebra=self.linear_algebra,
            freeze_tangent_after_iter=self.freeze_tangent_after_iter,
            line_search=self.line_search,
            verbose=False,
            near_nullspace=assembler.near_nullspace,
        )

    def solve(self) -> TransientResult:
        _log.info("--- INICIANDO SOLVER NEWMARK NO LINEAL (Newton dentro de Newmark) ---")

        # Disparar el build de topología del assembler antes de leer ndof.
        # `assembler.ndof` solo se rellena tras un primer ensamblaje.
        if self.assembler.domain.total_dofs == 0:
            self.assembler.domain.generate_equation_numbers()
        ndof = self.assembler.domain.total_dofs

        # Reducción por Dirichlet / MPC (ADR 0004): ``u = T·u_free + g``,
        # ``u̇ = T·u̇_free``, ``ü = T·ü_free``. Trabajar con la incógnita
        # reducida y reconstruir los vectores completos con ``T`` cubre a la
        # vez los apoyos (fila nula + g) y las restricciones lineales con
        # maestros (fila con los coeficientes α_si); actualizar sólo
        # ``free_dofs`` dejaba los esclavos congelados y el Newton no podía
        # anular el residuo con MPC (auditoría 2026-09-22).
        cs = self.assembler.constraint_set
        T_op, g_vec = cs.build(ndof)
        free_dofs = cs.free_dofs(ndof)
        n_free = T_op.shape[1]

        # Condiciones iniciales del usuario proyectadas sobre la variedad de
        # restricciones ANTES de cualquier ensamblaje: se conservan sus
        # componentes libres y los DOFs no libres se reconstruyen (apoyo
        # constante ⇒ g; esclavo ⇒ combinación de maestros). Ensamblar
        # ``F_int(u₀)`` antes de imponer ``g`` daba ü₀ inconsistente con
        # apoyos prescritos no nulos.
        u0_global = (np.zeros(ndof) if self.u0 is None
                     else np.asarray(self.u0, dtype=float).reshape(ndof))
        udot0_global = (np.zeros(ndof) if self.u0_dot is None
                        else np.asarray(self.u0_dot, dtype=float).reshape(ndof))
        u_total = np.asarray(T_op @ u0_global[free_dofs] + g_vec, dtype=float)
        udot_total = np.asarray(T_op @ udot0_global[free_dofs], dtype=float)

        # Masa (constante).
        M = self.assembler.assemble_mass_matrix(lumping=self.lumping)

        # Sistema no lineal en u_0 para obtener K_0 (Rayleigh) y F_int(u_0).
        K_0, F_int = self.assembler.assemble_non_linear_system(u_total)

        alpha_r, beta_r = resolve_rayleigh_config(
            self.rayleigh_cfg, source=type(self).__name__,
        )
        C = alpha_r * M + beta_r * K_0

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
            near_nullspace=self.assembler.near_nullspace,
        )
        M_solver = select_solver(props_M, override=self.linear_algebra)
        uddot_free = solve_mass_system(M_solver, M_red, rhs0, type(self).__name__)
        uddot_total = np.asarray(T_op @ uddot_free, dtype=float)

        # Historiales.
        dt = self.dt
        beta = self.beta
        gamma = self.gamma
        n_steps = number_of_steps(self.t_end, dt)
        t_history = np.linspace(0.0, n_steps * dt, n_steps + 1)
        u_history = np.zeros((ndof, n_steps + 1))
        udot_history = np.zeros((ndof, n_steps + 1))
        uddot_history = np.zeros((ndof, n_steps + 1))
        u_history[:, 0] = u_total
        forces_rec = ElementForcesRecorder(self.assembler.domain, self.record_internal_forces)
        forces_rec.record(u_history[:, 0])
        udot_history[:, 0] = udot_total
        uddot_history[:, 0] = uddot_total

        # Coeficientes Newmark precomputados.
        half_dt2 = 0.5 * dt * dt
        one_minus_2beta = 1.0 - 2.0 * beta
        one_minus_gamma = 1.0 - gamma
        beta_dt2 = beta * dt * dt
        gamma_dt = gamma * dt

        def residual(F_int_at, udot_at, uddot_at):
            return F_ext_next - F_int_at - C @ udot_at - M @ uddot_at

        for step in range(n_steps):
            t_next = t_history[step + 1]

            # ADR 0010 §5: hook de preparación de paso con el estado
            # convergido del paso anterior (activación de discontinuidades
            # embebidas). Mismo protocolo que los solvers estáticos.
            self.assembler.prepare_all_steps(u_total)

            # Predictores Newmark.
            u_pred = u_total + dt * udot_total + half_dt2 * one_minus_2beta * uddot_total
            udot_pred = udot_total + dt * one_minus_gamma * uddot_total

            F_ext_next = (np.zeros(ndof) if self.F_func is None
                           else np.asarray(self.F_func(t_next), dtype=float).reshape(ndof))

            # Corrector compartido (ADR 0015): un ensamblaje por iteración,
            # convergencia evaluada ANTES de resolver con el par
            # (‖R(u_k)‖, ‖δu_{k−1}‖) y commit del estado trial del ensamblaje
            # convergido. El Newton arranca en ü^(0) = ü_n (continuidad).
            problem = _DynamicNewtonProblem(
                self.assembler, T_op, free_dofs, M, C, M_red, C_red,
                u_pred, udot_pred, beta_dt2, gamma_dt, F_ext_next, residual,
            )
            res = self.corrector.run(problem, problem.initial_iterate(uddot_free))

            if not res.converged:
                # Clasificar el modo y lanzar excepción tipada (ADR 0011).
                raise res.divergence_error(
                    last_load_factor=t_next,
                    extra_message=(
                        f"NewtonNewmarkSolver: paso {step+1} (t={t_next:.4e}) no "
                        f"convergió en {self.max_iter} iteraciones; "
                        f"line_search={'on' if self.line_search else 'off'}; "
                        f"último α={res.last_alpha:.3f}."
                    ),
                )
            _log.info(
                f"  [PASO {step+1}/{n_steps}] t={t_next:.4e} | "
                f"iter={res.n_solves} | R/tol_F={res.last_conv.ratio_force:.2e}"
            )

            u_total, udot_total, uddot_total, uddot_free = res.x

            u_history[:, step + 1] = u_total
            forces_rec.record(u_history[:, step + 1])
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
            element_forces_history=forces_rec.result(),
        )


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
        record_internal_forces: bool = False,
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
            record_internal_forces=record_internal_forces,
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
            near_nullspace=self.assembler.near_nullspace,
        )
        M_solver = select_solver(props_M, override=self.linear_algebra)
        uddot_free = solve_mass_system(M_solver, M_red, rhs0, type(self).__name__)

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
            near_nullspace=self.assembler.near_nullspace,
        )
        A_solver = select_solver(props_A, override=self.linear_algebra)
        A_factor = A_solver.factorize(A_eff)

        n_steps = number_of_steps(self.t_end, dt)
        t_history = np.linspace(0.0, n_steps * dt, n_steps + 1)
        u_history = np.zeros((ndof, n_steps + 1))
        udot_history = np.zeros((ndof, n_steps + 1))
        uddot_history = np.zeros((ndof, n_steps + 1))
        u_history[:, 0] = T @ u_free + g
        forces_rec = ElementForcesRecorder(self.assembler.domain, self.record_internal_forces)
        forces_rec.record(u_history[:, 0])
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
            forces_rec.record(u_history[:, step + 1])
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
            element_forces_history=forces_rec.result(),
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
        record_internal_forces: bool = False,
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
            record_internal_forces=record_internal_forces,
        )
        self.alpha = float(alpha)

    def solve(self) -> TransientResult:
        _log.info(f"--- INICIANDO SOLVER HHT-α NO LINEAL (alpha={self.alpha:.4f}) ---")

        if self.assembler.domain.total_dofs == 0:
            self.assembler.domain.generate_equation_numbers()
        ndof = self.assembler.domain.total_dofs

        # Reducción por Dirichlet / MPC y proyección de las condiciones
        # iniciales sobre la variedad de restricciones antes de ensamblar.
        # Ver NewtonNewmarkSolver.solve para la justificación.
        cs = self.assembler.constraint_set
        T_op, g_vec = cs.build(ndof)
        free_dofs = cs.free_dofs(ndof)
        n_free = T_op.shape[1]

        u0_global = (np.zeros(ndof) if self.u0 is None
                     else np.asarray(self.u0, dtype=float).reshape(ndof))
        udot0_global = (np.zeros(ndof) if self.u0_dot is None
                        else np.asarray(self.u0_dot, dtype=float).reshape(ndof))
        u_total = np.asarray(T_op @ u0_global[free_dofs] + g_vec, dtype=float)
        udot_total = np.asarray(T_op @ udot0_global[free_dofs], dtype=float)

        M = self.assembler.assemble_mass_matrix(lumping=self.lumping)

        K_0, F_int = self.assembler.assemble_non_linear_system(u_total)

        alpha_r, beta_r = resolve_rayleigh_config(
            self.rayleigh_cfg, source=type(self).__name__,
        )
        C = alpha_r * M + beta_r * K_0

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
            near_nullspace=self.assembler.near_nullspace,
        )
        M_solver = select_solver(props_M, override=self.linear_algebra)
        uddot_free = solve_mass_system(M_solver, M_red, rhs0, type(self).__name__)
        uddot_total = np.asarray(T_op @ uddot_free, dtype=float)

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
        n_steps = number_of_steps(self.t_end, dt)
        t_history = np.linspace(0.0, n_steps * dt, n_steps + 1)
        u_history = np.zeros((ndof, n_steps + 1))
        udot_history = np.zeros((ndof, n_steps + 1))
        uddot_history = np.zeros((ndof, n_steps + 1))
        u_history[:, 0] = u_total
        forces_rec = ElementForcesRecorder(self.assembler.domain, self.record_internal_forces)
        forces_rec.record(u_history[:, 0])
        udot_history[:, 0] = udot_total
        uddot_history[:, 0] = uddot_total

        half_dt2 = 0.5 * dt * dt
        one_minus_2beta = 1.0 - 2.0 * beta
        one_minus_gamma = 1.0 - gamma
        beta_dt2 = beta * dt * dt
        gamma_dt = gamma * dt

        def residual(F_int_at, udot_at, uddot_at):
            # Residuo HHT-α no lineal:
            # R = (1+α)·F_{n+1} − α·F_n
            #     − [(1+α)·F_int(u_{n+1}) − α·F_int_n]
            #     − [(1+α)·C·u̇_{n+1} − α·C·u̇_n]
            #     − M·ü_{n+1}
            return (one_plus_alpha * F_ext_next - alpha_hht * F_prev_global
                    - one_plus_alpha * F_int_at + alpha_hht * F_int_prev
                    - one_plus_alpha * (C @ udot_at) + alpha_hht * (C @ udot_prev)
                    - M @ uddot_at)

        for step in range(n_steps):
            t_next = t_history[step + 1]

            # ADR 0010 §5: hook de preparación de paso.
            self.assembler.prepare_all_steps(u_total)

            u_pred = u_total + dt * udot_total + half_dt2 * one_minus_2beta * uddot_total
            udot_pred = udot_total + dt * one_minus_gamma * uddot_total

            F_ext_next = (np.zeros(ndof) if self.F_func is None
                           else np.asarray(self.F_func(t_next), dtype=float).reshape(ndof))

            # Corrector compartido (ADR 0015); jacobiano
            # J = M + (1+α)·γΔt·C + (1+α)·βΔt²·K_t.
            problem = _DynamicNewtonProblem(
                self.assembler, T_op, free_dofs, M, C, M_red, C_red,
                u_pred, udot_pred, beta_dt2, gamma_dt, F_ext_next, residual,
                jac_scale=one_plus_alpha,
            )
            res = self.corrector.run(problem, problem.initial_iterate(uddot_free))

            if not res.converged:
                raise res.divergence_error(
                    last_load_factor=t_next,
                    extra_message=(
                        f"NewtonHHTSolver (alpha={alpha_hht:.4f}): "
                        f"paso {step+1} (t={t_next:.4e}) no convergió en "
                        f"{self.max_iter} iteraciones; "
                        f"line_search={'on' if self.line_search else 'off'}."
                    ),
                )
            _log.info(
                f"  [PASO {step+1}/{n_steps}] t={t_next:.4e} | "
                f"iter={res.n_solves} | R/tol_F={res.last_conv.ratio_force:.2e}"
            )

            u_total, udot_total, uddot_total, uddot_free = res.x

            u_history[:, step + 1] = u_total
            forces_rec.record(u_history[:, step + 1])
            udot_history[:, step + 1] = udot_total
            uddot_history[:, step + 1] = uddot_total

            # Cache para el siguiente paso (términos α·X_n). F_int es el del
            # ensamblaje en u_{n+1} convergido.
            F_int_prev = res.state[1]
            udot_prev = udot_total.copy()
            F_prev_global = F_ext_next

        _log.info(f"  -> {n_steps} pasos completados (HHT-α no lineal). "
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
            element_forces_history=forces_rec.result(),
        )
