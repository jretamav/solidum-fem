# fenix_fem/fenix/math/solvers/arclength.py
"""``ArcLengthSolver`` — método de longitud de arco cilíndrico (Crisfield).
"""
from __future__ import annotations

import numpy as np

from fenix.constants import ARCLENGTH_MIN_DL_FACTOR, ZERO_TOL
from fenix.math.convergence import (
    ConvergenceCriterion,
    stiffness_diag_scale,
)
from fenix.math.linalg import LUSolver, StiffnessProperties, select_solver
from fenix.math.solvers._shared import (
    CholeskyNotPositiveDefiniteError,
    _log,
    domain_is_symmetric,
)
from fenix.registry import SolverRegistry


@SolverRegistry.register
class ArcLengthSolver:
    """
    Solucionador no lineal con Método de Longitud de Arco Cilíndrico (Crisfield).
    Permite trazar curvas de equilibrio con fenómenos de snap-through y snap-back
    variando simultáneamente los desplazamientos y la carga externa.
    """
    def __init__(self, assembler, convergence: ConvergenceCriterion | None = None, max_iter=20, max_lambda=1.0, initial_dl=0.1, max_steps=100,
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
            is_symmetric=domain_is_symmetric(self.assembler.domain),
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
                if final_step:
                    dU_iter = dU_iter + dU_update
                # En la rama no-final `dU_iter` ya quedó actualizado arriba.
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
                if dl < ARCLENGTH_MIN_DL_FACTOR * self.dl:
                    raise RuntimeError("Arc-Length fracasó irreparablemente.")

        return U_current
