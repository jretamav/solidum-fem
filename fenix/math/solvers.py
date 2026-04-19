# fenix_fem/fenix/math/solvers.py
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from fenix.constants import CONVERGENCE_TOL, ZERO_TOL
from fenix.registry import SolverRegistry


@SolverRegistry.register
class LinearSolver:
    """Solucionador de sistemas algebraicos lineales en un solo paso."""
    def __init__(self, assembler, penalty_value=1e15):
        self.assembler = assembler
        self.penalty = penalty_value

    def solve(self, F_ext_global: np.ndarray) -> np.ndarray:
        print("\n--- INICIANDO SOLVER LINEAL ---")
        self.assembler.assemble_system()
        K_global = self.assembler.K_global.copy()
        R = F_ext_global.copy()

        K_global, R = self.assembler.apply_bcs_to_system(K_global, R, penalty_value=self.penalty)
        
        U = spla.spsolve(K_global, R)
        print("  -> CONVERGENCIA ALCANZADA (1 Iteración).")
        return U

@SolverRegistry.register
class NonlinearSolver:
    """Solucionador incremental-iterativo de Newton-Raphson con paso adaptativo."""
    def __init__(self, assembler, tol=CONVERGENCE_TOL, max_iter=20, num_steps=10, adaptive=True, penalty_value=1e15, min_delta_lambda=1e-5):
        self.assembler = assembler
        self.tol = tol
        self.max_iter = max_iter
        self.num_steps = num_steps
        self.adaptive = adaptive
        self.penalty = penalty_value
        self.min_delta_lambda = min_delta_lambda

    def solve(self, F_ext_global: np.ndarray, step_callback=None) -> np.ndarray:
        domain = self.assembler.domain
        ndof = domain.total_dofs
        U_current = np.zeros(ndof)
        
        print("\n--- INICIANDO SOLVER NO LINEAL (CONTROL DE PASO ADAPTATIVO) ---")
        
        load_factor = 0.0
        target_load = 1.0
        delta_lambda = 1.0 / self.num_steps
        step = 0
        
        while load_factor < target_load - 1e-9:
            step += 1
            
            if load_factor + delta_lambda > target_load:
                delta_lambda = target_load - load_factor
                
            next_load_factor = load_factor + delta_lambda
            print(f"\n[PASO {step}] Intentando Factor de Carga: {next_load_factor:.4f} (Incremento: {delta_lambda:.4f})")
            
            F_ext_step = F_ext_global * next_load_factor
            U_iter = U_current.copy()
            converged = False
            
            for iteration in range(self.max_iter):
                K_global, F_int_global = self.assembler.assemble_non_linear_system(U_iter)
                R = F_ext_step - F_int_global
                
                K_global, R = self.assembler.apply_bcs_to_system(K_global, R, penalty_value=self.penalty, load_factor=next_load_factor, U_current=U_iter)
                
                try:
                    delta_U = spla.spsolve(K_global, R)
                except RuntimeError:
                    print("  -> Error: Matriz Singular detectada.")
                    break
                
                U_iter += delta_U

                # Criterio dual: norma de desplazamiento Y norma de residuo de fuerza.
                # Ambos deben converger para garantizar equilibrio.
                # En problemas controlados por desplazamientos F_ext_step ≈ 0; usamos
                # max(|F_ext|, |F_int|) como referencia para no dividir por ~ZERO_TOL.
                err_disp = np.linalg.norm(delta_U) / (np.linalg.norm(U_iter) + ZERO_TOL)
                ref_force = max(np.linalg.norm(F_ext_step), np.linalg.norm(F_int_global), ZERO_TOL)
                err_force = np.linalg.norm(R) / ref_force
                error = max(err_disp, err_force)
                print(f"  Iteración {iteration+1:2d} | Err_dU: {err_disp:.4e} | Err_R: {err_force:.4e}")
                
                if error < self.tol:
                    print("  -> CONVERGENCIA ALCANZADA.")
                    self.assembler.commit_all_states()

                    U_current = U_iter
                    load_factor = next_load_factor
                    converged = True
                    
                    if self.adaptive and iteration < 4 and delta_lambda < (1.0 / self.num_steps):
                        delta_lambda = min(delta_lambda * 1.5, 1.0 / self.num_steps)
                        print(f"  -> Acelerando el próximo incremento a {delta_lambda:.4f}")
                        
                    if step_callback:
                        step_callback(step, U_current, load_factor)
                        
                    break
                    
            if not converged:
                if self.adaptive:
                    delta_lambda /= 2.0
                    print(f"  -> NO CONVERGIÓ. Bisección: reduciendo incremento a {delta_lambda:.4f}")
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
    def __init__(self, assembler, tol=CONVERGENCE_TOL, max_iter=20, max_lambda=1.0, initial_dl=0.1, max_steps=100, penalty_value=1e15,
                 dl_grow_factor=1.5, dl_max_factor=5.0, dl_shrink_factor=0.6,
                 dl_grow_iter_threshold=4, dl_shrink_iter_threshold=8):
        self.assembler = assembler
        self.tol = tol
        self.max_iter = max_iter
        self.max_lambda = max_lambda
        self.dl = initial_dl
        self.max_steps = max_steps
        self.penalty = penalty_value
        # Factores de auto-ajuste de la longitud de arco:
        #   Si converge en < dl_grow_iter_threshold iter → ampliar dl × dl_grow_factor (max: initial_dl × dl_max_factor)
        #   Si converge en > dl_shrink_iter_threshold iter → reducir dl × dl_shrink_factor
        self.dl_grow_factor = dl_grow_factor
        self.dl_max_factor = dl_max_factor
        self.dl_shrink_factor = dl_shrink_factor
        self.dl_grow_iter_threshold = dl_grow_iter_threshold
        self.dl_shrink_iter_threshold = dl_shrink_iter_threshold

    def solve(self, F_ext_ref: np.ndarray, step_callback=None) -> np.ndarray:
        domain = self.assembler.domain
        ndof = domain.total_dofs
        U_current = np.zeros(ndof)
        
        lambda_curr = 0.0
        step = 0
        dl = self.dl
        
        delta_U_step = np.zeros(ndof)  # Historial del incremento del paso para guiar el arco
        
        print("\n--- INICIANDO SOLVER NO LINEAL (MÉTODO ARC-LENGTH) ---")
        
        while lambda_curr < self.max_lambda and step < self.max_steps:
            step += 1
            print(f"\n[PASO {step}] Longitud de Arco (dl): {dl:.4e}")
            
            U_iter = U_current.copy()
            lambda_iter = lambda_curr
            converged = False
            
            # --- 1. PREDICTOR ---
            K_global, F_int_global = self.assembler.assemble_non_linear_system(U_iter)
            
            K_t = K_global.copy()
            F_t = F_ext_ref.copy()
            
            K_t, F_t = self.assembler.apply_bcs_to_system(K_t, F_t, penalty_value=self.penalty)
            
            try:
                du_t = spla.spsolve(K_t, F_t)
            except RuntimeError:
                print("  -> Error: Matriz singular en predictor. Bisección de dl...")
                dl /= 2.0
                continue
                
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

                K_t = K_global.copy()
                F_t = F_ext_ref.copy()

                K_t, F_t = self.assembler.apply_bcs_to_system(K_t, F_t, penalty_value=self.penalty)
                K_global, R = self.assembler.apply_bcs_to_system(K_global, R, penalty_value=self.penalty, load_factor=lambda_iter, U_current=U_iter)

                try:
                    du_R = spla.spsolve(K_global, R)
                    du_t = spla.spsolve(K_t, F_t)
                except RuntimeError:
                    print("  -> Error: Matriz Singular en corrector.")
                    break

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
                        print("  -> Raíces imaginarias. La solución diverge del arco.")
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
                
                err_disp = np.linalg.norm(dU_update) / (np.linalg.norm(U_iter) + ZERO_TOL)
                ref_force = max(
                    np.linalg.norm(F_ext_ref) * abs(lambda_iter),
                    np.linalg.norm(F_int_global),
                    ZERO_TOL,
                )
                err_force = np.linalg.norm(R) / ref_force
                error = max(err_disp, err_force)
                print(f"  Iter. {iteration+1:2d} | lam={lambda_iter:.4f} | Err_dU: {err_disp:.4e} | Err_R: {err_force:.4e}")
                
                if error < self.tol:
                    print(f"  -> CONVERGENCIA. (Lambda alcanzado: {lambda_iter:.4f})")
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
                print(f"  -> Bisección: reduciendo longitud de arco a {dl:.4e}")
                if dl < 1e-6 * self.dl:
                    raise RuntimeError("Arc-Length fracasó irreparablemente.")
                    
        return U_current
