# fenix_fem/fenix/math/solvers/modal.py
"""``ModalSolver`` — problema de valores característicos ``K·φ = ω²·M·φ``
(ADR 0009 fase 1).
"""
import numpy as np

from fenix.math.linalg import EigenSolver
from fenix.math.solvers._shared import _log
from fenix.registry import SolverRegistry
from fenix.results import ModalResult


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
                1.0 / frequencies_hz,
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
