# fenix_fem/fenix/math/solvers/linear.py
"""``LinearSolver`` — un paso, sistema ``K·U = F`` con Dirichlet por eliminación.
"""
import numpy as np

from fenix.math.linalg import LUSolver, StiffnessProperties, select_solver
from fenix.math.solvers._shared import (
    CholeskyNotPositiveDefiniteError,
    _log,
    domain_is_symmetric,
)
from fenix.registry import SolverRegistry


@SolverRegistry.register
class LinearSolver:
    """Solucionador de sistemas algebraicos lineales en un solo paso."""

    PIPELINE_KIND = "static"

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
            is_symmetric=domain_is_symmetric(self.assembler.domain),
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
