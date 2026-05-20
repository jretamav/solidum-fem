# solidum_fem/solidum/math/solvers/linear.py
"""``LinearSolver`` — un paso, sistema ``K·U = F`` con Dirichlet por eliminación.

La factorización se cachea entre llamadas a :meth:`solve` (ADR 0003): la
primera llamada ensambla ``K``, reduce a DOFs libres (ADR 0004) y factoriza;
las siguientes reúsan el factor para resolver con un ``F`` distinto. El
patrón es el mismo que aplica :class:`NewmarkSolver` a ``A_eff`` y permite
barridos baratos sobre cargas (multiplicidad de cargas, FRF lineal,
parametric studies) sin pagar la factorización en cada llamada.

Si el usuario modifica el modelo entre llamadas (nuevos elementos, nuevas
BCs, cambios de materiales) debe invocar :meth:`invalidate_cache` para
forzar reensamblaje y refactorización.
"""
import numpy as np

from solidum.math.linalg import LUSolver, StiffnessProperties, select_solver
from solidum.math.solvers._shared import (
    CholeskyNotPositiveDefiniteError,
    _log,
    domain_is_symmetric,
)
from solidum.registry import SolverRegistry


@SolverRegistry.register
class LinearSolver:
    """Solucionador de sistemas algebraicos lineales en un solo paso."""

    PIPELINE_KIND = "static"

    def __init__(self, assembler, linear_algebra: str = "auto"):
        self.assembler = assembler
        self.linear_algebra = linear_algebra
        # Cache lazy — se llena en la 1ª llamada a `solve`.
        self._factor = None
        self._T = None
        self._g_full = None
        self._F_dir = None
        self._n_free = None

    def invalidate_cache(self) -> None:
        """Descarta la factorización cacheada.

        Necesario sólo si el modelo cambia entre llamadas a :meth:`solve`
        (nuevos elementos, nuevas BCs, cambios de materiales). Para
        barrido de cargas con el mismo modelo no hace falta — la cache
        sobrevive entre invocaciones de ``solve``.
        """
        self._factor = None
        self._T = None
        self._g_full = None
        self._F_dir = None
        self._n_free = None

    def _build_cache(self) -> None:
        """Ensambla ``K``, reduce a libres, factoriza y guarda el factor.

        Llamado lazy desde la 1ª :meth:`solve`. Maneja el fallback
        Cholesky → LU (ADR 0003 §5) si la matriz no es PD efectiva.
        """
        _log.info("--- LINEARSOLVER · ensamblando y factorizando (1ª llamada) ---")
        self.assembler.assemble_system()
        K_global = self.assembler.K_global.copy()

        # Reducir con F dummy: sólo necesitamos K_red, T y g_full aquí.
        K_red, _F_red_dummy, T, g_full = self.assembler.reduce(
            K_global, np.zeros(self.assembler.ndof),
        )

        # F_dir = T.T · K · g (contribución de Dirichlet no homogéneo).
        # Constante mientras K y g lo sean — se cachea junto al factor.
        F_dir = T.T @ (K_global @ g_full)

        # En estática lineal, si el dominio es simétrico (todos los
        # materiales con tangente simétrica) asumimos PD; si no, LU. El
        # fallback automático Cholesky→LU sigue cubriendo los casos donde
        # la simetría no implique PD (degenerados, casi-singulares).
        is_sym = domain_is_symmetric(self.assembler.domain)
        props = StiffnessProperties(
            is_symmetric=is_sym,
            is_positive_definite=is_sym,
            size=K_red.shape[0],
        )
        linalg = select_solver(props, override=self.linear_algebra)
        try:
            factor = linalg.factorize(K_red)
        except CholeskyNotPositiveDefiniteError:
            # Fallback automático SPD→LU (ADR 0003 §5).
            _log.warning("Cholesky reportó no-positividad. Degradando a LU.")
            factor = LUSolver().factorize(K_red)

        self._factor = factor
        self._T = T
        self._g_full = g_full
        self._F_dir = F_dir
        self._n_free = K_red.shape[0]

    def solve(self, F_ext_global: np.ndarray) -> np.ndarray:
        if self._factor is None:
            self._build_cache()
        else:
            _log.info("--- LINEARSOLVER · reusando factorización cacheada ---")

        # F_red = T.T · (F − K · g) = T.T · F − F_dir, con F_dir cacheado.
        F_red = self._T.T @ F_ext_global - self._F_dir
        u_red = self._factor.solve(F_red)
        U = self.assembler.expand(u_red, self._T, self._g_full)
        _log.info("  -> CONVERGENCIA ALCANZADA (1 Iteración).")
        return U
