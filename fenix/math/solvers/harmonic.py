# fenix_fem/fenix/math/solvers/harmonic.py
"""``HarmonicSolver`` — respuesta forzada armónica en el dominio de la
frecuencia (ADR 0009 fase 6).

Resuelve el sistema lineal complejo

.. math:: (-\\omega^2 \\mathbf{M} + i\\omega\\mathbf{C} + \\mathbf{K})\\,\\hat{\\mathbf{u}}(\\omega) = \\hat{\\mathbf{F}}

para un barrido de frecuencias ``ω``. La carga armónica ``F(t) = Re{F̂·e^{iωt}}``
y la respuesta estacionaria ``u(t) = Re{û(ω)·e^{iωt}}`` quedan codificadas
en la amplitud compleja ``û``: módulo = amplitud real, argumento = desfase.

Hipótesis:

- Sistema lineal: ``K`` y ``M`` constantes, ``C = α·M + β·K`` Rayleigh
  (mismo contrato que :class:`NewmarkSolver`).
- Régimen estacionario — el transitorio inicial no se modela (no
  asumiciones sobre condiciones iniciales; la solución dependería sólo
  de la excitación).
- Aritmética compleja completa: ``û``, ``F̂``, y la matriz efectiva
  ``Z(ω) = -ω²M + iωC + K`` son complejas. Cada ω requiere una
  factorización LU compleja distinta (la matriz cambia con ω; no hay
  caché de factorización reutilizable).

Apropiado para:

- Funciones de transferencia (FRFs) entrada→salida.
- Análisis sísmico en el dominio de la frecuencia.
- Diagnóstico de resonancias antes de un análisis transitorio caro.

Referencias:

- Bathe, K.-J. (2014). *Finite Element Procedures*, §9.5 (análisis en el
  dominio de la frecuencia).
- Clough, R. W., Penzien, J. (1993). *Dynamics of Structures*, §4.6.
- Chopra, A. K. (2017). *Dynamics of Structures*, cap. 4 (respuesta
  forzada armónica de sistemas SDOF y MDOF).
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from fenix.math.damping import rayleigh_from_modes
from fenix.math.solvers._shared import _log
from fenix.registry import SolverRegistry
from fenix.results import HarmonicResult


@SolverRegistry.register
class HarmonicSolver:
    """Análisis de respuesta forzada armónica (ADR 0009 fase 6).

    Parameters
    ----------
    assembler : Assembler
        Vinculación al modelo.
    omega : array-like or None, optional
        Vector explícito de frecuencias angulares (rad/s) a evaluar. Si
        se proporciona, ``omega_min``, ``omega_max``, ``n_omega`` y
        ``scale`` se ignoran. Forma ``(n_omega,)``.
    omega_min, omega_max : float, optional
        Límites del barrido en rad/s (incluyentes). Requeridos si
        ``omega`` no se pasa.
    n_omega : int, default 100
        Número de puntos del barrido.
    scale : {"linear", "log"}, default "linear"
        Espaciado del barrido. ``"log"`` usa ``np.geomspace`` —
        requerido por la convención de FRF en escala dB sobre rangos
        amplios. ``"linear"`` con ``np.linspace``.
    F_amplitude : np.ndarray or None
        Amplitud compleja de la carga, shape ``(ndof,)``. Constante a lo
        largo del barrido. Si ``None``, ceros (caso degenerado: ``û = 0``
        en todas las frecuencias).
    rayleigh : dict or None, default None
        Amortiguamiento Rayleigh ``C = α·M + β·K``. Mismo contrato que
        :class:`NewmarkSolver`. Sin amortiguamiento, los picos de
        resonancia son singulares teóricamente; numéricamente la
        factorización LU se vuelve ill-conditioned cerca de ω_n.
    lumping : str, default "consistent"
        Discretización de masa. ``"consistent"`` (default) o ``"lumped"``
        (ADR 0009 fase 2).
    linear_algebra : str, default "auto"
        Backend; aceptado por compatibilidad de firma (``scipy.sparse.linalg.spsolve``
        es la única ruta disponible para matrices complejas en esta
        versión inicial).
    """

    PIPELINE_KIND = "harmonic"

    def __init__(
        self,
        assembler,
        *,
        omega: np.ndarray | None = None,
        omega_min: float | None = None,
        omega_max: float | None = None,
        n_omega: int = 100,
        scale: str = "linear",
        F_amplitude: np.ndarray | None = None,
        rayleigh: dict | None = None,
        lumping: str = "consistent",
        linear_algebra: str = "auto",  # noqa: ARG002 — compat
    ):
        if omega is not None:
            omega_arr = np.asarray(omega, dtype=float).reshape(-1)
            if omega_arr.size == 0:
                raise ValueError("HarmonicSolver: omega vacío.")
            if np.any(omega_arr <= 0.0):
                raise ValueError(
                    "HarmonicSolver: todas las frecuencias deben ser "
                    "positivas (ω = 0 ⇒ análisis estático)."
                )
        else:
            if omega_min is None or omega_max is None:
                raise ValueError(
                    "HarmonicSolver: pasar `omega` o "
                    "(`omega_min`, `omega_max`, `n_omega`)."
                )
            if omega_min <= 0.0 or omega_max <= omega_min:
                raise ValueError(
                    f"HarmonicSolver: rango inválido "
                    f"omega_min={omega_min}, omega_max={omega_max}. "
                    f"Se requiere 0 < omega_min < omega_max."
                )
            if n_omega < 1:
                raise ValueError(f"HarmonicSolver: n_omega={n_omega} < 1.")
            if scale == "linear":
                omega_arr = np.linspace(omega_min, omega_max, n_omega)
            elif scale == "log":
                omega_arr = np.geomspace(omega_min, omega_max, n_omega)
            else:
                raise ValueError(
                    f"HarmonicSolver: scale='{scale}' no soportado. "
                    f"Admitidos: 'linear', 'log'."
                )

        self.assembler = assembler
        self.omega = omega_arr
        self.F_amplitude = F_amplitude
        self.rayleigh_cfg = rayleigh
        self.lumping = str(lumping)

    @staticmethod
    def _resolve_rayleigh(cfg: dict | None) -> tuple[float, float]:
        """Mismo contrato que :meth:`NewmarkSolver._resolve_rayleigh`."""
        if cfg is None:
            return 0.0, 0.0
        if not isinstance(cfg, dict):
            raise ValueError(
                f"HarmonicSolver.rayleigh: esperado dict o None, "
                f"recibido {type(cfg).__name__}."
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
            "HarmonicSolver.rayleigh: dict admitido es {'alpha', 'beta'} "
            "o {'xi1','omega1','xi2','omega2'}; recibido "
            + repr(set(cfg.keys()))
        )

    def solve(self) -> HarmonicResult:
        _log.info("--- INICIANDO SOLVER HARMONIC ---")

        # Ensamblar K y M (lineal). Para la respuesta armónica K es la
        # rigidez en u=0 (tangente inicial); el análisis es lineal.
        self.assembler.assemble_system()
        K = self.assembler.K_global
        M = self.assembler.assemble_mass_matrix(lumping=self.lumping)

        alpha_r, beta_r = self._resolve_rayleigh(self.rayleigh_cfg)
        C = alpha_r * M + beta_r * K

        # Reducción por Dirichlet. û en DOFs prescritos es cero (el apoyo
        # no oscila bajo la excitación armónica del problema homogéneo
        # asociado).
        cs = self.assembler.constraint_set
        T, _g = cs.build(self.assembler.ndof)

        K_red = (T.T @ K @ T).tocsr()
        M_red = (T.T @ M @ T).tocsr()
        C_red = (T.T @ C @ T).tocsr() if (alpha_r != 0.0 or beta_r != 0.0) \
                else sp.csr_matrix(K_red.shape, dtype=K_red.dtype)

        ndof = self.assembler.ndof
        if self.F_amplitude is None:
            F_amplitude = np.zeros(ndof, dtype=complex)
        else:
            F_amplitude = np.asarray(self.F_amplitude).reshape(ndof).astype(complex)
        F_red = T.T @ F_amplitude

        # Convertir a complejo desde el inicio — Z(ω) es siempre complejo
        # por el término iωC. Si C = 0 y ω real, la solución es real
        # numéricamente; mantenemos el tipo complejo por uniformidad.
        K_red_c = K_red.astype(complex)
        M_red_c = M_red.astype(complex)
        C_red_c = C_red.astype(complex)

        n_omega = self.omega.size
        u_complex = np.zeros((ndof, n_omega), dtype=complex)

        # Barrido. Cada ω requiere factorizar Z(ω) = -ω²M + iωC + K.
        # No hay cache: la matriz cambia con ω.
        for k, w in enumerate(self.omega):
            Z = (-w * w) * M_red_c + (1j * w) * C_red_c + K_red_c
            u_red_k = spla.spsolve(Z.tocsc(), F_red)
            u_complex[:, k] = T @ u_red_k

        _log.info(
            f"  -> {n_omega} frecuencias evaluadas en "
            f"[{self.omega[0]:.4e}, {self.omega[-1]:.4e}] rad/s. "
            f"Rayleigh: α={alpha_r:.4e}, β={beta_r:.4e}."
        )

        return HarmonicResult(
            omega=self.omega.copy(),
            u_complex=u_complex,
            F_complex=F_amplitude.copy(),
            alpha_rayleigh=alpha_r,
            beta_rayleigh=beta_r,
        )
