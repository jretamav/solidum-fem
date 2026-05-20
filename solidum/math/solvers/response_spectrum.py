# fenix_fem/solidum/math/solvers/response_spectrum.py
"""``ResponseSpectrumSolver`` — análisis de respuesta espectral por
combinación modal (ADR 0009 fase 7).

Orquesta:

1. Análisis modal interno vía :class:`ModalSolver` para obtener los
   primeros ``n_modes`` modos y frecuencias del sistema.
2. Cálculo de los factores de participación modal ``Γ_n = φ_nᵀ·M·r``
   donde ``r`` es el vector de excitación rígida (dirección sísmica).
3. Evaluación del espectro de respuesta ``S_d(ω_n)`` para cada modo.
4. Combinación modal: SRSS (default) o CQC (requiere amortiguamiento).
5. Devuelve :class:`ResponseSpectrumResult` con la respuesta máxima
   envolvente, la contribución modal individual, factores de
   participación y masas efectivas.

Limitación inherente del método: la respuesta es **envolvente máxima**,
no historia temporal. Si necesitas la historia completa, usa
transitorio (Newmark, HHT) con un acelerograma — el espectro de
respuesta es un descriptor estadístico de aceleraciones del suelo
agregadas sobre múltiples sismos.

Convención del espectro: ``spectrum_fn(omega)`` devuelve el
**desplazamiento espectral** ``S_d(ω)`` en las unidades del modelo. El
usuario convierte de ``S_a`` (aceleración espectral, lo más común en
normativas) a ``S_d`` vía ``S_d = S_a / ω²`` si su input es ``S_a``;
:func:`spectrum_from_sa` ofrece este wrapper.

Referencias:

- Chopra, A. K. (2017). *Dynamics of Structures*, cap. 13 (análisis
  espectral de sistemas MDOF).
- Wilson, E. L. (2002). *Three-Dimensional Static and Dynamic Analysis
  of Structures*, cap. 15 (CQC).
- Der Kiureghian, A. (1980). "Structural Response to Stationary
  Excitation". *J. Eng. Mech. Div. ASCE*, 106(EM6), 1195-1213.
"""
from __future__ import annotations

from typing import Callable

import numpy as np

from solidum.math.modal_response import (
    response_spectrum_cqc,
    response_spectrum_srss,
)
from solidum.math.solvers._shared import _log
from solidum.math.solvers.modal import ModalSolver
from solidum.registry import SolverRegistry
from solidum.results import ResponseSpectrumResult


def spectrum_from_sa(spectrum_sa_fn: Callable[[float], float]) -> Callable[[float], float]:
    """Convierte un espectro de aceleración ``S_a(ω)`` a desplazamiento
    ``S_d(ω) = S_a(ω) / ω²``.

    Aceleración espectral es lo más común en espectros normativos (NCh,
    ASCE 7, Eurocódigo 8, NTC); el solver internamente trabaja con
    desplazamiento espectral porque es lo que aparece en la fórmula
    ``u_n^max = Γ_n·φ_n·S_d(ω_n)``.

    Parameters
    ----------
    spectrum_sa_fn : callable
        ``spectrum_sa_fn(omega) -> S_a`` (aceleración espectral).

    Returns
    -------
    callable
        ``S_d(ω) = S_a(ω) / ω²`` listo para pasar a
        :class:`ResponseSpectrumSolver`.
    """
    def Sd(omega: float) -> float:
        if omega <= 0.0:
            return 0.0
        return float(spectrum_sa_fn(omega)) / (omega * omega)
    return Sd


def spectrum_tabulated(
    periods: np.ndarray,
    values: np.ndarray,
    *,
    kind: str = "Sd",
) -> Callable[[float], float]:
    """Construye un espectro como callable a partir de una tabla
    ``(periodo, valor)``.

    Interpolación lineal en el periodo. Fuera del rango tabulado,
    devuelve los extremos (clamp) — política conservadora; el usuario
    extiende la tabla si necesita extrapolación física.

    Parameters
    ----------
    periods : array-like
        Periodos en segundos, ordenados ascendentemente.
    values : array-like
        Valores del espectro en cada periodo.
    kind : {"Sd", "Sa"}, default "Sd"
        Si ``"Sa"``, se convierte internamente a ``S_d`` vía
        :func:`spectrum_from_sa`.

    Returns
    -------
    callable
        ``S_d(ω) -> float``.
    """
    periods = np.asarray(periods, dtype=float).ravel()
    values = np.asarray(values, dtype=float).ravel()
    if periods.size != values.size:
        raise ValueError(
            f"spectrum_tabulated: periods.size ({periods.size}) "
            f"!= values.size ({values.size})."
        )
    if np.any(np.diff(periods) <= 0.0):
        raise ValueError(
            "spectrum_tabulated: periods deben ser estrictamente crecientes."
        )
    if kind not in ("Sd", "Sa"):
        raise ValueError(
            f"spectrum_tabulated: kind='{kind}' no admitido. "
            f"Valores válidos: 'Sd', 'Sa'."
        )

    def interp_T(T: float) -> float:
        if T <= periods[0]:
            return float(values[0])
        if T >= periods[-1]:
            return float(values[-1])
        return float(np.interp(T, periods, values))

    if kind == "Sd":
        def Sd(omega: float) -> float:
            if omega <= 0.0:
                return 0.0
            T = 2.0 * np.pi / omega
            return interp_T(T)
        return Sd

    # kind == "Sa": convertir a Sd dentro del propio callable.
    def Sd_from_Sa(omega: float) -> float:
        if omega <= 0.0:
            return 0.0
        T = 2.0 * np.pi / omega
        Sa = interp_T(T)
        return Sa / (omega * omega)

    return Sd_from_Sa


@SolverRegistry.register
class ResponseSpectrumSolver:
    """Análisis de respuesta espectral por combinación modal (ADR 0009
    fase 7).

    Parameters
    ----------
    assembler : Assembler
        Vinculación al modelo.
    n_modes : int
        Número de modos a calcular y combinar. Recomendado lo suficiente
        para alcanzar masa efectiva acumulada ≥ 90% en la dirección de
        excitación (verificable a posteriori con
        :meth:`ResponseSpectrumResult.cumulative_effective_mass_ratio`).
    direction : array-like, shape (n_dof,) or dict
        Vector de excitación rígida ``r``. Puede pasarse como:

        - ``np.ndarray`` de shape ``(n_dof,)`` — vector directo.
        - ``dict`` con clave ``"dof_name"`` (e.g. ``"ux"``) — el solver
          construye el vector unitario en ese DOF en todos los nodos.
    spectrum : callable or dict
        Espectro de respuesta. Acepta:

        - ``callable`` ``spectrum_fn(omega) -> S_d`` directamente.
        - ``dict`` con forma ``{"type": "tabulated", "periods": [...],
          "values": [...], "kind": "Sd"|"Sa"}`` para construcción desde
          YAML (delega a :func:`spectrum_tabulated`).
        - ``dict`` con forma ``{"type": "constant_acceleration",
          "Sa": ...}`` para excitación aceleración constante.
    combination : {"SRSS", "CQC"}, default "SRSS"
        Método de combinación modal.
    damping : float, default 0.05
        Amortiguamiento modal ξ (mismo para todos los modos). Sólo
        usado por CQC para los coeficientes de correlación; SRSS lo
        ignora.
    sigma : float, default 0.0
        Shift para shift-invert en el cálculo modal interno. Mismo
        contrato que :class:`ModalSolver`.
    lumping : str, default "consistent"
        Pasado al ModalSolver interno.
    """

    PIPELINE_KIND = "spectrum"

    def __init__(
        self,
        assembler,
        *,
        n_modes: int,
        direction,
        spectrum,
        combination: str = "SRSS",
        damping: float = 0.05,
        sigma: float = 0.0,
        lumping: str = "consistent",
    ):
        if n_modes < 1:
            raise ValueError(
                f"ResponseSpectrumSolver: n_modes={n_modes} < 1."
            )
        if combination not in ("SRSS", "CQC"):
            raise ValueError(
                f"ResponseSpectrumSolver: combination='{combination}' "
                f"no admitido. Valores: 'SRSS', 'CQC'."
            )
        if damping < 0.0:
            raise ValueError(
                f"ResponseSpectrumSolver: damping={damping} negativo."
            )

        self.assembler = assembler
        self.n_modes = int(n_modes)
        self.direction_raw = direction
        self.spectrum_raw = spectrum
        self.combination = combination
        self.damping = float(damping)
        self.sigma = float(sigma)
        self.lumping = str(lumping)

    def _resolve_direction(self) -> np.ndarray:
        """Resuelve el vector de excitación rígida desde el formato de
        entrada (ndarray directo o dict ``{"dof_name": ...}``).
        """
        ndof = self.assembler.ndof
        raw = self.direction_raw
        if isinstance(raw, dict):
            dof_name = raw.get("dof_name")
            if dof_name is None:
                raise ValueError(
                    "ResponseSpectrumSolver.direction: dict requiere "
                    "clave 'dof_name' (e.g. 'ux'). Recibido: "
                    f"{list(raw.keys())}."
                )
            r = np.zeros(ndof)
            for node in self.assembler.domain.nodes.values():
                if dof_name in node.dofs:
                    r[node.dofs[dof_name]] = 1.0
            if np.linalg.norm(r) == 0.0:
                raise ValueError(
                    f"ResponseSpectrumSolver.direction: ningún nodo "
                    f"tiene DOF '{dof_name}'."
                )
            return r
        r = np.asarray(raw, dtype=float).ravel()
        if r.size != ndof:
            raise ValueError(
                f"ResponseSpectrumSolver.direction: shape {r.shape} "
                f"!= ndof={ndof}."
            )
        return r

    def _resolve_spectrum(self) -> Callable[[float], float]:
        """Resuelve el espectro desde el formato de entrada (callable o
        dict descriptivo).
        """
        raw = self.spectrum_raw
        if callable(raw):
            return raw
        if isinstance(raw, dict):
            stype = raw.get("type")
            if stype == "tabulated":
                return spectrum_tabulated(
                    periods=raw["periods"],
                    values=raw["values"],
                    kind=raw.get("kind", "Sd"),
                )
            if stype == "constant_acceleration":
                Sa_const = float(raw["Sa"])
                return spectrum_from_sa(lambda _w: Sa_const)
            raise ValueError(
                f"ResponseSpectrumSolver.spectrum: type='{stype}' "
                f"no admitido. Valores: 'tabulated', 'constant_acceleration'."
            )
        raise ValueError(
            f"ResponseSpectrumSolver.spectrum: tipo no admitido "
            f"({type(raw).__name__}). Debe ser callable o dict."
        )

    def solve(self) -> ResponseSpectrumResult:
        _log.info("--- INICIANDO SOLVER RESPONSE SPECTRUM ---")

        # Análisis modal interno.
        modal = ModalSolver(
            self.assembler,
            n_modes=self.n_modes,
            sigma=self.sigma,
            lumping=self.lumping,
        )
        modal_result = modal.solve()

        # Necesitamos M para los factores de participación y la
        # combinación. La masa ya está cacheada en el assembler tras
        # ModalSolver.solve().
        M = self.assembler.assemble_mass_matrix(lumping=self.lumping)

        direction = self._resolve_direction()
        spectrum_fn = self._resolve_spectrum()

        if self.combination == "SRSS":
            u_combined, u_per_mode, gamma = response_spectrum_srss(
                modal_result.modes,
                modal_result.frequencies_rad,
                M,
                direction,
                spectrum_fn,
            )
        else:  # CQC
            u_combined, u_per_mode, gamma = response_spectrum_cqc(
                modal_result.modes,
                modal_result.frequencies_rad,
                M,
                direction,
                spectrum_fn,
                damping=self.damping,
            )

        effective_masses = gamma * gamma

        _log.info(
            f"  -> {self.n_modes} modos combinados con {self.combination}. "
            f"ξ={self.damping:.3f}. "
            f"Masa efectiva acumulada: "
            f"{float(np.sum(effective_masses)):.3e}."
        )

        return ResponseSpectrumResult(
            u_combined=u_combined,
            u_per_mode=u_per_mode,
            frequencies_rad=modal_result.frequencies_rad.copy(),
            participation_factors=gamma,
            effective_masses=effective_masses,
            combination=self.combination,
            damping=self.damping if self.combination == "CQC" else 0.0,
            direction=direction.copy(),
        )
