# solidum_fem/solidum/math/convergence.py
"""Política de convergencia para solvers no lineales (ADR 0007).

Implementa el patrón ``atol + rtol·escala(estado)`` separado por criterio
(fuerza y desplazamiento). Es el único punto del proyecto donde vive la
fórmula del criterio de convergencia; los solvers no lineales la consumen
vía ``ConvergenceCriterion.evaluate(...)``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from solidum.constants import (
    CONVERGENCE_RTOL_FORCE,
    CONVERGENCE_RTOL_DISP,
    CONVERGENCE_ATOL_FORCE_FACTOR,
    CONVERGENCE_ATOL_DISP_FACTOR,
)


@dataclass(frozen=True)
class ConvergenceState:
    """Resultado de una evaluación del criterio en una iteración concreta.

    Attributes
    ----------
    converged : bool
        ``True`` si ambos criterios (fuerza y desplazamiento) se satisfacen
        simultáneamente (semántica AND, ADR 0007).
    ratio_force, ratio_disp : float
        Cocientes ``‖residuo‖/tol_force`` y ``‖δU‖/tol_disp``. Valor ``< 1``
        indica criterio cumplido. Útil para logging y diagnóstico de cuál
        de los dos criterios manda en cada iteración.
    err_force, err_disp : float
        Magnitudes crudas (no normalizadas): ``‖residuo‖`` y ``‖δU‖``.
    tol_force, tol_disp : float
        Tolerancias efectivas aplicadas, ``atol + rtol·escala``.
    """
    converged: bool
    ratio_force: float
    ratio_disp: float
    err_force: float
    err_disp: float
    tol_force: float
    tol_disp: float


class ConvergenceCriterion:
    """Criterio de convergencia ``atol + rtol·escala`` para Newton-Raphson y arc-length.

    Encapsula configuración (cuatro tolerancias) y estado calibrado (``atol``
    absolutos derivados del primer ensamblaje del problema). La calibración
    una sola vez por corrida garantiza invariancia completa bajo cambio de
    unidades: el usuario no necesita ajustar nada al pasar de N/m a kN/mm o
    a MPa/mm.

    Parameters
    ----------
    rtol_force, rtol_disp : float
        Tolerancias relativas dimensionalmente puras. Banda sobre la escala
        del estado corriente que se considera convergencia. Default ``1e-5``
        para ambos.
    atol_force_factor, atol_disp_factor : float
        Factores adimensionales del término absoluto. La tolerancia absoluta
        efectiva se autoderiva como ``factor · escala_inicial`` durante la
        calibración. Default ``1e-9`` para ambos — muy por debajo del rtol
        en problemas con escalas físicas normales, actúan solo como piso
        cuando la escala corriente colapsa.

    Notes
    -----
    El criterio se aplica como:

    .. code-block:: text

        ‖R[free_dofs]‖ ≤ atol_force + rtol_force · max(‖F_ext‖, ‖F_int‖)
        ‖δU‖           ≤ atol_disp  + rtol_disp  · ‖U_iter‖

    con semántica AND (ambos criterios simultáneamente). Ver ADR 0007 para
    la justificación de los defaults y de la elección AND vs OR.
    """

    def __init__(self,
                 rtol_force: float = CONVERGENCE_RTOL_FORCE,
                 rtol_disp: float = CONVERGENCE_RTOL_DISP,
                 atol_force_factor: float = CONVERGENCE_ATOL_FORCE_FACTOR,
                 atol_disp_factor: float = CONVERGENCE_ATOL_DISP_FACTOR):
        self.rtol_force = float(rtol_force)
        self.rtol_disp = float(rtol_disp)
        self.atol_force_factor = float(atol_force_factor)
        self.atol_disp_factor = float(atol_disp_factor)
        self._atol_force: float | None = None
        self._atol_disp: float | None = None

    @property
    def is_calibrated(self) -> bool:
        return self._atol_force is not None

    @property
    def atol_force(self) -> float:
        """Tolerancia absoluta efectiva de fuerza tras la calibración."""
        if self._atol_force is None:
            raise RuntimeError(
                "ConvergenceCriterion no calibrado: invocar calibrate() "
                "antes de leer atol_force."
            )
        return self._atol_force

    @property
    def atol_disp(self) -> float:
        """Tolerancia absoluta efectiva de desplazamiento tras la calibración."""
        if self._atol_disp is None:
            raise RuntimeError(
                "ConvergenceCriterion no calibrado: invocar calibrate() "
                "antes de leer atol_disp."
            )
        return self._atol_disp

    def calibrate(self, force_scale: float, disp_scale: float) -> None:
        """Fija los ``atol`` efectivos a partir de las escalas del problema.

        Pensado para invocarse **una vez** al inicio de ``solve()``, tras
        disponer del primer ensamblaje. Las escalas las computa el solver
        a partir de su propio dato:

        - ``force_scale``: típicamente ``max(‖F_ext_global‖, ‖F_int_initial‖)``.
          Con control en desplazamiento puro y sin prestress, ``F_int_initial``
          es cero y ``F_ext_global`` también; en ese caso pasar ``1.0`` como
          fallback (el rtol dominará en cuanto haya carga real).
        - ``disp_scale``: estimador de orden 1 del desplazamiento esperado,
          típicamente ``force_scale / max(|diag(K_initial)|)``.

        Una segunda llamada sobreescribe los atol (no es error). Útil si el
        solver decide recalibrar tras un cambio brusco de configuración del
        problema, aunque no es práctica recomendada.
        """
        self._atol_force = self.atol_force_factor * float(force_scale)
        self._atol_disp = self.atol_disp_factor * float(disp_scale)

    def evaluate(self,
                 residual_norm: float,
                 ref_force: float,
                 delta_u_norm: float,
                 u_norm: float) -> ConvergenceState:
        """Aplica el criterio ``atol + rtol·escala`` y devuelve el estado.

        Parameters
        ----------
        residual_norm : float
            ``‖R[free_dofs]‖`` — norma L2 del residuo en DOFs libres. En DOFs
            prescritos el residuo contiene la reacción física, no un fallo
            de equilibrio.
        ref_force : float
            Escala corriente de fuerza, típicamente
            ``max(‖F_ext_step‖, ‖F_int‖)``. Se reevalúa cada iteración para
            ser adaptativa al estado.
        delta_u_norm : float
            ``‖δU‖`` — norma L2 del incremento de la iteración.
        u_norm : float
            ``‖U_iter‖`` — norma L2 del desplazamiento corriente.

        Returns
        -------
        ConvergenceState
            Bool de convergencia + ratios + tolerancias efectivas para
            logging.
        """
        if not self.is_calibrated:
            raise RuntimeError(
                "ConvergenceCriterion.evaluate() invocado antes de calibrar. "
                "Llamar a calibrate(force_scale, disp_scale) al inicio de solve()."
            )

        tol_force = self._atol_force + self.rtol_force * float(ref_force)
        tol_disp = self._atol_disp + self.rtol_disp * float(u_norm)

        # tol_*  > 0 está garantizado por construcción: rtol·escala ≥ 0 y
        # atol > 0 si force_scale/disp_scale > 0; el solver pasa fallback 1.0
        # cuando las escalas reales colapsan.
        return ConvergenceState(
            converged=(residual_norm <= tol_force) and (delta_u_norm <= tol_disp),
            ratio_force=residual_norm / tol_force,
            ratio_disp=delta_u_norm / tol_disp,
            err_force=residual_norm,
            err_disp=delta_u_norm,
            tol_force=tol_force,
            tol_disp=tol_disp,
        )


def make_convergence_from_config(cfg: dict | None) -> ConvergenceCriterion:
    """Construye un ``ConvergenceCriterion`` desde un dict YAML.

    El dict acepta cualquier subconjunto de ``rtol_force``, ``rtol_disp``,
    ``atol_force_factor``, ``atol_disp_factor``. Claves desconocidas lanzan
    ``ValueError`` para detectar typos pronto.
    """
    if cfg is None:
        return ConvergenceCriterion()
    allowed = {"rtol_force", "rtol_disp", "atol_force_factor", "atol_disp_factor"}
    unknown = set(cfg.keys()) - allowed
    if unknown:
        raise ValueError(
            f"convergence: claves desconocidas {sorted(unknown)}. "
            f"Permitidas: {sorted(allowed)}."
        )
    return ConvergenceCriterion(**cfg)


def stiffness_diag_scale(K) -> float:
    """Escala característica de rigidez como ``max(|diag(K)|)``.

    Más barata que la norma matricial completa y representativa para mallas
    estructurales típicas. Funciona tanto con arrays densos como con sparse
    de scipy (todos exponen ``.diagonal()``). Garantiza ``> 0`` para evitar
    divisiones por cero aguas arriba.
    """
    diag = np.asarray(K.diagonal()).ravel()
    val = float(np.abs(diag).max()) if diag.size else 0.0
    return max(val, 1.0e-300)
