# solidum_fem/solidum/math/convergence.py
"""Política de convergencia para solvers no lineales (ADR 0007).

Implementa el patrón ``atol + rtol·escala(estado) + piso de redondeo``
separado por criterio (fuerza y desplazamiento). Es el único punto del
proyecto donde vive la fórmula del criterio de convergencia; los solvers no
lineales la consumen vía ``ConvergenceCriterion.evaluate(...)``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
from numba import njit

from solidum.constants import (
    CONVERGENCE_RTOL_FORCE,
    CONVERGENCE_RTOL_DISP,
    CONVERGENCE_ATOL_FORCE_FACTOR,
    CONVERGENCE_ATOL_DISP_FACTOR,
    CONVERGENCE_ROUNDOFF_FACTOR,
)

#: Épsilon de máquina de los cálculos en doble precisión.
_EPS = float(np.finfo(float).eps)


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
        de los dos criterios manda en cada iteración. Con tolerancia nula
        (problema sin fuerzas ni desplazamientos) vale 0 si la magnitud es
        nula e infinito si no.
    err_force, err_disp : float
        Magnitudes crudas (no normalizadas): ``‖residuo‖`` y ``‖δU‖``.
    tol_force, tol_disp : float
        Tolerancias efectivas aplicadas, ``atol + rtol·escala`` (más el piso
        de redondeo en la de fuerza).
    floor_force : float
        Parte de ``tol_force`` que aporta el piso de redondeo,
        ``roundoff_factor · roundoff_force``: si domina, el paso se acepta
        en el límite de la doble precisión.
    """
    converged: bool
    ratio_force: float
    ratio_disp: float
    err_force: float
    err_disp: float
    tol_force: float
    tol_disp: float
    floor_force: float


def _ratio(err: float, tol: float) -> float:
    """``err/tol``; con ``tol = 0``, 0 si ``err = 0`` e infinito si no."""
    if tol > 0.0:
        return err / tol
    return 0.0 if err == 0.0 else float("inf")


class ConvergenceCriterion:
    """Criterio de convergencia ``atol + rtol·escala`` para Newton-Raphson y arc-length.

    Encapsula configuración (cuatro tolerancias y el factor del piso de
    redondeo) y estado calibrado (``atol`` absolutos derivados del primer
    ensamblaje del problema). La calibración una sola vez por corrida y el
    piso de redondeo, calculado con el propio modelo, hacen el criterio
    invariante bajo cambio de unidades: el usuario no necesita ajustar nada
    al pasar de N/m a kN/mm o a MPa/mm.

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
    roundoff_factor : float
        Factor adimensional del piso de redondeo (ADR 0007, enmienda del
        2026-11-27): a la tolerancia de fuerza se suma
        ``roundoff_factor · roundoff_force``, con ``roundoff_force`` la escala
        del error de redondeo del residuo (:func:`roundoff_scale`). Un
        residuo por debajo no es un desequilibrio que el cálculo en doble
        precisión pueda resolver. Default ``10``, empírico: el ruido medido
        no pasa de 0.5 veces la escala en los elementos de cinemática
        lineal; no es una cota rigurosa.

    Notes
    -----
    El criterio se aplica como:

    .. code-block:: text

        ‖TᵀR‖ ≤ atol_force + rtol_force · max(‖F_ext‖, ‖F_int‖)
                + roundoff_factor · ε·‖|T|ᵀ·|K|·|U|‖
        ‖δU‖  ≤ atol_disp  + rtol_disp  · ‖U_iter‖

    con semántica AND (ambos criterios simultáneamente) y ``TᵀR`` el
    residuo reducido (ADR 0004), que con sólo condiciones de Dirichlet es
    ``R`` en los DOF libres. Ver ADR 0007 para la justificación de los
    defaults y de la elección AND vs OR.
    """

    def __init__(self,
                 rtol_force: float = CONVERGENCE_RTOL_FORCE,
                 rtol_disp: float = CONVERGENCE_RTOL_DISP,
                 atol_force_factor: float = CONVERGENCE_ATOL_FORCE_FACTOR,
                 atol_disp_factor: float = CONVERGENCE_ATOL_DISP_FACTOR,
                 roundoff_factor: float = CONVERGENCE_ROUNDOFF_FACTOR):
        values = {
            "rtol_force": rtol_force, "rtol_disp": rtol_disp,
            "atol_force_factor": atol_force_factor, "atol_disp_factor": atol_disp_factor,
            "roundoff_factor": roundoff_factor,
        }
        for name, value in values.items():
            _check_factor(name, value)
        self.rtol_force = float(rtol_force)
        self.rtol_disp = float(rtol_disp)
        self.atol_force_factor = float(atol_force_factor)
        self.atol_disp_factor = float(atol_disp_factor)
        self.roundoff_factor = float(roundoff_factor)
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
        a partir de su propio dato (``default_calibration_scales`` en
        ``solidum.math.solvers.corrector``):

        - ``force_scale``: ``max(‖F_ext_global‖, ‖F_int_initial‖)``. Con
          control en desplazamiento desde el reposo las dos son nulas y la
          escala también: el piso de fuerza lo da entonces sólo el término
          de redondeo, calculado con el modelo en cada iteración, y no una
          unidad de fuerza (enmienda del 2026-11-27).
        - ``disp_scale``: estimador del desplazamiento esperado,
          ``force_scale / max(|diag(K_initial)|)``, o la norma de los
          valores impuestos cuando no hay fuerzas.

        Las escalas nulas son válidas. Una segunda llamada sobreescribe los
        atol (no es error).
        """
        self._atol_force = self.atol_force_factor * float(force_scale)
        self._atol_disp = self.atol_disp_factor * float(disp_scale)

    def evaluate(self,
                 residual_norm: float,
                 ref_force: float,
                 delta_u_norm: float,
                 u_norm: float,
                 *,
                 roundoff_force: float) -> ConvergenceState:
        """Aplica el criterio ``atol + rtol·escala + piso`` y devuelve el estado.

        Parameters
        ----------
        residual_norm : float
            ``‖TᵀR‖`` — norma L2 del residuo reducido (ADR 0004). Excluye los
            DOF prescritos, donde el residuo contiene la reacción física y no
            un fallo de equilibrio, y suma a sus maestros el residuo de los
            esclavos de restricciones multipunto.
        ref_force : float
            Escala corriente de fuerza, típicamente
            ``max(‖F_ext_step‖, ‖F_int‖)``. Se reevalúa cada iteración para
            ser adaptativa al estado.
        delta_u_norm : float
            ``‖δU‖`` — norma L2 del incremento de la iteración.
        u_norm : float
            ``‖U_iter‖`` — norma L2 del desplazamiento corriente.
        roundoff_force : float
            Escala del error de redondeo del residuo en el iterado,
            :func:`roundoff_scale` (en unidades de fuerza). Obligatoria: un
            solver que la omitiera volvería a un piso que depende de las
            unidades.

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

        floor_force = self.roundoff_factor * float(roundoff_force)
        tol_force = self._atol_force + self.rtol_force * float(ref_force) + floor_force
        tol_disp = self._atol_disp + self.rtol_disp * float(u_norm)

        # Una tolerancia nula sólo ocurre en un problema sin fuerzas ni
        # desplazamientos (todo cero): converge si la magnitud es nula.
        return ConvergenceState(
            converged=(residual_norm <= tol_force) and (delta_u_norm <= tol_disp),
            ratio_force=_ratio(residual_norm, tol_force),
            ratio_disp=_ratio(delta_u_norm, tol_disp),
            err_force=residual_norm,
            err_disp=delta_u_norm,
            tol_force=tol_force,
            tol_disp=tol_disp,
            floor_force=floor_force,
        )


def make_convergence_from_config(cfg: dict | None) -> ConvergenceCriterion:
    """Construye un ``ConvergenceCriterion`` desde un dict YAML.

    El dict acepta cualquier subconjunto de ``rtol_force``, ``rtol_disp``,
    ``atol_force_factor``, ``atol_disp_factor`` y ``roundoff_factor``. Claves
    desconocidas lanzan ``ValueError`` para detectar typos pronto.
    """
    if cfg is None:
        return ConvergenceCriterion()
    allowed = {"rtol_force", "rtol_disp", "atol_force_factor", "atol_disp_factor",
               "roundoff_factor"}
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


@njit(cache=True)
def _abs_csr_matvec_add(indptr, indices, data, x_abs, out):  # pragma: no cover - numba
    """``out += |A|·x_abs`` para ``A`` en CSR, sin copiar ``|A|``."""
    for i in range(indptr.size - 1):
        s = 0.0
        for k in range(indptr[i], indptr[i + 1]):
            s += abs(data[k]) * x_abs[indices[k]]
        out[i] += s


def roundoff_vector(*terms, T=None) -> np.ndarray | None:
    """Escala, componente a componente, del error de redondeo del residuo:
    ``ε·|T|ᵀ·Σ|A|·|x|`` sobre los términos ``(A, x)`` que lo forman.

    ``(K, U)`` en estática; ``(K, u)``, ``(C, u̇)`` y ``(M, ü)`` en dinámica.
    Calcular ``A·x`` en doble precisión deja en cada componente un error del
    orden de ``ε·(|A|·|x|)_i``, que ninguna iteración de Newton reduce (ADR
    0007, enmienda del 2026-11-27). ``T`` (``u = T·u_libre + g``, ADR 0004)
    lleva la escala a las componentes del residuo reducido ``TᵀR``; sin
    ``T``, todas las componentes. Cada componente tiene las unidades de la
    componente del residuo (fuerza o momento), así que la escala es tan
    invariante ante las unidades como el propio residuo, giros incluidos.

    Vale para elementos de cinemática lineal, donde la deformación sale de
    ``U``. Los corotacionales forman ``X + u`` y su error escala con las
    coordenadas, que esta escala no cubre (limitación anotada en el ADR).

    Devuelve ``None`` si todos los vectores son nulos. Un producto ``|A|·|x|``
    por término, sin copiar la matriz (núcleo de Numba sobre CSR).
    """
    acc = None
    for A, x in terms:
        x_abs = np.abs(np.asarray(x, dtype=float)).ravel()
        if not x_abs.any():
            continue
        if acc is None:
            acc = np.zeros(A.shape[0])
        if sp.issparse(A):
            A = A.tocsr()
            _abs_csr_matvec_add(A.indptr, A.indices, A.data, x_abs, acc)
        else:
            acc += np.abs(np.asarray(A, dtype=float)) @ x_abs
    if acc is None:
        return None
    if T is not None:
        acc = abs(T).T @ acc
    return _EPS * np.asarray(acc, dtype=float).ravel()


def roundoff_scale(*terms, T=None) -> float:
    """Norma euclídea de :func:`roundoff_vector`, la escala del error de
    redondeo de ``‖TᵀR‖`` que el criterio de fuerza suma, multiplicada por
    ``roundoff_factor``, como piso. Medido: el residuo de un sistema resuelto
    no pasa de 0.5 veces esta escala (bloque con holgura de contacto, Hex8
    elástico y plástico, marcos con giros en m, mm y km, Newton-Newmark y
    Newton-HHT)."""
    v = roundoff_vector(*terms, T=T)
    return 0.0 if v is None else float(np.linalg.norm(v))


def _check_factor(name: str, value) -> None:
    """Un factor del criterio es un número real finito no negativo; el
    mensaje nombra la clave también si el valor no es un número."""
    ok = isinstance(value, (int, float, np.integer, np.floating)) and not isinstance(value, bool)
    if ok:
        v = float(value)
        ok = bool(np.isfinite(v)) and v >= 0.0
    if not ok:
        raise ValueError(f"convergence: {name} = {value!r} debe ser un número finito no "
                         f"negativo.")
