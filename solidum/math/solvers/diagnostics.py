# solidum_fem/solidum/math/solvers/diagnostics.py
"""Diagnóstico tipado de divergencia para solvers no lineales (ADR 0011).

Sustituye al ``RuntimeError`` genérico ``"el solver ha divergido"`` por una
familia de excepciones que distinguen los modos de fallo cualitativamente
distintos: oscilación del Newton fuera del pozo de convergencia, singularidad
real de la tangente, carga por encima de la capacidad plástica del modelo, o
modo desconocido cuando ninguno de los patrones anteriores aplica.

Cada excepción lleva tres tipos de información:

- **Métricas numéricas**: último residuo, último ``δU``, último factor de
  carga, número de bisecciones consumidas. Útil para post-mortem programático.
- **``hint`` textual**: una pista sobre el diagnóstico probable y la acción
  típica a tomar. No es prescriptiva — el usuario decide.
- **``mode`` simbólico**: cadena corta para logging/agrupación.

El módulo es consumido por ``NonlinearSolver`` y ``NewtonNewmarkSolver``; el
``ArcLengthSolver`` queda fuera del alcance del ADR 0011 (su geometría tiene
globalización implícita y no requiere line search; ver ADR §"Decisión").
"""
from __future__ import annotations


class SolverDivergedError(RuntimeError):
    """Base de los fallos de convergencia con diagnóstico tipado (ADR 0011).

    No se lanza directamente — se lanza alguna de sus subclases. Hereda de
    ``RuntimeError`` por compatibilidad con código que captura solo el tipo
    genérico (incluidos los tests de la fase A).

    Attributes
    ----------
    mode : str
        Identificador simbólico del modo de divergencia.
    last_residual : float
        ``‖R[free_dofs]‖`` en la última iteración antes de abandonar.
    last_delta : float
        ``‖δU‖`` de la última iteración (norma L2 del incremento de Newton).
    last_load_factor : float
        Factor de carga ``λ`` (o tiempo ``t`` en dinámica) del último paso
        intentado. Para arc-length sería ``λ_iter``; aquí no aplica.
    n_bisections : int
        Bisecciones del incremento de carga / paso temporal consumidas antes
        de abandonar.
    hint : str
        Pista textual sobre el diagnóstico probable.
    """

    mode: str = "unknown"
    hint: str = ""

    def __init__(self, *, last_residual: float, last_delta: float,
                 last_load_factor: float, n_bisections: int,
                 extra_message: str = ""):
        self.last_residual = float(last_residual)
        self.last_delta = float(last_delta)
        self.last_load_factor = float(last_load_factor)
        self.n_bisections = int(n_bisections)

        msg_parts = [
            f"[{self.mode}] Solver divergió.",
            f"  Último residuo: {self.last_residual:.3e}",
            f"  Último |δU|:    {self.last_delta:.3e}",
            f"  Factor de carga: {self.last_load_factor:.4f}",
            f"  Bisecciones: {self.n_bisections}",
        ]
        if self.hint:
            msg_parts.append(f"  Hint: {self.hint}")
        if extra_message:
            msg_parts.append(f"  Detalle: {extra_message}")
        super().__init__("\n".join(msg_parts))


class OscillatingNewtonError(SolverDivergedError):
    """Newton oscila entre dos pozos sin converger (ADR 0011).

    Detectado cuando el residuo crece (o no decrece monótonamente) entre
    iteraciones consecutivas durante ≥ 3 iteraciones, **con line search
    activado**. Sin line search, oscilación leve es esperable en regímenes
    plásticos marginales y no es diagnóstica.

    Patrón canónico: el incremento de Newton aterriza fuera del pozo de
    convergencia y rebota. Si el line search activado no rescató, la
    globalización Armijo no fue suficiente — el problema requiere
    globalización más fuerte (trust region) o tiene una patología más
    profunda (mal escalado severo, formulación inconsistente).
    """

    mode = "oscillating"
    hint = (
        "El residuo no decrece entre iteraciones consecutivas a pesar del "
        "line search activado. Considera reducir la carga, revisar el "
        "escalado del problema, o evaluar un solver con globalización más "
        "fuerte (trust region, no disponible hoy)."
    )


class SingularTangentError(SolverDivergedError):
    """La rigidez tangente perdió rango (ADR 0011).

    Detectado cuando el solver lineal lanza ``RuntimeError`` por matriz
    singular, o cuando Cholesky degrada a LU **y además** el residuo crece
    tras la degradación. Síntoma de bifurcación, punto límite o problema
    físico mal planteado.
    """

    mode = "singular_tangent"
    hint = (
        "La rigidez tangente perdió rango. Probable bifurcación, punto "
        "límite o problema mal planteado. Considera ArcLengthSolver para "
        "atravesar puntos críticos."
    )


class LoadExceedsCapacityError(SolverDivergedError):
    """La carga aplicada excede la capacidad plástica del modelo (ADR 0011).

    Detectado cuando el residuo se estabiliza por encima de la tolerancia
    durante ≥ 5 iteraciones (variación relativa < 1% entre iteraciones
    consecutivas), sin oscilación clara ni singularidad. Heurística: vale
    más en plasticidad perfecta que en daño.
    """

    mode = "load_exceeds_capacity"
    hint = (
        "El residuo se estabiliza por encima de la tolerancia sin oscilar "
        "ni diverger. Probable carga > capacidad plástica del modelo. "
        "Revisa la magnitud de la carga aplicada o usa endurecimiento "
        "positivo (H > 0)."
    )


class UnknownDivergenceError(SolverDivergedError):
    """Modo de divergencia que no encaja en los patrones conocidos.

    Fallback cuando el solver agota iteraciones y bisecciones sin caer en
    ningún patrón claro. Equivalente al ``RuntimeError`` previo al ADR 0011
    pero con métricas numéricas adjuntas para post-mortem.
    """

    mode = "unknown"
    hint = (
        "Modo de divergencia no clasificado. Revisa el log de iteraciones "
        "para diagnóstico manual."
    )


# Detector de modo de divergencia
# ============================================================================

def classify_divergence(residual_history: list[float],
                        delta_history: list[float],
                        *,
                        oscillation_window: int = 3,
                        stagnation_window: int = 5,
                        stagnation_rel_tol: float = 1.0e-2,
                        singular_tangent_detected: bool = False) -> type[SolverDivergedError]:
    """Clasifica el modo de divergencia a partir del historial de iteraciones.

    Parameters
    ----------
    residual_history : list of float
        ``‖R‖`` por iteración, en orden cronológico. Debe contener al menos
        ``oscillation_window`` valores para que la detección de oscilación
        sea aplicable.
    delta_history : list of float
        ``‖δU‖`` por iteración. Solo usado para clasificación futura;
        actualmente no determina el modo (reservado).
    oscillation_window : int, default 3
        Número mínimo de iteraciones donde el residuo debe **no decrecer**
        para clasificar como ``oscillating``. "No decrecer" = la diferencia
        relativa al mínimo del historial es positiva.
    stagnation_window : int, default 5
        Número mínimo de iteraciones donde el residuo debe estabilizarse
        (variación relativa < ``stagnation_rel_tol``) para clasificar como
        ``load_exceeds_capacity``.
    stagnation_rel_tol : float, default 1e-2
        Umbral de variación relativa para considerar el residuo estable.
    singular_tangent_detected : bool, default False
        Bandera externa: ``True`` si el solver lineal ya señaló
        singularidad. Prioritario sobre las demás clasificaciones.

    Returns
    -------
    type
        La clase de excepción a lanzar (subclase de ``SolverDivergedError``).
    """
    if singular_tangent_detected:
        return SingularTangentError

    n = len(residual_history)
    if n < 2:
        return UnknownDivergenceError

    # Oscillation: residuo crece (o no decrece monótonamente) en una ventana reciente.
    if n >= oscillation_window:
        window = residual_history[-oscillation_window:]
        # Contamos cuántas veces el residuo no decreció en la ventana.
        non_decreasing = sum(1 for i in range(1, len(window))
                              if window[i] >= window[i - 1])
        # Si ≥ floor(window/2) iteraciones no decrecen, es oscilación.
        if non_decreasing >= (oscillation_window - 1) // 2 + 1:
            return OscillatingNewtonError

    # Stagnation: residuo estable durante una ventana reciente.
    if n >= stagnation_window:
        window = residual_history[-stagnation_window:]
        ref = max(abs(window[0]), 1.0e-300)
        rel_range = (max(window) - min(window)) / ref
        if rel_range < stagnation_rel_tol:
            return LoadExceedsCapacityError

    return UnknownDivergenceError
