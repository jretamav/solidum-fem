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

from enum import Enum

from solidum.math.linalg.eigen import SemidefiniteMassError


class Termination(str, Enum):
    """Motivo por el que terminó un análisis incremental o en el tiempo
    (ADR 0025).

    Lo dejan ``NonlinearSolver``, la familia del arco y los solvers
    transitorios (mecánicos y térmico) en su atributo ``termination`` tras
    ``solve()``, también cuando termina con una excepción. ``solidum.run``
    lo copia en ``SolveResult.termination``, y los transitorios en el de su
    resultado. Es una cadena (``Termination.MAX_STEPS == "max_steps"``), así
    que se puede comparar sin importar la clase.

    - ``MAX_LAMBDA``: alcanzó el factor de carga final (``max_lambda``; en el
      ``NonlinearSolver``, la carga completa λ = 1).
    - ``CONTROL_TARGET``: alcanzó el valor objetivo de la magnitud controlada
      (``IndirectDisplacementSolver``, parámetro ``control_target``).
    - ``MAX_STEPS``: agotó ``max_steps`` antes de completar el trazado. El
      resultado es el del último paso convergido, **no** el análisis pedido.
    - ``T_END``: un análisis en el tiempo llegó al instante final.
    - ``DIVERGED``: un paso no convergió ni reduciendo su tamaño al mínimo,
      o la integración en el tiempo se volvió inestable; ``solve()`` lanzó
      una :class:`SolverDivergedError`.
    """

    MAX_LAMBDA = "max_lambda"
    CONTROL_TARGET = "control_target"
    MAX_STEPS = "max_steps"
    T_END = "t_end"
    DIVERGED = "diverged"

    def __str__(self) -> str:
        return self.value


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
        ``‖TᵀR‖`` (residuo reducido, ADR 0004) en la última iteración antes
        de abandonar.
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
        self.extra_message = str(extra_message)
        # Historial hasta el último paso convergido de un análisis en el
        # tiempo que fracasó (ADR 0025); lo fija el solver transitorio.
        self.partial_result = None

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
        "límite, un mecanismo (apoyos insuficientes, una rótula o nudo sin "
        "rigidez suficiente) o problema mal planteado. Considera "
        "ArcLengthSolver para atravesar puntos críticos."
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


class StepConstraintError(SolverDivergedError):
    """La restricción del paso de un solver de trayectoria no tiene solución
    con la longitud de paso pedida (ADR 0025).

    En el arco cilíndrico, las raíces de la ecuación de la restricción son
    imaginarias: el arco de radio Δl no corta la trayectoria de equilibrio.
    En el control indirecto, la magnitud controlada deja de responder a la
    carga. El solver reduce el paso y lo reintenta; esta excepción se lanza
    cuando ni el paso mínimo lo resuelve. La causa concreta va en el
    ``Detalle`` del mensaje.
    """

    mode = "step_constraint"
    hint = (
        "La restricción del paso no tiene solución aunque se reduzca la "
        "longitud de paso. Con el arco cilíndrico ocurre en un retroceso "
        "(snap-back) por localización o tras un punto límite agudo: el arco "
        "no corta la trayectoria. Considera IndirectDisplacementSolver con "
        "una magnitud que crezca de forma monótona (p. ej. el alargamiento "
        "de la zona que se ablanda) o un primer paso menor (initial_dlambda)."
    )


class TimeIntegrationInstabilityError(SolverDivergedError):
    """La integración en el tiempo se volvió inestable (ADR 0025).

    La solución crece sin límite paso a paso: el paso de tiempo supera el
    límite de estabilidad de un esquema condicionalmente estable (diferencia
    central, θ-método con θ < 1/2). No es un fallo de Newton: no hay residuo
    que reportar, y ``last_residual`` es ``nan``. ``last_delta`` es la norma
    de la solución en el paso en que se detectó y ``last_load_factor`` el
    instante.
    """

    mode = "unstable_time_integration"
    hint = (
        "El paso de tiempo supera el límite de estabilidad del esquema. "
        "Reduce dt o usa un esquema incondicionalmente estable (Newmark con "
        "β = 1/4, γ = 1/2; HHT-α; θ-método con θ ≥ 1/2)."
    )

    def __init__(self, *, time: float, solution_norm: float, threshold: float,
                 extra_message: str = ""):
        self.threshold = float(threshold)
        super().__init__(
            last_residual=float("nan"), last_delta=float(solution_norm),
            last_load_factor=float(time), n_bisections=0,
            extra_message=extra_message,
        )

    def __str__(self) -> str:
        parts = [
            f"[{self.mode}] La integración en el tiempo se volvió inestable.",
            f"  Instante: t = {self.last_load_factor:.6g}",
            f"  Norma de la solución: {self.last_delta:.3e} "
            f"(umbral {self.threshold:.3e})",
            f"  Hint: {self.hint}",
        ]
        if self.extra_message:
            parts.append(f"  Detalle: {self.extra_message}")
        return "\n".join(parts)


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


# ----------------------------------------------------------------------
# Errores de planteamiento del modelo (ADR 0019)
# ----------------------------------------------------------------------

class MechanismError(ValueError):
    """El modelo puede moverse como sólido rígido: los apoyos no bastan para
    un análisis estático (ADR 0019, capa 1).

    Se detecta **antes** de resolver, comprobando qué modos de cuerpo rígido
    dejan libres las restricciones. Es un error de modelado, no numérico, así
    que es ``ValueError`` y el mensaje describe el movimiento libre en
    términos del modelo (traslación, giro alrededor de un eje por un punto,
    campo escalar sin valor prescrito), no en términos de pivotes.

    Attributes
    ----------
    motions : list[str]
        Descripción de cada movimiento rígido libre.
    """

    def __init__(self, motions: list[str], solver: str = ""):
        self.motions = list(motions)
        quien = f" ({solver})" if solver else ""
        lista = "\n".join(f"  - {m}" for m in self.motions)
        super().__init__(
            f"El modelo no está suficientemente apoyado para un análisis "
            f"estático{quien}: puede moverse como sólido rígido sin "
            f"deformarse, así que la matriz de rigidez es singular y el "
            f"equilibrio no tiene solución única. Movimientos libres:\n"
            f"{lista}\n"
            f"Añade apoyos (o restricciones lineales) que impidan esos "
            f"movimientos. Un modelo libre sí es válido en un análisis modal "
            f"o dinámico, pero no en uno estático."
        )


class IllPosedSystemError(RuntimeError):
    """El sistema de un análisis estático no es fiable (ADR 0019).

    - **Estático lineal** (capas 2 y 3): la solución no resuelve el
      sistema (error hacia atrás normado, desde el 2026-10-20), el modelo
      tiene un mecanismo interno (desde el 2026-10-20
      se decide por la energía del modo de menor energía, y se localiza), o
      la factorización directa no pudo completarse por un pivote exactamente
      nulo (``factorization_error``; hasta el 2026-09-30 salía el
      ``RuntimeError`` del backend, *Factor is exactly singular*).
    - **Estático no lineal** (``initial_tangent=True``, 2026-09-30): la
      matriz tangente inicial (``U = 0``) tiene un mecanismo interno y ningún
      componente del modelo puede ganar rigidez al deformarse (sin elementos
      geométricamente no lineales ni materiales unilaterales), así que el
      mecanismo lo es durante todo el análisis.

    Cubre lo que la detección de mecanismos rígidos no ve: mecanismos
    **internos** (una rótula de más, una barra suelta, un nudo sin rigidez)
    y matrices tan mal condicionadas que el resultado carece de sentido.

    Attributes
    ----------
    backward_error : float | None
        Error hacia atrás normado de la solución (capa 2, desde el
        2026-10-20), si fue la causa.
    zero_pivots : int | None
        Número de mecanismos internos (capa 3): uno por cada modo que se
        mueve sin deformar ningún elemento, o el número de pivotes pequeños si
        su modo no pudo calcularse. ``None`` si la causa fue el equilibrio o
        si la factorización directa no pudo completarse por un pivote
        exactamente nulo.
    mechanisms : list[str]
        Cada mecanismo interno en palabras, con los elementos que mueve
        (vacía si no se pudieron localizar o la causa fue otra).
    initial_tangent : bool
        ``True`` si se detectó en la tangente inicial de un solver no lineal.
    factorization_error : str | None
        Mensaje del backend cuya factorización falló (estático lineal).
    """

    def __init__(self, *, backward_error: float | None = None,
                 zero_pivots: int | None = None, rtol: float | None = None,
                 solver: str = "", initial_tangent: bool = False,
                 factorization_error: str | None = None,
                 mechanisms: list[str] | None = None):
        self.backward_error = backward_error
        self.zero_pivots = zero_pivots
        self.mechanisms = list(mechanisms or [])
        self.initial_tangent = bool(initial_tangent)
        self.factorization_error = factorization_error
        if self.initial_tangent:
            super().__init__(_initial_tangent_message(zero_pivots, solver, self.mechanisms))
            return
        quien = f" ({solver})" if solver else ""
        if self.mechanisms:
            super().__init__(
                f"Resultado descartado{quien}: el modelo tiene "
                f"{zero_pivots or len(self.mechanisms)} mecanismo(s) interno(s), partes que "
                f"pueden moverse sin energía de deformación en ningún elemento "
                f"aunque los apoyos impidan el movimiento de sólido rígido del "
                f"conjunto:\n"
                f"{_lista(self.mechanisms)}\n"
                f"La matriz de rigidez es singular: si la carga no activa el "
                f"mecanismo, los desplazamientos incluirían una parte arbitraria "
                f"de él; si lo activa, no hay equilibrio. Causas habituales: una "
                f"rótula o articulación de más, una barra o elemento suelto, dos "
                f"partes unidas por un solo nodo o por una sola restricción. "
                f"Revisa la conectividad de esos elementos."
            )
            return
        if zero_pivots:
            causa = (f"el solver directo encontró {zero_pivots} pivote(s) "
                     f"pequeño(s) y no se pudo calcular su modo, así que la "
                     f"matriz es singular aunque los desplazamientos parezcan "
                     f"razonables")
        elif factorization_error is not None:
            causa = (f"el solver directo no pudo factorizar la matriz, porque "
                     f"tiene al menos un pivote exactamente nulo o una fila "
                     f"vacía (backend: {factorization_error})")
        else:
            causa = (f"la solución no resuelve el sistema (error hacia atrás "
                     f"normado {backward_error:.1e}, tolerancia {rtol:.0e}): el "
                     f"solver no pudo resolverlo con precisión")
        super().__init__(
            f"Resultado descartado{quien}: {causa}. La matriz de rigidez es singular "
            f"o está tan mal condicionada que los desplazamientos no son "
            f"fiables. Causas habituales: un mecanismo interno (una rótula o "
            f"articulación de más, una barra o elemento suelto, dos partes "
            f"unidas por un solo nodo que pueden girar entre sí), o rigideces "
            f"muy desproporcionadas entre partes del modelo (penalizaciones, "
            f"unidades mezcladas). Revisa la conectividad y los apoyos."
        )


class MasslessDofsError(SemidefiniteMassError):
    """Un análisis dinámico no tiene solución con los grados de libertad sin
    masa del modelo (ADR 0024).

    En el transitorio:

    - ``initial_acceleration``: los integradores implícitos (Newmark, HHT y
      sus variantes no lineales) aún no condensan los grados de libertad sin
      masa para la aceleración inicial (deuda #34);
    - ``explicit_integration``: la diferencia central exige masa en todos
      los grados de libertad libres.

    En el modal, los grados de libertad sin masa (masas sólo en algunos
    nudos, elementos con ``density = 0``) se condensan estáticamente, de
    forma exacta. Este error se lanza sólo cuando eso no basta:

    - ``too_many_modes``: se piden más modos que grados de libertad con
      masa, y sólo ellos tienen frecuencia finita;
    - ``no_mass``: ningún grado de libertad libre tiene masa;
    - ``massless_mechanism``: los grados de libertad sin masa pueden
      moverse sin rigidez ni inercia, y su frecuencia está indeterminada;
    - ``singular_mass`` / ``indefinite_mass``: la masa de los grados de
      libertad que sí la tienen es singular, o la matriz de masa no es
      semidefinida positiva.

    Es un error de modelado (``ValueError``). Hereda de
    :class:`~solidum.math.linalg.eigen.SemidefiniteMassError`, la versión
    algebraica con índices.

    Attributes
    ----------
    massless_dofs : list[str]
        Grados de libertad libres sin masa, como ``"nodo <id> (<dof>)"``.
    n_modes : int
        Modos pedidos.
    """

    _MAX_LISTED = 12

    def __init__(self, *, cause: str, massless, massless_dofs: list[str],
                 n_mass: int, n_modes: int = 0, solver: str = "ModalSolver"):
        self.massless_dofs = list(massless_dofs)
        self.n_modes = int(n_modes)
        n_zero = len(self.massless_dofs)
        muestra = ", ".join(self.massless_dofs[:self._MAX_LISTED])
        if n_zero > self._MAX_LISTED:
            muestra += f", … (+{n_zero - self._MAX_LISTED})"
        sin_masa = f"Grados de libertad sin masa ({n_zero}): {muestra}."
        dar_masa = ("densidad en los materiales (`density`) o masas de nudo "
                    "(`Domain.add_nodal_mass`, bloque YAML `nodal_masses`)")
        if cause == "too_many_modes":
            texto = (
                f"se pidieron {n_modes} modos, pero sólo {n_mass} grado(s) de "
                f"libertad libre(s) tienen masa, y sólo ellos tienen frecuencia "
                f"finita: los demás no tienen inercia y se condensan "
                f"estáticamente (ADR 0024). Pide n_modes ≤ {n_mass}, o añade "
                f"{dar_masa}. {sin_masa}"
            )
        elif cause == "no_mass":
            texto = (
                f"ningún grado de libertad libre tiene masa, así que no hay "
                f"frecuencias finitas. Declara {dar_masa}."
            )
        elif cause == "initial_acceleration":
            texto = (
                f"la aceleración inicial ü₀ = M⁻¹·(F(0) − F_int(u₀) − C·u̇₀) no "
                f"está definida con grados de libertad libres sin masa: la "
                f"integración implícita en el tiempo todavía no los condensa "
                f"(deuda #34 de docs/STATUS.md). Declara {dar_masa}. {sin_masa}"
            )
        elif cause == "explicit_integration":
            texto = (
                f"la integración explícita resuelve ü = M⁻¹·(…) en cada paso y "
                f"exige masa en todos los grados de libertad libres. Declara "
                f"{dar_masa}, o usa un integrador implícito. {sin_masa}"
            )
        elif cause == "massless_mechanism":
            texto = (
                f"los grados de libertad sin masa forman un mecanismo: pueden "
                f"moverse sin rigidez que los fije y sin inercia, así que su "
                f"frecuencia está indeterminada. Revisa la conectividad y los "
                f"apoyos de esos grados de libertad, o dales masa. {sin_masa}"
            )
        else:
            texto = (
                f"la matriz de masa de los grados de libertad que sí tienen "
                f"masa es singular, o la matriz de masa no es semidefinida "
                f"positiva (¿una densidad o masa negativa?). {sin_masa}"
            )
        super().__init__(f"{solver}: {texto}", cause=cause, massless=massless,
                         n_mass=n_mass)


def _lista(lineas: list[str]) -> str:
    return "\n".join(f"  - {linea}" for linea in lineas)


def _initial_tangent_message(zero_pivots: int | None, solver: str,
                             mechanisms: list[str] | None = None) -> str:
    """Mensaje de :class:`IllPosedSystemError` para la tangente inicial de un
    análisis estático no lineal (ADR 0019, 2026-09-30)."""
    quien = f" ({solver})" if solver else ""
    consecuencia = (
        f"Ningún componente del modelo puede ganar rigidez al deformarse (no "
        f"hay elementos geométricamente no lineales ni materiales "
        f"unilaterales, y la no linealidad del material sólo la reduce), así "
        f"que el mecanismo lo es durante todo el análisis: si la carga no lo "
        f"activa, el análisis podría converger, pero sus desplazamientos "
        f"incluirían una parte arbitraria del mecanismo; si lo activa, no hay "
        f"equilibrio."
    )
    if mechanisms:
        return (
            f"Análisis descartado{quien}: la matriz tangente inicial (U = 0) es "
            f"singular. El modelo tiene {zero_pivots or len(mechanisms)} mecanismo(s) "
            f"interno(s), partes que pueden moverse sin energía de deformación "
            f"en ningún elemento aunque los apoyos impidan el movimiento de "
            f"sólido rígido del conjunto:\n{_lista(mechanisms)}\n{consecuencia} Causas "
            f"habituales: una rótula o articulación de más, una barra o "
            f"elemento suelto, dos partes unidas por un solo nodo o por una sola "
            f"restricción. Revisa la conectividad de esos elementos."
        )
    if zero_pivots:
        pivotes = (f"tiene {zero_pivots} pivote(s) pequeño(s) y no se pudo "
                   f"calcular su modo, así que es singular")
    else:
        pivotes = ("no se pudo factorizar porque tiene al menos un pivote "
                   "exactamente nulo: es singular")
    return (
        f"Análisis descartado{quien}: la matriz tangente inicial (U = 0) "
        f"{pivotes}. Los apoyos impiden que el modelo se mueva como sólido "
        f"rígido, así que tiene un mecanismo interno: una parte que puede "
        f"moverse sin deformar ningún elemento (una rótula o articulación de "
        f"más, una barra o elemento suelto, dos partes unidas por un solo "
        f"nodo que pueden girar entre sí), o rigideces tan desproporcionadas "
        f"entre partes del modelo (penalizaciones, unidades mezcladas) que la "
        f"matriz es singular a efectos numéricos. {consecuencia} Para "
        f"localizarlo, haz un análisis modal del mismo modelo con densidad en "
        f"los materiales: cada mecanismo aparece como una frecuencia nula y su "
        f"forma modal muestra qué parte se mueve. Revisa la conectividad y los "
        f"apoyos."
    )
