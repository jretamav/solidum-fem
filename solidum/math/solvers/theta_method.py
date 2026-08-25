# solidum_fem/solidum/math/solvers/theta_method.py
"""``ThetaMethodSolver`` — integración temporal de primer orden (Etapa 8).

Resuelve el sistema semidiscreto de conducción de calor::

    C·Ṫ + K·T = F(t)

**No es una variante de** :class:`NewmarkSolver`. La familia Newmark integra
la ecuación de *segundo* orden ``M·ü + C·u̇ + K·u = F`` mediante hipótesis
sobre la aceleración dentro del paso. En conducción no existe segunda
derivada temporal: no hay aceleración, no hay inercia, y las hipótesis de
Newmark no tienen sobre qué aplicarse. Forzar Newmark con ``M = 0``
degenera el esquema. Es por tanto un solver de familia propia.

El régimen **estacionario** no lo usa: ``K·T = F`` lo resuelve el
``LinearSolver`` existente sin modificación.

Ver ``docs/specs/ThetaMethodSolver.md``.
"""
from __future__ import annotations

from typing import Callable

import numpy as np

from solidum.math.linalg import LUSolver, StiffnessProperties, select_solver
from solidum.math.solvers._shared import (
    CholeskyNotPositiveDefiniteError,
    _log,
)
from solidum.registry import SolverRegistry
from solidum.results import ThermalTransientResult

#: Tolerancia relativa para declarar que ``T_initial`` contradice un valor
#: Dirichlet impuesto en el mismo DOF. Escalada con la magnitud del campo
#: para ser invariante ante el cambio de unidades (K vs °C).
_DIRICHLET_MISMATCH_RTOL = 1e-9

#: Factor sobre el rango de los datos por encima del cual se declara
#: divergencia. Generoso a propósito: el objetivo es distinguir una
#: explosión numérica (θ < 0.5 fuera de su límite de estabilidad, donde el
#: campo crece órdenes de magnitud por paso) de una oscilación espuria
#: acotada, que es un resultado malo pero no un fallo del bucle.
_DIVERGENCE_FACTOR = 1e6


@SolverRegistry.register
class ThetaMethodSolver:
    """Integración temporal θ para ``C·Ṫ + K·T = F(t)`` (Etapa 8).

    Aproxima la ecuación en un punto intermedio del paso, ponderando el
    estado con ``θ ∈ [0, 1]``:

    .. math::

        \\mathbf C\\,\\frac{\\mathbf T_{n+1} - \\mathbf T_n}{\\Delta t}
        + \\mathbf K\\left[\\theta\\,\\mathbf T_{n+1}
        + (1-\\theta)\\,\\mathbf T_n\\right]
        = \\theta\\,\\mathbf F_{n+1} + (1-\\theta)\\,\\mathbf F_n

    lo que conduce al sistema por paso

    .. math::

        \\underbrace{(\\mathbf C + \\theta\\Delta t\\,\\mathbf K)}_{\\mathbf A}
        \\,\\mathbf T_{n+1}
        = \\left[\\mathbf C - (1-\\theta)\\Delta t\\,\\mathbf K\\right]\\mathbf T_n
        + \\Delta t\\left[\\theta\\mathbf F_{n+1} + (1-\\theta)\\mathbf F_n\\right]

    ``A`` es simétrica y definida positiva para ``θ > 0`` (suma de ``C``
    definida positiva y ``K`` semidefinida, con Dirichlet impuesto), luego
    el despachador algebraico (ADR 0003) elige Cholesky. Con ``Δt``
    constante **se factoriza una sola vez** y se reutiliza en todos los
    pasos: el problema es lineal y sin historia, así que ``K`` y ``C`` se
    ensamblan también una única vez. El coste por paso es una sustitución
    triangular.

    Casos particulares de ``θ``
    ---------------------------

    ==========  =========================  ======  ==========================
    ``θ``       Nombre                     Orden   Estabilidad
    ==========  =========================  ======  ==========================
    ``0``       Euler explícito            1       Condicional ``Δt ≤ 2/λ_max``
    ``1/2``     Crank-Nicolson             **2**   Incondicional (A-estable)
    ``2/3``     Galerkin                   1       Incondicional
    ``1``       Euler implícito            1       Incondicional, **L-estable**
    ==========  =========================  ======  ==========================

    Estabilidad no es ausencia de oscilaciones
    ------------------------------------------

    Crank-Nicolson es el más preciso de la tabla, pero A-estable y **no**
    L-estable: su factor de amplificación tiende a ``−1`` para los modos de
    frecuencia alta en vez de a ``0``. Ante un cambio brusco —un escalón de
    temperatura en la frontera, típico en el arranque— produce oscilaciones
    amortiguadas lentamente, con temperaturas que pueden salirse del rango
    de los datos y **violar el principio del máximo** de la ecuación de
    difusión: un resultado cualitativamente imposible, no sólo impreciso.

    Euler implícito es sólo de primer orden, pero L-estable: amortigua los
    modos altos por completo en un paso y nunca oscila. **Por eso el default
    es** ``theta = 1``. Es el mismo criterio que fijó ``lumped`` como default
    de la capacidad, aplicado al eje temporal en vez del espacial. Como
    contrapartida honesta, el solver **reporta el orden efectivo** del
    esquema (:attr:`order`, y en el ``ThermalTransientResult``) para que el
    coste en precisión del default robusto sea visible.

    Parameters
    ----------
    assembler : Assembler
        Vinculación al modelo. El dominio debe estar mallado con elementos
        térmicos y tener **al menos un Dirichlet** de temperatura.
    dt : float
        Paso temporal constante [s]. No limitado por estabilidad si
        ``θ ≥ 0.5``, sí por precisión: contrastar con el ``Δt``
        característico ``h²/α`` que se reporta al arrancar.
    n_steps : int
        Número de pasos a integrar.
    T_initial : float or array-like
        Condición inicial. Escalar ⇒ campo uniforme (el caso frecuente:
        cuerpo a temperatura ambiente); vector de tamaño ``ndof`` ⇒ valor
        por nodo (p. ej. el resultado de un estacionario previo). A
        diferencia del mecánico **no hay condición inicial de velocidad**:
        la ecuación es de primer orden y ``T_0`` la determina por completo.
    theta : float, default 1.0
        Peso del esquema ``∈ [0, 1]``. Ver tabla superior.
    lumping : {"lumped", "consistent"}, default "lumped"
        Forma de la matriz de capacidad. Default ``lumped`` por el mismo
        criterio físico que en los elementos térmicos.
    F_func : callable or None
        ``F_func(t) -> np.ndarray`` de shape ``(ndof,)``: vector global de
        cargas térmicas (fuente volumétrica y flujo de frontera ya
        ensamblados) en el instante ``t``. ``None`` ⇒ sin carga aplicada,
        el transitorio lo mueven las condiciones iniciales y Dirichlet.
    dirichlet_func : callable or None
        ``dirichlet_func(t) -> np.ndarray`` de shape ``(ndof,)``: vector de
        valores prescritos en el instante ``t``, en la misma disposición
        que el ``g`` del ``ConstraintSet``. Permite una temperatura impuesta
        variable en el tiempo (un ciclo térmico sobre una superficie).
        ``None`` ⇒ los valores declarados en el modelo, constantes.
    output_every : int, default 1
        Almacenar el campo cada ``N`` pasos. Limita la memoria en corridas
        largas; **no** afecta a la integración, que recorre todos los pasos.
        El instante inicial y el final se almacenan siempre.
    linear_algebra : str, default "auto"
        Backend algebraico (ADR 0003). ``A`` es SPD ⇒ Cholesky por defecto.

    Notes
    -----
    **Sin amortiguamiento de Rayleigh.** ``C = α·M + β·K`` es un modelo de
    disipación de la ecuación de segundo orden. En conducción la disipación
    es la conducción misma, ya contenida en ``K``. El parámetro no existe en
    este solver.
    """

    PIPELINE_KIND = "thermal_transient"

    def __init__(
        self,
        assembler,
        dt: float,
        n_steps: int,
        T_initial: float | np.ndarray,
        *,
        theta: float = 1.0,
        lumping: str = "lumped",
        F_func: Callable[[float], np.ndarray] | None = None,
        dirichlet_func: Callable[[float], np.ndarray] | None = None,
        output_every: int = 1,
        linear_algebra: str = "auto",
    ):
        cls = type(self).__name__
        if not (0.0 <= float(theta) <= 1.0):
            raise ValueError(
                f"{cls}: theta={theta} fuera de [0, 1]. Valores usuales: "
                "1.0 Euler implícito (L-estable, default), 0.5 Crank-Nicolson "
                "(2º orden), 2/3 Galerkin, 0.0 explícito."
            )
        if float(dt) <= 0.0:
            raise ValueError(
                f"{cls}: dt={dt} debe ser positivo. Es el paso temporal en "
                "segundos; contrástalo con el Δt característico h²/α que el "
                "solver reporta al arrancar."
            )
        if int(n_steps) < 1:
            raise ValueError(
                f"{cls}: n_steps={n_steps} debe ser ≥ 1."
            )
        if int(output_every) < 1:
            raise ValueError(
                f"{cls}: output_every={output_every} debe ser ≥ 1. Es cada "
                "cuántos pasos se almacena el campo; 1 almacena todos."
            )

        self.assembler = assembler
        self.dt = float(dt)
        self.n_steps = int(n_steps)
        self.T_initial = T_initial
        self.theta = float(theta)
        self.lumping = str(lumping)
        self.F_func = F_func
        self.dirichlet_func = dirichlet_func
        self.output_every = int(output_every)
        self.linear_algebra = str(linear_algebra)

    # ------------------------------------------------------------------
    # Diagnósticos declarativos
    # ------------------------------------------------------------------

    @property
    def order(self) -> int:
        """Orden de convergencia temporal efectivo del esquema.

        ``2`` sólo para Crank-Nicolson exacto (``θ = 1/2``), donde el
        término de error de truncación de primer orden se cancela; ``1``
        en cualquier otro ``θ``. Se expone como propiedad para que el coste
        en precisión del default robusto ``θ = 1`` sea consultable antes de
        correr el análisis, y no una penalización silenciosa.
        """
        return 2 if abs(self.theta - 0.5) < 1e-12 else 1

    @property
    def is_unconditionally_stable(self) -> bool:
        """``True`` si ``θ ≥ 1/2`` — ningún ``Δt`` produce divergencia."""
        return self.theta >= 0.5 - 1e-12

    @property
    def is_l_stable(self) -> bool:
        """``True`` sólo para ``θ = 1``.

        L-estabilidad ⇒ el factor de amplificación tiende a ``0`` para los
        modos altos, que quedan amortiguados por completo en un paso. Es la
        propiedad que garantiza monotonía ante un escalón, más fuerte que la
        A-estabilidad de Crank-Nicolson (factor → ``−1``, que oscila).
        """
        return abs(self.theta - 1.0) < 1e-12

    def characteristic_dt(self) -> float | None:
        """``Δt`` característico ``h²/α`` estimado sobre el modelo [s].

        ``α = k/(ρc)`` es la difusividad térmica y ``h`` el tamaño
        característico de elemento. Es el tiempo que tarda el frente térmico
        en atravesar un elemento: un ``Δt ≫ Δt_car`` es estable (con
        ``θ ≥ 1/2``) pero **no resuelve el transitorio** — salta por encima
        de la física que se quiere ver.

        Es **información, no restricción**: igual que el
        ``CentralDifferenceSolver`` no impone el ``Δt`` crítico, aquí no se
        limita el paso. Devuelve ``None`` si el modelo no permite estimarlo
        (material anisótropo, sin ``ρc`` declarada, o elementos sin
        geometría interrogable).
        """
        elements = list(self.assembler.domain.elements.values())
        if not elements:
            return None

        dt_min = np.inf
        for elem in elements:
            mat = elem.material
            try:
                alpha = mat.thermal_diffusivity
            except (AttributeError, ValueError):
                # Anisótropo (difusividad no escalar) o sin ρc declarada.
                return None
            if not np.isfinite(alpha) or alpha <= 0.0:
                return None

            h = self._element_size(elem)
            if h is None:
                return None
            dt_min = min(dt_min, h * h / alpha)

        return float(dt_min) if np.isfinite(dt_min) else None

    @staticmethod
    def _element_size(elem) -> float | None:
        """Tamaño característico ``h``: diagonal del bounding box del elemento.

        Estimador deliberadamente grosero. El diagnóstico busca el orden de
        magnitud del ``Δt`` característico, no un valor exacto: una medida
        más fina (raíz del volumen, mínima distancia internodal) cambiaría
        el resultado en un factor de orden 1, irrelevante frente a los
        órdenes de magnitud que separan un ``Δt`` razonable de uno que se
        salta el transitorio.
        """
        ndim = getattr(elem, "FLUX_DIM", None)
        if ndim is None:
            return None
        try:
            coords = elem.get_coordinate_matrix(ndim=ndim)
        except (AttributeError, TypeError, ValueError):
            return None
        extent = np.max(coords, axis=0) - np.min(coords, axis=0)
        h = float(np.linalg.norm(extent))
        return h if h > 0.0 else None

    # ------------------------------------------------------------------
    # Condición inicial
    # ------------------------------------------------------------------

    def _expand_initial(self, ndof: int) -> np.ndarray:
        """Normaliza ``T_initial`` a un vector global de tamaño ``ndof``."""
        T0 = self.T_initial
        if np.isscalar(T0) or (isinstance(T0, np.ndarray) and T0.ndim == 0):
            return np.full(ndof, float(T0))

        arr = np.asarray(T0, dtype=float).ravel()
        if arr.size != ndof:
            raise ValueError(
                f"{type(self).__name__}: T_initial es un vector de tamaño "
                f"{arr.size}, pero el modelo tiene {ndof} DOFs. Pasa un "
                "escalar para un campo inicial uniforme, o un vector con un "
                "valor por DOF de temperatura."
            )
        return arr

    def _warn_initial_dirichlet_mismatch(
        self, T0: np.ndarray, g: np.ndarray, prescribed: np.ndarray,
    ) -> None:
        """Avisa —sin abortar— si ``T_0`` contradice un Dirichlet impuesto.

        Es una discontinuidad en ``t = 0``, físicamente legítima: un choque
        térmico es exactamente eso, y es un caso de uso central. Pero es
        también la situación que más excita las oscilaciones espurias de un
        esquema A-estable no L-estable, así que merece señalarse.
        """
        if prescribed.size == 0:
            return
        diff = np.abs(T0[prescribed] - g[prescribed])
        scale = max(
            float(np.max(np.abs(g[prescribed]))),
            float(np.max(np.abs(T0))) if T0.size else 0.0,
            1.0,
        )
        bad = prescribed[diff > _DIRICHLET_MISMATCH_RTOL * scale]
        if bad.size == 0:
            return

        muestra = ", ".join(str(int(d)) for d in bad[:8])
        if bad.size > 8:
            muestra += f", … (+{bad.size - 8})"
        extra = ""
        if not self.is_l_stable:
            extra = (
                f" Con theta={self.theta:g} el esquema es A-estable pero no "
                "L-estable: este escalón puede producir oscilación amortiguada "
                "en los primeros pasos, con temperaturas fuera del rango de "
                "los datos. Usa theta=1.0 si necesitas monotonía garantizada."
            )
        _log.warning(
            "%s: T_initial contradice el valor Dirichlet impuesto en %d DOF(s) "
            "[%s]. Es una discontinuidad en t=0 — modelización legítima (choque "
            "térmico), no un error; el análisis continúa tomando el valor "
            "prescrito desde el primer paso.%s",
            type(self).__name__, int(bad.size), muestra, extra,
        )

    # ------------------------------------------------------------------
    # Integración
    # ------------------------------------------------------------------

    def solve(self) -> ThermalTransientResult:
        cls = type(self).__name__
        _log.info("--- INICIANDO SOLVER THETA-METHOD (térmico transitorio) ---")

        theta, dt, n_steps = self.theta, self.dt, self.n_steps
        nombre = {1.0: "Euler implícito", 0.5: "Crank-Nicolson",
                  0.0: "Euler explícito"}.get(round(theta, 12), f"θ={theta:g}")
        _log.info(
            "  theta=%g (%s) · orden efectivo %d · %s%s",
            theta, nombre, self.order,
            "incondicionalmente estable" if self.is_unconditionally_stable
            else "CONDICIONALMENTE estable",
            ", L-estable" if self.is_l_stable else "",
        )
        if not self.is_unconditionally_stable:
            _log.warning(
                "%s: theta=%g < 0.5 ⇒ estabilidad CONDICIONAL, el paso está "
                "limitado por Δt ≤ 2/λ_max. Con capacidad consistente λ_max es "
                "grande y el esquema resulta prácticamente inutilizable; usa "
                "lumping='lumped', o theta ≥ 0.5 para eliminar el límite.",
                cls, theta,
            )

        # --- Ensamblaje: una sola vez (problema lineal, sin historia) ---
        self.assembler.assemble_system()
        K = self.assembler.K_global
        C = self.assembler.assemble_mass_matrix(lumping=self.lumping)

        ndof = self.assembler.ndof
        cs = self.assembler.constraint_set
        T_op, g0 = cs.build(ndof)
        free_dofs = cs.free_dofs(ndof)
        n_free = T_op.shape[1]

        if n_free == ndof:
            raise ValueError(
                f"{cls}: el modelo no impone ningún Dirichlet de temperatura. "
                "La matriz de conductividad es semidefinida —tiene el modo de "
                "temperatura uniforme, que no produce gradiente ni flujo, "
                "análogo térmico de un modo de sólido rígido— y el problema "
                "queda indeterminado en el nivel absoluto de temperatura. "
                "Fija la temperatura en al menos un nodo (`fix_dof(\"T\", ...)`)."
            )

        prescribed = np.setdiff1d(np.arange(ndof), free_dofs)

        # Diagnóstico de Δt: información, no restricción.
        dt_car = self.characteristic_dt()
        if dt_car is not None:
            ratio = dt / dt_car
            msg = ("  Δt característico h²/α ≈ %.4e s · Δt usado %.4e s "
                    "(razón %.3g)")
            if ratio > 10.0:
                _log.warning(
                    msg + " — el paso es mucho mayor que el tiempo que tarda "
                    "el frente térmico en cruzar un elemento: el resultado "
                    "será estable pero puede saltarse el transitorio que se "
                    "quiere observar.", dt_car, dt, ratio,
                )
            else:
                _log.info(msg, dt_car, dt, ratio)
        else:
            _log.info("  Δt característico no estimable en este modelo "
                       "(material anisótropo o sin ρc declarada).")

        # --- Condición inicial ---
        T_full = self._expand_initial(ndof)
        self._warn_initial_dirichlet_mismatch(T_full, g0, prescribed)
        # El valor prescrito manda desde t=0: el DOF no es incógnita.
        T_free = T_full[free_dofs].copy()

        # --- Matrices reducidas y factorización única ---
        K_ff = (T_op.T @ K @ T_op).tocsr()
        C_ff = (T_op.T @ C @ T_op).tocsr()

        A = (C_ff + (theta * dt) * K_ff).tocsr()
        B = (C_ff - ((1.0 - theta) * dt) * K_ff).tocsr()

        props = StiffnessProperties(
            is_symmetric=True,
            # A = C + θΔt·K con C definida positiva y K semidefinida ⇒ SPD
            # para θ ≥ 0. Con θ = 0 se reduce a C, que sigue siendo SPD.
            is_positive_definite=True,
            size=n_free,
        )
        linalg = select_solver(props, override=self.linear_algebra)
        try:
            factor = linalg.factorize(A)
        except CholeskyNotPositiveDefiniteError:
            _log.warning(
                "%s: Cholesky reportó no-positividad sobre A = C + θΔt·K. "
                "Degradando a LU (ADR 0003 §5). Revisa que la capacidad "
                "volumétrica ρc sea positiva en todos los materiales.", cls,
            )
            factor = LUSolver().factorize(A)
        self._n_factorizations = 1   # contador para el test de reutilización

        # --- Historial: submuestreo por `output_every` ---
        stored = [k for k in range(n_steps + 1) if k % self.output_every == 0]
        if stored[-1] != n_steps:
            stored.append(n_steps)          # el instante final siempre entra
        store_at = set(stored)

        t_history = np.array([k * dt for k in stored], dtype=float)
        T_history = np.zeros((ndof, len(stored)))
        col_of = {k: i for i, k in enumerate(stored)}

        def _ensamblar_global(T_f: np.ndarray, g: np.ndarray) -> np.ndarray:
            return T_op @ T_f + g

        T_history[:, 0] = _ensamblar_global(T_free, g0)

        # --- Estado en t = 0 para el primer paso ---
        F_prev = self._F_at(0.0, ndof)
        g_prev = self._g_at(0.0, g0)
        converged = True
        # Rango de referencia para detectar explosión numérica.
        ref = max(float(np.max(np.abs(T_history[:, 0]))),
                   float(np.max(np.abs(g0))), 1.0)

        for step in range(n_steps):
            t_next = (step + 1) * dt
            F_next = self._F_at(t_next, ndof)
            g_next = self._g_at(t_next, g0)

            # rhs = B·T_n + Δt·[θ·F_{n+1} + (1−θ)·F_n] − acoplamiento Dirichlet.
            rhs = B @ T_free + dt * (
                T_op.T @ (theta * F_next + (1.0 - theta) * F_prev)
            )

            # Acoplamiento con los DOFs prescritos. Dos contribuciones:
            #   - conductividad: −Δt·K_fp·[θ·ḡ_{n+1} + (1−θ)·ḡ_n]
            #   - capacidad: −C_fp·(ḡ_{n+1} − ḡ_n), sólo si ḡ varía en el
            #     tiempo. Es el término que la mayoría de implementaciones
            #     olvida y que sesga el resultado cuando el Dirichlet es un
            #     ciclo térmico: la energía que entra por calentar el nodo
            #     prescrito no aparece si sólo se considera K.
            g_theta = theta * g_next + (1.0 - theta) * g_prev
            rhs -= dt * (T_op.T @ (K @ g_theta))
            if self.dirichlet_func is not None:
                rhs -= T_op.T @ (C @ (g_next - g_prev))

            T_free = factor.solve(rhs)

            if step + 1 in store_at:
                T_history[:, col_of[step + 1]] = _ensamblar_global(
                    T_free, g_next,
                )

            if not np.all(np.isfinite(T_free)) or (
                np.max(np.abs(T_free)) > _DIVERGENCE_FACTOR * ref
            ):
                converged = False
                _log.error(
                    "%s: divergencia detectada en el paso %d (t=%.4e s). "
                    "Con theta=%g el esquema es %s. %s",
                    cls, step + 1, t_next, theta,
                    "incondicionalmente estable — revisa los datos del modelo "
                    "(ρc positiva, conductividad bien condicionada)"
                    if self.is_unconditionally_stable else
                    "condicionalmente estable",
                    "" if self.is_unconditionally_stable else
                    f"Reduce Δt o usa theta ≥ 0.5 (Δt actual: {dt:g} s).",
                )
                break

            F_prev, g_prev = F_next, g_next

        n_done = step + 1 if not converged else n_steps
        _log.info("  -> %d pasos completados. T ∈ [%.6g, %.6g].",
                   n_done, float(np.min(T_history)), float(np.max(T_history)))

        return ThermalTransientResult(
            t_history=t_history,
            T_history=T_history,
            n_steps=n_steps,
            theta=theta,
            dt=dt,
            order=self.order,
            converged=converged,
        )

    # ------------------------------------------------------------------
    # Evaluación de cargas y Dirichlet en un instante
    # ------------------------------------------------------------------

    def _F_at(self, t: float, ndof: int) -> np.ndarray:
        if self.F_func is None:
            return np.zeros(ndof)
        F = np.asarray(self.F_func(t), dtype=float).ravel()
        if F.size != ndof:
            raise ValueError(
                f"{type(self).__name__}: F_func({t}) devolvió un vector de "
                f"tamaño {F.size}; se esperaban {ndof} (un valor por DOF)."
            )
        return F

    def _g_at(self, t: float, g0: np.ndarray) -> np.ndarray:
        if self.dirichlet_func is None:
            return g0
        g = np.asarray(self.dirichlet_func(t), dtype=float).ravel()
        if g.size != g0.size:
            raise ValueError(
                f"{type(self).__name__}: dirichlet_func({t}) devolvió un "
                f"vector de tamaño {g.size}; se esperaban {g0.size} (un valor "
                "por DOF, no nulo sólo en los DOFs prescritos)."
            )
        return g
