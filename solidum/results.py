"""API pública de resultados de una solución (ADR 0002).

Contiene los tipos de dato que exponen los resultados al exterior. Todos son
``frozen`` para impedir que un consumidor mute resultados y sorprenda a otro
posterior:

- :class:`ElementForces`: fuerzas internas (N, V, M, T) en ejes locales de un
  elemento estructural 1D, evaluadas en los nodos i, j. Convenciones de signo
  definidas en ``Reglas.md §5``.
- :class:`SolveResult`: agregado de un análisis estático — desplazamientos,
  reacciones de los apoyos, fuerzas de restricción de las uniones multipunto
  (ADR 0004 §3), fuerzas internas. Producido por ``LinearSolver``,
  ``NonlinearSolver``, ``ArcLengthSolver``.
- :class:`ModalResult`: frecuencias propias y modos M-ortonormales del problema
  ``K·φ = ω²·M·φ`` (ADR 0009 fase 1), con los factores de participación y las
  masas efectivas de cada traslación y los grados de libertad sin masa
  condensados. Expone ``free_vibration`` como wrapper delgado sobre
  ``solidum.math.modal_response``.
- :class:`TransientResult`: historiales temporales ``(t, u, u̇, ü)`` de un
  análisis dinámico transitorio (ADR 0009 fases 3, 4 y 5). Producido por
  ``NewmarkSolver``, ``NewtonNewmarkSolver``, ``HHTSolver``,
  ``NewtonHHTSolver``, ``CentralDifferenceSolver``.
- :class:`ThermalTransientResult`: historia ``(t, T)`` de un análisis térmico
  transitorio (Etapa 8). Producido por ``ThetaMethodSolver``. Dataclass propia
  y no ``TransientResult`` porque la ecuación térmica es de primer orden: no
  hay velocidad, aceleración ni amortiguamiento de Rayleigh.
- :class:`HarmonicResult`: amplitud compleja ``û(ω)`` por barrido en frecuencia
  (ADR 0009 fase 6). Producido por ``HarmonicSolver``. Métodos derivados
  ``.amplitude()`` y ``.phase()``.
- :class:`ResponseSpectrumResult`: envolvente máxima por combinación modal
  SRSS/CQC contra espectro de respuesta (ADR 0009 fase 7). Producido por
  ``ResponseSpectrumSolver``. Incluye factores de participación, masas
  efectivas y método ``.cumulative_effective_mass_ratio()``.

Cada solver retorna el ``Result`` específico de su pipeline; el `Domain`
cachea el último en ``domain.last_result`` independientemente del tipo.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

# ElementForces vive en ``solidum.core`` (parte del contrato de Element,
# ADR 0002). Re-exportado aquí para que la API pública
# ``from solidum.results import ElementForces`` siga funcionando.
from solidum.core.element_forces import ElementForces, ElementKind  # noqa: F401
from solidum.core.load_vector import LoadVector, as_load_vector, equivalent_loads_of

if TYPE_CHECKING:
    from solidum.core.domain import Domain
    from solidum.math.assembly import Assembler


@dataclass(frozen=True)
class ModalResult:
    """Resultado de un análisis modal (ADR 0009 fase 1). Inmutable (shallow).

    El ``frozen=True`` de la dataclass impide reasignar los atributos, pero
    los arrays NumPy y diccionarios que contiene **siguen siendo mutables
    por contenido** — un consumidor con la referencia puede hacer
    ``result.modes[0, 0] = 0`` y modificar el resultado in-place. No tratamos
    estos objetos como cápsulas opacas: el contrato es "no reasignes;
    tampoco modifiques los arrays interiores".

    Parameters
    ----------
    frequencies_rad
        Frecuencias naturales ``ω_n`` en rad/s, shape ``(n_modes,)``.
        Ordenadas ascendentemente. Modos de cuerpo rígido se reportan
        con ``ω ≈ 0`` (clip numérico aplicado).
    frequencies_hz
        ``f_n = ω_n / (2π)`` en Hz, shape ``(n_modes,)``. Mismo orden.
    periods
        ``T_n = 1/f_n`` en segundos. ``np.inf`` si ``f_n = 0``
        (modo de cuerpo rígido).
    modes
        Modos ``φ_n`` en columnas, shape ``(n_dof, n_modes)``,
        M-ortonormales (``ΦᵀMΦ = I``). En DOFs prescritos por Dirichlet
        la componente del modo es 0. Signo arbitrario por columna.
    n_modes
        Número de modos calculados. Coincide con ``frequencies_rad.size``.
    converged
        ``True`` si ARPACK convergió en los modos pedidos. ``False`` solo
        si el solver fue forzado a parar (max iteraciones); típicamente
        no debe ocurrir si las tolerancias del YAML están en el rango
        razonable por defecto.
    participation_factors
        ``{"ux": Γ, "uy": Γ, ...}``: para cada traslación del modelo, los
        factores de participación ``Γₙ = φₙᵀ·M·r``, shape ``(n_modes,)``, con
        ``r`` la traslación rígida unitaria en esa dirección (vale 1 también
        en los apoyos, que se mueven con el suelo; Chopra 2017, §13.1). Con
        modos M-ortonormales no hace falta dividir entre ``φₙᵀ·M·φₙ``. Su
        signo sigue el signo arbitrario del modo.
    effective_masses
        ``{"ux": Γ², ...}``: masas modales efectivas ``Mₙ* = Γₙ²`` en cada
        traslación, en unidades de masa. Para otro vector ``r`` (un giro
        alrededor de un punto, una dirección oblicua), usar
        :func:`solidum.math.modal_response.participation_factors` con
        :func:`solidum.math.rigid_body.rigid_body_vector`.
    condensed_dofs
        ``((id del nodo, grado de libertad), ...)``: los grados de libertad
        sin masa condensados estáticamente (ADR 0024); vacío si todos tienen
        masa.
    total_mass
        ``{"ux": rᵀ·M·r, ...}``: masa total del modelo en cada traslación,
        incluida la de los nudos con apoyo (para revisar la entrada de
        masas).
    participating_mass
        ``{"ux": bᵀ·M_LL⁻¹·b, ...}``: la masa que puede participar, suma de
        las masas efectivas de todos los modos; la referencia de
        :meth:`cumulative_effective_mass_ratio` (decisión del usuario del
        2026-11-09, ADR 0009 §8).
    eigenvalues
        Autovalores ``λₙ = ωₙ²`` con su signo, shape ``(n_modes,)``. Uno
        negativo indica que el estado sobre el que se calculó el modal (la
        rigidez tangente de una etapa, ADR 0031) es inestable en ese modo:
        su frecuencia y su periodo son ``NaN`` (:attr:`unstable_modes`).
    """

    frequencies_rad: np.ndarray
    frequencies_hz: np.ndarray
    periods: np.ndarray
    modes: np.ndarray
    n_modes: int
    converged: bool = True
    participation_factors: dict = field(default_factory=dict)
    effective_masses: dict = field(default_factory=dict)
    condensed_dofs: tuple = ()
    total_mass: dict = field(default_factory=dict)
    participating_mass: dict = field(default_factory=dict)
    eigenvalues: np.ndarray | None = None

    @property
    def unstable_modes(self) -> np.ndarray:
        """Índices (desde 0) de los modos con autovalor negativo: el estado
        sobre el que se calculó (ADR 0031) es inestable en ellos."""
        if self.eigenvalues is None:
            return np.empty(0, dtype=np.int64)
        return np.flatnonzero(np.asarray(self.eigenvalues) < 0.0)

    def cumulative_effective_mass_ratio(self, direction: str) -> np.ndarray:
        """Razón de masa efectiva acumulada en ``direction`` (``"ux"``…), shape
        ``(n_modes,)``: ``Σ_{n≤k} Mₙ* / M_part`` con ``M_part`` la masa que
        puede participar. Vale 1 con todos los modos; el criterio usual pide
        ≥ 0.9 con los calculados."""
        if direction not in self.effective_masses:
            raise ValueError(f"cumulative_effective_mass_ratio: dirección {direction!r} "
                             f"sin masas efectivas; hay {sorted(self.effective_masses)}.")
        reference = self.participating_mass.get(direction, 0.0)
        if reference <= 0.0:
            return np.zeros_like(self.effective_masses[direction])
        return np.cumsum(self.effective_masses[direction]) / reference

    def free_vibration(
        self,
        M,
        u0: np.ndarray,
        u0_dot: np.ndarray,
        t: np.ndarray,
    ) -> np.ndarray:
        """Respuesta temporal de vibración libre sin amortiguamiento por
        superposición modal.

        Wrapper delgado sobre :func:`solidum.math.modal_response.free_vibration`
        (ADR 0009 fase 7, regla D — el algoritmo vive en
        ``solidum/math/modal_response.py`` para componer con las funciones
        de combinación modal espectral). Preserva la API histórica.

        Parameters
        ----------
        M, u0, u0_dot, t
            Ver :func:`solidum.math.modal_response.free_vibration`.

        Returns
        -------
        np.ndarray, shape (n_dof, n_t)
        """
        from solidum.math.modal_response import free_vibration as _free_vibration
        return _free_vibration(
            self.modes, self.frequencies_rad, M, u0, u0_dot, t,
        )


@dataclass(frozen=True)
class TransientResult:
    """Resultado de un análisis dinámico transitorio (ADR 0009 fase 3). Inmutable (shallow).

    El ``frozen=True`` impide reasignar atributos; los historiales NumPy
    interiores siguen siendo mutables por contenido. Tratar los arrays
    como read-only por convención.

    Almacena las historias temporales completas del estado dinámico
    ``(u, u̇, ü)`` en cada paso del integrador. Mismo patrón que
    ``lambda_history`` en arc-length, generalizado a tres campos.

    Parameters
    ----------
    t_history
        Vector de instantes temporales evaluados, shape ``(n_steps + 1,)``.
        Incluye ``t = 0`` (condiciones iniciales).
    u_history, udot_history, uddot_history
        Historiales de desplazamiento, velocidad y aceleración en globales,
        shape ``(n_dof, n_steps + 1)``. La columna ``[:, k]`` corresponde
        a ``t_history[k]``. En DOFs prescritos por Dirichlet la componente
        es el valor constante del apoyo (u), o cero (u̇, ü).
    n_steps
        Número de pasos temporales realizados (``len(t_history) - 1``).
    alpha_rayleigh, beta_rayleigh
        Coeficientes Rayleigh efectivos usados en ``C = α·M + β·K``.
        ``(0.0, 0.0)`` si el análisis fue sin amortiguamiento.
    converged
        ``True`` si el integrador llegó al instante final. Un análisis que
        se detiene (Newton sin convergencia, integración inestable) lanza
        una ``SolverDivergedError`` cuyo ``partial_result`` es un
        ``TransientResult`` hasta el último paso convergido, con
        ``converged = False`` (ADR 0025).
    termination
        ``"t_end"`` (llegó al instante final) o ``"diverged"`` (resultado
        parcial de un análisis que se detuvo), ADR 0025.
    element_forces_history
        ``{elem_id: [ElementForces por paso]}`` registrado **durante** el
        análisis con el estado interno de cada paso, cuando el solver se
        construyó con ``record_internal_forces=True``; ``None`` si no se
        registró. Es lo que devuelve :meth:`internal_forces_history` cuando
        existe, y la única forma fiel con materiales con historia.
    F_constant
        Carga que actuó constante durante todo el análisis, además de la de
        ``F_func``: la del análisis estático del que continúa (ADR 0029).
        ``None`` si no continúa uno. Un análisis que continúa éste la hereda.
    reaction_labels
        ``((nodo, dof), ...)``: los grados de libertad de apoyo, una fila de
        ``R_history`` cada uno (los Dirichlet y los de las restricciones
        multipunto de apoyo, ADR 0004 §3).
    R_history
        Reacciones de los apoyos en cada instante, shape ``(n_apoyos,
        n_steps + 1)``: la ecuación de movimiento en esas filas,
        ``R = F_int(u) + M·ü + C·u̇ − F`` (decisión del usuario del
        2026-11-09): la fuerza de los elementos, la inercia de la masa del
        apoyo o acoplada a él y el amortiguamiento. Cumple el equilibrio de
        toda la estructura. ``None`` en un resultado construido sin ellas.
    R_damping_history
        La parte de amortiguamiento de ``R_history``, ``(C·u̇)`` en los
        apoyos: restarla da la reacción sin el amortiguamiento de Rayleigh
        (``damping=False`` en :meth:`total_reaction` y
        :meth:`reactions_by_node`).
    """

    t_history: np.ndarray
    u_history: np.ndarray
    udot_history: np.ndarray
    uddot_history: np.ndarray
    n_steps: int
    alpha_rayleigh: float = 0.0
    beta_rayleigh: float = 0.0
    converged: bool = True
    element_forces_history: dict[int, list["ElementForces"]] | None = None
    termination: str | None = None
    F_constant: np.ndarray | None = None
    reaction_labels: tuple = ()
    R_history: np.ndarray | None = None
    R_damping_history: np.ndarray | None = None

    def _reactions(self, damping: bool) -> np.ndarray:
        if self.R_history is None:
            raise ValueError("Este resultado no lleva reacciones (R_history = None).")
        return self.R_history if damping else self.R_history - self.R_damping_history

    def total_reaction(self, dof_name: str, *, damping: bool = True) -> np.ndarray:
        """Suma de las reacciones de todos los apoyos en ``dof_name`` en cada
        instante, shape ``(n_steps + 1,)``: el cortante basal con ``"ux"``.
        Con ``damping=False``, sin la parte de amortiguamiento."""
        rows = _rows_of(self.reaction_labels, dof_name)
        return self._reactions(damping)[rows].sum(axis=0)

    def reactions_by_node(self, step: int = -1, *, damping: bool = True) -> dict:
        """``{nodo: {dof: valor}}`` de las reacciones en el instante ``step``
        (por omisión, el último)."""
        return _by_node(self.reaction_labels, self._reactions(damping)[:, step])

    def internal_forces_history(
        self, domain: "Domain",
    ) -> dict[int, list["ElementForces"]]:
        """Fuerzas internas (N, V, M, T) por paso temporal para cada elemento.

        Si el solver registró las fuerzas durante el análisis
        (``record_internal_forces=True`` → :attr:`element_forces_history`),
        devuelve ese registro, fiel a la historia interna de cada paso.
        Si no, las reconstruye de forma *lazy* como se describe abajo.

        Para cada ``elem_id`` con ``internal_forces`` implementado (contrato
        ADR 0002 — barras/vigas/cables 2D y 3D), devuelve la lista de
        ``ElementForces`` evaluados con el campo ``u_history[:, k]`` en
        cada paso ``k``. Los elementos cuya ``internal_forces(U)`` devuelve
        ``None`` (sólidos 2D — deuda ADR 0002 documentada) se omiten del
        dict.

        Returns
        -------
        history
            ``{elem_id: [ElementForces_t_0, ..., ElementForces_t_n_steps]}``.
            Cada lista tiene longitud ``n_steps + 1`` (incluye ``t=0``).

        Notes
        -----
        Cierra la asimetría arquitectural detectada en la auditoría H-1.3:
        ``SolveResult.element_forces`` exponía las fuerzas internas estáticas pero
        ``TransientResult`` exigía al consumidor reconstruirlos manualmente.

        **Limitación — sólo exacto sin historia**: ``internal_forces``
        reevalúa cada paso desde ``u_history[:, k]`` con el estado interno
        **final** de los elementos (el committed al terminar el análisis).
        Es exacto para materiales sin variables internas y elementos
        geométricamente lineales (``Truss2D``, ``Frame2DEuler``,
        ``Frame2DTimoshenko``, ``Frame3D`` con material elástico). Con
        **cualquier** historia —plasticidad o daño (``ε^p``, ``α``, ``κ``
        del último paso aplicados a desplazamientos anteriores) o
        cinemática corotacional— las fuerzas de pasos intermedios son
        incorrectas, no sólo aproximadas (auditoría 2026-09-22). El método
        avisa por log cuando detecta materiales con ``PRIMARY_STATE_VAR``;
        la historia fiel exige que el solver almacene ``F_int`` por paso.

        Cálculo eager (recorre todos los pasos al llamar). Si necesitas
        múltiples consultas, almacena el dict devuelto externamente — no
        se cachea en la dataclass para no inflar la memoria del resultado
        cuando el consumidor no lo necesita.
        """
        if self.element_forces_history is not None:
            return {k: list(v) for k, v in self.element_forces_history.items()}
        with_history = sorted({
            type(m).__name__ for e in domain.elements.values()
            for m in e.materials()
            if getattr(m, 'PRIMARY_STATE_VAR', None)
        })
        if with_history:
            from solidum.logging import get_logger
            get_logger("results").warning(
                f"internal_forces_history: materiales con historia {with_history}; "
                f"las fuerzas de pasos intermedios se evalúan con el estado interno "
                f"final y NO son fieles a la historia. Construya el solver con "
                f"record_internal_forces=True para registrarlas durante el análisis."
            )
        history: dict[int, list[ElementForces]] = {}
        n_t = self.u_history.shape[1]
        for elem_id, elem in domain.elements.items():
            # Sondeo: si el primer paso devuelve ``None``, el elemento no
            # implementa el contrato — saltar todo el elemento.
            first = elem.internal_forces(self.u_history[:, 0])
            if first is None:
                continue
            per_step = [first]
            for k in range(1, n_t):
                per_step.append(elem.internal_forces(self.u_history[:, k]))
            history[elem_id] = per_step
        return history


@dataclass(frozen=True)
class ThermalTransientResult:
    """Resultado de un análisis térmico transitorio (Etapa 8). Inmutable (shallow).

    El ``frozen=True`` impide reasignar atributos; los historiales NumPy
    interiores siguen siendo mutables por contenido. Tratar como read-only.

    Producido por :class:`ThetaMethodSolver`. Almacena la historia del campo
    de temperatura al integrar ``C·Ṫ + K·T = F(t)``.

    **Por qué no reutiliza** :class:`TransientResult`. La ecuación térmica es
    de **primer orden** en el tiempo: no existen ``udot_history`` ni
    ``uddot_history`` (no hay velocidad ni aceleración de la temperatura con
    significado físico en el modelo), y ``alpha_rayleigh``/``beta_rayleigh``
    describen un amortiguamiento viscoso propio de la ecuación de segundo
    orden que aquí no aplica —la disipación es la conducción misma, ya
    contenida en ``K``. Rellenar esos campos con ceros produciría un objeto
    que miente sobre su contenido. Se sigue el criterio ya establecido con
    ``ModalResult``, ``HarmonicResult`` y ``ResponseSpectrumResult``: un tipo
    de análisis con semántica propia tiene su dataclass propia.

    Parameters
    ----------
    t_history
        Instantes almacenados, shape ``(n_stored,)``. Incluye siempre
        ``t = 0`` (condición inicial) y el instante final. Con
        ``output_every > 1`` los instantes intermedios se omiten, pero la
        integración recorre **todos** los pasos: el submuestreo afecta al
        almacenamiento, no al cálculo.
    T_history
        Campo de temperatura global, shape ``(n_dof, n_stored)``. La columna
        ``[:, k]`` corresponde a ``t_history[k]``. En DOFs prescritos por
        Dirichlet la componente es el valor impuesto en ese instante (que
        puede variar en el tiempo).
    n_steps
        Número de pasos temporales **integrados**, no de columnas
        almacenadas. Con ``output_every = 1`` se cumple
        ``len(t_history) == n_steps + 1``.
    theta
        Peso del esquema efectivamente utilizado.
    dt
        Paso temporal constante utilizado [s].
    order
        Orden de convergencia temporal efectivo del esquema: ``2`` si
        ``theta == 0.5`` (Crank-Nicolson), ``1`` en cualquier otro caso.
        Se expone —en vez de dejarlo implícito— para que el coste en
        precisión del default robusto ``theta = 1`` sea visible al
        consumidor y no una penalización silenciosa.
    converged
        ``True`` si el integrador llegó al instante final. Con
        ``theta < 0.5`` y un paso que excede el límite condicional de
        estabilidad, ``solve`` lanza ``TimeIntegrationInstabilityError``
        (ADR 0025), cuyo ``partial_result`` es un ``ThermalTransientResult``
        hasta el último instante almacenado antes de la divergencia, con
        ``converged = False``.
    termination
        ``"t_end"`` o ``"diverged"`` (resultado parcial), ADR 0025.
    F_constant
        Carga que actuó constante durante todo el análisis, además de la de
        ``F_func`` (fuentes y flujos): la del análisis estático del que continúa (ADR 0029).
        ``None`` si no continúa uno. Un análisis que continúa éste la hereda.
    """

    t_history: np.ndarray
    T_history: np.ndarray
    n_steps: int
    theta: float
    dt: float
    order: int
    converged: bool = True
    termination: str | None = None
    F_constant: np.ndarray | None = None

    @property
    def T_final(self) -> np.ndarray:
        """Campo de temperatura en el último instante almacenado."""
        return self.T_history[:, -1]

    def temperature_at(self, dof: int) -> np.ndarray:
        """Historia temporal de un DOF concreto, shape ``(n_stored,)``.

        Atajo sobre ``T_history[dof, :]`` para el caso frecuente de graficar
        la evolución en un punto de control.
        """
        return self.T_history[dof, :]

    def extremes(self) -> tuple[float, float]:
        """``(T_min, T_max)`` sobre todos los nodos e instantes almacenados.

        Herramienta de diagnóstico del **principio del máximo**: en un
        problema sin fuente, la solución exacta de la ecuación de difusión
        nunca excede los extremos de los datos iniciales y de frontera. Un
        resultado fuera de ese rango señala oscilación espuria del esquema
        —típicamente Crank-Nicolson ante un escalón (ver
        ``docs/specs/ThetaMethodSolver.md`` §3)— y no un fenómeno físico.
        """
        return float(np.min(self.T_history)), float(np.max(self.T_history))


@dataclass(frozen=True)
class HarmonicResult:
    """Resultado de un análisis de respuesta armónica (ADR 0009 fase 6). Inmutable (shallow).

    El ``frozen=True`` impide reasignar atributos; los arrays NumPy
    interiores siguen siendo mutables por contenido. Tratar como read-only.

    Almacena la amplitud compleja del desplazamiento ``û(ω)`` para cada
    frecuencia del barrido. La fase queda codificada en el argumento del
    complejo: ``|û|`` es la amplitud y ``arg(û)`` el desfase respecto a la
    excitación. El régimen estacionario es asumido (transitorio no
    incluido — el análisis es lineal en el dominio de la frecuencia).

    Parameters
    ----------
    omega
        Vector de frecuencias angulares evaluadas, shape ``(n_omega,)``,
        en rad/s. Las frecuencias en Hz se obtienen vía ``omega / (2π)``.
    u_complex
        Amplitudes complejas del desplazamiento global, shape
        ``(n_dof, n_omega)``. La columna ``[:, k]`` corresponde a
        ``omega[k]``. En DOFs prescritos por Dirichlet la componente es
        cero (el apoyo no oscila).
    F_complex
        Amplitudes complejas de la carga aplicada, shape ``(n_dof,)``.
        Constante a lo largo del barrido (la dependencia con ω está en la
        respuesta, no en la entrada). Puede tener parte imaginaria si la
        carga lleva fase.
    alpha_rayleigh, beta_rayleigh
        Coeficientes Rayleigh efectivos usados en ``C = α·M + β·K``.
        ``(0.0, 0.0)`` si el análisis fue sin amortiguamiento (caso
        singular en resonancia exacta — pico infinito teórico).
    reaction_labels
        ``((nodo, dof), ...)``: los grados de libertad de apoyo, una fila de
        ``R_complex`` cada uno.
    R_complex
        Amplitud compleja de las reacciones, shape ``(n_apoyos, n_omega)``:
        ``R̂ = (K − ω²·M + iω·C)·û − F̂`` en los apoyos, con las mismas partes
        que la del transitorio (decisión del usuario del 2026-11-09).
    R_damping_complex
        Su parte de amortiguamiento, ``iω·C·û`` en los apoyos.
    """

    omega: np.ndarray
    u_complex: np.ndarray
    F_complex: np.ndarray
    alpha_rayleigh: float = 0.0
    beta_rayleigh: float = 0.0
    reaction_labels: tuple = ()
    R_complex: np.ndarray | None = None
    R_damping_complex: np.ndarray | None = None

    def total_reaction(self, dof_name: str, *, damping: bool = True) -> np.ndarray:
        """Suma compleja de las reacciones de los apoyos en ``dof_name`` para
        cada frecuencia, shape ``(n_omega,)``; ``damping=False`` la da sin la
        parte de amortiguamiento."""
        if self.R_complex is None:
            raise ValueError("Este resultado no lleva reacciones (R_complex = None).")
        R = self.R_complex if damping else self.R_complex - self.R_damping_complex
        return R[_rows_of(self.reaction_labels, dof_name)].sum(axis=0)

    @property
    def frequencies_hz(self) -> np.ndarray:
        """Frecuencias del barrido en Hz."""
        return self.omega / (2.0 * np.pi)

    @property
    def n_omega(self) -> int:
        return int(self.omega.shape[0])

    def amplitude(self) -> np.ndarray:
        """``|û(ω)|`` shape ``(n_dof, n_omega)`` — espectro de amplitudes."""
        return np.abs(self.u_complex)

    def phase(self) -> np.ndarray:
        """``arg(û(ω))`` shape ``(n_dof, n_omega)`` en radianes."""
        return np.angle(self.u_complex)


@dataclass(frozen=True)
class ResponseSpectrumResult:
    """Resultado de un análisis de respuesta espectral (ADR 0009 fase 7).
    Inmutable (shallow): ``frozen=True`` impide reasignar atributos; los
    arrays NumPy interiores siguen siendo mutables por contenido. Tratar
    como read-only.

    Captura la respuesta máxima envolvente bajo una excitación
    caracterizada por un espectro de respuesta — no es una historia
    temporal, sino la combinación estadística de las respuestas modales
    máximas. La fase relativa entre modos se pierde por construcción;
    sólo se conservan amplitudes.

    Parameters
    ----------
    u_combined
        Respuesta máxima combinada (SRSS o CQC), shape ``(n_dof,)``.
        Positiva por construcción (raíz cuadrada de suma de cuadrados).
    u_per_mode
        Contribución modal individual ``Γ_n·φ_n·S_d(ω_n)`` (con signo),
        shape ``(n_dof, n_modes)``. Útil para diagnóstico de qué modos
        dominan la respuesta.
    frequencies_rad
        Frecuencias naturales utilizadas en la combinación, shape
        ``(n_modes,)``, rad/s.
    participation_factors
        Factores de participación ``Γ_n = φ_nᵀ·M·r``, shape
        ``(n_modes,)``.
    effective_masses
        Masas modales efectivas ``Γ_n²`` (M-ortonormalidad), shape
        ``(n_modes,)``. La suma sobre todos los modos converge a la
        masa total del modelo en la dirección de excitación cuando se
        incluyen suficientes modos.
    combination
        Identificador del método de combinación usado: ``"SRSS"`` o
        ``"CQC"``.
    damping
        Amortiguamiento modal asumido para el espectro y para la
        combinación CQC. ``0.0`` si SRSS (CQC requiere ``ξ > 0`` para
        coeficientes de correlación no triviales).
    direction
        Vector de excitación rígida ``r`` aplicado al modelo, shape
        ``(n_dof,)``. Para excitación sísmica unidireccional, sólo los
        DOFs traslacionales en esa dirección son no nulos.
    reaction_labels
        ``((nodo, dof), ...)``: los grados de libertad de apoyo, una fila de
        ``R_per_mode`` y de ``R_combined`` cada uno.
    R_per_mode
        Reacción de cada modo, ``(K − ωₙ²·M)·uₙ`` en los apoyos (con signo),
        shape ``(n_apoyos, n_modes)``.
    R_combined
        Reacción máxima de cada apoyo, combinada como los desplazamientos
        (SRSS o CQC), shape ``(n_apoyos,)``. Para la de varios apoyos juntos
        (el cortante basal) use :meth:`total_reaction`, que suma antes de
        combinar.
    """

    u_combined: np.ndarray
    u_per_mode: np.ndarray
    frequencies_rad: np.ndarray
    participation_factors: np.ndarray
    effective_masses: np.ndarray
    combination: str
    damping: float
    direction: np.ndarray
    # Masa total en la dirección de excitación, ``rᵀ·M·r``, incluida la de
    # los nudos con apoyo (para revisar la entrada de masas).
    total_mass: float | None = None
    # Masa que puede participar, ``bᵀ·M_LL⁻¹·b`` (suma de las masas efectivas
    # de todos los modos): la referencia de la razón acumulada (decisión del
    # usuario del 2026-11-09, ADR 0009 §8). ``None`` sólo en resultados
    # construidos a mano (se degrada a ``total_mass``).
    participating_mass: float | None = None
    reaction_labels: tuple = ()
    R_per_mode: np.ndarray | None = None
    R_combined: np.ndarray | None = None

    def total_reaction(self, dof_name: str) -> float:
        """Máximo de la suma de las reacciones de todos los apoyos en
        ``dof_name`` (el cortante basal con ``"ux"``): se suma en cada modo y
        después se combina (SRSS o CQC), porque la combinación no es lineal."""
        if self.R_per_mode is None:
            raise ValueError("Este resultado no lleva reacciones (R_per_mode = None).")
        from solidum.math.modal_response import combine_modal
        per_mode = self.R_per_mode[_rows_of(self.reaction_labels, dof_name)].sum(axis=0)
        return float(combine_modal(per_mode[None, :], self.frequencies_rad,
                                   self.combination, self.damping)[0])

    def reactions_by_node(self) -> dict:
        """``{nodo: {dof: valor}}`` de las reacciones máximas combinadas."""
        return _by_node(self.reaction_labels, self.R_combined)

    @property
    def n_modes(self) -> int:
        return int(self.frequencies_rad.shape[0])

    @property
    def frequencies_hz(self) -> np.ndarray:
        return self.frequencies_rad / (2.0 * np.pi)

    def cumulative_effective_mass_ratio(self) -> np.ndarray:
        """Razón de masa efectiva acumulada por modo, shape ``(n_modes,)``.

        ``Σ_{n≤k} Γ_n² / M_part``: la fracción de la masa que puede
        participar en la dirección de excitación (``participating_mass``,
        la suma de las masas efectivas de todos los modos) que capturan los
        primeros ``k`` modos; vale 1 con todos. Criterio normativo típico:
        ≥ 0.9 con los ``n_modes`` incluidos. Antes del 2026-11-09 se dividía
        entre ``rᵀ·M·r``, que incluye la masa de los apoyos y no se alcanzaba
        ni con todos los modos (0.875 en un marco con dos elementos por
        barra); antes del 2026-09-22, entre la suma de los modos calculados,
        que daba siempre 1.0.
        """
        total = (float(self.participating_mass) if self.participating_mass is not None
                 else float(self.total_mass) if self.total_mass is not None
                 else float(np.sum(self.effective_masses)))
        if total <= 0.0:
            return np.zeros_like(self.effective_masses)
        return np.cumsum(self.effective_masses) / total


@dataclass(frozen=True)
class StagedResult:
    """Resultados de un análisis por etapas declarado en el YAML (``stages``,
    ADR 0026 y 0029): uno por etapa, en orden.

    ``result[i]`` o ``result["nombre"]`` da el de una etapa; ``final``, el de
    la última; ``names``, los nombres. Cada etapa continúa la anterior (sus
    desplazamientos y el estado interno de los elementos; la carga ya
    aplicada se mantiene constante).
    """

    names: tuple
    results: tuple

    def __getitem__(self, key):
        if isinstance(key, str):
            try:
                return self.results[self.names.index(key)]
            except ValueError:
                raise KeyError(f"no hay una etapa {key!r}; etapas: {list(self.names)}") from None
        return self.results[key]

    def __len__(self) -> int:
        return len(self.results)

    def __iter__(self):
        return iter(self.results)

    @property
    def final(self):
        """Resultado de la última etapa."""
        return self.results[-1]


@dataclass(frozen=True)
class SolveResult:
    """Resultado agregado de una solución. Inmutable (shallow); calculado
    eager al final de ``solver.solve``.

    El ``frozen=True`` impide reasignar atributos; los arrays NumPy
    interiores siguen siendo mutables por contenido. Tratar como read-only.

    Parameters
    ----------
    U
        Desplazamientos globales, shape ``(n_dof,)``.
    F_applied
        Cargas aplicadas tal como las ve el sistema: nodales del usuario más
        equivalentes de elemento (distribuidas, térmicas) ya ensambladas a
        nodos, **escaladas por el factor de carga final** del solver
        (``λ_final · F_ref`` en arc-length; ``λ = 1`` en los demás). Shape
        ``(n_dof,)``. No incluye reacciones.
    R
        Reacciones globales, shape ``(n_dof,)``: la fuerza que los **apoyos**
        ejercen sobre la estructura (ADR 0004 §3, convención decidida el
        2026-09-30). Son apoyos los DOF con Dirichlet —que reportan su
        reacción completa, incluida la parte que les llega por nodos atados a
        ellos— y las restricciones multipunto (MPC) que impiden un movimiento
        de sólido rígido, como un rodillo oblicuo ``u_y = tanθ·u_x``: su
        fuerza entra en los DOF de su esclavo y de sus maestros (la reacción
        del rodillo es normal al plano de rodadura). Cero en los demás DOF.
        Las MPC de unión no aportan: su fuerza va a ``constraint_forces``. Se
        expone redundantemente con ``reactions_by_node`` por conveniencia
        (vector para GUIs, dict para scripts).
    reactions_by_node
        Vista ``{node_id: {dof_name: valor}}`` de ``R``: los nodos con algún
        Dirichlet (sus DOF prescritos) y los DOF de las MPC de apoyo.
    element_forces
        ``{elem_id: ElementForces}`` con las fuerzas internas de cada elemento
        que implementa ``internal_forces``. Elementos sólidos pueden omitirse.
    converged
        ``True`` si el solver completó el análisis: siempre para el lineal y
        el Newton-Raphson (que lanzan excepción si no convergen); en la
        familia del arco, si llegó a ``max_lambda`` o al objetivo de la
        magnitud controlada (``control_target``), y ``False`` cuando el
        trazado se detuvo por ``max_steps``.
    num_steps
        Número de pasos convergidos. ``1`` para solver lineal; ``N`` para
        Newton-Raphson / arc-length.
    constraint_forces
        Fuerzas de restricción de las MPC de **unión** (atado de mallas, nodo
        colgante, periodicidad, unión rígida), shape ``(n_dof,)``: la fuerza
        que cada unión ejerce sobre la estructura en los DOF de su esclavo y
        de sus maestros. Transmiten fuerza entre partes de la estructura, no
        la toman del exterior: no son reacciones, y por dirección suman cero.
        Una MPC es de unión si una traslación rígida de toda la estructura (o
        un campo escalar uniforme, como una temperatura uniforme) la cumple.
        Con ``R`` cumple ``F_int − F_applied = R + constraint_forces`` en todos
        los DOF, salvo el residuo de equilibrio del solver en los libres.
        ``build_solve_result`` lo rellena siempre (ceros sin uniones); ``None``
        sólo en un ``SolveResult`` construido a mano sin él.
    constraint_forces_by_node
        Vista ``{node_id: {dof_name: valor}}`` de ``constraint_forces`` en los
        DOF de las MPC de unión.
    condition_estimate
        Número de condición estimado ``κ̃`` de la rigidez reducida y escalada
        por la diagonal (ADR 0019, capa 4): ``log10(κ̃)`` acota las cifras que
        la solución puede haber perdido de las ~16 de la doble precisión, y
        por encima de ``1e12`` el solver avisa. Lo da el ``LinearSolver``
        (estática lineal y conducción estacionaria) con un backend directo y
        ``K`` simétrica; ``None`` en el resto.
    termination
        Motivo de la parada de un solver incremental (ADR 0025):
        ``"max_lambda"``, ``"control_target"`` o ``"max_steps"`` (un fracaso
        lanza excepción y no llega a construir el resultado). ``None`` en el
        ``LinearSolver``, que resuelve en un paso.
    load_factor
        Factor de carga final ``λ`` de la carga que crece en el análisis (o
        en la etapa): ``F_applied = F_constant + λ·F_ref``.
    F_constant
        Carga constante de un análisis por etapas (ADR 0026): la ya aplicada
        al empezar la etapa, que actúa completa durante toda ella. ``None``
        fuera de una etapa. ``F_variable`` da la otra parte, ``λ·F_ref``.
    perturbation
        ``True`` si es una perturbación lineal sobre una etapa (``LinearSolver``
        con ``previous``, ADR 0031): el estado de la etapa más el incremento
        lineal con su rigidez tangente, con reacciones y fuerzas internas
        lineales. No cambia el estado interno de los elementos, así que no se
        puede continuar: las etapas siguientes continúan la etapa base.
    equivalent_loads
        Cargas nodales equivalentes de la carga distribuida (peso propio,
        fuerza de cuerpo) por elemento, ya multiplicadas por su factor, que
        se restaron a las fuerzas nodales para dar ``element_forces``. Una
        etapa siguiente las hereda como parte de su carga constante.
    """

    U: np.ndarray
    F_applied: np.ndarray
    R: np.ndarray
    reactions_by_node: dict[int, dict[str, float]] = field(default_factory=dict)
    element_forces: dict[int, ElementForces] = field(default_factory=dict)
    converged: bool = True
    num_steps: int = 1
    constraint_forces: np.ndarray | None = None
    constraint_forces_by_node: dict[int, dict[str, float]] = field(default_factory=dict)
    condition_estimate: float | None = None
    termination: str | None = None
    load_factor: float = 1.0
    F_constant: np.ndarray | None = None
    equivalent_loads: dict[int, np.ndarray] = field(default_factory=dict)
    # Perturbación lineal sobre una etapa (ADR 0031): no cambió el estado
    # interno de los elementos; una etapa posterior no puede continuarla.
    perturbation: bool = False
    # Versión del registro de cargas equivalentes del ensamblador con que se
    # construyó el resultado (ADR 0026): una etapa siguiente la compara para
    # saber si las cargas registradas son suyas o de esta etapa.
    _equivalent_loads_version: int | None = field(default=None, repr=False)

    @property
    def F_variable(self) -> np.ndarray:
        """Parte de la carga que creció con λ: ``F_applied − F_constant``."""
        if self.F_constant is None:
            return self.F_applied
        return self.F_applied - self.F_constant


# Tolerancia relativa del criterio de unión: |C_k·b| ≤ tol·Σ_j |C_kj·b_j|.
_JOINT_RTOL = 1e-10


def _is_joint(constraint, B: np.ndarray) -> bool:
    """``True`` si la MPC es de unión: toda columna ``b`` de ``B`` (las
    traslaciones rígidas y los campos escalares uniformes de la base rígida)
    la cumple, ``|C_k·b| ≤ tol·Σ_j |C_kj·b_j|``.

    Un DOF de giro vale cero en una traslación, así que una unión rígida
    ``u_s = u_m + d·θ_m`` es de unión.
    """
    dofs = [constraint.slave_dof, *constraint.masters]
    coefs = np.array([1.0, *(-a for a in constraint.coefficients)])
    rows = B[dofs]                               # (len(dofs), m)
    value = coefs @ rows
    scale = np.abs(coefs) @ np.abs(rows)
    return bool(np.all(np.abs(value) <= _JOINT_RTOL * scale))


def _dof_labels_of_constraints(domain: "Domain") -> dict[int, tuple[int, str]]:
    """Número de ecuación → ``(node_id, dof_name)`` de los DOF que aparecen en
    las restricciones lineales del dominio. Recorre sólo las restricciones,
    no todos los nodos."""
    labels: dict[int, tuple[int, str]] = {}
    for spec in getattr(domain, "linear_constraints", []):
        for nid, dof_name in [tuple(spec["slave"]), *(tuple(m) for m in spec["masters"])]:
            labels[domain.nodes[nid].dofs[dof_name]] = (nid, dof_name)
    return labels


def _split_constraint_forces(assembler, residual):
    """Reparte el residuo entre reacciones de apoyo y fuerzas de unión.

    Devuelve ``(R, F_c, support_dofs, joint_dofs)``: los vectores de
    reacciones y de fuerzas de restricción y los DOF (con MPC) de las
    restricciones de apoyo y de unión. Ver ``build_solve_result``.
    """
    cs = assembler.constraint_set
    ndof = residual.shape[0]
    dofs, mu = cs.multiplier_arrays(residual)
    declared = cs.declared
    is_dirichlet = np.fromiter((c.is_dirichlet() for c in declared), dtype=bool,
                               count=len(declared))
    R = np.zeros(ndof)
    F_c = np.zeros(ndof)
    R[dofs[is_dirichlet]] = mu[is_dirichlet]

    basis = assembler.rigid_basis()
    cols = [j for j, (kind, _) in enumerate(basis.labels)
            if kind in ("translation", "scalar")]
    B = basis.B[:, cols]

    support_dofs: list[int] = []
    joint_dofs: list[int] = []
    for i in np.flatnonzero(~is_dirichlet):
        c = declared[i]
        if _is_joint(c, B):
            target, bucket = F_c, joint_dofs
        else:
            target, bucket = R, support_dofs
        m_k = mu[i]
        target[c.slave_dof] += m_k
        for m, alpha in zip(c.masters, c.coefficients):
            target[m] -= alpha * m_k
        bucket.append(c.slave_dof)
        bucket.extend(c.masters)
    return R, F_c, support_dofs, joint_dofs


def reaction_dofs(assembler) -> np.ndarray:
    """Grados de libertad en los que un análisis reporta reacciones, ordenados:
    los Dirichlet y los de las restricciones multipunto de apoyo (ADR 0004
    §3; las de unión no aportan)."""
    cs = assembler.constraint_set
    if not cs.has_multipoint:
        return np.unique(np.asarray(cs.slave_dofs, dtype=np.int64))
    ndof = assembler.domain.total_dofs
    _, _, support, _ = _split_constraint_forces(assembler, np.zeros(ndof))
    dirichlet = [c.slave_dof for c in cs.declared if c.is_dirichlet()]
    return np.unique(np.asarray(dirichlet + list(support), dtype=np.int64))


def reactions_from_residual(assembler, residual: np.ndarray) -> np.ndarray:
    """Reacciones de los apoyos, vector global, a partir del residuo de la
    ecuación de equilibrio (o de movimiento) ``r`` en todos los grados de
    libertad: el mismo reparto entre apoyos y uniones que ``build_solve_result``.
    Es lineal en ``r``: un residuo complejo se reparte por partes."""
    residual = np.asarray(residual)
    if np.iscomplexobj(residual):
        return (reactions_from_residual(assembler, residual.real)
                + 1j * reactions_from_residual(assembler, residual.imag))
    cs = assembler.constraint_set
    if not cs.has_multipoint:
        R = np.zeros(residual.shape[0])
        bc = cs.slave_dofs
        if len(bc):
            R[bc] = residual[bc]
        return R
    return _split_constraint_forces(assembler, residual)[0]


def reaction_labels(domain: "Domain", dofs) -> tuple:
    """``((nodo, dof), ...)`` de los grados de libertad globales ``dofs``."""
    index = {int(i): (node.id, name)
             for node in domain.nodes.values() for name, i in node.dofs.items()}
    return tuple(index[int(d)] for d in dofs)


def _rows_of(labels: tuple, dof_name: str) -> list[int]:
    rows = [i for i, (_, d) in enumerate(labels) if d == dof_name]
    if not rows:
        raise ValueError(f"ningún apoyo tiene reacción en {dof_name!r}; los apoyos "
                         f"son {sorted({d for _, d in labels})}.")
    return rows


def _by_node(labels: tuple, values) -> dict:
    out: dict[int, dict[str, float]] = {}
    for (node_id, dof), value in zip(labels, values):
        out.setdefault(node_id, {})[dof] = value.item() if hasattr(value, "item") else value
    return out


def legacy_equivalent_loads(assembler, load_factor: float) -> dict[int, np.ndarray]:
    """Cargas equivalentes de un vector aplicado que no las lleva (``np.ndarray``
    simple, ADR 0027): se supone que contiene la última carga de cuerpo
    ensamblada, escalada por ``load_factor``, como antes del ADR 0027.

    Esa suposición es ambigua si se ensamblaron varias cargas de cuerpo: se
    avisa, porque las fuerzas internas de extremo pueden salir mal. Las cargas
    que devuelven ``assemble_self_weight`` y ``assemble_body_load`` llevan las
    suyas; se pierden al convertirlas con ``np.asarray``/``np.array`` o al
    sumarlas sobre un arreglo de ceros con ``+=``.
    """
    recorded = getattr(assembler, "equivalent_loads", None) or {}
    if not recorded:
        return {}
    if getattr(assembler, "equivalent_loads_version", 0) > 1:
        from solidum.logging import get_logger
        get_logger("results").warning(
            "Fuerzas internas: la carga aplicada es un np.ndarray simple, sin "
            "sus cargas equivalentes por elemento, y se ensamblaron %d cargas de "
            "cuerpo. Se descuenta sólo la última (escalada por %.6g), que puede "
            "no ser la aplicada: los extremos de los elementos con carga "
            "distribuida pueden salir mal. Combine las cargas que devuelven "
            "assemble_self_weight/assemble_body_load con + y escalares (F = F1 + "
            "λ·F2) sin convertirlas con np.asarray ni sumarlas con += sobre un "
            "arreglo de ceros (ADR 0027).",
            assembler.equivalent_loads_version, load_factor,
        )
    return {eid: load_factor * np.asarray(f, dtype=float) for eid, f in recorded.items()}


def build_solve_result(
    domain: "Domain",
    assembler: "Assembler",
    U: np.ndarray,
    F_applied: np.ndarray,
    *,
    converged: bool = True,
    num_steps: int = 1,
    load_factor: float = 1.0,
    condition_estimate: float | None = None,
    termination: str | None = None,
    F_constant: np.ndarray | None = None,
    equivalent_loads: dict[int, np.ndarray] | None = None,
    F_int: np.ndarray | None = None,
    element_forces: dict | None = None,
    perturbation: bool = False,
) -> SolveResult:
    """Post-procesa una solución ``U`` para construir un ``SolveResult`` completo.

    Calcula:
    - Reacciones ``R`` y fuerzas de restricción ``constraint_forces`` a partir
      del residuo ``r = F_int(U) − F_applied`` (ADR 0004 §3). Con ``C`` las
      restricciones declaradas sin cerrar (fila ``e_dᵀ`` de un Dirichlet,
      ``(e_s − Σ α_i·e_mi)ᵀ`` de una MPC), los multiplicadores ``μ`` de
      ``Cᵀ·μ = r`` (:meth:`ConstraintSet.multipliers`) dan la fuerza
      ``μ_k·C_kᵀ`` de cada restricción sobre la estructura. ``R`` suma las de
      los Dirichlet y las de las MPC de apoyo; ``constraint_forces``, las de
      las MPC de unión. Una MPC es de unión si ``|C_k·b| ≤ 1e-10·Σ_j |C_kj·b_j|``
      para toda columna ``b`` de traslación o de campo escalar de la base
      rígida del dominio (``Assembler.rigid_basis``); si no, es de apoyo.
      Así ``F_int − F_applied = R + constraint_forces`` en todos los DOF,
      salvo el residuo de equilibrio del solver en los libres, y un apoyo
      cuyo nodo es maestro de una MPC recibe su reacción completa.
      Sin MPC, ``R = F_int − F_applied`` en los DOF prescritos, exactamente
      como antes de la convención (sin recorrer los DOF en Python).
    - ``reactions_by_node``: nodos con algún Dirichlet (sus DOF prescritos,
      también los declarados como restricción lineal sin maestros) más los DOF
      de las MPC de apoyo.
    - ``constraint_forces_by_node``: los DOF de las MPC de unión.
    - ``element_forces``: ``element.internal_forces(U)`` para cada elemento
      que implemente el contrato (los que devuelvan ``None`` se omiten).

    ``condition_estimate`` se copia tal cual: lo calcula el solver (capa 4
    del ADR 0019). Una perturbación lineal (ADR 0031) pasa ``F_int`` y
    ``element_forces`` ya calculados, lineales en el incremento, en vez de
    evaluar los elementos en ``U``, y marca ``perturbation``. El llamador (típicamente ``solidum.run``) asigna el
    resultado a ``domain.last_result``.
    """
    ndof = U.shape[0]

    if F_int is None:
        _, F_int = assembler.assemble_non_linear_system(U)

    cs = assembler.constraint_set
    support_dofs: list[int] = []
    joint_dofs: list[int] = []
    if not cs.has_multipoint:
        # Sólo Dirichlet: la reacción es el residuo en los DOF prescritos (en
        # equilibrio F_int[libres] = F_applied[libres]).
        bc_dofs = cs.slave_dofs
        R = np.zeros(ndof)
        if bc_dofs.size > 0:
            R[bc_dofs] = F_int[bc_dofs] - F_applied[bc_dofs]
        F_c = np.zeros(ndof)
    else:
        R, F_c, support_dofs, joint_dofs = _split_constraint_forces(
            assembler, F_int - F_applied)

    reactions_by_node: dict[int, dict[str, float]] = {}
    for node_id, node in domain.nodes.items():
        if not node.boundary_conditions:
            continue
        per_dof = {
            dof_name: float(R[node.dofs[dof_name]])
            for dof_name in node.boundary_conditions
            if dof_name in node.dofs
        }
        if per_dof:
            reactions_by_node[node_id] = per_dof

    constraint_forces_by_node: dict[int, dict[str, float]] = {}
    linear_constraints = getattr(domain, "linear_constraints", [])
    if linear_constraints:
        labels = _dof_labels_of_constraints(domain)
        # Un Dirichlet declarado como restricción lineal sin maestros es un
        # apoyo como cualquier otro.
        extra_support = [domain.nodes[spec["slave"][0]].dofs[spec["slave"][1]]
                         for spec in linear_constraints if not spec["masters"]]
        for dof in (*support_dofs, *extra_support):
            nid, dof_name = labels[dof]
            reactions_by_node.setdefault(nid, {})[dof_name] = float(R[dof])
        for dof in joint_dofs:
            nid, dof_name = labels[dof]
            constraint_forces_by_node.setdefault(nid, {})[dof_name] = float(F_c[dof])

    # Cargas nodales equivalentes de la carga distribuida (peso propio /
    # fuerza de cuerpo) por elemento: las fuerzas internas de extremo son
    # F_int − f_eq, no K·u. Las lleva la carga aplicada (``LoadVector``, ADR
    # 0027), combinadas con la aritmética de las cargas (también la de una
    # etapa, ADR 0026). Un ``np.ndarray`` simple las perdió: se recurre al
    # registro heredado del ensamblador (ver ``legacy_equivalent_loads``).
    if equivalent_loads is None:
        carried = equivalent_loads_of(F_applied)
        if isinstance(F_applied, LoadVector):
            equivalent_loads = {} if carried is None else carried.as_dict()
        else:
            equivalent_loads = legacy_equivalent_loads(assembler, load_factor)
    if element_forces is None:
        element_forces = {}
        for elem_id, elem in domain.elements.items():
            f_eq = equivalent_loads.get(elem_id)
            ef = elem.internal_forces(U, equivalent_load=f_eq)
            if ef is not None:
                element_forces[elem_id] = ef

    if not isinstance(F_applied, LoadVector):
        F_applied = as_load_vector(F_applied, equivalent_loads)
    return SolveResult(
        U=U,
        F_applied=F_applied,
        R=R,
        reactions_by_node=reactions_by_node,
        element_forces=element_forces,
        converged=converged,
        num_steps=num_steps,
        constraint_forces=F_c,
        constraint_forces_by_node=constraint_forces_by_node,
        condition_estimate=condition_estimate,
        termination=None if termination is None else str(termination),
        load_factor=float(load_factor),
        F_constant=(None if F_constant is None else
                    F_constant if isinstance(F_constant, LoadVector)
                    else np.asarray(F_constant, dtype=float)),
        equivalent_loads=equivalent_loads,
        perturbation=bool(perturbation),
        _equivalent_loads_version=getattr(assembler, "equivalent_loads_version", None),
    )
