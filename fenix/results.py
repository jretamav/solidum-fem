"""API pública de resultados de una solución (ADR 0002).

Contiene los tipos de dato que exponen los resultados al exterior. Todos son
``frozen`` para impedir que un consumidor mute resultados y sorprenda a otro
posterior:

- :class:`ElementForces`: esfuerzos internos en ejes locales de un elemento, en
  los nodos i, j. Convenciones de signo definidas en ``Reglas.md §5``.
- :class:`SolveResult`: agregado de un análisis estático — desplazamientos,
  reacciones, fuerzas internas. Producido por ``LinearSolver``,
  ``NonlinearSolver``, ``ArcLengthSolver``.
- :class:`ModalResult`: frecuencias propias y modos M-ortonormales del problema
  ``K·φ = ω²·M·φ`` (ADR 0009 fase 1). Expone ``free_vibration`` como wrapper
  delgado sobre ``fenix.math.modal_response``.
- :class:`TransientResult`: historiales temporales ``(t, u, u̇, ü)`` de un
  análisis dinámico transitorio (ADR 0009 fases 3, 4 y 5). Producido por
  ``NewmarkSolver``, ``NewtonNewmarkSolver``, ``HHTSolver``,
  ``NewtonHHTSolver``, ``CentralDifferenceSolver``.
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
from typing import TYPE_CHECKING, Literal

import numpy as np

if TYPE_CHECKING:
    from fenix.core.domain import Domain
    from fenix.math.assembly import Assembler


ElementKind = Literal["truss", "cable", "frame2d", "frame3d"]

# Claves válidas por familia. Un componente ausente del dict significa "no aplica"
# a ese tipo de elemento (p. ej. no hay M en un truss) — no se rellena con ceros.
_VALID_COMPONENTS: dict[ElementKind, frozenset[str]] = {
    "truss":   frozenset({"N"}),
    "cable":   frozenset({"N"}),
    "frame2d": frozenset({"N", "V", "M"}),
    "frame3d": frozenset({"N", "Vy", "Vz", "T", "My", "Mz"}),
}


@dataclass(frozen=True)
class ElementForces:
    """Esfuerzos internos en ejes locales, evaluados en los nodos i, j.

    Parameters
    ----------
    kind
        Familia del elemento. Determina qué claves pueden aparecer en ``components``.
    components
        Diccionario ``{nombre: array(2,)}`` con el valor en el nodo i (índice 0) y
        en el nodo j (índice 1). Las claves válidas dependen de ``kind``:

        - ``"truss"`` / ``"cable"``: ``{"N"}``.
        - ``"frame2d"``: ``{"N", "V", "M"}``. Convención de viga estructural
          (Reglas.md §5): ``V`` rota el diferencial en horario si es positivo,
          ``M`` positivo es sagging.
        - ``"frame3d"``: ``{"N", "Vy", "Vz", "T", "My", "Mz"}``. Convención
          stress-resultant / RHR pura (Reglas.md §5).
    """

    kind: ElementKind
    components: dict[str, np.ndarray]

    def __post_init__(self) -> None:
        valid = _VALID_COMPONENTS.get(self.kind)
        if valid is None:
            raise ValueError(f"ElementForces.kind desconocido: {self.kind!r}")

        extra = set(self.components) - valid
        if extra:
            raise ValueError(
                f"Componentes inválidos para kind={self.kind!r}: {sorted(extra)}. "
                f"Válidos: {sorted(valid)}"
            )

        for name, arr in self.components.items():
            if not isinstance(arr, np.ndarray) or arr.shape != (2,):
                raise ValueError(
                    f"components[{name!r}] debe ser np.ndarray shape (2,), "
                    f"recibido {type(arr).__name__} shape={getattr(arr, 'shape', None)}"
                )

    def at_node_i(self) -> dict[str, float]:
        """Valores en el nodo i como escalares."""
        return {k: float(v[0]) for k, v in self.components.items()}

    def at_node_j(self) -> dict[str, float]:
        """Valores en el nodo j como escalares."""
        return {k: float(v[1]) for k, v in self.components.items()}


@dataclass(frozen=True)
class ModalResult:
    """Resultado de un análisis modal (ADR 0009 fase 1). Inmutable.

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
    """

    frequencies_rad: np.ndarray
    frequencies_hz: np.ndarray
    periods: np.ndarray
    modes: np.ndarray
    n_modes: int
    converged: bool = True

    def free_vibration(
        self,
        M,
        u0: np.ndarray,
        u0_dot: np.ndarray,
        t: np.ndarray,
    ) -> np.ndarray:
        """Respuesta temporal de vibración libre sin amortiguamiento por
        superposición modal.

        Wrapper delgado sobre :func:`fenix.math.modal_response.free_vibration`
        (ADR 0009 fase 7, regla D — el algoritmo vive en
        ``fenix/math/modal_response.py`` para componer con las funciones
        de combinación modal espectral). Preserva la API histórica.

        Parameters
        ----------
        M, u0, u0_dot, t
            Ver :func:`fenix.math.modal_response.free_vibration`.

        Returns
        -------
        np.ndarray, shape (n_dof, n_t)
        """
        from fenix.math.modal_response import free_vibration as _free_vibration
        return _free_vibration(
            self.modes, self.frequencies_rad, M, u0, u0_dot, t,
        )


@dataclass(frozen=True)
class TransientResult:
    """Resultado de un análisis dinámico transitorio (ADR 0009 fase 3). Inmutable.

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
        ``True`` si el integrador completó todos los pasos. ``False`` solo
        si se detuvo prematuramente por inestabilidad numérica detectada
        (no aplica al Newmark incondicionalmente estable con default
        β=1/4, γ=1/2).
    """

    t_history: np.ndarray
    u_history: np.ndarray
    udot_history: np.ndarray
    uddot_history: np.ndarray
    n_steps: int
    alpha_rayleigh: float = 0.0
    beta_rayleigh: float = 0.0
    converged: bool = True


@dataclass(frozen=True)
class HarmonicResult:
    """Resultado de un análisis de respuesta armónica (ADR 0009 fase 6). Inmutable.

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
    """

    omega: np.ndarray
    u_complex: np.ndarray
    F_complex: np.ndarray
    alpha_rayleigh: float = 0.0
    beta_rayleigh: float = 0.0

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
    Inmutable.

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
    """

    u_combined: np.ndarray
    u_per_mode: np.ndarray
    frequencies_rad: np.ndarray
    participation_factors: np.ndarray
    effective_masses: np.ndarray
    combination: str
    damping: float
    direction: np.ndarray

    @property
    def n_modes(self) -> int:
        return int(self.frequencies_rad.shape[0])

    @property
    def frequencies_hz(self) -> np.ndarray:
        return self.frequencies_rad / (2.0 * np.pi)

    def cumulative_effective_mass_ratio(self) -> np.ndarray:
        """Razón de masa efectiva acumulada por modo, shape ``(n_modes,)``.

        Útil para verificar cuántos modos hacen falta para capturar la
        respuesta. Valor objetivo en códigos sísmicos: ≥ 0.9 con los
        ``n_modes`` incluidos.
        """
        total = float(np.sum(self.effective_masses))
        if total <= 0.0:
            return np.zeros_like(self.effective_masses)
        return np.cumsum(self.effective_masses) / total


@dataclass(frozen=True)
class SolveResult:
    """Resultado agregado de una solución. Inmutable; calculado eager al final
    de ``solver.solve``.

    Parameters
    ----------
    U
        Desplazamientos globales, shape ``(n_dof,)``.
    F_applied
        Cargas aplicadas tal como las ve el sistema: nodales del usuario más
        equivalentes de elemento (distribuidas, térmicas) ya ensambladas a
        nodos. Shape ``(n_dof,)``. No incluye reacciones.
    R
        Reacciones globales, shape ``(n_dof,)``. Ceros en DOF libres; valor no
        nulo solo en DOF restringidos. Se expone redundantemente con
        ``reactions_by_node`` por conveniencia (vector para GUIs, dict para scripts).
    reactions_by_node
        Vista ``{node_id: {dof_name: valor}}`` filtrada a nodos con al menos un
        DOF restringido.
    element_forces
        ``{elem_id: ElementForces}`` con los esfuerzos internos de cada elemento
        que implementa ``internal_forces``. Elementos sólidos pueden omitirse.
    converged
        ``True`` si el solver alcanzó el criterio de convergencia. Para solvers
        lineales, siempre ``True`` tras una resolución exitosa.
    num_steps
        Número de pasos/iteraciones realizados. ``1`` para solver lineal; ``N``
        para Newton-Raphson / arc-length.
    """

    U: np.ndarray
    F_applied: np.ndarray
    R: np.ndarray
    reactions_by_node: dict[int, dict[str, float]] = field(default_factory=dict)
    element_forces: dict[int, ElementForces] = field(default_factory=dict)
    converged: bool = True
    num_steps: int = 1


def build_solve_result(
    domain: "Domain",
    assembler: "Assembler",
    U: np.ndarray,
    F_applied: np.ndarray,
    *,
    converged: bool = True,
    num_steps: int = 1,
) -> SolveResult:
    """Post-procesa una solución ``U`` para construir un ``SolveResult`` completo.

    Calcula:
    - Reacciones ``R`` globales: ``F_int(U) − F_applied`` en DOFs restringidos,
      cero en DOFs libres (por equilibrio al convergir).
    - ``reactions_by_node``: vista filtrada a nodos con al menos un DOF con
      condición de Dirichlet.
    - ``element_forces``: ``element.internal_forces(U)`` para cada elemento
      que implemente el contrato (los que devuelvan ``None`` se omiten).

    El llamador (típicamente ``fenix.run``) asigna el resultado a
    ``domain.last_result``.
    """
    ndof = U.shape[0]

    _, F_int = assembler.assemble_non_linear_system(U)

    # Reacciones por residuo en DOFs prescritos (ADR 0004 §3). En equilibrio
    # F_int[free] = F_applied[free]; en DOFs restringidos la diferencia es la
    # reacción que el apoyo ejerce sobre el sistema.
    bc_dofs = assembler.constraint_set.slave_dofs

    R = np.zeros(ndof)
    if bc_dofs.size > 0:
        R[bc_dofs] = F_int[bc_dofs] - F_applied[bc_dofs]

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

    element_forces: dict[int, ElementForces] = {}
    for elem_id, elem in domain.elements.items():
        ef = elem.internal_forces(U)
        if ef is not None:
            element_forces[elem_id] = ef

    return SolveResult(
        U=U,
        F_applied=F_applied,
        R=R,
        reactions_by_node=reactions_by_node,
        element_forces=element_forces,
        converged=converged,
        num_steps=num_steps,
    )
