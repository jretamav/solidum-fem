# solidum_fem/solidum/math/solvers/_shared.py
"""Utilidades compartidas entre solvers del paquete ``solidum.math.solvers``.

- :data:`_log`: logger común ``"solvers"`` reutilizado por todos los módulos
  para mantener el nombre histórico en filtros y captura de tests.
- :func:`domain_is_symmetric`: agregación de los flags declarativos
  ``IS_SYMMETRIC`` / ``PRESERVES_SYMMETRY`` (ADR 0003 §2). Lo consumen
  ``LinearSolver``, ``NonlinearSolver``, ``ArcLengthSolver`` y ``NewmarkSolver``
  al construir las ``StiffnessProperties`` que entran al despachador.
- :class:`CholeskyNotPositiveDefiniteError`: alias del error real cuando
  ``scikit-sparse`` está instalado; placeholder inerte si no lo está
  (así los ``except`` no fallan en imports). ADR 0003 §5.
- :func:`warn_initial_conditions_off_constraints`: aviso de log de los
  solvers transitorios cuando ``u0`` o ``u̇0`` incumplen una restricción
  multipunto con maestros (ADR 0004).
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from solidum.constants import INITIAL_CONDITION_MPC_RTOL
from solidum.core.load_vector import LoadVector, as_load_vector, equivalent_loads_of
from solidum.logging import get_logger
from solidum.math.linalg.eigen import massless_dofs
from solidum.math.solvers.diagnostics import MasslessDofsError, Termination


_log = get_logger("solvers")


class TransientTerminationMixin:
    """Estado de terminación de un análisis en el tiempo (ADR 0025).

    Tras ``solve()`` —también si termina con una excepción— el solver deja:

    - ``termination``: ``"t_end"`` si llegó al instante final, ``"diverged"``
      si se detuvo (``Termination``); ``None`` antes de integrar o si
      ``solve`` falló antes del primer paso (modelo mal planteado).
    - ``completed``: ``True`` si llegó al instante final.
    - ``steps_done``, ``t_final``: último paso convergido y su instante.
    - ``U_final``: campo global (desplazamientos o temperaturas) de ese paso.
    - ``partial_result``: si se detuvo, el resultado hasta el último paso
      convergido (también en ``exc.partial_result``); ``None`` si terminó.

    ``solve`` lo reinicia al empezar con :meth:`_begin_trace`.

    Progreso: si ``step_callback`` no es ``None``, se llama
    ``step_callback(step, U, t)`` tras cada paso convergido (``step ≥ 1``),
    con el campo global (desplazamientos o temperaturas) y el instante: la
    misma forma que el ``callback(step, U, λ)`` de los solvers estáticos.

    Continuación (ADR 0029): :meth:`continue_from` hace que el análisis parta
    del estado final de otro —un estático (``SolveResult``) o un transitorio
    del mismo tipo—, con su instante (``t_start``) y su carga constante
    (``constant_load``), que se suma en todo instante a la de ``F_func``.
    """

    termination: Termination | None = None
    steps_done: int = 0
    t_final: float = 0.0
    U_final: np.ndarray | None = None
    partial_result = None
    #: ``callback(step, U, t)`` tras cada paso convergido, o ``None``.
    step_callback = None

    def _begin_trace(self) -> None:
        self.termination = None
        self.steps_done = 0
        self.t_final = 0.0
        self.U_final = None
        self.partial_result = None

    @property
    def completed(self) -> bool:
        return self.termination == Termination.T_END

    # ------------------------------------------------------------------
    # Continuación de un análisis anterior (ADR 0029)
    # ------------------------------------------------------------------

    #: Instante en que empieza el análisis (``t_history[0]``).
    t_start: float = 0.0
    #: Carga que actúa constante en todo instante, además de ``F_func(t)``.
    constant_load = None
    _continues_static: bool = False
    _continuing: bool = False
    _carried_uddot = None

    @property
    def _is_thermal(self) -> bool:
        return getattr(self, "PIPELINE_KIND", None) == "thermal_transient"

    def continue_from(self, previous) -> None:
        """Hace que este análisis continúe ``previous`` (ADR 0029).

        - Desde un **estático** (``SolveResult``): parte de sus desplazamientos
          (o temperaturas) en reposo, en ``t = 0``, y su carga aplicada se
          mantiene constante (con sus cargas distribuidas, ADR 0027); ``F_func``
          añade la que varía. Si el estado no está en equilibrio con esa carga,
          se avisa: la diferencia daría una aceleración inicial espuria.
        - Desde un **transitorio** del mismo tipo: parte de su estado final
          (desplazamientos y velocidades, o temperaturas) en su instante final,
          con la misma carga constante. ``F_func`` se evalúa en tiempo
          absoluto y ``t_end`` es el instante final absoluto.
        - Un resultado parcial (análisis detenido) se admite: el nuevo parte
          del último paso convergido, p. ej. con un paso de tiempo menor.

        El estado interno de los elementos es el que guarda el ensamblador:
        use el mismo. Las condiciones iniciales del constructor (``u0``,
        ``u0_dot``, ``T_initial``) y ``previous`` se excluyen.
        """
        from solidum.results import (
            SolveResult, ThermalTransientResult, TransientResult,
        )
        name = type(self).__name__
        rate = None
        if isinstance(previous, SolveResult):
            state = np.asarray(previous.U, dtype=float)
            t0 = 0.0
            constant = as_load_vector(previous.F_applied, previous.equivalent_loads)
            static = True
        elif isinstance(previous, TransientResult) and not self._is_thermal:
            state = previous.u_history[:, -1]
            rate = previous.udot_history[:, -1]
            t0 = float(previous.t_history[-1])
            constant = previous.F_constant
            static = False
            # La aceleración final del anterior, no la del equilibrio
            # instantáneo: con HHT-α no coinciden (la ecuación del método
            # promedia con α), y continuar debe ser igual que no parar.
            self._carried_uddot = np.array(previous.uddot_history[:, -1], dtype=float)
        elif isinstance(previous, ThermalTransientResult) and self._is_thermal:
            state = previous.T_history[:, -1]
            t0 = float(previous.t_history[-1])
            constant = previous.F_constant
            static = False
        else:
            kind = "térmico" if self._is_thermal else "mecánico"
            raise TypeError(
                f"{name}.continue_from: un análisis {kind} en el tiempo continúa un "
                f"estático (SolveResult) o un transitorio del mismo tipo; se recibió "
                f"{type(previous).__name__}.")
        if self._is_thermal:
            if getattr(self, "T_initial", None) is not None and not self._continuing:
                raise ValueError(
                    f"{name}: T_initial y el análisis anterior (previous) se excluyen: "
                    f"la temperatura inicial es la final de aquél.")
            self.T_initial = np.array(state, dtype=float)
        else:
            if (self.u0 is not None or self.u0_dot is not None) and not self._continuing:
                raise ValueError(
                    f"{name}: u0/u0_dot y el análisis anterior (previous) se excluyen: "
                    f"el estado inicial es el final de aquél.")
            self.u0 = np.array(state, dtype=float)
            self.u0_dot = None if rate is None else np.array(rate, dtype=float)
        self.t_start = t0
        self.constant_load = constant
        self._continues_static = static
        self._continuing = True
        if not getattr(previous, "converged", True):
            _log.warning(
                f"{name}: continúa un análisis que no terminó "
                f"({getattr(previous, 'termination', None)}); parte de su último "
                f"estado convergido, en t = {t0:.6g}.")

    def _carried_acceleration(self, free_dofs):
        """Aceleración inicial heredada de un transitorio anterior (en los
        grados de libertad libres), o ``None``: entonces se calcula del
        equilibrio en el instante inicial."""
        if self._carried_uddot is None:
            return None
        return self._carried_uddot[free_dofs].copy()

    def _external_load(self, t: float, ndof: int):
        """``(F_crudo, F)`` en ``t``: la carga constante más ``F_func(t)``. El
        crudo conserva las cargas equivalentes por elemento (ADR 0027)."""
        raw, F = external_load(self.F_func, t, ndof)
        if self.constant_load is None:
            return raw, F
        total = self.constant_load if raw is None else self.constant_load + raw
        return total, np.asarray(total, dtype=float).reshape(ndof)

    def _number_of_steps(self, dt: float) -> int:
        """Pasos de ``t_start`` a ``t_end`` (instante final absoluto)."""
        if not self.t_end > self.t_start:
            raise ValueError(
                f"{type(self).__name__}: t_end = {self.t_end:.6g} no es posterior al "
                f"instante en que empieza el análisis, t_start = {self.t_start:.6g}. "
                f"Al continuar un transitorio, t_end es el instante final absoluto.")
        return number_of_steps(self.t_end - self.t_start, dt)

    def _check_initial_equilibrium(self, T_op, F_int_red) -> None:
        """Al continuar un estático, el estado inicial debe equilibrar la carga
        constante: ``Tᵀ·F_c − F_int,red(u₀) ≈ 0`` (fuerzas internas reducidas,
        apoyos prescritos incluidos). Si no, avisa: la diferencia se convierte
        en una aceleración inicial espuria."""
        if not self._continues_static or self.constant_load is None:
            return
        Fc = T_op.T @ np.asarray(self.constant_load, dtype=float).ravel()
        Fi = np.asarray(F_int_red, dtype=float).ravel()
        r = Fc - Fi
        scale = max(float(np.max(np.abs(Fc), initial=0.0)),
                    float(np.max(np.abs(Fi), initial=0.0)), 1e-300)
        rel = float(np.max(np.abs(r), initial=0.0)) / scale
        if rel > 1e-6:
            _log.warning(
                f"{type(self).__name__}: el estado inicial no está en equilibrio con "
                f"la carga constante del análisis anterior (residuo relativo "
                f"{rel:.2e}): la diferencia actúa como una fuerza súbita en t = "
                f"{self.t_start:.6g}. ¿El análisis estático no convergió, o cambió "
                f"el modelo entre los dos?")


def external_load(F_func, t: float, ndof: int):
    """``(F_crudo, F)``: lo que devuelve ``F_func(t)`` (un ``LoadVector`` lleva
    sus cargas equivalentes, ADR 0027) y su vector global ``(ndof,)``;
    ``(None, 0)`` sin ``F_func``."""
    if F_func is None:
        return None, np.zeros(ndof)
    raw = F_func(t)
    return raw, np.asarray(raw, dtype=float).reshape(ndof)


def mechanical_trace(solver, t_history, ndof: int, alpha_r: float, beta_r: float):
    """:class:`TransientTrace` de un transitorio mecánico (``TransientResult``)."""
    from solidum.results import TransientResult

    def build(h, t, n_steps, termination, forces):
        return TransientResult(
            t_history=t, u_history=h["u"], udot_history=h["udot"],
            uddot_history=h["uddot"], n_steps=int(t.size - 1),
            alpha_rayleigh=alpha_r, beta_rayleigh=beta_r,
            converged=termination == Termination.T_END.value,
            termination=termination, element_forces_history=forces,
            F_constant=getattr(solver, "constant_load", None),
        )

    return TransientTrace(solver, t_history=t_history, n_rows=ndof,
                          fields=("u", "udot", "uddot"), build=build,
                          record_forces=solver.record_internal_forces)


class TransientTrace:
    """Historiales de un análisis transitorio y su terminación (ADR 0025).

    Reserva los historiales de ``n_steps + 1`` instantes, guarda cada paso
    convergido (:meth:`store`) y construye el resultado al llegar al final
    (:meth:`finish`) o hasta el último paso convergido si el análisis se
    detiene (:meth:`fail`). ``fields`` nombra los historiales (``"u"``,
    ``"udot"``, ``"uddot"`` en mecánica; ``"T"`` en el térmico) y
    ``build(histories, t, n_steps, termination, forces)`` crea el resultado.
    """

    def __init__(self, solver, *, t_history, n_rows: int, fields: tuple[str, ...],
                 build, record_forces: bool = False, store_columns=None):
        self.solver = solver
        self.t = np.asarray(t_history, dtype=float)
        self.fields = fields
        self.histories = {f: np.zeros((n_rows, self.t.size)) for f in fields}
        self.build = build
        self.forces = ElementForcesRecorder(
            solver.assembler.domain, record_forces, assembler=solver.assembler)
        self.last = -1          # última columna guardada
        self.step = 0           # último paso convergido
        # Columna de cada paso cuando sólo se guardan algunos (térmico con
        # ``output_every``); ``None`` = una columna por paso.
        self.store_columns = store_columns

    def converged(self, step: int, t: float, U: np.ndarray) -> None:
        """Marca el paso ``step`` (instante ``t``) como convergido y avisa al
        ``step_callback`` del solver, si lo tiene (no en el estado inicial)."""
        self.step = int(step)
        self.solver.steps_done = int(step)
        self.solver.t_final = float(t)
        self._U = U
        callback = getattr(self.solver, "step_callback", None)
        if callback is not None and step > 0:
            callback(int(step), U, float(t))

    def store(self, column: int, *values, F_ext=None) -> None:
        for name, value in zip(self.fields, values):
            self.histories[name][:, column] = value
        if self.fields[0] == "u":
            self.forces.record(self.histories["u"][:, column], F_ext)
        self.last = column

    def _result(self, termination: Termination):
        n = self.last + 1
        histories = {f: (h if n == h.shape[1] else h[:, :n].copy())
                     for f, h in self.histories.items()}
        t = self.t if n == self.t.size else self.t[:n].copy()
        return self.build(histories, t, self.step, str(termination), self.forces.result())

    def finish(self):
        s = self.solver
        s.termination = Termination.T_END
        s.U_final = np.array(self._U, dtype=float, copy=True)
        return self._result(Termination.T_END)

    def fail(self, error, t_attempt: float):
        """Deja el estado de un análisis que se detiene y devuelve la excepción
        tipada (con ``partial_result``), lista para lanzar."""
        s = self.solver
        s.termination = Termination.DIVERGED
        U = getattr(self, "_U", None)
        s.U_final = None if U is None else np.array(U, dtype=float, copy=True)
        partial = self._result(Termination.DIVERGED) if self.last >= 0 else None
        s.partial_result = partial
        error.partial_result = partial
        _log.error(
            f"{type(s).__name__}: análisis detenido ({error.mode}) en t = "
            f"{t_attempt:.6g}. Último paso convergido: {s.steps_done} "
            f"(t = {s.t_final:.6g}); el historial hasta ahí está en "
            f"partial_result.")
        return error


def domain_is_symmetric(domain) -> bool:
    """``True`` si todos los componentes del dominio preservan simetría de ``K``.

    Recorre elementos y materiales y agrega los flags declarativos
    ``PRESERVES_SYMMETRY`` y ``IS_SYMMETRIC``. Default en ambos es ``True`` —
    sólo plasticidad no asociada, follower loads y similares lo desactivan.

    Los flags se leen **a nivel de instancia** (no de clase) para permitir
    que un material derive su simetría dinámicamente del estado del
    constructor (p. ej. ``DruckerPrager2D`` asociado con ``ψ=φ`` tiene
    tangente simétrica; con ψ≠φ no). La instancia hereda del ClassVar
    cuando no hay override, así que el caso clásico sigue funcionando.

    Se recorren **todos** los materiales de cada elemento
    (``element.materials()``, ADR 0021 §4): un marco con sección por fibras
    tiene los de sus fibras.
    """
    for elem in domain.elements.values():
        if not getattr(elem, "PRESERVES_SYMMETRY", True):
            return False
        for material in elem.materials():
            if not getattr(material, "IS_SYMMETRIC", True):
                return False
    return True


try:
    from solidum.math.linalg.cholesky import CholeskyNotPositiveDefiniteError
except ImportError:
    class CholeskyNotPositiveDefiniteError(Exception):
        """Placeholder cuando ``scikit-sparse`` no está instalado (nunca se lanza)."""


# Tolerancia relativa para reconocer ``t_end = n·dt`` pese al redondeo
# binario (``3 * 0.1 / 0.1 = 3.0000000000000004``). Sin ella ``ceil`` añade
# un paso espurio en una fracción apreciable de las parejas ``(n, dt)``.
_STEP_COUNT_RTOL = 1.0e-9


def number_of_steps(t_end: float, dt: float) -> int:
    """Número de pasos temporales que cubren ``[0, t_end]`` con paso ``dt``.

    ``ceil(t_end/dt)`` salvo que el cociente esté a menos de
    ``_STEP_COUNT_RTOL`` (relativo) de un entero, en cuyo caso se toma ese
    entero: así ``t_end = k·T`` con ``dt = T/N`` produce exactamente ``k·N``
    pasos y el instante final es ``t_end`` y no ``t_end + dt``.
    """
    if dt <= 0.0:
        raise ValueError(f"dt={dt} debe ser > 0.")
    ratio = float(t_end) / float(dt)
    nearest = round(ratio)
    if abs(ratio - nearest) <= _STEP_COUNT_RTOL * max(1.0, abs(ratio)):
        return int(nearest)
    import math
    return int(math.ceil(ratio))


class ElementForcesRecorder:
    """Registra ``elem.internal_forces(U)`` en cada paso comiteado de un
    análisis transitorio, para ``TransientResult.element_forces_history``.

    Con materiales con historia (plasticidad, daño) o elementos
    corotacionales, las fuerzas internas de un paso sólo pueden evaluarse
    **durante** el análisis, con el estado interno de ese paso: reevaluar
    después desde ``u_history`` usa el estado final y no es fiel (deuda #15
    de ``STATUS.md``). El registro es opt-in (``record_internal_forces``)
    porque cuesta una evaluación por elemento 1D y paso.

    Las fuerzas de extremo descuentan las cargas equivalentes de la carga
    distribuida del instante (``F_int − f_eq``, ADR 0002), que lleva la carga
    externa ``F(t)`` cuando es un ``LoadVector`` (ADR 0027): con
    ``F(t) = F_peso + a(t)·F_viento``, ``f_eq = f_eq,peso + a(t)·f_eq,viento``.
    Si ``F(t)`` es un ``np.ndarray`` simple y el ensamblador registró cargas
    de cuerpo, no hay forma de saber con qué factor actúan en cada instante:
    se avisa una vez y no se descuentan.

    Sólo participan los elementos cuyo ``internal_forces`` no devuelve
    ``None`` (estructurales 1D, ADR 0002); se sondea en el primer registro.
    """

    def __init__(self, domain, enabled: bool, assembler=None):
        self.domain = domain
        self.enabled = bool(enabled)
        self.assembler = assembler
        self._history: dict[int, list] = {}
        self._elements: list | None = None
        self._warned_plain = False

    def _equivalent_loads(self, F_ext):
        carried = equivalent_loads_of(F_ext)
        if carried is not None or isinstance(F_ext, LoadVector):
            return carried
        recorded = getattr(self.assembler, "equivalent_loads", None)
        if recorded and not self._warned_plain:
            self._warned_plain = True
            _log.warning(
                "Historial de fuerzas internas: la carga externa F(t) es un "
                "np.ndarray simple y el ensamblador tiene cargas de cuerpo "
                "ensambladas; no se sabe con qué factor actúan en cada instante, "
                "así que los extremos de los elementos con carga distribuida NO "
                "descuentan sus cargas equivalentes. Construya F(t) combinando "
                "con + y escalares las cargas que devuelven assemble_self_weight / "
                "assemble_body_load (ADR 0027).")
        return None

    def record(self, U: np.ndarray, F_ext=None) -> None:
        """Añade el paso actual; llamar tras comitear el estado del paso.
        ``F_ext`` es la carga externa del instante (para ``f_eq``)."""
        if not self.enabled:
            return
        U = np.asarray(U, dtype=float)
        eq = self._equivalent_loads(F_ext)

        def forces_of(elem_id, elem):
            f_eq = None if eq is None else eq.get(elem_id)
            return elem.internal_forces(U, equivalent_load=f_eq)

        if self._elements is None:
            self._elements = []
            for elem_id, elem in self.domain.elements.items():
                forces = forces_of(elem_id, elem)
                if forces is None:
                    continue
                self._elements.append((elem_id, elem))
                self._history[elem_id] = [forces]
            return
        for elem_id, elem in self._elements:
            self._history[elem_id].append(forces_of(elem_id, elem))

    def result(self) -> dict[int, list] | None:
        """``{elem_id: [ElementForces por paso]}`` o ``None`` si no se registró."""
        if not self.enabled:
            return None
        return {k: list(v) for k, v in self._history.items()}


def _dof_index_map(domain) -> dict[int, tuple[int, str]]:
    """``índice global → (id del nodo, nombre del grado de libertad)``."""
    if domain is None:
        return {}
    return {int(idx): (node.id, dof_name)
            for node in domain.nodes.values() for dof_name, idx in node.dofs.items()}


def _dof_keys(domain, dofs) -> list[tuple[int, str]]:
    """``(id del nodo, nombre del grado de libertad)`` de los DOFs globales dados."""
    index = _dof_index_map(domain)
    return [index[int(d)] for d in dofs]


def _dof_labels(domain, dofs) -> list[str]:
    """Etiquetas ``"nodo <id> (<dof>)"`` para los DOFs globales dados."""
    index = _dof_index_map(domain)
    return [f"nodo {index[int(d)][0]} ({index[int(d)][1]})" if int(d) in index
            else f"DOF {int(d)}" for d in dofs]


def warn_initial_conditions_off_constraints(
    T, g, free_dofs, u0, u0_dot, solver_name: str, domain=None,
) -> None:
    """Avisa (log, no error) si ``u0`` o ``u̇0`` incumplen una MPC con maestros.

    Los solvers transitorios conservan sólo las componentes libres de las
    condiciones iniciales y reconstruyen los DOFs esclavos con el operador
    de eliminación directa (ADR 0004): ``u0 ← T·u0[libres] + g`` y
    ``u̇0 ← T·u̇0[libres]``. El valor que el usuario dé en un esclavo se
    descarta. Si es coherente con la restricción no cambia nada; si no lo es,
    la respuesta depende de qué DOF se declaró esclavo, y
    ``ModalResult.free_vibration`` —que proyecta en norma M— daría otra.
    Este aviso lo hace visible sin cambiar la conducta.

    Sólo se examinan las filas esclavas con maestros (al menos un
    coeficiente no nulo en ``T``). Las filas Dirichlet (sin maestros) quedan
    fuera: su proyección sobre el valor prescrito es deliberada y tiene test
    propio.

    Parameters
    ----------
    T, g : operador y término independiente de ``ConstraintSet.build``.
    free_dofs : DOFs libres (``ConstraintSet.free_dofs``).
    u0, u0_dot : condiciones iniciales globales dadas por el usuario, o
        ``None`` (no se examina).
    solver_name : nombre del solver para el mensaje.
    domain : opcional, para nombrar nodo y DOF en el mensaje.
    """
    if u0 is None and u0_dot is None:
        return
    T_csr = sp.csr_matrix(T, copy=True)
    T_csr.eliminate_zeros()
    ndof = T_csr.shape[0]
    slaves = np.setdiff1d(np.arange(ndof), np.asarray(free_dofs, dtype=np.int64))
    if slaves.size == 0:
        return
    with_masters = slaves[np.diff(T_csr.indptr)[slaves] > 0]
    if with_masters.size == 0:
        return
    free_dofs = np.asarray(free_dofs, dtype=np.int64)
    g = np.asarray(g, dtype=float)
    checks = (
        ("u0", u0, g, "T·u0[libres] + g"),
        ("u0_dot", u0_dot, np.zeros_like(g), "T·u0_dot[libres]"),
    )
    for label, vec, offset, formula in checks:
        if vec is None:
            continue
        vec = np.asarray(vec, dtype=float).reshape(ndof)
        recon = np.asarray(T_csr @ vec[free_dofs], dtype=float) + offset
        scale = max(float(np.max(np.abs(vec))), float(np.max(np.abs(recon))))
        if scale == 0.0:
            continue
        diff = np.abs(vec[with_masters] - recon[with_masters])
        bad = with_masters[diff > INITIAL_CONDITION_MPC_RTOL * scale]
        if bad.size == 0:
            continue
        labels = _dof_labels(domain, bad[:8])
        muestra = ", ".join(labels)
        if bad.size > 8:
            muestra += f", … (+{bad.size - 8})"
        _log.warning(
            "%s: la condición inicial %s incumple la restricción multipunto en "
            "%d DOF esclavo(s) [%s]. El valor dado en el esclavo se descarta y "
            "se reconstruye desde sus maestros como %s (ADR 0004): la "
            "respuesta depende de qué DOF se declaró esclavo, y "
            "ModalResult.free_vibration, que proyecta en norma M, daría otra. "
            "Da condiciones iniciales coherentes con las restricciones.",
            solver_name, label, int(bad.size), muestra, formula,
        )


def massless_error(M_red, solver_name: str, cause: str, *, domain=None,
                   free_dofs=None) -> MasslessDofsError | None:
    """``MasslessDofsError`` que nombra los grados de libertad libres sin masa
    (diagonal nula de ``M_red``) como ``"nodo <id> (<dof>)"``, o ``None`` si
    todos tienen masa. ``free_dofs`` traduce el índice reducido al global."""
    zero = massless_dofs(sp.csr_matrix(M_red))
    if zero.size == 0:
        return None
    global_idx = zero if free_dofs is None else np.asarray(free_dofs)[zero]
    return MasslessDofsError(
        cause=cause, massless=zero,
        massless_dofs=_dof_labels(domain, global_idx),
        n_mass=int(M_red.shape[0] - zero.size), solver=solver_name,
    )


def solve_mass_system(M_solver, M_red, rhs, solver_name: str, *, domain=None,
                      free_dofs=None):
    """``M_red · x = rhs`` con diagnóstico accionable si ``M_red`` es singular.

    Para una ``M_red`` con masa en todos los grados de libertad (los que no
    la tienen los trata :func:`initial_acceleration`). Una ``M_red`` singular
    por otra causa da ``ValueError`` en vez de propagar ``NaN`` en silencio.
    """
    try:
        x = M_solver.solve(M_red, rhs)
    except RuntimeError as exc:
        raise ValueError(
            f"{solver_name}: la matriz de masa reducida es singular "
            f"({exc}). Algún DOF libre carece de masa: revise materiales con "
            f"density = 0.0 o secciones sin inercia rotacional en los "
            f"elementos que lo comparten."
        ) from exc
    if not np.all(np.isfinite(x)):
        raise ValueError(
            f"{solver_name}: la matriz de masa reducida es singular o casi "
            f"singular (la aceleración inicial no es finita). Algún DOF libre "
            f"carece de masa: revise materiales con density = 0.0 o secciones "
            f"sin inercia rotacional en los elementos que lo comparten."
        )
    return x


# Residuo de las filas sin masa en t₀, relativo a la escala de las fuerzas,
# por encima del cual el estado inicial no las equilibra.
_MASSLESS_RESIDUAL_RTOL = 1.0e-8


def initial_acceleration(M_solver, M_red, K_red, rhs, solver_name: str, *,
                         domain=None, free_dofs=None, force_scale: float = 0.0):
    """Aceleración inicial ``ü₀`` de ``M·ü₀ = rhs`` en un integrador implícito,
    con grados de libertad sin masa (decisión del usuario, 2026-11-09; ADR 0024).

    ``rhs = F(t₀) − F_int(u₀) − C·u̇₀`` reducido. Sin grados de libertad sin
    masa es ``M⁻¹·rhs``. Con ellos (diagonal nula de ``M``; con ``M``
    semidefinida, fila y columna nulas), con ``m`` los que tienen masa y ``0``
    los demás:

    - ``üₘ = Mₘₘ⁻¹·rhsₘ``: las filas con masa la determinan;
    - ``ü₀ = −K₀₀⁻¹·K₀ₘ·üₘ``: la aceleración **consistente**. Las filas sin
      masa no tienen inercia, ``K₀·u = F₀`` en todo instante (con Rayleigh,
      ``C₀ = a₁·K₀`` en esas filas), y derivarla da esa relación: los grados
      de libertad sin masa siguen de forma cuasiestática a los que la tienen.
      Es la condensación del análisis modal.

    La elección de ``ü₀`` no cambia los desplazamientos ni las velocidades de
    la integración (la ecuación algebraica fija ``u₀`` en cada paso); sí la
    aceleración que se reporta en los grados de libertad sin masa: con
    ``ü₀ = 0`` arrastraría un error que Newmark no amortigua (medido: 50 % de
    la aceleración de piso con una carga aplicada desde ``t₀``).

    ``K_red`` es la rigidez (tangente) reducida en ``u₀``. Si el estado inicial
    no equilibra las filas sin masa (una carga aplicada de golpe en un nudo
    sin masa), se avisa: esos grados de libertad se ajustan en el primer paso.
    Se lanza :class:`~solidum.math.solvers.diagnostics.MasslessDofsError` si
    ningún grado de libertad tiene masa, si los que no la tienen forman un
    mecanismo (``K₀₀`` singular) o si la masa de los demás es singular.
    """
    M_red = sp.csr_matrix(M_red)
    zero = massless_dofs(M_red)
    if zero.size == 0:
        return solve_mass_system(M_solver, M_red, rhs, solver_name,
                                 domain=domain, free_dofs=free_dofs)
    n = M_red.shape[0]
    mass = np.setdiff1d(np.arange(n), zero)
    global_idx = zero if free_dofs is None else np.asarray(free_dofs)[zero]

    def fail(cause: str) -> MasslessDofsError:
        return MasslessDofsError(cause=cause, massless=zero,
                                 massless_dofs=_dof_labels(domain, global_idx),
                                 n_mass=int(mass.size), solver=solver_name)

    rows = M_red[zero]
    if rows.nnz and np.any(rows.data != 0.0):
        raise fail("indefinite_mass")
    if mass.size == 0:
        raise fail("no_mass")
    K_red = sp.csr_matrix(K_red)
    try:
        lu_zz = spla.splu(K_red[zero][:, zero].tocsc())
    except RuntimeError as exc:
        raise fail("massless_mechanism") from exc
    try:
        a_m = M_solver.solve(M_red[mass][:, mass].tocsr(), rhs[mass])
    except RuntimeError as exc:
        raise fail("singular_mass") from exc
    if not np.all(np.isfinite(a_m)):
        raise fail("singular_mass")
    a = np.empty(n)
    a[mass] = a_m
    a[zero] = -lu_zz.solve(np.asarray(K_red[zero][:, mass] @ a_m).ravel())

    residual = float(np.linalg.norm(rhs[zero]))
    scale = max(float(force_scale), float(np.linalg.norm(rhs)))
    if residual > _MASSLESS_RESIDUAL_RTOL * scale and residual > 0.0:
        bad = zero[np.abs(rhs[zero]) > _MASSLESS_RESIDUAL_RTOL * scale]
        labels = _dof_labels(domain, bad if free_dofs is None else np.asarray(free_dofs)[bad])
        muestra = ", ".join(labels[:8]) + (f", … (+{len(labels) - 8})" if len(labels) > 8 else "")
        _log.warning(
            "%s: en t₀ la carga sobre grados de libertad sin masa no está "
            "equilibrada por el estado inicial [%s]. No tienen inercia: se "
            "ajustan de golpe en el primer paso, y el estado reportado en t₀ "
            "no es de equilibrio en ellos.", solver_name, muestra)
    _log.info("  %d grado(s) de libertad sin masa: aceleración inicial consistente "
              "ü₀ = −K₀₀⁻¹·K₀ₘ·üₘ (ADR 0024).", zero.size)
    return a
