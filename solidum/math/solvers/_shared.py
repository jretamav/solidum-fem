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

from solidum.constants import INITIAL_CONDITION_MPC_RTOL
from solidum.core.load_vector import LoadVector, equivalent_loads_of
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
    """

    termination: Termination | None = None
    steps_done: int = 0
    t_final: float = 0.0
    U_final: np.ndarray | None = None
    partial_result = None

    def _begin_trace(self) -> None:
        self.termination = None
        self.steps_done = 0
        self.t_final = 0.0
        self.U_final = None
        self.partial_result = None

    @property
    def completed(self) -> bool:
        return self.termination == Termination.T_END


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
        """Marca el paso ``step`` (instante ``t``) como convergido."""
        self.step = int(step)
        self.solver.steps_done = int(step)
        self.solver.t_final = float(t)
        self._U = U

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


def _dof_labels(domain, dofs) -> list[str]:
    """Etiquetas ``"nodo <id> (<dof>)"`` para los DOFs globales dados."""
    wanted = {int(d) for d in dofs}
    names: dict[int, str] = {}
    if domain is not None:
        for node in domain.nodes.values():
            for dof_name, idx in node.dofs.items():
                if int(idx) in wanted:
                    names[int(idx)] = f"nodo {node.id} ({dof_name})"
    return [names.get(int(d), f"DOF {int(d)}") for d in dofs]


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

    Un grado de libertad libre sin masa (diagonal nula: sin densidad en los
    elementos que lo comparten ni masa de nudo, o sin inercia rotacional)
    deja la aceleración inicial indefinida: se lanza ``MasslessDofsError``
    (``cause = "initial_acceleration"``) con esos grados de libertad
    nombrados, en vez de propagar un historial de ``NaN`` en silencio. Una
    ``M_red`` singular por otra causa da ``ValueError``.
    """
    error = massless_error(M_red, solver_name, "initial_acceleration",
                           domain=domain, free_dofs=free_dofs)
    if error is not None:
        raise error
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
