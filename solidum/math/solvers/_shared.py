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
from solidum.logging import get_logger


_log = get_logger("solvers")


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

    Sólo participan los elementos cuyo ``internal_forces`` no devuelve
    ``None`` (estructurales 1D, ADR 0002); se sondea en el primer registro.
    """

    def __init__(self, domain, enabled: bool):
        self.domain = domain
        self.enabled = bool(enabled)
        self._history: dict[int, list] = {}
        self._elements: list | None = None

    def record(self, U: np.ndarray) -> None:
        """Añade el paso actual; llamar tras comitear el estado del paso."""
        if not self.enabled:
            return
        U = np.asarray(U, dtype=float)
        if self._elements is None:
            self._elements = []
            for elem_id, elem in self.domain.elements.items():
                forces = elem.internal_forces(U)
                if forces is None:
                    continue
                self._elements.append((elem_id, elem))
                self._history[elem_id] = [forces]
            return
        for elem_id, elem in self._elements:
            self._history[elem_id].append(elem.internal_forces(U))

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


def solve_mass_system(M_solver, M_red, rhs, solver_name: str):
    """``M_red · x = rhs`` con diagnóstico accionable si ``M_red`` es singular.

    Una fila nula en ``M_red`` (DOF libre sin masa: material con
    ``density = 0.0`` en todos los elementos que lo comparten, o inercia
    rotacional nula) hace que el backend devuelva ``NaN``/``inf`` o lance
    ``RuntimeError``. Aquí se convierte en ``ValueError`` que apunta a la
    causa, en vez de propagar un historial de ``NaN`` en silencio.
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
