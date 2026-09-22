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
"""
import numpy as np

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
    """
    for elem in domain.elements.values():
        if not getattr(elem, "PRESERVES_SYMMETRY", True):
            return False
        material = getattr(elem, "material", None)
        if material is not None and not getattr(material, "IS_SYMMETRIC", True):
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
