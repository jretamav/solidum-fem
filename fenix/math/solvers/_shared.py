# fenix_fem/fenix/math/solvers/_shared.py
"""Utilidades compartidas entre solvers del paquete ``fenix.math.solvers``.

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
from fenix.logging import get_logger


_log = get_logger("solvers")


def domain_is_symmetric(domain) -> bool:
    """``True`` si todos los componentes del dominio preservan simetría de ``K``.

    Recorre elementos y materiales y agrega los flags declarativos
    ``PRESERVES_SYMMETRY`` y ``IS_SYMMETRIC``. Default en ambos es ``True`` —
    sólo plasticidad no asociada, follower loads y similares lo desactivan.
    """
    for elem in domain.elements.values():
        if not getattr(type(elem), "PRESERVES_SYMMETRY", True):
            return False
        material = getattr(elem, "material", None)
        if material is not None and not getattr(type(material), "IS_SYMMETRIC", True):
            return False
    return True


try:
    from fenix.math.linalg.cholesky import CholeskyNotPositiveDefiniteError
except ImportError:
    class CholeskyNotPositiveDefiniteError(Exception):
        """Placeholder cuando ``scikit-sparse`` no está instalado (nunca se lanza)."""
