"""Imposición de condiciones de frontera por eliminación directa (ADR 0004).

Este subpaquete implementa el modelo afín ``u = T·u_libre + g``: toda
restricción lineal sobre los DOFs se reduce a la forma ``u_s = g_s + Σ α_si·u_mi``
y se elimina del sistema antes de resolverlo.

En la fase 1 sólo se expone Dirichlet (sin maestros). MPC lineales y cierre
transitivo entrarán en la fase 2 cuando aparezca el primer caso de uso (ver
``docs/adr/0004-imposicion-de-condiciones-de-dirichlet.md``).
"""

from fenix.bc.constraints import AffineConstraint, ConstraintSet

__all__ = ["AffineConstraint", "ConstraintSet"]
