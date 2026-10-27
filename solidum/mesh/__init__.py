"""Generación de mallas dentro de Solidum (ADR 0022).

:func:`structured_mesh` malla un bloque transformado con cualquier elemento
que declare ``REFERENCE_CELL`` y ``REFERENCE_NODES``. Las geometrías reales se
mallan con Gmsh y se leen con ``solidum.utils.gmsh_parser``.
"""
from solidum.mesh.structured import SIDES_2D, SIDES_3D, StructuredMesh, structured_mesh

__all__ = ["SIDES_2D", "SIDES_3D", "StructuredMesh", "structured_mesh"]
