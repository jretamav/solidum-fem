"""Generación de mallas dentro de Solidum (ADR 0022 y 0037).

:func:`structured_mesh` malla un bloque transformado con cualquier elemento
que declare ``REFERENCE_CELL`` y ``REFERENCE_NODES``, con líneas de rejilla
por tramos y gradación (:func:`axis_lines`), celdas excluidas, partición de
triángulos y fusión con bloques anteriores; desde el YAML,
``mesh: {structured: ...}`` (:mod:`solidum.mesh.structured_blocks`).
:func:`mesh_quality` mide la calidad de los elementos y :func:`free_boundary`
da el contorno libre (lados o caras de un solo elemento; ADR 0042). Las geometrías reales
se mallan con Gmsh y se leen con ``solidum.utils.gmsh_parser``; desde el YAML,
``mesh: {geometry: archivo.geo, ...}`` llama a Gmsh con
:mod:`solidum.mesh.gmsh_geometry` (ADR 0035).
"""
from solidum.mesh.boundary import free_boundary
from solidum.mesh.quality import MeshQuality, mesh_quality
from solidum.mesh.structured import (
    DIAGONALS,
    SIDES_1D,
    SIDES_2D,
    SIDES_3D,
    StructuredMesh,
    axis_lines,
    structured_mesh,
)

__all__ = ["DIAGONALS", "MeshQuality", "SIDES_1D", "SIDES_2D", "SIDES_3D", "StructuredMesh",
           "axis_lines", "free_boundary", "mesh_quality", "structured_mesh"]
