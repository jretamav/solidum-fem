# solidum_fem/solidum/core/dof_names.py
"""Nombres de los grados de libertad mecánicos: convención única del programa.

Los elementos declaran sus grados de libertad por nombre (``DOF_NAMES``) y el
resto del programa los reconoce por ese nombre: los modos de cuerpo rígido
(ADR 0018), la descripción de mecanismos (ADR 0019), la exportación VTK y las
masas de nudo (ADR 0024). Antes cada uno guardaba su propia copia de estas
tuplas; aquí viven una sola vez.

Cualquier otro nombre (``T`` en el problema térmico) es un campo escalar sin
significado mecánico.
"""

#: Traslaciones según ``x``, ``y``, ``z``.
TRANSLATIONS: tuple[str, ...] = ("ux", "uy", "uz")
#: Giros alrededor de ``x``, ``y``, ``z`` (regla de la mano derecha).
ROTATIONS: tuple[str, ...] = ("rx", "ry", "rz")
#: Todos los grados de libertad mecánicos.
MECHANICAL: frozenset[str] = frozenset(TRANSLATIONS + ROTATIONS)
