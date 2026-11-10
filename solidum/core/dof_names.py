# solidum_fem/solidum/core/dof_names.py
"""Nombres de los grados de libertad mecánicos: convención única del programa.

Los elementos declaran sus grados de libertad por nombre (``DOF_NAMES``) y el
resto del programa los reconoce por ese nombre: los modos de cuerpo rígido
(ADR 0018), la descripción de mecanismos (ADR 0019), la exportación VTK y las
masas de nudo (ADR 0024). Antes cada uno guardaba su propia copia de estas
tuplas; aquí viven una sola vez.

Cualquier otro nombre es un campo escalar sin significado mecánico; el del
problema térmico es :data:`TEMPERATURE`.
"""

#: Traslaciones según ``x``, ``y``, ``z``.
TRANSLATIONS: tuple[str, ...] = ("ux", "uy", "uz")
#: Giros alrededor de ``x``, ``y``, ``z`` (regla de la mano derecha).
ROTATIONS: tuple[str, ...] = ("rx", "ry", "rz")
#: Todos los grados de libertad mecánicos.
MECHANICAL: frozenset[str] = frozenset(TRANSLATIONS + ROTATIONS)
#: Orden de las seis componentes de un movimiento (o de una fuerza) en 3D:
#: traslaciones y giros según x, y, z.
MECHANICAL_ORDER: tuple[str, ...] = TRANSLATIONS + ROTATIONS
#: Grado de libertad del problema térmico (Etapa 8).
TEMPERATURE: str = "T"


def translation_dims(nodes) -> int:
    """Dimensión del espacio de traslaciones de un conjunto de nodos: 3 si
    alguno tiene ``uz``, 2 en otro caso."""
    return 3 if any(TRANSLATIONS[2] in n.dofs for n in nodes) else 2
