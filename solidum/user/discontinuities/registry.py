"""Registro de la familia de materiales cohesivos (ADR 0010, ADR 0020).

Familia propia, paralela a los materiales del programa principal: un
cohesivo relaciona la tracción ``t`` con el salto ``[[u]]`` sobre ``Γ_d``, no
el esfuerzo con la deformación en Voigt. Declarar ``YAML_SECTION`` la
convierte en una familia de material que el lector YAML del programa
principal recorre sin conocerla (ADR 0020, P1-P2).

Sólo declara atributos públicos de :class:`solidum.registry.Registry`: el
almacén de clases lo crea la base al definir la subclase.

La familia vive aquí porque hoy su único consumidor es ``CST_Embedded2D``. Una
ley tracción-salto genérica con elementos de interfaz es estándar según el
criterio de admisión (Reglas.md §1); si entra al programa principal un
elemento de interfaz que la use, la familia pasará con él.
"""
from __future__ import annotations

from solidum.registry import Registry


class CohesiveMaterialRegistry(Registry):
    SPEC_KIND = "cohesive_material"
    YAML_SECTION = "cohesive_materials"
    YAML_LABEL = "material cohesivo"
