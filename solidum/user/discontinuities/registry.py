"""Registro de la familia de materiales cohesivos (ADR 0010, ADR 0020).

Familia propia, paralela a los materiales del programa principal: un
cohesivo relaciona la tracción ``t`` con el salto ``[[u]]`` sobre ``Γ_d``, no
el esfuerzo con la deformación en Voigt. Declarar ``YAML_SECTION`` la
convierte en una familia de material que el lector YAML del programa
principal recorre sin conocerla (ADR 0020, P1-P2).
"""
from __future__ import annotations

from typing import Dict, Type

from solidum.registry import Registry


class CohesiveMaterialRegistry(Registry):
    _items: Dict[str, Type] = {}
    _kind = "MaterialCohesivo"
    SPEC_KIND = "cohesive_material"
    YAML_SECTION = "cohesive_materials"
    YAML_LABEL = "material cohesivo"
