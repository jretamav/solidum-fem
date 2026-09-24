"""Módulo de usuario: discontinuidades interiores embebidas (ADR 0010, ADR 0020).

Formulación de J. Retama (2010), *Formulation and Approximation to Problems in
Solids by Embedded Discontinuity Models*, UNAM: triángulo CST con un salto de
desplazamientos embebido en su interior (formulación KOS), condensado a nivel
de elemento, y una ley cohesiva *tracción-salto* que gobierna la superficie de
discontinuidad. No es un elemento finito estándar, así que vive fuera del
programa principal, como un elemento de usuario de FEAP.

Contenido:

- ``CST_Embedded2D`` — el elemento (``embedded_cst.py``).
- ``CohesiveMaterial`` — contrato de la familia cohesiva
  (``cohesive_material.py``); ``CohesiveMaterialRegistry`` — su registro,
  con la sección YAML ``cohesive_materials`` (``registry.py``).
- ``CohesiveDamageIsotropic`` — ley cohesiva de daño isótropo con
  penalización (``damage_isotropic.py``).
- ``DiscontinuityState`` — geometría y estado del salto de un elemento
  agrietado (``discontinuity_state.py``).

Uso: en YAML, ``user_modules: [discontinuities]``; en Python,
``solidum.load_user_module("discontinuities")`` o importar de aquí.
"""
from solidum.user.discontinuities.registry import CohesiveMaterialRegistry
from solidum.user.discontinuities.cohesive_material import CohesiveMaterial
from solidum.user.discontinuities.discontinuity_state import DiscontinuityState
from solidum.user.discontinuities.damage_isotropic import CohesiveDamageIsotropic
from solidum.user.discontinuities.embedded_cst import CST_Embedded2D

__all__ = [
    "CST_Embedded2D",
    "CohesiveDamageIsotropic",
    "CohesiveMaterial",
    "CohesiveMaterialRegistry",
    "DiscontinuityState",
]
