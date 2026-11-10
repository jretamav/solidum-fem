"""Módulo de usuario: discontinuidades interiores embebidas (ADR 0010, ADR 0020).

Formulación de J. Retama (2010), *Formulation and Approximation to Problems in
Solids by Embedded Discontinuity Models*, UNAM: triángulo CST con un salto de
desplazamientos embebido en su interior (formulación KOS), condensado a nivel
de elemento, y una ley cohesiva *tracción-salto* que gobierna la superficie de
discontinuidad. No es un elemento finito estándar, así que vive fuera del
programa principal, como un elemento de usuario de FEAP.

Contenido:

- ``CST_Embedded2D`` — el elemento (``embedded_cst.py``).
- La ley cohesiva de la superficie es un material de interfaz del programa
  principal (sección YAML ``interface_materials``; p. ej.
  ``solidum.materials.interface.CohesiveDamageIsotropic``; ADR 0033).
- ``DiscontinuityState`` — geometría y estado del salto de un elemento
  agrietado (``discontinuity_state.py``).
- ``GlobalCrackTracking`` — seguimiento global de la trayectoria de grieta,
  que decide qué elementos se agrietan y su nodo solitario
  (``crack_tracking.py``); ``CrackTrackerRegistry`` — su registro, con la
  sección YAML ``crack_trackers``. ``PrescribedCrackPath`` — variante con la
  trayectoria dada (grieta prescrita).

Uso: en YAML, ``user_modules: [discontinuities]``; en Python,
``solidum.load_user_module("discontinuities")`` o importar de aquí.
"""
from solidum.user.discontinuities.registry import CrackTrackerRegistry
from solidum.user.discontinuities.discontinuity_state import DiscontinuityState
from solidum.user.discontinuities.crack_tracking import (
    Crack,
    GlobalCrackTracking,
    PrescribedCrackPath,
)
from solidum.user.discontinuities.embedded_cst import CST_Embedded2D

__all__ = [
    "CST_Embedded2D",
    "Crack",
    "CrackTrackerRegistry",
    "DiscontinuityState",
    "GlobalCrackTracking",
    "PrescribedCrackPath",
]
