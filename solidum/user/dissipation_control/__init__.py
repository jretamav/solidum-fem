"""Módulo de usuario: control del paso por disipación de energía (ADR 0020).

Arco con restricción de disipación de Gutiérrez (2004), con conmutación
automática cilíndrico↔disipación a la manera de Verhoosel, Remmers y Gutiérrez
(2009): la magnitud que el paso mantiene constante es la energía disipada
incremental, no la longitud de arco euclídea. Es una variante del
``ArcLengthSolver`` del programa principal que sólo cambia la restricción del
paso, y se conecta a él por la clase de solver y por el problema de Newton por
paso (``ArcProblem``) que el programa principal ofrece como punto de extensión.

Por qué vive fuera del programa principal (criterio de admisión, Reglas.md
§1): ningún código general establecido implementa esta restricción de forma
nativa. DIANA ofrece como restricciones del arco el *Spherical Path* y el
*Updated Normal Plane*, y su método «energy based» es una regla de tamaño de
paso; el ``PILOTAGE`` de Code_Aster no la incluye; en Abaqus sólo existen
implementaciones de investigación mediante subrutinas de usuario. Como un
macro de usuario de FEAP (``umacr``), vive en su propio módulo.

Contenido:

- ``DissipationArcLengthSolver`` — el solver (``dissipation_arclength.py``).

Uso: en YAML, ``user_modules: [dissipation_control]`` y
``solver: {type: DissipationArcLengthSolver, initial_tau: ...}``; en Python,
``solidum.load_user_module("dissipation_control")`` o importar de aquí.
"""
from solidum.user.dissipation_control.dissipation_arclength import DissipationArcLengthSolver

__all__ = [
    "DissipationArcLengthSolver",
]
