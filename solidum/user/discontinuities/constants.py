"""Constantes del módulo de discontinuidades interiores (ADR 0010).

Newton local del salto cohesivo en ``CST_Embedded2D``: tolerancia relativa
sobre el residuo del salto ``[[u]]`` y máximo de iteraciones. Más estricta
que la del Newton global porque típicamente converge en 1-3 iteraciones
gracias a la rigidez de penalización ``K_e``.

Seguimiento global de la grieta (spec ``GlobalCrackTracking``):
regularización adimensional ``ε`` del problema auxiliar (§2) y magnitud
mínima, relativa al tamaño del elemento, de los valores de Dirichlet del
campo de seguimiento (§2 y §5).
"""
LOCAL_JUMP_RTOL = 1.0e-10
LOCAL_JUMP_MAX_ITER = 30

TRACKING_EPSILON = 1.0e-6
TRACKING_MIN_DISTANCE = 1.0e-3
