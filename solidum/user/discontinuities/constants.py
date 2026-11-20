"""Constantes del módulo de discontinuidades interiores (ADR 0010).

Newton local del salto cohesivo en ``CST_Embedded2D``: tolerancia relativa
sobre el residuo del salto ``[[u]]`` y máximo de iteraciones. Más estricta
que la del Newton global porque típicamente converge en 1-3 iteraciones
gracias a la rigidez de penalización ``K_e``.

Seguimiento de la grieta (spec ``CrackTracking``): magnitud mínima, relativa
al tamaño del elemento ``h_e = √(2A_e)``, de los valores nodales de la grieta
(§1–§3 y §5); y la inclinación mínima con que la prolongación debe dirigirse a
un lado del borde libre para salir por él: d·n_ext > ``FREE_EDGE_MIN_APPROACH``,
el seno del ángulo entre la prolongación y el lado (0.05 ≈ 2.9°; §3, decisión
del usuario del 2026-11-20). Con > 0, una grieta paralela al borde salía o no
según el redondeo de d·n_ext ≈ 0.
"""
LOCAL_JUMP_RTOL = 1.0e-10
LOCAL_JUMP_MAX_ITER = 30

TRACKING_MIN_DISTANCE = 1.0e-3
FREE_EDGE_MIN_APPROACH = 0.05
