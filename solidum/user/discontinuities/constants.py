"""Constantes del módulo de discontinuidades interiores (ADR 0010).

Newton local del salto cohesivo en ``CST_Embedded2D``: tolerancia relativa
sobre el residuo del salto ``[[u]]`` y máximo de iteraciones. Más estricta
que la del Newton global porque típicamente converge en 1-3 iteraciones
gracias a la rigidez de penalización ``K_e``.
"""
LOCAL_JUMP_RTOL = 1.0e-10
LOCAL_JUMP_MAX_ITER = 30
