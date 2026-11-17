"""Configuración común de las pruebas de Solidum FEM.

Las pruebas corren repartidas entre los núcleos (pytest-xdist, ``-n auto`` en
``pyproject.toml``). Cada proceso usa pocos hilos internos: con 24 procesos y
los hilos por omisión de MKL, OpenMP y Numba compiten entre sí, y la suite
tardaba 4 min 22 s frente a 2 min 45 s con este límite (medido el 2026-11-17,
3550 pruebas, los mismos resultados). Sin xdist (``-n 0``) no se toca nada.

``PYTEST_XDIST_WORKER`` lo fija xdist en cada proceso antes de leer este
archivo, y NumPy y Numba leen estas variables al importarse, después.
"""
import os

if os.environ.get("PYTEST_XDIST_WORKER"):
    for _var, _n in (("OMP_NUM_THREADS", "1"), ("MKL_NUM_THREADS", "1"),
                     ("NUMBA_NUM_THREADS", "2")):
        os.environ.setdefault(_var, _n)
