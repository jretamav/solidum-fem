"""Benchmarks de validación contra referencias externas (NAFEMS, Lamé, MacNeal-Harder slender beam, etc.).

Esta carpeta contiene tests que validan Fenix FEM contra soluciones analíticas
clásicas y benchmarks publicados, complementando los tests unitarios y de
integración del directorio padre.

Convención: cada archivo ``test_<benchmark>.py`` documenta su referencia
bibliográfica en el módulo-docstring, el valor objetivo, y la tolerancia
justificada.
"""
