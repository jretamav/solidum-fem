"""Elementos térmicos — conducción de calor (Etapa 8).

Familia paralela a los sólidos mecánicos: un DOF escalar ``T`` por nodo,
matriz de conductividad en vez de rigidez y capacidad calorífica en vez de
masa. Comparten los kernels de forma y jacobiano de sus gemelos mecánicos.
"""
from solidum.elements.thermal.quad4_thermal import Quad4Thermal

__all__ = ["Quad4Thermal"]
