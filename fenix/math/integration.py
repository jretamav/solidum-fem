# fenix_fem/fenix/math/integration.py
import numpy as np
from fenix.registry import QuadratureRegistry

class GaussQuadrature:
    @staticmethod
    def get_points_2d_1x1():
        points = [(0.0, 0.0)]
        weights = [4.0]
        return points, weights

    @staticmethod
    def get_points_2d_2x2():
        pt = 1.0 / np.sqrt(3.0)
        points = [(-pt, -pt), ( pt, -pt), ( pt,  pt), (-pt,  pt)]
        weights = [1.0, 1.0, 1.0, 1.0]
        return points, weights

# Registrar reglas automáticamente
QuadratureRegistry.register("1x1", *GaussQuadrature.get_points_2d_1x1())
QuadratureRegistry.register("2x2", *GaussQuadrature.get_points_2d_2x2())
