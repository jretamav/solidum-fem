# fenix_fem/fenix/__init__.py

# 0. Logging configurable (ADR 0005). Inicializa el logger raíz "fenix" antes
#    de cualquier otro import, para que los módulos pidan sus loggers hijos
#    sobre una configuración ya consistente.
from fenix.logging import get_logger, set_log_level

# 1. Auto-descubrimiento: importa todos los materiales/elementos/solvers/cuadraturas,
#    activando los decoradores @Registry.register que los inscriben automáticamente.
#    Añadir un archivo nuevo en fenix/materials/ o fenix/elements/ es suficiente:
#    NO hace falta tocar este archivo ni mantener listas de imports manuales.
from fenix.autodiscover import initialize as _initialize_registries
_initialize_registries()

from fenix.registry import MaterialRegistry, ElementRegistry, SolverRegistry, QuadratureRegistry

# 2. Re-exports de API pública (conveniencia para `from fenix import Foo`).
#    El registro YAML funciona aunque un material/elemento nuevo NO se exporte aquí;
#    esta lista es solo azúcar para uso programático directo.
from fenix.core.domain import Domain
from fenix.core.node import Node
from fenix.results import ElementForces, SolveResult, build_solve_result
from fenix.entry import run, run_yaml

from fenix.materials.elastic import Elastic1D
from fenix.materials.elastic_2d import Elastic2D
from fenix.materials.von_mises_2d import VonMises2D
from fenix.materials.plastic_1d import Elastoplastic1D
from fenix.materials.damage_1d import IsotropicDamage1D
from fenix.materials.damage_2d import IsotropicDamage2D
from fenix.materials.cable_1d import CableMaterial1D

from fenix.elements.solid_2d import Quad4, Tri3
from fenix.elements.truss import Truss2D, Truss2DCorot, Truss3D, Truss3DCorot
from fenix.elements.cable import Cable2DCorot, Cable3DCorot
from fenix.elements.frame import Frame2DEuler, Frame2DEulerCorot, Frame2DTimoshenko
from fenix.elements.frame3d import Frame3D

from fenix.math.assembly import Assembler
from fenix.math.solvers import LinearSolver, NonlinearSolver, ArcLengthSolver

from fenix.utils.yaml_parser import YamlParser
from fenix.utils.vtk_exporter import VtkExporter
from fenix.utils.gmsh_parser import GmshParser
