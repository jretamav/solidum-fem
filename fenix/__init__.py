# fenix_fem/fenix/__init__.py

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

from fenix.materials.elastic import Elastic1D
from fenix.materials.elastic_2d import Elastic2D
from fenix.materials.von_mises_2d import VonMises2D
from fenix.materials.plastic_1d import Elastoplastic1D
from fenix.materials.damage_1d import IsotropicDamage1D
from fenix.materials.damage_2d import IsotropicDamage2D

from fenix.elements.solid_2d import Quad4, Tri3
from fenix.elements.structural import Truss2D, Truss3D, Frame2DEuler, Frame2DTimoshenko

from fenix.math.assembly import Assembler
from fenix.math.solvers import LinearSolver, NonlinearSolver, ArcLengthSolver

from fenix.utils.yaml_parser import YamlParser
from fenix.utils.vtk_exporter import VtkExporter
from fenix.utils.gmsh_parser import GmshParser
