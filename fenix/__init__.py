# fenix_fem/fenix/__init__.py

# 1. Inicializar los registros automáticamente al importar la librería
from fenix import registry_initialization
from fenix.registry import MaterialRegistry, ElementRegistry, SolverRegistry

# 2. Exponer el Core
from fenix.core.domain import Domain
from fenix.core.node import Node

# 3. Exponer Materiales
from fenix.materials.elastic import Elastic1D
from fenix.materials.elastic_2d import Elastic2D
from fenix.materials.von_mises_2d import VonMises2D
from fenix.materials.plastic_1d import Elastoplastic1D
from fenix.materials.damage_1d import IsotropicDamage1D
from fenix.materials.damage_2d import IsotropicDamage2D

# 4. Exponer Elementos
from fenix.elements.solid_2d import Quad4, Tri3
from fenix.elements.structural import Truss2D, Truss3D, Frame2DEuler, Frame2DTimoshenko

# 5. Exponer Matemáticas y Solvers
from fenix.math.assembly import Assembler
from fenix.math.solvers import LinearSolver, NonlinearSolver, ArcLengthSolver

# 6. Exponer Utilidades
from fenix.utils.yaml_parser import YamlParser
from fenix.utils.vtk_exporter import VtkExporter
from fenix.utils.gmsh_parser import GmshParser