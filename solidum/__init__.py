# fenix_fem/solidum/__init__.py

# 0. Logging configurable (ADR 0005). Inicializa el logger raíz "solidum" antes
#    de cualquier otro import, para que los módulos pidan sus loggers hijos
#    sobre una configuración ya consistente.
from solidum.logging import get_logger, set_log_level

# 1. Auto-descubrimiento: importa todos los materiales/elementos/solvers/cuadraturas,
#    activando los decoradores @Registry.register que los inscriben automáticamente.
#    Añadir un archivo nuevo en solidum/materials/ o solidum/elements/ es suficiente:
#    NO hace falta tocar este archivo ni mantener listas de imports manuales.
from solidum.autodiscover import initialize as _initialize_registries
_initialize_registries()

from solidum.registry import (
    MaterialRegistry,
    CohesiveMaterialRegistry,
    ElementRegistry,
    SolverRegistry,
    QuadratureRegistry,
)

# 2. Re-exports de API pública (conveniencia para `from solidum import Foo`).
#    El registro YAML funciona aunque un material/elemento nuevo NO se exporte aquí;
#    esta lista es solo azúcar para uso programático directo.
from solidum.core.domain import Domain
from solidum.core.node import Node
from solidum.results import ElementForces, SolveResult, build_solve_result
from solidum.entry import run, run_yaml

from solidum.materials.elastic import Elastic1D
from solidum.materials.elastic_2d import Elastic2D
from solidum.materials.von_mises_2d import VonMises2D
from solidum.materials.plastic_1d import Elastoplastic1D
from solidum.materials.damage_1d import IsotropicDamage1D
from solidum.materials.damage_2d import IsotropicDamage2D
from solidum.materials.cable_1d import CableMaterial1D

from solidum.cohesive_materials.damage_isotropic import CohesiveDamageIsotropic

from solidum.elements.solid_2d import Quad4, Quad8, Quad9, Tri3, Tri6
from solidum.elements.solid_2d.embedded_cst import CST_Embedded2D
from solidum.elements.truss import Truss2D, Truss2DCorot, Truss3D, Truss3DCorot
from solidum.elements.cable import Cable2DCorot, Cable3DCorot
from solidum.elements.frame import Frame2DEuler, Frame2DEulerCorot, Frame2DTimoshenko
from solidum.elements.frame3d import Frame3D

from solidum.math.assembly import Assembler
from solidum.math.solvers import LinearSolver, NonlinearSolver, ArcLengthSolver

from solidum.utils.yaml_parser import YamlParser
from solidum.utils.vtk_exporter import VtkExporter
from solidum.utils.gmsh_parser import GmshParser
