from fenix.registry import MaterialRegistry, ElementRegistry, SolverRegistry

# Materiales
from fenix.materials.von_mises_2d import VonMises2D
from fenix.materials.elastic_2d import Elastic2D
from fenix.materials.elastic import Elastic1D
from fenix.materials.plastic_1d import Elastoplastic1D
from fenix.materials.damage_2d import IsotropicDamage2D
from fenix.materials.damage_1d import IsotropicDamage1D

MaterialRegistry.register('VonMises2D', VonMises2D)
MaterialRegistry.register('Elastic2D', Elastic2D)
MaterialRegistry.register('Elastic1D', Elastic1D)
MaterialRegistry.register('Elastoplastic1D', Elastoplastic1D)
MaterialRegistry.register('IsotropicDamage2D', IsotropicDamage2D)
MaterialRegistry.register('IsotropicDamage1D', IsotropicDamage1D)

# Elementos
from fenix.elements.solid_2d import Quad4, Tri3
from fenix.elements.structural import Frame2DEuler, Frame2DTimoshenko, Truss2D, Truss3D

ElementRegistry.register('Quad4', Quad4)
ElementRegistry.register('Tri3', Tri3)
ElementRegistry.register('Frame2DEuler', Frame2DEuler)
ElementRegistry.register('Frame2DTimoshenko', Frame2DTimoshenko)
ElementRegistry.register('Truss2D', Truss2D)
ElementRegistry.register('Truss3D', Truss3D)

# Solucionadores
from fenix.math.solvers import LinearSolver, NonlinearSolver

SolverRegistry.register('LinearSolver', LinearSolver)
SolverRegistry.register('NonlinearSolver', NonlinearSolver)