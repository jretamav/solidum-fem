"""Cubo Lamé 3D — solución cerrada de elasticidad lineal isótropa 3D.

Referencias
-----------
Timoshenko, S.P.; Goodier, J.N. (1970). *Theory of Elasticity*, 3rd ed.,
McGraw-Hill, §6 (elasticidad 3D). Cubo unitario bajo tracción uniaxial y
compresión hidrostática son las dos pruebas analíticas más básicas.

Concepto
--------
Cubo unitario de material elástico isótropo con dos escenarios:

1. **Tracción uniaxial**. Cara ``x=1`` con tracción ``(p, 0, 0)``; cara ``x=0``
   con ``u_x = 0`` (y restricciones mínimas que evitan movimiento rígido).
   Solución cerrada:

   .. math::

       u_x = \\frac{p\\,x}{E},\\quad u_y = -\\frac{\\nu p\\,y}{E},\\quad u_z = -\\frac{\\nu p\\,z}{E}

       \\sigma_{xx} = p,\\ \\sigma_{yy} = \\sigma_{zz} = 0,\\ \\sigma_{xy} = \\sigma_{yz} = \\sigma_{xz} = 0

2. **Compresión hidrostática**. Tracción normal ``-p`` en las seis caras
   (con BCs simétricas mínimas). Solución:

   .. math::

       \\sigma_{xx} = \\sigma_{yy} = \\sigma_{zz} = -p,\\quad \\varepsilon_{xx} = \\varepsilon_{yy} = \\varepsilon_{zz} = -\\frac{p}{3K}

   con ``K = E / [3(1-2ν)]`` el módulo de compresibilidad.

En ambos escenarios el campo de desplazamientos es **lineal en (x, y, z)**.
Tanto ``Hex8`` (trilineal Q₁) como ``Tet4`` (lineal P₁) reproducen campos
lineales **exactamente**, así que la solución FEM debe coincidir con la
analítica a **precisión máquina** independientemente del refinamiento.
Esta propiedad sale del patch test (Reglas.md §6, "tests blindan física").

Geometría y carga
-----------------
- Cubo unitario [0, 1]³.
- E = 1000, ν = 0.3, p = 10.
- Caso A (uniaxial Hex8): un Hex8; cara ``x=0`` con ``u_x = 0``; nodo ``(0,0,0)``
  con ``u_y = u_z = 0``; nodo ``(0,1,0)`` con ``u_z = 0``. Tracción ``(p, 0, 0)``
  en cara ``+ξ`` (face index 3).
- Caso B (uniaxial Tet4, malla 5 tetraedros). Mismas BCs.
- Caso C (hidrostática Hex8). Símetria por planos: ``u_x = 0`` en ``x=0``,
  ``u_y = 0`` en ``y=0``, ``u_z = 0`` en ``z=0``. Tracción normal ``-p``
  en cada una de las 3 caras opuestas (face 1, 3, 4).
"""
import os
import sys

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from fenix.core.domain import Domain
from fenix.elements.solid_3d import Hex8, Tet4
from fenix.materials.elastic_3d import Elastic3D
from fenix.math.assembly import Assembler
from fenix.math.solvers import LinearSolver


E_YOUNG = 1000.0
NU = 0.3
PRESSURE = 10.0

# División canónica de un cubo unitario en 5 tetraedros con volumen total = 1.
# Verificada analíticamente: T3 tiene V=1/3, los otros cuatro V=1/6 cada uno.
TET5_CONNECTIVITY = (
    (0, 1, 2, 5),
    (0, 2, 3, 7),
    (0, 2, 7, 5),
    (0, 5, 7, 4),
    (2, 7, 5, 6),
)
HEX8_VTK_COORDS = [
    (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
]


def _apply_face_traction(domain, elem, face, t_vec, F):
    f_elem = elem.compute_face_traction(face, t_vec)
    for k_local, k_global in enumerate(elem.get_global_dof_indices()):
        if k_global >= 0:
            F[k_global] += f_elem[k_local]


def _solve_static(domain):
    domain.generate_equation_numbers(verbose=False)
    assembler = Assembler(domain)
    F = np.zeros(domain.total_dofs)
    return assembler, F


def _node_id_by_coord(domain, x, y, z, tol=1e-12):
    for n in domain.nodes.values():
        if (abs(n.coordinates[0] - x) < tol
                and abs(n.coordinates[1] - y) < tol
                and abs(n.coordinates[2] - z) < tol):
            return n.id
    raise KeyError(f"Nodo con coordenadas ({x},{y},{z}) no encontrado.")


# =============================================================================
# Helpers — construcción de mallas y BCs.
# =============================================================================

def _build_cube_hex8():
    domain = Domain()
    for i, c in enumerate(HEX8_VTK_COORDS):
        domain.add_node(i + 1, list(c))
    nodes = [domain.nodes[i + 1] for i in range(8)]
    material = Elastic3D(E=E_YOUNG, nu=NU)
    elem = Hex8(1, nodes, material)
    domain.add_element(elem)
    return domain, elem


def _build_cube_tet5():
    domain = Domain()
    for i, c in enumerate(HEX8_VTK_COORDS):
        domain.add_node(i + 1, list(c))
    material = Elastic3D(E=E_YOUNG, nu=NU)
    elems = []
    for k, conn in enumerate(TET5_CONNECTIVITY):
        nodes = [domain.nodes[i + 1] for i in conn]
        e = Tet4(k + 1, nodes, material)
        domain.add_element(e)
        elems.append(e)
    return domain, elems


def _apply_uniaxial_bcs(domain):
    """Restringe la cara ``x=0`` (4 nodos) en u_x; nodo (0,0,0) en (u_y, u_z);
    nodo (0,1,0) en u_z. Suficiente para fijar los 6 modos rígidos sin
    sobre-restringir."""
    for n in domain.nodes.values():
        if abs(n.coordinates[0]) < 1e-12:
            n.fix_dof('ux', 0.0)
    domain.nodes[1].fix_dof('uy', 0.0)
    domain.nodes[1].fix_dof('uz', 0.0)
    domain.nodes[4].fix_dof('uz', 0.0)  # nodo (0,1,0) = id 4


def _apply_symmetry_bcs(domain):
    """Para hidrostática: simetría por planos coordenados.

    u_x = 0 en x=0, u_y = 0 en y=0, u_z = 0 en z=0.
    """
    for n in domain.nodes.values():
        if abs(n.coordinates[0]) < 1e-12:
            n.fix_dof('ux', 0.0)
        if abs(n.coordinates[1]) < 1e-12:
            n.fix_dof('uy', 0.0)
        if abs(n.coordinates[2]) < 1e-12:
            n.fix_dof('uz', 0.0)


# =============================================================================
# Test 1 — Tracción uniaxial Hex8 (1 elemento)
# =============================================================================

def test_uniaxial_traction_hex8_exact():
    """Un Hex8 reproduce la solución de tracción uniaxial exacta a precisión máquina."""
    domain, elem = _build_cube_hex8()
    _apply_uniaxial_bcs(domain)
    assembler, F = _solve_static(domain)
    _apply_face_traction(domain, elem, face=3, t_vec=np.array([PRESSURE, 0.0, 0.0]), F=F)

    U = LinearSolver(assembler).solve(F)

    # Verificación nodal: u_x(x) = p·x/E, u_y(y) = -νp·y/E, u_z(z) = -νp·z/E.
    for n in domain.nodes.values():
        x, y, z = n.coordinates
        ix, iy, iz = n.dofs['ux'], n.dofs['uy'], n.dofs['uz']
        ux = 0.0 if ix < 0 else U[ix]
        uy = 0.0 if iy < 0 else U[iy]
        uz = 0.0 if iz < 0 else U[iz]
        np.testing.assert_allclose(ux, PRESSURE * x / E_YOUNG, atol=1.0e-10)
        np.testing.assert_allclose(uy, -NU * PRESSURE * y / E_YOUNG, atol=1.0e-10)
        np.testing.assert_allclose(uz, -NU * PRESSURE * z / E_YOUNG, atol=1.0e-10)

    # Verificación de esfuerzos en puntos de Gauss.
    gs = elem.compute_gauss_state(U)
    expected_sigma = np.array([PRESSURE, 0.0, 0.0, 0.0, 0.0, 0.0])
    for k in range(gs['stress'].shape[0]):
        np.testing.assert_allclose(gs['stress'][k], expected_sigma, atol=1.0e-9)


# =============================================================================
# Test 2 — Tracción uniaxial Tet4 (malla de 5 tetraedros)
# =============================================================================

def test_uniaxial_traction_tet4_mesh_exact():
    """Cubo en 5 tetraedros bajo tracción uniaxial — solución exacta nodal.

    Los Tet4 son CST: reproducen campos lineales exactamente. La tracción en
    la cara ``x=1`` se aplica recorriendo cada tetraedro cuya cara coincida
    con el plano ``x=1`` (cara con los 3 nodos en x=1).
    """
    domain, elems = _build_cube_tet5()
    _apply_uniaxial_bcs(domain)
    assembler, F = _solve_static(domain)

    # Detectar caras que coincidan con x=1 y aplicar tracción.
    for e in elems:
        for face_idx, face_local in enumerate(e.FACE_NODES):
            xs = [e.nodes[i].coordinates[0] for i in face_local]
            if all(abs(x - 1.0) < 1e-12 for x in xs):
                _apply_face_traction(
                    domain, e, face=face_idx,
                    t_vec=np.array([PRESSURE, 0.0, 0.0]), F=F,
                )

    U = LinearSolver(assembler).solve(F)

    for n in domain.nodes.values():
        x, y, z = n.coordinates
        ix, iy, iz = n.dofs['ux'], n.dofs['uy'], n.dofs['uz']
        ux = 0.0 if ix < 0 else U[ix]
        uy = 0.0 if iy < 0 else U[iy]
        uz = 0.0 if iz < 0 else U[iz]
        np.testing.assert_allclose(ux, PRESSURE * x / E_YOUNG, atol=1.0e-10)
        np.testing.assert_allclose(uy, -NU * PRESSURE * y / E_YOUNG, atol=1.0e-10)
        np.testing.assert_allclose(uz, -NU * PRESSURE * z / E_YOUNG, atol=1.0e-10)

    # Esfuerzos: σ_xx = p, resto ≈ 0 en todos los elementos.
    expected_sigma = np.array([PRESSURE, 0.0, 0.0, 0.0, 0.0, 0.0])
    for e in elems:
        gs = e.compute_gauss_state(U)
        np.testing.assert_allclose(gs['stress'][0], expected_sigma, atol=1.0e-8)


# =============================================================================
# Test 3 — Compresión hidrostática Hex8
# =============================================================================

def test_hydrostatic_compression_hex8():
    """Tracción normal -p en las tres caras opuestas a los planos de simetría.

    Solución analítica:
        σ_ii = -p,  ε_ii = -p / (3K),  con K = E/[3(1-2ν)].
    """
    domain, elem = _build_cube_hex8()
    _apply_symmetry_bcs(domain)
    assembler, F = _solve_static(domain)

    # Tracción normal -p en cara +ξ (3), +η (4), +ζ (1).
    _apply_face_traction(domain, elem, face=3, t_vec=np.array([-PRESSURE, 0.0, 0.0]), F=F)
    _apply_face_traction(domain, elem, face=4, t_vec=np.array([0.0, -PRESSURE, 0.0]), F=F)
    _apply_face_traction(domain, elem, face=1, t_vec=np.array([0.0, 0.0, -PRESSURE]), F=F)

    U = LinearSolver(assembler).solve(F)

    # ε_ii = -p · (1-2ν) / E  (cada componente)
    eps_axial = -PRESSURE * (1.0 - 2.0 * NU) / E_YOUNG
    for n in domain.nodes.values():
        x, y, z = n.coordinates
        ix, iy, iz = n.dofs['ux'], n.dofs['uy'], n.dofs['uz']
        ux = 0.0 if ix < 0 else U[ix]
        uy = 0.0 if iy < 0 else U[iy]
        uz = 0.0 if iz < 0 else U[iz]
        # u_i = ε_axial · x_i por simetría
        np.testing.assert_allclose(ux, eps_axial * x, atol=1.0e-10)
        np.testing.assert_allclose(uy, eps_axial * y, atol=1.0e-10)
        np.testing.assert_allclose(uz, eps_axial * z, atol=1.0e-10)

    gs = elem.compute_gauss_state(U)
    expected_sigma = np.array([-PRESSURE, -PRESSURE, -PRESSURE, 0.0, 0.0, 0.0])
    for k in range(gs['stress'].shape[0]):
        np.testing.assert_allclose(gs['stress'][k], expected_sigma, atol=1.0e-8)
