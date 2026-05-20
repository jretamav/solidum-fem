"""NAFEMS LE1 — Elliptic Membrane (plane stress).

Referencia
----------
NAFEMS (1990). *The Standard NAFEMS Benchmarks*, Issue 1, Revision 3,
test 11.1 ("Elliptic Membrane"). Glasgow, UK.

Definición
----------
Cuadrante de placa entre dos elipses concéntricas, plane stress:

- Elipse exterior: x²/3.25² + y²/2.75² = 1
- Elipse interior: x²/2²    + y²/1²    = 1
- Espesor t = 0.1
- Material: E = 210 000, ν = 0.3

Puntos canónicos del contorno:

    A = (0, 1)        — extremo superior de la elipse interior
    B = (3.25, 0)     — extremo derecho de la elipse exterior
    C = (0, 2.75)     — extremo superior de la elipse exterior
    D = (2, 0)        — extremo derecho de la elipse interior

Condiciones de borde
--------------------
- Borde AD (x = 0, y ∈ [1, 2.75])  → u_x = 0  (simetría)
- Borde DC (y = 0, x ∈ [2, 3.25])  → u_y = 0  (simetría)
- Borde BC (elipse exterior)        → tracción normal saliente σ_n = 10
  (presión exterior de 10 unidades aplicada como tracción del material en la
  dirección de la normal exterior a la elipse exterior).
- Borde AB (elipse interior)        → libre.

Cantidad objetivo
-----------------
**σ_yy en el punto D = (2, 0) = 92.7** (valor canónico NAFEMS).

Es la tensión tangencial sobre la elipse interior — en D la tangente a la
elipse interior es vertical, así que la "tangencial" coincide con la
componente y global. Esta tensión es relativamente sensible a:

1. Aplicación correcta de la presión normal variable sobre el borde curvo
   exterior.
2. Cinemática del elemento ante un campo de tensiones no uniforme.
3. Calidad de la malla en la zona próxima al punto interior D.

Stress recovery en el punto D
-----------------------------
σ_yy(D) es un **valor extremo en un nodo de esquina**. El "valor numérico"
de σ en el punto exacto D no es lo que devuelve directamente el FEM —
la formulación entrega σ en los puntos de Gauss. Hay tres estrategias de
recovery posibles:

1. **Promedio de los Gauss del elemento que contiene a D** — subestima
   sistemáticamente (σ_yy varía de 35 a 77 dentro del elemento en una
   malla coarse; el promedio queda en ~48 vs 92.7 esperado).
2. **σ en el Gauss más cercano a D** — converge limpio para Q4/Q8/Q9 con
   refinamiento (el Gauss se acerca al nodo a tasa h). Para Tri3/Tri6 el
   "Gauss más cercano" está en el centroide del triángulo, no cerca del
   vértice, y converge mucho más lento.
3. **Extrapolación de Barlow / SPR de Zienkiewicz** — recovery
   "superconvergente" que devolvería ≈92.7 con malla coarse. **No está
   implementado en post-proceso de Solidum.** Sería el siguiente paso si el
   benchmark se quiere afinar.

Este test usa **opción 2** (σ en el Gauss más cercano) y verifica:

- **Cota cuantitativa** por elemento y malla, calibrada empíricamente.
- **Convergencia h**: σ_yy(Gauss_más_cercano) → 92.7 al refinar la malla.

Las cotas son honestas con la limitación del recovery: Q4 32×32 da
≈90 (3% error); Q8 16×16 da ≈89 (4%); Tri3 32×32 da ≈78 (16%);
Tri6 16×16 da ≈78 (15%). La convergencia es monótona en todos los casos.

**Implementar SPR como deuda técnica está documentado en el README de
validación.**
"""
import os
import sys

import numpy as np
import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from solidum.core.domain import Domain
from solidum.elements.solid_2d import Quad4, Quad8, Quad9, Tri3, Tri6
from solidum.materials.elastic_2d import Elastic2D
from solidum.math.assembly import Assembler
from solidum.math.solvers import LinearSolver


# Parámetros físicos del benchmark.
A_INNER, B_INNER = 2.0, 1.0      # semiejes de la elipse interior
A_OUTER, B_OUTER = 3.25, 2.75    # semiejes de la elipse exterior
SIGMA_N = 10.0                    # presión normal exterior (tracción del material)
THICKNESS = 0.1
E_YOUNG = 210_000.0
NU = 0.3

SIGMA_YY_AT_D_REF = 92.7          # σ_yy(D) — valor de referencia NAFEMS


def _elliptic_node(s: float, t: float) -> tuple[float, float]:
    """Mapeo (s, t) → (x, y) en el cuadrante elíptico.

    s ∈ [0, π/2]: ángulo paramétrico (s=0 → eje x, s=π/2 → eje y).
    t ∈ [0, 1]:   coordenada radial paramétrica (0 = elipse interior,
                   1 = elipse exterior).
    """
    a = (1.0 - t) * A_INNER + t * A_OUTER
    b = (1.0 - t) * B_INNER + t * B_OUTER
    return a * np.cos(s), b * np.sin(s)


def _outward_normal_at_outer(s: float) -> np.ndarray:
    """Normal unitaria exterior a la elipse exterior en el punto (s, t=1).

    La gradiente de F(x,y) = x²/a² + y²/b² − 1 es (2x/a², 2y/b²); apunta
    hacia el exterior. Se normaliza.
    """
    x = A_OUTER * np.cos(s)
    y = B_OUTER * np.sin(s)
    nx = x / A_OUTER**2
    ny = y / B_OUTER**2
    norm = np.hypot(nx, ny)
    return np.array([nx / norm, ny / norm])


# =============================================================================
# Malla estructurada elíptica.
# =============================================================================

def _build_elliptic_mesh(elem_cls, nr: int, nt: int, material):
    """Construye una malla estructurada nr×nt en (t, s) para uno de los cinco
    elementos 2D del catálogo.

    Returns
    -------
    domain : Domain
    outer_edges : list[(element, edge_idx)]
        Pares (elemento, edge_idx) cuya arista coincide con la elipse exterior
        (t = 1). Para aplicar la presión normal saliente.
    nid_map : dict[(I, J), Node]
    nI, nJ : int
    point_D : Node
        Nodo en (2, 0) — para evaluar σ_yy(D) por interpolación nodal.
    """
    is_quadratic = elem_cls in (Quad8, Quad9, Tri6)
    sub = 2 if is_quadratic else 1
    nI = nr * sub + 1    # índice radial (t)
    nJ = nt * sub + 1    # índice angular (s)

    domain = Domain()
    nid_map: dict[tuple[int, int], object] = {}
    nid = 0
    for J in range(nJ):
        for I in range(nI):
            is_center = is_quadratic and (I % 2 == 1) and (J % 2 == 1)
            if is_center and elem_cls is Quad8:
                continue
            nid += 1
            t = I / (nI - 1)
            s = (np.pi / 2.0) * J / (nJ - 1)
            x, y = _elliptic_node(s, t)
            nid_map[(I, J)] = domain.add_node(nid, [x, y])

    outer_edges: list[tuple[object, int]] = []
    eid = 0

    if elem_cls is Quad4:
        for jy in range(nt):
            for ix in range(nr):
                n0 = nid_map[(ix,   jy)]
                n1 = nid_map[(ix+1, jy)]
                n2 = nid_map[(ix+1, jy+1)]
                n3 = nid_map[(ix,   jy+1)]
                eid += 1
                elem = Quad4(eid, [n0, n1, n2, n3], material, thickness=THICKNESS)
                domain.add_element(elem)
                if ix == nr - 1:
                    # Arista exterior (t=1) entre n1 y n2 = edge_idx 1.
                    outer_edges.append((elem, 1))

    elif elem_cls in (Quad8, Quad9):
        for jy in range(nt):
            for ix in range(nr):
                I0, J0 = 2 * ix, 2 * jy
                c1 = nid_map[(I0,   J0)]
                c2 = nid_map[(I0+2, J0)]
                c3 = nid_map[(I0+2, J0+2)]
                c4 = nid_map[(I0,   J0+2)]
                m12 = nid_map[(I0+1, J0)]
                m23 = nid_map[(I0+2, J0+1)]
                m34 = nid_map[(I0+1, J0+2)]
                m41 = nid_map[(I0,   J0+1)]
                nodes = [c1, c2, c3, c4, m12, m23, m34, m41]
                if elem_cls is Quad9:
                    nodes.append(nid_map[(I0+1, J0+1)])
                eid += 1
                elem = elem_cls(eid, nodes, material, thickness=THICKNESS)
                domain.add_element(elem)
                if ix == nr - 1:
                    outer_edges.append((elem, 1))

    elif elem_cls is Tri3:
        for jy in range(nt):
            for ix in range(nr):
                c1 = nid_map[(ix,   jy)]
                c2 = nid_map[(ix+1, jy)]
                c3 = nid_map[(ix+1, jy+1)]
                c4 = nid_map[(ix,   jy+1)]
                eid += 1
                tri_a = Tri3(eid, [c1, c2, c3], material, thickness=THICKNESS)
                domain.add_element(tri_a)
                eid += 1
                domain.add_element(Tri3(eid, [c1, c3, c4], material, thickness=THICKNESS))
                if ix == nr - 1:
                    # Tri A (c1, c2, c3): arista entre c2 (vertex 1) y c3 (vertex 2) = edge 1.
                    outer_edges.append((tri_a, 1))

    elif elem_cls is Tri6:
        for jy in range(nt):
            for ix in range(nr):
                I0, J0 = 2 * ix, 2 * jy
                c1 = nid_map[(I0,   J0)]
                c2 = nid_map[(I0+2, J0)]
                c3 = nid_map[(I0+2, J0+2)]
                c4 = nid_map[(I0,   J0+2)]
                m12 = nid_map[(I0+1, J0)]
                m23 = nid_map[(I0+2, J0+1)]
                m34 = nid_map[(I0+1, J0+2)]
                m41 = nid_map[(I0,   J0+1)]
                center = nid_map[(I0+1, J0+1)]
                eid += 1
                tri_a = Tri6(eid, [c1, c2, c3, m12, m23, center],
                             material, thickness=THICKNESS)
                domain.add_element(tri_a)
                eid += 1
                domain.add_element(Tri6(eid, [c1, c3, c4, center, m34, m41],
                                        material, thickness=THICKNESS))
                if ix == nr - 1:
                    outer_edges.append((tri_a, 1))

    else:
        raise ValueError(f"Elemento no soportado: {elem_cls.__name__}")

    # BCs simétricas:
    #   s = 0 (J = 0)   → y = 0       → u_y = 0
    #   s = π/2 (J = nJ-1) → x = 0     → u_x = 0
    for I in range(nI):
        if (I, 0) in nid_map:
            nid_map[(I, 0)].fix_dof('uy', 0.0)
        if (I, nJ - 1) in nid_map:
            nid_map[(I, nJ - 1)].fix_dof('ux', 0.0)

    domain.generate_equation_numbers(verbose=False)
    # El punto D = (2, 0) corresponde a (I=0, J=0) — esquina inferior-izquierda
    # del grid paramétrico (t=0, s=0).
    point_D = nid_map[(0, 0)]
    return domain, outer_edges, nid_map, nI, nJ, point_D


def _apply_outer_pressure(domain: Domain,
                           outer_edges: list[tuple[object, int]]) -> np.ndarray:
    """Tracción normal saliente σ_n·n̂ sobre la elipse exterior.

    Como el borde no es de curvatura constante en s (la normal varía con s),
    se aplica como tracción uniforme **por elemento** evaluada en el centro
    paramétrico (s_mid) del borde. Error O(h²) que se anula con refinamiento.
    """
    F = np.zeros(domain.total_dofs)
    for elem, edge_idx in outer_edges:
        # Extremos de la arista — primer y último nodo de la tupla.
        edge_tuple = elem.EDGE_NODES[edge_idx]
        a, c = edge_tuple[0], edge_tuple[-1]
        xa = np.asarray(elem.nodes[a].coordinates[:2], dtype=np.float64)
        xc = np.asarray(elem.nodes[c].coordinates[:2], dtype=np.float64)
        # s en el centro de la arista: invertir la parametrización
        # (x/A_OUTER)² + (y/B_OUTER)² = 1; usar arctan del punto medio.
        xmid = 0.5 * (xa + xc)
        # arctan2(y/B_OUTER, x/A_OUTER) da el parámetro angular.
        s_mid = np.arctan2(xmid[1] / B_OUTER, xmid[0] / A_OUTER)
        n_hat = _outward_normal_at_outer(s_mid)
        t_vec = SIGMA_N * n_hat
        f_elem = elem.compute_edge_traction(edge_idx, t_vec)
        for k_local, k_global in enumerate(elem.get_global_dof_indices()):
            if k_global >= 0:
                F[k_global] += f_elem[k_local]
    return F


def _sigma_yy_at_nearest_gauss(node, U: np.ndarray, domain: Domain) -> float:
    """σ_yy en el punto de Gauss más próximo al nodo dado.

    Para nodos en una esquina del dominio (como D), este es el "best
    sampling point" disponible sin extrapolación nodal — converge a la
    solución exacta al refinar la malla a tasa h (Quad) o más lenta (Tri,
    cuyo Gauss está en el centroide del elemento).
    """
    target = np.asarray(node.coordinates[:2], dtype=np.float64)
    best_dist = np.inf
    best_sigma_yy = 0.0
    for elem in domain.elements.values():
        if node not in elem.nodes:
            continue
        gs = elem.compute_gauss_state(U)
        for (xy, sig) in zip(gs['points_global'], gs['stress']):
            d = float(np.hypot(xy[0] - target[0], xy[1] - target[1]))
            if d < best_dist:
                best_dist = d
                best_sigma_yy = float(sig[1])  # componente 1 de Voigt 2D
    if best_dist == np.inf:
        raise RuntimeError("Nodo no pertenece a ningún elemento del dominio.")
    return best_sigma_yy


# =============================================================================
# Tests parametrizados por elemento.
# =============================================================================

# (elem_cls, nr, nt, tol_rel_sigma_yy_D) — calibradas empíricamente con
# σ_yy en el Gauss más cercano. Reflejan la limitación del recovery sin SPR.
PARAMS = [
    pytest.param(Quad4, 32, 32, 0.05, id='Quad4_32x32'),
    pytest.param(Quad8,  8,  8, 0.10, id='Quad8_8x8'),
    pytest.param(Quad9,  8,  8, 0.10, id='Quad9_8x8'),
    pytest.param(Tri3,  32, 32, 0.20, id='Tri3_32x32'),
    pytest.param(Tri6,  16, 16, 0.20, id='Tri6_16x16'),
]


@pytest.mark.parametrize("elem_cls, nr, nt, tol_rel", PARAMS)
def test_le1_sigma_yy_at_point_D(elem_cls, nr, nt, tol_rel):
    """σ_yy en el Gauss más cercano a D = (2, 0) coincide con NAFEMS (92.7).

    Tolerancia documenta la limitación del recovery sin SPR; ver módulo
    docstring para discusión.
    """
    material = Elastic2D(E=E_YOUNG, nu=NU, hypothesis='plane_stress')
    domain, outer_edges, _, _, _, point_D = _build_elliptic_mesh(
        elem_cls, nr, nt, material,
    )
    F = _apply_outer_pressure(domain, outer_edges)
    U = LinearSolver(Assembler(domain)).solve(F)
    sigma_yy = _sigma_yy_at_nearest_gauss(point_D, U, domain)
    rel_err = abs(sigma_yy - SIGMA_YY_AT_D_REF) / SIGMA_YY_AT_D_REF
    assert rel_err < tol_rel, (
        f"{elem_cls.__name__} {nr}×{nt}: σ_yy(D, nearest Gauss)={sigma_yy:.3f} "
        f"vs ref={SIGMA_YY_AT_D_REF}, err_rel={rel_err:.4%} > {tol_rel:.2%}"
    )


@pytest.mark.parametrize("elem_cls, meshes", [
    pytest.param(Quad4, (8, 16, 32), id='Quad4'),
    pytest.param(Quad8, (2, 4, 8),   id='Quad8'),
    pytest.param(Quad9, (2, 4, 8),   id='Quad9'),
    pytest.param(Tri3,  (8, 16, 32), id='Tri3'),
    pytest.param(Tri6,  (4, 8, 16),  id='Tri6'),
])
def test_le1_h_refinement_converges(elem_cls, meshes):
    """Convergencia h: σ_yy(Gauss más cercano a D) → 92.7 monotónamente."""
    material = Elastic2D(E=E_YOUNG, nu=NU, hypothesis='plane_stress')
    errors = []
    for n in meshes:
        domain, outer_edges, _, _, _, point_D = _build_elliptic_mesh(
            elem_cls, nr=n, nt=n, material=material,
        )
        F = _apply_outer_pressure(domain, outer_edges)
        U = LinearSolver(Assembler(domain)).solve(F)
        sigma_yy = _sigma_yy_at_nearest_gauss(point_D, U, domain)
        errors.append(abs(sigma_yy - SIGMA_YY_AT_D_REF))
    assert errors[1] < errors[0] and errors[2] < errors[1], (
        f"{elem_cls.__name__} LE1: error no decrece monótonamente "
        f"con mallas {meshes}: errores={errors}"
    )


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
