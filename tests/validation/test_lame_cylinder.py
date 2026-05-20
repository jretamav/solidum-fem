"""Cilindro grueso bajo presión interna — solución de Lamé (plane strain).

Referencia
----------
Timoshenko, S.P.; Goodier, J.N. (1970). *Theory of Elasticity*, 3rd ed.,
McGraw-Hill, §28 ("Pure shear. Stresses in a hollow circular cylinder").

Concepto
--------
Cuadrante de corona circular sometido a presión interna uniforme. La
solución de Lamé en coordenadas polares es:

    σ_rr(r) = A − B/r²
    σ_θθ(r) = A + B/r²
    u_r(r)  = (1+ν)/E · [(1−2ν)·A·r + B/r]      (plane strain)

con A = p·Rᵢ²/(Rₑ²−Rᵢ²) y B = p·Rᵢ²·Rₑ²/(Rₑ²−Rᵢ²).

Es un benchmark canónico de:

1. Anisotropía radial-circunferencial (σ_rr ≠ σ_θθ): caza errores de signo
   o transposición en la matriz B / matriz constitutiva C.
2. Hipótesis plane strain: σ_zz ≠ 0 internamente, aunque no se observa en
   2D; el módulo efectivo cambia (1−ν²/(1−ν))·E respecto a plane stress.
3. Convergencia bajo refinamiento h con malla estructurada en (r, θ).

Geometría y carga
-----------------
- Rᵢ = 1, Rₑ = 2, p = 1, E = 1000, ν = 0.3, espesor = 1.
- BCs simétricas: u_y=0 en y=0; u_x=0 en x=0.
- Presión interna aplicada como tracción radial uniforme sobre cada
  elemento del borde Rᵢ, evaluada en el centro angular del borde
  (aproximación constante por elemento que converge como O(h²)).

Convención de evaluación
------------------------
Los esfuerzos se evalúan en los puntos de Gauss y se transforman a polares:

    σ_rr = σ_xx·c² + σ_yy·s² + 2·σ_xy·c·s
    σ_θθ = σ_xx·s² + σ_yy·c² − 2·σ_xy·c·s

donde (c, s) = (cosθ, sinθ) del punto de Gauss en globales.

Métrica de error y tolerancias
------------------------------
El test usa **error L²-relativo global** sobre los puntos de Gauss del
interior del dominio (descartando los anillos próximos al borde, donde la
aproximación "tracción constante por arista" introduce error de frontera
O(h)):

    e_rel(σ) = sqrt( Σ_g (σ_num − σ_ref)² ) / sqrt( Σ_g σ_ref² )

Es la convención estándar en validación FEM (Bathe §4.3, Hughes §3.10) y
amortigua el ruido por punto/anillo. Las tolerancias por elemento son
fruto de medida empírica con un margen, no de adivinación a priori — quedan
documentadas como límites de aceptación.

Adicionalmente, un test de convergencia h con Quad4 verifica que el error
decrece monótonamente al refinar, lo cual es la propiedad cualitativa
más fuerte que se pide a una formulación.
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
R_INNER = 1.0
R_OUTER = 2.0
P_INNER = 1.0
E_YOUNG = 1000.0
NU = 0.3


def _lame_stress_polar(r: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """σ_rr(r), σ_θθ(r) según Lamé."""
    A = P_INNER * R_INNER**2 / (R_OUTER**2 - R_INNER**2)
    B = P_INNER * R_INNER**2 * R_OUTER**2 / (R_OUTER**2 - R_INNER**2)
    return A - B / r**2, A + B / r**2


def _lame_radial_displacement(r: np.ndarray) -> np.ndarray:
    """u_r(r) en plane strain."""
    A = P_INNER * R_INNER**2 / (R_OUTER**2 - R_INNER**2)
    B = P_INNER * R_INNER**2 * R_OUTER**2 / (R_OUTER**2 - R_INNER**2)
    return (1.0 + NU) / E_YOUNG * ((1.0 - 2.0 * NU) * A * r + B / r)


def _cart_to_polar_stress(sigma_xx: np.ndarray, sigma_yy: np.ndarray,
                          sigma_xy: np.ndarray,
                          x: np.ndarray, y: np.ndarray
                          ) -> tuple[np.ndarray, np.ndarray]:
    """Transforma σ globales a (σ_rr, σ_θθ) usando el ángulo del punto."""
    r = np.hypot(x, y)
    c = x / r
    s = y / r
    sigma_rr = sigma_xx * c**2 + sigma_yy * s**2 + 2.0 * sigma_xy * c * s
    sigma_tt = sigma_xx * s**2 + sigma_yy * c**2 - 2.0 * sigma_xy * c * s
    return sigma_rr, sigma_tt


# =============================================================================
# Generador de malla polar estructurada.
# =============================================================================

def _polar_node(I: int, J: int, nI: int, nJ: int) -> tuple[float, float]:
    r = R_INNER + (R_OUTER - R_INNER) * I / (nI - 1)
    theta = (np.pi / 2.0) * J / (nJ - 1)
    return r * np.cos(theta), r * np.sin(theta)


def _build_polar_mesh(elem_cls, nr: int, nt: int, material):
    """Construye un dominio con malla estructurada en (r, θ) para uno de los
    cinco elementos 2D del catálogo.

    Returns
    -------
    domain : Domain
    inner_edges : list[(element, edge_idx)]
        Pares de (elemento, índice de arista) cuya arista coincide con r=Rᵢ.
        Para aplicar la presión interna como tracción radial.
    """
    # Quad4/Tri3 → 1 subdivisión por celda; Quad8/Quad9/Tri6 → 2.
    is_quadratic = elem_cls in (Quad8, Quad9, Tri6)
    sub = 2 if is_quadratic else 1
    nI = nr * sub + 1
    nJ = nt * sub + 1

    domain = Domain()
    nid_map: dict[tuple[int, int], object] = {}
    nid = 0
    for J in range(nJ):
        for I in range(nI):
            # Quad8 omite el nodo central de la celda (I,J ambos impares).
            is_center = is_quadratic and (I % 2 == 1) and (J % 2 == 1)
            if is_center and elem_cls is Quad8:
                continue
            nid += 1
            x, y = _polar_node(I, J, nI, nJ)
            nid_map[(I, J)] = domain.add_node(nid, [x, y])

    inner_edges: list[tuple[object, int]] = []
    eid = 0

    if elem_cls is Quad4:
        for jy in range(nt):
            for ix in range(nr):
                n0 = nid_map[(ix,   jy)]      # (r_i, θ_j)
                n1 = nid_map[(ix+1, jy)]      # (r_{i+1}, θ_j)
                n2 = nid_map[(ix+1, jy+1)]    # (r_{i+1}, θ_{j+1})
                n3 = nid_map[(ix,   jy+1)]    # (r_i, θ_{j+1})
                eid += 1
                elem = Quad4(eid, [n0, n1, n2, n3], material, thickness=1.0)
                domain.add_element(elem)
                if ix == 0:
                    # Arista interior (r=Rᵢ) entre n3 y n0 = edge_idx 3.
                    inner_edges.append((elem, 3))

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
                elem = elem_cls(eid, nodes, material, thickness=1.0)
                domain.add_element(elem)
                if ix == 0:
                    inner_edges.append((elem, 3))

    elif elem_cls is Tri3:
        # Cada celda → dos Tri3 a lo largo de la diagonal c1→c3.
        for jy in range(nt):
            for ix in range(nr):
                c1 = nid_map[(ix,   jy)]
                c2 = nid_map[(ix+1, jy)]
                c3 = nid_map[(ix+1, jy+1)]
                c4 = nid_map[(ix,   jy+1)]
                eid += 1
                # Tri A: (c1, c2, c3) — no toca r=Rᵢ.
                domain.add_element(Tri3(eid, [c1, c2, c3], material, thickness=1.0))
                eid += 1
                # Tri B: (c1, c3, c4) — arista entre c4 (vertex 2) y c1 (vertex 0) = edge 2.
                tri_b = Tri3(eid, [c1, c3, c4], material, thickness=1.0)
                domain.add_element(tri_b)
                if ix == 0:
                    inner_edges.append((tri_b, 2))

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
                domain.add_element(Tri6(eid, [c1, c2, c3, m12, m23, center],
                                        material, thickness=1.0))
                eid += 1
                tri_b = Tri6(eid, [c1, c3, c4, center, m34, m41],
                             material, thickness=1.0)
                domain.add_element(tri_b)
                if ix == 0:
                    inner_edges.append((tri_b, 2))

    else:  # pragma: no cover — guard contra elementos no soportados.
        raise ValueError(f"Elemento no soportado: {elem_cls.__name__}")

    # BCs de simetría: u_y=0 en y=0 (J=0); u_x=0 en x=0 (J=nJ-1).
    for I in range(nI):
        if (I, 0) in nid_map:
            nid_map[(I, 0)].fix_dof('uy', 0.0)
        if (I, nJ - 1) in nid_map:
            nid_map[(I, nJ - 1)].fix_dof('ux', 0.0)

    domain.generate_equation_numbers(verbose=False)
    return domain, inner_edges, nid_map, nI, nJ


def _edge_endpoints(elem, edge_idx: int) -> tuple[int, int]:
    """Extremos de la arista ``edge_idx`` — primer y último nodo de la tupla.

    ``EDGE_NODES`` para Quad4/Tri3 es ``(n_a, n_c)`` y para Quad8/Quad9/Tri6
    ``(n_a, n_mid, n_c)``. En ambos casos ``[0]`` y ``[-1]`` son los extremos.
    """
    edge_tuple = elem.EDGE_NODES[edge_idx]
    return edge_tuple[0], edge_tuple[-1]


def _apply_internal_pressure(domain: Domain,
                              inner_edges: list[tuple[object, int]]) -> np.ndarray:
    """Tracción radial p·r̂ uniforme por elemento, evaluada en el centro del borde."""
    F = np.zeros(domain.total_dofs)
    for elem, edge_idx in inner_edges:
        a, c = _edge_endpoints(elem, edge_idx)
        xa = np.asarray(elem.nodes[a].coordinates[:2], dtype=np.float64)
        xc = np.asarray(elem.nodes[c].coordinates[:2], dtype=np.float64)
        xmid = 0.5 * (xa + xc)
        r_mid = np.hypot(xmid[0], xmid[1])
        n_radial = xmid / r_mid                          # r̂ en el centro
        t_vec = P_INNER * n_radial                       # tracción aplicada al material
        f_elem = elem.compute_edge_traction(edge_idx, t_vec)
        for k_local, k_global in enumerate(elem.get_global_dof_indices()):
            if k_global >= 0:
                F[k_global] += f_elem[k_local]
    return F


def _solve_and_collect_gauss_stress(domain: Domain) -> tuple[np.ndarray, np.ndarray,
                                                              np.ndarray, np.ndarray,
                                                              np.ndarray]:
    """Resuelve el problema y devuelve (x_g, y_g, σ_rr, σ_θθ, r_g) en todos los Gauss."""
    assembler = Assembler(domain)
    material = next(iter(domain.elements.values())).material
    F = _apply_internal_pressure(domain,
                                  [(e, _inner_edge(e)) for e in domain.elements.values()
                                   if _inner_edge(e) is not None])
    U = LinearSolver(assembler).solve(F)

    xs, ys, srrs, stts = [], [], [], []
    for elem in domain.elements.values():
        gs = elem.compute_gauss_state(U)
        xg = gs['points_global'][:, 0]
        yg = gs['points_global'][:, 1]
        sxx = gs['stress'][:, 0]
        syy = gs['stress'][:, 1]
        sxy = gs['stress'][:, 2]
        srr, stt = _cart_to_polar_stress(sxx, syy, sxy, xg, yg)
        xs.append(xg); ys.append(yg); srrs.append(srr); stts.append(stt)
    xg = np.concatenate(xs); yg = np.concatenate(ys)
    srr_g = np.concatenate(srrs); stt_g = np.concatenate(stts)
    rg = np.hypot(xg, yg)
    return xg, yg, srr_g, stt_g, rg, U


def _inner_edge(elem) -> int | None:
    """Devuelve el edge_idx del borde r=Rᵢ si el elemento lo tiene, o None.

    Detecta el borde comparando radios de los nodos del elemento: la arista
    interna tiene ambos nodos extremos en r ≈ Rᵢ.
    """
    for edge_idx in range(len(elem.EDGE_NODES)):
        a, c = _edge_endpoints(elem, edge_idx)
        xa = np.asarray(elem.nodes[a].coordinates[:2], dtype=np.float64)
        xc = np.asarray(elem.nodes[c].coordinates[:2], dtype=np.float64)
        if (abs(np.hypot(*xa) - R_INNER) < 1e-9
                and abs(np.hypot(*xc) - R_INNER) < 1e-9):
            return edge_idx
    return None


# =============================================================================
# Tests parametrizados por elemento — métrica L²-relativa global.
# =============================================================================

# Error L²-relativo global del campo de esfuerzos sobre Gauss interiores.
# Tolerancias calibradas empíricamente con margen razonable; reflejan el
# orden de convergencia teórico de cada formulación:
#   - Quad4/Tri3: σ ∈ P0 por elemento → O(h¹) en σ.
#   - Quad8/Quad9/Tri6: σ ∈ P1 por elemento → O(h²) en σ.
# Para Quad4 8×8 el error σ_rr es ≈18%; con malla 16×16 baja a ≈9%, lo cual
# valida la tasa de convergencia y la formulación, pero impone un techo en
# la malla 8×8. El test de h-refinement complementa la verificación
# cualitativa (la propiedad arquitecturalmente más fuerte).
PARAMS = [
    pytest.param(Quad4, 8, 8, 0.22, 0.05,  id='Quad4_8x8'),
    pytest.param(Tri3,  8, 8, 0.22, 0.11,  id='Tri3_8x8'),
    pytest.param(Quad8, 4, 4, 0.06, 0.02,  id='Quad8_4x4'),
    pytest.param(Quad9, 4, 4, 0.06, 0.02,  id='Quad9_4x4'),
    pytest.param(Tri6,  4, 4, 0.06, 0.02,  id='Tri6_4x4'),
]


def _l2_relative_error(num: np.ndarray, ref: np.ndarray) -> float:
    """Error L² relativo discreto entre dos campos en los mismos puntos."""
    return float(np.sqrt(np.mean((num - ref) ** 2))
                 / max(np.sqrt(np.mean(ref ** 2)), 1e-12))


@pytest.mark.parametrize("elem_cls, nr, nt, tol_rr, tol_tt", PARAMS)
def test_lame_stress_l2_error(elem_cls, nr, nt, tol_rr, tol_tt):
    """Error L²-relativo de σ_rr y σ_θθ en Gauss interiores ≤ tolerancia."""
    material = Elastic2D(E=E_YOUNG, nu=NU, hypothesis='plane_strain')
    domain, _, _, _, _ = _build_polar_mesh(elem_cls, nr, nt, material)
    xg, yg, srr_g, stt_g, rg, _ = _solve_and_collect_gauss_stress(domain)

    # Descartar Gauss próximos a los bordes Rᵢ (BC tracción aproximada) y Rₑ
    # (σ_rr → 0, denominador degenerado en error relativo).
    mask = (rg > R_INNER + 0.15) & (rg < R_OUTER - 0.15)
    assert mask.sum() > 0, f"Sin Gauss interiores en malla {elem_cls.__name__}"

    srr_ref, stt_ref = _lame_stress_polar(rg[mask])
    e_rr = _l2_relative_error(srr_g[mask], srr_ref)
    e_tt = _l2_relative_error(stt_g[mask], stt_ref)

    assert e_rr < tol_rr, (
        f"{elem_cls.__name__} {nr}×{nt}: e_L2(σ_rr)={e_rr:.4%} > {tol_rr:.2%}"
    )
    assert e_tt < tol_tt, (
        f"{elem_cls.__name__} {nr}×{nt}: e_L2(σ_θθ)={e_tt:.4%} > {tol_tt:.2%}"
    )


@pytest.mark.parametrize("elem_cls, nr, nt, tol_rr, tol_tt", PARAMS)
def test_lame_radial_displacement(elem_cls, nr, nt, tol_rr, tol_tt):
    """u_r en los nodos del borde interno y externo coincide con Lamé.

    El desplazamiento converge **más rápido** que el esfuerzo (σ pierde un
    orden por la derivación de B), por eso usamos una tolerancia más estricta:
    igual a la del hoop stress, que también es un campo "regular".
    """
    material = Elastic2D(E=E_YOUNG, nu=NU, hypothesis='plane_strain')
    domain, _, nid_map, nI, nJ = _build_polar_mesh(elem_cls, nr, nt, material)
    _, _, _, _, _, U = _solve_and_collect_gauss_stress(domain)

    tol_u = tol_tt
    for I in (0, nI - 1):  # bordes interno y externo
        r_ref = R_INNER if I == 0 else R_OUTER
        ur_ref = _lame_radial_displacement(np.array([r_ref]))[0]
        for J in range(nJ):
            if (I, J) not in nid_map:
                continue
            node = nid_map[(I, J)]
            ux = U[node.dofs['ux']]
            uy = U[node.dofs['uy']]
            x, y = node.coordinates[:2]
            r = np.hypot(x, y)
            ur = (ux * x + uy * y) / r
            rel_err = abs(ur - ur_ref) / max(abs(ur_ref), 1e-12)
            assert rel_err < tol_u, (
                f"{elem_cls.__name__}: u_r(r={r:.4f}, "
                f"θ={np.degrees(np.arctan2(y, x)):.1f}°): "
                f"num={ur:.6e}, ref={ur_ref:.6e}, err={rel_err:.4%} > {tol_u:.2%}"
            )


def test_lame_h_refinement_quad4_converges():
    """Refinamiento h con Quad4: error L²-relativo de σ_rr y σ_θθ decrece.

    Convergencia esperada O(h¹) para Quad4 (σ ∈ P0 por elemento).
    """
    material = Elastic2D(E=E_YOUNG, nu=NU, hypothesis='plane_strain')
    errors_rr = []
    errors_tt = []
    for n in (4, 8, 16):
        domain, _, _, _, _ = _build_polar_mesh(Quad4, nr=n, nt=n, material=material)
        xg, yg, srr_g, stt_g, rg, _ = _solve_and_collect_gauss_stress(domain)
        mask = (rg > R_INNER + 0.15) & (rg < R_OUTER - 0.15)
        srr_ref, stt_ref = _lame_stress_polar(rg[mask])
        errors_rr.append(_l2_relative_error(srr_g[mask], srr_ref))
        errors_tt.append(_l2_relative_error(stt_g[mask], stt_ref))
    assert errors_rr[1] < errors_rr[0] and errors_rr[2] < errors_rr[1], (
        f"σ_rr no decrece monótonamente: {errors_rr}"
    )
    assert errors_tt[1] < errors_tt[0] and errors_tt[2] < errors_tt[1], (
        f"σ_θθ no decrece monótonamente: {errors_tt}"
    )


def test_lame_h_refinement_quad8_quadratic_order():
    """Quad8 4×4 → 8×8 reduce el error σ_rr más rápido que orden 1.

    Verifica convergencia cuadrática: pasar de h a h/2 debe reducir el error
    en factor > 2 (idealmente ~4). Toleramos > 2.5 para márgenes.
    """
    material = Elastic2D(E=E_YOUNG, nu=NU, hypothesis='plane_strain')
    errs = []
    for n in (4, 8):
        domain, _, _, _, _ = _build_polar_mesh(Quad8, nr=n, nt=n, material=material)
        xg, yg, srr_g, _, rg, _ = _solve_and_collect_gauss_stress(domain)
        mask = (rg > R_INNER + 0.15) & (rg < R_OUTER - 0.15)
        srr_ref, _ = _lame_stress_polar(rg[mask])
        errs.append(_l2_relative_error(srr_g[mask], srr_ref))
    ratio = errs[0] / errs[1]
    assert ratio > 2.5, (
        f"Quad8: ratio de mejora h→h/2 = {ratio:.2f} < 2.5 (esperado ~4 para "
        f"orden cuadrático). Errores: {errs}"
    )


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
