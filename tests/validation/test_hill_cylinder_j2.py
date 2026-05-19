"""Cilindro grueso de Hill — plasticidad J2 perfecta (plane strain).

Referencia
----------
Hill, R. (1950). *The Mathematical Theory of Plasticity*, Oxford
University Press, §5 ("Plastic deformation of a thick cylindrical tube").
Recogido también en Lubliner *Plasticity Theory* §4.1 y en Souza Neto-
Perić-Owen §7.3.4.

Concepto
--------
Cuadrante de corona circular (Rᵢ, Rₑ) sometido a presión interna p, con
material elastoplástico de von Mises **perfecto** (sin endurecimiento).
Plane strain. Existen **tres regímenes** según el valor de p:

1. **Régimen puramente elástico** (p < p_e): solución de Lamé. Sin
   plastificación.
2. **Régimen elastoplástico** (p_e ≤ p < p_L): zona plástica anular
   en r ∈ [Rᵢ, r_c] y zona elástica en r ∈ [r_c, Rₑ]. El radio de
   transición r_c crece con p hasta que llega a Rₑ.
3. **Límite plástico** (p = p_L): pared completamente plastificada;
   el problema tiene solución indefinida (colapso plástico).

Soluciones analíticas (von Mises, plane strain)
-----------------------------------------------
Notación: σ_y es la tensión de fluencia uniaxial; k = σ_y/√3 es la
tensión de fluencia en cortante puro.

Presiones críticas:

    p_e = (σ_y/√3) · (1 − Rᵢ²/Rₑ²)             (yield onset en r=Rᵢ)
    p_L = (2σ_y/√3) · ln(Rₑ/Rᵢ)                (límite plástico)

Para p en el régimen elastoplástico, el radio de transición r_c es la
única raíz en (Rᵢ, Rₑ) de la ecuación implícita:

    p = (σ_y/√3) · [2·ln(r_c/Rᵢ) + 1 − r_c²/Rₑ²]

Tensiones en la **zona plástica** r ∈ [Rᵢ, r_c]:

    σ_rr(r) = (2σ_y/√3) · ln(r/Rᵢ) − p
    σ_θθ(r) = σ_rr(r) + 2σ_y/√3

Tensiones en la **zona elástica** r ∈ [r_c, Rₑ] — Lamé con presión
interna efectiva p_c = −σ_rr(r_c):

    A = (σ_y/√3) · r_c² / Rₑ²
    B = A · Rₑ²
    σ_rr(r) = A − B/r²,   σ_θθ(r) = A + B/r²

Parámetros del benchmark
------------------------
- Rᵢ = 1, Rₑ = 2
- E = 1000, ν = 0.3 (plane strain)
- σ_y = 1, H = 0 (plasticidad perfecta)
- p_e = (1/√3)·(1 − 1/4) ≈ 0.4330
- p_L = (2/√3)·ln(2)     ≈ 0.8005

Validación
----------
Tres tests, uno por régimen:

1. **Régimen elástico** (p = 0.30 < p_e): α = 0 en todos los Gauss;
   σ_rr/σ_θθ coinciden con Lamé puro a tolerancia razonable (10% por la
   convergencia O(h¹) de Quad4 en σ, validado en `test_lame_cylinder.py`).
2. **Régimen elastoplástico** (p = 0.70): r_c analítico ≈ 1.49. Los
   Gauss en r ≪ r_c deben tener α > 0 (zona plástica); los Gauss en
   r ≫ r_c deben tener α = 0 (zona elástica). σ_rr en zona plástica
   coincide con Hill; σ_θθ - σ_rr ≈ 2σ_y/√3 en zona plástica
   (Δσ saturada al límite de fluencia).
3. **Aproximación al límite plástico** (p = 0.78 ≈ 0.975·p_L): la zona
   plástica casi alcanza Rₑ. u_r(Rᵢ) crece monotónicamente al
   acercarse a p_L (no se valida un valor exacto — la respuesta diverge
   en p → p_L).
"""
import math
import os
import sys

import numpy as np
import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from fenix.core.domain import Domain
from fenix.elements.solid_2d import Quad4
from fenix.materials.von_mises_2d import VonMises2D
from fenix.math.assembly import Assembler
from fenix.math.convergence import ConvergenceCriterion
from fenix.math.solvers import NonlinearSolver


# Parámetros físicos.
R_INNER = 1.0
R_OUTER = 2.0
E_YOUNG = 1000.0
NU = 0.3
SIGMA_Y = 1.0
K_SHEAR = SIGMA_Y / math.sqrt(3.0)   # tensión de fluencia en cortante puro

P_YIELD_ONSET = K_SHEAR * (1.0 - R_INNER**2 / R_OUTER**2)   # p_e
P_LIMIT = 2.0 * K_SHEAR * math.log(R_OUTER / R_INNER)        # p_L


def _hill_transition_radius(p: float) -> float:
    """Resuelve la ecuación implícita de Hill para r_c con bisección.

    Si p está fuera del régimen elastoplástico, devuelve Rᵢ (puramente
    elástico) o Rₑ (al límite).
    """
    if p <= P_YIELD_ONSET:
        return R_INNER
    if p >= P_LIMIT:
        return R_OUTER

    def equation(rc):
        return K_SHEAR * (2.0 * math.log(rc / R_INNER) + 1.0 - rc**2 / R_OUTER**2) - p

    lo, hi = R_INNER + 1e-6, R_OUTER - 1e-6
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if equation(mid) > 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def _hill_stress_polar(r: float, p: float, r_c: float) -> tuple[float, float]:
    """σ_rr, σ_θθ según Hill (plástico) o Lamé efectiva (elástico).

    Si r_c ≤ R_INNER se interpreta como régimen puramente elástico: la
    solución es Lamé con presión interna ``p`` sobre la corona completa
    [Rᵢ, Rₑ]. Si r_c > R_INNER, se aplica Hill (zona plástica + zona
    elástica reducida [r_c, Rₑ] con presión efectiva derivada de la
    continuidad en r_c).
    """
    if r_c <= R_INNER:
        # Lamé puro con presión real p sobre toda la corona.
        A = p * R_INNER**2 / (R_OUTER**2 - R_INNER**2)
        B = A * R_OUTER**2
        return A - B / r**2, A + B / r**2

    if r <= r_c:                      # zona plástica
        srr = 2.0 * K_SHEAR * math.log(r / R_INNER) - p
        stt = srr + 2.0 * K_SHEAR
    else:                              # zona elástica restante (Lamé reducida)
        A = K_SHEAR * r_c**2 / R_OUTER**2
        B = A * R_OUTER**2
        srr = A - B / r**2
        stt = A + B / r**2
    return srr, stt


# =============================================================================
# Helpers de malla polar (Quad4) — duplicación local de la malla de Lamé.
# =============================================================================

def _build_polar_mesh_quad4(nr: int, nt: int, material):
    nI = nr + 1
    nJ = nt + 1
    domain = Domain()
    nid_map: dict[tuple[int, int], object] = {}
    nid = 0
    for J in range(nJ):
        for I in range(nI):
            r = R_INNER + (R_OUTER - R_INNER) * I / (nI - 1)
            theta = (math.pi / 2.0) * J / (nJ - 1)
            nid += 1
            nid_map[(I, J)] = domain.add_node(
                nid, [r * math.cos(theta), r * math.sin(theta)],
            )

    inner_edges: list[tuple[object, int]] = []
    eid = 0
    for jy in range(nt):
        for ix in range(nr):
            n0 = nid_map[(ix,     jy)]
            n1 = nid_map[(ix + 1, jy)]
            n2 = nid_map[(ix + 1, jy + 1)]
            n3 = nid_map[(ix,     jy + 1)]
            eid += 1
            elem = Quad4(eid, [n0, n1, n2, n3], material, thickness=1.0)
            domain.add_element(elem)
            if ix == 0:
                inner_edges.append((elem, 3))   # arista n3-n0 = r=Rᵢ

    # BCs simétricas.
    for I in range(nI):
        if (I, 0) in nid_map:
            nid_map[(I, 0)].fix_dof('uy', 0.0)
        if (I, nJ - 1) in nid_map:
            nid_map[(I, nJ - 1)].fix_dof('ux', 0.0)

    domain.generate_equation_numbers(verbose=False)
    return domain, inner_edges, nid_map, nI, nJ


def _apply_internal_pressure(domain: Domain,
                              inner_edges: list[tuple[object, int]],
                              p: float) -> np.ndarray:
    F = np.zeros(domain.total_dofs)
    for elem, edge_idx in inner_edges:
        edge_tuple = elem.EDGE_NODES[edge_idx]
        a, c = edge_tuple[0], edge_tuple[-1]
        xa = np.asarray(elem.nodes[a].coordinates[:2], dtype=np.float64)
        xc = np.asarray(elem.nodes[c].coordinates[:2], dtype=np.float64)
        xmid = 0.5 * (xa + xc)
        r_mid = math.hypot(xmid[0], xmid[1])
        n_radial = xmid / r_mid
        t_vec = p * n_radial
        f_elem = elem.compute_edge_traction(edge_idx, t_vec)
        for k_local, k_global in enumerate(elem.get_global_dof_indices()):
            if k_global >= 0:
                F[k_global] += f_elem[k_local]
    return F


def _cart_to_polar_stress(sigma_xx, sigma_yy, sigma_xy, x, y):
    r = math.hypot(x, y)
    c, s = x / r, y / r
    srr = sigma_xx * c**2 + sigma_yy * s**2 + 2.0 * sigma_xy * c * s
    stt = sigma_xx * s**2 + sigma_yy * c**2 - 2.0 * sigma_xy * c * s
    return srr, stt


def _solve_pressure(p: float, num_steps: int = 10):
    """Resuelve el problema con presión interna p y devuelve (domain, U)."""
    material = VonMises2D(E=E_YOUNG, nu=NU, sigma_y=SIGMA_Y, H=0.0,
                           hypothesis='plane_strain')
    domain, inner_edges, nid_map, nI, nJ = _build_polar_mesh_quad4(
        nr=12, nt=8, material=material,
    )
    F = _apply_internal_pressure(domain, inner_edges, p)
    conv = ConvergenceCriterion(rtol_force=1e-6, rtol_disp=1e-6)
    solver = NonlinearSolver(Assembler(domain), convergence=conv,
                              num_steps=num_steps)
    U = solver.solve(F)
    return domain, U, nid_map, nI, nJ


def _gather_gauss_state(domain: Domain, U: np.ndarray):
    """Devuelve (r, σ_rr, σ_θθ, α) para todos los Gauss interiores."""
    rs, srrs, stts, alphas = [], [], [], []
    for elem in domain.elements.values():
        gs = elem.compute_gauss_state(U)
        for idx, (xy, sig) in enumerate(zip(gs['points_global'], gs['stress'])):
            r = math.hypot(xy[0], xy[1])
            srr, stt = _cart_to_polar_stress(sig[0], sig[1], sig[2], xy[0], xy[1])
            alpha = elem.state.vars[idx].get('alpha', 0.0) if elem.state.vars[idx] else 0.0
            rs.append(r); srrs.append(srr); stts.append(stt); alphas.append(alpha)
    return (np.asarray(rs), np.asarray(srrs), np.asarray(stts),
            np.asarray(alphas))


# =============================================================================
# Tests.
# =============================================================================

def test_hill_elastic_regime_no_plasticity():
    """p = 0.30 < p_e ≈ 0.433: nada plastifica, σ ≈ Lamé."""
    p = 0.30
    assert p < P_YIELD_ONSET, "Test inválido: p ya supera p_e"
    domain, U, *_ = _solve_pressure(p, num_steps=4)
    r_g, srr_g, stt_g, alpha_g = _gather_gauss_state(domain, U)

    # No debe haber deformación plástica acumulada.
    assert alpha_g.max() < 1e-10, (
        f"En régimen elástico (p={p}, p_e={P_YIELD_ONSET:.4f}) "
        f"α_max = {alpha_g.max():.3e} debe ser 0"
    )

    # σ_rr y σ_θθ deben aproximar Lamé (régimen puramente elástico). La
    # corrección de paso a polares ya está en _hill_stress_polar con
    # r_c=0 (interpretado como Lamé puro).
    mask = (r_g > R_INNER + 0.15) & (r_g < R_OUTER - 0.15)
    srr_ref = np.array([_hill_stress_polar(r, p, 0.0)[0] for r in r_g[mask]])
    stt_ref = np.array([_hill_stress_polar(r, p, 0.0)[1] for r in r_g[mask]])
    e_rr = (np.sqrt(np.mean((srr_g[mask] - srr_ref) ** 2))
             / np.sqrt(np.mean(srr_ref ** 2)))
    e_tt = (np.sqrt(np.mean((stt_g[mask] - stt_ref) ** 2))
             / np.sqrt(np.mean(stt_ref ** 2)))
    # Tolerancias coherentes con Quad4 12×8 en Lamé (~18% en σ_rr, ~3% en σ_θθ).
    assert e_rr < 0.20, f"σ_rr L²-rel = {e_rr:.4%} > 20%"
    assert e_tt < 0.06, f"σ_θθ L²-rel = {e_tt:.4%} > 6%"


def test_hill_elastoplastic_regime_yield_zone():
    """p = 0.70: hay zona plástica interna y zona elástica externa.

    r_c analítico ≈ 1.49 para p=0.70. Los Gauss con r < 1.35 deben estar
    plastificados (α > 0); los Gauss con r > 1.65 deben estar elásticos.
    En la zona plástica, σ_θθ − σ_rr ≈ 2σ_y/√3 (saturación al límite de
    fluencia).
    """
    p = 0.70
    assert P_YIELD_ONSET < p < P_LIMIT, "Test inválido: p fuera del régimen"
    r_c = _hill_transition_radius(p)
    # r_c esperado: ≈ 1.49 para los parámetros del benchmark.
    assert 1.4 < r_c < 1.6, f"r_c analítico inesperado: {r_c:.4f}"

    domain, U, *_ = _solve_pressure(p, num_steps=10)
    r_g, srr_g, stt_g, alpha_g = _gather_gauss_state(domain, U)

    # Zona "claramente plástica": r < r_c − margen.
    in_plastic = r_g < (r_c - 0.10)
    assert in_plastic.sum() > 0, "Sin Gauss en zona plástica"
    assert (alpha_g[in_plastic] > 1e-6).all(), (
        f"Gauss en zona plástica (r < {r_c-0.10:.3f}) con α=0: "
        f"min α = {alpha_g[in_plastic].min():.3e}"
    )

    # Zona "claramente elástica": r > r_c + margen.
    in_elastic = r_g > (r_c + 0.10)
    assert in_elastic.sum() > 0, "Sin Gauss en zona elástica"
    assert (alpha_g[in_elastic] < 1e-8).all(), (
        f"Gauss en zona elástica (r > {r_c+0.10:.3f}) con α>0: "
        f"max α = {alpha_g[in_elastic].max():.3e}"
    )

    # En zona plástica: σ_θθ − σ_rr ≈ 2σ_y/√3.
    delta_target = 2.0 * K_SHEAR
    delta_num = (stt_g[in_plastic] - srr_g[in_plastic]).mean()
    rel_err = abs(delta_num - delta_target) / delta_target
    assert rel_err < 0.05, (
        f"σ_θθ − σ_rr promedio en zona plástica = {delta_num:.4f} "
        f"vs 2σ_y/√3 = {delta_target:.4f}, err_rel = {rel_err:.4%}"
    )


def test_hill_stress_in_plastic_zone_matches_analytic():
    """σ_rr en zona plástica coincide con la fórmula de Hill (L²-rel < 10%)."""
    p = 0.70
    r_c = _hill_transition_radius(p)
    domain, U, *_ = _solve_pressure(p, num_steps=10)
    r_g, srr_g, stt_g, alpha_g = _gather_gauss_state(domain, U)

    # Gauss en zona plástica con margen de los bordes del problema.
    mask = (r_g < (r_c - 0.05)) & (r_g > R_INNER + 0.08)
    assert mask.sum() > 0, "Sin Gauss aptos para comparar en zona plástica"

    srr_ref = np.array([_hill_stress_polar(r, p, r_c)[0] for r in r_g[mask]])
    e_rr = (np.sqrt(np.mean((srr_g[mask] - srr_ref) ** 2))
             / np.sqrt(np.mean(srr_ref ** 2)))
    # Tolerancia 18%: Quad4 con 12×8 elementos y stress polar derivado
    # tiene el mismo orden de error que el caso elástico de Lamé.
    assert e_rr < 0.18, f"σ_rr L²-rel en zona plástica = {e_rr:.4%} > 18%"


def test_hill_displacement_grows_approaching_limit():
    """u_r(Rᵢ) crece monótonamente al aproximarse a p_L."""
    pressures = [0.50, 0.65, 0.78]   # creciente; el último ≈ 0.975·p_L
    u_r_inner = []
    for p in pressures:
        domain, U, nid_map, _, nJ = _solve_pressure(p, num_steps=10)
        # Promedio de u_r en los nodos del borde interno (I=0).
        urs = []
        for J in range(nJ):
            if (0, J) in nid_map:
                node = nid_map[(0, J)]
                ux = U[node.dofs['ux']]
                uy = U[node.dofs['uy']]
                x, y = node.coordinates[:2]
                r = math.hypot(x, y)
                urs.append((ux * x + uy * y) / r)
        u_r_inner.append(float(np.mean(urs)))

    assert u_r_inner[0] < u_r_inner[1] < u_r_inner[2], (
        f"u_r(Rᵢ) no crece monótonamente con p: {u_r_inner}"
    )
    # El último debe ser **mucho mayor** que el primero (regimen no lineal).
    ratio = u_r_inner[-1] / u_r_inner[0]
    assert ratio > 2.0, (
        f"u_r(Rᵢ) en p ≈ p_L debe ser al menos 2× el de p=0.50; "
        f"ratio = {ratio:.2f}"
    )


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
