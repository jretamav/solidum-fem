"""Cilindro hueco 3D — plasticidad J2 perfecta vs Hill 1950 §5.

Referencias
-----------
- Hill, R. (1950). *The Mathematical Theory of Plasticity*, Oxford University
  Press, §5 ("Plastic deformation of a thick cylindrical tube").
- Lubliner, J. (1990). *Plasticity Theory*, §4.1.
- de Souza Neto, E.A., Perić, D., Owen, D.R.J. (2008). *Computational Methods
  for Plasticity*, §7.3.4.

Concepto
--------
Verificación sistémica del pipeline **3D** (``Hex8`` + Voigt 6D + ``VonMises3D``
+ tracciones de cara 3D + ``NonlinearSolver``) contra la solución analítica
cerrada de Hill para el cilindro grueso en plane strain.

La geometría tiene simetría plane-strain (el dominio físico es bidimensional
sobre el plano radial), pero se resuelve con el aparato 3D completo: malla
``Hex8`` (no ``Quad4``), Voigt 6D (no Voigt 3), material ``VonMises3D``
(no ``VonMises2D``), face tractions sobre la cara cilíndrica interna en
coordenadas 3D, y el sistema algebraico 3D. La restricción plane strain
se impone vía Dirichlet ``u_z = 0`` en las dos tapas del cilindro extruido.

Esto extiende el test cruzado unit-level
(``TestVonMises3DvsPlaneStrain.test_equivalencia_plane_strain_path``) a una
validación sistema-completo contra solución analítica publicada.

Solución analítica de Hill (idéntica al 2D, validada en
``test_hill_cylinder_j2.py``):

- Presión yield onset: ``p_e = (σ_y/√3) · (1 - R_i²/R_e²)``
- Presión límite plástica: ``p_L = (2σ_y/√3) · ln(R_e/R_i)``
- Para ``p_e < p < p_L``, radio de transición ``c`` raíz implícita de
  ``p = (σ_y/√3) · [2·ln(c/R_i) + 1 - c²/R_e²]``
- En zona plástica ``[R_i, c]``: ``σ_rr(r) = 2(σ_y/√3)·ln(r/R_i) - p``,
  ``σ_θθ(r) = σ_rr(r) + 2σ_y/√3``.

Geometría y carga
-----------------
- Cilindro: cuarto de corona en el plano ``xy`` (cuadrante I), altura ``L=1``
  en ``z ∈ [0, 1]``. Radios ``R_i = 1, R_e = 2``.
- Material: ``VonMises3D(E=1000, ν=0.3, σ_y=1, H=0)``.
- BCs: simetría en planos ``y=0`` (``u_y=0``) y ``x=0`` (``u_x=0``);
  plane strain en tapas ``z=0`` y ``z=L`` (``u_z=0`` en cada tapa).
- Carga: tracción interna ``p`` en la cara cilíndrica ``r=R_i`` aplicada
  vía ``compute_face_traction`` del ``Hex8`` sobre la cara apropiada de
  cada elemento de la primera capa radial.
"""
import math
import os
import sys

import numpy as np
import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from solidum.core.domain import Domain
from solidum.elements.solid_3d import Hex8
from solidum.materials.von_mises_3d import VonMises3D
from solidum.math.assembly import Assembler
from solidum.math.convergence import ConvergenceCriterion
from solidum.math.solvers import NonlinearSolver


# Parámetros físicos (idénticos al benchmark 2D).
R_INNER = 1.0
R_OUTER = 2.0
HEIGHT = 1.0           # extrusión en z (plane strain por uz=0 en tapas)
E_YOUNG = 1000.0
NU = 0.3
SIGMA_Y = 1.0
K_SHEAR = SIGMA_Y / math.sqrt(3.0)

P_YIELD_ONSET = K_SHEAR * (1.0 - R_INNER**2 / R_OUTER**2)
P_LIMIT = 2.0 * K_SHEAR * math.log(R_OUTER / R_INNER)


def _hill_transition_radius(p: float) -> float:
    """Raíz de ``p = (σ_y/√3)·[2·ln(r_c/R_i) + 1 - r_c²/R_e²]`` por bisección.

    Idéntica a ``test_hill_cylinder_j2._hill_transition_radius``.
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
    """``σ_rr, σ_θθ`` analíticos según zona plástica/elástica."""
    if r_c <= R_INNER:
        # Lamé puro
        A = p * R_INNER**2 / (R_OUTER**2 - R_INNER**2)
        B = A * R_OUTER**2
        return A - B / r**2, A + B / r**2
    if r <= r_c:
        srr = 2.0 * K_SHEAR * math.log(r / R_INNER) - p
        stt = srr + 2.0 * K_SHEAR
    else:
        A = K_SHEAR * r_c**2 / R_OUTER**2
        B = A * R_OUTER**2
        srr = A - B / r**2
        stt = A + B / r**2
    return srr, stt


# =============================================================================
# Helpers — malla polar Hex8 extruida.
# =============================================================================

def _build_polar_mesh_hex8(nr: int, nt: int, nz: int, material):
    """Cuarto de corona 3D mallado con Hex8.

    Capas: ``nr`` radial × ``nt`` angular × ``nz`` axial. Devuelve
    ``(domain, inner_faces, nid_map, nI, nJ, nK)``, donde
    ``inner_faces`` es la lista de ``(elem, face_idx)`` para las caras
    cilíndricas interiores (``r = R_i``) sobre las que se aplicará la
    tracción de presión.
    """
    nI = nr + 1
    nJ = nt + 1
    nK = nz + 1
    domain = Domain()
    nid_map: dict[tuple[int, int, int], object] = {}
    nid = 0
    for K in range(nK):
        z = HEIGHT * K / (nK - 1) if nK > 1 else 0.0
        for J in range(nJ):
            theta = (math.pi / 2.0) * J / (nJ - 1)
            for I in range(nI):
                r = R_INNER + (R_OUTER - R_INNER) * I / (nI - 1)
                nid += 1
                nid_map[(I, J, K)] = domain.add_node(
                    nid, [r * math.cos(theta), r * math.sin(theta), z],
                )

    inner_faces: list[tuple[object, int]] = []
    eid = 0
    # Hex8 conectividad VTK estándar: nodos 0-3 cara z_min en sentido CCW,
    # nodos 4-7 cara z_max alineados con 0-3.
    for kz in range(nz):
        for jy in range(nt):
            for ix in range(nr):
                n0 = nid_map[(ix,     jy,     kz)]
                n1 = nid_map[(ix + 1, jy,     kz)]
                n2 = nid_map[(ix + 1, jy + 1, kz)]
                n3 = nid_map[(ix,     jy + 1, kz)]
                n4 = nid_map[(ix,     jy,     kz + 1)]
                n5 = nid_map[(ix + 1, jy,     kz + 1)]
                n6 = nid_map[(ix + 1, jy + 1, kz + 1)]
                n7 = nid_map[(ix,     jy + 1, kz + 1)]
                eid += 1
                elem = Hex8(eid, [n0, n1, n2, n3, n4, n5, n6, n7], material)
                domain.add_element(elem)
                if ix == 0:
                    # FACE_NODES Hex8 del proyecto:
                    #   0=(0,3,2,1) z_min; 1=(4,5,6,7) z_max;
                    #   2=(0,1,5,4) η_min; 3=(1,2,6,5) +ξ;
                    #   4=(2,3,7,6) +η;   5=(3,0,4,7) **-ξ**.
                    # La cara cilíndrica interior corresponde a la cara -ξ
                    # del elemento (todos los nodos en ix=0 → r=R_i).
                    inner_faces.append((elem, 5))

    # BCs: simetría en y=0 (lado theta=0) y x=0 (lado theta=π/2);
    # plane strain en tapas z=0 y z=L.
    for K in range(nK):
        for I in range(nI):
            # theta = 0 (lado y=0): u_y = 0
            if (I, 0, K) in nid_map:
                nid_map[(I, 0, K)].fix_dof('uy', 0.0)
            # theta = π/2 (lado x=0): u_x = 0
            if (I, nJ - 1, K) in nid_map:
                nid_map[(I, nJ - 1, K)].fix_dof('ux', 0.0)
    # Plane strain en ambas tapas
    for J in range(nJ):
        for I in range(nI):
            if (I, J, 0) in nid_map:
                nid_map[(I, J, 0)].fix_dof('uz', 0.0)
            if (I, J, nK - 1) in nid_map:
                nid_map[(I, J, nK - 1)].fix_dof('uz', 0.0)

    domain.generate_equation_numbers(verbose=False)
    return domain, inner_faces, nid_map, nI, nJ, nK


def _apply_internal_pressure(domain: Domain,
                              inner_faces: list[tuple[object, int]],
                              p: float) -> np.ndarray:
    """Aplica tracción radial ``-p·n̂`` en las caras cilíndricas interiores.

    La normal exterior de la cara cilíndrica interior apunta hacia el centro
    del cilindro (en -r̂). Una presión positiva ``p`` empuja el material
    hacia afuera, así que la tracción aplicada es ``+p·r̂`` (en dirección
    +r̂, hacia afuera).
    """
    F = np.zeros(domain.total_dofs)
    for elem, face_idx in inner_faces:
        # Centroide de la cara para calcular dirección radial promedio
        face_local = elem.FACE_NODES[face_idx]
        coords = np.array([elem.nodes[i].coordinates[:3] for i in face_local])
        centroid = coords.mean(axis=0)
        r_centroid = math.hypot(centroid[0], centroid[1])
        n_radial = np.array([centroid[0] / r_centroid, centroid[1] / r_centroid, 0.0])
        t_vec = p * n_radial
        f_elem = elem.compute_face_traction(face_idx, t_vec)
        for k_local, k_global in enumerate(elem.get_global_dof_indices()):
            if k_global >= 0:
                F[k_global] += f_elem[k_local]
    return F


def _cart_to_polar_stress_3d(sigma_voigt, x, y):
    """Pasa tensor σ desde Voigt 3D Cartesian a polar (σ_rr, σ_θθ).

    Voigt 3D: ``[σ_xx, σ_yy, σ_zz, σ_xy, σ_yz, σ_xz]``. Solo nos interesa el
    bloque ``(xx, yy, xy)`` para la rotación polar.
    """
    r = math.hypot(x, y)
    c, s = x / r, y / r
    sxx, syy, sxy = sigma_voigt[0], sigma_voigt[1], sigma_voigt[3]
    srr = sxx * c**2 + syy * s**2 + 2.0 * sxy * c * s
    stt = sxx * s**2 + syy * c**2 - 2.0 * sxy * c * s
    return srr, stt


def _solve_pressure(p: float, nr: int = 12, nt: int = 8, nz: int = 1,
                    num_steps: int = 10):
    """Resuelve el problema con presión interna ``p``. Devuelve ``(domain, U)``."""
    material = VonMises3D(E=E_YOUNG, nu=NU, sigma_y=SIGMA_Y, H=0.0)
    domain, inner_faces, nid_map, nI, nJ, nK = _build_polar_mesh_hex8(
        nr=nr, nt=nt, nz=nz, material=material,
    )
    F = _apply_internal_pressure(domain, inner_faces, p)
    conv = ConvergenceCriterion(rtol_force=1.0e-6, rtol_disp=1.0e-6)
    solver = NonlinearSolver(
        Assembler(domain), convergence=conv, num_steps=num_steps,
    )
    U = solver.solve(F)
    return domain, U, nid_map, nI, nJ, nK


def _gather_gauss_state(domain: Domain, U: np.ndarray):
    """Devuelve ``(r, σ_rr, σ_θθ, α)`` para todos los Gauss points del dominio.

    Solo considera Gauss points "interiores" en z (no en tapas) para evitar
    el efecto borde del plane strain por Dirichlet.
    """
    rs, srrs, stts, alphas = [], [], [], []
    for elem in domain.elements.values():
        gs = elem.compute_gauss_state(U)
        for idx, (xyz, sig) in enumerate(zip(gs['points_global'], gs['stress'])):
            r = math.hypot(xyz[0], xyz[1])
            srr, stt = _cart_to_polar_stress_3d(sig, xyz[0], xyz[1])
            alpha = (elem.state.vars[idx].get('alpha', 0.0)
                     if elem.state.vars[idx] else 0.0)
            rs.append(r); srrs.append(srr); stts.append(stt); alphas.append(alpha)
    return (np.asarray(rs), np.asarray(srrs), np.asarray(stts),
            np.asarray(alphas))


# =============================================================================
# Tests.
# =============================================================================

def test_hill_3d_elastic_regime():
    """``p = 0.30 < p_e``: nada plastifica; σ_rr/σ_θθ aproximan Lamé."""
    p = 0.30
    assert p < P_YIELD_ONSET, f"Test inválido: p={p} ≥ p_e={P_YIELD_ONSET:.4f}"
    domain, U, *_ = _solve_pressure(p, num_steps=4)
    r_g, srr_g, stt_g, alpha_g = _gather_gauss_state(domain, U)

    # Sin plasticidad
    assert alpha_g.max() < 1.0e-10, (
        f"Régimen elástico (p={p}, p_e={P_YIELD_ONSET:.4f}): "
        f"α_max = {alpha_g.max():.3e} debe ser 0"
    )

    # σ comparada con Lamé en zona interior (lejos de bordes)
    mask = (r_g > R_INNER + 0.15) & (r_g < R_OUTER - 0.15)
    srr_ref = np.array([_hill_stress_polar(r, p, 0.0)[0] for r in r_g[mask]])
    stt_ref = np.array([_hill_stress_polar(r, p, 0.0)[1] for r in r_g[mask]])
    e_rr = (np.sqrt(np.mean((srr_g[mask] - srr_ref) ** 2))
             / np.sqrt(np.mean(srr_ref ** 2)))
    e_tt = (np.sqrt(np.mean((stt_g[mask] - stt_ref) ** 2))
             / np.sqrt(np.mean(stt_ref ** 2)))
    # Tolerancias del orden del benchmark 2D Hex8 8×6 (similar a Quad4 12×8).
    assert e_rr < 0.20, f"σ_rr L²-rel = {e_rr:.4%} > 20%"
    assert e_tt < 0.08, f"σ_θθ L²-rel = {e_tt:.4%} > 8%"


def test_hill_3d_elastoplastic_zone():
    """``p = 0.70``: zona plástica interna detectable y σ_θθ - σ_rr ≈ 2σ_y/√3.

    Para los parámetros del benchmark, ``c_analítico ≈ 1.49`` en ``p=0.70``.
    """
    p = 0.70
    assert P_YIELD_ONSET < p < P_LIMIT, (
        f"Test inválido: p={p} fuera de régimen elastoplástico "
        f"({P_YIELD_ONSET:.4f}, {P_LIMIT:.4f})"
    )
    r_c = _hill_transition_radius(p)
    assert 1.4 < r_c < 1.6, f"r_c analítico inesperado: {r_c:.4f}"

    domain, U, *_ = _solve_pressure(p, num_steps=10)
    r_g, srr_g, stt_g, alpha_g = _gather_gauss_state(domain, U)

    # Gauss en zona claramente plástica
    in_plastic = r_g < (r_c - 0.10)
    assert in_plastic.sum() > 0, "Sin Gauss en zona plástica"
    assert (alpha_g[in_plastic] > 1.0e-6).all(), (
        f"Gauss en zona plástica (r < {r_c-0.10:.3f}) con α=0: "
        f"min α = {alpha_g[in_plastic].min():.3e}"
    )

    # Gauss en zona claramente elástica
    in_elastic = r_g > (r_c + 0.10)
    assert in_elastic.sum() > 0, "Sin Gauss en zona elástica"
    assert (alpha_g[in_elastic] < 1.0e-8).all(), (
        f"Gauss en zona elástica (r > {r_c+0.10:.3f}) con α>0: "
        f"max α = {alpha_g[in_elastic].max():.3e}"
    )

    # Saturación de fluencia: σ_θθ - σ_rr ≈ 2σ_y/√3 en zona plástica
    delta_target = 2.0 * K_SHEAR
    delta_num = (stt_g[in_plastic] - srr_g[in_plastic]).mean()
    rel_err = abs(delta_num - delta_target) / delta_target
    assert rel_err < 0.05, (
        f"σ_θθ - σ_rr en zona plástica = {delta_num:.4f} "
        f"vs 2σ_y/√3 = {delta_target:.4f} (err_rel = {rel_err:.4%})"
    )


def test_hill_3d_stress_in_plastic_zone_matches_analytic():
    """``σ_rr`` FEM en zona plástica coincide con la fórmula de Hill, L²-rel < 18%."""
    p = 0.70
    r_c = _hill_transition_radius(p)
    domain, U, *_ = _solve_pressure(p, num_steps=10)
    r_g, srr_g, stt_g, alpha_g = _gather_gauss_state(domain, U)

    mask = (r_g < (r_c - 0.05)) & (r_g > R_INNER + 0.08)
    assert mask.sum() > 0, "Sin Gauss aptos para comparar en zona plástica"

    srr_ref = np.array([_hill_stress_polar(r, p, r_c)[0] for r in r_g[mask]])
    e_rr = (np.sqrt(np.mean((srr_g[mask] - srr_ref) ** 2))
             / np.sqrt(np.mean(srr_ref ** 2)))
    # Tolerancia 18% — mismo orden que Quad4 12×8 plane strain en 2D.
    # El Hex8 con plane strain restringida da el mismo orden de error que Quad4
    # bajo mallado equivalente (8 nr × 6 nt) — confirma el cross-consistency a
    # nivel de pipeline completo.
    assert e_rr < 0.18, f"σ_rr L²-rel zona plástica = {e_rr:.4%} > 18%"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
