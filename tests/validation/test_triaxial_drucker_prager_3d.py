"""Triaxial Drucker-Prager 3D — verificación cuantitativa de la superficie del cono.

Referencias
-----------
- Drucker D.C., Prager W. (1952). Soil mechanics and plastic analysis or limit
  design. Quart. Appl. Math. 10, 157-165.
- Chen W.F., Han D.J. (1988). *Plasticity for Structural Engineers*, Springer,
  §7 (cohesivo-friccionales).
- de Souza Neto E.A., Perić D., Owen D.R.J. (2008). *Computational Methods for
  Plasticity*, Wiley, §8.4 (calibración 3D outer/inner cone vs Mohr-Coulomb).

Concepto
--------
1 elemento ``Hex8`` con tres planos de simetría como rodillos y carga aplicada
en las caras opuestas. El estado de esfuerzos es uniforme; el cubo tiene caras
de área unidad, así que las tracciones aplicadas corresponden directamente a
las componentes diagonales del tensor.

Esto permite ejercer trayectorias controladas en el espacio de invariantes
``(I_1, √J_2)`` y verificar que el modelo numérico cruza la superficie del
cono Drucker-Prager exactamente cuando

.. math::

   f(σ, α) = √J_2 + η_f·I_1 - k(α) = 0

con ``η_f`` y ``k_0`` derivados de los parámetros físicos ``(c_0, φ)`` según
la variante de calibración elegida.

Decisión de control (force vs displacement)
-------------------------------------------
**Force control** (tracción aplicada) en los tests sobre la rama regular del
cono. Requiere ``H > 0`` para que la capacidad plástica crezca y el solver
pueda equilibrar cargas ligeramente por encima del yield inicial; con
``H = 0`` el cono pone un techo duro a la capacidad tensil y el solver no
encuentra equilibrio (limitación física, no del algoritmo).

**Displacement control** en el test del ápice. El tangente algorítmico al
ápice es rank-1 (solo volumétrico); las modos cortantes libres del cubo
quedan en el kernel de la matriz tangente → sistema globalmente singular
en force control. Prescribir todos los DOFs vía Dirichlet evita el solve
del sistema singular y permite leer σ del estado committed para compararlo
contra la fórmula analítica del ápice.

Validaciones
------------
1. **Yield onset uniaxial tensión** (force control, H > 0): yield exacto en
   σ_1 = k_0 / (1/√3 + η_f). Para outer_cone, φ=30°, c_0=10: σ_1_yield ≈ 14.85.

2. **Apex pressure formula** (displacement control, H > 0): aplicar expansión
   hidrostática prescrita; verificar σ_avg = k(α)/(3·η_f) tras converger.

3. **Calibración outer vs inner** (force control, H > 0): misma carga uniaxial
   tensil entre los dos yields predicchos; outer permanece elástico, inner
   plastifica. Cada yield onset coincide con su fórmula del cono.

4. **Compresión hidrostática sin yield** (force control, H = 0): σ_1 = σ_2 =
   σ_3 = -p (compresión isotrópica) **no produce yield** independientemente
   de la magnitud porque η_f·I_1 < 0 ayuda a cumplir f ≤ 0.

5. **Dilatancia plástica** (force control, H > 0): tras alcanzar yield bajo
   carga uniaxial tensil, ``tr(ε^p) = 3·η_g·α`` en cada Gauss point —
   invariante cinemática de la regla de flujo no asociada.
"""
import math
import os
import sys

import numpy as np
import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from solidum.core.domain import Domain
from solidum.elements.solid_3d import Hex8
from solidum.materials.drucker_prager_3d import DruckerPrager3D
from solidum.math.assembly import Assembler
from solidum.math.convergence import ConvergenceCriterion
from solidum.math.solvers import NonlinearSolver


# Parámetros físicos del benchmark.
E_YOUNG = 2.0e4
NU = 0.2
COHESION = 10.0
PHI_DEG = 30.0
PSI_DEG = 10.0   # no asociada (η_g < η_f)
HARDENING = 1.0e3   # H_k > 0 para que la capacidad plástica crezca con α

# Cubo unitario Hex8.
HEX8_UNIT_COORDS = [
    (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
]
HEX8_FACE_XMIN = (0, 3, 4, 7)
HEX8_FACE_XMAX = (1, 2, 5, 6)
HEX8_FACE_YMIN = (0, 1, 4, 5)
HEX8_FACE_YMAX = (2, 3, 6, 7)
HEX8_FACE_ZMIN = (0, 1, 2, 3)
HEX8_FACE_ZMAX = (4, 5, 6, 7)


def _build_triaxial_cube(material):
    """Cubo unitario Hex8 con tres rodillos simétricos (x=0, y=0, z=0)."""
    dom = Domain()
    nodes = [dom.add_node(i + 1, list(c)) for i, c in enumerate(HEX8_UNIT_COORDS)]
    elem = Hex8(1, nodes, material)
    dom.add_element(elem)

    for i in HEX8_FACE_XMIN:
        nodes[i].fix_dof('ux', 0.0)
    for i in HEX8_FACE_YMIN:
        nodes[i].fix_dof('uy', 0.0)
    for i in HEX8_FACE_ZMIN:
        nodes[i].fix_dof('uz', 0.0)

    dom.generate_equation_numbers(verbose=False)
    return dom, elem, nodes


def _apply_triaxial_traction(dom, nodes, sigma_xx, sigma_yy, sigma_zz):
    """Aplica tracción uniforme en las tres caras +x, +y, +z (área unidad)."""
    F = np.zeros(dom.total_dofs)
    for i in HEX8_FACE_XMAX:
        F[nodes[i].dofs['ux']] = sigma_xx / 4.0
    for i in HEX8_FACE_YMAX:
        F[nodes[i].dofs['uy']] = sigma_yy / 4.0
    for i in HEX8_FACE_ZMAX:
        F[nodes[i].dofs['uz']] = sigma_zz / 4.0
    return F


def _solve_triaxial_force(material, sigma_xx, sigma_yy, sigma_zz, num_steps=5):
    """Resuelve el problema force-controlled triaxial. Devuelve (elem, U)."""
    dom, elem, nodes = _build_triaxial_cube(material)
    F = _apply_triaxial_traction(dom, nodes, sigma_xx, sigma_yy, sigma_zz)
    conv = ConvergenceCriterion(rtol_force=1.0e-9, rtol_disp=1.0e-9)
    solver = NonlinearSolver(
        Assembler(dom), convergence=conv, num_steps=num_steps,
    )
    U = solver.solve(F)
    return elem, U


def _build_apex_displacement_cube(material, eps_iso):
    """Cubo unitario con expansión hidrostática prescrita (Dirichlet en todos los DOFs).

    Cada nodo se desplaza desde el centro (0.5, 0.5, 0.5) en una cantidad
    proporcional a (pos - center) por el factor eps_iso. Bajo deformación
    uniforme isotrópica, el estado de esfuerzos resultante es hidrostático
    σ_xx = σ_yy = σ_zz = σ_iso. Si eps_iso es suficientemente grande, el
    estado cae en la rama de ápice del cono DP.
    """
    dom = Domain()
    nodes = [dom.add_node(i + 1, list(c)) for i, c in enumerate(HEX8_UNIT_COORDS)]
    elem = Hex8(1, nodes, material)
    dom.add_element(elem)

    center = np.array([0.5, 0.5, 0.5])
    for n in nodes:
        pos = np.array(n.coordinates[:3], dtype=float)
        u = eps_iso * (pos - center)
        n.fix_dof('ux', u[0])
        n.fix_dof('uy', u[1])
        n.fix_dof('uz', u[2])

    dom.generate_equation_numbers(verbose=False)
    return dom, elem


# =============================================================================
# Test 1 — Yield onset uniaxial tensión, outer_cone.
# =============================================================================

def test_yield_onset_uniaxial_tension_outer_cone():
    """σ_1 > 0, σ_2 = σ_3 = 0: yield exacto en σ_1 = k_0 / (1/√3 + η_f).

    Bajo este path I_1 = σ_1, J_2 = σ_1²/3 ⇒ √J_2 = σ_1/√3.
    Criterio f = σ_1·(1/√3 + η_f) - k_0 = 0 ⇒ σ_1_yield = k_0/(1/√3 + η_f).
    Con ``H > 0``, la capacidad plástica crece con α y el solver puede
    equilibrar cargas ligeramente por encima del yield inicial.
    """
    mat = DruckerPrager3D(
        E=E_YOUNG, nu=NU, cohesion=COHESION,
        phi_deg=PHI_DEG, psi_deg=PSI_DEG, H=HARDENING,
        variant='outer_cone',
    )
    sigma_yield_analytical = mat.k0 / (1.0 / math.sqrt(3.0) + mat.eta_f)
    # Sanity check del valor analítico para los parámetros del benchmark
    assert abs(sigma_yield_analytical - 14.85) < 0.5, (
        f"σ_yield analítico esperado ≈ 14.85, obtenido {sigma_yield_analytical:.4f}"
    )

    # 0.95·σ_yield → régimen elástico inequívoco
    elem_below, _ = _solve_triaxial_force(
        mat, 0.95 * sigma_yield_analytical, 0.0, 0.0,
    )
    for k, state in enumerate(elem_below.state.vars):
        assert state['alpha'] == 0.0, (
            f"GP {k} debajo del yield (σ_1 = 0.95·σ_yield = {0.95*sigma_yield_analytical:.4f}): "
            f"α = {state['alpha']:.3e} ≠ 0"
        )

    # 1.05·σ_yield → régimen plástico (capacidad ajustada por H > 0)
    elem_above, _ = _solve_triaxial_force(
        mat, 1.05 * sigma_yield_analytical, 0.0, 0.0,
    )
    for k, state in enumerate(elem_above.state.vars):
        assert state['alpha'] > 0.0, (
            f"GP {k} por encima del yield (σ_1 = 1.05·σ_yield = {1.05*sigma_yield_analytical:.4f}): "
            f"α = {state['alpha']:.3e} debe ser > 0"
        )

    # Estimación analítica de α tras equilibrar 1.05·σ_yield:
    # σ_1 = (k_0 + H·α)/(1/√3 + η_f) ⇒ α = [σ_1·(1/√3+η_f) - k_0] / H
    alpha_expected = (1.05 * sigma_yield_analytical
                      * (1.0 / math.sqrt(3.0) + mat.eta_f) - mat.k0) / mat.H
    alpha_fem = np.mean([s['alpha'] for s in elem_above.state.vars])
    rel_err = abs(alpha_fem - alpha_expected) / alpha_expected
    assert rel_err < 0.05, (
        f"α FEM = {alpha_fem:.4e} vs analítico = {alpha_expected:.4e} "
        f"(err_rel = {rel_err:.4%} > 5%)"
    )


# =============================================================================
# Test 2 — Apex pressure formula vía displacement control.
# =============================================================================

def test_apex_pressure_formula_via_displacement_control():
    """Expansión hidrostática prescrita: σ_avg = k(α)/(3·η_f) tras converger.

    Con prescripción de TODOS los DOFs por Dirichlet, el sistema no requiere
    solve del K global (singular en apex por tangente rank-1). El kernel
    devuelve el σ correcto desde la formulación del return mapping al ápice.
    Bajo expansión uniforme ε_iso > p_apex_init/(3·K), el estado cae en
    ápice y σ_avg debe coincidir con k(α_committed)/(3·η_f).
    """
    mat = DruckerPrager3D(
        E=E_YOUNG, nu=NU, cohesion=COHESION,
        phi_deg=PHI_DEG, psi_deg=PSI_DEG, H=HARDENING,
        variant='outer_cone',
    )

    # Expansión grande para asegurar entrada al ápice
    eps_iso = 5.0e-3
    dom, elem = _build_apex_displacement_cube(mat, eps_iso)
    assembler = Assembler(dom)

    F_ext = np.zeros(dom.total_dofs)
    conv = ConvergenceCriterion(rtol_force=1.0e-7, rtol_disp=1.0e-7)
    solver = NonlinearSolver(assembler, convergence=conv, num_steps=6)
    U = solver.solve(F_ext)

    # Estado committed: leer σ_avg y α_avg desde los 8 Gauss points
    gs = elem.compute_gauss_state(U)
    for k in range(8):
        sigma_gp = gs['stress'][k]
        # σ_xx ≈ σ_yy ≈ σ_zz (hidrostático)
        avg_diag = (sigma_gp[0] + sigma_gp[1] + sigma_gp[2]) / 3.0
        for j in range(3):
            rel_err = abs(sigma_gp[j] - avg_diag) / abs(avg_diag)
            assert rel_err < 1.0e-6, (
                f"GP {k} no hidrostático: σ_{['xx','yy','zz'][j]} = {sigma_gp[j]:.4e}, "
                f"avg = {avg_diag:.4e}, err_rel = {rel_err:.3e}"
            )
        # Cortantes ≈ 0
        for j in range(3, 6):
            self_norm = max(abs(avg_diag), 1.0)
            assert abs(sigma_gp[j]) / self_norm < 1.0e-6, (
                f"GP {k}: cortante σ_{j} = {sigma_gp[j]:.4e} debería ser ≈ 0"
            )

    # Verificar fórmula del ápice: σ_avg = k(α)/(3·η_f)
    for k, state in enumerate(elem.state.vars):
        alpha = state['alpha']
        assert alpha > 0.0, f"GP {k}: α = {alpha:.3e} debe ser > 0 en apex"
        k_curr = mat.k0 + mat.H * alpha
        sigma_apex_expected = k_curr / (3.0 * mat.eta_f)
        sigma_gp_avg = (gs['stress'][k][0] + gs['stress'][k][1] + gs['stress'][k][2]) / 3.0
        rel_err = abs(sigma_gp_avg - sigma_apex_expected) / sigma_apex_expected
        assert rel_err < 1.0e-6, (
            f"GP {k}: σ_avg FEM = {sigma_gp_avg:.6e} vs k(α)/(3η_f) = {sigma_apex_expected:.6e} "
            f"(α = {alpha:.4e}, err_rel = {rel_err:.3e})"
        )


# =============================================================================
# Test 3 — Outer cone vs inner cone: yields distintos.
# =============================================================================

def test_outer_vs_inner_cone_yield_difference():
    """Carga uniaxial intermedia entre los dos yields: inner plastifica, outer no.

    Outer cone tiene η_f y k_0 mayores → σ_yield mayor (cono más ancho en
    tensión). Inner cone es geométricamente inscrito → cono más conservador.
    """
    mat_outer = DruckerPrager3D(
        E=E_YOUNG, nu=NU, cohesion=COHESION,
        phi_deg=PHI_DEG, psi_deg=PSI_DEG, H=HARDENING,
        variant='outer_cone',
    )
    mat_inner = DruckerPrager3D(
        E=E_YOUNG, nu=NU, cohesion=COHESION,
        phi_deg=PHI_DEG, psi_deg=PSI_DEG, H=HARDENING,
        variant='inner_cone',
    )

    sigma_yield_outer = mat_outer.k0 / (1.0 / math.sqrt(3.0) + mat_outer.eta_f)
    sigma_yield_inner = mat_inner.k0 / (1.0 / math.sqrt(3.0) + mat_inner.eta_f)
    assert sigma_yield_outer > sigma_yield_inner, (
        f"σ_yield_outer ({sigma_yield_outer:.4f}) debe ser > "
        f"σ_yield_inner ({sigma_yield_inner:.4f})"
    )

    # Carga intermedia con margen amplio respecto a ambos yields para
    # absorber cualquier ruido del integrador 8-Gauss.
    sigma_load = 0.5 * (sigma_yield_outer + sigma_yield_inner)
    assert sigma_load > 1.10 * sigma_yield_inner, "Margen mínimo del inner yield insuficiente"
    assert sigma_load < 0.90 * sigma_yield_outer, "Margen mínimo del outer yield insuficiente"

    # Inner cone: σ_load > σ_yield_inner → plástico
    elem_inner, _ = _solve_triaxial_force(mat_inner, sigma_load, 0.0, 0.0)
    for k, state in enumerate(elem_inner.state.vars):
        assert state['alpha'] > 0.0, (
            f"Inner cone GP {k}: σ_load = {sigma_load:.4f} > σ_yield_inner = "
            f"{sigma_yield_inner:.4f}, α = {state['alpha']:.3e} debe ser > 0"
        )

    # Outer cone: σ_load < σ_yield_outer → elástico
    elem_outer, _ = _solve_triaxial_force(mat_outer, sigma_load, 0.0, 0.0)
    for k, state in enumerate(elem_outer.state.vars):
        assert state['alpha'] == 0.0, (
            f"Outer cone GP {k}: σ_load = {sigma_load:.4f} < σ_yield_outer = "
            f"{sigma_yield_outer:.4f}, α = {state['alpha']:.3e} debe ser 0"
        )


# =============================================================================
# Test 4 — Compresión hidrostática sin yield.
# =============================================================================

def test_hydrostatic_compression_no_yield():
    """σ_1 = σ_2 = σ_3 = -p (compresión isótropa): no produce yield para
    ningún p > 0, porque I_1 < 0 ⇒ η_f·I_1 < 0 ⇒ f < 0 sin cortante.

    Es propiedad fundamental del cono DP: la compresión hidrostática es
    siempre admisible (en ausencia de cap). El material es "infinitamente
    resistente" a la compresión pura isotrópica.
    """
    mat = DruckerPrager3D(
        E=E_YOUNG, nu=NU, cohesion=COHESION,
        phi_deg=PHI_DEG, psi_deg=PSI_DEG, H=0.0,   # H=0 OK: no fuerza yield
        variant='outer_cone',
    )

    # Compresión muy grande (10× el k_0 esperado en magnitud)
    p_huge = 100.0 * mat.k0
    elem, _ = _solve_triaxial_force(
        mat, -p_huge, -p_huge, -p_huge, num_steps=4,
    )
    for k, state in enumerate(elem.state.vars):
        assert state['alpha'] == 0.0, (
            f"GP {k} bajo compresión hidrostática p={p_huge:.1f}: "
            f"α = {state['alpha']:.3e} ≠ 0; el cono DP no debería yield "
            f"bajo compresión isótropa pura"
        )


# =============================================================================
# Test 5 — Invariante dilatancia tras yield.
# =============================================================================

def test_dilatancy_invariant_after_yield():
    """Tras yield bajo carga uniaxial tensil, tr(ε^p) = 3·η_g·α en el Gauss point.

    Invariante kinemática de la regla de flujo no asociada del DP:
    ``∂g/∂σ`` tiene parte volumétrica ``η_g·I`` que se acumula en
    ``tr(ε^p)``. Cada incremento ``Δγ`` añade ``3·η_g·Δγ`` a la traza y
    suma ``Δγ`` a α; por tanto ``tr(ε^p) = 3·η_g·α`` en todo estado plástico.
    """
    mat = DruckerPrager3D(
        E=E_YOUNG, nu=NU, cohesion=COHESION,
        phi_deg=PHI_DEG, psi_deg=PSI_DEG, H=HARDENING,
        variant='outer_cone',
    )
    sigma_yield = mat.k0 / (1.0 / math.sqrt(3.0) + mat.eta_f)
    sigma_load = 1.10 * sigma_yield   # plásticamente por encima del yield

    elem, _ = _solve_triaxial_force(mat, sigma_load, 0.0, 0.0, num_steps=8)

    for k, state in enumerate(elem.state.vars):
        assert state['alpha'] > 0.0, f"GP {k}: α = {state['alpha']:.3e} debe ser > 0"
        eps_p = state['eps_p']
        tr_eps_p = eps_p[0] + eps_p[1] + eps_p[2]
        expected = 3.0 * mat.eta_g * state['alpha']
        rel_err = abs(tr_eps_p - expected) / (abs(expected) + 1.0e-30)
        assert rel_err < 1.0e-8, (
            f"GP {k}: tr(ε^p) = {tr_eps_p:.6e} ≠ 3·η_g·α = {expected:.6e} "
            f"(err_rel = {rel_err:.3e})"
        )


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
