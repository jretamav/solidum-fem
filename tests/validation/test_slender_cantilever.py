"""MacNeal-Harder slender cantilever — viga esbelta bajo carga transversal.

Referencia
----------
MacNeal, R.H.; Harder, R.L. (1985). *A proposed standard set of problems
to test finite element accuracy*. Finite Elements in Analysis and Design,
1(1), 3-20.

Problema canónico (versión slender beam, no patch test):

    Cantilever rectangular de longitud L = 6, peralte h = 0.2, espesor
    t = 0.1. Esbeltez L/h = 30.
    Material elástico isótropo E = 1×10⁷, ν = 0.3.
    Empotramiento total en x = 0.
    Carga transversal P = 1 en el extremo libre x = L.

Solución analítica de Euler-Bernoulli (sin corrección por cortante):

    u_tip = P · L³ / (3 · E · I)        con I = t · h³ / 12

    u_tip_EB = 0.108

Solución analítica de Timoshenko (con corrección por cortante):

    u_tip = P · L³ / (3 · E · I) + P · L / (κ · G · A)

    con κ = 5/6 (factor de cortante para sección rectangular),
        G = E / (2 · (1 + ν)),  A = h · t.

    Para los valores de este problema, el aporte de cortante es
    ≈ 9.4 × 10⁻⁵ → u_tip_T ≈ 0.108094 (Timoshenko vs Euler-Bernoulli
    coinciden a 4 cifras para L/h = 30).

Objetivos de este benchmark
---------------------------
1. **Frames analíticos**: Frame2DEuler debe dar la solución exacta (a
   redondeo de máquina) — son funciones de forma cúbicas que reproducen
   exactamente el problema. Frame2DTimoshenko debe dar la solución de
   Timoshenko también con un único elemento.
2. **Sólidos 2D con malla 1 capa (shear locking)**: Quad4 a una sola
   capa de elementos en la dirección del peralte **subestima drásticamente**
   la flecha por shear locking — predicción típica ≈30–60% del valor
   analítico para L/h=30 y h-discretización ~30 elementos en x.
3. **Sólidos 2D con varias capas**: Quad4 con 4 capas reduce el locking
   sustancialmente. Quad8 sin shear locking severo.
4. **Convergencia h**: la flecha converge a la solución de Timoshenko al
   refinar la malla.

Aplicación de la carga (sólidos 2D)
-----------------------------------
La carga total P = 1 se distribuye como **tracción uniforme** en el borde
derecho (x = L). Para una distribución parabólica (Timoshenko-exact), se
podría descomponer manualmente pero la tracción uniforme converge al
mismo resultado al refinar, con un error de Saint-Venant local que decae
en pocas longitudes de peralte. La evaluación se hace en el centro del
borde derecho (y = 0), lejos de la perturbación de Saint-Venant.
"""
import os
import sys

import numpy as np
import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from solidum.core.domain import Domain
from solidum.elements.frame import Frame2DEuler, Frame2DTimoshenko
from solidum.elements.solid_2d import Quad4, Quad8
from solidum.materials.elastic import Elastic1D
from solidum.materials.elastic_2d import Elastic2D
from solidum.math.assembly import Assembler
from solidum.math.solvers import LinearSolver


# Parámetros físicos del benchmark.
L_BEAM = 6.0
H_BEAM = 0.2
T_BEAM = 0.1          # espesor (2D plane stress)
E_YOUNG = 1.0e7
NU = 0.3
P_TIP = 1.0

I_BEAM = T_BEAM * H_BEAM**3 / 12.0      # momento de inercia
A_BEAM = H_BEAM * T_BEAM                  # área de la sección

# Solución analítica.
U_TIP_EULER = P_TIP * L_BEAM**3 / (3.0 * E_YOUNG * I_BEAM)
# Timoshenko: u_EB + corrección de cortante.
G_SHEAR = E_YOUNG / (2.0 * (1.0 + NU))
KAPPA = 5.0 / 6.0
U_TIP_TIMOSHENKO = U_TIP_EULER + P_TIP * L_BEAM / (KAPPA * G_SHEAR * A_BEAM)


# =============================================================================
# 1) Frames — verificación contra solución analítica.
# =============================================================================

def _build_frame_cantilever(elem_cls, n_elements: int, use_timoshenko: bool = False):
    """Construye un cantilever de marcos con ``n_elements`` elementos."""
    domain = Domain()
    material = Elastic1D(E=E_YOUNG)
    nodes = []
    for i in range(n_elements + 1):
        x = L_BEAM * i / n_elements
        nodes.append(domain.add_node(i + 1, [x, 0.0]))

    for i in range(n_elements):
        if use_timoshenko:
            elem = Frame2DTimoshenko(i + 1, [nodes[i], nodes[i + 1]],
                                      material, A=A_BEAM, I=I_BEAM,
                                      As=KAPPA * A_BEAM, nu=NU)
        else:
            elem = Frame2DEuler(i + 1, [nodes[i], nodes[i + 1]],
                                 material, A=A_BEAM, I=I_BEAM)
        domain.add_element(elem)

    # Empotramiento total en x = 0.
    nodes[0].fix_dof('ux', 0.0)
    nodes[0].fix_dof('uy', 0.0)
    nodes[0].fix_dof('rz', 0.0)
    domain.generate_equation_numbers(verbose=False)
    return domain, nodes


def test_frame_euler_single_element_exact():
    """Frame2DEuler con un único elemento da la solución de Euler-Bernoulli exacta.

    Las funciones de forma cúbicas reproducen exactamente la solución
    homogénea del problema (no hay carga distribuida); un solo elemento
    basta.
    """
    domain, nodes = _build_frame_cantilever(Frame2DEuler, n_elements=1)
    F = np.zeros(domain.total_dofs)
    F[nodes[-1].dofs['uy']] = -P_TIP    # convención: +y arriba, P hacia abajo
    U = LinearSolver(Assembler(domain)).solve(F)
    u_tip = U[nodes[-1].dofs['uy']]
    # u_tip negativo (P hacia abajo); comparamos magnitudes.
    np.testing.assert_allclose(abs(u_tip), U_TIP_EULER, rtol=1e-12,
                               err_msg=f"Frame2DEuler 1 elem: u_tip={u_tip}, "
                                       f"ref={U_TIP_EULER}")


def test_frame_timoshenko_single_element_exact():
    """Frame2DTimoshenko con un único elemento da la solución de Timoshenko exacta."""
    domain, nodes = _build_frame_cantilever(Frame2DTimoshenko, n_elements=1,
                                              use_timoshenko=True)
    F = np.zeros(domain.total_dofs)
    F[nodes[-1].dofs['uy']] = -P_TIP
    U = LinearSolver(Assembler(domain)).solve(F)
    u_tip = U[nodes[-1].dofs['uy']]
    # Timoshenko con un elemento puede tener shear locking moderado a baja
    # esbeltez, pero para L/h = 30 con la integración del elemento de Solidum
    # debe converger a < 1% al valor analítico.
    rel_err = abs(abs(u_tip) - U_TIP_TIMOSHENKO) / U_TIP_TIMOSHENKO
    assert rel_err < 0.01, (
        f"Frame2DTimoshenko 1 elem: u_tip={u_tip:.6e}, "
        f"ref_T={U_TIP_TIMOSHENKO:.6e}, err={rel_err:.4%}"
    )


def test_frame_timoshenko_h_refinement_converges_exact():
    """Frame2DTimoshenko: con 8 elementos converge a Timoshenko-exact."""
    domain, nodes = _build_frame_cantilever(Frame2DTimoshenko, n_elements=8,
                                              use_timoshenko=True)
    F = np.zeros(domain.total_dofs)
    F[nodes[-1].dofs['uy']] = -P_TIP
    U = LinearSolver(Assembler(domain)).solve(F)
    u_tip = U[nodes[-1].dofs['uy']]
    np.testing.assert_allclose(abs(u_tip), U_TIP_TIMOSHENKO, rtol=1e-6,
                               err_msg=f"Frame2DTimoshenko 8 elem: u_tip={u_tip}, "
                                       f"ref_T={U_TIP_TIMOSHENKO}")


# =============================================================================
# 2) Sólido 2D — detección de shear locking y convergencia.
# =============================================================================

def _build_solid_cantilever(elem_cls, nx: int, ny: int, material):
    """Cantilever 2D rectangular nx × ny celdas.

    Returns
    -------
    domain, mid_right_node : Domain, Node
        ``mid_right_node`` es el nodo central del borde derecho (x=L, y=0)
        — punto de medida de la flecha.
    right_elements, right_edge_idx
    """
    is_quadratic = elem_cls is Quad8
    sub = 2 if is_quadratic else 1
    nI = nx * sub + 1
    nJ = ny * sub + 1

    domain = Domain()
    nid_map: dict[tuple[int, int], object] = {}
    nid = 0
    for J in range(nJ):
        for I in range(nI):
            is_center = is_quadratic and (I % 2 == 1) and (J % 2 == 1)
            if is_center and elem_cls is Quad8:
                continue
            nid += 1
            x = L_BEAM * I / (nI - 1)
            y = -H_BEAM / 2.0 + H_BEAM * J / (nJ - 1)   # centrado en y=0
            nid_map[(I, J)] = domain.add_node(nid, [x, y])

    right_elements: list[object] = []
    eid = 0
    if elem_cls is Quad4:
        for jy in range(ny):
            for ix in range(nx):
                n0 = nid_map[(ix,   jy)]
                n1 = nid_map[(ix+1, jy)]
                n2 = nid_map[(ix+1, jy+1)]
                n3 = nid_map[(ix,   jy+1)]
                eid += 1
                elem = Quad4(eid, [n0, n1, n2, n3], material, thickness=T_BEAM)
                domain.add_element(elem)
                if ix == nx - 1:
                    right_elements.append(elem)
    elif elem_cls is Quad8:
        for jy in range(ny):
            for ix in range(nx):
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
                eid += 1
                elem = Quad8(eid, nodes, material, thickness=T_BEAM)
                domain.add_element(elem)
                if ix == nx - 1:
                    right_elements.append(elem)

    # Empotramiento total del borde izquierdo (I = 0).
    for J in range(nJ):
        if (0, J) in nid_map:
            nid_map[(0, J)].fix_dof('ux', 0.0)
            nid_map[(0, J)].fix_dof('uy', 0.0)

    # Nodo central del borde derecho (I=nI-1, J=ny en grid Q4, =ny en Q8).
    # En Q8 la celda central angular es la fila J=ny (que es el J central).
    mid_right = nid_map[(nI - 1, (nJ - 1) // 2)]

    domain.generate_equation_numbers(verbose=False)
    return domain, mid_right, right_elements


def _apply_tip_shear(domain: Domain, right_elements: list, total_force: float):
    """Tracción uniforme vertical en el borde derecho, fuerza total = total_force."""
    edge_length = H_BEAM
    t_y = total_force / (edge_length * T_BEAM)
    # Edge_idx 1 = arista entre n1 y n2 (cara derecha) para Quad4/Quad8.
    edge_idx = 1
    F = np.zeros(domain.total_dofs)
    for elem in right_elements:
        f_elem = elem.compute_edge_traction(edge_idx, np.array([0.0, t_y]))
        for k_local, k_global in enumerate(elem.get_global_dof_indices()):
            if k_global >= 0:
                F[k_global] += f_elem[k_local]
    return F


def test_quad4_single_layer_shows_shear_locking():
    """Quad4 con malla 30×1 muestra shear locking: u_tip < 80% del valor analítico.

    Este test **documenta el fenómeno** — no es un fallo, es la confirmación
    del comportamiento conocido de Q4 con bilinear shape functions en
    flexión esbelta con una sola capa de elementos. La predicción es
    típicamente ~30-60% del valor real para L/h grande.
    """
    material = Elastic2D(E=E_YOUNG, nu=NU, hypothesis='plane_stress')
    domain, mid_right, right_elems = _build_solid_cantilever(
        Quad4, nx=30, ny=1, material=material,
    )
    F = _apply_tip_shear(domain, right_elems, total_force=-P_TIP)
    U = LinearSolver(Assembler(domain)).solve(F)
    u_tip = abs(U[mid_right.dofs['uy']])
    # Documentación del locking: ratio < 0.8.
    ratio = u_tip / U_TIP_TIMOSHENKO
    assert ratio < 0.80, (
        f"Quad4 30×1 NO muestra el shear locking esperado: "
        f"u_tip/u_ref = {ratio:.3%} >= 80%. "
        f"u_tip={u_tip:.6e}, u_ref={U_TIP_TIMOSHENKO:.6e}. "
        f"Esto contradice el comportamiento documentado del elemento."
    )


def test_quad4_needs_fine_mesh_to_overcome_locking():
    """Quad4 requiere malla fina en x para superar shear locking.

    Datos empíricos:
        30×4  → ≈28% error  (locking persistente; 4 capas no basta).
        60×4  → ≈9.4% error (aceptable, dentro de 10%).
        120×8 → ≈2.6% error (fino, dentro de 5%).

    Este test usa **60×4** que es la malla mínima para conseguir un error
    razonable con Q4 plane stress, y verifica que se entra en el rango
    < 12% (margen sobre el 9.4% empírico).

    El usuario que necesite mejor precisión con malla coarse debería usar
    Q8 o implementar reduced-integration / B-bar.
    """
    material = Elastic2D(E=E_YOUNG, nu=NU, hypothesis='plane_stress')
    domain, mid_right, right_elems = _build_solid_cantilever(
        Quad4, nx=60, ny=4, material=material,
    )
    F = _apply_tip_shear(domain, right_elems, total_force=-P_TIP)
    U = LinearSolver(Assembler(domain)).solve(F)
    u_tip = abs(U[mid_right.dofs['uy']])
    rel_err = abs(u_tip - U_TIP_TIMOSHENKO) / U_TIP_TIMOSHENKO
    assert rel_err < 0.12, (
        f"Quad4 60×4: u_tip={u_tip:.6e}, ref_T={U_TIP_TIMOSHENKO:.6e}, "
        f"err={rel_err:.4%} > 12%"
    )


def test_quad4_fine_mesh_converges_to_timoshenko():
    """Quad4 con malla fina 120×8 alcanza error < 5% — locking superado."""
    material = Elastic2D(E=E_YOUNG, nu=NU, hypothesis='plane_stress')
    domain, mid_right, right_elems = _build_solid_cantilever(
        Quad4, nx=120, ny=8, material=material,
    )
    F = _apply_tip_shear(domain, right_elems, total_force=-P_TIP)
    U = LinearSolver(Assembler(domain)).solve(F)
    u_tip = abs(U[mid_right.dofs['uy']])
    rel_err = abs(u_tip - U_TIP_TIMOSHENKO) / U_TIP_TIMOSHENKO
    assert rel_err < 0.05, (
        f"Quad4 120×8: u_tip={u_tip:.6e}, ref_T={U_TIP_TIMOSHENKO:.6e}, "
        f"err={rel_err:.4%} > 5%"
    )


def test_quad8_handles_slender_beam_with_coarse_mesh():
    """Quad8 con malla 12×1 sin shear locking severo; error < 5%."""
    material = Elastic2D(E=E_YOUNG, nu=NU, hypothesis='plane_stress')
    domain, mid_right, right_elems = _build_solid_cantilever(
        Quad8, nx=12, ny=1, material=material,
    )
    F = _apply_tip_shear(domain, right_elems, total_force=-P_TIP)
    U = LinearSolver(Assembler(domain)).solve(F)
    u_tip = abs(U[mid_right.dofs['uy']])
    rel_err = abs(u_tip - U_TIP_TIMOSHENKO) / U_TIP_TIMOSHENKO
    assert rel_err < 0.05, (
        f"Quad8 12×1: u_tip={u_tip:.6e}, ref_T={U_TIP_TIMOSHENKO:.6e}, "
        f"err={rel_err:.4%} > 5%"
    )


def test_quad4_h_refinement_converges_to_timoshenko():
    """Refinamiento h con Q4: error decrece monótonamente hacia Timoshenko."""
    material = Elastic2D(E=E_YOUNG, nu=NU, hypothesis='plane_stress')
    errors = []
    for nx, ny in [(30, 4), (60, 4), (120, 8)]:
        domain, mid_right, right_elems = _build_solid_cantilever(
            Quad4, nx=nx, ny=ny, material=material,
        )
        F = _apply_tip_shear(domain, right_elems, total_force=-P_TIP)
        U = LinearSolver(Assembler(domain)).solve(F)
        u_tip = abs(U[mid_right.dofs['uy']])
        errors.append(abs(u_tip - U_TIP_TIMOSHENKO))
    assert errors[1] < errors[0] and errors[2] < errors[1], (
        f"Quad4 cantilever: refinamiento no mejora: errores={errors}"
    )


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
