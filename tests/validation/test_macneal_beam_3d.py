"""MacNeal-Harder slender cantilever — extensión 3D con Hex8.

Referencia
----------
MacNeal, R.H.; Harder, R.L. (1985). *A proposed standard set of problems
to test finite element accuracy*. Finite Elements in Analysis and Design,
1(1), 3-20.

Problema
--------
Cantilever rectangular de longitud ``L = 6``, peralte ``h = 0.2`` (eje y),
espesor ``t = 0.1`` (eje z). Esbeltez ``L/h = 30``. Material elástico
isótropo ``E = 1×10⁷``, ``ν = 0.3``. Empotramiento total en ``x = 0``
(``u_x = u_y = u_z = 0`` en los 4 nodos de la cara). Carga transversal
``P = 1`` aplicada como tracción uniforme ``(0, P/A_cara, 0)`` en la cara
``x = L``, con ``A_cara = h · t``.

Solución analítica de Euler-Bernoulli (``L/h = 30`` ⇒ aporte de cortante
< 0.1%):

.. math::

    u_y(\\text{tip}) = \\frac{P \\cdot L^3}{3 \\cdot E \\cdot I},
    \\quad I = \\frac{t \\cdot h^3}{12}

Para los valores del problema: ``I = 0.1 · 0.008 / 12 = 6.6̄ × 10⁻⁵``,
y ``u_tip_EB = 1 · 216 / (3 · 1e7 · 6.6̄e-5) = 0.108``.

Objetivos del test
------------------
1. **Hex8 1×1 en la sección (12 elementos a lo largo)**: subestima
   severamente la flecha por shear locking — el espejo natural del
   Quad4 1 capa en 2D. Documenta la limitación.
2. **Convergencia con refinamiento en la sección (2×2 y a lo largo más
   denso)**: la flecha crece monótonamente hacia el valor analítico al
   reducir el shear locking. La propiedad cualitativa importante es la
   **monotonía**, no alcanzar el valor exacto con malla coarse — eso
   pediría Hex20.
"""
import os
import sys

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from solidum.core.domain import Domain
from solidum.elements.solid_3d import Hex8, Hex20, Hex27
from solidum.materials.elastic_3d import Elastic3D
from solidum.math.assembly import Assembler
from solidum.math.solvers import LinearSolver


L_BEAM = 6.0
H_BEAM = 0.2
T_BEAM = 0.1
E_YOUNG = 1.0e7
NU = 0.3
P_LOAD = 1.0


def _u_tip_euler_bernoulli():
    I = T_BEAM * H_BEAM ** 3 / 12.0
    return P_LOAD * L_BEAM ** 3 / (3.0 * E_YOUNG * I)


def _build_hex_beam(nx: int, ny: int, nz: int):
    """Cantilever de Hex8 con nx × ny × nz elementos."""
    domain = Domain()
    material = Elastic3D(E=E_YOUNG, nu=NU)

    nIx, nIy, nIz = nx + 1, ny + 1, nz + 1
    node_id = 0
    nid = {}
    for k in range(nIz):
        for j in range(nIy):
            for i in range(nIx):
                node_id += 1
                x = L_BEAM * i / nx
                y = -0.5 * H_BEAM + H_BEAM * j / ny
                z = -0.5 * T_BEAM + T_BEAM * k / nz
                nid[(i, j, k)] = domain.add_node(node_id, [x, y, z])

    eid = 0
    elems = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                # VTK_HEXAHEDRON order: bottom face (z=k) antihorario visto desde dentro;
                # top face (z=k+1) antihorario visto desde fuera.
                local = [
                    nid[(i,   j,   k)],
                    nid[(i+1, j,   k)],
                    nid[(i+1, j+1, k)],
                    nid[(i,   j+1, k)],
                    nid[(i,   j,   k+1)],
                    nid[(i+1, j,   k+1)],
                    nid[(i+1, j+1, k+1)],
                    nid[(i,   j+1, k+1)],
                ]
                eid += 1
                e = Hex8(eid, local, material)
                domain.add_element(e)
                elems.append((e, i, j, k))

    # BCs: empotramiento total en x=0.
    for n in domain.nodes.values():
        if abs(n.coordinates[0]) < 1e-12:
            n.fix_dof('ux', 0.0)
            n.fix_dof('uy', 0.0)
            n.fix_dof('uz', 0.0)

    return domain, elems, nid, nx, ny, nz


def _apply_tip_load(domain, elems, F, nx):
    """Tracción uniforme (0, P/A, 0) sobre la cara x=L de los elementos del extremo."""
    A_face = H_BEAM * T_BEAM
    t_vec = np.array([0.0, P_LOAD / A_face, 0.0])
    for e, i, j, k in elems:
        if i == nx - 1:
            # Cara +ξ del Hex8 = face index 3.
            f_elem = e.compute_face_traction(3, t_vec)
            for k_local, k_global in enumerate(e.get_global_dof_indices()):
                if k_global >= 0:
                    F[k_global] += f_elem[k_local]


def _solve_and_get_tip_uy(domain, elems, nx, ny, nz, nid):
    domain.generate_equation_numbers(verbose=False)
    assembler = Assembler(domain)
    F = np.zeros(domain.total_dofs)
    _apply_tip_load(domain, elems, F, nx)
    U = LinearSolver(assembler).solve(F)

    # Tip uy = promedio de uy en los 4 nodos centrados (axis (i=nx, j∈{ny/2}, k∈{nz/2}))
    # Para simplificar promediamos los 4 nodos de la cara x=L del centro de sección.
    if ny % 2 == 0 and nz % 2 == 0:
        j_mid = [ny // 2]
        k_mid = [nz // 2]
    else:
        j_mid = [ny // 2, ny // 2 + 1] if ny > 1 else [ny // 2]
        k_mid = [nz // 2, nz // 2 + 1] if nz > 1 else [nz // 2]
    samples = []
    for jm in j_mid:
        for km in k_mid:
            if (nx, jm, km) in nid:
                idx = nid[(nx, jm, km)].dofs['uy']
                if idx >= 0:
                    samples.append(U[idx])
    return float(np.mean(samples))


def test_macneal_3d_locking_coarse():
    """Hex8 12×1×1 documenta shear locking severo: u_tip < 40% del valor EB."""
    domain, elems, nid, nx, ny, nz = _build_hex_beam(nx=12, ny=1, nz=1)
    u_tip = _solve_and_get_tip_uy(domain, elems, nx, ny, nz, nid)
    u_analytical = _u_tip_euler_bernoulli()
    ratio = u_tip / u_analytical
    # Hex8 con 1 capa en sección y 12 a lo largo lockea fuerte. El valor
    # esperado es del orden 0.2-0.4 (similar al Quad4 1 capa en 2D).
    assert 0.10 < ratio < 0.55, (
        f"Hex8 12x1x1 lockea: ratio u_tip / u_EB = {ratio:.3f} fuera de "
        f"[0.10, 0.55]; u_tip={u_tip:.4e}, u_EB={u_analytical:.4e}."
    )


def test_macneal_3d_h_refinement_reduces_locking():
    """Refinar la malla reduce el locking monotonamente."""
    cases = [(12, 1, 1), (24, 2, 2), (48, 4, 4)]
    ratios = []
    for nx, ny, nz in cases:
        domain, elems, nid, nx_, ny_, nz_ = _build_hex_beam(nx, ny, nz)
        u_tip = _solve_and_get_tip_uy(domain, elems, nx_, ny_, nz_, nid)
        ratios.append(u_tip / _u_tip_euler_bernoulli())
    assert ratios[0] < ratios[1] < ratios[2], (
        f"Refinar no reduce locking monotonamente: ratios={ratios}"
    )
    # Con la malla más fina ya debemos estar razonablemente cerca de EB.
    assert ratios[-1] > 0.85, (
        f"Hex8 48×4×4 sigue lockeando demasiado: ratio={ratios[-1]:.3f}; "
        f"se esperaba >0.85 para tener confianza en la formulación."
    )


# =============================================================================
# Hex20 — convergencia del cuadrático sobre la misma viga.
# =============================================================================
#
# El cuadrático serendípito mitiga drásticamente el shear locking del Hex8 en
# este benchmark. Comparativa cuantitativa típica:
#   - Hex8  12×1×1 → ratio < 0.55  (lockea severo)
#   - Hex20 6×1×1  → ratio ≈ 0.97  (con la mitad de elementos, sin locking)
#   - Hex20 12×1×1 → ratio ≈ 0.99
# Esta es la motivación práctica de la sub-etapa A.ter: dominios con flexión
# o geometría curva donde el Hex8 obliga a mallas excesivamente densas.

def _build_hex20_beam(nx: int, ny: int, nz: int):
    """Cantilever Hex20 con nx × ny × nz elementos en una malla estructurada.

    Cada Hex20 tiene 20 nodos: 8 vértices + 12 medios de arista. La malla se
    construye usando un índice ``(2i, 2j, 2k)`` para vértices y ``(2i±1, ...)``
    para medios — sólo allocando posiciones con como mucho un índice impar
    (vértice + arista, descartando centros de cara e interior).
    """
    domain = Domain()
    material = Elastic3D(E=E_YOUNG, nu=NU)

    nIx = 2 * nx + 1
    nIy = 2 * ny + 1
    nIz = 2 * nz + 1
    nid = {}
    node_id = 0
    for k in range(nIz):
        for j in range(nIy):
            for i in range(nIx):
                n_odd = (i % 2) + (j % 2) + (k % 2)
                if n_odd <= 1:  # vértice (0) o medio de arista (1)
                    node_id += 1
                    x = L_BEAM * i / (2 * nx)
                    y = -0.5 * H_BEAM + H_BEAM * j / (2 * ny)
                    z = -0.5 * T_BEAM + T_BEAM * k / (2 * nz)
                    nid[(i, j, k)] = domain.add_node(node_id, [x, y, z])

    elems = []
    eid = 0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                ic, jc, kc = 2 * i, 2 * j, 2 * k
                vertices = [
                    nid[(ic,     jc,     kc)],
                    nid[(ic + 2, jc,     kc)],
                    nid[(ic + 2, jc + 2, kc)],
                    nid[(ic,     jc + 2, kc)],
                    nid[(ic,     jc,     kc + 2)],
                    nid[(ic + 2, jc,     kc + 2)],
                    nid[(ic + 2, jc + 2, kc + 2)],
                    nid[(ic,     jc + 2, kc + 2)],
                ]
                midsides = [
                    nid[(ic + 1, jc,     kc)],       # 8: 0-1
                    nid[(ic + 2, jc + 1, kc)],       # 9: 1-2
                    nid[(ic + 1, jc + 2, kc)],       # 10: 2-3
                    nid[(ic,     jc + 1, kc)],       # 11: 3-0
                    nid[(ic + 1, jc,     kc + 2)],   # 12: 4-5
                    nid[(ic + 2, jc + 1, kc + 2)],   # 13: 5-6
                    nid[(ic + 1, jc + 2, kc + 2)],   # 14: 6-7
                    nid[(ic,     jc + 1, kc + 2)],   # 15: 7-4
                    nid[(ic,     jc,     kc + 1)],   # 16: 0-4
                    nid[(ic + 2, jc,     kc + 1)],   # 17: 1-5
                    nid[(ic + 2, jc + 2, kc + 1)],   # 18: 2-6
                    nid[(ic,     jc + 2, kc + 1)],   # 19: 3-7
                ]
                eid += 1
                e = Hex20(eid, vertices + midsides, material)
                domain.add_element(e)
                elems.append((e, i, j, k))

    for n in domain.nodes.values():
        if abs(n.coordinates[0]) < 1e-12:
            n.fix_dof('ux', 0.0)
            n.fix_dof('uy', 0.0)
            n.fix_dof('uz', 0.0)

    return domain, elems, nid, nx, ny, nz


def _solve_hex20_and_get_tip_uy(domain, elems, nx, ny, nz, nid):
    domain.generate_equation_numbers(verbose=False)
    assembler = Assembler(domain)
    F = np.zeros(domain.total_dofs)
    A_face = H_BEAM * T_BEAM
    t_vec = np.array([0.0, P_LOAD / A_face, 0.0])
    for e, i, j, k in elems:
        if i == nx - 1:
            f_elem = e.compute_face_traction(3, t_vec)
            for k_local, k_global in enumerate(e.get_global_dof_indices()):
                if k_global >= 0:
                    F[k_global] += f_elem[k_local]
    U = LinearSolver(assembler).solve(F)
    # Tip uy en el nodo central de la cara x=L (vértice (2*nx, ny, nz)
    # si ny,nz son pares; en otro caso promedio sobre los 4 vértices).
    if (2 * nx, ny, nz) in nid:
        return float(U[nid[(2 * nx, ny, nz)].dofs['uy']])
    samples = []
    for jm in (0, 2 * ny):
        for km in (0, 2 * nz):
            if (2 * nx, jm, km) in nid:
                samples.append(U[nid[(2 * nx, jm, km)].dofs['uy']])
    return float(np.mean(samples))


def test_macneal_3d_hex20_coarse_already_accurate():
    """Hex20 6×1×1 con 1 capa en sección alcanza > 90% de la flecha EB.

    Comparativa cuantitativa de la mitigación del locking respecto al Hex8 con
    1 capa en sección — la motivación práctica de la sub-etapa A.ter.
    """
    domain, elems, nid, nx, ny, nz = _build_hex20_beam(nx=6, ny=1, nz=1)
    u_tip = _solve_hex20_and_get_tip_uy(domain, elems, nx, ny, nz, nid)
    ratio = u_tip / _u_tip_euler_bernoulli()
    assert ratio > 0.90, (
        f"Hex20 6×1×1 NO alcanza > 90% de u_EB: ratio = {ratio:.4f}; "
        f"se esperaba > 0.90. u_tip={u_tip:.4e}, u_EB={_u_tip_euler_bernoulli():.4e}."
    )
    assert ratio < 1.05, (
        f"Hex20 6×1×1 sobrepasa significativamente u_EB: ratio = {ratio:.4f}; "
        f"se esperaba < 1.05."
    )


def test_macneal_3d_hex20_dominates_hex8():
    """Hex20 6×1×1 alcanza > 90% de u_EB con la mitad de elementos que el
    Hex8 12×1×1 (que se queda en < 55%). Ratio cuantitativo Hex20/Hex8 > 1.7.
    """
    # Hex20 6×1×1
    domH20, elH20, nidH20, _, _, _ = _build_hex20_beam(nx=6, ny=1, nz=1)
    u_tip_h20 = _solve_hex20_and_get_tip_uy(domH20, elH20, 6, 1, 1, nidH20)
    # Hex8 12×1×1 (mismo número de nodos a lo largo, doble de elementos)
    domH8, elH8, nidH8, _, _, _ = _build_hex_beam(nx=12, ny=1, nz=1)
    u_tip_h8 = _solve_and_get_tip_uy(domH8, elH8, 12, 1, 1, nidH8)
    ratio = u_tip_h20 / u_tip_h8
    assert ratio > 1.7, (
        f"Hex20 6×1×1 no domina cuantitativamente a Hex8 12×1×1: "
        f"u_h20={u_tip_h20:.4e}, u_h8={u_tip_h8:.4e}, ratio={ratio:.3f} "
        f"(se esperaba > 1.7)."
    )


def test_macneal_3d_hex20_h_refinement_monotonic():
    """Refinar la malla Hex20 acerca monotonamente al valor analítico EB."""
    cases = [(6, 1, 1), (12, 1, 1), (12, 2, 2)]
    ratios = []
    for nx, ny, nz in cases:
        domain, elems, nid, _, _, _ = _build_hex20_beam(nx, ny, nz)
        u_tip = _solve_hex20_and_get_tip_uy(domain, elems, nx, ny, nz, nid)
        ratios.append(u_tip / _u_tip_euler_bernoulli())
    assert ratios[0] <= ratios[1] <= ratios[2], (
        f"Refinar no acerca monotonamente al valor analítico: ratios={ratios}"
    )
    assert ratios[-1] > 0.98, (
        f"Hex20 12×2×2 debería estar a > 98% de u_EB: ratio={ratios[-1]:.4f}."
    )


# =============================================================================
# Hex27 — espejo del Hex20 sobre el Lagrangiano completo.
# =============================================================================

def _build_hex27_beam(nx: int, ny: int, nz: int):
    """Cantilever Hex27. Cada elemento tiene 27 nodos: 20 del Hex20 + 6 centros
    de cara + 1 centro del cuerpo. La malla allocá **todas** las posiciones
    de la rejilla (2*nx+1)×(2*ny+1)×(2*nz+1) sin filtrar por ``n_odd``.
    """
    domain = Domain()
    material = Elastic3D(E=E_YOUNG, nu=NU)

    nIx = 2 * nx + 1
    nIy = 2 * ny + 1
    nIz = 2 * nz + 1
    nid = {}
    node_id = 0
    for k in range(nIz):
        for j in range(nIy):
            for i in range(nIx):
                node_id += 1
                x = L_BEAM * i / (2 * nx)
                y = -0.5 * H_BEAM + H_BEAM * j / (2 * ny)
                z = -0.5 * T_BEAM + T_BEAM * k / (2 * nz)
                nid[(i, j, k)] = domain.add_node(node_id, [x, y, z])

    elems = []
    eid = 0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                ic, jc, kc = 2 * i, 2 * j, 2 * k
                vertices = [
                    nid[(ic,     jc,     kc)],
                    nid[(ic + 2, jc,     kc)],
                    nid[(ic + 2, jc + 2, kc)],
                    nid[(ic,     jc + 2, kc)],
                    nid[(ic,     jc,     kc + 2)],
                    nid[(ic + 2, jc,     kc + 2)],
                    nid[(ic + 2, jc + 2, kc + 2)],
                    nid[(ic,     jc + 2, kc + 2)],
                ]
                midsides = [
                    nid[(ic + 1, jc,     kc)],       # 8: 0-1
                    nid[(ic + 2, jc + 1, kc)],       # 9: 1-2
                    nid[(ic + 1, jc + 2, kc)],       # 10: 2-3
                    nid[(ic,     jc + 1, kc)],       # 11: 3-0
                    nid[(ic + 1, jc,     kc + 2)],   # 12: 4-5
                    nid[(ic + 2, jc + 1, kc + 2)],   # 13: 5-6
                    nid[(ic + 1, jc + 2, kc + 2)],   # 14: 6-7
                    nid[(ic,     jc + 1, kc + 2)],   # 15: 7-4
                    nid[(ic,     jc,     kc + 1)],   # 16: 0-4
                    nid[(ic + 2, jc,     kc + 1)],   # 17: 1-5
                    nid[(ic + 2, jc + 2, kc + 1)],   # 18: 2-6
                    nid[(ic,     jc + 2, kc + 1)],   # 19: 3-7
                ]
                face_centers = [
                    nid[(ic,     jc + 1, kc + 1)],   # 20: face -x
                    nid[(ic + 2, jc + 1, kc + 1)],   # 21: face +x
                    nid[(ic + 1, jc,     kc + 1)],   # 22: face -y
                    nid[(ic + 1, jc + 2, kc + 1)],   # 23: face +y
                    nid[(ic + 1, jc + 1, kc)],       # 24: face -z
                    nid[(ic + 1, jc + 1, kc + 2)],   # 25: face +z
                ]
                body_center = [nid[(ic + 1, jc + 1, kc + 1)]]  # 26
                all_nodes = vertices + midsides + face_centers + body_center
                eid += 1
                e = Hex27(eid, all_nodes, material)
                domain.add_element(e)
                elems.append((e, i, j, k))

    for n in domain.nodes.values():
        if abs(n.coordinates[0]) < 1e-12:
            n.fix_dof('ux', 0.0)
            n.fix_dof('uy', 0.0)
            n.fix_dof('uz', 0.0)

    return domain, elems, nid, nx, ny, nz


def _solve_hex27_and_get_tip_uy(domain, elems, nx, ny, nz, nid):
    domain.generate_equation_numbers(verbose=False)
    assembler = Assembler(domain)
    F = np.zeros(domain.total_dofs)
    A_face = H_BEAM * T_BEAM
    t_vec = np.array([0.0, P_LOAD / A_face, 0.0])
    for e, i, j, k in elems:
        if i == nx - 1:
            f_elem = e.compute_face_traction(3, t_vec)
            for k_local, k_global in enumerate(e.get_global_dof_indices()):
                if k_global >= 0:
                    F[k_global] += f_elem[k_local]
    U = LinearSolver(assembler).solve(F)
    if (2 * nx, ny, nz) in nid:
        return float(U[nid[(2 * nx, ny, nz)].dofs['uy']])
    samples = []
    for jm in (0, 2 * ny):
        for km in (0, 2 * nz):
            if (2 * nx, jm, km) in nid:
                samples.append(U[nid[(2 * nx, jm, km)].dofs['uy']])
    return float(np.mean(samples))


def test_macneal_3d_hex27_coarse_already_accurate():
    """Hex27 6×1×1 alcanza > 90% de u_EB (mismo orden que el Hex20)."""
    domain, elems, nid, nx, ny, nz = _build_hex27_beam(nx=6, ny=1, nz=1)
    u_tip = _solve_hex27_and_get_tip_uy(domain, elems, nx, ny, nz, nid)
    ratio = u_tip / _u_tip_euler_bernoulli()
    assert ratio > 0.90, (
        f"Hex27 6×1×1 NO alcanza > 90% de u_EB: ratio={ratio:.4f}; "
        f"u_tip={u_tip:.4e}, u_EB={_u_tip_euler_bernoulli():.4e}."
    )
    assert ratio < 1.05, (
        f"Hex27 6×1×1 sobrepasa u_EB: ratio={ratio:.4f}."
    )


def test_macneal_3d_hex27_comparable_to_hex20():
    """Hex27 y Hex20 dan resultados muy similares en flexión simple — el
    espacio extra triquadrático del Hex27 rara vez gobierna en problemas
    donde la flexión Bernoulli domina (Cook §6.6).
    """
    domH27, elH27, nidH27, _, _, _ = _build_hex27_beam(nx=6, ny=1, nz=1)
    u_h27 = _solve_hex27_and_get_tip_uy(domH27, elH27, 6, 1, 1, nidH27)
    domH20, elH20, nidH20, _, _, _ = _build_hex20_beam(nx=6, ny=1, nz=1)
    u_h20 = _solve_hex20_and_get_tip_uy(domH20, elH20, 6, 1, 1, nidH20)
    # Ratio dentro de 3% — ambos elementos atrapan el mismo régimen flexional.
    rel_diff = abs(u_h27 - u_h20) / abs(u_h20)
    assert rel_diff < 0.03, (
        f"Hex27 y Hex20 difieren más del 3% en flexión simple: "
        f"u_h27={u_h27:.4e}, u_h20={u_h20:.4e}, rel_diff={rel_diff:.4f}."
    )
