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

from fenix.core.domain import Domain
from fenix.elements.solid_3d import Hex8
from fenix.materials.elastic_3d import Elastic3D
from fenix.math.assembly import Assembler
from fenix.math.solvers import LinearSolver


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
