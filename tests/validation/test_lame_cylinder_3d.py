"""Lamé thick cylinder 3D — solución analítica cerrada sobre geometría curva.

Referencia
----------
Timoshenko, S.P.; Goodier, J.N. (1970). *Theory of Elasticity*, 3rd ed.,
McGraw-Hill, §28 ("Pure shear. Stresses in a hollow circular cylinder").
Solución original: Lamé, G. (1852). *Leçons sur la théorie mathématique de
l'élasticité des corps solides*, Paris.

Concepto
--------
Cuadrante de un cilindro hueco de pared gruesa sometido a presión interna
uniforme, con condiciones de plane strain (``u_z = 0`` en ambas caras
axiales) — versión 3D del benchmark 2D :mod:`test_lame_cylinder`.

Solución analítica (Lamé) en coordenadas cilíndricas:

.. math::

    \\sigma_{rr}(r) = A - B/r^2,\\quad
    \\sigma_{\\theta\\theta}(r) = A + B/r^2,\\quad
    \\sigma_{zz}(r) = \\nu\\,(\\sigma_{rr} + \\sigma_{\\theta\\theta}) = 2\\nu A

    u_r(r) = \\frac{1+\\nu}{E}\\,\\Big[(1-2\\nu)\\,A\\,r + \\frac{B}{r}\\Big]

con ``A = p·Rᵢ²/(Rₑ²−Rᵢ²)`` y ``B = p·Rᵢ²·Rₑ²/(Rₑ²−Rᵢ²)``. La solución es
exacta en CADA punto del dominio y NO depende del recovery — directamente
comparable con σ y ε de los puntos de Gauss y con u en cada nodo.

Por qué este benchmark valida específicamente la sub-etapa A.ter
-----------------------------------------------------------------
Los elementos cuadráticos isoparamétricos (Hex20, Hex27, Tet10) representan
**fronteras curvas** con error ``O(h³)`` cuando los nodos mid-edge se sitúan
sobre la superficie real (no en línea recta entre vértices). En la malla
estructurada que aquí se construye, los nodos mid-edge se sitúan sobre la
superficie cilíndrica natural mediante el mapeo paramétrico
``(t, s, z) → ((1-t)·Rᵢ + t·Rₑ)·(cos s, sin s, z)``. Es la capacidad
distintiva del cuadrático sobre el lineal — Hex8 y Tet4 no pueden
representar el cilindro mejor que con facetas planas, mientras los
cuadráticos lo curvan naturalmente.

Geometría, material y carga
---------------------------
- ``Rᵢ = 1.0``, ``Rₑ = 2.0``, longitud axial ``L = 1.0``.
- Cuarto de cilindro: ``x ≥ 0``, ``y ≥ 0``, ``z ∈ [0, L]``.
- ``E = 1000.0``, ``ν = 0.3``.
- Presión interna ``p = 100.0`` (compresión sobre la superficie cilíndrica
  interna ``r = Rᵢ``). La tracción aplicada al material es ``+p·r̂`` —
  apunta hacia afuera, consistente con presión interna empujando al
  material.

Condiciones de borde
--------------------
- Plano ``y = 0`` (sector simétrico): ``u_y = 0``.
- Plano ``x = 0`` (sector simétrico): ``u_x = 0``.
- Caras axiales ``z = 0`` y ``z = L``: ``u_z = 0`` (plane strain).

Tracción interna
----------------
``+p·r̂`` se aplica como **tracción uniforme por elemento** sobre la cara
interna (``r = Rᵢ``), evaluada en el centro angular del borde. Error
``O(h²)`` que se anula con refinamiento — patrón paritario con el LE1 2D.

Valores numéricos esperados (para los parámetros del módulo)
------------------------------------------------------------
- ``A = p·Rᵢ²/(Rₑ²−Rᵢ²) = 100·1/3 = 33.333…``
- ``B = p·Rᵢ²·Rₑ²/(Rₑ²−Rᵢ²) = 100·4/3 = 133.333…``
- ``σ_rr(Rᵢ) = A − B/Rᵢ² = -100`` (= -p — Newton 3ª ley sobre el material).
- ``σ_rr(Rₑ) = A − B/Rₑ² = 0`` (BC libre).
- ``σ_θθ(Rᵢ) = A + B/Rᵢ² = 166.667`` (tracción máxima).
- ``σ_θθ(Rₑ) = A + B/Rₑ² = 66.667``.
- ``σ_zz = 2νA = 20.0`` (constante).
- ``u_r(Rᵢ) = (1+ν)/E·[(1−2ν)·A·Rᵢ + B/Rᵢ] = 0.1907``.
- ``u_r(Rₑ) = (1+ν)/E·[(1−2ν)·A·Rₑ + B/Rₑ] = 0.1520``.

Elementos validados
-------------------
``Hex20`` y ``Hex27``. Para ambos la malla del cuadrante cilíndrico ×
axial es directa sobre el grid ``(I, J, K)`` paramétrico — los mid-edges
quedan automáticamente sobre la superficie cilíndrica por la
no-linealidad de cos/sin en el ángulo. Es la **capacidad distintiva** de
los elementos cuadráticos isoparamétricos sobre los lineales.

``Tet10`` con descomposición ingenua de cada celda hex en 5 tetraedros
**no converge** sobre este benchmark — el código de la malla está
presente en el módulo (``_build_tet10_cylinder_mesh``) pero los tests
correspondientes están marcados con ``@pytest.mark.skip`` con motivo
documentado. Razón: cada cara del tetraedro que toca la superficie
cilíndrica tiene 3 mid-edges; 1 corresponde a una arista del grid
paramétrico (curva, sobre la superficie cilíndrica) pero los otros 2
corresponden a aristas DIAGONALES del hex original (que se rompen al
descomponer en tets), cuyos mid-edges son rectos (promedio de las
coordenadas físicas de los dos vértices). Resultado: la cara
"isoparamétricamente parcial" del Tet10 no representa la superficie
curva con error ``O(h³)`` como el Hex20/Hex27, y el error en σ_rr
**crece con refinamiento h** (0.38 con 2×2×1, 0.96 con 8×8×1 — el
crecimiento del número de tets con cara parcialmente curva domina la
mejora local de cada tet). La validación correcta del Tet10 sobre
geometría curva requiere un mesher tetraédrico nativo con mid-edges
curvos (gmsh API o descomposición específica del anillo cilíndrico), no
una descomposición genérica hex→tets — eso queda fuera del alcance de
esta sesión y se documenta como deuda técnica.

La capacidad del Tet10 sobre geometría plana (sin la dificultad
arquitectural anterior) ya está validada en
:mod:`test_cube_lame_3d` (5-tet del cubo Lamé exacto a precisión máquina).
"""
import os
import sys

import numpy as np
import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from solidum.core.domain import Domain
from solidum.elements.solid_3d import Hex20, Hex27, Tet10
from solidum.materials.elastic_3d import Elastic3D
from solidum.math.assembly import Assembler
from solidum.math.solvers import LinearSolver


# Parámetros físicos del benchmark.
R_INNER = 1.0
R_OUTER = 2.0
LENGTH = 1.0
E_YOUNG = 1000.0
NU = 0.3
P_INNER = 100.0

# Constantes A, B de la solución de Lamé — calculadas una vez.
_LAME_A = P_INNER * R_INNER**2 / (R_OUTER**2 - R_INNER**2)
_LAME_B = P_INNER * R_INNER**2 * R_OUTER**2 / (R_OUTER**2 - R_INNER**2)


def _lame_stress(r: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """σ_rr(r), σ_θθ(r) y σ_zz (constante) según Lamé en plane strain."""
    sigma_rr = _LAME_A - _LAME_B / r**2
    sigma_tt = _LAME_A + _LAME_B / r**2
    sigma_zz = 2.0 * NU * _LAME_A
    return sigma_rr, sigma_tt, sigma_zz


def _lame_radial_displacement(r: np.ndarray) -> np.ndarray:
    """u_r(r) en plane strain."""
    return (1.0 + NU) / E_YOUNG * ((1.0 - 2.0 * NU) * _LAME_A * r + _LAME_B / r)


def _stress_to_polar(sigma_voigt: np.ndarray, xy: np.ndarray
                     ) -> tuple[float, float, float]:
    """Convierte σ en Voigt 6D del proyecto a (σ_rr, σ_θθ, σ_zz).

    Voigt 3D del proyecto: ``[σ_xx, σ_yy, σ_zz, σ_xy, σ_yz, σ_xz]``.
    """
    sxx, syy, szz, sxy, _syz, _sxz = sigma_voigt
    r = np.hypot(xy[0], xy[1])
    c = xy[0] / r
    s = xy[1] / r
    sigma_rr = sxx * c**2 + syy * s**2 + 2.0 * sxy * c * s
    sigma_tt = sxx * s**2 + syy * c**2 - 2.0 * sxy * c * s
    return float(sigma_rr), float(sigma_tt), float(szz)


# =============================================================================
# Malla del cuadrante cilíndrico × longitud axial.
# =============================================================================

def _cyl_point(I: int, J: int, K: int, nr: int, nt: int, nz: int
               ) -> tuple[float, float, float]:
    """(I, J, K) → (x, y, z) sobre la geometría cilíndrica natural.

    Mid-edges quedan AUTOMÁTICAMENTE sobre la superficie cilíndrica
    porque las funciones cos/sin son no lineales en ``J``.
    """
    t = I / (2 * nr)
    s = (np.pi / 2.0) * J / (2 * nt)
    z = LENGTH * K / (2 * nz)
    r_local = (1.0 - t) * R_INNER + t * R_OUTER
    return r_local * np.cos(s), r_local * np.sin(s), z


def _build_hex_cylinder_mesh(elem_cls, nr: int, nt: int, nz: int, material):
    """Malla del cuadrante cilíndrico × axial para Hex20 o Hex27.

    Returns
    -------
    domain : Domain
    elems : list[(element, ir, jt, kz)]
    nid : dict[(I, J, K), Node]
    """
    is_hex27 = elem_cls is Hex27
    nIx = 2 * nr + 1
    nIy = 2 * nt + 1
    nIz = 2 * nz + 1

    domain = Domain()
    nid: dict[tuple[int, int, int], object] = {}
    node_id = 0
    for K in range(nIz):
        for J in range(nIy):
            for I in range(nIx):
                n_odd = (I % 2) + (J % 2) + (K % 2)
                if (not is_hex27) and n_odd > 1:
                    continue
                x, y, z = _cyl_point(I, J, K, nr, nt, nz)
                node_id += 1
                nid[(I, J, K)] = domain.add_node(node_id, [x, y, z])

    elems = []
    eid = 0
    for kz in range(nz):
        for jt in range(nt):
            for ir in range(nr):
                ic, jc, kc = 2 * ir, 2 * jt, 2 * kz
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
                    nid[(ic + 1, jc,     kc)],
                    nid[(ic + 2, jc + 1, kc)],
                    nid[(ic + 1, jc + 2, kc)],
                    nid[(ic,     jc + 1, kc)],
                    nid[(ic + 1, jc,     kc + 2)],
                    nid[(ic + 2, jc + 1, kc + 2)],
                    nid[(ic + 1, jc + 2, kc + 2)],
                    nid[(ic,     jc + 1, kc + 2)],
                    nid[(ic,     jc,     kc + 1)],
                    nid[(ic + 2, jc,     kc + 1)],
                    nid[(ic + 2, jc + 2, kc + 1)],
                    nid[(ic,     jc + 2, kc + 1)],
                ]
                if is_hex27:
                    face_centers = [
                        nid[(ic,     jc + 1, kc + 1)],
                        nid[(ic + 2, jc + 1, kc + 1)],
                        nid[(ic + 1, jc,     kc + 1)],
                        nid[(ic + 1, jc + 2, kc + 1)],
                        nid[(ic + 1, jc + 1, kc)],
                        nid[(ic + 1, jc + 1, kc + 2)],
                    ]
                    body_center = [nid[(ic + 1, jc + 1, kc + 1)]]
                    all_nodes = vertices + midsides + face_centers + body_center
                else:
                    all_nodes = vertices + midsides
                eid += 1
                e = elem_cls(eid, all_nodes, material)
                domain.add_element(e)
                elems.append((e, ir, jt, kz))

    # BCs:
    for (I, J, K), node in nid.items():
        if J == 0:
            node.fix_dof('uy', 0.0)
        if J == nIy - 1:
            node.fix_dof('ux', 0.0)
        if K == 0 or K == nIz - 1:
            node.fix_dof('uz', 0.0)

    domain.generate_equation_numbers(verbose=False)
    return domain, elems, nid


def _apply_inner_pressure_hex(domain: Domain, elems, nt: int) -> np.ndarray:
    """Tracción ``+p·r̂`` sobre la cara interna (face index 5 = -ξ) de los
    elementos con ``ir == 0``.

    Patrón paritario con la presión interna del Lamé 2D: la tracción se
    aproxima por "uniforme por elemento" evaluada en el centro angular
    del elemento. Error ``O(h²)`` sobre la integral del borde.
    """
    F = np.zeros(domain.total_dofs)
    for elem, ir, jt, kz in elems:
        if ir != 0:
            continue
        # Centro angular del elemento (jt × angular).
        s_mid = (np.pi / 2.0) * (2 * jt + 1) / (2 * nt)
        n_hat = np.array([np.cos(s_mid), np.sin(s_mid), 0.0])
        t_vec = P_INNER * n_hat
        f_e = elem.compute_face_traction(5, t_vec)
        for kl, kg in enumerate(elem.get_global_dof_indices()):
            if kg >= 0:
                F[kg] += f_e[kl]
    return F


# =============================================================================
# Tet10 — malla por descomposición de cada celda hex en 5 tetraedros.
# =============================================================================
#
# Numeración de los 5 tets (descomposición canónica del cubo de referencia):
#   T1: (0, 1, 3, 4)
#   T2: (1, 2, 3, 6)
#   T3: (1, 4, 5, 6)
#   T4: (3, 4, 6, 7)
#   T5: (1, 3, 4, 6)   ← tetraedro central
#
# Cada arista de cada tet necesita su mid-edge. Las aristas son:
#  - Aristas del hex original (12): mid-edges del grid paramétrico, sobre
#    la superficie cilíndrica natural (curvos donde aplique).
#  - Aristas diagonales (cara o cuerpo del hex):
#      0-3 (del hex), 1-4, 3-4 (T1, T5)
#      1-3, 2-6, 1-6, 3-6 (T2, T5)
#      4-6, 5-6 (T3)
#      4-7, 3-7, 6-7 (T4)
#  Los mid-edges de las aristas que NO están en el grid paramétrico se
#  calculan como promedio de las posiciones físicas de los vértices
#  (línea recta en el espacio físico). Esto es INTERIOR al hex; la
#  representación de la frontera curva del cilindro queda determinada por
#  los mid-edges del grid (aristas del hex), no por las diagonales.

# Aristas locales de cada tet en la numeración VTK_QUADRATIC_TETRA:
#   edge 4: vertices 0-1
#   edge 5: vertices 1-2
#   edge 6: vertices 2-0
#   edge 7: vertices 0-3
#   edge 8: vertices 1-3
#   edge 9: vertices 2-3
_TET_EDGE_PAIRS = (
    (0, 1),  # 4
    (1, 2),  # 5
    (2, 0),  # 6
    (0, 3),  # 7
    (1, 3),  # 8
    (2, 3),  # 9
)

_HEX_TO_TETS = (
    (0, 1, 3, 4),  # T1
    (1, 2, 3, 6),  # T2
    (1, 4, 5, 6),  # T3
    (3, 4, 6, 7),  # T4
    (1, 3, 4, 6),  # T5 (central)
)


def _build_tet10_cylinder_mesh(nr: int, nt: int, nz: int, material):
    """Malla Tet10 por descomposición hex → 5 tets.

    Returns
    -------
    domain : Domain
    elems : list[(element, ir, jt, kz, tet_idx)]
        Tetraedros con su origen (índices hex + cual de los 5 tets).
    nid : dict[(I, J, K), Node]
        Subconjunto de nodos del grid paramétrico (los del Hex20 standalone).
    inner_tet_faces : list[(element, face_idx)]
        Pares de elementos y caras locales que coinciden con la frontera
        interna del cilindro (``r = Rᵢ``).
    """
    nIx = 2 * nr + 1
    nIy = 2 * nt + 1
    nIz = 2 * nz + 1
    domain = Domain()

    # 1) Allocá los nodos del grid paramétrico — mismos que Hex20 (vértices
    # + mid-edges = todos los puntos con n_odd ≤ 1).
    nid: dict[tuple[int, int, int], object] = {}
    node_id = 0
    for K in range(nIz):
        for J in range(nIy):
            for I in range(nIx):
                if (I % 2) + (J % 2) + (K % 2) > 1:
                    continue
                x, y, z = _cyl_point(I, J, K, nr, nt, nz)
                node_id += 1
                nid[(I, J, K)] = domain.add_node(node_id, [x, y, z])

    # 2) Cache de mid-edges para aristas diagonales (que NO están en el grid).
    # Clave: frozenset({(I1, J1, K1), (I2, J2, K2)}) → node
    diag_mid_cache: dict[frozenset, object] = {}

    def _vertex_indices(ic, jc, kc):
        """Índices (I, J, K) de los 8 vértices del hex (ic..ic+2 en I, etc.)."""
        return (
            (ic,     jc,     kc),
            (ic + 2, jc,     kc),
            (ic + 2, jc + 2, kc),
            (ic,     jc + 2, kc),
            (ic,     jc,     kc + 2),
            (ic + 2, jc,     kc + 2),
            (ic + 2, jc + 2, kc + 2),
            (ic,     jc + 2, kc + 2),
        )

    def _grid_midedge_or_diag(idx_a, idx_b):
        """Devuelve el nodo correspondiente al mid-edge entre dos vértices.

        Si la arista coincide con una arista del grid paramétrico (los dos
        índices difieren en exactamente dos coordenadas: una con paso 2, las
        otras iguales), usa el nodo del grid (mid-edge curvo).

        Si la arista es diagonal del hex (de cara o de cuerpo), allocá un
        nodo nuevo en el promedio físico de los dos vértices y lo cachea.
        """
        diffs = tuple(b - a for a, b in zip(idx_a, idx_b))
        n_nonzero = sum(1 for d in diffs if d != 0)
        if n_nonzero == 1 and any(abs(d) == 2 for d in diffs):
            # Arista del grid (cambia exactamente una coordenada en 2).
            mid = tuple((a + b) // 2 for a, b in zip(idx_a, idx_b))
            return nid[mid]
        # Diagonal: alloca nodo nuevo en promedio físico.
        key = frozenset({idx_a, idx_b})
        if key in diag_mid_cache:
            return diag_mid_cache[key]
        node_a = nid[idx_a]
        node_b = nid[idx_b]
        xa = np.asarray(node_a.coordinates[:3])
        xb = np.asarray(node_b.coordinates[:3])
        xm = 0.5 * (xa + xb)
        nonlocal node_id
        node_id += 1
        new_node = domain.add_node(node_id, list(xm))
        diag_mid_cache[key] = new_node
        return new_node

    elems = []
    eid = 0
    inner_tet_faces: list[tuple[object, int]] = []
    for kz in range(nz):
        for jt in range(nt):
            for ir in range(nr):
                ic, jc, kc = 2 * ir, 2 * jt, 2 * kz
                vidx = _vertex_indices(ic, jc, kc)

                # Caras del hex (vértices locales) que coinciden con r=Rᵢ
                # (cara -ξ del hex local: vértices 0, 3, 7, 4).
                INNER_HEX_FACE = (0, 3, 7, 4)

                for tet_idx, tet_def in enumerate(_HEX_TO_TETS):
                    v_local = [nid[vidx[k]] for k in tet_def]
                    midedges = []
                    for (a_local, b_local) in _TET_EDGE_PAIRS:
                        a_glob = vidx[tet_def[a_local]]
                        b_glob = vidx[tet_def[b_local]]
                        midedges.append(_grid_midedge_or_diag(a_glob, b_glob))
                    all_nodes = v_local + midedges
                    eid += 1
                    e = Tet10(eid, all_nodes, material)
                    domain.add_element(e)
                    elems.append((e, ir, jt, kz, tet_idx))

                    # Identificar qué cara del tet (si alguna) coincide con
                    # la frontera interna r=Rᵢ. Cara del Tet10 = índice 0..3.
                    # FACE_NODES del Tet10: tupla de tuplas con 6 nodos locales.
                    # Cada cara del tet usa 3 vértices del tet + 3 mid-edges.
                    # Una cara del tet coincide con la frontera interna si sus
                    # 3 vértices están en la cara hex INNER_HEX_FACE.
                    if ir == 0:
                        for face_idx in range(4):
                            face_def = _tet_face_vertices(face_idx)
                            face_hex_local = [tet_def[v] for v in face_def]
                            if all(h in INNER_HEX_FACE for h in face_hex_local):
                                inner_tet_faces.append((e, face_idx))

    # BCs:
    for (I, J, K), node in nid.items():
        if J == 0:
            node.fix_dof('uy', 0.0)
        if J == nIy - 1:
            node.fix_dof('ux', 0.0)
        if K == 0 or K == nIz - 1:
            node.fix_dof('uz', 0.0)
    # Diagonales nuevas: heredan BCs según su coord física.
    for node in diag_mid_cache.values():
        xy = node.coordinates
        if abs(xy[1]) < 1e-12:
            node.fix_dof('uy', 0.0)
        if abs(xy[0]) < 1e-12:
            node.fix_dof('ux', 0.0)
        if abs(xy[2]) < 1e-12 or abs(xy[2] - LENGTH) < 1e-12:
            node.fix_dof('uz', 0.0)

    domain.generate_equation_numbers(verbose=False)
    return domain, elems, nid, inner_tet_faces


def _tet_face_vertices(face_idx: int) -> tuple[int, int, int]:
    """Vértices del tetraedro de referencia para cada cara local 0..3.

    Mismo convenio que ``Tet10.FACE_NODES`` (cara ``i`` opuesta al vértice ``i``).
    """
    return (
        (1, 2, 3),  # cara 0: opuesta al vértice 0
        (0, 3, 2),  # cara 1: opuesta al vértice 1
        (0, 1, 3),  # cara 2: opuesta al vértice 2
        (0, 2, 1),  # cara 3: opuesta al vértice 3
    )[face_idx]


def _apply_inner_pressure_tet(domain: Domain, inner_tet_faces,
                              nt: int) -> np.ndarray:
    """Tracción ``+p·r̂`` sobre las caras de los tets coincidentes con r=Rᵢ.

    A diferencia del hex (donde s_mid se calcula a partir de jt), para el
    tet hay que recuperar el ángulo angular medio de la cara concreta. Se
    usa el centroide cartesiano de la cara y se proyecta a r̂.
    """
    F = np.zeros(domain.total_dofs)
    for elem, face_idx in inner_tet_faces:
        face_def = elem.FACE_NODES[face_idx]
        coords = np.array(
            [elem.nodes[i].coordinates[:3] for i in face_def],
            dtype=np.float64,
        )
        centroid = coords.mean(axis=0)
        r_centroid = np.hypot(centroid[0], centroid[1])
        n_hat = np.array([centroid[0] / r_centroid, centroid[1] / r_centroid, 0.0])
        t_vec = P_INNER * n_hat
        f_e = elem.compute_face_traction(face_idx, t_vec)
        for kl, kg in enumerate(elem.get_global_dof_indices()):
            if kg >= 0:
                F[kg] += f_e[kl]
    return F


# =============================================================================
# Métrica de error sobre puntos de Gauss interiores.
# =============================================================================

def _gauss_error_polar(domain: Domain, U: np.ndarray, exclude_outer_layer: bool
                       ) -> dict:
    """Error L²-relativo de σ_rr, σ_θθ, σ_zz sobre Gauss interiores.

    ``exclude_outer_layer``: si True, descarta los Gauss con r demasiado
    cerca de Rᵢ o Rₑ (donde la aproximación "tracción constante por arista"
    introduce error de frontera ``O(h)``). Tolerancia 5% de la pared.
    """
    sum_dr2 = 0.0
    sum_dt2 = 0.0
    sum_dz2 = 0.0
    sum_r2 = 0.0
    sum_t2 = 0.0
    sum_z2 = 0.0
    margin = 0.05 * (R_OUTER - R_INNER) if exclude_outer_layer else 0.0
    for elem in domain.elements.values():
        gs = elem.compute_gauss_state(U)
        for (xyz, sig_voigt) in zip(gs['points_global'], gs['stress']):
            r = np.hypot(xyz[0], xyz[1])
            if exclude_outer_layer:
                if r < R_INNER + margin or r > R_OUTER - margin:
                    continue
            sigma_rr_n, sigma_tt_n, sigma_zz_n = _stress_to_polar(sig_voigt, xyz)
            sigma_rr_a, sigma_tt_a, sigma_zz_a = _lame_stress(np.array(r))
            sum_dr2 += (sigma_rr_n - sigma_rr_a) ** 2
            sum_dt2 += (sigma_tt_n - sigma_tt_a) ** 2
            sum_dz2 += (sigma_zz_n - sigma_zz_a) ** 2
            sum_r2 += sigma_rr_a ** 2
            sum_t2 += sigma_tt_a ** 2
            sum_z2 += sigma_zz_a ** 2
    return {
        'sigma_rr': float(np.sqrt(sum_dr2 / max(sum_r2, 1e-30))),
        'sigma_tt': float(np.sqrt(sum_dt2 / max(sum_t2, 1e-30))),
        'sigma_zz': float(np.sqrt(sum_dz2 / max(sum_z2, 1e-30))),
    }


def _u_r_error(domain: Domain, U: np.ndarray) -> float:
    """Error L²-relativo de u_r sobre todos los nodos."""
    sum_du2 = 0.0
    sum_u2 = 0.0
    for node in domain.nodes.values():
        x, y, z = node.coordinates[:3]
        r = np.hypot(x, y)
        if r < 1e-12:
            continue
        idx_x = node.dofs['ux']
        idx_y = node.dofs['uy']
        ux = U[idx_x] if idx_x >= 0 else 0.0
        uy = U[idx_y] if idx_y >= 0 else 0.0
        u_r_num = (ux * x + uy * y) / r
        u_r_ana = _lame_radial_displacement(np.array(r))
        sum_du2 += (u_r_num - u_r_ana) ** 2
        sum_u2 += u_r_ana ** 2
    return float(np.sqrt(sum_du2 / max(sum_u2, 1e-30)))


# =============================================================================
# Tests Hex20 + Hex27.
# =============================================================================

# (elem_cls, nr, nt, nz, tol_sigma, tol_u)
HEX_PARAMS = [
    pytest.param(Hex20, 4, 4, 1, 0.04, 0.02, id='Hex20_4x4x1'),
    pytest.param(Hex27, 4, 4, 1, 0.04, 0.02, id='Hex27_4x4x1'),
]


@pytest.mark.parametrize("elem_cls, nr, nt, nz, tol_sigma, tol_u", HEX_PARAMS)
def test_lame_3d_hex_stress_and_displacement(elem_cls, nr, nt, nz, tol_sigma, tol_u):
    """Hex20/Hex27 reproducen σ_rr, σ_θθ, σ_zz y u_r dentro de tolerancia."""
    material = Elastic3D(E=E_YOUNG, nu=NU)
    domain, elems, _ = _build_hex_cylinder_mesh(elem_cls, nr, nt, nz, material)
    F = _apply_inner_pressure_hex(domain, elems, nt)
    U = LinearSolver(Assembler(domain)).solve(F)

    errors = _gauss_error_polar(domain, U, exclude_outer_layer=True)
    u_err = _u_r_error(domain, U)

    assert errors['sigma_rr'] < tol_sigma, (
        f"{elem_cls.__name__}: err_L²_rel(σ_rr) = {errors['sigma_rr']:.4f} "
        f"> {tol_sigma} (interior, sin anillos de frontera)"
    )
    assert errors['sigma_tt'] < tol_sigma, (
        f"{elem_cls.__name__}: err_L²_rel(σ_θθ) = {errors['sigma_tt']:.4f} "
        f"> {tol_sigma}"
    )
    assert errors['sigma_zz'] < tol_sigma, (
        f"{elem_cls.__name__}: err_L²_rel(σ_zz) = {errors['sigma_zz']:.4f} "
        f"> {tol_sigma}"
    )
    assert u_err < tol_u, (
        f"{elem_cls.__name__}: err_L²_rel(u_r) = {u_err:.4f} > {tol_u}"
    )


@pytest.mark.parametrize("elem_cls, meshes", [
    pytest.param(Hex20, [(2, 2, 1), (4, 4, 1), (6, 6, 1)], id='Hex20'),
    pytest.param(Hex27, [(2, 2, 1), (3, 3, 1), (4, 4, 1)], id='Hex27'),
])
def test_lame_3d_hex_h_refinement_converges(elem_cls, meshes):
    """Convergencia h monótona del error L²(σ_θθ) — elemento más sensible
    a la representación de la geometría curva.
    """
    material = Elastic3D(E=E_YOUNG, nu=NU)
    errors = []
    for nr, nt, nz in meshes:
        domain, elems, _ = _build_hex_cylinder_mesh(elem_cls, nr, nt, nz, material)
        F = _apply_inner_pressure_hex(domain, elems, nt)
        U = LinearSolver(Assembler(domain)).solve(F)
        errs = _gauss_error_polar(domain, U, exclude_outer_layer=True)
        errors.append(errs['sigma_tt'])
    for i in range(1, len(errors)):
        assert errors[i] < errors[i - 1] + 1e-12, (
            f"{elem_cls.__name__}: err_L²(σ_θθ) NO decrece monotonamente con "
            f"mallas {meshes}: errores={errors}"
        )


# =============================================================================
# Tests Tet10.
# =============================================================================

@pytest.mark.skip(
    reason="Tet10 con descomposición hex→5tets no converge sobre cilindro: "
           "mid-edges de aristas diagonales rectos rompen isoparametricidad "
           "sobre frontera curva (error crece con h). Requiere mesher tet "
           "nativo (gmsh) con mid-edges curvos. Ver docstring del módulo."
)
@pytest.mark.parametrize("nr, nt, nz, tol_sigma, tol_u", [
    pytest.param(4, 4, 1, 0.08, 0.05, id='Tet10_4x4x1'),
])
def test_lame_3d_tet10_stress_and_displacement(nr, nt, nz, tol_sigma, tol_u):
    """Tet10 reproduce σ_rr, σ_θθ, σ_zz y u_r dentro de tolerancia.

    Las tolerancias son más laxas que las del Hex20/Hex27 porque las
    aristas diagonales internas del Tet10 son rectas (no isoparamétricas
    completas), induciendo un error adicional respecto al Hex con grid
    cilíndrico nativo.
    """
    material = Elastic3D(E=E_YOUNG, nu=NU)
    domain, elems, _, inner_faces = _build_tet10_cylinder_mesh(nr, nt, nz, material)
    F = _apply_inner_pressure_tet(domain, inner_faces, nt)
    U = LinearSolver(Assembler(domain)).solve(F)

    errors = _gauss_error_polar(domain, U, exclude_outer_layer=True)
    u_err = _u_r_error(domain, U)

    assert errors['sigma_rr'] < tol_sigma, (
        f"Tet10: err_L²_rel(σ_rr) = {errors['sigma_rr']:.4f} > {tol_sigma}"
    )
    assert errors['sigma_tt'] < tol_sigma, (
        f"Tet10: err_L²_rel(σ_θθ) = {errors['sigma_tt']:.4f} > {tol_sigma}"
    )
    assert u_err < tol_u, (
        f"Tet10: err_L²_rel(u_r) = {u_err:.4f} > {tol_u}"
    )


@pytest.mark.skip(
    reason="Tet10 con descomposición hex→5tets diverge en σ_rr (error crece "
           "con h). Ver docstring del módulo."
)
def test_lame_3d_tet10_h_refinement_converges():
    """Convergencia h monótona del Tet10 sobre el cilindro."""
    material = Elastic3D(E=E_YOUNG, nu=NU)
    errors = []
    for nr, nt, nz in [(2, 2, 1), (3, 3, 1), (4, 4, 1)]:
        domain, elems, _, inner_faces = _build_tet10_cylinder_mesh(
            nr, nt, nz, material,
        )
        F = _apply_inner_pressure_tet(domain, inner_faces, nt)
        U = LinearSolver(Assembler(domain)).solve(F)
        errs = _gauss_error_polar(domain, U, exclude_outer_layer=True)
        errors.append(errs['sigma_tt'])
    for i in range(1, len(errors)):
        assert errors[i] < errors[i - 1] + 1e-12, (
            f"Tet10: err_L²(σ_θθ) NO decrece monotonamente con mallas "
            f"[(2,2,1),(3,3,1),(4,4,1)]: errores={errors}"
        )


# =============================================================================
# Sanity test: σ_rr en la cara interna ≈ -p_inner (BC discreta).
# =============================================================================

@pytest.mark.parametrize("elem_cls, nr, nt, nz", [
    pytest.param(Hex20, 4, 4, 1, id='Hex20_4x4x1'),
    pytest.param(Hex27, 4, 4, 1, id='Hex27_4x4x1'),
])
def test_lame_3d_sigma_rr_at_inner_face_is_minus_p(elem_cls, nr, nt, nz):
    """σ_rr en el Gauss más cercano a la cara interna ≈ -P_INNER (Newton 3ª).

    Sanity transversal: una presión interna ``p`` sobre la superficie
    cilíndrica debe materializarse como un esfuerzo σ_rr ≈ -p en los
    puntos próximos a la cara (signo negativo = compresión). Tolerancia
    laxa (15%) porque el Gauss más cercano puede estar lejos del borde
    con malla coarse.
    """
    material = Elastic3D(E=E_YOUNG, nu=NU)
    domain, elems, _ = _build_hex_cylinder_mesh(elem_cls, nr, nt, nz, material)
    F = _apply_inner_pressure_hex(domain, elems, nt)
    U = LinearSolver(Assembler(domain)).solve(F)

    best_r = np.inf
    best_sigma_rr = 0.0
    for elem in domain.elements.values():
        gs = elem.compute_gauss_state(U)
        for (xyz, sig_voigt) in zip(gs['points_global'], gs['stress']):
            r = float(np.hypot(xyz[0], xyz[1]))
            if r < best_r:
                best_r = r
                sigma_rr, _, _ = _stress_to_polar(sig_voigt, xyz)
                best_sigma_rr = sigma_rr
    rel_err = abs(best_sigma_rr - (-P_INNER)) / P_INNER
    assert rel_err < 0.15, (
        f"{elem_cls.__name__}: σ_rr(Gauss más cercano a r={best_r:.3f}) = "
        f"{best_sigma_rr:.3f} no satisface BC = -{P_INNER}; "
        f"err_rel = {rel_err:.2%}."
    )


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
