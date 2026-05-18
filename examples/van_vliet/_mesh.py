"""Generador de malla estructurada para el benchmark Van Vliet (§8.1 Retama 2010).

Probeta dog-bone en tracción uniaxial con cabezales de acero, fiel a la
geometría de la fig 8.1 de la tesis y aclarada por el usuario. Unidades SI
(m, N, Pa).

Geometría (confirmada con el autor de la tesis, 2026-05-18)
-----------------------------------------------------------
- **Cuerpo (concreto)** dividido en tres tramos verticales:
  - Zona recta superior: 50 mm × 200 mm.
  - Cuello curvo central: 200 mm de altura, arco de circunferencia radio
    ``r = 145 mm``, ancho mínimo 120 mm en ``z = 0`` y ancho 200 mm en
    los extremos del arco (``z = ±100 mm``). Centro del arco en
    ``x_c = 0.060 + r = 0.205 m`` (off-body, a la derecha del lado derecho).
  - Zona recta inferior: 50 mm × 200 mm.
  - Altura total del cuerpo: 1.5 · D = 300 mm.
- **Cabezales (placas de acero)**: dos placas rectangulares, una encima y
  una debajo del cuerpo, **25 mm** de altura cada una. Ancho real de la
  placa = 250 mm (25 mm de alerón por cada lado del cuerpo), pero el
  alerón está **fuera de la línea de transmisión de carga** entre la
  rótula (x = +e, eje vertical) y el cuerpo, por lo que aquí se modela
  el cabezal con ancho 200 mm (igual al cuerpo en su zona recta).
  Documentado como simplificación en el ADR 0010.
- **Profundidad** (out-of-plane): 100 mm. Estado plane stress.

Condiciones de contorno (rotating boundary conditions, Van Vliet)
----------------------------------------------------------------
La carga entra por una rótula (hinge) en el cabezal superior y la
reacción por otra rótula en el cabezal inferior. Ambas rótulas
desplazadas lateralmente ``e = D/50 = 4 mm`` respecto al eje de simetría
vertical del cuerpo (eccentric pure tension, line of action en ``x = +e``).

Implementación MEF:
- ``load_node``: nodo único en el borde superior del cabezal superior, en
  ``x ≈ +e``. Se le aplica fuerza ``(0, +P)``.
- ``fix_node``: nodo único en el borde inferior del cabezal inferior, en
  ``x ≈ +e``. ``ux = uy = 0``.
- ``antirot_node``: con load_node y fix_node alineados en ``x = +e``, el
  sistema admite rotación rígida alrededor del eje vertical ``x = +e``
  (auto-valor cero de K_t). Sin momento driver no se excita, pero la
  matriz queda singular; se añade ``ux = 0`` en una esquina del cabezal
  inferior (``x ≈ -CABEZAL_WIDTH/2``, fuera del eje de carga) que rompe
  esa libertad sin introducir reacción física esperada.

Crack prescrita
---------------
Una fila central de elementos en ``z ≈ 0`` se construye con
``CST_Embedded2D`` y se pre-activa con ``n = (0, 1)``. El resto del
cuerpo y los cabezales son ``Tri3`` estándar.
"""
from __future__ import annotations

import numpy as np

from fenix.cohesive_materials.damage_isotropic import CohesiveDamageIsotropic
from fenix.core.domain import Domain
from fenix.elements.solid_2d.embedded_cst import CST_Embedded2D
from fenix.elements.solid_2d.tri3 import Tri3
from fenix.materials.elastic_2d import Elastic2D


# --- Parámetros físicos de §8.1 (Retama 2010) -----------------------------
D = 0.200                # m — dimensión característica
R = 0.145                # m — radio del arco del cuello
DEPTH = 0.100            # m — espesor (out-of-plane)

E_CONCRETE = 39.8e9      # Pa
NU_CONCRETE = 0.2
SIGMA_T0 = 2.57e6        # Pa
G_F = 121.9              # N/m  (= 0.1219 N/mm)

E_STEEL = 210.0e9        # Pa  (acero estructural estándar)
NU_STEEL = 0.3

E_ECC = D / 50.0         # m — excentricidad de la rótula (4 mm)
L_S = 0.6 * D            # m — distancia entre puntos de referencia (120 mm)

# --- Geometría de la probeta (concreto) -----------------------------------
BODY_STRAIGHT_H = 0.050  # m — altura de cada zona recta (50 mm)
BODY_CURVE_H = D         # m — altura del cuello curvo (200 mm)
BODY_TOTAL_H = BODY_CURVE_H + 2 * BODY_STRAIGHT_H  # 0.300 m = 1.5·D
BODY_MAX_WIDTH = D       # m — 200 mm (zonas rectas y extremos del cuello)
BODY_MIN_WIDTH = 0.6 * D # m — 120 mm (centro del cuello)
X_C_ARC = BODY_MIN_WIDTH / 2.0 + R  # m — centro del arco, off-body (0.205)

# --- Cabezales (acero) ----------------------------------------------------
CABEZAL_H = 0.025        # m — altura de cada placa (25 mm)
# Ancho real 250 mm (alerón de 25 mm por lado); aquí simplificado a 200 mm,
# justificación en el docstring del módulo.
CABEZAL_WIDTH = D        # m — 200 mm

# --- Penalty cohesivo -----------------------------------------------------
K_E_COHESIVE = 1.0e15    # N/m³


def half_width(z: float) -> float:
    """Mitad del ancho a altura ``z`` (z=0 en el centro del cuello).

    Cinco regiones (con z creciendo hacia arriba):
    - ``z ∈ [-Z_TOT, -Z_BODY]``: cabezal inferior, ancho ``CABEZAL_WIDTH``.
    - ``z ∈ [-Z_BODY, -Z_CURVE]``: zona recta inferior, ancho ``BODY_MAX_WIDTH``.
    - ``z ∈ [-Z_CURVE, +Z_CURVE]``: cuello curvo, ancho variable según el arco.
    - ``z ∈ [+Z_CURVE, +Z_BODY]``: zona recta superior, ancho ``BODY_MAX_WIDTH``.
    - ``z ∈ [+Z_BODY, +Z_TOT]``: cabezal superior, ancho ``CABEZAL_WIDTH``.
    """
    z_abs = abs(z)
    z_body = BODY_TOTAL_H / 2.0
    z_curve = BODY_CURVE_H / 2.0
    if z_abs > z_body:
        return CABEZAL_WIDTH / 2.0
    if z_abs > z_curve:
        return BODY_MAX_WIDTH / 2.0
    return X_C_ARC - np.sqrt(R**2 - z_abs**2)


def build_van_vliet_domain(
    n_x: int = 21,
    n_z_curve: int = 41,
    n_z_straight: int = 5,
    n_z_cabezal: int = 2,
    softening: str = 'linear',
):
    """Construye dominio completo y devuelve ``(domain, info)``.

    Parameters
    ----------
    n_x
        Nodos por fila horizontal (mapeo ``s ∈ [-1, +1]`` al medio-ancho local).
    n_z_curve
        Celdas verticales en la zona curva. **Debe ser impar** para colocar
        una fila central con centroide en ``z = 0``.
    n_z_straight
        Celdas verticales por zona recta del cuerpo (top, bottom).
    n_z_cabezal
        Celdas verticales por cabezal (top, bottom).
    softening
        ``'linear'`` o ``'exponential'``.

    Returns
    -------
    domain : Domain
    info : dict
        ``load_node``, ``fix_node``, ``antirot_node``, ``ref_top_node``,
        ``ref_bot_node``, ``crack_elements`` (pre-activados), ``cohesive``.
    """
    if n_z_curve % 2 != 1:
        raise ValueError(f"n_z_curve debe ser impar; recibido {n_z_curve}.")
    if softening not in ('linear', 'exponential'):
        raise ValueError(f"softening debe ser 'linear' o 'exponential'; recibido {softening!r}.")

    domain = Domain()

    bulk_conc = Elastic2D(E=E_CONCRETE, nu=NU_CONCRETE, hypothesis='plane_stress')
    bulk_steel = Elastic2D(E=E_STEEL, nu=NU_STEEL, hypothesis='plane_stress')
    cohesive = CohesiveDamageIsotropic(
        sigma_t0=SIGMA_T0, G_f=G_F, K_e=K_E_COHESIVE, softening=softening,
    )

    # === 1. Generar niveles de z (cinco regiones concatenadas) =============
    z_body = BODY_TOTAL_H / 2.0
    z_curve = BODY_CURVE_H / 2.0
    z_tot = z_body + CABEZAL_H

    z_cab_bot = np.linspace(-z_tot, -z_body, n_z_cabezal + 1)
    z_str_bot = np.linspace(-z_body, -z_curve, n_z_straight + 1)
    z_curve_levels = np.linspace(-z_curve, +z_curve, n_z_curve + 1)
    z_str_top = np.linspace(+z_curve, +z_body, n_z_straight + 1)
    z_cab_top = np.linspace(+z_body, +z_tot, n_z_cabezal + 1)

    z_all = np.concatenate([
        z_cab_bot, z_str_bot[1:], z_curve_levels[1:], z_str_top[1:], z_cab_top[1:],
    ])
    n_z_levels = len(z_all)  # (n_z_cabezal+1) + n_z_straight + n_z_curve + n_z_straight + n_z_cabezal

    # === 2. Generación de nodos (grid s × z con s mapeando a half_width(z)) =
    s_grid = np.linspace(-1.0, +1.0, n_x)
    node_id = 1
    nodes: dict[tuple[int, int], object] = {}

    for i_z, z in enumerate(z_all):
        w_half = half_width(z)
        for i_s, s in enumerate(s_grid):
            x = s * w_half
            nodes[(i_z, i_s)] = domain.add_node(node_id, [x, z])
            node_id += 1

    # === 3. Asignación celdas → región ====================================
    # Las celdas se indexan por i_z_cell ∈ [0, n_z_levels - 2].
    n_cells_cab_bot = n_z_cabezal
    n_cells_str_bot = n_z_straight
    n_cells_curve = n_z_curve
    n_cells_str_top = n_z_straight
    n_cells_total = n_z_levels - 1

    # Fila central del cuello (con centroide en z = 0): celda de rango
    # (n_cells_cab_bot + n_cells_str_bot + (n_cells_curve - 1) // 2).
    crack_row = n_cells_cab_bot + n_cells_str_bot + (n_cells_curve - 1) // 2

    def material_for_cell(i_z_cell: int):
        # Cabezal inferior
        if i_z_cell < n_cells_cab_bot:
            return bulk_steel
        # Cuerpo (recto inf + curva + recto sup)
        if i_z_cell < n_cells_cab_bot + n_cells_str_bot + n_cells_curve + n_cells_str_top:
            return bulk_conc
        # Cabezal superior
        return bulk_steel

    # === 4. Generación de elementos =======================================
    elem_id = 1
    crack_elements = []

    for i_z_cell in range(n_cells_total):
        is_crack_row = (i_z_cell == crack_row)
        material = material_for_cell(i_z_cell)
        for i_s_cell in range(n_x - 1):
            ll = nodes[(i_z_cell, i_s_cell)]
            lr = nodes[(i_z_cell, i_s_cell + 1)]
            ul = nodes[(i_z_cell + 1, i_s_cell)]
            ur = nodes[(i_z_cell + 1, i_s_cell + 1)]
            for tri_nodes in ([ll, lr, ur], [ll, ur, ul]):
                if is_crack_row:
                    elem = CST_Embedded2D(
                        elem_id, tri_nodes, material, cohesive, thickness=DEPTH,
                    )
                    crack_elements.append(elem)
                else:
                    elem = Tri3(elem_id, tri_nodes, material, thickness=DEPTH)
                domain.add_element(elem)
                elem_id += 1

    # === 5. Identificación de nodos especiales ============================
    def _closest_in_row(i_z: int, target_x: float):
        best = None
        best_d = np.inf
        for i_s in range(n_x):
            nd = nodes[(i_z, i_s)]
            d = abs(nd.coordinates[0] - target_x)
            if d < best_d:
                best_d = d
                best = nd
        return best

    load_node = _closest_in_row(n_z_levels - 1, +E_ECC)  # borde superior del cabezal sup.
    fix_node = _closest_in_row(0, +E_ECC)                # borde inferior del cabezal inf.
    antirot_node = _closest_in_row(0, -CABEZAL_WIDTH / 2.0)  # esquina inf-izq del cabezal inf.

    def _closest_centerline(target_z: float):
        best = None
        best_score = np.inf
        for (i_z, i_s), nd in nodes.items():
            score = abs(nd.coordinates[1] - target_z) + 1e-3 * abs(nd.coordinates[0])
            if score < best_score:
                best_score = score
                best = nd
        return best

    ref_top_node = _closest_centerline(+L_S / 2.0)
    ref_bot_node = _closest_centerline(-L_S / 2.0)

    # === 6. BCs y numeración ==============================================
    load_node.add_dof('ux'); load_node.add_dof('uy')
    fix_node.fix_dof('ux', 0.0)
    fix_node.fix_dof('uy', 0.0)
    antirot_node.fix_dof('ux', 0.0)

    domain.generate_equation_numbers(verbose=False)

    # === 7. Pre-activación de la fila de crack ============================
    n_dir = np.array([0.0, 1.0])
    for elem in crack_elements:
        coords = elem.get_coordinate_matrix(ndim=2)
        elem._activate(n_dir, coords)

    info = {
        'load_node': load_node,
        'fix_node': fix_node,
        'antirot_node': antirot_node,
        'ref_top_node': ref_top_node,
        'ref_bot_node': ref_bot_node,
        'crack_elements': crack_elements,
        'cohesive': cohesive,
        'BODY_TOTAL_H': BODY_TOTAL_H,
        'CABEZAL_WIDTH': CABEZAL_WIDTH,
    }
    return domain, info
