"""Utilidades compartidas por los elementos del paquete ``solidum.elements.solid_3d`` (ADR 0012).

Cinco grupos:

1. **Kinematics Numba** para los elementos lineales (Hex8, Tet4):
   :func:`_compute_kinematics_hex8`, :func:`_compute_kinematics_tet4`,
   :func:`_compute_integrands_3d`, :func:`_shape_functions_hex8`,
   :func:`_det_jacobian_hex8`, :func:`_shape_functions_tet4`.
2. **Funciones de forma y derivadas en numpy puro** para los elementos
   de orden superior 3D (sub-etapa A.ter): :func:`_N_hex20`,
   :func:`_dN_hex20`, :func:`_N_hex27`, :func:`_dN_hex27` (y, en sub-fase
   sucesiva, sus equivalentes para Tet10).
3. **Kinematics genérico de orden superior 3D**:
   :func:`_kinematics_higher_order_3d` — calcula ``B`` (6×(3·n_nodes)) y
   ``detJ`` a partir de un par ``(shape_fn, grad_fn)`` arbitrario y la
   matriz de coordenadas globales.
4. **Expansor escalar→bloque 3D** para masa traslacional:
   :func:`_expand_scalar_mass_3d`.
5. **Base interna** :class:`_HigherOrderSolid3D` que comparten los
   elementos sólidos 3D isoparamétricos de orden superior: Hex20 y Hex27
   (y Tet10 en sub-fase 3). Comparte los bucles de Gauss para ``K``,
   ``F_int``, ``gauss_state``, ``body_load``, ``face_traction`` y masa
   consistente. Cada subclase declara las funciones de forma, la
   cuadratura por defecto, la cuadratura específica de masa (cuando la
   del elemento subintegra), las funciones de forma 2D de cara y la
   cuadratura de cara.

Convenciones de caras (numeración y nodos locales) viven en cada
elemento como atributo de clase, no aquí.

Convención Voigt 3D del proyecto (ADR 0012, ``Reglas.md §5``):
``[ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz]`` con ``γ_ij = 2·ε_ij`` *engineering*.

Todos los símbolos llevan guion bajo: son privados al paquete.
"""
from typing import List

import numpy as np
from numba import njit

from solidum.constants import ZERO_JACOBIAN_TOL
from solidum.core.element import Element, validate_lumping_kwarg
from solidum.core.material import Material
from solidum.core.node import Node
from solidum.math.mass_lumping import lump_hrz
from solidum.registry import QuadratureRegistry


# ---------------------------------------------------------------------------
# Hex8 — hexaedro trilineal isoparamétrico.
# ---------------------------------------------------------------------------

@njit(cache=True)
def _shape_functions_hex8(xi: float, eta: float, zeta: float) -> np.ndarray:
    """Funciones de forma trilineales del Hex8 en (ξ, η, ζ) ∈ [-1, 1]³.

    Orden de nodos VTK_HEXAHEDRON:
        0: (-1, -1, -1)   4: (-1, -1, +1)
        1: (+1, -1, -1)   5: (+1, -1, +1)
        2: (+1, +1, -1)   6: (+1, +1, +1)
        3: (-1, +1, -1)   7: (-1, +1, +1)
    """
    N = np.zeros(8, dtype=np.float64)
    xi_sign = (-1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0)
    eta_sign = (-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0)
    zeta_sign = (-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0)
    for i in range(8):
        N[i] = 0.125 * (1.0 + xi_sign[i] * xi) \
                     * (1.0 + eta_sign[i] * eta) \
                     * (1.0 + zeta_sign[i] * zeta)
    return N


@njit(cache=True)
def _compute_kinematics_hex8(xi: float, eta: float, zeta: float,
                             coords: np.ndarray):
    """Calcula B (6×24) y det(J) del Hex8 en (ξ, η, ζ).

    ``coords`` es la matriz 8×3 de coordenadas globales de los nodos
    en orden VTK. Devuelve la matriz B en convención Voigt 6D del
    proyecto (ADR 0012) y el determinante del Jacobiano.
    """
    xi_sign = (-1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0)
    eta_sign = (-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0)
    zeta_sign = (-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0)

    dN_dxi = np.zeros((3, 8), dtype=np.float64)
    for i in range(8):
        sx = xi_sign[i]
        sy = eta_sign[i]
        sz = zeta_sign[i]
        dN_dxi[0, i] = 0.125 * sx * (1.0 + sy * eta) * (1.0 + sz * zeta)
        dN_dxi[1, i] = 0.125 * sy * (1.0 + sx * xi) * (1.0 + sz * zeta)
        dN_dxi[2, i] = 0.125 * sz * (1.0 + sx * xi) * (1.0 + sy * eta)

    J = dN_dxi @ coords  # 3×3

    # det(J) por expansión de la primera fila.
    detJ = (
        J[0, 0] * (J[1, 1] * J[2, 2] - J[1, 2] * J[2, 1])
        - J[0, 1] * (J[1, 0] * J[2, 2] - J[1, 2] * J[2, 0])
        + J[0, 2] * (J[1, 0] * J[2, 1] - J[1, 1] * J[2, 0])
    )

    if detJ <= ZERO_JACOBIAN_TOL:
        raise ValueError(
            "Jacobiano negativo o cero detectado en Hex8. Revisa la "
            "conectividad o la distorsión."
        )

    # Inversa 3×3 por cofactores.
    invJ = np.zeros((3, 3), dtype=np.float64)
    invJ[0, 0] = (J[1, 1] * J[2, 2] - J[1, 2] * J[2, 1]) / detJ
    invJ[0, 1] = -(J[0, 1] * J[2, 2] - J[0, 2] * J[2, 1]) / detJ
    invJ[0, 2] = (J[0, 1] * J[1, 2] - J[0, 2] * J[1, 1]) / detJ
    invJ[1, 0] = -(J[1, 0] * J[2, 2] - J[1, 2] * J[2, 0]) / detJ
    invJ[1, 1] = (J[0, 0] * J[2, 2] - J[0, 2] * J[2, 0]) / detJ
    invJ[1, 2] = -(J[0, 0] * J[1, 2] - J[0, 2] * J[1, 0]) / detJ
    invJ[2, 0] = (J[1, 0] * J[2, 1] - J[1, 1] * J[2, 0]) / detJ
    invJ[2, 1] = -(J[0, 0] * J[2, 1] - J[0, 1] * J[2, 0]) / detJ
    invJ[2, 2] = (J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0]) / detJ

    # Derivadas en globales (3×8).
    dN_dx = invJ @ dN_dxi

    # B en Voigt 6D del proyecto: [xx, yy, zz, xy, yz, xz].
    B = np.zeros((6, 24), dtype=np.float64)
    for i in range(8):
        bx = dN_dx[0, i]
        by = dN_dx[1, i]
        bz = dN_dx[2, i]
        c0 = 3 * i
        B[0, c0]     = bx
        B[1, c0 + 1] = by
        B[2, c0 + 2] = bz
        # γ_xy
        B[3, c0]     = by
        B[3, c0 + 1] = bx
        # γ_yz
        B[4, c0 + 1] = bz
        B[4, c0 + 2] = by
        # γ_xz
        B[5, c0]     = bz
        B[5, c0 + 2] = bx
    return B, detJ


@njit(cache=True)
def _det_jacobian_hex8(xi: float, eta: float, zeta: float,
                       coords: np.ndarray) -> float:
    """Solo det(J) del Hex8 (para body load, cuando B no se necesita)."""
    xi_sign = (-1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0)
    eta_sign = (-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0)
    zeta_sign = (-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0)

    dN_dxi = np.zeros((3, 8), dtype=np.float64)
    for i in range(8):
        sx = xi_sign[i]
        sy = eta_sign[i]
        sz = zeta_sign[i]
        dN_dxi[0, i] = 0.125 * sx * (1.0 + sy * eta) * (1.0 + sz * zeta)
        dN_dxi[1, i] = 0.125 * sy * (1.0 + sx * xi) * (1.0 + sz * zeta)
        dN_dxi[2, i] = 0.125 * sz * (1.0 + sx * xi) * (1.0 + sy * eta)

    J = dN_dxi @ coords
    detJ = (
        J[0, 0] * (J[1, 1] * J[2, 2] - J[1, 2] * J[2, 1])
        - J[0, 1] * (J[1, 0] * J[2, 2] - J[1, 2] * J[2, 0])
        + J[0, 2] * (J[1, 0] * J[2, 1] - J[1, 1] * J[2, 0])
    )
    return detJ


# ---------------------------------------------------------------------------
# Tet4 — tetraedro lineal (CST 3D).
# ---------------------------------------------------------------------------

@njit(cache=True)
def _shape_functions_tet4(xi: float, eta: float, zeta: float) -> np.ndarray:
    """Funciones de forma lineales del Tet4 en coordenadas baricéntricas.

    Vértices del tetraedro de referencia:
        0: (0, 0, 0)
        1: (1, 0, 0)
        2: (0, 1, 0)
        3: (0, 0, 1)
    """
    N = np.zeros(4, dtype=np.float64)
    N[0] = 1.0 - xi - eta - zeta
    N[1] = xi
    N[2] = eta
    N[3] = zeta
    return N


@njit(cache=True)
def _compute_kinematics_tet4(coords: np.ndarray):
    """Calcula B (6×12) y det(J) del Tet4. Ambos constantes sobre el elemento.

    ``coords`` es la matriz 4×3 de coordenadas globales de los nodos.
    Devuelve B en convención Voigt 6D y det(J) = 6·V_e con V_e el volumen
    del tetraedro.
    """
    # Derivadas constantes de N_i en coordenadas naturales:
    # ∂N_1/∂(ξ,η,ζ) = (-1, -1, -1)
    # ∂N_2/∂(ξ,η,ζ) = ( 1,  0,  0)
    # ∂N_3/∂(ξ,η,ζ) = ( 0,  1,  0)
    # ∂N_4/∂(ξ,η,ζ) = ( 0,  0,  1)
    dN_dxi = np.zeros((3, 4), dtype=np.float64)
    dN_dxi[0, 0] = -1.0
    dN_dxi[1, 0] = -1.0
    dN_dxi[2, 0] = -1.0
    dN_dxi[0, 1] = 1.0
    dN_dxi[1, 2] = 1.0
    dN_dxi[2, 3] = 1.0

    J = dN_dxi @ coords  # 3×3, constante sobre el elemento

    detJ = (
        J[0, 0] * (J[1, 1] * J[2, 2] - J[1, 2] * J[2, 1])
        - J[0, 1] * (J[1, 0] * J[2, 2] - J[1, 2] * J[2, 0])
        + J[0, 2] * (J[1, 0] * J[2, 1] - J[1, 1] * J[2, 0])
    )

    if detJ <= ZERO_JACOBIAN_TOL:
        raise ValueError(
            "Jacobiano negativo o cero detectado en Tet4. Cuatro nodos "
            "coplanares o en orden invertido."
        )

    invJ = np.zeros((3, 3), dtype=np.float64)
    invJ[0, 0] = (J[1, 1] * J[2, 2] - J[1, 2] * J[2, 1]) / detJ
    invJ[0, 1] = -(J[0, 1] * J[2, 2] - J[0, 2] * J[2, 1]) / detJ
    invJ[0, 2] = (J[0, 1] * J[1, 2] - J[0, 2] * J[1, 1]) / detJ
    invJ[1, 0] = -(J[1, 0] * J[2, 2] - J[1, 2] * J[2, 0]) / detJ
    invJ[1, 1] = (J[0, 0] * J[2, 2] - J[0, 2] * J[2, 0]) / detJ
    invJ[1, 2] = -(J[0, 0] * J[1, 2] - J[0, 2] * J[1, 0]) / detJ
    invJ[2, 0] = (J[1, 0] * J[2, 1] - J[1, 1] * J[2, 0]) / detJ
    invJ[2, 1] = -(J[0, 0] * J[2, 1] - J[0, 1] * J[2, 0]) / detJ
    invJ[2, 2] = (J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0]) / detJ

    dN_dx = invJ @ dN_dxi  # 3×4, constante

    B = np.zeros((6, 12), dtype=np.float64)
    for i in range(4):
        bx = dN_dx[0, i]
        by = dN_dx[1, i]
        bz = dN_dx[2, i]
        c0 = 3 * i
        B[0, c0]     = bx
        B[1, c0 + 1] = by
        B[2, c0 + 2] = bz
        B[3, c0]     = by
        B[3, c0 + 1] = bx
        B[4, c0 + 1] = bz
        B[4, c0 + 2] = by
        B[5, c0]     = bz
        B[5, c0 + 2] = bx
    return B, detJ


# ---------------------------------------------------------------------------
# Integrandos 3D — comunes a Hex8, Tet4.
# ---------------------------------------------------------------------------

@njit(cache=True)
def _compute_integrands_3d(B: np.ndarray, C: np.ndarray, sigma: np.ndarray,
                           detJ: float, w: float):
    """Contribuciones K_contrib y F_contrib de un punto de Gauss 3D.

    Sin thickness: el volumen ya está incluido en detJ·w.
    """
    factor = detJ * w
    K_contrib = B.T @ C @ B * factor
    F_contrib = B.T @ sigma * factor
    return K_contrib, F_contrib


# ---------------------------------------------------------------------------
# Hex20 — hexaedro serendípito 3D de orden 2 (20 nodos).
#
# Numeración VTK_QUADRATIC_HEXAHEDRON:
#   0-7   : vértices (idéntico Hex8)
#   8-11  : medios de aristas de la cara inferior (z=-1)
#             8 entre (0,1), 9 entre (1,2), 10 entre (2,3), 11 entre (3,0)
#   12-15 : medios de aristas de la cara superior (z=+1)
#             12 entre (4,5), 13 entre (5,6), 14 entre (6,7), 15 entre (7,4)
#   16-19 : medios de aristas verticales (entre caras inferior y superior)
#             16 entre (0,4), 17 entre (1,5), 18 entre (2,6), 19 entre (3,7)
#
# Funciones de forma serendípitas (Cook-Malkus-Plesha-Witt §6.5):
#   - Vértice i con (xi_i, eta_i, zeta_i) ∈ {±1}³:
#       N_i = (1/8)(1+u)(1+v)(1+w)(u+v+w-2)
#       con u=xi_i·xi, v=eta_i·eta, w=zeta_i·zeta.
#   - Medio paralelo a xi (xi_i=0): N = (1/4)(1-xi²)(1+v)(1+w).
#   - Medio paralelo a eta (eta_i=0): N = (1/4)(1+u)(1-eta²)(1+w).
#   - Medio paralelo a zeta (zeta_i=0): N = (1/4)(1+u)(1+v)(1-zeta²).
# ---------------------------------------------------------------------------

# Signos de las coordenadas naturales por nodo en el orden VTK_QUADRATIC_HEXAHEDRON.
# Para vértices: ±1; para medios: 0 en el eje del cual son punto medio.
_HEX20_XI = (
    -1.0,  1.0,  1.0, -1.0, -1.0,  1.0,  1.0, -1.0,  # 0-7 vértices
     0.0,  1.0,  0.0, -1.0,                          # 8-11 cara inferior
     0.0,  1.0,  0.0, -1.0,                          # 12-15 cara superior
    -1.0,  1.0,  1.0, -1.0,                          # 16-19 verticales
)
_HEX20_ETA = (
    -1.0, -1.0,  1.0,  1.0, -1.0, -1.0,  1.0,  1.0,
    -1.0,  0.0,  1.0,  0.0,
    -1.0,  0.0,  1.0,  0.0,
    -1.0, -1.0,  1.0,  1.0,
)
_HEX20_ZETA = (
    -1.0, -1.0, -1.0, -1.0,  1.0,  1.0,  1.0,  1.0,
    -1.0, -1.0, -1.0, -1.0,
     1.0,  1.0,  1.0,  1.0,
     0.0,  0.0,  0.0,  0.0,
)
# Tipo de nodo: 0=vértice, 1=medio paralelo a xi, 2=medio paralelo a eta,
# 3=medio paralelo a zeta. Se infiere del cero en la terna (xi_i, eta_i, zeta_i).
_HEX20_KIND = (
    0, 0, 0, 0, 0, 0, 0, 0,         # vértices
    1, 2, 1, 2,                     # 8(xi), 9(eta), 10(xi), 11(eta) — cara inferior
    1, 2, 1, 2,                     # 12(xi), 13(eta), 14(xi), 15(eta) — cara superior
    3, 3, 3, 3,                     # 16-19 verticales
)


def _N_hex20(xi: float, eta: float, zeta: float) -> np.ndarray:
    """Funciones de forma serendípitas del Hex20 en (ξ, η, ζ) ∈ [-1, 1]³.

    Devuelve un array de longitud 20 con ``N[i]`` evaluado en el punto.
    Convención VTK_QUADRATIC_HEXAHEDRON.
    """
    N = np.zeros(20, dtype=np.float64)
    for i in range(20):
        kind = _HEX20_KIND[i]
        xi_i = _HEX20_XI[i]
        eta_i = _HEX20_ETA[i]
        zeta_i = _HEX20_ZETA[i]
        if kind == 0:  # vértice
            u = xi_i * xi
            v = eta_i * eta
            w = zeta_i * zeta
            N[i] = 0.125 * (1.0 + u) * (1.0 + v) * (1.0 + w) * (u + v + w - 2.0)
        elif kind == 1:  # medio paralelo a ξ (xi_i = 0)
            N[i] = 0.25 * (1.0 - xi * xi) * (1.0 + eta_i * eta) * (1.0 + zeta_i * zeta)
        elif kind == 2:  # medio paralelo a η (eta_i = 0)
            N[i] = 0.25 * (1.0 + xi_i * xi) * (1.0 - eta * eta) * (1.0 + zeta_i * zeta)
        else:  # kind == 3: medio paralelo a ζ (zeta_i = 0)
            N[i] = 0.25 * (1.0 + xi_i * xi) * (1.0 + eta_i * eta) * (1.0 - zeta * zeta)
    return N


def _dN_hex20(xi: float, eta: float, zeta: float) -> np.ndarray:
    """Derivadas de las funciones de forma del Hex20 en (ξ, η, ζ).

    Devuelve una matriz 3×20 con ``dN[0,i] = ∂N_i/∂ξ``,
    ``dN[1,i] = ∂N_i/∂η``, ``dN[2,i] = ∂N_i/∂ζ``.
    """
    dN = np.zeros((3, 20), dtype=np.float64)
    for i in range(20):
        kind = _HEX20_KIND[i]
        xi_i = _HEX20_XI[i]
        eta_i = _HEX20_ETA[i]
        zeta_i = _HEX20_ZETA[i]
        if kind == 0:  # vértice
            u = xi_i * xi
            v = eta_i * eta
            w = zeta_i * zeta
            dN[0, i] = 0.125 * xi_i * (1.0 + v) * (1.0 + w) * (2.0 * u + v + w - 1.0)
            dN[1, i] = 0.125 * eta_i * (1.0 + u) * (1.0 + w) * (u + 2.0 * v + w - 1.0)
            dN[2, i] = 0.125 * zeta_i * (1.0 + u) * (1.0 + v) * (u + v + 2.0 * w - 1.0)
        elif kind == 1:  # medio paralelo a ξ
            dN[0, i] = -0.5 * xi * (1.0 + eta_i * eta) * (1.0 + zeta_i * zeta)
            dN[1, i] = 0.25 * eta_i * (1.0 - xi * xi) * (1.0 + zeta_i * zeta)
            dN[2, i] = 0.25 * zeta_i * (1.0 - xi * xi) * (1.0 + eta_i * eta)
        elif kind == 2:  # medio paralelo a η
            dN[0, i] = 0.25 * xi_i * (1.0 - eta * eta) * (1.0 + zeta_i * zeta)
            dN[1, i] = -0.5 * eta * (1.0 + xi_i * xi) * (1.0 + zeta_i * zeta)
            dN[2, i] = 0.25 * zeta_i * (1.0 + xi_i * xi) * (1.0 - eta * eta)
        else:  # kind == 3: medio paralelo a ζ
            dN[0, i] = 0.25 * xi_i * (1.0 + eta_i * eta) * (1.0 - zeta * zeta)
            dN[1, i] = 0.25 * eta_i * (1.0 + xi_i * xi) * (1.0 - zeta * zeta)
            dN[2, i] = -0.5 * zeta * (1.0 + xi_i * xi) * (1.0 + eta_i * eta)
    return dN


# ---------------------------------------------------------------------------
# Hex27 — hexaedro Lagrangiano 3D triquadrático (27 nodos, sub-etapa A.ter).
#
# Numeración VTK_TRIQUADRATIC_HEXAHEDRON:
#   0-7   : vértices (idéntico Hex8/Hex20)
#   8-19  : medios de arista (idéntico Hex20)
#   20    : centro de cara en ξ = -1 (cara 5 de la convención Hex8: -x)
#   21    : centro de cara en ξ = +1 (cara 3: +x)
#   22    : centro de cara en η = -1 (cara 2: -y)
#   23    : centro de cara en η = +1 (cara 4: +y)
#   24    : centro de cara en ζ = -1 (cara 0: -z)
#   25    : centro de cara en ζ = +1 (cara 1: +z)
#   26    : centro del cuerpo (ξ = η = ζ = 0)
#
# Funciones de forma: producto tensorial de Lagrange cuadráticos:
#   L_{-1}(x) = x(x-1)/2,  L_0(x) = 1 - x²,  L_{+1}(x) = x(x+1)/2.
#   N_{ijk}(ξ, η, ζ) = L_i(ξ) · L_j(η) · L_k(ζ)  con i, j, k ∈ {-1, 0, +1}.
# ---------------------------------------------------------------------------

# (i, j, k) signs per node — orden VTK_TRIQUADRATIC_HEXAHEDRON.
_HEX27_I = (
    -1, 1, 1, -1, -1, 1, 1, -1,   # 0-7 vértices
     0, 1, 0, -1,                  # 8-11 bottom edge mids
     0, 1, 0, -1,                  # 12-15 top edge mids
    -1, 1, 1, -1,                  # 16-19 vertical edge mids
    -1, 1,                         # 20-21 face centers -x, +x
     0, 0,                         # 22-23 face centers -y, +y
     0, 0,                         # 24-25 face centers -z, +z
     0,                            # 26 body center
)
_HEX27_J = (
    -1, -1, 1, 1, -1, -1, 1, 1,
    -1, 0, 1, 0,
    -1, 0, 1, 0,
    -1, -1, 1, 1,
     0, 0,
    -1, 1,
     0, 0,
     0,
)
_HEX27_K = (
    -1, -1, -1, -1, 1, 1, 1, 1,
    -1, -1, -1, -1,
     1, 1, 1, 1,
     0, 0, 0, 0,
     0, 0,
     0, 0,
    -1, 1,
     0,
)


def _L_lagrange_quad(idx: int, x: float) -> float:
    """Lagrange cuadrático 1D evaluado en ``x`` para índice ``idx ∈ {-1, 0, +1}``."""
    if idx == -1:
        return 0.5 * x * (x - 1.0)
    if idx == 0:
        return 1.0 - x * x
    return 0.5 * x * (x + 1.0)


def _dL_lagrange_quad(idx: int, x: float) -> float:
    """Derivada del Lagrange cuadrático 1D."""
    if idx == -1:
        return x - 0.5
    if idx == 0:
        return -2.0 * x
    return x + 0.5


def _N_hex27(xi: float, eta: float, zeta: float) -> np.ndarray:
    """Funciones de forma triquadráticas Lagrange del Hex27.

    Devuelve un array de longitud 27 con ``N[i]`` en el punto natural
    ``(ξ, η, ζ) ∈ [-1, 1]³``. Convención VTK_TRIQUADRATIC_HEXAHEDRON.
    """
    N = np.zeros(27, dtype=np.float64)
    for n in range(27):
        N[n] = (
            _L_lagrange_quad(_HEX27_I[n], xi)
            * _L_lagrange_quad(_HEX27_J[n], eta)
            * _L_lagrange_quad(_HEX27_K[n], zeta)
        )
    return N


def _dN_hex27(xi: float, eta: float, zeta: float) -> np.ndarray:
    """Derivadas de las funciones de forma del Hex27 en (ξ, η, ζ).

    Devuelve una matriz 3×27 con ``dN[0,i] = ∂N_i/∂ξ``, etc.
    """
    dN = np.zeros((3, 27), dtype=np.float64)
    for n in range(27):
        i_idx = _HEX27_I[n]
        j_idx = _HEX27_J[n]
        k_idx = _HEX27_K[n]
        Lx = _L_lagrange_quad(i_idx, xi)
        Ly = _L_lagrange_quad(j_idx, eta)
        Lz = _L_lagrange_quad(k_idx, zeta)
        dLx = _dL_lagrange_quad(i_idx, xi)
        dLy = _dL_lagrange_quad(j_idx, eta)
        dLz = _dL_lagrange_quad(k_idx, zeta)
        dN[0, n] = dLx * Ly * Lz
        dN[1, n] = Lx * dLy * Lz
        dN[2, n] = Lx * Ly * dLz
    return dN


# ---------------------------------------------------------------------------
# Tet10 — tetraedro cuadrático isoparamétrico (10 nodos, sub-etapa A.ter).
#
# Numeración VTK_QUADRATIC_TETRA:
#   0-3 : vértices (idéntico Tet4)
#       0: (0, 0, 0),  1: (1, 0, 0),  2: (0, 1, 0),  3: (0, 0, 1)
#   4-9 : medios de arista
#       4: medio arista 0-1  →  (0.5, 0, 0)
#       5: medio arista 1-2  →  (0.5, 0.5, 0)
#       6: medio arista 2-0  →  (0, 0.5, 0)
#       7: medio arista 0-3  →  (0, 0, 0.5)
#       8: medio arista 1-3  →  (0.5, 0, 0.5)
#       9: medio arista 2-3  →  (0, 0.5, 0.5)
#
# Coordenadas baricéntricas en (ξ, η, ζ):
#   L_1 = 1 − ξ − η − ζ,  L_2 = ξ,  L_3 = η,  L_4 = ζ
#
# Funciones de forma cuadráticas:
#   - Vértices: N_i = L_i (2 L_i - 1)  con i ∈ {1, 2, 3, 4}.
#   - Medios:   N_ij = 4 L_i L_j       con (i, j) la arista correspondiente.
# ---------------------------------------------------------------------------

# Mapeo nodo local → barycentric flags: para cada nodo, la lista de los
# L_k que aparecen multiplicados. Vértice 0 ↔ L_1, vértice 1 ↔ L_2, etc.
# Para medios, dos L_k.
_TET10_VERTEX_L = (1, 2, 3, 4)  # nodo i → índice barycentric (1..4)
_TET10_EDGE_LS = (
    (1, 2),  # 4: medio 0-1
    (2, 3),  # 5: medio 1-2
    (3, 1),  # 6: medio 2-0  (= 0-2)
    (1, 4),  # 7: medio 0-3
    (2, 4),  # 8: medio 1-3
    (3, 4),  # 9: medio 2-3
)


def _N_tet10(xi: float, eta: float, zeta: float) -> np.ndarray:
    """Funciones de forma cuadráticas del Tet10 en coordenadas naturales
    ``(ξ, η, ζ)`` (= ``(L_2, L_3, L_4)`` baricéntricas).

    Convención VTK_QUADRATIC_TETRA.
    """
    L = (1.0 - xi - eta - zeta, xi, eta, zeta)  # L_1..L_4 indexable 0..3
    N = np.zeros(10, dtype=np.float64)
    # Vértices: N_i = L_i (2 L_i - 1).
    for i in range(4):
        N[i] = L[i] * (2.0 * L[i] - 1.0)
    # Medios: N_ij = 4 L_i L_j.
    for k, (i_idx, j_idx) in enumerate(_TET10_EDGE_LS):
        # _TET10_EDGE_LS usa índices 1..4; convertimos a 0..3.
        N[4 + k] = 4.0 * L[i_idx - 1] * L[j_idx - 1]
    return N


def _dN_tet10(xi: float, eta: float, zeta: float) -> np.ndarray:
    """Derivadas de las funciones de forma del Tet10 en (ξ, η, ζ).

    Devuelve matriz 3×10 con ``dN[0,i] = ∂N_i/∂ξ``, etc.

    Recordando que en (ξ, η, ζ):
        ∂L_1/∂(ξ,η,ζ) = (-1, -1, -1)
        ∂L_2/∂(ξ,η,ζ) = ( 1,  0,  0)
        ∂L_3/∂(ξ,η,ζ) = ( 0,  1,  0)
        ∂L_4/∂(ξ,η,ζ) = ( 0,  0,  1)
    """
    L = (1.0 - xi - eta - zeta, xi, eta, zeta)
    # dL_k/d(ξ, η, ζ) por k = 0..3:
    dL = (
        (-1.0, -1.0, -1.0),   # L_1
        ( 1.0,  0.0,  0.0),   # L_2
        ( 0.0,  1.0,  0.0),   # L_3
        ( 0.0,  0.0,  1.0),   # L_4
    )
    dN = np.zeros((3, 10), dtype=np.float64)
    # Vértices: N_i = L_i(2L_i-1) → dN_i/dx_k = (4L_i - 1) · dL_i/dx_k.
    for i in range(4):
        factor = 4.0 * L[i] - 1.0
        for axis in range(3):
            dN[axis, i] = factor * dL[i][axis]
    # Medios: N_ij = 4 L_i L_j → dN_ij/dx_k = 4 (L_j dL_i/dx_k + L_i dL_j/dx_k).
    for k, (i_idx, j_idx) in enumerate(_TET10_EDGE_LS):
        i_ = i_idx - 1
        j_ = j_idx - 1
        for axis in range(3):
            dN[axis, 4 + k] = 4.0 * (L[j_] * dL[i_][axis] + L[i_] * dL[j_][axis])
    return dN


# ---------------------------------------------------------------------------
# Kinematics genérico de orden superior 3D (sub-etapa A.ter).
# ---------------------------------------------------------------------------

def _kinematics_higher_order_3d(grad_fn, xi: float, eta: float, zeta: float,
                                coords: np.ndarray, n_nodes: int):
    """Calcula B (6×3·n_nodes) y det(J) para un sólido 3D de orden superior.

    Paritario con :func:`solidum.elements.solid_2d._shared._kinematics_higher_order`,
    pero en 3D y con la convención Voigt 6D del proyecto (ADR 0012):
    ``[ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz]``.

    Parameters
    ----------
    grad_fn : callable
        Función ``(xi, eta, zeta) -> ndarray (3, n_nodes)`` con las derivadas
        en coordenadas naturales de las funciones de forma del elemento
        concreto (e.g. :func:`_dN_hex20`).
    xi, eta, zeta : float
        Coordenadas naturales del punto de Gauss.
    coords : ndarray (n_nodes, 3)
        Coordenadas globales de los nodos del elemento.
    n_nodes : int
        Número de nodos del elemento (necesario para dimensionar B sin
        depender del shape de ``coords`` exclusivamente).

    Returns
    -------
    B : ndarray (6, 3·n_nodes)
        Matriz deformación-desplazamiento en convención Voigt 6D del proyecto.
    detJ : float
        Determinante del Jacobiano 3×3.
    """
    dN_dxi = grad_fn(xi, eta, zeta)  # (3, n_nodes)
    J = dN_dxi @ coords              # (3, 3)
    detJ = (
        J[0, 0] * (J[1, 1] * J[2, 2] - J[1, 2] * J[2, 1])
        - J[0, 1] * (J[1, 0] * J[2, 2] - J[1, 2] * J[2, 0])
        + J[0, 2] * (J[1, 0] * J[2, 1] - J[1, 1] * J[2, 0])
    )
    if detJ <= ZERO_JACOBIAN_TOL:
        raise ValueError(
            "Jacobiano negativo o cero detectado en sólido 3D de orden "
            "superior. Revisa la conectividad o la distorsión."
        )
    # Inversa 3×3 por cofactores.
    invJ = np.empty((3, 3), dtype=np.float64)
    invJ[0, 0] = (J[1, 1] * J[2, 2] - J[1, 2] * J[2, 1]) / detJ
    invJ[0, 1] = -(J[0, 1] * J[2, 2] - J[0, 2] * J[2, 1]) / detJ
    invJ[0, 2] = (J[0, 1] * J[1, 2] - J[0, 2] * J[1, 1]) / detJ
    invJ[1, 0] = -(J[1, 0] * J[2, 2] - J[1, 2] * J[2, 0]) / detJ
    invJ[1, 1] = (J[0, 0] * J[2, 2] - J[0, 2] * J[2, 0]) / detJ
    invJ[1, 2] = -(J[0, 0] * J[1, 2] - J[0, 2] * J[1, 0]) / detJ
    invJ[2, 0] = (J[1, 0] * J[2, 1] - J[1, 1] * J[2, 0]) / detJ
    invJ[2, 1] = -(J[0, 0] * J[2, 1] - J[0, 1] * J[2, 0]) / detJ
    invJ[2, 2] = (J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0]) / detJ
    dN_dx = invJ @ dN_dxi  # (3, n_nodes) en globales

    B = np.zeros((6, 3 * n_nodes), dtype=np.float64)
    for i in range(n_nodes):
        bx = dN_dx[0, i]
        by = dN_dx[1, i]
        bz = dN_dx[2, i]
        c0 = 3 * i
        B[0, c0]     = bx
        B[1, c0 + 1] = by
        B[2, c0 + 2] = bz
        # γ_xy
        B[3, c0]     = by
        B[3, c0 + 1] = bx
        # γ_yz
        B[4, c0 + 1] = bz
        B[4, c0 + 2] = by
        # γ_xz
        B[5, c0]     = bz
        B[5, c0 + 2] = bx
    return B, detJ


# ---------------------------------------------------------------------------
# Expansor escalar → bloque 3D para masa traslacional.
# ---------------------------------------------------------------------------

def _expand_scalar_mass_3d(M_scalar: np.ndarray) -> np.ndarray:
    """Expande una matriz escalar (n_nodes × n_nodes) a una matriz 3D
    traslacional (3·n_nodes × 3·n_nodes) por producto de Kronecker con I_3.

    Análogo 3D de :func:`solidum.elements.solid_2d._shared._expand_scalar_mass`.
    Usado por ``Hex8.compute_mass_matrix`` y ``Tet4.compute_mass_matrix``.
    """
    n = M_scalar.shape[0]
    M = np.zeros((3 * n, 3 * n), dtype=np.float64)
    for i in range(n):
        for j in range(n):
            v = M_scalar[i, j]
            M[3 * i,     3 * j]     = v
            M[3 * i + 1, 3 * j + 1] = v
            M[3 * i + 2, 3 * j + 2] = v
    return M


# ---------------------------------------------------------------------------
# Base interna compartida por los sólidos 3D isoparamétricos de orden
# superior: Hex20, Hex27 (sub-etapa A.ter), y Tet10 (sub-fase 3).
# ---------------------------------------------------------------------------

class _HigherOrderSolid3D(Element):
    """Base interna para sólidos 3D isoparamétricos de orden superior.

    Comparte los bucles de Gauss para ``K``, ``F_int``, ``gauss_state``,
    ``body_load``, ``face_traction`` y masa consistente. Cada subclase
    declara:

    - ``_SHAPE_FN`` / ``_GRAD_FN`` — funciones de forma y derivadas 3D.
    - ``_DEFAULT_QUADRATURE`` — nombre de la cuadratura por defecto en
      ``QuadratureRegistry`` (e.g. ``"hex_3x3x3"``).
    - ``_MASS_QUADRATURE`` — cuadratura específica para la matriz de masa
      (independiente de la elegida para ``K``). ``None`` = usa la del
      elemento.
    - ``FACE_NODES`` — tupla de tuplas con los nodos locales por cara.
      Cada tupla tiene tantos elementos como nodos por cara (8 en
      Hex20, 9 en Hex27, 6 en Tet10).
    - ``_FACE_N_FN`` / ``_FACE_DN_FN`` — funciones de forma 2D de cara
      y su gradiente (e.g. Quad8 para Hex20, Quad9 para Hex27, Tri6
      para Tet10).
    - ``_FACE_QUADRATURE`` — nombre de la cuadratura 2D de cara en
      ``QuadratureRegistry`` (e.g. ``"3x3"`` para cuadriláteros,
      ``"tri_6"`` para triángulos cuadráticos).

    El número de nodos se infiere de ``len(self.nodes)``: la base no
    presupone forma (cubo, tetraedro, etc.).

    Convención Voigt 6D del proyecto (ADR 0012):
    ``[ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz]`` con ``γ_ij = 2·ε_ij``.
    """

    DOF_NAMES = ['ux', 'uy', 'uz']
    STRAIN_DIM = 6

    # Atributos a sobreescribir por cada subclase concreta:
    _SHAPE_FN = staticmethod(lambda xi, eta, zeta: None)
    _GRAD_FN = staticmethod(lambda xi, eta, zeta: None)
    _DEFAULT_QUADRATURE: str = ""
    _MASS_QUADRATURE: str | None = None
    _FACE_N_FN = staticmethod(lambda s, t: None)
    _FACE_DN_FN = staticmethod(lambda s, t: None)
    _FACE_QUADRATURE: str = ""
    FACE_NODES: tuple = ()

    def __init__(self, element_id: int, nodes: List[Node], material: Material,
                 quadrature: str | None = None):
        if quadrature is None:
            quadrature = self._DEFAULT_QUADRATURE
        pts, ws = QuadratureRegistry.get(quadrature)
        self.points = pts
        self.weights = ws
        self.N_INTEGRATION_POINTS = len(self.points)
        super().__init__(element_id, nodes, material)

    @property
    def _n_nodes(self) -> int:
        return len(self.nodes)

    # ------------------------------------------------------------------
    # Contrato principal
    # ------------------------------------------------------------------

    def compute_element_state(self, u_e: np.ndarray):
        n = self._n_nodes
        ndof = 3 * n
        K_e = np.zeros((ndof, ndof))
        F_int_e = np.zeros(ndof)
        coords = self.get_coordinate_matrix(ndim=3)
        for idx, ((xi, eta, zeta), w) in enumerate(zip(self.points, self.weights)):
            B, detJ = _kinematics_higher_order_3d(
                self._GRAD_FN, xi, eta, zeta, coords, n
            )
            strain = B @ u_e
            sigma, C_tangent, new_state = self.material.compute_state(
                strain, self.state.vars[idx]
            )
            self.state.vars_trial[idx] = new_state
            self.state.stresses_trial[idx] = sigma
            factor = detJ * w
            K_e += (B.T @ C_tangent @ B) * factor
            F_int_e += (B.T @ sigma) * factor
        return K_e, F_int_e

    # ------------------------------------------------------------------
    # Salida por punto de Gauss (ADR 0012)
    # ------------------------------------------------------------------

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        """Promedio del campo σ y ε sobre los puntos de Gauss del elemento.

        Para acceso no promediado por punto de integración úsese
        :meth:`compute_gauss_state`.
        """
        gs = self.compute_gauss_state(U_global)
        return {
            'stress': gs['stress'].mean(axis=0),
            'strain': gs['strain'].mean(axis=0),
        }

    def compute_gauss_state(self, U_global: np.ndarray) -> dict:
        """Estado por punto de Gauss para sólidos 3D de orden superior."""
        u_e = self.get_local_displacements(U_global)
        coords = self.get_coordinate_matrix(ndim=3)
        n = self._n_nodes
        n_g = self.N_INTEGRATION_POINTS

        nat = np.asarray(self.points, dtype=np.float64).reshape(n_g, 3)
        glb = np.zeros((n_g, 3))
        eps = np.zeros((n_g, 6))
        sig = np.zeros((n_g, 6))

        for idx, (xi, eta, zeta) in enumerate(self.points):
            N = self._SHAPE_FN(xi, eta, zeta)
            glb[idx] = N @ coords
            B, _ = _kinematics_higher_order_3d(
                self._GRAD_FN, xi, eta, zeta, coords, n
            )
            strain = B @ u_e
            stress, _, _ = self.material.compute_state(strain, self.state.vars[idx])
            eps[idx] = strain
            sig[idx] = stress
        return {
            'points_natural': nat,
            'points_global': glb,
            'strain': eps,
            'stress': sig,
        }

    # ------------------------------------------------------------------
    # Cargas distribuidas consistentes
    # ------------------------------------------------------------------

    def compute_body_load(self, b: np.ndarray) -> np.ndarray:
        """Vector consistente con fuerza de cuerpo uniforme.

        Integra ``∫ N^T b dV`` con la cuadratura del elemento. ``b`` en
        globales ``(b_x, b_y, b_z)``. Reparto entre nodos no uniforme
        para elementos de orden superior; suma global ``b·V_e`` exacta
        para ``b`` uniforme y geometría sin distorsión severa.
        """
        b = np.asarray(b, dtype=np.float64).reshape(3)
        coords = self.get_coordinate_matrix(ndim=3)
        n = self._n_nodes
        f = np.zeros(3 * n)
        for (xi, eta, zeta), w in zip(self.points, self.weights):
            N = self._SHAPE_FN(xi, eta, zeta)
            _, detJ = _kinematics_higher_order_3d(
                self._GRAD_FN, xi, eta, zeta, coords, n
            )
            factor = detJ * w
            for i in range(n):
                f[3 * i]     += N[i] * b[0] * factor
                f[3 * i + 1] += N[i] * b[1] * factor
                f[3 * i + 2] += N[i] * b[2] * factor
        return f

    def compute_face_traction(self, face: int, t_vec: np.ndarray) -> np.ndarray:
        """Tracción uniforme sobre la cara ``face`` integrada con la cuadratura
        2D de cara (``_FACE_QUADRATURE``).

        Cada subclase declara:
        - ``FACE_NODES[face]`` — índices locales de los nodos de la cara.
        - ``_FACE_N_FN(s, t)``, ``_FACE_DN_FN(s, t)`` — funciones de forma
          2D de cara y su gradiente (vector dirección a lo largo de los
          dos parámetros).

        ``t_vec`` se asume **constante** sobre la cara y se especifica en
        coordenadas globales ``(t_x, t_y, t_z)``.
        """
        if face < 0 or face >= len(self.FACE_NODES):
            raise ValueError(
                f"face={face} fuera de rango para "
                f"{type(self).__name__} (0..{len(self.FACE_NODES) - 1})."
            )
        t_vec = np.asarray(t_vec, dtype=np.float64).reshape(3)

        face_local = self.FACE_NODES[face]
        n_face = len(face_local)
        X = np.array(
            [self.nodes[i].coordinates[:3] for i in face_local],
            dtype=np.float64,
        )  # n_face × 3

        face_pts, face_wts = QuadratureRegistry.get(self._FACE_QUADRATURE)

        n = self._n_nodes
        f = np.zeros(3 * n)
        for (sp, tp), w in zip(face_pts, face_wts):
            N_face = self._FACE_N_FN(sp, tp)        # (n_face,)
            dN_face = self._FACE_DN_FN(sp, tp)      # (2, n_face)
            g_s = dN_face[0] @ X                    # (3,)
            g_t = dN_face[1] @ X                    # (3,)
            dA = float(np.linalg.norm(np.cross(g_s, g_t)))
            factor = dA * w
            for k, node_local in enumerate(face_local):
                base = 3 * node_local
                contrib = N_face[k] * factor
                f[base]     += contrib * t_vec[0]
                f[base + 1] += contrib * t_vec[1]
                f[base + 2] += contrib * t_vec[2]
        return f

    # ------------------------------------------------------------------
    # Masa elemental (ADR 0009)
    # ------------------------------------------------------------------

    def compute_mass_matrix(self, lumping: str = "consistent") -> np.ndarray:
        """Masa del sólido 3D de orden superior.

        **Consistente** (default): cuando la subclase declara
        ``_MASS_QUADRATURE`` (no None), se usa esa cuadratura
        independientemente de la elegida para ``K``. Si es ``None`` se
        usa la del elemento — el caso típico cuando la cuadratura del
        elemento ya integra exactamente la masa.

        **Lumped**: HRZ canónico (ADR 0009 fase 2). Para elementos de
        orden superior con nodos intermedios/face/body, HRZ es la única
        opción razonable: row-sum produciría masas negativas en
        vértices (Bathe FEP §9.2.4).

        Returns
        -------
        np.ndarray, shape (3·n_nodes, 3·n_nodes)
        """
        rho = self.material.density
        coords = self.get_coordinate_matrix(ndim=3)
        n = self._n_nodes

        if self._MASS_QUADRATURE is None:
            mass_pts, mass_wts = self.points, self.weights
        else:
            mass_pts, mass_wts = QuadratureRegistry.get(self._MASS_QUADRATURE)

        M_s = np.zeros((n, n))
        m_total = 0.0
        for (xi, eta, zeta), w in zip(mass_pts, mass_wts):
            N = self._SHAPE_FN(xi, eta, zeta)
            _, detJ = _kinematics_higher_order_3d(
                self._GRAD_FN, xi, eta, zeta, coords, n
            )
            weight = detJ * w
            M_s += rho * np.outer(N, N) * weight
            m_total += rho * weight
        M_consistent = _expand_scalar_mass_3d(M_s)

        validate_lumping_kwarg(lumping, type(self).__name__)
        if lumping == "lumped":
            return lump_hrz(M_consistent, total_mass=m_total, n_translational_dirs=3)
        return M_consistent
