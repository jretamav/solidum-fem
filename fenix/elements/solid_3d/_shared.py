"""Utilidades compartidas por los elementos del paquete ``fenix.elements.solid_3d`` (ADR 0012).

Tres grupos:

1. **Kinematics Numba** para los elementos lineales (Hex8, Tet4):
   :func:`_compute_kinematics_hex8`, :func:`_compute_kinematics_tet4`,
   :func:`_compute_integrands_3d`, :func:`_shape_functions_hex8`,
   :func:`_det_jacobian_hex8`, :func:`_shape_functions_tet4`.
2. **Expansor escalar→bloque 3D** para masa traslacional:
   :func:`_expand_scalar_mass_3d`.
3. **Convenciones de caras** (numeración y nodos locales) que ambos
   elementos usan para ``compute_face_traction``: viven en cada elemento
   como atributo de clase, no aquí.

Convención Voigt 3D del proyecto (ADR 0012, ``Reglas.md §5``):
``[ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz]`` con ``γ_ij = 2·ε_ij`` *engineering*.

Todos los símbolos llevan guion bajo: son privados al paquete.
"""
import numpy as np
from numba import njit

from fenix.constants import ZERO_JACOBIAN_TOL


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
# Expansor escalar → bloque 3D para masa traslacional.
# ---------------------------------------------------------------------------

def _expand_scalar_mass_3d(M_scalar: np.ndarray) -> np.ndarray:
    """Expande una matriz escalar (n_nodes × n_nodes) a una matriz 3D
    traslacional (3·n_nodes × 3·n_nodes) por producto de Kronecker con I_3.

    Análogo 3D de :func:`fenix.elements.solid_2d._shared._expand_scalar_mass`.
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
