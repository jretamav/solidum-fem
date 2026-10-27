"""``Hex8`` — hexaedro trilineal 3D isoparamétrico (ADR 0012, Etapa 7)."""
from typing import List

import numpy as np

from solidum.core.element import Element, validate_lumping_kwarg
from solidum.core.material import Material
from solidum.core.node import Node
from solidum.elements.solid_2d._shared import _dN_quad4, _shape_functions_quad4
from solidum.elements.solid_3d._shared import (
    _SolidLoads3D,
    _batch_kin_hex8,
    _compute_integrands_3d,
    _compute_kinematics_hex8,
    _det_jacobian_hex8,
    _expand_scalar_mass_3d,
    _shape_functions_hex8,
)
from solidum.math.integration import resolve_quadrature
from solidum.math.mass_lumping import lump_hrz
from solidum.registry import ElementRegistry, QuadratureRegistry


@ElementRegistry.register
class Hex8(_SolidLoads3D, Element):
    """Elemento hexaédrico trilineal 3D isoparamétrico.

    Parameters
    ----------
    element_id : int
    nodes : List[Node]
        8 nodos en orden VTK_HEXAHEDRON (0-3 cara inferior antihorario
        visto desde dentro, 4-7 cara superior antihorario visto desde fuera).
    material : Material
        Material 3D (STRAIN_DIM=6).
    quadrature : str, optional
        Nombre de cuadratura desde ``QuadratureRegistry``. Default ``hex_2x2x2``
        (8 puntos). Alternativas: ``hex_3x3x3`` (27 puntos, no lineales severos),
        ``hex_1x1x1`` (1 punto, reducido — propenso a hourglass).

    Notes
    -----
    Locking volumétrico con ν → 0.5 y hourglass con integración reducida son
    limitaciones declaradas en la spec; sin mitigación implementada.
    """

    DOF_NAMES = ['ux', 'uy', 'uz']
    STRAIN_DIM = 6
    N_INTEGRATION_POINTS = 8  # default Gauss 2×2×2
    # Camino por lotes (ADR 0014): misma cinemática compilada que
    # compute_element_state; sin escala (el volumen va en detJ·w).
    BATCH_KINEMATICS = _batch_kin_hex8
    BATCH_SHAPE_FUNCTIONS = staticmethod(_shape_functions_hex8)
    REFERENCE_CELL = "hex"
    REFERENCE_NODES = (
        (-1.0, -1.0, -1.0), (1.0, -1.0, -1.0), (1.0, 1.0, -1.0), (-1.0, 1.0, -1.0),
        (-1.0, -1.0, 1.0), (1.0, -1.0, 1.0), (1.0, 1.0, 1.0), (-1.0, 1.0, 1.0),
    )

    # ADR 0012 — 6 caras con normal saliente. Cada tupla son los 4 nodos
    # locales de la cara en orden antihorario vistos desde fuera, es decir
    # con (b−a)×(c−a) apuntando hacia fuera: la normal que recibe un campo
    # de tracción t(x, n) (ADR 0022).
    FACE_NODES = (
        (0, 3, 2, 1),  # 0: −ζ (inferior)
        (4, 5, 6, 7),  # 1: +ζ (superior)
        (0, 1, 5, 4),  # 2: −η (frontal)
        (1, 2, 6, 5),  # 3: +ξ (derecha)
        (2, 3, 7, 6),  # 4: +η (trasera)
        (3, 0, 4, 7),  # 5: −ξ (izquierda)
    )

    def __init__(self, element_id: int, nodes: List[Node], material: Material,
                 quadrature="hex_2x2x2"):
        self.points, self.weights, self.quadrature_key = resolve_quadrature(
            quadrature, "hex_2x2x2", family="hex",
        )
        self.N_INTEGRATION_POINTS = len(self.points)
        super().__init__(element_id, nodes, material)

    # ------------------------------------------------------------------
    # Contrato principal
    # ------------------------------------------------------------------

    def compute_element_state(self, u_e: np.ndarray):
        K_e = np.zeros((24, 24))
        F_int_e = np.zeros(24)

        coords = self.get_coordinate_matrix(ndim=3)

        for idx, (p, w) in enumerate(zip(self.points, self.weights)):
            xi, eta, zeta = p
            B, detJ = _compute_kinematics_hex8(xi, eta, zeta, coords)
            strain = B @ u_e

            sigma, C_tangent, new_state = self.material.compute_state(
                strain, self.state.vars[idx]
            )
            self.state.vars_trial[idx] = new_state
            self.state.stresses_trial[idx] = sigma

            K_contrib, F_contrib = _compute_integrands_3d(
                B, C_tangent, sigma, detJ, w
            )
            K_e += K_contrib
            F_int_e += F_contrib

        return K_e, F_int_e

    # ------------------------------------------------------------------
    # Salida por punto de Gauss (ADR 0012: sólidos no exponen internal_forces)
    # ------------------------------------------------------------------

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        """Promedio del campo σ y ε sobre los puntos de Gauss (utilidad rápida).

        Para visualizar el campo no promediado úsese ``compute_gauss_state``.
        """
        gs = self.compute_gauss_state(U_global)
        return {
            'stress': gs['stress'].mean(axis=0),
            'strain': gs['strain'].mean(axis=0),
        }

    def compute_gauss_state(self, U_global: np.ndarray) -> dict:
        """Estado en cada punto de Gauss del Hex8.

        Returns
        -------
        dict con claves:
        - ``points_natural`` : ndarray (n_g, 3) con (ξ, η, ζ).
        - ``points_global``  : ndarray (n_g, 3) con (x, y, z) en globales.
        - ``strain`` : ndarray (n_g, 6) — Voigt 3D [ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz].
        - ``stress`` : ndarray (n_g, 6) — Voigt 3D [σ_xx, σ_yy, σ_zz, σ_xy, σ_yz, σ_xz].
        """
        u_e = self.get_local_displacements(U_global)
        coords = self.get_coordinate_matrix(ndim=3)
        n_g = self.N_INTEGRATION_POINTS

        nat = np.asarray(self.points, dtype=np.float64).reshape(n_g, 3)
        glb = np.zeros((n_g, 3))
        eps = np.zeros((n_g, 6))
        sig = np.zeros((n_g, 6))

        for idx, (xi, eta, zeta) in enumerate(self.points):
            N = _shape_functions_hex8(xi, eta, zeta)
            glb[idx] = N @ coords
            B, _ = _compute_kinematics_hex8(xi, eta, zeta, coords)
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
    # Cargas distribuidas consistentes: compute_body_load y
    # compute_face_traction de _SolidLoads3D (ADR 0022), con la regla
    # 2×2×2 de la masa aunque la rigidez sea 1×1×1, y caras bilineales.
    # ------------------------------------------------------------------

    _LOAD_QUADRATURE = "hex_2x2x2"
    _FACE_N_FN = staticmethod(_shape_functions_quad4)
    _FACE_DN_FN = staticmethod(_dN_quad4)
    _FACE_QUADRATURE = "2x2"

    # ------------------------------------------------------------------
    # Masa elemental (ADR 0009)
    # ------------------------------------------------------------------

    def compute_mass_matrix(self, lumping: str = "consistent") -> np.ndarray:
        """Masa del Hex8 por cuadratura del propio elemento.

        **Consistente** (default): ``M_e = ∫ρ·N^T·N dV``. Para Hex8 trilineal
        con cuadratura 2×2×2 la integral es exacta.

        **Lumped**: HRZ canónico (ADR 0009 fase 2). Para Hex8 regular coincide
        con row-sum por simetría — ``ρ·V_e/8`` en cada DOF traslacional.

        Returns
        -------
        np.ndarray, shape (24, 24)
        """
        rho = self.material.density
        coords = self.get_coordinate_matrix(ndim=3)
        M_s = np.zeros((8, 8))
        m_total = 0.0
        # Masa siempre con la regla completa 2×2×2 (exacta para el producto
        # trilineal×trilineal); con hex_1x1x1 M sería de rango 3.
        for (xi, eta, zeta), w in zip(*QuadratureRegistry.get("hex_2x2x2")):
            N = _shape_functions_hex8(xi, eta, zeta)
            detJ = _det_jacobian_hex8(xi, eta, zeta, coords)
            if detJ <= 0.0:
                raise ValueError(f"Jacobiano no positivo en Hex8 id={self.id}.")
            weight = detJ * w
            M_s += rho * np.outer(N, N) * weight
            m_total += rho * weight
        M_consistent = _expand_scalar_mass_3d(M_s)
        validate_lumping_kwarg(lumping, type(self).__name__)
        if lumping == "lumped":
            return lump_hrz(M_consistent, total_mass=m_total, n_translational_dirs=3)
        return M_consistent
