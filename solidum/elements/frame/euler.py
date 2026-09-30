"""``Frame2DEuler`` — viga 2D basada en Euler-Bernoulli (vigas esbeltas).
"""
from typing import List

import numpy as np

from solidum.core.element import Element
from solidum.core.material import Material
from solidum.core.node import Node
from solidum.core.section import Section
from solidum.elements.frame._section import (
    DEFAULT_N_IP,
    FrameSectionMixin,
    hermite_second_derivatives,
)
from solidum.elements.frame._shared import (
    _frame2d_consistent_body_load,
    _frame2d_consistent_mass_local,
    _frame2d_forces_from_local,
    _frame2d_lumped_mass_local,
    build_geometry_2d,
)
from solidum.registry import ElementRegistry, MaterialRegistry, SectionRegistry
from solidum.results import ElementForces


@ElementRegistry.register
class Frame2DEuler(FrameSectionMixin, Element):
    """Marco/viga 2D basado en Euler-Bernoulli.

    Dos nodos rígidamente conectados, 3 DOFs por nodo (ux, uy, rz). Transmite
    fuerza axial, cortante y momento flector. Régimen de linealidad
    geométrica y régimen de validez Euler-Bernoulli (vigas esbeltas,
    L/h ≳ 10).

    Dos caminos excluyentes para la constitutiva (ADR 0021):

    - ``material`` + ``A``, ``I``: forma cerrada elástica. Sólo materiales
      elásticos lineales (``IS_LINEAR_ELASTIC = True``).
    - ``section`` (+ ``n_ip``): elemento en desplazamientos que integra la
      sección (p. ej. :class:`~solidum.sections.FiberSection`) en ``n_ip``
      puntos de Gauss-Legendre (spec ``FiberSection`` §5). Axial lineal y
      transversal cúbica de Hermite: ``ε₀ = (u₂ − u₁)/L`` constante y
      ``κ_z(x) = N₁''v₁ + N₂''θ₁ + N₃''v₂ + N₄''θ₂`` lineal.

    Parameters
    ----------
    element_id : int
    nodes : List[Node]
        Exactamente 2 nodos extremos.
    material : Material, optional
        Material 1D elástico lineal (``STRAIN_DIM = 1``,
        ``IS_LINEAR_ELASTIC = True``).
    A : float, optional
        Área de la sección transversal (camino con material).
    I : float, optional
        Momento de inercia respecto al eje perpendicular al plano (z)
        (camino con material).
    section : Section, optional
        Sección de marco; excluyente con ``material``, ``A`` e ``I``.
    n_ip : int, default 3
        Puntos de Gauss-Legendre a lo largo del elemento (≥ 2); sólo con
        ``section``.
    """

    DOF_NAMES = ['ux', 'uy', 'rz']
    STRAIN_DIM = 1
    N_INTEGRATION_POINTS = 1
    REFERENCE_KWARGS = {"material": MaterialRegistry, "section": SectionRegistry}
    SECTION_DIM = 2
    SECTION_RESPONSE = "response_2d"

    def __init__(self, element_id: int, nodes: List[Node], material: Material | None = None,
                 A: float | None = None, I: float | None = None, *,
                 section: Section | None = None, n_ip: int = DEFAULT_N_IP):
        if len(nodes) != 2:
            raise ValueError("El elemento Frame2DEuler requiere exactamente 2 nodos.")

        self._setup_constitution(element_id, material, section, n_ip,
                                 geometry={"A": A, "I": I})
        self.A = A
        self.I = I
        super().__init__(element_id, nodes, material)

        self.L0, self.c, self.s, self.T = build_geometry_2d(self.nodes)
        if self.section is not None:
            self._B = self._build_section_B()
            self._wL = self._gauss_w * self.L0

    def _build_section_B(self) -> np.ndarray:
        """``B_g`` (2×6) en ejes locales: ``e_g = [ε₀, κ_z(x_g)] = B_g·u_local``."""
        L = self.L0
        B = np.zeros((self.n_ip, 2, 6))
        for g, xi in enumerate(self._gauss_xi):
            d2 = hermite_second_derivatives(xi, L)
            B[g, 0, 0] = -1.0 / L
            B[g, 0, 3] = 1.0 / L
            B[g, 1, [1, 2, 4, 5]] = d2
        return B

    def compute_element_state(self, u_e: np.ndarray):
        if self.section is not None:
            return self._section_element_state(u_e)
        L = self.L0
        u_local = self.T @ u_e

        epsilon = (u_local[3] - u_local[0]) / L
        sigma, E_t, new_state = self.material.compute_state(epsilon, self.state.vars[0])
        self.state.vars_trial[0] = new_state
        self.state.stresses_trial[0] = sigma

        EA_L = E_t * self.A / L
        EI_L = E_t * self.I / L
        EI_L2 = EI_L / L
        EI_L3 = EI_L2 / L

        K_local = np.array([
            [ EA_L,         0,         0, -EA_L,         0,         0],
            [    0,  12*EI_L3,   6*EI_L2,     0, -12*EI_L3,   6*EI_L2],
            [    0,   6*EI_L2,    4*EI_L,     0,  -6*EI_L2,    2*EI_L],
            [-EA_L,         0,         0,  EA_L,         0,         0],
            [    0, -12*EI_L3,  -6*EI_L2,     0,  12*EI_L3,  -6*EI_L2],
            [    0,   6*EI_L2,    2*EI_L,     0,  -6*EI_L2,    4*EI_L],
        ])

        F_int_local = K_local @ u_local
        # Componente axial con el σ del material. Con material elástico lineal
        # (única opción admitida en este camino, ADR 0021) coincide con la
        # fila axial de ``K_local @ u_local``.
        F_int_local[0] = -sigma * self.A
        F_int_local[3] =  sigma * self.A

        K_global = self.T.T @ K_local @ self.T
        F_int_e = self.T.T @ F_int_local
        return K_global, F_int_e

    def _section_element_state(self, u_e: np.ndarray):
        """Camino con sección: ``K_local = Σ w_g L B_gᵀ k_s B_g``,
        ``F_local = Σ w_g L B_gᵀ s_g`` (spec §5)."""
        u_local = self.T @ u_e
        K_local, F_int_local = self._integrate_section(self._B, self._wL, u_local)
        return self.T.T @ K_local @ self.T, self.T.T @ F_int_local

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        """Diccionario legado de post-proceso (no es el contrato del ADR 0002).

        Con material: fuerzas de extremo, ``stress`` y ``strain`` axiales del
        material. Con sección: fuerzas de extremo, ``strain`` (``ε₀`` del
        eje), ``section_deformations`` ``(n_ip, 2)`` = ``[ε₀, κ_z]`` y
        ``section_forces`` ``(n_ip, 2)`` = ``[N, M_z]`` en los puntos de Gauss.
        """
        u_e = self.get_local_displacements(U_global)
        _, F_int = self._evaluate_without_trace(u_e)
        # ``compute_element_state`` aplica T.T para rotar F_int_local → F_int (globales).
        # Aquí lo rotamos de vuelta a locales para extraer N/V/M. La doble rotación
        # ``self.T @ (self.T.T @ F_int_local)`` es identidad porque T es ortogonal
        # (rotación pura) — pagamos un producto 6×6·6 cosmético por mantener el
        # contrato uniforme de ``compute_element_state`` devolviendo en globales
        # (auditoría H-2.8, severidad bajo).
        F_local = self.T @ F_int

        u_local = self.T @ u_e
        epsilon = (u_local[3] - u_local[0]) / self.L0
        out = {
            'axial_force': F_local[0],
            'shear_force': F_local[1],
            'moment_i': F_local[2],
            'moment_j': F_local[5],
            'strain': epsilon,
        }
        if self.section is None:
            sigma, _, _ = self.material.compute_state(epsilon, self.state.vars[0])
            out['stress'] = sigma
        else:
            e, s = self._section_forces_at(self._B, u_local)
            out['section_deformations'] = e
            out['section_forces'] = s
        return out

    def internal_forces(self, U_global: np.ndarray,
                        equivalent_load: np.ndarray | None = None) -> ElementForces:
        """API pública (ADR 0002): N, V, M en nodos i, j, convención §5."""
        u_e = self.get_local_displacements(U_global)
        _, F_int = self._evaluate_without_trace(u_e)
        if equivalent_load is not None:
            # Fuerzas internas de extremo = fuerzas nodales − carga equivalente.
            F_int = F_int - np.asarray(equivalent_load, dtype=float)
        F_local = self.T @ F_int
        return _frame2d_forces_from_local(F_local)

    def compute_body_load(self, b: np.ndarray) -> np.ndarray:
        """Vector nodal consistente con una fuerza de cuerpo uniforme ``b``
        (por unidad de volumen): carga por unidad de longitud ``A·b``, con
        ``A`` el área bruta (``Σ A_i`` con sección). Ver
        :func:`_frame2d_consistent_body_load`.
        """
        return _frame2d_consistent_body_load(b, self._axial_area(), self.L0, self.T)

    def compute_self_weight(self, g: np.ndarray) -> np.ndarray:
        """Peso propio: carga por unidad de longitud ``ρA·g``. Con sección,
        ``ρA = Σ ρ_i A_i`` (spec §6); con material, ``ρ·A``."""
        if self.section is None:
            return super().compute_self_weight(g)
        rhoA = self.section.mass_properties()["rhoA"]
        return _frame2d_consistent_body_load(g, rhoA, self.L0, self.T)

    def _mass_arguments(self):
        """``(ρ, A, I)`` para las funciones de masa de ``_shared``: con sección,
        ``(1, ρA, ρI_z)`` de la sección (spec §6)."""
        if self.section is None:
            return self.material.density, self.A, self.I
        mp = self.section.mass_properties()
        return 1.0, mp["rhoA"], mp["rhoIz"]

    def compute_mass_matrix(self, lumping: str = "consistent") -> np.ndarray:
        """Matriz de masa del Frame2DEuler en ejes globales (ADR 0009).

        **Consistente** (default, fase 1): axial lineal + flexional
        Hermitiana cúbica + inercia rotacional de sección ``ρI``. Ver
        :func:`_frame2d_consistent_mass_local`.

        **Lumped** (fase 2): nodal directo ``ρAL/2`` traslacional + ``ρIL/2``
        rotacional por nodo. Ver :func:`_frame2d_lumped_mass_local`.

        Con sección, ``ρA`` y ``ρI`` son los de ``section.mass_properties()``.

        Ambos casos pasan por ``T^T · M_local · T``. Para el lumped, T es
        ortogonal por bloques compatibles con la diagonalidad: el resultado
        global también es diagonal.

        Returns
        -------
        np.ndarray, shape (6, 6)
        """
        if lumping not in ("lumped", "consistent"):
            raise NotImplementedError(
                f"Frame2DEuler.compute_mass_matrix: lumping='{lumping}' no "
                f"soportado. Valores admitidos: 'consistent', 'lumped'."
            )
        rho, A, I = self._mass_arguments()
        if lumping == "lumped":
            M_local = _frame2d_lumped_mass_local(rho, A, I, self.L0)
        else:
            M_local = _frame2d_consistent_mass_local(rho, A, I, self.L0)
        return self.T.T @ M_local @ self.T
