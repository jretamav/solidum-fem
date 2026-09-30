"""``Frame2DEulerCorot`` — viga 2D Euler-Bernoulli corotacional
(Updated Lagrangian, grandes desplazamientos y rotaciones rígidas).
"""
import math
from typing import List

import numpy as np

from solidum.core.element import Element
from solidum.core.material import Material
from solidum.core.node import Node
from solidum.core.section import Section
from solidum.elements.frame._section import DEFAULT_N_IP, FrameSectionMixin
from solidum.elements.frame._shared import (
    _frame2d_consistent_body_load,
    _frame2d_consistent_mass_local,
    _frame2d_forces_from_local,
    _frame2d_lumped_mass_local,
)
from solidum.registry import ElementRegistry, MaterialRegistry, SectionRegistry
from solidum.results import ElementForces


@ElementRegistry.register
class Frame2DEulerCorot(FrameSectionMixin, Element):
    """Marco/viga 2D Euler-Bernoulli corotacional (Updated Lagrangian).

    Dos nodos, 3 DOFs por nodo (ux, uy, rz). Captura grandes desplazamientos
    y grandes rotaciones rígidas del elemento con pequeñas deformaciones
    axiales. Rotaciones deformacionales nodales moderadas (|θ̄| ≲ 30°).

    Implementación autónoma: hereda directamente de Element, sin relación
    con Frame2DEuler ni con ningún otro elemento. La construcción de la
    cinemática corotacional se hace íntegramente dentro de la clase.

    La ley del elemento vive en el **sistema básico**
    ``v = [Δl, θ̄₁, θ̄₂]`` con fuerzas básicas ``q = [N, M₁, M₂]`` y su
    tangente ``K_ll``; la cinemática corotacional las transforma a globales
    y añade la rigidez geométrica con ``N``, ``M₁``, ``M₂``. Dos caminos
    excluyentes para la ley básica (ADR 0021):

    - ``material`` + ``A``, ``I``: forma cerrada elástica
      (``N = σA``, ``M = (EI/L₀)·[[4, 2], [2, 4]]·θ̄``). Sólo materiales
      elásticos lineales.
    - ``section`` (+ ``n_ip``): la sección se integra en ``n_ip`` puntos de
      Gauss-Legendre con ``ε₀ = Δl/L₀`` y
      ``κ_z(x) = [(6x/L₀ − 4)·θ̄₁ + (6x/L₀ − 2)·θ̄₂]/L₀`` (spec
      ``FiberSection`` §5).

    Parameters
    ----------
    element_id : int
    nodes : List[Node]
        Exactamente 2 nodos extremos.
    material : Material, optional
        Material 1D elástico lineal (``STRAIN_DIM = 1``).
    A : float, optional
        Área de la sección transversal (camino con material).
    I : float, optional
        Momento de inercia respecto al eje perpendicular al plano (camino
        con material).
    section : Section, optional
        Sección de marco; excluyente con ``material``, ``A`` e ``I``.
    n_ip : int, default 3
        Puntos de Gauss-Legendre (≥ 2); sólo con ``section``.
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
            raise ValueError("El elemento Frame2DEulerCorot requiere exactamente 2 nodos.")

        self._setup_constitution(element_id, material, section, n_ip,
                                 geometry={"A": A, "I": I})
        self.A = A
        self.I = I
        super().__init__(element_id, nodes, material)

        c1 = self.nodes[0].coordinates
        c2 = self.nodes[1].coordinates
        dx0, dy0 = c2[0] - c1[0], c2[1] - c1[1]
        self.L0 = math.sqrt(dx0**2 + dy0**2)
        if self.L0 == 0.0:
            raise ValueError("La longitud del elemento no puede ser cero.")
        self.alpha0 = math.atan2(dy0, dx0)
        if self.section is not None:
            self._B_basic = self._build_basic_B()
            self._wL = self._gauss_w * self.L0

    def _build_basic_B(self) -> np.ndarray:
        """``B_g`` (2×3) del sistema básico: ``[ε₀, κ_z(x_g)] = B_g·[Δl, θ̄₁, θ̄₂]``."""
        L0 = self.L0
        B = np.zeros((self.n_ip, 2, 3))
        for g, xi in enumerate(self._gauss_xi):
            B[g, 0, 0] = 1.0 / L0
            B[g, 1, 1] = (6.0 * xi - 4.0) / L0
            B[g, 1, 2] = (6.0 * xi - 2.0) / L0
        return B

    @staticmethod
    def _unwrap(angle: float) -> float:
        """Lleva un ángulo al intervalo (-π, π]."""
        return math.atan2(math.sin(angle), math.cos(angle))

    def _current_geometry(self, u_e):
        c1 = self.nodes[0].coordinates
        c2 = self.nodes[1].coordinates
        x1 = c1[0] + u_e[0]; y1 = c1[1] + u_e[1]
        x2 = c2[0] + u_e[3]; y2 = c2[1] + u_e[4]
        dx, dy = x2 - x1, y2 - y1
        l = math.sqrt(dx**2 + dy**2)
        if l == 0.0:
            raise ValueError("Longitud corriente nula en Frame2DEulerCorot.")
        alpha = math.atan2(dy, dx)
        alpha_e = self._unwrap(alpha - self.alpha0)
        return l, dx / l, dy / l, alpha_e

    def _basic_deformations(self, u_e):
        """Geometría corriente y deformaciones básicas ``v = [Δl, θ̄₁, θ̄₂]``."""
        l, c, s, alpha_e = self._current_geometry(u_e)
        u_elong = l - self.L0
        theta_bar_1 = self._unwrap(u_e[2] - alpha_e)
        theta_bar_2 = self._unwrap(u_e[5] - alpha_e)
        return l, c, s, u_elong, theta_bar_1, theta_bar_2

    def _basic_response(self, u_elong, theta_bar_1, theta_bar_2):
        """Fuerzas básicas ``q = [N, M₁, M₂]`` y su tangente ``K_ll`` (3×3);
        actualiza el trial."""
        if self.section is not None:
            v = np.array([u_elong, theta_bar_1, theta_bar_2])
            K_ll, F_l = self._integrate_section(self._B_basic, self._wL, v)
            return F_l, K_ll

        epsilon = u_elong / self.L0
        sigma, E_t, new_state = self.material.compute_state(epsilon, self.state.vars[0])
        self.state.vars_trial[0] = new_state
        self.state.stresses_trial[0] = sigma

        EA_L0 = E_t * self.A / self.L0
        EI_L0 = E_t * self.I / self.L0
        K_ll = np.array([
            [EA_L0,        0,        0],
            [    0, 4*EI_L0, 2*EI_L0],
            [    0, 2*EI_L0, 4*EI_L0],
        ])

        N  = sigma * self.A
        M1 = K_ll[1, 1] * theta_bar_1 + K_ll[1, 2] * theta_bar_2
        M2 = K_ll[2, 1] * theta_bar_1 + K_ll[2, 2] * theta_bar_2
        return np.array([N, M1, M2]), K_ll

    def compute_element_state(self, u_e):
        l, c, s, u_elong, theta_bar_1, theta_bar_2 = self._basic_deformations(u_e)
        F_l, K_ll = self._basic_response(u_elong, theta_bar_1, theta_bar_2)
        N, M1, M2 = F_l[0], F_l[1], F_l[2]

        B = np.array([
            [  -c,   -s, 0,    c,    s, 0],
            [-s/l,  c/l, 1,  s/l, -c/l, 0],
            [-s/l,  c/l, 0,  s/l, -c/l, 1],
        ])

        F_int_e = B.T @ F_l

        K_material = B.T @ K_ll @ B

        r = np.array([  -c,   -s, 0,    c,    s, 0])
        z = np.array([  -s,    c, 0,    s,   -c, 0])
        K_sigma = (N / l) * np.outer(z, z) \
                - ((M1 + M2) / (l**2)) * (np.outer(r, z) + np.outer(z, r))

        K_T = K_material + K_sigma
        return K_T, F_int_e

    def compute_internal_forces(self, U_global):
        """Diccionario legado de post-proceso (no es el contrato del ADR 0002).

        Fuerzas de extremo proyectadas en el sistema corotado. Con material,
        además ``stress`` y ``strain`` axiales; con sección, ``strain``
        (``ε₀ = Δl/L₀``), ``section_deformations`` y ``section_forces``
        ``(n_ip, 2)`` en los puntos de Gauss.
        """
        u_e = self.get_local_displacements(U_global)
        l, c, s, u_elong, theta_bar_1, theta_bar_2 = self._basic_deformations(u_e)
        epsilon = u_elong / self.L0

        _, F_int = self._evaluate_without_trace(u_e)
        # Proyección a sistema local corotacional para separar axial/cortante
        lam2 = np.array([[c, s], [-s, c]])
        F1 = lam2 @ F_int[0:2]
        out = {
            'axial_force': F1[0],
            'shear_force': F1[1],
            'moment_i': F_int[2],
            'moment_j': F_int[5],
            'strain': epsilon,
        }
        if self.section is None:
            sigma, _, _ = self.material.compute_state(epsilon, self.state.vars[0])
            out['stress'] = sigma
        else:
            v = np.array([u_elong, theta_bar_1, theta_bar_2])
            e, q = self._section_forces_at(self._B_basic, v)
            out['section_deformations'] = e
            out['section_forces'] = q
        return out

    def internal_forces(self, U_global: np.ndarray,
                        equivalent_load: np.ndarray | None = None) -> ElementForces:
        """API pública (ADR 0002): N, V, M en nodos i, j, convención §5.

        Usa los cosenos directores de la configuración corriente (corotacional);
        el ``state`` interno ya está comprometido tras ``solve()``.
        """
        u_e = self.get_local_displacements(U_global)
        _, c, s, _ = self._current_geometry(u_e)
        _, F_int = self._evaluate_without_trace(u_e)
        if equivalent_load is not None:
            # Fuerzas internas de extremo = fuerzas nodales − carga equivalente.
            F_int = F_int - np.asarray(equivalent_load, dtype=float)
        lam2 = np.array([[c, s], [-s, c]])
        F_local = np.empty(6)
        F_local[0:2] = lam2 @ F_int[0:2]
        F_local[2]   = F_int[2]
        F_local[3:5] = lam2 @ F_int[3:5]
        F_local[5]   = F_int[5]
        return _frame2d_forces_from_local(F_local)

    def compute_body_load(self, b: np.ndarray) -> np.ndarray:
        """Vector nodal consistente con una fuerza de cuerpo uniforme ``b``,
        evaluado en la configuración de referencia (carga por unidad de
        longitud ``A·b``; ``A = Σ A_i`` con sección).

        Para grandes rotaciones con cargas conservadoras (gravedad), la
        diferencia frente a la integración sobre la geometría corriente es
        de segundo orden y se acepta como aproximación estándar. Ver
        :func:`_frame2d_consistent_body_load`.
        """
        T = self._reference_transform()
        return _frame2d_consistent_body_load(b, self._axial_area(), self.L0, T)

    def compute_self_weight(self, g: np.ndarray) -> np.ndarray:
        """Peso propio en la configuración de referencia: ``ρA·g`` por unidad
        de longitud (``ρA = Σ ρ_i A_i`` con sección, spec §6)."""
        if self.section is None:
            return super().compute_self_weight(g)
        rhoA = self.section.mass_properties()["rhoA"]
        return _frame2d_consistent_body_load(g, rhoA, self.L0, self._reference_transform())

    def _reference_transform(self) -> np.ndarray:
        """Matriz de transformación 6×6 en la configuración de referencia."""
        c, s = math.cos(self.alpha0), math.sin(self.alpha0)
        return np.array([
            [ c,  s, 0, 0, 0, 0],
            [-s,  c, 0, 0, 0, 0],
            [ 0,  0, 1, 0, 0, 0],
            [ 0,  0, 0, c, s, 0],
            [ 0,  0, 0,-s, c, 0],
            [ 0,  0, 0, 0, 0, 1],
        ])

    def _mass_arguments(self):
        """``(ρ, A, I)`` para las funciones de masa de ``_shared``: con sección,
        ``(1, ρA, ρI_z)`` de la sección (spec §6)."""
        if self.section is None:
            return self.material.density, self.A, self.I
        mp = self.section.mass_properties()
        return 1.0, mp["rhoA"], mp["rhoIz"]

    def compute_mass_matrix(self, lumping: str = "consistent") -> np.ndarray:
        """Matriz de masa del Frame2DEulerCorot en ejes globales (ADR 0009).

        Evaluada en la configuración de referencia (Lagrangeano total).
        Análoga a Frame2DEuler pero la matriz de transformación se
        reconstruye desde ``alpha0`` ya que el corotacional no almacena
        ``self.T`` precalculada. Con sección, ``ρA`` y ``ρI`` son los de
        ``section.mass_properties()``.

        **Consistente**: :func:`_frame2d_consistent_mass_local`.

        **Lumped** (fase 2): :func:`_frame2d_lumped_mass_local`.
        """
        if lumping not in ("lumped", "consistent"):
            raise NotImplementedError(
                f"Frame2DEulerCorot.compute_mass_matrix: lumping='{lumping}' "
                f"no soportado. Valores admitidos: 'consistent', 'lumped'."
            )
        rho, A, I = self._mass_arguments()
        if lumping == "lumped":
            M_local = _frame2d_lumped_mass_local(rho, A, I, self.L0)
        else:
            M_local = _frame2d_consistent_mass_local(rho, A, I, self.L0)
        T = self._reference_transform()
        return T.T @ M_local @ T
