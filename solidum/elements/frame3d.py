# solidum_fem/solidum/elements/frame3d.py
"""Elemento de marco/viga 3D Euler-Bernoulli.

Implementación **independiente** del resto de elementos: no hereda de
ningún otro elemento; la construcción de ejes locales, la matriz de
transformación 12×12 y la rigidez local se generan íntegramente dentro de
la clase. Comparte con los marcos planos la constitutiva genérica de
``solidum.elements.frame._section`` (ADR 0021: forma cerrada con material
elástico lineal o integración de una sección) y los bloques de masa de la
flexión de ``solidum.elements.frame._shared`` (masa hermítica y masa
rotacional HRZ). Ver docs/specs/Frame3D.md.
"""
from __future__ import annotations

import math
from typing import List

import numpy as np

from solidum.core.element import Element
from solidum.core.material import Material
from solidum.core.node import Node
from solidum.core.section import Section
from solidum.elements._line import LineElementMixin, hermite_functions, linear_tangent
from solidum.elements.frame._section import (
    DEFAULT_N_IP,
    FrameSectionMixin,
    hermite_second_derivatives,
)
from solidum.elements.frame._shared import hermite_consistent_mass, hrz_rotational_mass
from solidum.logging import get_logger
from solidum.registry import ElementRegistry, MaterialRegistry, SectionRegistry
from solidum.results import ElementForces

_log = get_logger("elements.frame3d")


@ElementRegistry.register
class Frame3D(LineElementMixin, FrameSectionMixin, Element):
    """Marco/viga 3D basado en Euler-Bernoulli.

    Dos nodos rígidamente conectados, 6 DOFs por nodo (ux, uy, uz, rx, ry, rz).
    Transmite axial + cortante y flexión en dos planos + torsión. Régimen de
    linealidad geométrica.

    Dos caminos excluyentes para la constitutiva (ADR 0021):

    - ``material`` + ``A``, ``Iy``, ``Iz``, ``J`` (+ ``nu``): forma cerrada
      elástica. Sólo materiales elásticos lineales.
    - ``section`` (+ ``n_ip``): elemento en desplazamientos que integra la
      sección en ``n_ip`` puntos de Gauss-Legendre (spec ``FiberSection``
      §5) con ``e = [ε₀, κ_y, κ_z]``: ``ε₀ = (u₂ − u₁)/L``,
      ``κ_z = N₁''v₁ + N₂''θ_z1 + N₃''v₂ + N₄''θ_z2`` y, como ``w' = −θ_y``
      por la regla de la mano derecha,
      ``κ_y = −w'' = −N₁''w₁ + N₂''θ_y1 − N₃''w₂ + N₄''θ_y2``. Torsión
      elástica desacoplada ``T = GJ·θₓ'`` con el ``GJ`` de la sección
      (obligatorio).

    Parameters
    ----------
    element_id : int
    nodes : List[Node]
        Exactamente 2 nodos extremos. Acepta nodos con 3 coordenadas.
    material : Material, optional
        Material 1D elástico lineal (``STRAIN_DIM = 1``). Si expone `nu`, se
        usa ese valor.
    A : float, optional
        Área de la sección transversal (camino con material).
    Iy, Iz : float, optional
        Momentos de inercia respecto a los ejes locales y, z de la sección
        (camino con material).
    J : float, optional
        Constante torsional de Saint-Venant (camino con material).
    nu : float, optional
        Coeficiente de Poisson si el material no lo expone. Default 0.3
        (camino con material).
    ref_vector : list[float] | None, optional
        Vector de referencia 3D para fijar la orientación de los ejes
        locales y, z. Default `[0, 0, 1]`; con fallback `[1, 0, 0]` si
        la barra es cercanamente vertical.
    section : Section, optional
        Sección de marco con ``GJ`` declarado; excluyente con ``material``,
        ``A``, ``Iy``, ``Iz``, ``J`` y ``nu``.
    n_ip : int, default 3
        Puntos de Gauss-Legendre (≥ 2); sólo con ``section``.
    """

    DOF_NAMES = ['ux', 'uy', 'uz', 'rx', 'ry', 'rz']
    STRAIN_DIM = 1
    N_INTEGRATION_POINTS = 1
    # Campo interior en ejes locales (ADR 0022, fase 4).
    LINE_FIELDS = ("u", "v", "w", "theta_x", "theta_y", "theta_z")
    LINE_STRAINS = ("eps", "twist", "kappa_y", "kappa_z")
    REFERENCE_KWARGS = {"material": MaterialRegistry, "section": SectionRegistry}
    SECTION_DIM = 3
    SECTION_RESPONSE = "response_3d"

    _VERTICAL_COS_TOL = 0.99  # umbral para detectar barras casi-verticales
    # Flag de clase para emitir el warning de fallback ref_vector sólo una vez
    # por sesión (auditoría H-2.9): en estructuras con muchas columnas
    # verticales el aviso por cada elemento es ruido sin valor — la primera
    # emisión basta para alertar al usuario.
    _vertical_warning_emitted: bool = False

    def __init__(self, element_id: int, nodes: List[Node], material: Material | None = None,
                 A: float | None = None, Iy: float | None = None,
                 Iz: float | None = None, J: float | None = None,
                 nu: float | None = None,
                 ref_vector: List[float] | None = None, *,
                 section: Section | None = None, n_ip: int = DEFAULT_N_IP):
        if len(nodes) != 2:
            raise ValueError("El elemento Frame3D requiere exactamente 2 nodos.")

        self._setup_constitution(element_id, material, section, n_ip,
                                 geometry={"A": A, "Iy": Iy, "Iz": Iz, "J": J},
                                 material_only={"nu": nu})
        if self.section is not None and self.section.GJ is None:
            raise ValueError(
                f"Frame3D(id={element_id}): la sección no declara GJ. Frame3D "
                f"necesita la rigidez torsional elástica de la sección (torsión "
                f"desacoplada T = GJ·θx'); declárela con GJ en la sección."
            )
        self.A = A
        self.Iy = Iy
        self.Iz = Iz
        self.J = J
        super().__init__(element_id, nodes, material)

        if self.section is None:
            # Poisson: material.nu si existe, si no el parámetro. `None`
            # distingue "no declarado" de un 0.3 declarado a propósito: antes
            # se comparaba contra 0.3 y el aviso salía también cuando el
            # usuario lo escribía.
            if hasattr(material, 'nu'):
                self.nu = material.nu
            else:
                if nu is None:
                    nu = 0.3
                    _log.warning(
                        f"Frame3D (id={element_id}): el material no expone 'nu'. "
                        f"Se usará nu={nu} (default). Especifique 'nu' en el YAML si esto es incorrecto."
                    )
                self.nu = float(nu)
        else:
            # La torsión la da el GJ de la sección: no hace falta Poisson.
            self.nu = None

        self.L0, self.lam, self.T = self._build_local_frame(self.nodes, ref_vector,
                                                             element_id)
        if self.section is not None:
            self._B = self._build_section_B()
            self._wL = self._gauss_w * self.L0

    def _build_section_B(self) -> np.ndarray:
        """``B_g`` (3×12) en ejes locales: ``e_g = [ε₀, κ_y, κ_z](x_g) = B_g·u_local``.

        Layout local ``[u₁, v₁, w₁, θx₁, θy₁, θz₁, u₂, v₂, w₂, θx₂, θy₂, θz₂]``.
        """
        L = self.L0
        B = np.zeros((self.n_ip, 3, 12))
        for g, xi in enumerate(self._gauss_xi):
            d1, d2, d3, d4 = hermite_second_derivatives(xi, L)
            B[g, 0, 0] = -1.0 / L
            B[g, 0, 6] = 1.0 / L
            # κ_y = −w'' (plano xz, w' = −θ_y).
            B[g, 1, 2] = -d1
            B[g, 1, 4] = d2
            B[g, 1, 8] = -d3
            B[g, 1, 10] = d4
            # κ_z = v'' (plano xy, v' = θ_z).
            B[g, 2, 1] = d1
            B[g, 2, 5] = d2
            B[g, 2, 7] = d3
            B[g, 2, 11] = d4
        return B

    @classmethod
    def _build_local_frame(cls, nodes: List[Node],
                            ref_vector: List[float] | None,
                            element_id: int):
        """Construye longitud, matriz λ (3×3) y T (12×12). Una sola vez."""
        c1 = nodes[0].coordinates
        c2 = nodes[1].coordinates
        z1 = c1[2] if len(c1) > 2 else 0.0
        z2 = c2[2] if len(c2) > 2 else 0.0
        dx = c2[0] - c1[0]
        dy = c2[1] - c1[1]
        dz = z2 - z1
        L0 = math.sqrt(dx**2 + dy**2 + dz**2)
        if L0 == 0.0:
            raise ValueError("La longitud del elemento no puede ser cero.")

        x_local = np.array([dx, dy, dz]) / L0

        # Vector de referencia: usuario > default con fallback
        if ref_vector is None:
            if abs(x_local[2]) > cls._VERTICAL_COS_TOL:
                v_ref = np.array([1.0, 0.0, 0.0])
                # Aviso emitido una sola vez por sesión (H-2.9): estructuras
                # con N columnas verticales sin ref_vector explícito generaban
                # N warnings idénticos.
                if not cls._vertical_warning_emitted:
                    _log.warning(
                        f"Frame3D (id={element_id}): barra cercanamente vertical sin ref_vector explícito. "
                        f"Se usa fallback [1, 0, 0]. Especifique ref_vector si la orientación de la sección "
                        f"debe ser otra. (Este aviso se emite una vez por sesión; otras barras verticales "
                        f"con el mismo fallback no lo repetirán.)"
                    )
                    cls._vertical_warning_emitted = True
            else:
                v_ref = np.array([0.0, 0.0, 1.0])
        else:
            v_ref = np.asarray(ref_vector, dtype=float)
            if v_ref.shape != (3,):
                raise ValueError(
                    f"Frame3D (id={element_id}): ref_vector debe ser una lista "
                    f"de 3 componentes; se recibió shape {v_ref.shape}."
                )

        # Proyección de v_ref sobre el plano perpendicular a x_local
        y_tilde = v_ref - np.dot(v_ref, x_local) * x_local
        ny = np.linalg.norm(y_tilde)
        if ny < 1e-10:
            raise ValueError(
                f"Frame3D (id={element_id}): ref_vector es paralelo al eje de "
                f"la barra; no puede definir la orientación de la sección. "
                f"Elija otro vector."
            )
        y_local = y_tilde / ny
        z_local = np.cross(x_local, y_local)

        lam = np.array([x_local, y_local, z_local])  # 3×3
        T = np.zeros((12, 12))
        for block in range(4):
            i = 3 * block
            T[i:i+3, i:i+3] = lam

        return L0, lam, T

    def _build_local_stiffness(self, E_t: float, G: float) -> np.ndarray:
        """Matriz K_local 12×12 en el sistema de ejes de la barra."""
        L = self.L0
        EA_L = E_t * self.A / L
        GJ_L = G * self.J / L

        EIz = E_t * self.Iz
        EIy = E_t * self.Iy
        az = 12.0 * EIz / L**3
        bz = 6.0 * EIz / L**2
        cz = 4.0 * EIz / L
        dz = 2.0 * EIz / L
        ay = 12.0 * EIy / L**3
        by = 6.0 * EIy / L**2
        cy = 4.0 * EIy / L
        dy = 2.0 * EIy / L

        K = np.zeros((12, 12))

        # Axial (DOFs 0 ↔ 6)
        K[0, 0]  =  EA_L;  K[0, 6]  = -EA_L
        K[6, 0]  = -EA_L;  K[6, 6]  =  EA_L

        # Torsión (DOFs 3 ↔ 9)
        K[3, 3]  =  GJ_L;  K[3, 9]  = -GJ_L
        K[9, 3]  = -GJ_L;  K[9, 9]  =  GJ_L

        # Flexión en plano xy (DOFs 1, 5, 7, 11) con Iz
        K[1,  1]  =  az;  K[1,  5]  =  bz;  K[1,  7]  = -az;  K[1,  11] =  bz
        K[5,  1]  =  bz;  K[5,  5]  =  cz;  K[5,  7]  = -bz;  K[5,  11] =  dz
        K[7,  1]  = -az;  K[7,  5]  = -bz;  K[7,  7]  =  az;  K[7,  11] = -bz
        K[11, 1]  =  bz;  K[11, 5]  =  dz;  K[11, 7]  = -bz;  K[11, 11] =  cz

        # Flexión en plano xz (DOFs 2, 4, 8, 10) con Iy (signos invertidos
        # respecto a xy por orientación dextrógira uz ↔ -ry)
        K[2,  2]  =  ay;  K[2,  4]  = -by;  K[2,  8]  = -ay;  K[2,  10] = -by
        K[4,  2]  = -by;  K[4,  4]  =  cy;  K[4,  8]  =  by;  K[4,  10] =  dy
        K[8,  2]  = -ay;  K[8,  4]  =  by;  K[8,  8]  =  ay;  K[8,  10] =  by
        K[10, 2]  = -by;  K[10, 4]  =  dy;  K[10, 8]  =  by;  K[10, 10] =  cy

        return K

    def compute_element_state(self, u_e: np.ndarray):
        if self.section is not None:
            return self._section_element_state(u_e)
        u_local = self.T @ u_e
        L = self.L0

        epsilon = (u_local[6] - u_local[0]) / L
        sigma, E_t, new_state = self.material.compute_state(epsilon, self.state.vars[0])
        self.state.vars_trial[0] = new_state
        self.state.stresses_trial[0] = sigma

        G = E_t / (2.0 * (1.0 + self.nu))
        K_local = self._build_local_stiffness(E_t, G)

        F_int_local = K_local @ u_local
        # Componente axial con el σ del material. Con material elástico lineal
        # (única opción admitida en este camino, ADR 0021) coincide con la
        # fila axial de ``K_local @ u_local``.
        F_int_local[0] = -sigma * self.A
        F_int_local[6] =  sigma * self.A

        K_global = self.T.T @ K_local @ self.T
        F_int_e = self.T.T @ F_int_local
        return K_global, F_int_e

    def _section_element_state(self, u_e: np.ndarray):
        """Camino con sección (spec §5): flexión y axial integrados en los
        puntos de Gauss más la torsión elástica ``GJ/L`` en los DOFs θx."""
        u_local = self.T @ u_e
        K_local, F_int_local = self._integrate_section(self._B, self._wL, u_local)
        GJ_L = self.section.GJ / self.L0
        K_local[3, 3] += GJ_L
        K_local[3, 9] -= GJ_L
        K_local[9, 3] -= GJ_L
        K_local[9, 9] += GJ_L
        torque = GJ_L * (u_local[9] - u_local[3])     # T = GJ·θx'
        F_int_local[3] -= torque
        F_int_local[9] += torque
        return self.T.T @ K_local @ self.T, self.T.T @ F_int_local

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        """Diccionario legado de post-proceso (no es el contrato del ADR 0002).

        Con material, además ``stress`` y ``strain`` axiales del material;
        con sección, ``strain`` (``ε₀`` del eje), ``section_deformations``
        ``(n_ip, 3)`` = ``[ε₀, κ_y, κ_z]`` y ``section_forces`` ``(n_ip, 3)``
        = ``[N, M_y, M_z]`` en los puntos de Gauss.
        """
        u_e = self.get_local_displacements(U_global)
        _, F_int = self._evaluate_without_trace(u_e)
        F_local = self.T @ F_int

        u_local = self.T @ u_e
        epsilon = (u_local[6] - u_local[0]) / self.L0
        out = {
            'axial_force': F_local[0],
            'shear_y': F_local[1],
            'shear_z': F_local[2],
            'torsion': F_local[3],
            'moment_y_i': F_local[4],
            'moment_z_i': F_local[5],
            'moment_y_j': F_local[10],
            'moment_z_j': F_local[11],
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
        """API pública (ADR 0002): N, Vy, Vz, T, My, Mz en nodos i, j.

        Convención stress-resultant / RHR pura (Reglas.md §5). F_local tiene
        semántica ``K·u = F_ext`` en ejes locales, layout
        ``[Fi, Mi, Fj, Mj]`` (cada bloque 3-componente). El nodo i vive en la
        cara con normal saliente ``−x_local`` → signo invertido para todas las
        componentes; el nodo j en cara ``+x_local`` → signo directo.

        Verificado con: tracción pura (N=[+F,+F]), tip-load y (Vy=[-F,-F],
        Mz=[-F,0]), tip-load z (Vz=[-F,-F], My=[+F,0]), torsión pura
        (T=[+T,+T]).
        """
        u_e = self.get_local_displacements(U_global)
        _, F_int = self._evaluate_without_trace(u_e)
        if equivalent_load is not None:
            # Fuerzas internas de extremo = fuerzas nodales − carga equivalente.
            F_int = F_int - np.asarray(equivalent_load, dtype=float)
        F_local = self.T @ F_int
        return ElementForces(
            kind="frame3d",
            components={
                "N":  np.array([-F_local[0],  F_local[6]]),
                "Vy": np.array([-F_local[1],  F_local[7]]),
                "Vz": np.array([-F_local[2],  F_local[8]]),
                "T":  np.array([-F_local[3],  F_local[9]]),
                "My": np.array([-F_local[4],  F_local[10]]),
                "Mz": np.array([-F_local[5],  F_local[11]]),
            },
        )

    # Signos de la flexión en el plano xz: con (w, θ_y) y w' = −θ_y, la masa
    # hermítica del plano xy (v, θ_z con v' = θ_z) cambia de signo en los
    # términos que llevan exactamente un giro: H_xz = S·H_xy·S.
    _XZ_SIGNS = np.diag([1.0, -1.0, 1.0, -1.0])
    _IDX_XY = [1, 5, 7, 11]     # (uy_i, rz_i, uy_j, rz_j)
    _IDX_XZ = [2, 4, 8, 10]     # (uz_i, ry_i, uz_j, ry_j)

    def _build_local_mass(self, rhoA: float, rhoJp: float) -> np.ndarray:
        """Matriz de masa consistente 12×12 en ejes locales del elemento.

        Recibe las propiedades másicas por unidad de longitud ``ρA`` y
        ``ρJ_p`` (ver :meth:`_mass_properties`). Contribuciones desacopladas
        en el layout local
        ``[ux_i, uy_i, uz_i, rx_i, ry_i, rz_i, ux_j, uy_j, uz_j, rx_j, ry_j, rz_j]``:

        - **Axial** ``(ux_i, ux_j)``: ``(ρAL/6)·[[2,1],[1,2]]``.
        - **Torsión** ``(rx_i, rx_j)``: ``(ρJ_pL/6)·[[2,1],[1,2]]`` con
          ``J_p = I_y + I_z`` (momento polar geométrico de área, no la
          constante torsional de Saint-Venant ``J``, que rige la rigidez).
        - **Flexión xy** ``(uy_i, rz_i, uy_j, rz_j)``: masa hermítica
          ``ρAL/420·H`` (:func:`hermite_consistent_mass`) con ``v' = θ_z``.
        - **Flexión xz** ``(uz_i, ry_i, uz_j, ry_j)``: la misma con
          ``w' = −θ_y``: ``S·H·S`` con ``S = diag(1, −1, 1, −1)``, es decir,
          los términos de acoplamiento traslación-giro cambian de signo.
          Coherente con la rigidez ``K_local`` de este elemento.

        Sin inercia rotatoria de flexión: Euler-Bernoulli la desprecia por
        hipótesis (antes de 2026-09-29 se sumaba ``ρI·L/6·[[2,1],[1,2]]`` en
        ``ry`` y ``rz``, cinemática de Timoshenko).
        """
        L = self.L0
        c_ax = rhoA * L / 6.0
        c_tor = rhoJp * L / 6.0

        M = np.zeros((12, 12))

        # Axial (0 ↔ 6)
        M[0, 0]  = 2.0 * c_ax;  M[0, 6]  = 1.0 * c_ax
        M[6, 0]  = 1.0 * c_ax;  M[6, 6]  = 2.0 * c_ax

        # Torsión (3 ↔ 9)
        M[3, 3]  = 2.0 * c_tor; M[3, 9]  = 1.0 * c_tor
        M[9, 3]  = 1.0 * c_tor; M[9, 9]  = 2.0 * c_tor

        # Flexión en los dos planos.
        H_xy = hermite_consistent_mass(rhoA, L)
        S = self._XZ_SIGNS
        M[np.ix_(self._IDX_XY, self._IDX_XY)] = H_xy
        M[np.ix_(self._IDX_XZ, self._IDX_XZ)] = S @ H_xy @ S
        return M

    def _build_local_mass_lumped(self, rhoA: float, rhoJp: float) -> np.ndarray:
        """Matriz de masa concentrada 12×12 en ejes locales (ADR 0009).

        Mismas propiedades másicas por unidad de longitud que
        :meth:`_build_local_mass`. Por nodo:

        - Traslacional ``m_t = ρAL/2`` en ``ux``, ``uy``, ``uz``.
        - Torsional ``m_rx = ρJ_pL/2`` con ``J_p = I_y + I_z``.
        - Flexional HRZ ``m_ry = m_rz = ρAL³/78`` (:func:`hrz_rotational_mass`),
          sin inercia rotatoria de la sección (Euler-Bernoulli).

        Diagonal en locales. **En globales la diagonalidad sólo se conserva
        cuando el eje del elemento coincide con un eje global.** En
        orientación oblicua el bloque traslacional 3×3 por nodo sigue siendo
        diagonal (``m_t·I₃`` es invariante ante SO(3)), pero el rotacional
        ``diag(m_rx, m_θ, m_θ)`` no es escalar: ``m_rx = ρJ_pL/2`` y
        ``m_θ = ρAL³/78`` coincidirían sólo si ``J_p/A = L²/39``, imposible
        en una barra esbelta. Tras rotar queda lleno y ``M_global`` es
        **diagonal por bloques** por nodo: la limitación estándar de la masa
        concentrada en marcos 3D, aceptable para Newmark/HHT;
        ``CentralDifferenceSolver`` la rechaza.

        Returns
        -------
        np.ndarray, shape (12, 12)
        """
        L = self.L0
        m_t = rhoA * L / 2.0
        m_rx = rhoJp * L / 2.0
        m_th = hrz_rotational_mass(rhoA, L)
        return np.diag([m_t, m_t, m_t, m_rx, m_th, m_th,
                        m_t, m_t, m_t, m_rx, m_th, m_th])

    def compute_mass_matrix(self, lumping: str = "consistent") -> np.ndarray:
        """Matriz de masa del Frame3D en ejes globales (ADR 0009).

        **Consistente**: :meth:`_build_local_mass` (axial, torsión ``ρJ_p``
        y flexión hermítica en los dos planos, sin inercia rotatoria de
        flexión).

        **Concentrada** (``"lumped"``): :meth:`_build_local_mass_lumped`
        (HRZ en la flexión). Ver allí la nota sobre la diagonalidad por
        bloques en orientación oblicua.

        Con sección, las propiedades másicas salen de
        ``section.mass_properties()`` (spec ``FiberSection`` §6).

        Returns
        -------
        np.ndarray, shape (12, 12)
        """
        if lumping not in ("lumped", "consistent"):
            raise NotImplementedError(
                f"Frame3D.compute_mass_matrix: lumping='{lumping}' no "
                f"soportado. Valores admitidos: 'consistent', 'lumped'."
            )
        props = self._mass_properties()
        if lumping == "lumped":
            M_local = self._build_local_mass_lumped(*props)
        else:
            M_local = self._build_local_mass(*props)
        return self.T.T @ M_local @ self.T

    def _mass_properties(self):
        """``(ρA, ρJ_p)`` por unidad de longitud. Con material, ``ρ·A`` y
        ``ρ·(I_y + I_z)``; con sección, ``ρA`` y ``ρI_y + ρI_z`` de
        ``section.mass_properties()``. ``ρJ_p`` sólo entra en la torsión."""
        if self.section is None:
            rho = self.material.density
            return rho * self.A, rho * (self.Iy + self.Iz)
        mp = self.section.mass_properties()
        return mp["rhoA"], mp["rhoIy"] + mp["rhoIz"]

    def compute_body_load(self, b: np.ndarray) -> np.ndarray:
        """Vector nodal consistente con ∫NᵀbA dx para marco 3D (peso propio).

        La carga distribuida por unidad de longitud ``q = A · b`` (``A`` el
        área bruta: ``Σ A_i`` con sección) se reparte con
        :meth:`_line_load_vector` si es constante. Un campo ``b(x)`` (ADR
        0022) se proyecta a ejes locales en 4 puntos de Gauss y se reparte
        con las funciones de la rigidez (axial lineal, flexión de Hermite en
        los dos planos); centrado en el eje, no produce torsión.

        Parameters
        ----------
        b : np.ndarray, shape (3,), o callable
            Fuerza de cuerpo por unidad de volumen en ejes globales (e.g.
            ``b = (0, 0, -ρ·g)`` para peso propio con ``+z`` hacia arriba), o
            un campo ``b(x)`` de la posición.

        Returns
        -------
        np.ndarray, shape (12,)
            Vector consistente en ejes globales, layout estándar de Frame3D:
            ``[Fx_i, Fy_i, Fz_i, Mx_i, My_i, Mz_i, Fx_j, ..., Mz_j]``.
        """
        if callable(b):
            return self._field_line_load(b, self._axial_area())
        b = np.asarray(b, dtype=float).reshape(3)
        return self._line_load_vector(self._axial_area() * b)

    # ------------------------------------------------------------------
    # Campo interior (ADR 0022, fase 4)
    # ------------------------------------------------------------------

    def _to_local(self, u_e):
        return self.T @ u_e

    def _to_global(self, f_local):
        return self.T.T @ f_local

    def _load_axes(self):
        return self.lam

    def _line_interpolation(self, xi):
        """Axial y torsión lineales; flexión de Hermite en los dos planos:
        ``v = N·[v₁, θz₁, v₂, θz₂]`` (``θz = v'``) y, por la regla de la mano
        derecha, ``w = N·[w₁, −θy₁, w₂, −θy₂]`` (``θy = −w'``). Deformaciones
        ``[ε, θx', κ_y = θy' = −w'', κ_z = v'']``, las de
        :meth:`_build_section_B`."""
        L = self.L0
        n = len(xi)
        H, dH, d2H = hermite_functions(xi, L)
        sign = np.array([1.0, -1.0, 1.0, -1.0])
        N = np.zeros((n, 6, 12))
        B = np.zeros((n, 4, 12))
        N[:, 0, 0], N[:, 0, 6] = 1.0 - xi, xi
        N[:, 1, [1, 5, 7, 11]] = H
        N[:, 2, [2, 4, 8, 10]] = H * sign
        N[:, 3, 3], N[:, 3, 9] = 1.0 - xi, xi
        N[:, 4, [2, 4, 8, 10]] = -dH * sign
        N[:, 5, [1, 5, 7, 11]] = dH
        B[:, 0, 0], B[:, 0, 6] = -1.0 / L, 1.0 / L
        B[:, 1, 3], B[:, 1, 9] = -1.0 / L, 1.0 / L
        B[:, 2, [2, 4, 8, 10]] = -d2H * sign
        B[:, 3, [1, 5, 7, 11]] = d2H
        return N, B

    def section_rigidities(self) -> np.ndarray:
        """``[EA, GJ, EI_y, EI_z]`` con ``G = E/(2(1 + ν))``; sólo con
        material elástico lineal."""
        if self.section is not None:
            raise ValueError(f"Frame3D(id={self.id}): con sección por fibras las "
                             "rigideces de la sección no son constantes.")
        E = linear_tangent(self.material)
        G = E / (2.0 * (1.0 + self.nu))
        return np.array([E * self.A, G * self.J, E * self.Iy, E * self.Iz])

    def compute_self_weight(self, g: np.ndarray) -> np.ndarray:
        """Peso propio: carga por unidad de longitud ``ρA·g`` (``ρA = Σ ρ_i A_i``
        con sección, spec §6)."""
        if self.section is None:
            return super().compute_self_weight(g)
        g = np.asarray(g, dtype=float).reshape(3)
        return self._line_load_vector(self.section.mass_properties()["rhoA"] * g)

    def _line_load_vector(self, q_global: np.ndarray) -> np.ndarray:
        """Vector nodal consistente de una carga uniforme por unidad de
        longitud ``q_global`` (ejes globales).

        ``q`` se proyecta sobre los ejes locales mediante la matriz λ (3×3).
        En local, cada componente reparte (forma Hermite cúbica para
        transversales, lineal para axial):

        - Axial ``q_x_local``: mitad-mitad por nodo, sin momentos.
        - Transversal ``q_y_local``: ``q_y · L/2`` por nodo en F_y;
          ``+q_y · L²/12`` en Mz_i y ``-q_y · L²/12`` en Mz_j.
        - Transversal ``q_z_local``: ``q_z · L/2`` por nodo en F_z;
          ``-q_z · L²/12`` en My_i y ``+q_z · L²/12`` en My_j (signo
          opuesto a ``M_z`` por la regla de la mano derecha aplicada al
          eje de flexión y_local).
        - Sin torsión (sin excentricidad de la carga respecto al centro de
          cortadura).

        Los signos de los momentos siguen la convención stress-resultant /
        RHR pura del proyecto (Reglas.md §5) y son consistentes con la
        matriz ``K_local`` del elemento.
        """
        # Proyección a ejes locales: λ es la rotación 3×3 global→local.
        q_local = self.lam @ q_global  # (q_axial, q_y_local, q_z_local)
        qa, qy, qz = q_local[0], q_local[1], q_local[2]
        L = self.L0
        m_y = qz * L * L / 12.0
        m_z = qy * L * L / 12.0
        f_local = np.array([
            0.5 * qa * L,   # Fx_i
            0.5 * qy * L,   # Fy_i
            0.5 * qz * L,   # Fz_i
            0.0,            # Mx_i  (sin torsión por carga uniforme centrada)
            -m_y,           # My_i
            +m_z,           # Mz_i
            0.5 * qa * L,   # Fx_j
            0.5 * qy * L,   # Fy_j
            0.5 * qz * L,   # Fz_j
            0.0,            # Mx_j
            +m_y,           # My_j
            -m_z,           # Mz_j
        ])
        return self.T.T @ f_local
