# solidum_fem/solidum/elements/frame3d.py
"""Elemento de marco/viga 3D Euler-Bernoulli.

Implementación **independiente** del resto de elementos: no hereda de
ningún otro elemento; la construcción de ejes locales, la matriz de
transformación 12×12 y la rigidez local se generan íntegramente dentro de
la clase. Comparte con los marcos planos sólo la constitutiva genérica de
``solidum.elements.frame._section`` (ADR 0021: forma cerrada con material
elástico lineal o integración de una sección). Ver docs/specs/Frame3D.md.
"""
from __future__ import annotations

import math
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
from solidum.logging import get_logger
from solidum.registry import ElementRegistry, MaterialRegistry, SectionRegistry
from solidum.results import ElementForces

_log = get_logger("elements.frame3d")


@ElementRegistry.register
class Frame3D(FrameSectionMixin, Element):
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

    def _build_local_mass(self, rhoA: float, rhoIy: float, rhoIz: float,
                          rhoJp: float) -> np.ndarray:
        """Matriz de masa consistente 12×12 en ejes locales del elemento.

        Recibe las propiedades másicas por unidad de longitud: ``ρA``,
        ``ρI_y``, ``ρI_z`` y ``ρJ_p`` (con material, ``ρ`` por las
        propiedades geométricas; con sección, las de
        ``section.mass_properties()`` y ``ρJ_p = ρI_y + ρI_z``).

        Contribuciones desacopladas, ensambladas en el layout local
        ``[ux_i, uy_i, uz_i, rx_i, ry_i, rz_i, ux_j, uy_j, uz_j, rx_j, ry_j, rz_j]``:

        - **Axial** (ux_i, ux_j): ``(ρAL/6) · [[2,1],[1,2]]``.
        - **Torsión** (rx_i, rx_j): ``(ρ·Jp·L/6) · [[2,1],[1,2]]`` con
          ``Jp = Iy + Iz`` (momento polar geométrico de área, no la constante
          torsional de Saint-Venant ``J`` — esa rige la rigidez, no la inercia).
        - **Flexión xy** (uy_i, rz_i, uy_j, rz_j): masa Hermitiana cúbica
          ``(ρAL/420) · H`` con la convención sign-aware ``θ_z = +dw_y/dx``.
        - **Flexión xz** (uz_i, ry_i, uz_j, ry_j): masa Hermitiana con
          ``θ_y = −dw_z/dx`` (signo opuesto al de la flexión xy en todos los
          términos que llevan una rotación una sola vez). Coherente con la
          construcción de la rigidez ``K_local`` en este mismo elemento.
        - **Inercia rotacional propia de sección** en flexión (ADR 0009 §1):
          ``ρ·Iz·L/6 · [[2,1],[1,2]]`` para los DOFs ``rz_i, rz_j`` y
          ``ρ·Iy·L/6 · [[2,1],[1,2]]`` para los DOFs ``ry_i, ry_j``. Recoge
          la energía cinética de rotación de la sección alrededor de cada
          eje local. La componente torsional ``rx`` ya está cubierta por el
          término torsional (con ``Jp = Iy + Iz``).
        """
        L = self.L0
        L2 = L * L
        c_ax = rhoA * L / 6.0
        c_tr = rhoA * L / 420.0
        c_tor = rhoJp * L / 6.0
        c_rot_y = rhoIy * L / 6.0
        c_rot_z = rhoIz * L / 6.0

        M = np.zeros((12, 12))

        # Axial (0 ↔ 6)
        M[0, 0]  = 2.0 * c_ax;  M[0, 6]  = 1.0 * c_ax
        M[6, 0]  = 1.0 * c_ax;  M[6, 6]  = 2.0 * c_ax

        # Torsión (3 ↔ 9)
        M[3, 3]  = 2.0 * c_tor; M[3, 9]  = 1.0 * c_tor
        M[9, 3]  = 1.0 * c_tor; M[9, 9]  = 2.0 * c_tor

        # Flexión xy: idx_xy = [uy_i=1, rz_i=5, uy_j=7, rz_j=11]
        H_xy = c_tr * np.array([
            [156.0,   22.0 * L,   54.0,  -13.0 * L],
            [22.0 * L,  4.0 * L2,  13.0 * L, -3.0 * L2],
            [54.0,    13.0 * L,  156.0,  -22.0 * L],
            [-13.0 * L, -3.0 * L2, -22.0 * L,  4.0 * L2],
        ])
        idx_xy = [1, 5, 7, 11]
        for i, gi in enumerate(idx_xy):
            for j, gj in enumerate(idx_xy):
                M[gi, gj] = H_xy[i, j]

        # Flexión xz: idx_xz = [uz_i=2, ry_i=4, uz_j=8, ry_j=10]. Mismos
        # coeficientes que H_xy, pero con signo invertido en los términos
        # que involucran exactamente una rotación (filas/columnas 1 y 3).
        S = np.diag([1.0, -1.0, 1.0, -1.0])
        H_xz = S @ H_xy @ S
        idx_xz = [2, 4, 8, 10]
        for i, gi in enumerate(idx_xz):
            for j, gj in enumerate(idx_xz):
                M[gi, gj] = H_xz[i, j]

        # Inercia rotacional propia de sección sobre ry (5,11) y rz (4,10)
        # son los DOFs de rotación; idx_ry = [4, 10], idx_rz = [5, 11].
        M[5, 5]  += 2.0 * c_rot_z; M[5, 11] += 1.0 * c_rot_z
        M[11, 5] += 1.0 * c_rot_z; M[11, 11] += 2.0 * c_rot_z
        M[4, 4]  += 2.0 * c_rot_y; M[4, 10] += 1.0 * c_rot_y
        M[10, 4] += 1.0 * c_rot_y; M[10, 10] += 2.0 * c_rot_y

        return M

    def _build_local_mass_lumped(self, rhoA: float, rhoIy: float, rhoIz: float,
                                 rhoJp: float) -> np.ndarray:
        """Matriz de masa lumped 12×12 en ejes locales (ADR 0009 fase 2).

        Mismas propiedades másicas por unidad de longitud que
        :meth:`_build_local_mass`.

        Esquema nodal directo:

        - Traslacional ``m_t = ρAL/2`` por DOF (``ux``, ``uy``, ``uz``).
        - Rotacional torsional ``m_rx = ρ·Jp·L/2`` con ``Jp = Iy + Iz``
          (momento polar geométrico).
        - Rotacional flexional ``m_ry = ρ·Iy·L/2``, ``m_rz = ρ·Iz·L/2``.

        Diagonal en local. **En global la diagonalidad se preserva sólo
        cuando el eje del elemento coincide con un eje global** (la
        rotación T no mezcla DOFs locales con globales en los grados
        rotacionales). En orientación oblicua, el bloque traslacional 3×3
        por nodo permanece diagonal porque ``m_t·I_3`` es invariante a
        SO(3); el bloque rotacional 3×3 por nodo no es escalar porque
        ``m_rx = ρJp·L/2``, ``m_ry = ρIy·L/2``, ``m_rz = ρIz·L/2`` son
        valores siempre distintos para una sección real (la torsional usa
        el momento polar geométrico ``Jp = Iy + Iz``, las flexionales sólo
        ``Iy`` o ``Iz``), de modo que tras rotar T el bloque rotacional
        queda lleno. M_global resulta entonces **bloque-diagonal** por
        nodo. Es la limitación estándar del lumping en frames 3D
        (Cook-Malkus-Plesha §11.4): aceptable para Newmark/HHT
        (factorización requerida igualmente) y para diferencias centradas
        con inversión local de bloques 6×6 baratos por nodo.

        Layout local: ``[ux_i, uy_i, uz_i, rx_i, ry_i, rz_i,
        ux_j, uy_j, uz_j, rx_j, ry_j, rz_j]``.

        Returns
        -------
        np.ndarray, shape (12, 12)
        """
        L = self.L0
        m_t = rhoA * L / 2.0
        m_rx = rhoJp * L / 2.0
        m_ry = rhoIy * L / 2.0
        m_rz = rhoIz * L / 2.0
        return np.diag([m_t, m_t, m_t, m_rx, m_ry, m_rz,
                        m_t, m_t, m_t, m_rx, m_ry, m_rz])

    def compute_mass_matrix(self, lumping: str = "consistent") -> np.ndarray:
        """Matriz de masa del Frame3D en ejes globales (ADR 0009).

        **Consistente** (fase 1): :meth:`_build_local_mass` (axial +
        Hermitiana cúbica en los dos planos de flexión + inercia
        rotacional de sección ``ρI`` en ambos ejes locales + torsión
        ``ρJp``).

        **Lumped** (fase 2): :meth:`_build_local_mass_lumped`. Ver allí la
        nota sobre Iy ≠ Iz y bloque-diagonalidad en global.

        Con sección, las propiedades másicas son las de
        ``section.mass_properties()`` (spec ``FiberSection`` §6), con la
        inercia torsional ``ρJ_p = ρI_y + ρI_z``.

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
        """``(ρA, ρI_y, ρI_z, ρJ_p)`` por unidad de longitud. Con material,
        ``ρ`` por ``A``, ``I_y``, ``I_z`` y el momento polar geométrico
        ``J_p = I_y + I_z``; con sección, ``section.mass_properties()``."""
        if self.section is None:
            rho = self.material.density
            return (rho * self.A, rho * self.Iy, rho * self.Iz,
                    rho * (self.Iy + self.Iz))
        mp = self.section.mass_properties()
        return mp["rhoA"], mp["rhoIy"], mp["rhoIz"], mp["rhoIy"] + mp["rhoIz"]

    def compute_body_load(self, b: np.ndarray) -> np.ndarray:
        """Vector nodal consistente con ∫NᵀbA dx para marco 3D (peso propio).

        La carga distribuida por unidad de longitud ``q = A · b`` (``A`` el
        área bruta: ``Σ A_i`` con sección) se reparte con
        :meth:`_line_load_vector`.

        Parameters
        ----------
        b : np.ndarray, shape (3,)
            Fuerza de cuerpo por unidad de volumen en ejes globales (e.g.
            ``b = (0, 0, -ρ·g)`` para peso propio con ``+z`` hacia arriba).

        Returns
        -------
        np.ndarray, shape (12,)
            Vector consistente en ejes globales, layout estándar de Frame3D:
            ``[Fx_i, Fy_i, Fz_i, Mx_i, My_i, Mz_i, Fx_j, ..., Mz_j]``.
        """
        b = np.asarray(b, dtype=float).reshape(3)
        return self._line_load_vector(self._axial_area() * b)

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
