# fenix_fem/fenix/elements/frame3d.py
"""Elemento de marco/viga 3D Euler-Bernoulli.

Implementación **independiente** del resto de elementos. No hereda de
ningún otro elemento ni comparte helpers: la construcción de ejes
locales, la matriz de transformación 12×12 y la rigidez local se generan
íntegramente dentro de la clase. Ver docs/specs/Frame3D.md.
"""
import math
from typing import List, Optional

import numpy as np

from fenix.core.element import Element
from fenix.core.material import Material
from fenix.core.node import Node
from fenix.logging import get_logger
from fenix.registry import ElementRegistry
from fenix.results import ElementForces

_log = get_logger("elements.frame3d")


@ElementRegistry.register
class Frame3D(Element):
    """Marco/viga 3D basado en Euler-Bernoulli.

    Dos nodos rígidamente conectados, 6 DOFs por nodo (ux, uy, uz, rx, ry, rz).
    Transmite axial + cortante y flexión en dos planos + torsión. Régimen de
    linealidad geométrica.

    Parameters
    ----------
    element_id : int
    nodes : List[Node]
        Exactamente 2 nodos extremos. Acepta nodos con 3 coordenadas.
    material : Material
        Material 1D (STRAIN_DIM=1). Si expone `nu`, se usa ese valor.
    A : float
        Área de la sección transversal.
    Iy, Iz : float
        Momentos de inercia respecto a los ejes locales y, z de la sección.
    J : float
        Constante torsional de Saint-Venant.
    nu : float, optional
        Coeficiente de Poisson si el material no lo expone. Default 0.3.
    ref_vector : Optional[List[float]], optional
        Vector de referencia 3D para fijar la orientación de los ejes
        locales y, z. Default `[0, 0, 1]`; con fallback `[1, 0, 0]` si
        la barra es cercanamente vertical.
    """

    DOF_NAMES = ['ux', 'uy', 'uz', 'rx', 'ry', 'rz']
    STRAIN_DIM = 1
    N_INTEGRATION_POINTS = 1

    _VERTICAL_COS_TOL = 0.99  # umbral para detectar barras casi-verticales

    def __init__(self, element_id: int, nodes: List[Node], material: Material,
                 A: float, Iy: float, Iz: float, J: float,
                 nu: float = 0.3,
                 ref_vector: Optional[List[float]] = None):
        if len(nodes) != 2:
            raise ValueError("El elemento Frame3D requiere exactamente 2 nodos.")

        self.A = A
        self.Iy = Iy
        self.Iz = Iz
        self.J = J
        super().__init__(element_id, nodes, material)

        # Poisson: material.nu si existe, si no el parámetro
        if hasattr(material, 'nu'):
            self.nu = material.nu
        else:
            if nu == 0.3:
                _log.warning(
                    f"Frame3D (id={element_id}): el material no expone 'nu'. "
                    f"Se usará nu={nu} (default). Especifique 'nu' en el YAML si esto es incorrecto."
                )
            self.nu = nu

        self.L0, self.lam, self.T = self._build_local_frame(self.nodes, ref_vector,
                                                             element_id)

    @classmethod
    def _build_local_frame(cls, nodes: List[Node],
                            ref_vector: Optional[List[float]],
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
                _log.warning(
                    f"Frame3D (id={element_id}): barra cercanamente vertical sin ref_vector explícito. "
                    f"Se usa fallback [1, 0, 0]. Especifique ref_vector si la orientación de la sección debe ser otra."
                )
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
        u_local = self.T @ u_e
        L = self.L0

        epsilon = (u_local[6] - u_local[0]) / L
        sigma, E_t, new_state = self.material.compute_state(epsilon, self.state.vars[0])
        self.state.vars_trial[0] = new_state
        self.state.stresses_trial[0] = sigma

        G = E_t / (2.0 * (1.0 + self.nu))
        K_local = self._build_local_stiffness(E_t, G)

        F_int_local = K_local @ u_local
        # Componente axial: usar sigma directo del material
        F_int_local[0] = -sigma * self.A
        F_int_local[6] =  sigma * self.A

        K_global = self.T.T @ K_local @ self.T
        F_int_e = self.T.T @ F_int_local
        return K_global, F_int_e

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        u_e = self.get_local_displacements(U_global)
        _, F_int = self.compute_element_state(u_e)
        F_local = self.T @ F_int

        u_local = self.T @ u_e
        epsilon = (u_local[6] - u_local[0]) / self.L0
        sigma, _, _ = self.material.compute_state(epsilon, self.state.vars[0])

        return {
            'axial_force': F_local[0],
            'shear_y': F_local[1],
            'shear_z': F_local[2],
            'torsion': F_local[3],
            'moment_y_i': F_local[4],
            'moment_z_i': F_local[5],
            'moment_y_j': F_local[10],
            'moment_z_j': F_local[11],
            'stress': sigma,
            'strain': epsilon,
        }

    def internal_forces(self, U_global: np.ndarray) -> ElementForces:
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
        _, F_int = self.compute_element_state(u_e)
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

    def _build_local_mass(self, rho: float) -> np.ndarray:
        """Matriz de masa consistente 12×12 en ejes locales del elemento.

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
        c_ax = rho * self.A * L / 6.0
        c_tr = rho * self.A * L / 420.0
        Jp = self.Iy + self.Iz
        c_tor = rho * Jp * L / 6.0
        c_rot_y = rho * self.Iy * L / 6.0
        c_rot_z = rho * self.Iz * L / 6.0

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

    def compute_mass_matrix(self, lumping: str = "consistent") -> np.ndarray:
        """Matriz de masa consistente del Frame3D en ejes globales.

        Ensambla :meth:`_build_local_mass` y la rota con ``T^T · M_local · T``,
        donde ``T`` es la matriz 12×12 ya construida en ``__init__``. La
        inercia rotacional propia de sección (``ρI``) se omite (Bernoulli
        esbelta — ver ADR 0009 §1 y :meth:`_build_local_mass`).

        Returns
        -------
        np.ndarray, shape (12, 12)
        """
        if lumping != "consistent":
            raise NotImplementedError(
                f"Frame3D.compute_mass_matrix: lumping='{lumping}' no "
                f"implementado. Fase 1 (ADR 0009) solo admite 'consistent'."
            )
        M_local = self._build_local_mass(self.material.density)
        return self.T.T @ M_local @ self.T

    def compute_body_load(self, b: np.ndarray) -> np.ndarray:
        """Vector nodal consistente con ∫NᵀbA dx para marco 3D (peso propio).

        La carga distribuida por unidad de longitud ``q = A · b`` se proyecta
        sobre los ejes locales mediante la matriz λ (3×3). En local, cada
        componente reparte (forma Hermite cúbica para transversales, lineal
        para axial):

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
        # Carga distribuida por unidad de longitud en globales.
        q_global = self.A * b
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
