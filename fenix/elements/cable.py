# fenix_fem/fenix/elements/cable.py
"""Elementos de cable.

Implementación **independiente** de los elementos de armadura (`Truss*`).
La duplicación de la maquinaria cinemática respecto a `Truss2DCorot` es
deliberada: si la formulación de armaduras cambia, los cables no se ven
arrastrados, y viceversa. Ver docs/specs/Cable2DCorot.md.
"""
import math
from typing import List

import numpy as np

from fenix.core.element import Element, validate_lumping_kwarg
from fenix.core.material import Material
from fenix.core.node import Node
from fenix.math.geometry import perpendicular_projector
from fenix.registry import ElementRegistry
from fenix.results import ElementForces


@ElementRegistry.register
class Cable2DCorot(Element):
    """Cable 2D corotacional — tensión unilateral, grandes rotaciones.

    Dos nodos articulados en el plano, 2 DOFs por nodo (ux, uy). La
    cinemática corotacional (Updated Lagrangian) recalcula longitud y
    cosenos directores en configuración corriente en cada evaluación.

    La propiedad de unilateralidad no vive en esta clase — el elemento
    delega la relación constitutiva al material. Para obtener respuesta
    físicamente correcta de cable, emparejar con un material unilateral
    como `CableMaterial1D`.

    Parameters
    ----------
    element_id : int
    nodes : List[Node]
        Exactamente 2 nodos extremos.
    material : Material
        Material 1D (STRAIN_DIM=1). Para comportamiento de cable genuino,
        debe ser unilateral (p. ej. CableMaterial1D).
    A : float
        Área de la sección transversal del cable.
    """

    DOF_NAMES = ['ux', 'uy']
    STRAIN_DIM = 1
    N_INTEGRATION_POINTS = 1
    ACCEPTS_UNILATERAL = True

    def __init__(self, element_id: int, nodes: List[Node], material: Material,
                 A: float):
        if len(nodes) != 2:
            raise ValueError("El elemento Cable2DCorot requiere exactamente 2 nodos.")

        self.A = A
        super().__init__(element_id, nodes, material)

        coords1 = self.nodes[0].coordinates
        coords2 = self.nodes[1].coordinates
        dx0, dy0 = coords2[0] - coords1[0], coords2[1] - coords1[1]
        self.L0 = math.sqrt(dx0**2 + dy0**2)
        if self.L0 == 0.0:
            raise ValueError("La longitud del cable no puede ser cero.")

    def _current_geometry(self, u_e: np.ndarray):
        """Longitud y cosenos directores en configuración corriente."""
        x1 = self.nodes[0].coordinates[0] + u_e[0]
        y1 = self.nodes[0].coordinates[1] + u_e[1]
        x2 = self.nodes[1].coordinates[0] + u_e[2]
        y2 = self.nodes[1].coordinates[1] + u_e[3]
        dx, dy = x2 - x1, y2 - y1
        l = math.sqrt(dx**2 + dy**2)
        if l == 0.0:
            raise ValueError("Longitud corriente nula en Cable2DCorot.")
        return l, dx / l, dy / l

    def compute_element_state(self, u_e: np.ndarray):
        l, c_t, s_t = self._current_geometry(u_e)

        epsilon = (l - self.L0) / self.L0
        sigma, E_t, new_state = self.material.compute_state(epsilon, self.state.vars[0])
        self.state.vars_trial[0] = new_state
        self.state.stresses_trial[0] = sigma

        N = sigma * self.A
        d = np.array([-c_t, -s_t, c_t, s_t])
        n = np.array([-s_t, c_t, s_t, -c_t])

        K_M = ((E_t * self.A) / self.L0) * np.outer(d, d)
        K_G = (N / l) * np.outer(n, n)
        F_int_e = N * d

        return K_M + K_G, F_int_e

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        u_e = self.get_local_displacements(U_global)
        l, _, _ = self._current_geometry(u_e)
        epsilon = (l - self.L0) / self.L0
        sigma, _, _ = self.material.compute_state(epsilon, self.state.vars[0])
        N = self.A * sigma
        return {'axial_force': N, 'stress': sigma, 'strain': epsilon}

    def internal_forces(self, U_global: np.ndarray) -> ElementForces:
        """API pública (ADR 0002): N en ejes locales, tracción positiva (§5).

        La unilateralidad (N≥0, cero cuando slack) está delegada al material;
        si se emparejó con CableMaterial1D, el valor ya refleja la rama correcta.
        """
        N = self.compute_internal_forces(U_global)['axial_force']
        return ElementForces(kind="cable", components={"N": np.array([N, N])})

    def compute_body_load(self, b: np.ndarray) -> np.ndarray:
        """Vector nodal consistente con ∫NᵀbA dx en configuración de referencia.

        Reparto exacto ``b · A · L₀ / 2`` por nodo, en cada componente global,
        evaluado sobre la geometría inicial. Aproximación estándar para
        cargas de cuerpo en elementos corotacionales: para grandes rotaciones
        con cargas conservadoras (gravedad), la diferencia frente a la
        integración sobre la geometría corriente es de segundo orden.

        Parameters
        ----------
        b : np.ndarray, shape (2,)
            Fuerza de cuerpo por unidad de volumen en ejes globales.

        Returns
        -------
        np.ndarray, shape (4,)
            ``[Fx_i, Fy_i, Fx_j, Fy_j]`` en ejes globales.
        """
        b = np.asarray(b, dtype=float).reshape(2)
        half = 0.5 * self.A * self.L0
        return np.array([half * b[0], half * b[1], half * b[0], half * b[1]])

    def compute_mass_matrix(self, lumping: str = "consistent") -> np.ndarray:
        """Matriz de masa elemental (ADR 0009 §1, fase 2 con ``lumping="lumped"``).

        **Consistente**: idéntica en forma a Truss2D
        ``M = (ρAL₀/6) · kron([[2,1],[1,2]], I₂)``, evaluada en configuración
        de referencia (Lagrangeano total). La unilateralidad del cable
        afecta a la rigidez, no a la inercia: la masa total siempre es física.

        **Lumped**: masa nodal uniforme ``ρAL₀/2`` por DOF traslacional
        (HRZ = nodal directo en barras articuladas).

        Returns
        -------
        np.ndarray, shape (4, 4)
        """
        validate_lumping_kwarg(lumping, type(self).__name__)
        m_e = self.material.density * self.A * self.L0
        if lumping == "lumped":
            return np.diag([0.5 * m_e] * 4)
        coef = m_e / 6.0
        return coef * np.array([
            [2.0, 0.0, 1.0, 0.0],
            [0.0, 2.0, 0.0, 1.0],
            [1.0, 0.0, 2.0, 0.0],
            [0.0, 1.0, 0.0, 2.0],
        ])


@ElementRegistry.register
class Cable3DCorot(Element):
    """Cable 3D corotacional — tensión unilateral, grandes rotaciones en el espacio.

    Dos nodos articulados, 3 DOFs por nodo (ux, uy, uz). Cinemática
    corotacional 3D: longitud y cosenos directores se recalculan en
    configuración corriente. La respuesta física de cable (unilateralidad)
    vive en el material; el elemento delega sin filtrar por signo.

    Clase totalmente autónoma: hereda directamente de Element, sin relación
    con Cable2DCorot ni con los elementos de armadura.

    Parameters
    ----------
    element_id : int
    nodes : List[Node]
        Exactamente 2 nodos extremos. Acepta nodos con 2 ó 3 coordenadas
        (completa con z=0 si falta).
    material : Material
        Material 1D (STRAIN_DIM=1). Para cable genuino, unilateral
        (p. ej. CableMaterial1D).
    A : float
        Área de la sección transversal del cable.
    """

    DOF_NAMES = ['ux', 'uy', 'uz']
    STRAIN_DIM = 1
    N_INTEGRATION_POINTS = 1
    ACCEPTS_UNILATERAL = True

    def __init__(self, element_id: int, nodes: List[Node], material: Material,
                 A: float):
        if len(nodes) != 2:
            raise ValueError("El elemento Cable3DCorot requiere exactamente 2 nodos.")

        self.A = A
        super().__init__(element_id, nodes, material)

        c1 = self.nodes[0].coordinates
        c2 = self.nodes[1].coordinates
        z1_0 = c1[2] if len(c1) > 2 else 0.0
        z2_0 = c2[2] if len(c2) > 2 else 0.0
        dx0, dy0, dz0 = c2[0] - c1[0], c2[1] - c1[1], z2_0 - z1_0
        self.L0 = math.sqrt(dx0**2 + dy0**2 + dz0**2)
        if self.L0 == 0.0:
            raise ValueError("La longitud del cable no puede ser cero.")

    def _current_geometry(self, u_e: np.ndarray):
        """Longitud y cosenos directores en configuración corriente."""
        c1 = self.nodes[0].coordinates
        c2 = self.nodes[1].coordinates
        z1 = c1[2] if len(c1) > 2 else 0.0
        z2 = c2[2] if len(c2) > 2 else 0.0
        x1 = c1[0] + u_e[0]; y1 = c1[1] + u_e[1]; zn1 = z1 + u_e[2]
        x2 = c2[0] + u_e[3]; y2 = c2[1] + u_e[4]; zn2 = z2 + u_e[5]
        dx, dy, dz = x2 - x1, y2 - y1, zn2 - zn1
        l = math.sqrt(dx**2 + dy**2 + dz**2)
        if l == 0.0:
            raise ValueError("Longitud corriente nula en Cable3DCorot.")
        return l, dx / l, dy / l, dz / l

    def compute_element_state(self, u_e: np.ndarray):
        l, cx, cy, cz = self._current_geometry(u_e)

        epsilon = (l - self.L0) / self.L0
        sigma, E_t, new_state = self.material.compute_state(epsilon, self.state.vars[0])
        self.state.vars_trial[0] = new_state
        self.state.stresses_trial[0] = sigma

        N = sigma * self.A
        d = np.array([-cx, -cy, -cz, cx, cy, cz])
        P = perpendicular_projector(np.array([cx, cy, cz]))

        K_M = ((E_t * self.A) / self.L0) * np.outer(d, d)

        K_G = np.zeros((6, 6))
        K_G[:3, :3] = P
        K_G[3:, 3:] = P
        K_G[:3, 3:] = -P
        K_G[3:, :3] = -P
        K_G *= (N / l)

        F_int_e = N * d
        return K_M + K_G, F_int_e

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        u_e = self.get_local_displacements(U_global)
        l, _, _, _ = self._current_geometry(u_e)
        epsilon = (l - self.L0) / self.L0
        sigma, _, _ = self.material.compute_state(epsilon, self.state.vars[0])
        N = self.A * sigma
        return {'axial_force': N, 'stress': sigma, 'strain': epsilon}

    def internal_forces(self, U_global: np.ndarray) -> ElementForces:
        """API pública (ADR 0002): N en ejes locales, tracción positiva (§5).

        La unilateralidad está delegada al material (ver Cable2DCorot).
        """
        N = self.compute_internal_forces(U_global)['axial_force']
        return ElementForces(kind="cable", components={"N": np.array([N, N])})

    def compute_body_load(self, b: np.ndarray) -> np.ndarray:
        """Vector nodal consistente con ∫NᵀbA dx para cable 3D.

        Reparto exacto ``b · A · L₀ / 2`` por nodo, en cada componente global,
        sobre la geometría inicial. Ver Cable2DCorot para la justificación
        del uso de la configuración de referencia.

        Parameters
        ----------
        b : np.ndarray, shape (3,)
            Fuerza de cuerpo por unidad de volumen en ejes globales.

        Returns
        -------
        np.ndarray, shape (6,)
            ``[Fx_i, Fy_i, Fz_i, Fx_j, Fy_j, Fz_j]`` en ejes globales.
        """
        b = np.asarray(b, dtype=float).reshape(3)
        half = 0.5 * self.A * self.L0
        return np.array([half * b[0], half * b[1], half * b[2],
                         half * b[0], half * b[1], half * b[2]])

    def compute_mass_matrix(self, lumping: str = "consistent") -> np.ndarray:
        """Matriz de masa elemental (ADR 0009 §1, fase 2 con ``lumping="lumped"``).

        **Consistente**: análoga a Truss3D
        ``M = (ρAL₀/6) · kron([[2,1],[1,2]], I₃)``.

        **Lumped**: masa nodal uniforme ``ρAL₀/2`` por DOF traslacional.

        Returns
        -------
        np.ndarray, shape (6, 6)
        """
        validate_lumping_kwarg(lumping, type(self).__name__)
        m_e = self.material.density * self.A * self.L0
        if lumping == "lumped":
            return np.diag([0.5 * m_e] * 6)
        coef = m_e / 6.0
        I3 = np.eye(3)
        M = np.block([[2.0 * I3, I3       ],
                      [I3,       2.0 * I3]])
        return coef * M
