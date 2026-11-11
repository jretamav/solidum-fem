"""Armaduras rectas de segundo orden: ``Truss2D3`` y ``Truss3D3``.

Barra de tres nodos ``[i, j, m]`` (extremos y nodo intermedio en el punto
medio) con desplazamiento axial cuadrático y deformación lineal; transmite
sólo fuerza axial. Specs ``docs/specs/Truss2D3.md`` y ``Truss3D3.md``
(Bathe 2014, §6.3.3; Abaqus T2D3/T3D3).

Con ``ξ = s/L ∈ [0, 1]``:

    N = [(1 − ξ)(1 − 2ξ),  ξ(2ξ − 1),  4ξ(1 − ξ)]
    B_loc = [4ξ − 3,  4ξ − 1,  4 − 8ξ] / L

Gauss de 2 puntos, exacto para la rigidez con ``E_t·A`` constante; cada
punto lleva su estado y su longitud característica ``w_g·L = L/2``
(ADR 0028).
"""
from __future__ import annotations

from typing import ClassVar, List

import numpy as np

from solidum.core.element import Element, validate_lumping_kwarg
from solidum.core.material import Material
from solidum.core.node import Node
from solidum.elements._line import _LINE_W, _LINE_XI, LineElementMixin, linear_tangent
from solidum.math.distributed_loads import field_values
from solidum.registry import ElementRegistry
from solidum.results import ElementForces

# Gauss-Legendre de 2 puntos en [0, 1].
_GAUSS_XI = np.array([0.5 - 0.5 / np.sqrt(3.0), 0.5 + 0.5 / np.sqrt(3.0)])
_GAUSS_W = np.array([0.5, 0.5])
# Extrapolación lineal de los puntos de Gauss a los extremos (ξ = 0 y ξ = 1):
# N(ξ) = N_1 + (N_2 − N_1)(ξ − ξ_1)/(ξ_2 − ξ_1).
_TO_ENDS = np.array([[1.0 + np.sqrt(3.0), 1.0 - np.sqrt(3.0)],
                     [1.0 - np.sqrt(3.0), 1.0 + np.sqrt(3.0)]]) / 2.0
# F_int·e en los extremos en función de (N_i, N_j) con N lineal:
# F_i·e = −(5N_i + N_j)/6,  F_j·e = (N_i + 5N_j)/6.
_ENDS_FROM_NODAL = np.linalg.inv(np.array([[-5.0, -1.0], [1.0, 5.0]]) / 6.0)
# Masa consistente por componente, orden [i, j, m], en unidades de ρAL.
_MASS = np.array([[4.0, -1.0, 2.0], [-1.0, 4.0, 2.0], [2.0, 2.0, 16.0]]) / 30.0
_MASS_LUMPED = np.array([1.0, 1.0, 4.0]) / 6.0
# Desviación admitida del nodo intermedio respecto del punto medio, relativa
# a la longitud del elemento.
_MIDPOINT_TOL = 1.0e-8


def quadratic_shape_functions(xi) -> tuple[np.ndarray, np.ndarray]:
    """``N`` y ``dN/dξ`` de la línea cuadrática en ``ξ ∈ [0, 1]``, cada una
    ``(n, 3)`` en el orden de nodos ``[i, j, m]``."""
    xi = np.atleast_1d(np.asarray(xi, dtype=np.float64))
    N = np.column_stack([(1.0 - xi) * (1.0 - 2.0 * xi), xi * (2.0 * xi - 1.0),
                         4.0 * xi * (1.0 - xi)])
    dN = np.column_stack([4.0 * xi - 3.0, 4.0 * xi - 1.0, 4.0 - 8.0 * xi])
    return N, dN


class _QuadraticTruss(LineElementMixin, Element):
    """Base de las armaduras de tres nodos: la dimensión la fija
    ``DOF_NAMES``."""
    STRAIN_DIM = 1
    SEEDS_CHARACTERISTIC_LENGTH = True  # siembra h con _seed_characteristic_lengths
    N_INTEGRATION_POINTS = 2
    REFERENCE_NODES: ClassVar[tuple] = ((-1.0,), (1.0,), (0.0,))
    VTK_CELL_TYPE: ClassVar[str] = "line3"
    LINE_FIELDS = ("u",)
    LINE_STRAINS = ("eps",)

    def __init__(self, element_id: int, nodes: List[Node], material: Material, A: float):
        name = type(self).__name__
        if len(nodes) != 3:
            raise ValueError(f"El elemento {name} requiere exactamente 3 nodos [i, j, m].")
        if not A > 0.0:
            raise ValueError(f"{name}(id={element_id}): el área A debe ser > 0 (recibido {A}).")
        self.A = A
        super().__init__(element_id, nodes, material)

        X = self._node_coordinates()
        axis = X[1] - X[0]
        self.L0 = float(np.linalg.norm(axis))
        if self.L0 == 0.0:
            raise ValueError("La longitud del elemento no puede ser cero.")
        offset = float(np.linalg.norm(X[2] - 0.5 * (X[0] + X[1])))
        if offset > _MIDPOINT_TOL * self.L0:
            raise ValueError(
                f"{name}(id={element_id}): el nodo intermedio debe estar en el punto "
                f"medio de los extremos (está a {offset:.3e}, L = {self.L0:.3e}); el "
                f"elemento es recto y no admite geometría curva.")
        self.e = axis / self.L0
        _, dN = quadratic_shape_functions(_GAUSS_XI)
        # B global en cada punto de Gauss: (2, 3·d), B_loc ⊗ eᵀ.
        self._B = np.stack([np.kron(dN[g] / self.L0, self.e) for g in range(2)])
        self._seed_characteristic_lengths(_GAUSS_W * self.L0)

    def _node_coordinates(self) -> np.ndarray:
        """Coordenadas ``(3, d)`` de los nodos con ``d`` las componentes de
        desplazamiento (las que faltan, cero)."""
        d = len(self.DOF_NAMES)
        X = np.zeros((3, d))
        for k, node in enumerate(self.nodes):
            c = np.asarray(node.coordinates, dtype=np.float64)[:d]
            X[k, :len(c)] = c
        return X

    def _line_ends(self) -> np.ndarray:
        return self._node_coordinates()[:2]

    # --- formulación --------------------------------------------------------

    def compute_element_state(self, u_e: np.ndarray):
        u_e = np.asarray(u_e, dtype=np.float64)
        ndof = self._B.shape[1]
        K_e = np.zeros((ndof, ndof))
        F_int = np.zeros(ndof)
        for g in range(2):
            B = self._B[g]
            sigma, E_t, new_state = self.material.compute_state(float(B @ u_e),
                                                                self.state.vars[g])
            self.state.vars_trial[g] = new_state
            self.state.stresses_trial[g] = sigma
            wL = _GAUSS_W[g] * self.L0
            K_e += (wL * E_t * self.A) * np.outer(B, B)
            F_int += (wL * sigma * self.A) * B
        return K_e, F_int

    def _gauss_axial(self, U_global: np.ndarray):
        """``(ε_g, σ_g)`` en los dos puntos de Gauss para ``U_global``, desde el
        estado confirmado del material (sin tocar el de prueba)."""
        u_e = np.asarray(self.get_local_displacements(U_global), dtype=np.float64)
        eps = self._B @ u_e
        sig = np.array([self.material.compute_state(float(eps[g]), self.state.vars[g])[0]
                        for g in range(2)])
        return eps, sig

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        """Utilidad de post-proceso: axial, esfuerzo y deformación en los dos
        puntos de Gauss."""
        eps, sig = self._gauss_axial(U_global)
        return {'axial_force': self.A * sig, 'stress': sig, 'strain': eps}

    def internal_forces(self, U_global: np.ndarray,
                        equivalent_load: np.ndarray | None = None) -> ElementForces:
        """API pública (ADR 0002): ``N`` en los extremos, tracción positiva
        (Reglas §5), por extrapolación lineal desde los puntos de Gauss.

        El axial del elemento es lineal, así que lo determinan sus dos valores
        de Gauss; con una carga distribuida el campo ya la incluye, por lo que
        ``equivalent_load`` no interviene (se acepta por el contrato)."""
        _, sig = self._gauss_axial(U_global)
        return ElementForces(kind="truss", components={"N": _TO_ENDS @ (self.A * sig)})

    def internal_forces_from_nodal(self, F_e, U_global, equivalent_load=None) -> ElementForces:
        """``N`` en los extremos a partir de unas fuerzas internas nodales
        ``F_e = ∫Bᵀ·N ds`` (ADR 0031): con ``N`` lineal, las proyecciones sobre
        el eje en los extremos son ``−(5N_i + N_j)/6`` y ``(N_i + 5N_j)/6``."""
        F = np.asarray(F_e, dtype=np.float64)
        d = len(self.e)
        ends = np.array([F[:d] @ self.e, F[d:2 * d] @ self.e])
        return ElementForces(kind="truss", components={"N": _ENDS_FROM_NODAL @ ends})

    def section_rigidities(self) -> np.ndarray:
        return np.array([linear_tangent(self.material) * self.A])

    def _line_interpolation(self, xi):
        N, dN = quadratic_shape_functions(xi)
        n = len(N)
        return N.reshape(n, 1, 3), (dN / self.L0).reshape(n, 1, 3)

    def _to_local(self, u_e):
        d = len(self.e)
        u = np.asarray(u_e, dtype=np.float64).reshape(3, d)
        return u @ self.e

    # --- cargas y masa ------------------------------------------------------

    def compute_body_load(self, b) -> np.ndarray:
        """Vector nodal consistente ``∫ Nᵀ·A·b dx`` en cada componente global.

        Uniforme: ``b·A·L·(1/6, 1/6, 2/3)``. Un campo ``b(x)`` se integra con
        Gauss de 4 puntos (exacto con ``b`` de grado ≤ 5), como las armaduras
        de dos nodos (ADR 0022)."""
        d = len(self.e)
        if callable(b):
            X = self._line_ends()
            x = X[0] + _LINE_XI[:, None] * (X[1] - X[0])
            values = field_values(b, x, d, "de cuerpo")
            N, _ = quadratic_shape_functions(_LINE_XI)
            return (self.A * self.L0 * (N * _LINE_W[:, None]).T @ values).reshape(-1)
        b = np.asarray(b, dtype=np.float64).reshape(d)
        return np.kron(self.A * self.L0 * _MASS_LUMPED, b)

    def compute_mass_matrix(self, lumping: str = "consistent") -> np.ndarray:
        """Masa elemental (ADR 0009 §1): consistente ``ρAL/30·[[4, −1, 2],
        [−1, 4, 2], [2, 2, 16]] ⊗ I``; concentrada ``ρAL·(1/6, 1/6, 2/3)`` por
        componente (Lobatto de 3 puntos y HRZ coinciden)."""
        validate_lumping_kwarg(lumping, type(self).__name__)
        m_e = self.material.density * self.A * self.L0
        I = np.eye(len(self.e))
        if lumping == "lumped":
            return np.kron(np.diag(m_e * _MASS_LUMPED), I)
        return np.kron(m_e * _MASS, I)


@ElementRegistry.register
class Truss2D3(_QuadraticTruss):
    """Armadura 2D recta de segundo orden (3 nodos ``[i, j, m]``), pequeños
    desplazamientos. Spec ``docs/specs/Truss2D3.md``.

    Parameters
    ----------
    element_id : int
    nodes : List[Node]
        Extremos ``i``, ``j`` y nodo intermedio ``m`` en el punto medio.
    material : Material
        Material 1D (``STRAIN_DIM = 1``).
    A : float
        Área de la sección transversal.
    """
    DOF_NAMES = ['ux', 'uy']


@ElementRegistry.register
class Truss3D3(_QuadraticTruss):
    """Armadura 3D recta de segundo orden (3 nodos ``[i, j, m]``), pequeños
    desplazamientos. Spec ``docs/specs/Truss3D3.md``. Acepta nodos con 2 ó 3
    coordenadas (``z = 0`` si falta)."""
    DOF_NAMES = ['ux', 'uy', 'uz']
