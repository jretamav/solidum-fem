"""Interfaz de adherencia-deslizamiento entre un tramo de barra embebida y su
anfitrión (spec ``EmbeddedReinforcement`` §12, ADR 0032 y ADR 0033).

La genera el refuerzo embebido, una por tramo, cuando se le da ``bond``; no se
registra ni se cita en el YAML por ``type``. Une los nodos del elemento de
barra (2 o 3) con los del anfitrión. En el punto ``ξ`` del tramo, el
desplazamiento relativo (barra menos concreto) es

    Δu(ξ) = Σ_a N^b_a(ξ)·u^b_a − Σ_i N^h_i(ξ_h(ξ))·u^h_i,

y en el sistema local del tramo, normales primero y el eje ``t`` al final,
``Δ = R·Δu``. La ley de interfaz da la tracción ``t(Δ)``; la normal actúa
sobre el ancho ``d_b`` y la de adherencia sobre el perímetro ``p``:

    F_int = Σ_g w_g·L·B_gᵀ·Rᵀ·W·t_g,   K = Σ_g w_g·L·B_gᵀ·Rᵀ·W·D_g·R·B_g,

con ``B_g = [N^b ⊗ I, −N^h ⊗ I]`` y ``W = diag(d_b, [d_b], p)``.
Integración nodal (Lobatto en los nodos del tramo) o de Gauss.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

import numpy as np

from solidum.core.element import Element, validate_lumping_kwarg
from solidum.core.interface_material import InterfaceMaterial
from solidum.core.dof_names import TRANSLATIONS
from solidum.elements.truss_quadratic import quadratic_shape_functions
from solidum.mesh.point_location import natural_coordinates, shape_values

INTEGRATIONS = ("nodal", "gauss")

# Reglas sobre el tramo, ξ ∈ [0, 1], en el orden de nodos de la barra.
_RULES = {
    (2, "nodal"): (np.array([0.0, 1.0]), np.array([0.5, 0.5])),
    (2, "gauss"): (np.array([0.5 - 0.5 / np.sqrt(3.0), 0.5 + 0.5 / np.sqrt(3.0)]),
                   np.array([0.5, 0.5])),
    (3, "nodal"): (np.array([0.0, 1.0, 0.5]), np.array([1.0, 1.0, 4.0]) / 6.0),
    (3, "gauss"): (np.array([0.5 - 0.5 * np.sqrt(0.6), 0.5, 0.5 + 0.5 * np.sqrt(0.6)]),
                   np.array([5.0, 8.0, 5.0]) / 18.0),
}


def bar_shape_functions(n_nodes: int, xi: np.ndarray) -> np.ndarray:
    """``N^b`` ``(n, n_nodes)`` de la barra: lineal (2 nodos) o cuadrática
    (3 nodos ``[i, j, m]``)."""
    xi = np.atleast_1d(np.asarray(xi, dtype=np.float64))
    if n_nodes == 2:
        return np.column_stack([1.0 - xi, xi])
    return quadratic_shape_functions(xi)[0]


def local_axes(t: np.ndarray) -> np.ndarray:
    """``R`` ``(d, d)``: filas ``n`` (2D) o ``n1, n2`` (3D) y después ``t``."""
    t = np.asarray(t, dtype=np.float64)
    if t.size == 2:
        return np.array([[-t[1], t[0]], [t[0], t[1]]])
    ref = np.eye(3)[int(np.argmin(np.abs(t)))]
    n1 = ref - (ref @ t) * t
    n1 /= np.linalg.norm(n1)
    return np.array([n1, np.cross(t, n1), t])


@dataclass(frozen=True)
class InterfacePoint:
    """Estado de un punto de integración de la interfaz."""
    position: np.ndarray        # coordenadas del punto
    slip: float                 # s: barra respecto del concreto a lo largo de t [m]
    bond_stress: float          # τ [Pa]
    normal_opening: np.ndarray  # Δu_n [m]
    normal_traction: np.ndarray  # t_n [Pa]
    bond_force: float           # p·τ, fuerza de adherencia por unidad de longitud [N/m]


class BondSlipInterface(Element):
    """Interfaz distribuida barra–sólido de un tramo (ver el módulo).

    Parameters
    ----------
    element_id : int
    bar : Element
        Elemento de barra del tramo (2 o 3 nodos, recto).
    host : Element
        Elemento de sólido que contiene el tramo.
    law : InterfaceMaterial
        Ley de interfaz de línea (``GEOMETRY = "line"``).
    diameter : float
        ``d_b``: ancho de la tracción normal.
    perimeter : float, optional
        Perímetro de la tracción de adherencia; por omisión ``π·d_b``.
    integration : {"nodal", "gauss"}
    """

    # Sin celda VTK propia: sus resultados se leen con interface_state().
    VTK_EXPORT: ClassVar[bool] = False

    def __init__(self, element_id: int, bar: Element, host: Element, law: InterfaceMaterial,
                 diameter: float, perimeter: float | None = None,
                 integration: str = "nodal"):
        name = type(self).__name__
        n_bar = len(bar.nodes)
        if n_bar not in (2, 3):
            raise ValueError(f"{name}(id={element_id}): la barra debe tener 2 o 3 nodos.")
        if integration not in INTEGRATIONS:
            raise ValueError(
                f"{name}(id={element_id}): integration debe ser uno de {INTEGRATIONS} "
                f"(recibido {integration!r}).")
        dim = sum(1 for d in bar.DOF_NAMES if d in TRANSLATIONS)
        if not isinstance(law, InterfaceMaterial) or not law.accepts("line", dim):
            raise ValueError(
                f"{name}(id={element_id}): la ley {type(law).__name__} no es un material "
                f"de interfaz de línea con {dim} componentes (GEOMETRY='line').")
        if not diameter > 0.0:
            raise ValueError(f"{name}(id={element_id}): diameter debe ser > 0.")
        perimeter = np.pi * diameter if perimeter is None else perimeter
        if not perimeter > 0.0:
            raise ValueError(f"{name}(id={element_id}): perimeter debe ser > 0.")

        self.DOF_NAMES = list(TRANSLATIONS[:dim])
        self.STRAIN_DIM = dim
        xi, w = _RULES[(n_bar, integration)]
        self.N_INTEGRATION_POINTS = len(xi)
        self.bar_element_id, self.host_element_id = bar.id, host.id
        self.diameter, self.perimeter = float(diameter), float(perimeter)
        self.integration = integration
        self.n_bar = n_bar
        super().__init__(element_id, list(bar.nodes) + list(host.nodes), law)

        X = np.array([np.pad(np.asarray(n.coordinates, float)[:dim],
                             (0, max(0, dim - len(n.coordinates))))
                      for n in bar.nodes[:2]])
        axis = X[1] - X[0]
        self.length = float(np.linalg.norm(axis))
        self.R = local_axes(axis / self.length)
        self.W = np.array([self.diameter] * (dim - 1) + [self.perimeter])
        self.points = X[0] + xi[:, None] * axis
        self.weights = w * self.length
        Nb = bar_shape_functions(n_bar, xi)
        I = np.eye(dim)
        B = []
        for g, x in enumerate(self.points):
            xi_h = natural_coordinates(host, x)
            if xi_h is None:
                raise ValueError(
                    f"{name}(id={element_id}): el punto {np.round(x, 12).tolist()} del "
                    f"tramo no se localiza en el anfitrión {host.id}.")
            Nh = shape_values(host, xi_h)
            B.append(np.hstack([np.kron(Nb[g], I), -np.kron(Nh, I)]))
        self._B = np.array(B)                       # (n_g, d, ndof)
        self._RB = np.einsum("ij,gjk->gik", self.R, self._B)

    # --- formulación ----------------------------------------------------------

    def compute_element_state(self, u_e: np.ndarray):
        u_e = np.asarray(u_e, dtype=np.float64)
        ndof = u_e.size
        K = np.zeros((ndof, ndof))
        F = np.zeros(ndof)
        for g in range(self.N_INTEGRATION_POINTS):
            RB = self._RB[g]
            t, D, new_state = self.material.compute_traction(RB @ u_e, self.state.vars[g])
            self.state.vars_trial[g] = new_state
            self.state.stresses_trial[g] = t
            WRB = self.W[:, None] * RB
            F += self.weights[g] * (WRB.T @ t)
            K += self.weights[g] * (WRB.T @ D @ RB)
        return K, F

    def compute_mass_matrix(self, lumping: str = "consistent") -> np.ndarray:
        """La interfaz no tiene masa: la de la barra y la del sólido están en
        sus elementos."""
        validate_lumping_kwarg(lumping, type(self).__name__)
        n = len(self.nodes) * len(self.DOF_NAMES)
        return np.zeros((n, n))

    # --- resultados ---------------------------------------------------------

    def interface_state(self, U_global: np.ndarray) -> list[InterfacePoint]:
        """Deslizamiento, tracciones y fuerza de adherencia en cada punto de
        integración para ``U_global``, desde el estado confirmado de la ley."""
        u_e = np.asarray(self.get_local_displacements(U_global), dtype=np.float64)
        out = []
        for g in range(self.N_INTEGRATION_POINTS):
            delta = self._RB[g] @ u_e
            t, _, _ = self.material.compute_traction(delta, self.state.vars[g])
            out.append(InterfacePoint(position=self.points[g].copy(), slip=float(delta[-1]),
                                      bond_stress=float(t[-1]),
                                      normal_opening=delta[:-1].copy(),
                                      normal_traction=t[:-1].copy(),
                                      bond_force=float(self.perimeter * t[-1])))
        return out
