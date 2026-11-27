"""Bases de los elementos de interfaz (ADR 0043).

``LocalJumpElement`` es lo común a todo elemento que entrega a una ley de
interfaz el salto en ejes locales (§2): bucle por punto, anchos por
componente, masa nula y registro por punto; la heredan
``InterfaceElement`` (espesor nulo, abajo) y la adherencia barra–sólido
(``solidum.bc.bond_interface.BondSlipInterface``).

Un elemento de interfaz de espesor nulo une dos caras A y B por pares de nodos y entrega a un
material de interfaz de superficie (ADR 0033) el salto en los ejes locales de
la cara A. Formulación: spec ``docs/specs/Interface2D4.md``.

Orden de nodos de los elementos cohesivos de Abaqus en los lineales: primero
la cara A, después la B; en los cuadráticos, primero las esquinas y después
los intermedios. ``A_NODES`` da las posiciones de los nodos de la cara A en el
orden de sus funciones de forma y ``B_NODES[i]``, la del nodo de B que está
frente al ``i``-ésimo de A (en ``Interface2D4``, 4 frente a 1 y 3 frente a 2).

En cada punto ``ξ`` de la cara:

    g = g_0·e_1 + R·Δu,   Δu = Σ_i N_i·(u_B,par(i) − u_A,i),
    g_0 = n·(x_B − x_A),
    f_e = Σ_g m_g·Bᵀ·Rᵀ·t,   K_e = Σ_g m_g·Bᵀ·Rᵀ·T·R·B,

con ``R`` las filas ``n, s`` (2D) o ``n, s1, s2`` (3D) sobre la cara A y
``m_g`` la medida del punto: peso × jacobiano de la cara A (× espesor en 2D).
En la regla nodal, ``m_i = ∫N_i J`` (área tributaria consistente; decisión
del usuario del 2026-11-26) en las caras lineales; en la línea cuadrática,
``m_i = w_i·J(ξ_i)`` con los pesos de Simpson ``(1/3, 1/3, 4/3)`` (en un lado
recto es lo mismo; en uno curvo es la norma del vector de área tributaria
``∫N_i n J`` y reproduce una presión uniforme: decisión del usuario del
2026-11-26, Diálogo Q1 de la spec ``Interface2D6``).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, ClassVar

import numpy as np

from solidum.core.element import Element, validate_lumping_kwarg
from solidum.core.element_state import ElementState
from solidum.core.interface_material import InterfaceMaterial, component_names
from solidum.registry import InterfaceMaterialRegistry

INTEGRATIONS = ("nodal", "gauss")

# Tolerancia geométrica relativa al tamaño de la cara (decisión del usuario del
# 2026-11-26): una holgura inicial |g_0| ≤ tol·h se lleva a 0 (redondeo); una
# interpenetración inicial mayor es error. La misma tolerancia acota el
# desfase tangencial de un par.
GAP_TOLERANCE = 1.0e-6


# -- caras: funciones de forma, derivadas y reglas ---------------------------------

def _line(xi):
    x = xi[0]
    return np.array([0.5 * (1.0 - x), 0.5 * (1.0 + x)]), np.array([[-0.5], [0.5]])


def _line3(xi):
    # Nodos en ξ = −1, +1, 0 (esquinas primero, como la cara A de Interface2D6).
    x = xi[0]
    return (np.array([0.5 * x * (x - 1.0), 0.5 * x * (x + 1.0), 1.0 - x * x]),
            np.array([[x - 0.5], [x + 0.5], [-2.0 * x]]))


def _tri(xi):
    x, y = xi
    return (np.array([1.0 - x - y, x, y]),
            np.array([[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]]))


_QUAD_CORNERS = np.array([[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]])


def _quad(xi):
    x, y = xi
    c = _QUAD_CORNERS
    N = 0.25 * (1.0 + c[:, 0] * x) * (1.0 + c[:, 1] * y)
    dN = np.column_stack([0.25 * c[:, 0] * (1.0 + c[:, 1] * y),
                          0.25 * c[:, 1] * (1.0 + c[:, 0] * x)])
    return N, dN


def _gauss_1d(n):
    x, w = np.polynomial.legendre.leggauss(n)
    return x, w


def _gauss_quad(n):
    x, w = _gauss_1d(n)
    return (np.array([[a, b] for b in x for a in x]),
            np.array([wa * wb for wb in w for wa in w]))


def _line3_degenerate(XA) -> bool:
    """Un intermedio con el que la tangente se anula o se invierte (el nodo a
    un cuarto, p. ej.) deja un extremo sin peso o con peso negativo (spec
    ``Interface2D6``, «Advertencias»)."""
    a = [(_line3(np.array([x]))[1].T @ XA)[0] for x in (-1.0, 0.0, 1.0)]
    return not (a[0] @ a[1] > 0.0 and a[2] @ a[1] > 0.0)


@dataclass(frozen=True)
class FaceCell:
    """Tipo de cara de un elemento de interfaz de espesor nulo: todo lo que
    la base necesita de ella (ADR 0043 §1). Un tipo de cara nuevo es una fila
    de :data:`FACE_CELLS`; un elemento nuevo con un tipo existente sólo
    declara ``FACE_CELL`` y sus nodos.

    - ``shape(ξ) → (N, dN/dξ)``, con los nodos de la cara A en el orden de sus
      funciones de forma.
    - Regla de Gauss: ``gauss_points`` y ``gauss_weights`` (la medida del
      punto es ``w·J``).
    - Regla nodal: los puntos son los nodos (``nodes``) y cada uno ve sólo su
      par. Ejes en ``nodal_axes_points``; medida ``w_i·J`` en esos puntos con
      ``nodal_weights``, o, sin ellos, el área tributaria ``∫N_i J`` con la
      regla exacta ``tributary_rule``.
    - ``degenerate(X_A) → bool``: comprobación propia de la cara (además de
      «medida ≤ 0»).
    """
    shape: Callable
    n_nodes: int
    dim: int                     # dimensión del modelo (JUMP_DIM): 2 en líneas, 3 en superficies
    gauss_points: np.ndarray
    gauss_weights: np.ndarray
    nodes: np.ndarray
    nodal_axes_points: np.ndarray
    nodal_weights: np.ndarray | None = None
    tributary_rule: tuple | None = None
    degenerate: Callable | None = None
    degenerate_message: str = ""

    def points(self, integration: str) -> np.ndarray:
        return self.nodes if integration == "nodal" else self.gauss_points


_TRI_NODES = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
# Regla de 6 puntos de grado 4 (Strang y Fix) para ∫N_i J en el triángulo.
_TA, _TB = 0.445948490915965, 0.091576213509771
_TRI_FINE = (np.array([[_TA, _TA], [1 - 2 * _TA, _TA], [_TA, 1 - 2 * _TA],
                       [_TB, _TB], [1 - 2 * _TB, _TB], [_TB, 1 - 2 * _TB]]),
             np.array([0.223381589678011 / 2] * 3 + [0.109951743655322 / 2] * 3))
_LINE_NODES = np.array([[-1.0], [1.0]])
_LINE3_NODES = np.array([[-1.0], [1.0], [0.0]])

#: Tipos de cara (``FACE_CELL``) de los elementos de interfaz de espesor nulo.
FACE_CELLS: dict[str, FaceCell] = {
    # Línea de 2 nodos: ∫N_i J = L/2 (exacta con Gauss de 3 puntos).
    "line": FaceCell(_line, 2, 2, _gauss_1d(2)[0][:, None], _gauss_1d(2)[1], _LINE_NODES,
                     _LINE_NODES, tributary_rule=(_gauss_1d(3)[0][:, None], _gauss_1d(3)[1])),
    # Línea de 3 nodos (esquinas primero), posiblemente curva: Simpson con el
    # jacobiano del nodo, la norma del vector de área tributaria ∫N_i n J que
    # reproduce una presión uniforme (Diálogo Q1 de la spec Interface2D6).
    "line3": FaceCell(_line3, 3, 2, _gauss_1d(3)[0][:, None], _gauss_1d(3)[1], _LINE3_NODES,
                      _LINE3_NODES, nodal_weights=np.array([1 / 3, 1 / 3, 4 / 3]),
                      degenerate=_line3_degenerate,
                      degenerate_message="su tangente se anula o se invierte entre los "
                                         "extremos y el nodo intermedio"),
    # Triángulo (plano): ∫N_i J = A/3.
    "tri": FaceCell(_tri, 3, 3, np.array([[1 / 6, 1 / 6], [2 / 3, 1 / 6], [1 / 6, 2 / 3]]),
                    np.full(3, 1 / 6), _TRI_NODES, _TRI_NODES, tributary_rule=_TRI_FINE),
    # Cuadrilátero bilineal, posiblemente alabeado: el nodo i toma la normal y
    # la medida del vector de área tributaria consistente ∫N_i (a_ξ × a_η)
    # dξdη, que es exactamente (a_ξ × a_η) en (ξ_i/3, η_i/3). En una cara plana
    # es ∫N_i J (J lineal); en una alabeada, la normal que equilibra a los
    # sólidos (verificación del 2026-11-26, decisión del usuario; spec
    # Interface3D8).
    "quad": FaceCell(_quad, 4, 3, *_gauss_quad(2), _QUAD_CORNERS, _QUAD_CORNERS / 3.0,
                     nodal_weights=np.ones(4)),
}


@dataclass(frozen=True)
class InterfacePointState:
    """Estado de un punto de integración de un elemento de interfaz."""
    position: np.ndarray        # coordenadas del punto (sobre la cara A, o sobre la barra)
    axes: np.ndarray            # R: filas en el orden de la ley (n, s / n, s1, s2 / n1, n2, s)
    jump: np.ndarray            # g en ejes locales (g_n con la holgura inicial)
    traction: np.ndarray        # t en ejes locales [Pa]
    state: dict                 # estado de la ley (confirmado)


class LocalJumpElement(Element):
    """Base común de los elementos cuya magnitud por punto es el salto en ejes
    locales que reciben de una ley de interfaz (ADR 0043 §2): los de espesor
    nulo (:class:`InterfaceElement`) y la adherencia barra–sólido
    (``BondSlipInterface``).

    Lo común: el bucle por punto, los ejes ``R`` y la matriz ``R·B`` de cada
    punto, los anchos ``W`` por componente, la masa nula y un registro
    genérico por punto. En cada punto ``g``, con la medida ``m_g`` (peso ×
    jacobiano, o longitud tributaria) y ``W = diag(w_1, …, w_d)``::

        f_e = Σ_g m_g·(W·R·B)ᵀ·t,   K_e = Σ_g m_g·(W·R·B)ᵀ·T·R·B.

    ``W`` es la identidad en una superficie (la tracción actúa sobre el área
    de la cara) y ``diag(d_b, [d_b], p)`` en una barra (la normal sobre el
    diámetro, la de adherencia sobre el perímetro).

    Contrato de las subclases: ``JUMP_DIM`` (de clase o de instancia),
    ``N_INTEGRATION_POINTS``, la geometría por punto con
    :meth:`_set_point_geometry` y el salto con :meth:`jumps`.
    """

    STRAIN_DIM = None
    JUMP_DIM: ClassVar[int | None] = None

    def _init_state(self) -> None:
        self.state = ElementState(self.N_INTEGRATION_POINTS, init_stress=np.zeros(self.JUMP_DIM))

    def _set_point_geometry(self, R, RB, measures, points, widths=None) -> None:
        """Fija la geometría de los puntos: ejes ``R`` ``(n_g, d, d)``, ``R·B``
        ``(n_g, d, ndof)``, medidas ``(n_g,)``, posiciones ``(n_g, d)`` y
        anchos ``(d,)`` (``None``: todos 1)."""
        self.R = np.asarray(R, dtype=np.float64)
        self._RB = np.asarray(RB, dtype=np.float64)
        self.measures = np.asarray(measures, dtype=np.float64)
        # Posiciones físicas (en los sólidos, ``points`` son coordenadas
        # naturales: otro nombre para no confundirlas).
        self.positions = np.asarray(points, dtype=np.float64)
        self._WRB = (self._RB if widths is None
                     else np.asarray(widths, dtype=np.float64)[None, :, None] * self._RB)

    def jumps(self, u_e: np.ndarray) -> np.ndarray:
        """Salto local ``(n_g, d)`` en cada punto para ``u_e``."""
        raise NotImplementedError(f"{type(self).__name__} no define jumps().")

    def jump_components(self) -> tuple[str, ...]:
        """Nombres de las componentes locales del salto y de la tracción, en
        su orden (``("n", "s")``, ``("n", "s1", "s2")`` o ``("n1", "n2",
        "s")``): los de la salida (ADR 0043 §7)."""
        return component_names(self.material.GEOMETRY, self.JUMP_DIM)

    def compute_element_state(self, u_e: np.ndarray):
        u_e = np.asarray(u_e, dtype=np.float64)
        ndof = u_e.size
        K = np.zeros((ndof, ndof))
        F = np.zeros(ndof)
        for g, jump in enumerate(self.jumps(u_e)):
            t, T, new_state = self.material.compute_traction(jump, self.state.vars[g])
            self.state.vars_trial[g] = new_state
            self.state.stresses_trial[g] = t
            WRB = self._WRB[g]
            F += self.measures[g] * (WRB.T @ t)
            K += self.measures[g] * (WRB.T @ T @ self._RB[g])
        return K, F

    def compute_mass_matrix(self, lumping: str = "consistent") -> np.ndarray:
        """La interfaz no tiene masa: la de cada lado está en sus elementos."""
        validate_lumping_kwarg(lumping, type(self).__name__)
        n = len(self.nodes) * len(self.DOF_NAMES)
        return np.zeros((n, n))

    def force_matrix(self) -> np.ndarray:
        """``G_e`` ``(ndof, n_g·d)`` con ``f_e = G_e·[t_1; …; t_ng]``: el
        bloque del punto ``g`` es ``m_g·(W·R·B)_gᵀ``. Fija en desplazamientos
        pequeños; la usan las energías de la salida (ADR 0038 y 0043 §5)."""
        n_g, d, ndof = self._WRB.shape
        return (self.measures[:, None, None] * self._WRB).transpose(2, 0, 1).reshape(ndof, n_g * d)

    def nodal_forces(self, trial: bool = False) -> np.ndarray:
        """Fuerzas nodales ``f_e`` de las tracciones del estado confirmado (o
        de prueba), sin evaluar la ley."""
        s = self.state.stresses_trial if trial else self.state.stresses
        return self.force_matrix() @ np.concatenate([np.atleast_1d(np.asarray(t, dtype=float))
                                                     for t in s])

    def interface_points(self, U_global: np.ndarray) -> list[InterfacePointState]:
        """Salto, tracción y estado de la ley en cada punto, desde el estado
        confirmado."""
        u_e = np.asarray(self.get_local_displacements(U_global), dtype=np.float64)
        out = []
        for g, jump in enumerate(self.jumps(u_e)):
            t, _, _ = self.material.compute_traction(jump, self.state.vars[g])
            out.append(InterfacePointState(position=self.positions[g].copy(), axes=self.R[g].copy(),
                                           jump=np.array(jump), traction=np.array(t),
                                           state=dict(self.state.vars[g] or {})))
        return out


class InterfaceElement(LocalJumpElement):
    """Base de los elementos de interfaz de espesor nulo (ADR 0043).

    Contrato de las subclases: ``DOF_NAMES``, ``JUMP_DIM`` (2 o 3),
    ``FACE_CELL`` (un tipo de :data:`FACE_CELLS`: ``"line"``, ``"line3"``,
    ``"tri"``, ``"quad"``),
    ``N_FACE_NODES``, ``A_NODES`` y ``B_NODES`` (posiciones, en la lista de
    nodos, de los nodos de la cara A en el orden de sus funciones de forma y
    de su par en B), ``N_INTEGRATION_POINTS`` y ``VTK_CELL_TYPE``.

    Parameters
    ----------
    element_id : int
    nodes : list[Node]
        Cara A y después cara B, en el orden de la clase.
    material : InterfaceMaterial
        Ley de interfaz de superficie con ``JUMP_DIM`` componentes.
    integration : {"nodal", "gauss"}
        Por omisión ``"nodal"`` (decisión del usuario del 2026-11-26).
    thickness : float
        Espesor fuera del plano (sólo 2D).
    """

    JUMP_DIM: ClassVar[int]
    FACE_CELL: ClassVar[str]
    N_FACE_NODES: ClassVar[int]
    A_NODES: ClassVar[tuple]
    B_NODES: ClassVar[tuple]
    VTK_CELL_TYPE: ClassVar[str]
    ACCEPTS_UNILATERAL = True
    REFERENCE_KWARGS: ClassVar[dict] = {"material": InterfaceMaterialRegistry}

    def __init_subclass__(cls, **kwargs):
        """Contrato de la subclase, comprobado al definirla: un ``FACE_CELL``
        mal escrito o unos nodos que no casan con la cara fallan al importar,
        no al construir el primer elemento."""
        super().__init_subclass__(**kwargs)
        cell = getattr(cls, "FACE_CELL", None)
        if cell is None:
            return
        name = cls.__name__
        if cell not in FACE_CELLS:
            raise TypeError(f"{name}: FACE_CELL={cell!r} no es un tipo de cara conocido "
                            f"({sorted(FACE_CELLS)}).")
        face = FACE_CELLS[cell]
        if getattr(cls, "JUMP_DIM", None) != face.dim:
            raise TypeError(f"{name}: JUMP_DIM={getattr(cls, 'JUMP_DIM', None)!r}, pero una cara "
                            f"{cell!r} es de un modelo {face.dim}D.")
        n = getattr(cls, "N_FACE_NODES", None)
        if n != face.n_nodes:
            raise TypeError(f"{name}: N_FACE_NODES={n!r}, pero una cara {cell!r} tiene "
                            f"{face.n_nodes} nodos.")
        a, b = tuple(getattr(cls, "A_NODES", ())), tuple(getattr(cls, "B_NODES", ()))
        if len(a) != n or len(b) != n or sorted(a + b) != list(range(2 * n)):
            raise TypeError(f"{name}: A_NODES={a} y B_NODES={b} deben repartir las "
                            f"posiciones 0…{2 * n - 1} entre las dos caras, {n} en cada una.")
        n_g = {len(face.points(i)) for i in INTEGRATIONS}
        if n_g != {getattr(cls, "N_INTEGRATION_POINTS", None)}:
            raise TypeError(f"{name}: N_INTEGRATION_POINTS={getattr(cls, 'N_INTEGRATION_POINTS', None)!r}, "
                            f"pero las reglas de la cara {cell!r} tienen {sorted(n_g)} puntos.")

    def __init__(self, element_id: int, nodes, material: InterfaceMaterial,
                 integration: str = "nodal", thickness: float = 1.0):
        name = type(self).__name__
        dim = self.JUMP_DIM
        if len(nodes) != 2 * self.N_FACE_NODES:
            raise ValueError(f"{name}(id={element_id}): necesita {2 * self.N_FACE_NODES} "
                             f"nodos (cara A y cara B); tiene {len(nodes)}.")
        if not isinstance(material, InterfaceMaterial) or not material.accepts("surface", dim):
            raise ValueError(
                f"{name}(id={element_id}): la ley {type(material).__name__} no es un "
                f"material de interfaz de superficie con {dim} componentes "
                f"(interface_materials con GEOMETRY='surface').")
        if integration not in INTEGRATIONS:
            raise ValueError(f"{name}(id={element_id}): integration debe ser uno de "
                             f"{INTEGRATIONS} (recibido {integration!r}).")
        if dim == 2 and not thickness > 0.0:
            raise ValueError(f"{name}(id={element_id}): thickness debe ser > 0.")
        self.integration = integration
        self.thickness = float(thickness) if dim == 2 else 1.0
        super().__init__(element_id, list(nodes), material)
        self._build_geometry()

    # --- geometría ------------------------------------------------------------

    def _face_coordinates(self, X=None):
        """Coordenadas ``(X_A, X_B)`` de las dos caras, la B en el orden de A."""
        d = self.JUMP_DIM
        if X is None:
            X = np.array([np.pad(np.asarray(n.coordinates, float)[:d],
                                 (0, max(0, d - len(n.coordinates)))) for n in self.nodes])
        return X[list(self.A_NODES)], X[list(self.B_NODES)]

    def local_axes(self, XA, xi) -> tuple[np.ndarray, float]:
        """Ejes locales ``R`` y jacobiano de la cara A en ``ξ``. Recibe las
        coordenadas para que una variante con desplazamientos grandes lo
        llame con ``X + u`` (ADR 0043 §1)."""
        _, dN = FACE_CELLS[self.FACE_CELL].shape(xi)
        a = dN.T @ XA                           # (d_cara, d)
        if self.JUMP_DIM == 2:
            J = float(np.linalg.norm(a[0]))
            s = a[0] / J
            n = np.array([-s[1], s[0]])         # e_z × s
            return np.array([n, s]), J
        cross = np.cross(a[0], a[1])
        J = float(np.linalg.norm(cross))
        n = cross / J
        s1 = a[0] / np.linalg.norm(a[0])
        return np.array([n, s1, np.cross(n, s1)]), J

    def _build_geometry(self) -> None:
        name = f"{type(self).__name__}(id={self.id})"
        d, nf = self.JUMP_DIM, self.N_FACE_NODES
        face = FACE_CELLS[self.FACE_CELL]
        XA, XB = self._face_coordinates()
        shape = face.shape
        if face.degenerate is not None and face.degenerate(XA):
            raise ValueError(f"{name}: la cara A es degenerada ({face.degenerate_message}).")
        pts = face.points(self.integration)

        # Ejes y medida de cada punto: en Gauss, en el punto, con w·J; en la
        # regla nodal, donde y como diga el tipo de cara (FaceCell).
        if self.integration == "nodal":
            axes_pts = face.nodal_axes_points
            if face.nodal_weights is not None:
                measure = np.array([wq * self.local_axes(XA, xi)[1]
                                    for xi, wq in zip(axes_pts, face.nodal_weights)])
            else:
                measure = np.zeros(nf)
                for xi, wq in zip(*face.tributary_rule):
                    N, _ = shape(xi)
                    measure += wq * self.local_axes(XA, xi)[1] * N
        else:
            axes_pts = pts
            measure = np.array([wq * self.local_axes(XA, xi)[1]
                                for xi, wq in zip(pts, face.gauss_weights)])
        if not np.all(measure > 0.0):
            raise ValueError(f"{name}: la cara A es degenerada (medida ≤ 0).")
        self.size = float(np.sum(measure)) ** (1.0 / (d - 1))   # h: longitud o √área

        tol = GAP_TOLERANCE * self.size
        R_all, N_all, g0_all, pos = [], [], [], []
        for xi, xi_axes in zip(pts, axes_pts):
            N, _ = shape(xi)
            R, _ = self.local_axes(XA, xi_axes)
            gap = R @ (N @ XB - N @ XA)
            if np.any(np.abs(gap[1:]) > tol):
                raise ValueError(
                    f"{name}: el par no está sobre la normal de la cara A (desfase "
                    f"tangencial {np.max(np.abs(gap[1:])):.3g} > {tol:.3g}).")
            g0 = float(gap[0])
            if g0 < -tol:
                raise ValueError(
                    f"{name}: las caras se interpenetran al inicio (holgura inicial "
                    f"{g0:.6g}, tolerancia {tol:.3g}).")
            if abs(g0) <= tol:
                g0 = 0.0
            R_all.append(R)
            N_all.append(N)
            g0_all.append(g0)
            pos.append(N @ XA)
        self.g0 = np.array(g0_all)
        if np.any(self.g0 != 0.0) and not self.material.ACCEPTS_INITIAL_GAP:
            raise ValueError(
                f"{name}: la ley {type(self.material).__name__} no admite holgura inicial "
                f"y las caras están separadas (holgura máxima {np.max(self.g0):.6g}).")
        R_all = np.array(R_all)                             # (n_g, d, d)

        # B: (n_g, d, ndof) con u_e = [u_A..., u_B...].
        ndof = 2 * nf * d
        I = np.eye(d)
        B = np.zeros((len(pts), d, ndof))
        for g, N in enumerate(N_all):
            for i in range(nf):
                a = self.A_NODES[i] * d
                b = self.B_NODES[i] * d
                B[g, :, a:a + d] -= N[i] * I
                B[g, :, b:b + d] += N[i] * I
        self._set_point_geometry(R_all, np.einsum("gij,gjk->gik", R_all, B),
                                 measure * self.thickness, np.array(pos))
        # Para el salto, diferencias por par antes de interpolar: una traslación
        # rígida da salto exactamente nulo, sin cancelación de desplazamientos
        # grandes frente a holguras pequeñas.
        self._N = np.array(N_all)                           # (n_g, n_cara)
        dofs_A = np.asarray(self.A_NODES)[:, None] * d + np.arange(d)
        dofs_B = np.asarray(self.B_NODES)[:, None] * d + np.arange(d)
        self._pair_dofs = (dofs_A, dofs_B)

    # --- formulación ----------------------------------------------------------

    def jumps(self, u_e: np.ndarray) -> np.ndarray:
        """Salto local ``g`` ``(n_g, d)`` en cada punto para ``u_e``."""
        u_e = np.asarray(u_e, dtype=np.float64)
        dofs_A, dofs_B = self._pair_dofs
        du = u_e[dofs_B] - u_e[dofs_A]                      # (n_cara, d)
        g = np.einsum("gij,gj->gi", self.R, self._N @ du)
        g[:, 0] += self.g0
        return g

    # --- resultados -------------------------------------------------------------

    def side_a_displacements(self, u_e: np.ndarray) -> np.ndarray:
        """Desplazamiento ``(n_g, d)`` de la cara A en cada punto (los puntos
        de ``points`` están sobre ella): el que deforma la nube de puntos de
        la salida."""
        u_e = np.asarray(u_e, dtype=np.float64)
        return self._N @ u_e[self._pair_dofs[0]]

    def stored_energy(self, U_global: np.ndarray, trial: bool = False) -> float:
        """Energía almacenada ``Σ_g m_g·ψ_g`` con el estado confirmado (o de
        prueba) (ADR 0043 §5). ``ψ`` es la de la ley, con la holgura actual en
        ``g_n``: con holgura inicial o con fricción no es ``½·u_e·f_e``. Lanza
        ``NotImplementedError`` si la ley no la conoce."""
        u_e = np.asarray(self.get_local_displacements(U_global), dtype=np.float64)
        vars_ = self.state.vars_trial if trial else self.state.vars
        return float(sum(self.measures[g] * self.material.stored_energy(jump, vars_[g])
                         for g, jump in enumerate(self.jumps(u_e))))
