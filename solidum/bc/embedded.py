"""Refuerzo embebido con adherencia perfecta (spec ``EmbeddedReinforcement``,
ADR 0032).

La barra es una cadena de elementos 1D ordinarios (``Truss2D``,
``Truss3D``, sus variantes corrotacionales, cables) con su propia malla,
definida por una polilínea. Cada nodo de la barra queda atado al
desplazamiento interpolado del elemento de sólido que lo contiene,

    u(X_b) = Σ_i N_i(ξ_b)·u_i,

una restricción afín por traslación (ADR 0004). La polilínea se corta en los
bordes de los anfitriones para que cada tramo quede dentro de uno solo: con
anfitriones de interpolación lineal a lo largo de una recta la deformación de
la barra es la del sólido proyectada sobre su eje.

El área de concreto que ocupa la barra no se descuenta (decisión del
usuario, 2026-11-10; de Borst et al. 2012, ec. 6.136): la barra se superpone
al sólido.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Optional, Sequence

import numpy as np

from solidum.core.dof_names import TRANSLATIONS
from solidum.mesh.point_location import PointLocator, is_locatable, shape_values
from solidum.registry import ConstraintRegistry, ElementRegistry, MaterialRegistry

# Elemento 1D por omisión según la dimensión del sólido.
_DEFAULT_ELEMENT = {2: "Truss2D", 3: "Truss3D"}
# Un coeficiente de interpolación menor que esto es un cero de redondeo.
_ZERO_COEFF = 1.0e-14
# Tramos más cortos que esta fracción del tamaño del anfitrión más pequeño se
# funden con el vecino (un roce con un vértice).
_MIN_PIECE = 1.0e-9
# Avance, relativo al anfitrión más pequeño, con que se entra en el anfitrión
# siguiente al salir de uno.
_PROBE = 1.0e-6


@dataclass(frozen=True)
class EmbeddedBar:
    """Lo que creó un refuerzo embebido en el dominio.

    Attributes
    ----------
    node_ids
        Nodos de la barra, en el orden de la polilínea.
    element_ids
        Elementos 1D de la barra (``element_ids[k]`` une ``node_ids[k]`` y
        ``node_ids[k + 1]``).
    host_ids
        Elemento anfitrión al que está atado cada nodo.
    natural_coordinates
        ``ξ`` de cada nodo en su anfitrión, shape ``(n_nodos, ndim)``.
    length
        Longitud de la barra en la configuración de referencia.
    dofs
        Grados de libertad atados en cada nodo (``("ux", "uy"[, "uz"])``).
    """

    node_ids: tuple
    element_ids: tuple
    host_ids: tuple
    natural_coordinates: np.ndarray
    length: float
    dofs: tuple

    def bond_forces(self, result) -> np.ndarray:
        """Fuerzas de adherencia concentradas que el sólido ejerce sobre la
        barra en cada uno de sus nodos, shape ``(n_nodos, ndim)``.

        Son las fuerzas de las restricciones de embebido en los nodos de la
        barra (``result.constraint_forces_by_node``, ADR 0004 §3): con ellas y
        las cargas aplicadas en sus nodos, la barra está en equilibrio.
        """
        by_node = getattr(result, "constraint_forces_by_node", None)
        if by_node is None:
            raise ValueError(
                "bond_forces: el resultado no tiene constraint_forces_by_node "
                "(sólo lo da un análisis estático, SolveResult)."
            )
        return np.array([[by_node.get(n, {}).get(d, 0.0) for d in self.dofs]
                         for n in self.node_ids])


@ConstraintRegistry.register
class EmbeddedReinforcement:
    """Barra de refuerzo embebida en elementos de sólido 2D/3D.

    Parameters
    ----------
    points
        Polilínea de la barra en la configuración de referencia: al menos
        dos puntos de ``ndim`` coordenadas, todos dentro de los anfitriones.
    area
        Área de la barra ``A_s``.
    material
        Material 1D de la barra (instancia).
    element
        Nombre del elemento 1D registrado (``Truss2D`` por omisión en 2D,
        ``Truss3D`` en 3D); cualquier elemento de dos nodos con firma
        ``(id, nodos, material, A, **element_kwargs)``.
    hosts
        Ids de los elementos anfitriones admitidos. Por omisión, todos los
        sólidos isoparamétricos de la dimensión de la barra.
    subdivisions
        Partes iguales en que se divide cada tramo (1 por omisión).
    element_kwargs
        Argumentos adicionales del elemento 1D (p. ej. ``L0`` de un cable).
    """

    def __init__(self, points: Sequence[Sequence[float]], area: float, material,
                 element: Optional[str] = None, hosts: Optional[Iterable[int]] = None,
                 subdivisions: int = 1, element_kwargs: Optional[dict] = None):
        pts = np.asarray(points, dtype=float)
        if pts.ndim != 2 or pts.shape[0] < 2 or pts.shape[1] not in (2, 3):
            raise ValueError(
                "EmbeddedReinforcement: 'points' debe ser una lista de al menos dos "
                f"puntos de 2 o 3 coordenadas; se recibió la forma {pts.shape}."
            )
        if np.any(np.linalg.norm(np.diff(pts, axis=0), axis=1) == 0.0):
            raise ValueError("EmbeddedReinforcement: 'points' tiene puntos consecutivos repetidos.")
        if not (isinstance(area, (int, float)) and area > 0):
            raise ValueError(f"EmbeddedReinforcement: 'area' debe ser positiva; se recibió {area!r}.")
        if not (isinstance(subdivisions, int) and subdivisions >= 1):
            raise ValueError(
                f"EmbeddedReinforcement: 'subdivisions' debe ser un entero ≥ 1; se recibió "
                f"{subdivisions!r}.")
        self.points = pts
        self.ndim = pts.shape[1]
        self.area = float(area)
        self.material = material
        self.element = element or _DEFAULT_ELEMENT[self.ndim]
        ElementRegistry.get(self.element)          # error temprano si no existe
        self.hosts = None if hosts is None else [int(h) for h in hosts]
        self.subdivisions = subdivisions
        self.element_kwargs = dict(element_kwargs or {})

    @classmethod
    def from_yaml(cls, data: dict, resolve) -> "EmbeddedReinforcement":
        """Construye el refuerzo desde su entrada del bloque ``constraints``::

            constraints:
              - type: EmbeddedReinforcement
                points: [[0.0, 0.04], [2.0, 0.04]]
                area: 4.0e-4
                material: 2          # id de un material 1D de 'materials'
                element: Truss2D     # opcional
                hosts: [1, 2, 3]     # opcional: ids de los anfitriones
                subdivisions: 1      # opcional
        """
        data = dict(data)
        allowed = {"points", "area", "material", "element", "hosts", "subdivisions",
                   "element_kwargs"}
        unknown = sorted(set(data) - allowed)
        if unknown:
            raise ValueError(
                f"EmbeddedReinforcement: parámetros no aceptados {unknown}. "
                f"Admitidos: {sorted(allowed)}.")
        for key in ("points", "area", "material"):
            if key not in data:
                raise ValueError(f"EmbeddedReinforcement: falta '{key}'.")
        data["material"] = resolve(MaterialRegistry, data["material"], "material")
        return cls(**data)

    # ------------------------------------------------------------------
    # Aplicación al dominio
    # ------------------------------------------------------------------

    def apply(self, domain) -> EmbeddedBar:
        """Crea en ``domain`` los nodos, los elementos 1D y las restricciones
        de la barra. Si el dominio ya estaba numerado, lo renumera."""
        hosts = self._host_elements(domain)
        locator = PointLocator(hosts, self.ndim)
        h_min = float(locator.sizes.min())
        positions = self._bar_positions(locator, h_min)

        located = []
        for x in positions:
            hit = locator.locate(x)
            if hit is None:
                raise ValueError(
                    f"EmbeddedReinforcement: el punto {np.round(x, 12).tolist()} de la "
                    f"barra no está dentro de ningún elemento anfitrión.")
            located.append(hit)

        dofs = TRANSLATIONS[:self.ndim]
        element_cls = ElementRegistry.get(self.element)
        next_node = max(domain.nodes, default=0) + 1
        next_element = max(domain.elements, default=0) + 1
        node_ids = []
        for x in positions:
            node = domain.add_node(next_node, [float(c) for c in x])
            node_ids.append(node.id)
            next_node += 1
        element_ids = []
        for a, b in zip(node_ids[:-1], node_ids[1:]):
            el = element_cls(next_element, [domain.nodes[a], domain.nodes[b]], self.material,
                             A=self.area, **self.element_kwargs)
            missing = set(dofs) - set(el.DOF_NAMES)
            if missing or len(el.nodes) != 2:
                raise ValueError(
                    f"EmbeddedReinforcement: {self.element} no es un elemento 1D de dos "
                    f"nodos con traslaciones {list(dofs)}.")
            domain.add_element(el)
            element_ids.append(el.id)
            next_element += 1

        host_ids, xis = [], []
        for nid, (index, xi) in zip(node_ids, located):
            host = locator.elements[index]
            N = shape_values(host, xi)
            keep = np.abs(N) > _ZERO_COEFF
            masters = [host.nodes[i].id for i in np.nonzero(keep)[0]]
            coeffs = [float(c) for c in N[keep]]
            for d in dofs:
                domain.add_linear_constraint((nid, d), [(m, d) for m in masters], coeffs)
            host_ids.append(host.id)
            xis.append(xi)
        domain.generated_nodes.update(node_ids)
        if domain.total_dofs:
            domain.generate_equation_numbers()
        length = float(np.sum(np.linalg.norm(np.diff(self.points, axis=0), axis=1)))
        return EmbeddedBar(tuple(node_ids), tuple(element_ids), tuple(host_ids),
                           np.array(xis), length, dofs)

    def _host_elements(self, domain) -> list:
        if self.hosts is None:
            hosts = [e for e in domain.elements.values() if is_locatable(e, self.ndim)]
            if not hosts:
                raise ValueError(
                    f"EmbeddedReinforcement: el dominio no tiene elementos de sólido "
                    f"{self.ndim}D donde embeber la barra.")
            return hosts
        missing = [h for h in self.hosts if h not in domain.elements]
        if missing:
            raise ValueError(f"EmbeddedReinforcement: anfitriones inexistentes {missing}.")
        wrong = [h for h in self.hosts if not is_locatable(domain.elements[h], self.ndim)]
        if wrong:
            raise ValueError(
                f"EmbeddedReinforcement: los elementos {wrong} no son sólidos "
                f"isoparamétricos {self.ndim}D con funciones de forma.")
        return [domain.elements[h] for h in self.hosts]

    def _bar_positions(self, locator: PointLocator, h_min: float) -> np.ndarray:
        """Nodos de la barra: los vértices de la polilínea, sus cortes con los
        bordes de los anfitriones y las subdivisiones de cada tramo."""
        positions = [self.points[0]]
        for p, q in zip(self.points[:-1], self.points[1:]):
            params = self._crossings(locator, p, q, h_min)
            for t0, t1 in zip(params[:-1], params[1:]):
                for k in range(1, self.subdivisions + 1):
                    t = t0 + (t1 - t0) * k / self.subdivisions
                    positions.append(p + t * (q - p))
        return np.array(positions)

    def _crossings(self, locator: PointLocator, p, q, h_min: float) -> list:
        """Parámetros ``t ∈ [0, 1]`` de ``p + t·(q − p)`` donde la recta cambia
        de anfitrión, con ``0`` y ``1``."""
        d = q - p
        L = float(np.linalg.norm(d))
        probe = _PROBE * h_min / L
        min_piece = _MIN_PIECE * h_min / L
        params = [0.0]
        t = 0.0
        while t < 1.0 - min_piece:
            hit = locator.locate(p + min(t + probe, 1.0) * d)
            if hit is None:
                x = p + min(t + probe, 1.0) * d
                raise ValueError(
                    f"EmbeddedReinforcement: la barra sale de los anfitriones cerca de "
                    f"{np.round(x, 9).tolist()}.")
            index = hit[0]
            if locator.natural_coordinates(index, q) is not None:
                t_exit = 1.0
            else:
                lo, hi = min(t + probe, 1.0), 1.0
                while (hi - lo) * L > 1.0e-12 * h_min:
                    mid = 0.5 * (lo + hi)
                    if locator.natural_coordinates(index, p + mid * d) is not None:
                        lo = mid
                    else:
                        hi = mid
                t_exit = 0.5 * (lo + hi)
            if t_exit - params[-1] > min_piece:
                params.append(t_exit)
            t = t_exit
        params[-1] = 1.0
        return params


def embed_reinforcement(domain, points, area: float, material, *, element: Optional[str] = None,
                        hosts: Optional[Iterable[int]] = None, subdivisions: int = 1,
                        element_kwargs: Optional[dict] = None) -> EmbeddedBar:
    """Embebe una barra de refuerzo con adherencia perfecta en los sólidos de
    ``domain`` (spec ``EmbeddedReinforcement``). Ver
    :class:`EmbeddedReinforcement` para los parámetros; devuelve el
    :class:`EmbeddedBar` con los nodos, elementos y anfitriones creados."""
    return EmbeddedReinforcement(points, area, material, element=element, hosts=hosts,
                                 subdivisions=subdivisions,
                                 element_kwargs=element_kwargs).apply(domain)
