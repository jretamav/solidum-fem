"""Elementos de línea rectos de dos nodos: celda de referencia, campo interior y
fuerza de cuerpo variable (ADR 0022, fase 4).

Lo comparten las armaduras y los marcos. Cada elemento aporta su
interpolación en ejes locales (``_line_interpolation``), el paso de sus DOFs
globales a los locales (``_to_local``) y, si reparte la carga en ejes
locales, los ejes (``_load_axes``) y el paso inverso (``_to_global``).
:class:`LineElementMixin` da con eso:

- ``REFERENCE_CELL = "line"`` y ``REFERENCE_NODES``, que usa el generador de
  mallas estructuradas;
- :meth:`~LineElementMixin.line_field`: desplazamientos, giros y deformaciones
  generalizadas en cualquier punto del elemento, en ejes locales; con
  :meth:`~LineElementMixin.section_rigidities`, la energía de deformación
  ``½ ∫ eᵀ·D·e dx`` que la rigidez del elemento integra exactamente;
- :meth:`~LineElementMixin._field_line_load`: el vector nodal consistente de
  una fuerza de cuerpo dada como campo ``b(x)``;
- la carga por unidad de longitud (ADR 0044): sus **muestras** en los puntos
  ``LINE_LOAD_POINTS`` (en ejes globales), el vector nodal consistente de unas
  muestras y las muestras de la fuerza de cuerpo y del peso propio, que viajan
  con la carga (``LoadVector``) para dar las fuerzas internas a lo largo del
  elemento (:func:`solidum.core.element_forces.forces_along`);
- las capacidades declaradas (``LINE_LOADS``, ``LINE_LOAD_LOCAL_AXES``,
  ``INTERNAL_FORCES_ALONG``) que leen el lector del YAML y la salida, y los ejes
  con que se expresan las fuerzas internas (:meth:`internal_force_axes`).
"""
from __future__ import annotations

from typing import ClassVar

import numpy as np

from solidum.math.distributed_loads import field_values

# Gauss-Legendre de 4 puntos en [0, 1]: exacto hasta grado 7, es decir
# N·q con N cúbica (Hermite) y una carga q de grado ≤ 4.
_LINE_XI, _LINE_W = np.polynomial.legendre.leggauss(4)
_LINE_XI = 0.5 * (_LINE_XI + 1.0)
_LINE_W = 0.5 * _LINE_W


def hermite_functions(xi: np.ndarray, L: float):
    """Funciones cúbicas de Hermite ``v = N₁v₁ + N₂θ₁ + N₃v₂ + N₄θ₂`` en
    ``x = ξ·L``, con sus derivadas primera y segunda respecto de ``x``, cada
    una ``(n, 4)``."""
    xi = np.asarray(xi, dtype=np.float64)
    N = np.column_stack([1.0 - 3.0 * xi**2 + 2.0 * xi**3, L * (xi - 2.0 * xi**2 + xi**3),
                         3.0 * xi**2 - 2.0 * xi**3, L * (xi**3 - xi**2)])
    dN = np.column_stack([6.0 * (xi**2 - xi) / L, 1.0 - 4.0 * xi + 3.0 * xi**2,
                          6.0 * (xi - xi**2) / L, 3.0 * xi**2 - 2.0 * xi])
    d2N = np.column_stack([(12.0 * xi - 6.0) / L**2, (6.0 * xi - 4.0) / L,
                           (6.0 - 12.0 * xi) / L**2, (6.0 * xi - 2.0) / L])
    return N, dN, d2N


def linear_tangent(material) -> float:
    """Módulo ``E`` de un material 1D elástico lineal: su tangente en el
    origen."""
    return float(material.compute_state(0.0, None)[1])


class LineElementMixin:
    """Elemento de línea recto de dos nodos (ADR 0022, fase 4).

    Contrato de las subclases
    -------------------------
    ``LINE_FIELDS`` : nombres de los desplazamientos y giros locales que da
        :meth:`line_field` (``"u"`` axial, ``"v"``, ``"w"`` transversales,
        ``"theta_x"``, ``"theta_y"``, ``"theta_z"`` giros); los
        desplazamientos van primero.
    ``LINE_STRAINS`` : nombres de las deformaciones generalizadas
        (``"eps"`` axial, ``"twist"`` = θx', ``"kappa_y"``, ``"kappa_z"``,
        ``"gamma_y"``), en el orden de :meth:`section_rigidities`.
    ``_line_interpolation(xi)`` : ``(N, B)`` de forma ``(n, k, m)`` y
        ``(n, s, m)``, con ``m`` los DOFs locales: campos y deformaciones en
        los puntos ``xi ∈ [0, 1]``.
    ``_to_local(u_e)`` : DOFs locales a partir del vector global del elemento.
    ``section_rigidities()`` : rigideces de la sección por deformación.

    Carga por unidad de longitud (ADR 0044)
    ---------------------------------------
    ``LINE_LOADS``, ``LINE_LOAD_LOCAL_AXES``, ``INTERNAL_FORCES_ALONG`` :
        capacidades declaradas: recibe ``line_loads``; la admite en ejes
        locales; da fuerzas internas a lo largo del elemento.
    ``LINE_LOAD_POINTS``, ``LINE_LOAD_WEIGHTS`` : abscisas ``ξ ∈ [0, 1]`` y
        pesos de la regla con que se muestrea e integra la carga (Gauss-Legendre
        de 4 puntos: exacta con la carga de grado ≤ 4 y Hermite).
    ``LINE_LOAD_BY_COMPONENT`` : ``True`` en las armaduras, que reparten cada
        componente global con sus funciones axiales; ``False`` en los marcos,
        que proyectan la carga a ``_load_axes()`` y la reparten con
        ``_line_interpolation``.
    ``_load_area()``, ``_linear_density()`` : el área de la fuerza de cuerpo y
        la densidad lineal del peso propio, declaradas una vez por elemento; de
        ellas salen las muestras (:meth:`body_load_samples`,
        :meth:`self_weight_samples`), que ningún elemento redefine.
    ``internal_force_axes(U)`` : ejes locales (filas, en globales) y
        longitudes de referencia y corriente de las fuerzas internas.
    """
    REFERENCE_CELL: ClassVar[str] = "line"
    REFERENCE_NODES: ClassVar[tuple] = ((-1.0,), (1.0,))
    LINE_FIELDS: ClassVar[tuple] = ()
    LINE_STRAINS: ClassVar[tuple] = ()
    LINE_LOADS: ClassVar[bool] = True
    LINE_LOAD_LOCAL_AXES: ClassVar[bool] = True
    INTERNAL_FORCES_ALONG: ClassVar[bool] = True
    LINE_LOAD_BY_COMPONENT: ClassVar[bool] = False
    LINE_LOAD_POINTS: ClassVar[tuple] = tuple(float(v) for v in _LINE_XI)
    LINE_LOAD_WEIGHTS: ClassVar[tuple] = tuple(float(v) for v in _LINE_W)

    def line_field(self, u_e: np.ndarray, xi) -> tuple[np.ndarray, np.ndarray]:
        """Campos ``(n, len(LINE_FIELDS))`` y deformaciones generalizadas
        ``(n, len(LINE_STRAINS))`` en los puntos ``xi ∈ [0, 1]`` (``ξ = 0`` en
        el nodo i), en ejes locales, con la interpolación del elemento y su
        vector global de DOFs ``u_e``."""
        N, B = self._line_interpolation(np.atleast_1d(np.asarray(xi, dtype=np.float64)))
        u = self._to_local(np.asarray(u_e, dtype=np.float64))
        return N @ u, B @ u

    def section_rigidities(self) -> np.ndarray:
        raise NotImplementedError

    def _line_ends(self) -> np.ndarray:
        """Coordenadas ``(2, d)`` de los nodos, con ``d`` las componentes de
        desplazamiento del elemento (las que faltan, cero)."""
        d = sum(1 for n in self.DOF_NAMES if n in ("ux", "uy", "uz"))
        X = np.zeros((2, d))
        for k, node in enumerate(self.nodes):
            c = np.asarray(node.coordinates, dtype=np.float64)[:d]
            X[k, :len(c)] = c
        return X

    # ------------------------------------------------------------------
    # Carga por unidad de longitud (ADR 0022 y 0044)
    # ------------------------------------------------------------------

    def _load_area(self) -> float:
        """Área por la que la fuerza de cuerpo ``b`` da la carga por unidad de
        longitud ``A·b``. Por omisión, ``A`` del elemento."""
        return float(self.A)

    def _linear_density(self) -> float:
        """Densidad lineal ``ρA`` del peso propio. Por omisión, ``ρ·A``."""
        return float(self.material.density) * self._load_area()

    def _sample_points(self) -> np.ndarray:
        """Posiciones ``(n, d)`` de ``LINE_LOAD_POINTS`` sobre la geometría de
        referencia."""
        X = self._line_ends()
        xi = np.asarray(self.LINE_LOAD_POINTS, dtype=np.float64)
        return X[0] + xi[:, None] * (X[1] - X[0])

    def _reference_axis(self) -> np.ndarray:
        """Eje unitario ``i → j`` de la geometría de referencia."""
        X = self._line_ends()
        return (X[1] - X[0]) / self.L0

    def _field_samples(self, value, name: str) -> np.ndarray:
        """Valores ``(n, d)`` de una constante o un campo ``value(x)`` en los
        puntos de muestreo."""
        x = self._sample_points()
        return np.array(field_values(value, x, x.shape[1], name), dtype=np.float64)

    def line_load_samples(self, q, *, local: bool = False, projected: bool = False) -> np.ndarray:
        """Carga por unidad de longitud del elemento en ``LINE_LOAD_POINTS``, en
        ejes globales, ``(n, d)``.

        ``q`` es constante (``d`` componentes) o un campo ``q(x)`` de la
        posición sobre la geometría de referencia. Con ``local=True`` sus
        componentes van en los ejes locales del elemento (``_load_axes``, los
        de la referencia: carga muerta). Con ``projected=True`` cada componente
        global ``q_k`` se da por unidad de longitud proyectada sobre el plano
        perpendicular a ``e_k`` y actúa ``q_k·√(1 − t_k²)`` por unidad de
        longitud del elemento (spec Frame2DEuler §14)."""
        if local and projected:
            raise ValueError(f"{type(self).__name__}(id={self.id}): la carga proyectada se da "
                             "en ejes globales, no en locales.")
        if local and not self.LINE_LOAD_LOCAL_AXES:
            raise ValueError(f"{type(self).__name__}(id={self.id}): no admite la carga en ejes "
                             "locales (el elemento no define ejes transversales).")
        values = self._field_samples(q, "de línea")
        if local:
            values = values @ np.asarray(self._load_axes(), dtype=np.float64)
        if projected:
            t = self._reference_axis()
            values = values * np.sqrt(np.clip(1.0 - t * t, 0.0, None))
        return values

    def body_load_samples(self, b) -> np.ndarray:
        """Muestras de la fuerza de cuerpo: ``A·b(x)`` (``_load_area``)."""
        return self._load_area() * self._field_samples(b, "de cuerpo")

    def self_weight_samples(self, g) -> np.ndarray:
        """Muestras del peso propio: ``ρA·g`` (``_linear_density``)."""
        return self._linear_density() * self._field_samples(g, "de peso propio")

    def line_load_vector(self, samples) -> np.ndarray:
        """Vector nodal consistente ``∫ Nᵀ·q ds`` de la carga representada por
        sus muestras ``(n, d)`` en ejes globales (la cúbica que las interpola
        con 4 puntos, spec Frame2DEuler §14), en ejes globales.

        Armaduras (``LINE_LOAD_BY_COMPONENT``): cada componente global con
        las funciones axiales de ``_line_interpolation`` (la componente
        perpendicular al eje es fuerza nodal pura). Marcos: la carga se
        proyecta sobre los ejes locales (``_load_axes``) y se reparte con la
        interpolación de los desplazamientos del elemento."""
        q_global = np.asarray(samples, dtype=np.float64)
        xi = np.asarray(self.LINE_LOAD_POINTS, dtype=np.float64)
        w = np.asarray(self.LINE_LOAD_WEIGHTS, dtype=np.float64)
        N, _ = self._line_interpolation(xi)
        if self.LINE_LOAD_BY_COMPONENT:
            Na = N[:, 0, :]
            return (self.L0 * (Na * w[:, None]).T @ q_global).reshape(-1)
        R = self._load_axes()
        q = q_global @ R.T
        f_local = np.einsum("g,gki,gk->i", w * self.L0, N[:, :R.shape[0], :], q)
        return self._to_global(f_local)

    def compute_line_load(self, q, *, local: bool = False, projected: bool = False) -> np.ndarray:
        """Vector nodal consistente de una carga por unidad de longitud
        (ADR 0044): ``line_load_vector(line_load_samples(q, ...))``."""
        return self.line_load_vector(self.line_load_samples(q, local=local, projected=projected))

    def _field_load_by_component(self, b, area: float) -> np.ndarray:
        """Vector nodal ``∫ N_a·A·b(x) dx`` de las armaduras con la fuerza de
        cuerpo como campo: :meth:`line_load_vector` de ``A·b(x_g)``."""
        return self.line_load_vector(area * self._field_samples(b, "de cuerpo"))

    def _field_line_load(self, b, area: float) -> np.ndarray:
        """Vector nodal ``∫ Nᵀ·q dx`` de la carga por unidad de longitud
        ``q = A·b(x)`` de los marcos: :meth:`line_load_vector` de ``A·b(x_g)``
        (mismas operaciones que antes del ADR 0044: bit a bit)."""
        return self.line_load_vector(area * self._field_samples(b, "de cuerpo"))

    # ------------------------------------------------------------------
    # Ejes de las fuerzas internas (ADR 0044 §1)
    # ------------------------------------------------------------------

    def internal_force_axes(self, U_global) -> tuple[np.ndarray, float, float]:
        """``(R, L0, l)``: ejes locales en que se expresan las fuerzas internas
        (filas, en globales; los que usa ``internal_forces_from_nodal``) y las
        longitudes de referencia y corriente. Por omisión, los ejes de carga de
        la referencia y ``l = L0`` (elementos geométricamente lineales)."""
        return np.asarray(self._load_axes(), dtype=np.float64), float(self.L0), float(self.L0)

    def internal_forces_along(self, forces, U_global, xi, samples=None):
        """Fuerzas internas en ``xi ∈ [0, 1]`` (ADR 0044 §2) por estática del
        tramo, a partir de las fuerzas de extremo ya calculadas ``forces``
        (:class:`~solidum.core.element_forces.ElementForces`, las del paso; no
        se reevalúa el elemento) y de las muestras de su carga en globales
        (``None``: sin carga a lo largo del elemento). ``U_global`` sólo da la
        geometría (los ejes de la cuerda en un corotacional)."""
        from solidum.core.element_forces import forces_along
        if not self.INTERNAL_FORCES_ALONG:
            raise ValueError(f"{type(self).__name__}(id={self.id}): no da fuerzas internas a lo "
                             "largo del elemento.")
        R, L0, l = self.internal_force_axes(U_global)
        local = None if samples is None else np.asarray(samples, dtype=np.float64) @ R.T
        return forces_along(forces, self.LINE_LOAD_POINTS, local, L0, l, xi)

    def perturbation_forces_along(self, forces, base_forces, U_base, dU, xi,
                                  samples=None, base_samples=None):
        """Fuerzas internas en ``xi`` del resultado de una perturbación lineal
        (ADR 0031, ADR 0044 §2): la linealización de la estática del tramo
        alrededor del estado base ``U_base``. ``forces``: las fuerzas de
        extremo de la perturbación (totales, en los ejes del estado base);
        ``samples``: la carga total; ``base_forces`` y ``base_samples``: las
        del estado base; ``dU``: el incremento. Ejes y longitud del estado
        base, más :meth:`_perturbation_geometric_terms`."""
        along = self.internal_forces_along(forces, U_base, xi, samples)
        terms = self._perturbation_geometric_terms(base_forces, U_base, dU, along.xi, base_samples)
        if not terms:
            return along
        comps = {c: v + terms[c] if c in terms else v for c, v in along.components.items()}
        return type(along)(kind=along.kind, xi=along.xi, s=along.s, components=comps)

    def _perturbation_geometric_terms(self, base_forces, U_base, dU, xi, base_samples) -> dict:
        """Términos geométricos de primer orden del diagrama de una
        perturbación, ``{componente: (m,)}``. Ninguno por omisión: en un
        elemento geométricamente lineal los ejes no dependen del estado, y la
        ``N`` de una armadura corrotacional, componente sobre la cuerda del
        estado base, no lleva término geométrico (spec Frame2DEulerCorot §11)."""
        return {}
