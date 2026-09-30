# solidum_fem/solidum/sections/fiber.py
"""Sección por fibras para los marcos Euler-Bernoulli (ADR 0021, spec
``docs/specs/FiberSection.md``).

La sección se discretiza en **fibras** ``(y_i, z_i, A_i, material_i)``,
con coordenadas en los ejes locales del elemento medidas desde el **eje del
elemento** (la línea de nodos), no desde el centroide. Cada fibra es un punto
material uniaxial con su propio material 1D del catálogo (``STRAIN_DIM = 1``)
y su propia historia, así que la plastificación o el daño avanzan fibra a
fibra. Una sección puede mezclar materiales (concreto y acero).

Cinemática (spec §2, Reglas.md §5: ``M_z > 0`` tracciona ``y < 0``;
``M_y > 0`` tracciona ``z > 0``)::

    2D:  ε_i = ε₀ − y_i·κ_z                 e = [ε₀, κ_z]
    3D:  ε_i = ε₀ − y_i·κ_z + z_i·κ_y       e = [ε₀, κ_y, κ_z]

Fuerzas internas y tangente (spec §4), regla del punto medio sobre las
fibras::

    s   = Σ σ_i A_i w_i          w_i = [1, −y_i]        (2D)
    k_s = Σ E_t,i A_i w_i w_iᵀ   w_i = [1, z_i, −y_i]   (3D)

``k_s = ∂s/∂e`` es exacta respecto de esa aproximación y simétrica para
cualquier ``E_t,i``. En forma matricial, con ``W`` la matriz de filas
``w_iᵀ``: ``ε = W·e``, ``s = Wᵀ(σ∘A)``, ``k_s = Wᵀ diag(E_t∘A) W``.

Primitivas de discretización (con la API de OpenSees: *patch*, *layer*,
*fiber*):

- :class:`RectangularPatch` — rejilla ``ny × nz`` de celdas; una fibra en
  el centro de cada celda con el área de la celda.
- :class:`CircularPatch` — corona circular (maciza con ``r_inner = 0``)
  dividida en ``n_rad × n_circ`` sectores anulares; una fibra en el
  centroide de cada sector con su área exacta. Los ángulos se miden desde
  el eje ``+y`` local hacia ``+z``; el primer sector empieza en el ángulo 0.
- :class:`BarLayer` — ``n`` barras de área ``area`` equiespaciadas de
  ``start`` a ``end`` inclusive (``n = 1``: una barra en ``start``).
- :class:`Fiber` — una fibra suelta.

La sección evalúa sus fibras **agrupadas por material** (ADR 0021 §6): el
bucle sobre las fibras de un mismo material reutiliza el método ligado del
material y deja la puerta abierta a una evaluación vectorizada por grupo
cuando la medición lo justifique.
"""
from __future__ import annotations

import inspect
import math
from typing import Tuple

import numpy as np

from solidum.core.section import Section
from solidum.registry import MaterialRegistry, SectionRegistry


# ---------------------------------------------------------------------------
# Validación de datos de entrada (mensajes que nombran la primitiva)
# ---------------------------------------------------------------------------

def _check_material(material, where: str) -> None:
    """Material de fibra: un objeto con ``compute_state`` y ``STRAIN_DIM = 1``."""
    if material is None or not callable(getattr(material, "compute_state", None)):
        raise ValueError(
            f"{where}: 'material' debe ser un material con compute_state; se "
            f"recibió {type(material).__name__}."
        )
    dim = getattr(material, "STRAIN_DIM", None)
    if dim != 1:
        raise ValueError(
            f"{where}: el material {type(material).__name__} declara "
            f"STRAIN_DIM={dim}; una fibra es uniaxial (sólo σ_xx) y requiere un "
            f"material 1D (STRAIN_DIM = 1), p. ej. Elastic1D, Elastoplastic1D o "
            f"IsotropicDamage1D."
        )


def _as_count(value, name: str, where: str) -> int:
    """Número de divisiones: entero ≥ 1."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise ValueError(f"{where}: {name}={value!r} debe ser un entero ≥ 1.")
    if value < 1:
        raise ValueError(f"{where}: {name}={value} debe ser ≥ 1.")
    return int(value)


def _as_number(value, name: str, where: str) -> float:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{where}: {name}={value!r} debe ser un número.")
    try:
        v = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{where}: {name}={value!r} debe ser un número.") from None
    if not math.isfinite(v):
        raise ValueError(f"{where}: {name}={value!r} debe ser finito.")
    return v


def _as_positive(value, name: str, where: str) -> float:
    v = _as_number(value, name, where)
    if v <= 0.0:
        raise ValueError(f"{where}: {name}={value!r} debe ser > 0.")
    return v


def _as_point(value, name: str, where: str) -> Tuple[float, float]:
    """Par ``(y, z)`` de números finitos."""
    try:
        seq = list(value)
    except TypeError:
        raise ValueError(f"{where}: {name}={value!r} debe ser un par (y, z).") from None
    if len(seq) != 2:
        raise ValueError(f"{where}: {name}={value!r} debe ser un par (y, z).")
    return (_as_number(seq[0], f"{name}[0]", where),
            _as_number(seq[1], f"{name}[1]", where))


def _as_interval(value, name: str, where: str) -> Tuple[float, float]:
    """Par ``(mín, máx)`` con ``mín < máx``."""
    try:
        seq = list(value)
    except TypeError:
        raise ValueError(f"{where}: {name}={value!r} debe ser un par (mín, máx).") from None
    if len(seq) != 2:
        raise ValueError(f"{where}: {name}={value!r} debe ser un par (mín, máx).")
    lo = _as_number(seq[0], f"{name}[0]", where)
    hi = _as_number(seq[1], f"{name}[1]", where)
    if not hi > lo:
        raise ValueError(f"{where}: {name}=({lo}, {hi}) requiere {name}[1] > {name}[0].")
    return lo, hi


# ---------------------------------------------------------------------------
# Primitivas de discretización
# ---------------------------------------------------------------------------

class RectangularPatch:
    """Parche rectangular ``[y_min, y_max] × [z_min, z_max]`` dividido en una
    rejilla de ``ny × nz`` celdas iguales; una fibra en el centro de cada
    celda, con el área de la celda.

    Con ``nz = 1`` (o ``ny = 1``) el parche es una pila de capas: basta para
    la flexión en un solo plano. Una sección rectangular de ``n`` capas
    integra ``EI = E·b·h³/12·(1 − 1/n²)`` (regla del punto medio).
    """

    def __init__(self, material, y, z, ny, nz):
        where = "RectangularPatch"
        _check_material(material, where)
        self.material = material
        self.y = _as_interval(y, "y", where)
        self.z = _as_interval(z, "z", where)
        self.ny = _as_count(ny, "ny", where)
        self.nz = _as_count(nz, "nz", where)

    def discretize(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """``(y, z, A)`` de las fibras del parche."""
        (y0, y1), (z0, z1) = self.y, self.z
        dy = (y1 - y0) / self.ny
        dz = (z1 - z0) / self.nz
        yc = y0 + (np.arange(self.ny) + 0.5) * dy
        zc = z0 + (np.arange(self.nz) + 0.5) * dz
        Y, Z = np.meshgrid(yc, zc, indexing="ij")
        n = self.ny * self.nz
        return Y.ravel(), Z.ravel(), np.full(n, dy * dz)


class CircularPatch:
    """Parche circular: corona ``r_inner ≤ r ≤ r_outer`` (maciza con
    ``r_inner = 0``) centrada en ``center``, dividida en ``n_rad`` anillos
    de igual espesor y ``n_circ`` sectores de igual ángulo.

    Una fibra por sector anular, en su **centroide** y con su **área
    exacta** ``A = Δθ/2·(r₂² − r₁²)``. El centroide dista del centro
    ``r_c = (2/3)·(sen α/α)·(r₂³ − r₁³)/(r₂² − r₁²)``, con ``α = Δθ/2``, a lo
    largo de la bisectriz del sector. Los ángulos se miden desde el eje
    ``+y`` local hacia ``+z``, y el primer sector empieza en el ángulo 0.
    """

    def __init__(self, material, r_outer, r_inner=0.0, center=(0.0, 0.0),
                 n_rad=None, n_circ=None):
        where = "CircularPatch"
        _check_material(material, where)
        self.material = material
        self.r_outer = _as_positive(r_outer, "r_outer", where)
        self.r_inner = _as_number(r_inner, "r_inner", where)
        if self.r_inner < 0.0 or not self.r_outer > self.r_inner:
            raise ValueError(
                f"{where}: se requiere r_outer > r_inner ≥ 0 (r_outer="
                f"{self.r_outer}, r_inner={self.r_inner})."
            )
        self.center = _as_point(center, "center", where)
        if n_rad is None or n_circ is None:
            faltan = [k for k, v in (("n_rad", n_rad), ("n_circ", n_circ)) if v is None]
            raise ValueError(f"{where}: faltan las divisiones {faltan} (enteros ≥ 1).")
        self.n_rad = _as_count(n_rad, "n_rad", where)
        self.n_circ = _as_count(n_circ, "n_circ", where)

    def discretize(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """``(y, z, A)`` de las fibras del parche."""
        r = np.linspace(self.r_inner, self.r_outer, self.n_rad + 1)
        r1, r2 = r[:-1], r[1:]
        dtheta = 2.0 * math.pi / self.n_circ
        alpha = 0.5 * dtheta
        # sen α/α con α = π (n_circ = 1) vale 0: la corona entera tiene el
        # centroide en el centro, que es lo correcto.
        ratio = math.sin(alpha) / alpha
        r_c = (2.0 / 3.0) * ratio * (r2**3 - r1**3) / (r2**2 - r1**2)
        area = alpha * (r2**2 - r1**2)
        theta = (np.arange(self.n_circ) + 0.5) * dtheta
        RC, TH = np.meshgrid(r_c, theta, indexing="ij")
        A = np.repeat(area, self.n_circ)
        cy, cz = self.center
        return (cy + (RC * np.cos(TH)).ravel(),
                cz + (RC * np.sin(TH)).ravel(),
                A)


class BarLayer:
    """Capa de ``n`` barras iguales de área ``area`` equiespaciadas sobre el
    segmento de ``start`` a ``end``, ambos extremos incluidos. Con ``n = 1``
    hay una sola barra, en ``start``."""

    def __init__(self, material, start, end, n, area):
        where = "BarLayer"
        _check_material(material, where)
        self.material = material
        self.start = _as_point(start, "start", where)
        self.end = _as_point(end, "end", where)
        self.n = _as_count(n, "n", where)
        self.area = _as_positive(area, "area", where)

    def discretize(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """``(y, z, A)`` de las barras."""
        t = np.zeros(1) if self.n == 1 else np.linspace(0.0, 1.0, self.n)
        (y0, z0), (y1, z1) = self.start, self.end
        return y0 + t * (y1 - y0), z0 + t * (z1 - z0), np.full(self.n, self.area)


class Fiber:
    """Fibra suelta en ``(y, z)`` con área ``area``."""

    def __init__(self, y, z, area, material):
        where = "Fiber"
        _check_material(material, where)
        self.material = material
        self.y = _as_number(y, "y", where)
        self.z = _as_number(z, "z", where)
        self.area = _as_positive(area, "area", where)

    def discretize(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """``(y, z, A)`` de la fibra."""
        return np.array([self.y]), np.array([self.z]), np.array([self.area])


# Forma de cada parche en el YAML (clave ``shape``).
_PATCH_SHAPES = {
    "rectangle": RectangularPatch,
    "circle": CircularPatch,
    "layer": BarLayer,
}


def _accepted_kwargs(cls) -> set:
    """Parámetros del constructor de una primitiva, salvo ``material``."""
    params = inspect.signature(cls.__init__).parameters
    return {n for n in params if n not in ("self", "material")}


# ---------------------------------------------------------------------------
# Sección por fibras
# ---------------------------------------------------------------------------

@SectionRegistry.register
class FiberSection(Section):
    """Sección por fibras (ADR 0021, spec ``FiberSection``).

    Parameters
    ----------
    patches : list, optional
        Primitivas de discretización (:class:`RectangularPatch`,
        :class:`CircularPatch`, :class:`BarLayer`).
    fibers : list, optional
        Fibras sueltas (:class:`Fiber`).
    GJ : float, optional
        Rigidez torsional elástica de la sección (> 0). Obligatoria sólo si
        la sección la usa un ``Frame3D``.

    Atributos
    ---------
    y, z, area : np.ndarray (n_fibers,)
        Coordenadas locales (desde el eje del elemento) y área de cada fibra,
        en el orden de declaración: primero los parches y después las fibras
        sueltas.
    n_fibers : int
    GJ : float | None

    Notes
    -----
    Sin estado (ADR 0021 §1): el estado de cada fibra en cada punto de
    Gauss lo guarda el elemento. ``initial_states()`` devuelve una lista de
    ``None`` (estado virgen de cada material) y ``response_*`` devuelve la
    lista de estados trial sin modificar la recibida.
    """

    def __init__(self, patches=None, fibers=None, GJ=None):
        primitives = []
        for label, group in (("patches", patches), ("fibers", fibers)):
            if group is None:
                continue
            if not isinstance(group, (list, tuple)):
                raise ValueError(
                    f"FiberSection: '{label}' debe ser una lista; se recibió "
                    f"{type(group).__name__}."
                )
            for i, prim in enumerate(group):
                if not callable(getattr(prim, "discretize", None)):
                    raise ValueError(
                        f"FiberSection: {label}[{i}] no es una primitiva de "
                        f"fibras (RectangularPatch, CircularPatch, BarLayer o "
                        f"Fiber); se recibió {type(prim).__name__}."
                    )
                primitives.append(prim)

        ys, zs, areas, fiber_materials = [], [], [], []
        for prim in primitives:
            y, z, a = prim.discretize()
            ys.append(np.asarray(y, dtype=float))
            zs.append(np.asarray(z, dtype=float))
            areas.append(np.asarray(a, dtype=float))
            fiber_materials.extend([prim.material] * len(y))
        if not fiber_materials:
            raise ValueError(
                "FiberSection: la sección no tiene fibras; declare al menos un "
                "parche ('patches') o una fibra ('fibers')."
            )

        self.y = np.ascontiguousarray(np.concatenate(ys))
        self.z = np.ascontiguousarray(np.concatenate(zs))
        self.area = np.ascontiguousarray(np.concatenate(areas))
        if not np.all(self.area > 0.0):
            raise ValueError("FiberSection: todas las fibras deben tener área > 0.")
        self._fiber_materials = tuple(fiber_materials)

        # Grupos por material (identidad), en orden de aparición.
        position: dict[int, int] = {}
        groups: list[tuple] = []
        for i, mat in enumerate(fiber_materials):
            key = id(mat)
            if key not in position:
                position[key] = len(groups)
                groups.append((mat, []))
            groups[position[key]][1].append(i)
        self._groups = [(mat, list(idx)) for mat, idx in groups]
        self._materials = [mat for mat, _ in groups]

        # Matrices W de la spec §4: ε = W·e, s = Wᵀ(σ∘A), k_s = Wᵀ diag(E_t∘A) W.
        ones = np.ones(self.n_fibers)
        self._W2 = np.ascontiguousarray(np.column_stack([ones, -self.y]))
        self._W3 = np.ascontiguousarray(np.column_stack([ones, self.z, -self.y]))

        if GJ is not None:
            GJ = _as_positive(GJ, "GJ", "FiberSection")
        self.GJ = GJ

    # ------------------------------------------------------------------
    # Construcción desde el YAML (ADR 0021 §2)
    # ------------------------------------------------------------------

    @classmethod
    def from_yaml(cls, data: dict, resolve) -> "FiberSection":
        """Construye la sección desde su entrada del bloque ``sections``.

        ``data`` son los campos de la entrada salvo ``id`` y ``type``;
        ``resolve(registro, id, contexto)`` devuelve el material de id dado
        (o lanza el error del lector, que nombra el contexto)::

            patches:
              - {shape: rectangle, material: 1, y: [-0.25, 0.25], z: [-0.15, 0.15], ny: 50, nz: 1}
              - {shape: circle, material: 1, r_outer: 0.2, r_inner: 0.0, center: [0.0, 0.0], n_rad: 20, n_circ: 32}
              - {shape: layer, material: 2, start: [-0.2, -0.1], end: [-0.2, 0.1], n: 3, area: 3.1e-4}
            fibers:
              - {y: 0.1, z: 0.0, area: 1.0e-4, material: 2}
            GJ: 1.0e6
        """
        data = dict(data)
        patches_data = data.pop("patches", None) or []
        fibers_data = data.pop("fibers", None) or []
        GJ = data.pop("GJ", None)
        if data:
            raise ValueError(
                f"FiberSection: parámetros no aceptados {sorted(data)}. "
                f"Admitidos: ['GJ', 'fibers', 'patches']."
            )
        if not isinstance(patches_data, list) or not isinstance(fibers_data, list):
            raise ValueError("FiberSection: 'patches' y 'fibers' deben ser listas.")

        patches = []
        for i, entry in enumerate(patches_data):
            ctx = f"patches[{i}]"
            if not isinstance(entry, dict):
                raise ValueError(f"{ctx}: cada parche debe ser un diccionario.")
            entry = dict(entry)
            shape = entry.pop("shape", None)
            patch_cls = _PATCH_SHAPES.get(shape)
            if patch_cls is None:
                raise ValueError(
                    f"{ctx}: shape={shape!r} desconocida. Admitidas: "
                    f"{sorted(_PATCH_SHAPES)}."
                )
            patches.append(cls._build_primitive(patch_cls, entry, resolve, ctx))

        fibers = []
        for i, entry in enumerate(fibers_data):
            ctx = f"fibers[{i}]"
            if not isinstance(entry, dict):
                raise ValueError(f"{ctx}: cada fibra debe ser un diccionario.")
            fibers.append(cls._build_primitive(Fiber, dict(entry), resolve, ctx))

        return cls(patches=patches, fibers=fibers, GJ=GJ)

    @staticmethod
    def _build_primitive(prim_cls, entry: dict, resolve, ctx: str):
        if "material" not in entry:
            raise ValueError(f"{ctx}: falta el campo obligatorio 'material'.")
        material = resolve(MaterialRegistry, entry.pop("material"), f"{ctx}.material")
        accepted = _accepted_kwargs(prim_cls)
        unknown = sorted(set(entry) - accepted)
        if unknown:
            raise ValueError(
                f"{ctx}: parámetros {unknown} no aceptados por "
                f"{prim_cls.__name__}. Admitidos: {sorted(accepted | {'material'})}."
            )
        try:
            return prim_cls(material=material, **entry)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{ctx}: {exc}") from exc

    # ------------------------------------------------------------------
    # Contrato de Section
    # ------------------------------------------------------------------

    @property
    def n_fibers(self) -> int:
        return int(self.area.shape[0])

    def materials(self) -> list:
        """Materiales distintos de la sección, en orden de aparición."""
        return list(self._materials)

    def initial_states(self) -> list:
        """Estado virgen de cada fibra (``None``: el material parte de su
        estado natural en la primera evaluación)."""
        return [None] * self.n_fibers

    def total_area(self) -> float:
        return float(self.area.sum())

    def strains_2d(self, e) -> np.ndarray:
        """``ε_i = ε₀ − y_i·κ_z`` para ``e = [ε₀, κ_z]``."""
        return self._W2 @ self._check_e(e, 2)

    def strains_3d(self, e) -> np.ndarray:
        """``ε_i = ε₀ − y_i·κ_z + z_i·κ_y`` para ``e = [ε₀, κ_y, κ_z]``."""
        return self._W3 @ self._check_e(e, 3)

    def response_2d(self, e, states):
        """``s = [N, M_z]``, ``k_s`` (2×2) y estados trial para ``e = [ε₀, κ_z]``."""
        return self._response(self._W2, self._check_e(e, 2), states)

    def response_3d(self, e, states):
        """``s = [N, M_y, M_z]``, ``k_s`` (3×3) y estados trial para
        ``e = [ε₀, κ_y, κ_z]``."""
        return self._response(self._W3, self._check_e(e, 3), states)

    def mass_properties(self) -> dict:
        """``ρA = Σρ_i A_i``, ``ρI_y = Σρ_i z_i² A_i``, ``ρI_z = Σρ_i y_i² A_i``
        (spec §6), respecto del eje del elemento."""
        missing = sorted({type(m).__name__ for m in self._materials
                          if getattr(m, "density", None) is None})
        if missing:
            raise ValueError(
                f"FiberSection.mass_properties: la masa de la sección no puede "
                f"calcularse; materiales de fibra sin `density` declarada (ADR "
                f"0008): {missing}. Declare `density` en cada uno (0.0 si el "
                f"material no tiene masa por diseño)."
            )
        rho_a = np.array([float(m.density) for m in self._fiber_materials]) * self.area
        return {
            "rhoA": float(rho_a.sum()),
            "rhoIy": float((rho_a * self.z * self.z).sum()),
            "rhoIz": float((rho_a * self.y * self.y).sum()),
        }

    # ------------------------------------------------------------------
    # Evaluación
    # ------------------------------------------------------------------

    @staticmethod
    def _check_e(e, dim: int) -> np.ndarray:
        e = np.asarray(e, dtype=float).reshape(-1)
        if e.shape[0] != dim:
            raise ValueError(
                f"FiberSection: se esperaban {dim} deformaciones generalizadas; "
                f"se recibieron {e.shape[0]}."
            )
        return e

    def _response(self, W: np.ndarray, e: np.ndarray, states):
        n = self.n_fibers
        if states is None:
            states = [None] * n
        elif len(states) != n:
            raise ValueError(
                f"FiberSection: {len(states)} estados de fibra para "
                f"{n} fibras."
            )
        eps = (W @ e).tolist()
        sigma = np.empty(n)
        tangent = np.empty(n)
        new_states = [None] * n
        for material, idx in self._groups:
            compute = material.compute_state
            for i in idx:
                s_i, et_i, st_i = compute(eps[i], states[i])
                sigma[i] = s_i
                tangent[i] = et_i
                new_states[i] = st_i
        s = W.T @ (sigma * self.area)
        k = W.T @ (W * (tangent * self.area)[:, None])
        # Simetría exacta (el producto matricial puede sumar en distinto orden
        # las entradas (i, j) y (j, i)).
        k = 0.5 * (k + k.T)
        return s, k, new_states

    def __repr__(self) -> str:
        return (f"FiberSection(n_fibers={self.n_fibers}, "
                f"materials={[type(m).__name__ for m in self._materials]}, GJ={self.GJ})")
