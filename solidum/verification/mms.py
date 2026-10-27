"""Verificación por soluciones manufacturadas (MMS) (ADR 0022).

Se inventa una solución suave y SymPy deriva, de forma exacta, las cargas que
la hacen exacta:

- elasticidad lineal (:class:`ManufacturedElasticity`), con un desplazamiento
  ``u*(x)`` en 2D o en 3D::

      ε* = ∇ˢu*  (Voigt del proyecto, γ_ij = 2·ε_ij),   σ* = C·ε*,
      b* = −∇·σ*,   t*(x, n) = σ*·n;

- conducción de calor (:class:`ManufacturedConduction`), con una temperatura
  ``T*(x)`` o ``T*(x, t)`` y conductividad ``k`` anisótropa o no::

      q* = −k·∇T*,   Q* = ρc·∂T*/∂t + ∇·q*,   q̄*(x, n) = q*·n  (> 0 saliente);

- barras y marcos rectos (:class:`ManufacturedFrame`), con los campos
  locales de la línea (axial, flexión de Navier-Bernoulli o de Timoshenko,
  torsión lineal): la carga de línea como fuerza de cuerpo y las fuerzas
  internas como cargas nodales en los extremos. En 1D los valores nodales son
  exactos (Tong 1969) y los órdenes van por componente.

Con la carga de volumen, la de contorno en la parte de Neumann y el campo
impuesto en los nodos de la parte de Dirichlet, la solución inventada es la
exacta. La carga de contorno se evalúa con la normal del borde o de la cara
discretos, así que la solución inventada resuelve exactamente el problema en
el dominio de la propia malla, también con contornos curvos: el error medido
es sólo de discretización.

Al refinar una malla asintóticamente afín, con ``p`` el grado del polinomio
completo que reproducen las funciones de forma (:func:`interpolation_degree`),

    ‖e‖_L2 ~ h^(p+1),    ‖e‖_E ~ h^p,

con ``‖e‖_E² = ∫ (ε*−ε_h)ᵀ C (ε*−ε_h) dΩ`` en elasticidad y
``∫ ∇e·k·∇e dΩ`` en conducción. En el transitorio, el error en el instante
final suma el de la integración en el tiempo con el método θ, de orden 2 con
Crank-Nicolson (θ = ½) y 1 con cualquier otro θ: :func:`convergence_study`
refina la malla con ``Δt`` acoplado a ``h`` y :func:`time_convergence_study`
refina ``Δt`` sobre una malla fija con un ``T*`` lineal en ``x``, que la malla
reproduce exactamente.

Requiere SymPy (``pip install solidum-fem[verification]``). Sirve para
cualquier elemento que declare ``REFERENCE_CELL``, ``REFERENCE_NODES``,
``BATCH_KINEMATICS`` y funciones de forma por lotes, incluidos los de un
módulo de usuario: sólidos de desplazamientos ``(ux, uy[, uz])`` con un
material elástico lineal y elementos térmicos ``(T)`` con un material de
conducción lineal.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from solidum.math.assembly import Assembler
from solidum.math.solvers import LinearSolver
from solidum.math.solvers.theta_method import ThetaMethodSolver
from solidum.mesh import structured_mesh

_DOFS = ("ux", "uy", "uz")
_CELL_DIM = {"quad": 2, "tri": 2, "hex": 3, "tet": 3}


def _sympy():
    try:
        import sympy
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "solidum.verification necesita SymPy: pip install solidum-fem[verification]"
        ) from exc
    return sympy


def _parse(expressions, symbols, names: str) -> list:
    """Expresiones de SymPy a partir de cadenas o expresiones, con los nombres
    ``names`` ligados a ``symbols``; rechaza cualquier otro símbolo."""
    sp = _sympy()
    loc = dict(zip(names, symbols))
    out = [sp.sympify(e, locals=loc) for e in expressions]
    extra = set().union(*(e.free_symbols for e in out)) - set(symbols)
    if extra:
        raise ValueError(f"Símbolos desconocidos en la solución manufacturada: "
                         f"{', '.join(sorted(map(str, extra)))}; se admiten "
                         f"{', '.join(names)}.")
    return out


def linear_elastic_matrix(material, rtol: float = 1e-10) -> np.ndarray:
    """Matriz constitutiva ``C`` (Voigt 2D 3×3 o 3D 6×6) de un material
    elástico lineal, leída de ``compute_state``. Lanza ``ValueError`` si
    ``σ(ε)`` no es ``C·ε`` (material no lineal: la MMS de este módulo es
    lineal)."""
    m = int(getattr(material, "STRAIN_DIM", 3))
    _, C, _ = material.compute_state(np.zeros(m), None)
    C = np.array(C, dtype=np.float64)
    eps = 1e-4 * np.linspace(0.7, -1.1, m)
    sigma, _, _ = material.compute_state(eps, None)
    if np.linalg.norm(np.asarray(sigma) - C @ eps) > rtol * np.linalg.norm(C @ eps):
        raise ValueError(f"{type(material).__name__} no es elástico lineal: "
                         "la verificación por soluciones manufacturadas de este "
                         "módulo exige σ = C·ε.")
    return C


class _ManufacturedSolution:
    """Base de las soluciones manufacturadas: evalúa las expresiones de SymPy,
    compiladas, en un punto o en muchos.

    Cada subclase fija ``dim``, ``dofs`` (los DOFs del campo primario),
    ``energy_matrix`` (``C`` o ``k``), ``expressions``, los nombres
    ``_PRIMARY`` y ``_GRADIENT`` de las expresiones del campo y de la medida
    de su gradiente que entra en la norma de energía (deformación Voigt o
    ``∇T``), :meth:`_check_material` y :meth:`_loads`.
    """
    time_dependent = False
    _PRIMARY = _GRADIENT = ""

    @property
    def mesh_dim(self) -> int:
        """Dimensión de la malla (de la celda de referencia): la del
        espacio, salvo en los elementos de línea."""
        return self.dim

    def _compile(self, symbols) -> None:
        sp = _sympy()
        # Dos versiones: con math para un punto (las cargas se evalúan punto a
        # punto, y NumPy sobre escalares es ~10 veces más lento) y con NumPy
        # para muchos a la vez (normas del error).
        self._f = {k: sp.lambdify(symbols, list(v), "numpy", cse=True)
                   for k, v in self.expressions.items()}
        self._f_point = {k: sp.lambdify(symbols, list(v), "math", cse=True)
                         for k, v in self.expressions.items()}

    def _args(self, x: np.ndarray) -> tuple:
        """Argumentos espaciales de las expresiones compiladas: las
        coordenadas ``x``, ``y`` (``z``) de un punto ``(d,)`` o de varios
        ``(n, d)``."""
        return tuple(x[..., k] for k in range(self.dim))

    def _eval(self, name, x, t: float = 0.0) -> np.ndarray:
        """Valores en un punto ``(d,)`` → ``(m,)`` o en varios ``(n, d)`` →
        ``(n, m)``, en el instante ``t`` si la solución depende del tiempo."""
        x = np.asarray(x, dtype=np.float64)
        args = self._args(x) + ((float(t),) if self.time_dependent else ())
        if x.ndim == 1:
            return np.array(self._f_point[name](*(float(a) for a in args)), dtype=np.float64)
        vals = self._f[name](*args)
        return np.column_stack([np.broadcast_to(np.asarray(v, dtype=np.float64), (len(x),))
                                for v in vals])

    def _dof_values(self, x) -> np.ndarray:
        """Valores de los DOFs (``dofs``) del campo exacto en el nodo ``x``,
        los que se imponen en la parte de Dirichlet."""
        return self._eval(self._PRIMARY, x)

    def _constrain(self, mesh, fixed: dict) -> None:
        """Restricciones adicionales antes de numerar las ecuaciones; ninguna
        por omisión. ``fixed``: los nodos de Dirichlet por id."""

    def _check_mesh(self, mesh) -> None:
        """Rechaza una malla cuyos elementos no tengan la constitutiva de la
        solución (con otra, la solución manufacturada no sería la exacta)."""
        self._check_material(next(iter(mesh.domain.elements.values())).material)

    def _check_material(self, material) -> None:
        raise NotImplementedError

    def _loads(self, mesh, assembler, dirichlet, t: float = 0.0) -> np.ndarray:
        """Vector global de cargas en el instante ``t``: la de volumen en todos
        los elementos y la de contorno en los lados que no son ``dirichlet``."""
        raise NotImplementedError


class ManufacturedElasticity(_ManufacturedSolution):
    """Solución manufacturada de la elasticidad lineal con ``C`` constante.

    Parameters
    ----------
    displacement : 2 ó 3 expresiones en ``x``, ``y`` (``z``)
        Las componentes de ``u*`` como cadenas (``"1e-3*sin(x + 2*y)"``) o
        expresiones de SymPy; su número fija la dimensión.
    C : (3, 3) en 2D o (6, 6) en 3D
        Matriz constitutiva en el Voigt del proyecto; :meth:`from_material`
        la toma del material.
    """
    _PRIMARY, _GRADIENT = "displacement", "strain"

    def __init__(self, displacement, C):
        sp = _sympy()
        self.dim = d = len(displacement)
        if d not in (2, 3):
            raise ValueError("displacement: se esperan 2 ó 3 componentes.")
        self.dofs = _DOFS[:d]
        xs = sp.symbols("x y z", real=True)[:d]
        u = _parse(displacement, xs, "xyz"[:d])
        m = 3 if d == 2 else 6
        self.C = np.array(C, dtype=np.float64).reshape(m, m)
        g = [[sp.diff(u[i], xs[j]) for j in range(d)] for i in range(d)]
        # Voigt del proyecto: normales y después γ_xy[, γ_yz, γ_xz].
        pairs = ((0, 1),) if d == 2 else ((0, 1), (1, 2), (0, 2))
        eps = sp.Matrix([g[i][i] for i in range(d)] + [g[i][j] + g[j][i] for i, j in pairs])
        sig = sp.Matrix(self.C.tolist()) * eps
        # Índice Voigt de cada componente del tensor de esfuerzos.
        self._voigt = np.empty((d, d), dtype=int)
        for i in range(d):
            self._voigt[i, i] = i
        for k, (i, j) in enumerate(pairs):
            self._voigt[i, j] = self._voigt[j, i] = d + k
        b = sp.Matrix([-sum(sp.diff(sig[self._voigt[i, j]], xs[j]) for j in range(d))
                       for i in range(d)])
        self.expressions = {"displacement": sp.Matrix(u), "strain": eps,
                            "stress": sig, "body_force": b}
        self._compile(xs)

    @classmethod
    def from_material(cls, displacement, material) -> "ManufacturedElasticity":
        return cls(displacement, linear_elastic_matrix(material))

    @property
    def energy_matrix(self) -> np.ndarray:
        return self.C

    def displacement(self, x) -> np.ndarray:
        return self._eval("displacement", x)

    def strain(self, x) -> np.ndarray:
        return self._eval("strain", x)

    def stress(self, x) -> np.ndarray:
        return self._eval("stress", x)

    def body_force(self, x) -> np.ndarray:
        """``b*(x)``: el campo que se pasa a ``compute_body_load``."""
        return self._eval("body_force", x)

    def traction(self, x, n) -> np.ndarray:
        """``t*(x, n) = σ*(x)·n``: el campo que se pasa a
        ``compute_edge_traction`` (2D) o ``compute_face_traction`` (3D)."""
        return self.stress(x)[self._voigt] @ np.asarray(n, dtype=np.float64)[:self.dim]

    def _check_material(self, material) -> None:
        if not np.allclose(linear_elastic_matrix(material), self.C, rtol=1e-12, atol=0.0):
            raise ValueError("La C del material no es la de la solución manufacturada.")

    def _loads(self, mesh, assembler, dirichlet, t: float = 0.0) -> np.ndarray:
        F = assembler.assemble_body_load(self.body_force)
        for side, entities in mesh.boundary.items():
            if side in dirichlet:
                continue
            for el, k in entities:
                load = el.compute_edge_traction if mesh.dim == 2 else el.compute_face_traction
                F[el.get_global_dof_indices()] += load(k, self.traction)
        return F


class ManufacturedConduction(_ManufacturedSolution):
    """Solución manufacturada de la conducción de calor con ``k`` y ``ρc``
    constantes.

    Parameters
    ----------
    temperature : expresión en ``x``, ``y`` (``z``) y, si es transitoria, ``t``
        ``T*`` como cadena (``"sin(x + 2*y)*exp(-t)"``) o expresión de SymPy.
        Si depende de ``t`` el problema es transitorio.
    k : (2, 2) o (3, 3)
        Tensor de conductividad, anisótropo o no; su tamaño fija la
        dimensión. :meth:`from_material` lo toma del material.
    rho_c : float, opcional
        Capacidad volumétrica ``ρ·c``; obligatoria si ``T*`` depende de ``t``.
    """
    dofs = ("T",)
    _PRIMARY, _GRADIENT = "temperature", "gradient"

    def __init__(self, temperature, k, rho_c: float | None = None):
        sp = _sympy()
        self.k = np.array(k, dtype=np.float64)
        if self.k.ndim != 2 or self.k.shape[0] != self.k.shape[1] or self.k.shape[0] not in (2, 3):
            raise ValueError("k: se espera el tensor de conductividad, 2×2 ó 3×3.")
        self.dim = d = self.k.shape[0]
        xs = sp.symbols("x y z", real=True)[:d]
        t = sp.Symbol("t", real=True)
        (T,) = _parse([temperature], (*xs, t), "xyz"[:d] + "t")
        self.time_dependent = t in T.free_symbols
        if self.time_dependent and rho_c is None:
            raise ValueError("T* depende de t: hace falta rho_c (ρ·c) para la fuente ρc·∂T*/∂t.")
        self.rho_c = None if rho_c is None else float(rho_c)
        grad = sp.Matrix([sp.diff(T, s) for s in xs])
        q = -sp.Matrix(self.k.tolist()) * grad
        Q = sum(sp.diff(q[i], xs[i]) for i in range(d))
        if self.time_dependent:
            Q += self.rho_c * sp.diff(T, t)
        # T* lineal en x (cualquier dependencia de t): lo reproduce todo
        # elemento isoparamétrico, en cualquier malla.
        self.affine_in_space = all(sp.simplify(sp.diff(T, a, b)) == 0 for a in xs for b in xs)
        self.expressions = {"temperature": sp.Matrix([T]), "gradient": grad,
                            "flux": q, "source": sp.Matrix([Q])}
        self._compile((*xs, t) if self.time_dependent else xs)

    @classmethod
    def from_material(cls, temperature, material) -> "ManufacturedConduction":
        """Con la conductividad del material y, si ``T*`` depende de ``t``,
        su capacidad volumétrica ``ρ·c``."""
        rho_c = None
        if material.density is not None and material.specific_heat is not None:
            rho_c = material.volumetric_capacity()
        return cls(temperature, material.conductivity, rho_c)

    @property
    def energy_matrix(self) -> np.ndarray:
        return self.k

    def _scalar(self, name, x, t):
        v = self._eval(name, x, t)
        return float(v[0]) if v.ndim == 1 else v[:, 0]

    def temperature(self, x, t: float = 0.0):
        """``T*`` en un punto (escalar) o en varios ``(n,)``."""
        return self._scalar("temperature", x, t)

    def gradient(self, x, t: float = 0.0) -> np.ndarray:
        return self._eval("gradient", x, t)

    def flux(self, x, t: float = 0.0) -> np.ndarray:
        """Flujo de Fourier ``q* = −k·∇T*``."""
        return self._eval("flux", x, t)

    def source(self, x, t: float = 0.0):
        """``Q*(x, t)``: el campo que se pasa a ``compute_body_source``."""
        return self._scalar("source", x, t)

    def boundary_flux(self, x, n, t: float = 0.0) -> float:
        """``q̄*(x, n) = q*(x)·n``, positivo saliente: el campo que se pasa a
        ``compute_edge_flux`` (2D) o ``compute_face_flux`` (3D)."""
        return float(self.flux(x, t) @ np.asarray(n, dtype=np.float64)[:self.dim])

    def _check_material(self, material) -> None:
        if not np.allclose(material.conductivity, self.k, rtol=1e-12, atol=0.0):
            raise ValueError("La conductividad del material no es la de la solución "
                             "manufacturada.")
        if self.time_dependent and not np.isclose(
                material.volumetric_capacity(), self.rho_c, rtol=1e-12, atol=0.0):
            raise ValueError("La capacidad ρ·c del material no es la de la solución "
                             "manufacturada.")

    def _loads(self, mesh, assembler, dirichlet, t: float = 0.0) -> np.ndarray:
        F = np.zeros(mesh.domain.total_dofs)
        for el in mesh.domain.elements.values():
            F[el.get_global_dof_indices()] += el.compute_body_source(lambda x: self.source(x, t))
        for side, entities in mesh.boundary.items():
            if side in dirichlet:
                continue
            for el, k in entities:
                load = el.compute_edge_flux if mesh.dim == 2 else el.compute_face_flux
                F[el.get_global_dof_indices()] += load(
                    k, lambda x, n: self.boundary_flux(x, n, t))
        return F


class ManufacturedFrame(_ManufacturedSolution):
    """Solución manufacturada de barras y marcos rectos, elásticos y lineales
    (fase 4).

    Los campos se dan en los ejes locales de la línea, como expresiones de la
    abscisa ``s = (x − origin)·e_x`` a lo largo del eje:

    - ``"u"``: desplazamiento axial; ``ε = u'``, ``N = EA·ε``;
    - ``"v"``: flexión en el plano ``xy`` (Navier-Bernoulli): ``θz = v'``,
      ``κz = v''``, ``Mz = EI_z·κz``;
    - ``"w"``: flexión en el plano ``xz`` (3D): ``θy = −w'`` por la regla de la
      mano derecha, ``κy = θy'``, ``My = EI_y·κy``;
    - ``"psi"``: en Timoshenko, la función de la que salen ``θz = ψ'`` y
      ``v = ψ − (EI/(G·A_s))·ψ''``, la pareja que no necesita momento
      distribuido: ``Vy = G·A_s·(v' − θz) = −EI·ψ'''``; con cortante rígido
      ``v = ψ``;
    - ``"theta_x"``: giro de torsión (3D), lineal, porque no hay torsor
      distribuido: ``T = GJ·θx'``.

    De la forma débil de cada elemento salen las cargas que hacen exacta la
    solución: las de línea ``p = −N'``, ``q_y = −Vy'`` y ``q_z = −Vz'``, con
    ``Vy = −Mz'`` y ``Vz = My'`` en Navier-Bernoulli (``Mz' = −Vy`` y
    ``My' = Vz``, convención stress-resultant, Reglas §5), que se aplican como
    la fuerza de cuerpo ``b = λᵀ·q/A``; y en un extremo libre, las fuerzas
    internas de la sección por la normal exterior ``±e_x`` como cargas
    nodales.

    Parameters
    ----------
    fields : dict
        Expresiones en ``s`` (cadena o SymPy) de los campos que el elemento
        admite: ``{"u"}`` en una armadura, ``{"u", "v"}`` en un marco 2D,
        ``{"u", "psi"}`` en el Timoshenko, ``{"u", "v", "w", "theta_x"}`` en
        el marco 3D; los que falten valen cero.
    rigidities : dict
        Rigidez de la sección por deformación generalizada, con los nombres
        de ``LINE_STRAINS`` del elemento (``"eps"``: EA, ``"kappa_z"``: EI_z,
        ``"kappa_y"``: EI_y, ``"gamma_y"``: G·A_s, ``"twist"``: GJ).
    area : float
        Área con la que el elemento convierte la fuerza de cuerpo en carga por
        unidad de longitud.
    origin : (d,)
        Punto de ``s = 0``.
    axes : (d, d)
        Ejes locales en globales, por filas (``e_x``, ``e_y`` [, ``e_z``]).
    dofs : sequence
        ``DOF_NAMES`` del elemento.

    :meth:`from_element` toma todo salvo los campos y el origen de un
    elemento de la malla.
    """

    _PRIMARY, _GRADIENT = "state", ""
    _KINDS = {"bar": {"u"}, "frame2d": {"u", "v"}, "timoshenko": {"u", "psi"},
              "frame3d": {"u", "v", "w", "theta_x"}}

    def __init__(self, fields: dict, rigidities: dict, area: float, origin, axes, dofs):
        sp = _sympy()
        self.axes = np.array(axes, dtype=np.float64)
        self.dim = d = self.axes.shape[0]
        self.origin = np.asarray(origin, dtype=np.float64).reshape(d)
        self.dofs = tuple(dofs)
        self.rigidities = {k: float(v) for k, v in rigidities.items()}
        self.area = float(area)
        kind = ("timoshenko" if "gamma_y" in rigidities else "frame3d" if "kappa_y" in rigidities
                else "frame2d" if "kappa_z" in rigidities else "bar")
        self.kind = kind
        unknown = set(fields) - self._KINDS[kind]
        if unknown:
            raise ValueError(f"Campos {sorted(unknown)} no admitidos en un elemento de tipo "
                             f"{kind}; se admiten {sorted(self._KINDS[kind])}.")
        s = sp.Symbol("s", real=True)
        f = {k: _parse([v], (s,), "s")[0] for k, v in fields.items()}
        zero = sp.Integer(0)
        R = self.rigidities
        u = f.get("u", zero)
        eps = sp.diff(u, s)
        N = R["eps"] * eps
        v = w = tx = ty = tz = kz = ky = gy = twist = zero
        Vy = Vz = T = My = Mz = zero
        if kind == "timoshenko":
            psi = f.get("psi", zero)
            tz = sp.diff(psi, s)
            v = psi - R["kappa_z"] / R["gamma_y"] * sp.diff(psi, s, 2)
            kz = sp.diff(tz, s)
            gy = sp.diff(v, s) - tz
            Mz, Vy = R["kappa_z"] * kz, R["gamma_y"] * gy
        elif kind in ("frame2d", "frame3d"):
            v = f.get("v", zero)
            tz = sp.diff(v, s)
            kz = sp.diff(tz, s)
            Mz = R["kappa_z"] * kz
            Vy = -sp.diff(Mz, s)
        if kind == "frame3d":
            w = f.get("w", zero)
            ty = -sp.diff(w, s)
            ky = sp.diff(ty, s)
            My = R["kappa_y"] * ky
            Vz = sp.diff(My, s)
            tx = f.get("theta_x", zero)
            if sp.simplify(sp.diff(tx, s, 2)) != 0:
                raise ValueError("theta_x debe ser lineal en s: no hay torsor distribuido.")
            twist = sp.diff(tx, s)
            T = R["twist"] * twist
        strains = {"eps": eps, "twist": twist, "kappa_y": ky, "kappa_z": kz, "gamma_y": gy}
        self.expressions = {
            "state": sp.Matrix([u, v, w, tx, ty, tz]),
            "strains": sp.Matrix([strains[k] for k in rigidities]),
            "end_forces": sp.Matrix([N, Vy, Vz, T, My, Mz]),
            "line_load": sp.Matrix([-sp.diff(N, s), -sp.diff(Vy, s), -sp.diff(Vz, s)]),
        }
        self._strain_names = tuple(rigidities)
        self._compile((s,))
        # Ejes en 3D (en 2D, el plano xy con z = e_z) para pasar de locales a
        # globales desplazamientos, giros, fuerzas y momentos.
        self._axes3 = np.eye(3)
        self._axes3[:d, :d] = self.axes

    @classmethod
    def from_element(cls, fields: dict, element, origin) -> "ManufacturedFrame":
        """Con las rigideces, el área, los ejes y los DOFs del elemento."""
        rig = dict(zip(element.LINE_STRAINS, element.section_rigidities()))
        area = element._axial_area() if hasattr(element, "_axial_area") else element.A
        return cls(fields, rig, area, origin, _element_axes(element), element.DOF_NAMES)

    @property
    def mesh_dim(self) -> int:
        return 1

    @property
    def energy_matrix(self) -> np.ndarray:
        return np.diag([self.rigidities[k] for k in self._strain_names])

    def _args(self, x: np.ndarray) -> tuple:
        return ((x[..., :self.dim] - self.origin) @ self.axes[0],)

    def local(self, names, x) -> np.ndarray:
        """Campos (``"u"``, ``"v"``, ``"w"``, ``"theta_x"``, ``"theta_y"``,
        ``"theta_z"``) o deformaciones (``LINE_STRAINS``) locales en los puntos
        ``x`` ``(n, d)``, ``(n, len(names))``."""
        state = ("u", "v", "w", "theta_x", "theta_y", "theta_z")
        out = []
        for name in names:
            if name in state:
                out.append(self._eval("state", x)[:, state.index(name)])
            else:
                out.append(self._eval("strains", x)[:, self._strain_names.index(name)])
        return np.column_stack(out)

    def _to_dofs(self, local6: np.ndarray) -> np.ndarray:
        """``[u, v, w, θx, θy, θz]`` (o fuerzas y momentos) locales → valores
        en los ``dofs`` globales del elemento."""
        g = np.concatenate([self._axes3.T @ local6[:3], self._axes3.T @ local6[3:]])
        index = {"ux": 0, "uy": 1, "uz": 2, "rx": 3, "ry": 4, "rz": 5}
        return np.array([g[index[n]] for n in self.dofs])

    def _dof_values(self, x) -> np.ndarray:
        return self._to_dofs(self._eval("state", x))

    def body_force(self, x) -> np.ndarray:
        """``b(x) = λᵀ·q(s)/A``: el campo que se pasa a ``compute_body_load``."""
        return (self._axes3.T @ self._eval("line_load", x))[:self.dim] / self.area

    def end_loads(self, x, outward: float) -> np.ndarray:
        """Cargas nodales en un extremo de normal exterior ``outward·e_x``:
        las fuerzas internas de la sección, en los ``dofs`` globales."""
        return outward * self._to_dofs(self._eval("end_forces", x))

    def _constrain(self, mesh, fixed: dict) -> None:
        """Una barra no tiene rigidez transversal: cada nodo libre se mueve sólo
        a lo largo del eje, ``u_j = (e_j/e_k)·u_k`` con ``k`` la componente
        mayor de ``e_x``."""
        if self.kind != "bar":
            return
        e = self.axes[0]
        k = int(np.argmax(np.abs(e)))
        for node in mesh.domain.nodes.values():
            if node.id in fixed:
                continue
            for j in range(self.dim):
                if j != k:
                    mesh.domain.add_linear_constraint((node.id, self.dofs[j]),
                                                      [(node.id, self.dofs[k])], [e[j] / e[k]])

    def _check_mesh(self, mesh) -> None:
        for el in mesh.domain.elements.values():
            rig = dict(zip(el.LINE_STRAINS, el.section_rigidities()))
            if set(rig) != set(self.rigidities) or not all(
                    np.isclose(rig[k], self.rigidities[k], rtol=1e-12, atol=0.0) for k in rig):
                raise ValueError("Las rigideces de la sección no son las de la solución "
                                 "manufacturada.")
            if not np.allclose(_element_axes(el), self.axes, rtol=0.0, atol=1e-12):
                raise ValueError("Los ejes locales del elemento no son los de la solución "
                                 "manufacturada.")

    def _loads(self, mesh, assembler, dirichlet, t: float = 0.0) -> np.ndarray:
        F = assembler.assemble_body_load(self.body_force)
        for side, ends in mesh.boundary.items():
            if side in dirichlet:
                continue
            for el, k in ends:
                node = el.nodes[k]
                loads = self.end_loads(node.coordinates, 1.0 if side == "xmax" else -1.0)
                F[[node.dofs[n] for n in self.dofs]] += loads
        return F


def _element_axes(element) -> np.ndarray:
    """Ejes locales del elemento por filas: los de su carga en los marcos; en
    una armadura, ``e_x`` y una base ortonormal cualquiera que lo complete."""
    if hasattr(element, "_load_axes"):
        return np.asarray(element._load_axes(), dtype=np.float64)
    X = element._line_ends()
    e = (X[1] - X[0]) / element.L0
    if len(e) == 2:
        return np.array([e, [-e[1], e[0]]])
    ref = np.eye(3)[int(np.argmin(np.abs(e)))]
    y = ref - (ref @ e) * e
    y /= np.linalg.norm(y)
    return np.array([e, y, np.cross(e, y)])


def _line_error_norms(domain, U: np.ndarray, solution: ManufacturedFrame,
                      n_points: int | None = None) -> dict:
    """Normas del error en elementos de línea, por componente local:
    ``"L2:<campo>"`` = ‖f* − f_h‖_L2 y ``"energy:<deformación>"`` =
    (∫ D·(e* − e_h)² dx)^½, con sus totales ``"L2"`` (desplazamientos) y
    ``"energy"``, y ``"nodal"``, el máximo error en los DOFs de los nodos."""
    xi, w = np.polynomial.legendre.leggauss(8 if n_points is None else n_points)
    xi, w = 0.5 * (xi + 1.0), 0.5 * w
    D = solution.rigidities
    sq: dict = {}
    for el in domain.elements.values():
        u_e = np.asarray(el.get_local_displacements(U), dtype=np.float64)
        f_h, e_h = el.line_field(u_e, xi)
        X = el._line_ends()
        x = X[0] + xi[:, None] * (X[1] - X[0])
        dx = w * el.L0
        for k, name in enumerate(el.LINE_FIELDS):
            df = solution.local([name], x)[:, 0] - f_h[:, k]
            sq["L2:" + name] = sq.get("L2:" + name, 0.0) + dx @ df**2
        for k, name in enumerate(el.LINE_STRAINS):
            de = solution.local([name], x)[:, 0] - e_h[:, k]
            sq["energy:" + name] = sq.get("energy:" + name, 0.0) + D[name] * (dx @ de**2)
    out = {k: float(np.sqrt(v)) for k, v in sq.items()}
    out["L2"] = float(np.sqrt(sum(v for k, v in sq.items() if k in ("L2:u", "L2:v", "L2:w"))))
    out["energy"] = float(np.sqrt(sum(v for k, v in sq.items() if k.startswith("energy:"))))
    nodal = 0.0
    for node in domain.nodes.values():
        exact = solution._dof_values(node.coordinates)
        fe = U[[node.dofs[n] for n in solution.dofs]]
        nodal = max(nodal, float(np.max(np.abs(fe - exact))))
    out["nodal"] = nodal
    return out


def _dirichlet_sides(mesh, dirichlet) -> tuple:
    return (tuple(s for s in mesh.sides if s.endswith("min")) if dirichlet is None
            else tuple(dirichlet))


def _impose_dirichlet(mesh, solution, dirichlet) -> list:
    """Fija el campo de ``solution`` (en ``t = 0``) en los nodos de los lados
    ``dirichlet``, numera las ecuaciones y devuelve los nodos fijados."""
    nodes = {node.id: node for side in dirichlet for node in mesh.nodes_on(side)}
    for node in nodes.values():
        for dof, value in zip(solution.dofs, solution._dof_values(node.coordinates)):
            node.fix_dof(dof, float(value))
    solution._constrain(mesh, nodes)
    mesh.domain.generate_equation_numbers(verbose=False)
    return list(nodes.values())


def solve_manufactured(mesh, solution, dirichlet=None, linear_algebra=None) -> np.ndarray:
    """Resuelve la malla con las cargas de ``solution``: el campo exacto en los
    nodos de los lados ``dirichlet`` (por omisión, los lados ``min``) y la
    carga de contorno exacta en los demás lados del bloque. Devuelve el vector
    global de desplazamientos o de temperaturas."""
    if solution.time_dependent:
        raise ValueError("La solución manufacturada depende del tiempo: "
                         "solve_manufactured_transient.")
    dirichlet = _dirichlet_sides(mesh, dirichlet)
    _impose_dirichlet(mesh, solution, dirichlet)
    assembler = Assembler(mesh.domain)
    F = solution._loads(mesh, assembler, dirichlet)
    kwargs = {} if linear_algebra is None else {"linear_algebra": linear_algebra}
    return LinearSolver(assembler, **kwargs).solve(F)


def solve_manufactured_transient(mesh, solution: ManufacturedConduction, t_end: float,
                                 n_steps: int, *, theta: float = 1.0, lumping: str = "lumped",
                                 dirichlet=None, linear_algebra=None):
    """Integra la malla con el método θ (:class:`ThetaMethodSolver`) de
    ``t = 0`` a ``t_end`` en ``n_steps`` pasos: ``T*(x, 0)`` como condición
    inicial en todos los nodos, ``T*(x, t)`` en los nodos de los lados
    ``dirichlet`` (por omisión, los ``min``) y ``Q*`` y ``q̄*`` en cada
    instante. Devuelve el ``ThermalTransientResult``, con el campo final en
    ``T_history[:, -1]``."""
    if not solution.time_dependent:
        raise ValueError("La solución manufacturada no depende del tiempo: solve_manufactured.")
    dirichlet = _dirichlet_sides(mesh, dirichlet)
    fixed = _impose_dirichlet(mesh, solution, dirichlet)
    dom = mesh.domain
    ndof = dom.total_dofs
    nodes = list(dom.nodes.values())
    T0 = np.zeros(ndof)
    T0[[nd.dofs["T"] for nd in nodes]] = solution.temperature(
        np.array([nd.coordinates for nd in nodes]), 0.0)
    fixed_dofs = [nd.dofs["T"] for nd in fixed]
    fixed_x = np.array([nd.coordinates for nd in fixed])

    def prescribed(t):
        g = np.zeros(ndof)
        g[fixed_dofs] = solution.temperature(fixed_x, t)
        return g

    assembler = Assembler(dom)
    kwargs = {} if linear_algebra is None else {"linear_algebra": linear_algebra}
    result = ThetaMethodSolver(
        assembler, t_end / n_steps, n_steps, T0, theta=theta, lumping=lumping,
        F_func=lambda t: solution._loads(mesh, assembler, dirichlet, t),
        dirichlet_func=prescribed, output_every=n_steps, **kwargs).solve()
    if not result.converged:
        raise RuntimeError(f"El método θ (theta={theta}) divergió con Δt = {t_end / n_steps}.")
    return result


def _reference_rule(cell: str, n: int):
    """Regla de Gauss de ``n`` puntos por dirección en la celda de referencia:
    producto tensorial en ``[-1, 1]^d``; en el triángulo y el tetraedro,
    producto colapsado (Duffy), exacto hasta grado ``2n − 2`` y ``2n − 3``."""
    a, w = np.polynomial.legendre.leggauss(n)
    if cell in ("quad", "hex"):
        d = _CELL_DIM[cell]
        grids = np.meshgrid(*([a] * d), indexing="ij")
        W = np.ones([n] * d)
        for k in range(d):
            W = W * np.expand_dims(w, [j for j in range(d) if j != k])
        return np.column_stack([g.ravel() for g in grids]), W.ravel()
    s, ws = (a + 1.0) / 2.0, w / 2.0
    if cell == "tri":
        S, T = np.meshgrid(s, s, indexing="ij")
        W = np.outer(ws, ws) * (1.0 - S)
        return np.column_stack([S.ravel(), (T * (1.0 - S)).ravel()]), W.ravel()
    if cell == "tet":
        S, T, R = np.meshgrid(s, s, s, indexing="ij")
        W = (ws[:, None, None] * ws[None, :, None] * ws[None, None, :]
             * (1.0 - S) ** 2 * (1.0 - T))
        pts = np.column_stack([S.ravel(), (T * (1.0 - S)).ravel(),
                               (R * (1.0 - S) * (1.0 - T)).ravel()])
        return pts, W.ravel()
    raise ValueError(f"Celda de referencia {cell!r} no admitida (quad, tri, hex, tet).")


def error_norms(domain, U: np.ndarray, solution, n_points: int | None = None,
                t: float = 0.0) -> dict:
    """Normas del error entre la solución manufacturada (en el instante ``t``
    si depende del tiempo) y la de elementos finitos ``U``:
    ``{"L2": ‖e‖_L2, "energy": ‖e‖_E}`` (en 2D por unidad de espesor),
    integradas en cada elemento con ``n_points`` puntos de Gauss por
    dirección, más que la regla del elemento porque el error no es
    polinómico. Por omisión ``p + 3``, con ``p`` el grado de la
    interpolación.

    En los elementos de línea (:class:`ManufacturedFrame`) las normas van por
    componente local, con Gauss de ``n_points`` (8 por omisión) puntos por
    elemento: ver :func:`_line_error_norms`."""
    if isinstance(solution, ManufacturedFrame):
        return _line_error_norms(domain, U, solution, n_points)
    d = solution.dim
    m = len(solution.dofs)
    D = solution.energy_matrix
    e_l2 = 0.0
    e_energy = 0.0
    rules = {}
    for el in domain.elements.values():
        cls = type(el)
        if cls not in rules:
            n = interpolation_degree(el) + 3 if n_points is None else n_points
            pts, wts = _reference_rule(cls.REFERENCE_CELL, n)
            rules[cls] = (pts, wts, el.batch_shape_functions(pts))
        pts, wts, N = rules[cls]
        coords = np.ascontiguousarray(el.get_coordinate_matrix(ndim=d))
        u_e = np.asarray(el.get_local_displacements(U), dtype=np.float64)
        x = N @ coords
        du = solution._eval(solution._PRIMARY, x, t) - N @ u_e.reshape(-1, m)
        de = solution._eval(solution._GRADIENT, x, t)
        B = np.empty((de.shape[1], u_e.size))
        dV = np.empty(len(wts))
        for g in range(len(wts)):
            dV[g] = cls.BATCH_KINEMATICS(pts[g], coords, B) * wts[g]
            de[g] -= B @ u_e
        e_l2 += np.sum(du * du, axis=1) @ dV
        e_energy += np.einsum("gi,ij,gj->g", de, D, de) @ dV
    return {"L2": float(np.sqrt(e_l2)), "energy": float(np.sqrt(e_energy))}


def interpolation_degree(element, max_degree: int = 6) -> int:
    """Grado ``p`` del polinomio completo que reproducen exactamente las
    funciones de forma en la celda de referencia (``Σ N_i(ξ)·m(ξ_i) = m(ξ)``
    para todo monomio ``m`` de grado ≤ p). Quad4, Tri3, Hex8, Tet4,
    Quad4Thermal y Hex8Thermal: 1; Quad8, Quad9, Tri6, Hex20, Hex27 y Tet10:
    2. Es el que fija los órdenes de convergencia teóricos."""
    cls = type(element)
    nodes = np.asarray(cls.REFERENCE_NODES, dtype=np.float64)
    d = nodes.shape[1]
    rng = np.random.default_rng(0)
    if cls.REFERENCE_CELL in ("quad", "hex"):
        pts = rng.uniform(-1.0, 1.0, (16, d))
    else:
        pts = rng.dirichlet(np.ones(d + 1), 16)[:, 1:]
    N = element.batch_shape_functions(pts)
    degree = 0
    for k in range(1, max_degree + 1):
        for powers in np.ndindex(*([k + 1] * d)):
            if sum(powers) != k:
                continue
            exact = np.prod(pts ** np.array(powers), axis=1)
            interp = N @ np.prod(nodes ** np.array(powers), axis=1)
            if not np.allclose(interp, exact, rtol=0.0, atol=1e-10):
                return degree
        degree = k
    return degree


@dataclass
class ConvergenceStudy:
    """Resultado de :func:`convergence_study` o de :func:`time_convergence_study`.

    ``h`` es el parámetro que se refina (``parameter``): el tamaño relativo de
    la malla, ``1/divisions``, o el paso de tiempo ``Δt``, con ``divisions``
    entonces el número de pasos de cada nivel. ``errors``, las normas del
    error por nivel; ``expected``, los órdenes teóricos por norma cuando no
    son los de la interpolación (``p + 1`` en L² y ``p`` en energía), como en
    el estudio en el tiempo o en los elementos de línea, cuyos órdenes van por
    componente (``degree`` es entonces ``None``)."""
    element: str
    degree: int | None
    divisions: tuple
    h: np.ndarray
    errors: dict
    expected: dict | None = None
    parameter: str = "h"

    def expected_order(self, norm: str) -> int:
        if self.expected is not None:
            return self.expected[norm]
        return self.degree + 1 if norm == "L2" else self.degree

    def orders(self, norm: str) -> np.ndarray:
        """Órdenes observados entre niveles consecutivos."""
        e = self.errors[norm]
        return np.log(e[:-1] / e[1:]) / np.log(self.h[:-1] / self.h[1:])

    def table(self) -> str:
        """Tabla de texto: nivel, parámetro y, por norma (las de ``expected``,
        o L2 y energía), errores y órdenes observados."""
        norms = list(self.expected) if self.expected is not None else ["L2", "energy"]
        title = f"{self.element} (p = {self.degree}; esperado " + ", ".join(
            f"{k} {self.expected_order(k)}" for k in norms) + ")"
        rows = [title, f"{'n':>4} {self.parameter:>9} " + " ".join(
            f"{('e_' + k)[:15]:>15} {'orden':>6}" for k in norms)]
        orders = {k: self.orders(k) for k in norms}
        for i, n in enumerate(self.divisions):
            cells = []
            for k in norms:
                o = f"{orders[k][i - 1]:6.3f}" if i else " " * 6
                cells.append(f"{self.errors[k][i]:15.4e} {o}")
            rows.append(f"{n:4d} {self.h[i]:9.5f} " + " ".join(cells))
        return "\n".join(rows)


def _name(element) -> str:
    return element if isinstance(element, str) else element.__name__


def convergence_study(element, solution, divisions, *, mapping=None, size=None, origin=None,
                      dirichlet=None, n_points: int | None = None, linear_algebra=None,
                      t_end: float | None = None, n_steps=None, theta: float = 1.0,
                      lumping: str = "lumped", expected: dict | None = None,
                      **element_kwargs) -> ConvergenceStudy:
    """Estudio de convergencia en mallas estructuradas de ``n`` celdas por
    dirección, ``n`` en ``divisions``, sobre el bloque transformado por
    ``mapping``.

    Si ``solution`` depende del tiempo, cada nivel se integra hasta ``t_end``
    con el método θ (``theta``, capacidad ``lumping``) en ``n_steps`` pasos,
    un entero o uno por nivel, y el error se mide en ``t_end``. Los órdenes
    esperados son los del espacio: ``Δt`` debe bajar lo bastante deprisa para
    que el error en el tiempo no domine (con Crank-Nicolson y ``p = 1``,
    ``Δt ∝ h`` basta).

    En los elementos de línea la malla es una línea de ``n`` elementos y los
    órdenes esperados, por componente, se pasan en ``expected``
    (``{"L2:v": 4, "energy:kappa_z": 2, ...}``).

    ``element_kwargs`` va al constructor del elemento (``material``,
    ``thickness``...). La constitutiva del elemento debe ser la de
    ``solution`` (``C``; ``k`` y, en el transitorio, ``ρ·c``; las rigideces
    de la sección y los ejes en una línea): con otra, la solución manufacturada
    no sería la exacta."""
    if solution.time_dependent:
        if t_end is None or n_steps is None:
            raise ValueError("Solución transitoria: hacen falta t_end y n_steps.")
        steps = ((int(n_steps),) * len(divisions) if np.isscalar(n_steps)
                 else tuple(int(s) for s in n_steps))
        if len(steps) != len(divisions):
            raise ValueError("n_steps: un número de pasos por nivel de divisions.")
    errors = {"L2": [], "energy": []}
    degree = None
    for k, n in enumerate(divisions):
        mesh = structured_mesh(element, (n,) * solution.mesh_dim, size=size, origin=origin,
                               mapping=mapping, **element_kwargs)
        if k == 0:
            solution._check_mesh(mesh)
        if solution.time_dependent:
            U = solve_manufactured_transient(mesh, solution, t_end, steps[k], theta=theta,
                                             lumping=lumping, dirichlet=dirichlet,
                                             linear_algebra=linear_algebra).T_history[:, -1]
        else:
            U = solve_manufactured(mesh, solution, dirichlet, linear_algebra)
        for key, v in error_norms(mesh.domain, U, solution, n_points, t=t_end or 0.0).items():
            errors.setdefault(key, []).append(v)
        el = next(iter(mesh.domain.elements.values()))
        if degree is None and el.REFERENCE_CELL != "line":
            degree = interpolation_degree(el)
    return ConvergenceStudy(_name(element), degree, tuple(divisions),
                            1.0 / np.asarray(divisions, dtype=np.float64),
                            {k: np.asarray(v) for k, v in errors.items()}, expected=expected)


def time_convergence_study(element, solution: ManufacturedConduction, divisions: int,
                           n_steps, *, t_end: float, theta: float = 1.0, mapping=None,
                           size=None, origin=None, dirichlet=None, n_points: int | None = None,
                           linear_algebra=None, **element_kwargs) -> ConvergenceStudy:
    """Estudio de convergencia en el tiempo: una malla fija de ``divisions``
    celdas por dirección, integrada hasta ``t_end`` con ``n`` pasos del método
    θ, ``n`` en ``n_steps``.

    Mide sólo el error del esquema en el tiempo. Para ello exige un ``T*``
    lineal en ``x`` (con cualquier dependencia de ``t``), que todo elemento
    isoparamétrico reproduce, y usa la capacidad consistente: con la misma
    regla en la capacidad y en la fuente, la semidiscretización en el espacio
    es exacta (``C·Ṫ* + K·T* = F`` en los nodos) y el error en ``t_end`` es
    sólo el del método θ, de orden 2 con Crank-Nicolson (θ = ½) y 1 con
    cualquier otro θ. La capacidad concentrada (``lumped``) no deja exacta la
    semidiscretización: su error en el espacio pondría un suelo al estudio,
    así que se verifica en :func:`convergence_study`."""
    if not solution.time_dependent:
        raise ValueError("La solución manufacturada no depende del tiempo.")
    if not solution.affine_in_space:
        raise ValueError("El estudio en el tiempo exige un T* lineal en x, que la malla "
                         "reproduce exactamente; con otro, el error en el espacio se mezcla "
                         "con el del tiempo.")
    material = element_kwargs.get("material")
    if material is not None:
        solution._check_material(material)
    errors = {"L2": [], "energy": []}
    degree = order = None
    for m in n_steps:
        mesh = structured_mesh(element, (int(divisions),) * solution.dim, size=size,
                               origin=origin, mapping=mapping, **element_kwargs)
        result = solve_manufactured_transient(mesh, solution, t_end, int(m), theta=theta,
                                              lumping="consistent", dirichlet=dirichlet,
                                              linear_algebra=linear_algebra)
        for key, v in error_norms(mesh.domain, result.T_history[:, -1], solution, n_points,
                                  t=t_end).items():
            errors[key].append(v)
        if degree is None:
            degree = interpolation_degree(next(iter(mesh.domain.elements.values())))
            order = result.order
    return ConvergenceStudy(_name(element), degree, tuple(int(m) for m in n_steps),
                            t_end / np.asarray(n_steps, dtype=np.float64),
                            {k: np.asarray(v) for k, v in errors.items()},
                            expected={"L2": order, "energy": order}, parameter="dt")
