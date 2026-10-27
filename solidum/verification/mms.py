"""Verificación por soluciones manufacturadas (MMS) de la elasticidad lineal (ADR 0022).

Se inventa un desplazamiento suave ``u*(x)``, en 2D o en 3D, y SymPy deriva,
de forma exacta,

    ε* = ∇ˢu*  (Voigt del proyecto, γ_ij = 2·ε_ij),   σ* = C·ε*,
    b* = −∇·σ*,   t*(x, n) = σ*·n.

Con la fuerza de cuerpo ``b*``, la tracción ``t*`` en la parte de Neumann del
contorno y ``u*`` impuesto en los nodos de la parte de Dirichlet, ``u*`` es la
solución exacta. La tracción se evalúa con la normal del borde o de la cara
discretos, así que ``u*`` resuelve exactamente el problema en el dominio de la
propia malla, también con contornos curvos: el error medido es sólo de
discretización.

Al refinar una malla asintóticamente afine, con ``p`` el grado del polinomio
completo que reproducen las funciones de forma (:func:`interpolation_degree`),

    ‖u* − u_h‖_L2 ~ h^(p+1),    ‖u* − u_h‖_E = (∫ (ε*−ε_h)ᵀ C (ε*−ε_h) dΩ)^½ ~ h^p.

Requiere SymPy (``pip install solidum-fem[verification]``). Sirve para
cualquier sólido de desplazamientos ``(ux, uy[, uz])`` que declare
``REFERENCE_CELL``, ``REFERENCE_NODES``, ``BATCH_KINEMATICS`` y funciones de
forma por lotes, incluidos los de un módulo de usuario, con un material
elástico lineal.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from solidum.math.assembly import Assembler
from solidum.math.solvers import LinearSolver
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


class ManufacturedElasticity:
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

    def __init__(self, displacement, C):
        sp = _sympy()
        self.dim = d = len(displacement)
        if d not in (2, 3):
            raise ValueError("displacement: se esperan 2 ó 3 componentes.")
        xs = sp.symbols("x y z", real=True)[:d]
        u = [sp.sympify(e, locals=dict(zip("xyz", xs))) for e in displacement]
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
        # Dos versiones: con math para un punto (las cargas se evalúan punto a
        # punto, y NumPy sobre escalares es ~10 veces más lento) y con NumPy
        # para muchos a la vez (normas del error).
        self._f = {k: sp.lambdify(xs, list(v), "numpy", cse=True)
                   for k, v in self.expressions.items()}
        self._f_point = {k: sp.lambdify(xs, list(v), "math", cse=True)
                         for k, v in self.expressions.items()}

    @classmethod
    def from_material(cls, displacement, material) -> "ManufacturedElasticity":
        return cls(displacement, linear_elastic_matrix(material))

    def _eval(self, name, x) -> np.ndarray:
        """Valores en un punto ``(d,)`` → ``(m,)`` o en varios ``(n, d)`` →
        ``(n, m)``."""
        x = np.asarray(x, dtype=np.float64)
        if x.ndim == 1:
            return np.array(self._f_point[name](*(float(c) for c in x[:self.dim])),
                            dtype=np.float64)
        pts = x[:, :self.dim]
        vals = self._f[name](*pts.T)
        out = np.column_stack([np.broadcast_to(np.asarray(v, dtype=np.float64), (len(pts),))
                               for v in vals])
        return out

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


def solve_manufactured(mesh, solution: ManufacturedElasticity, dirichlet=None,
                       linear_algebra=None) -> np.ndarray:
    """Resuelve la malla con las cargas de ``solution``: ``u*`` en los nodos de
    los lados ``dirichlet`` (por omisión, los lados ``min``) y ``t*`` en los
    demás lados del bloque. Devuelve el vector global de desplazamientos."""
    dom = mesh.domain
    dirichlet = tuple(s for s in mesh.sides if s.endswith("min")) if dirichlet is None else dirichlet
    dofs = _DOFS[:mesh.dim]
    for side in dirichlet:
        for node in mesh.nodes_on(side):
            for dof, value in zip(dofs, solution.displacement(node.coordinates)):
                node.fix_dof(dof, float(value))
    dom.generate_equation_numbers(verbose=False)
    assembler = Assembler(dom)
    F = assembler.assemble_body_load(solution.body_force)
    for side, entities in mesh.boundary.items():
        if side in dirichlet:
            continue
        for el, k in entities:
            load = el.compute_edge_traction if mesh.dim == 2 else el.compute_face_traction
            F[el.get_global_dof_indices()] += load(k, solution.traction)
    kwargs = {} if linear_algebra is None else {"linear_algebra": linear_algebra}
    return LinearSolver(assembler, **kwargs).solve(F)


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


def error_norms(domain, U: np.ndarray, solution: ManufacturedElasticity,
                n_points: int | None = None) -> dict:
    """Normas del error ``u* − u_h``: ``{"L2": ‖e‖_L2, "energy": ‖e‖_E}``
    (en 2D por unidad de espesor), integradas en cada elemento con
    ``n_points`` puntos de Gauss por dirección, más que la regla del elemento
    porque el error no es polinómico. Por omisión ``p + 3``, con ``p`` el
    grado de la interpolación."""
    d = solution.dim
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
        du = solution.displacement(x) - N @ u_e.reshape(-1, d)
        de = solution.strain(x)
        B = np.empty((de.shape[1], u_e.size))
        dV = np.empty(len(wts))
        for g in range(len(wts)):
            dV[g] = cls.BATCH_KINEMATICS(pts[g], coords, B) * wts[g]
            de[g] -= B @ u_e
        e_l2 += np.sum(du * du, axis=1) @ dV
        e_energy += np.einsum("gi,ij,gj->g", de, solution.C, de) @ dV
    return {"L2": float(np.sqrt(e_l2)), "energy": float(np.sqrt(e_energy))}


def interpolation_degree(element, max_degree: int = 6) -> int:
    """Grado ``p`` del polinomio completo que reproducen exactamente las
    funciones de forma en la celda de referencia (``Σ N_i(ξ)·m(ξ_i) = m(ξ)``
    para todo monomio ``m`` de grado ≤ p). Quad4, Tri3, Hex8 y Tet4: 1;
    Quad8, Quad9, Tri6, Hex20, Hex27 y Tet10: 2. Es el que fija los órdenes de
    convergencia teóricos."""
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
    """Resultado de :func:`convergence_study`.

    ``h`` es el tamaño relativo de la malla (``1/divisions``); ``errors``,
    las normas del error por nivel."""
    element: str
    degree: int
    divisions: tuple
    h: np.ndarray
    errors: dict

    def expected_order(self, norm: str) -> int:
        return self.degree + 1 if norm == "L2" else self.degree

    def orders(self, norm: str) -> np.ndarray:
        """Órdenes observados entre niveles consecutivos."""
        e = self.errors[norm]
        return np.log(e[:-1] / e[1:]) / np.log(self.h[:-1] / self.h[1:])

    def table(self) -> str:
        """Tabla de texto: nivel, h, errores y órdenes observados."""
        rows = [f"{self.element} (p = {self.degree})",
                f"{'n':>4} {'h':>9} {'e_L2':>11} {'orden':>6} {'e_E':>11} {'orden':>6}"]
        o_l2, o_en = self.orders("L2"), self.orders("energy")
        for k, n in enumerate(self.divisions):
            a = f"{o_l2[k - 1]:6.3f}" if k else " " * 6
            b = f"{o_en[k - 1]:6.3f}" if k else " " * 6
            rows.append(f"{n:4d} {self.h[k]:9.5f} {self.errors['L2'][k]:11.4e} {a} "
                        f"{self.errors['energy'][k]:11.4e} {b}")
        return "\n".join(rows)


def convergence_study(element, solution: ManufacturedElasticity, divisions, *,
                      mapping=None, size=None, origin=None, dirichlet=None,
                      n_points: int | None = None, linear_algebra=None,
                      **element_kwargs) -> ConvergenceStudy:
    """Estudio de convergencia en mallas estructuradas de ``n`` celdas por
    dirección, ``n`` en ``divisions``, sobre el bloque transformado por
    ``mapping``.

    ``element_kwargs`` va al constructor del elemento (``material``,
    ``thickness``...). Si incluye ``material``, su ``C`` debe ser la de
    ``solution``: con otra, ``u*`` no sería la solución exacta."""
    material = element_kwargs.get("material")
    if material is not None and not np.allclose(linear_elastic_matrix(material), solution.C,
                                                rtol=1e-12, atol=0.0):
        raise ValueError("La C del material no es la de la solución manufacturada.")
    errors = {"L2": [], "energy": []}
    degree = None
    for n in divisions:
        mesh = structured_mesh(element, (n,) * solution.dim, size=size, origin=origin,
                               mapping=mapping, **element_kwargs)
        U = solve_manufactured(mesh, solution, dirichlet, linear_algebra)
        for k, v in error_norms(mesh.domain, U, solution, n_points).items():
            errors[k].append(v)
        if degree is None:
            degree = interpolation_degree(next(iter(mesh.domain.elements.values())))
    name = element if isinstance(element, str) else element.__name__
    return ConvergenceStudy(name, degree, tuple(divisions),
                            1.0 / np.asarray(divisions, dtype=np.float64),
                            {k: np.asarray(v) for k, v in errors.items()})
