"""Verificación por soluciones manufacturadas (MMS) de la elasticidad plana (ADR 0022).

Se inventa un desplazamiento suave ``u*(x)`` y SymPy deriva, de forma exacta,

    ε* = ∇ˢu*  (Voigt del proyecto, γ_xy = 2·ε_xy),   σ* = C·ε*,
    b* = −∇·σ*,   t*(x, n) = σ*·n.

Con la fuerza de cuerpo ``b*``, la tracción ``t*`` en la parte de Neumann del
contorno y ``u*`` impuesto en los nodos de la parte de Dirichlet, ``u*`` es la
solución exacta. La tracción se evalúa con la normal del borde discreto, así
que ``u*`` resuelve exactamente el problema en el dominio de la propia malla,
también con contornos curvos: el error medido es sólo de discretización.

Al refinar una malla asintóticamente afine, con ``p`` el grado del polinomio
completo que reproducen las funciones de forma (:func:`interpolation_degree`),

    ‖u* − u_h‖_L2 ~ h^(p+1),    ‖u* − u_h‖_E = (∫ (ε*−ε_h)ᵀ C (ε*−ε_h) dA)^½ ~ h^p.

Requiere SymPy (``pip install solidum-fem[verification]``). Sirve para
cualquier elemento plano de desplazamientos ``(ux, uy)`` que declare
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


def _sympy():
    try:
        import sympy
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "solidum.verification necesita SymPy: pip install solidum-fem[verification]"
        ) from exc
    return sympy


def linear_elastic_matrix(material, rtol: float = 1e-10) -> np.ndarray:
    """Matriz constitutiva ``C`` (3×3, Voigt 2D) de un material elástico
    lineal, leída de ``compute_state``. Lanza ``ValueError`` si ``σ(ε)`` no es
    ``C·ε`` (material no lineal: la MMS de este módulo es lineal)."""
    _, C, _ = material.compute_state(np.zeros(3), None)
    C = np.array(C, dtype=np.float64)
    eps = 1e-4 * np.array([0.7, -0.4, 1.1])
    sigma, _, _ = material.compute_state(eps, None)
    if np.linalg.norm(np.asarray(sigma) - C @ eps) > rtol * np.linalg.norm(C @ eps):
        raise ValueError(f"{type(material).__name__} no es elástico lineal: "
                         "la verificación por soluciones manufacturadas de este "
                         "módulo exige σ = C·ε.")
    return C


class ManufacturedElasticity2D:
    """Solución manufacturada de la elasticidad plana con ``C`` constante.

    Parameters
    ----------
    displacement : par de expresiones en ``x``, ``y``
        ``u*_x`` y ``u*_y`` como cadenas (``"1e-3*sin(x + 2*y)"``) o
        expresiones de SymPy.
    C : (3, 3)
        Matriz constitutiva en Voigt 2D; :meth:`from_material` la toma del
        material.
    """

    def __init__(self, displacement, C):
        sp = _sympy()
        x, y = sp.symbols("x y", real=True)
        u = [sp.sympify(e, locals={"x": x, "y": y}) for e in displacement]
        if len(u) != 2:
            raise ValueError("displacement: se esperan dos componentes (u_x, u_y).")
        self.C = np.array(C, dtype=np.float64).reshape(3, 3)
        eps = sp.Matrix([sp.diff(u[0], x), sp.diff(u[1], y),
                         sp.diff(u[0], y) + sp.diff(u[1], x)])
        sig = sp.Matrix(self.C.tolist()) * eps
        b = sp.Matrix([-(sp.diff(sig[0], x) + sp.diff(sig[2], y)),
                       -(sp.diff(sig[2], x) + sp.diff(sig[1], y))])
        self.expressions = {"displacement": sp.Matrix(u), "strain": eps,
                            "stress": sig, "body_force": b}
        self._f = {k: sp.lambdify((x, y), list(v), "numpy")
                   for k, v in self.expressions.items()}

    @classmethod
    def from_material(cls, displacement, material) -> "ManufacturedElasticity2D":
        return cls(displacement, linear_elastic_matrix(material))

    def _eval(self, name, x) -> np.ndarray:
        return np.array(self._f[name](float(x[0]), float(x[1])), dtype=np.float64)

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
        """``t*(x, n) = σ*(x)·n``: el campo que se pasa a ``compute_edge_traction``."""
        s = self.stress(x)
        return np.array([s[0] * n[0] + s[2] * n[1], s[2] * n[0] + s[1] * n[1]])


def solve_manufactured(mesh, solution: ManufacturedElasticity2D,
                       dirichlet=("xmin", "ymin"), linear_algebra=None) -> np.ndarray:
    """Resuelve la malla con las cargas de ``solution``: ``u*`` en los nodos de
    los lados ``dirichlet`` y ``t*`` en los demás lados del bloque. Devuelve
    el vector global de desplazamientos."""
    dom = mesh.domain
    for side in dirichlet:
        for node in mesh.nodes_on(side):
            ux, uy = solution.displacement(node.coordinates)
            node.fix_dof("ux", float(ux))
            node.fix_dof("uy", float(uy))
    dom.generate_equation_numbers(verbose=False)
    assembler = Assembler(dom)
    F = assembler.assemble_body_load(solution.body_force)
    for side, edges in mesh.boundary_edges.items():
        if side in dirichlet:
            continue
        for el, edge in edges:
            F[el.get_global_dof_indices()] += el.compute_edge_traction(edge, solution.traction)
    kwargs = {} if linear_algebra is None else {"linear_algebra": linear_algebra}
    return LinearSolver(assembler, **kwargs).solve(F)


def _reference_rule(cell: str, n: int):
    """Regla de Gauss de ``n`` puntos por dirección en la celda de referencia:
    producto tensorial en ``[-1, 1]²``; en el triángulo, producto colapsado
    (Duffy), exacto hasta grado ``2n − 2``."""
    a, w = np.polynomial.legendre.leggauss(n)
    if cell == "quad":
        X, Y = np.meshgrid(a, a, indexing="ij")
        return np.column_stack([X.ravel(), Y.ravel()]), np.outer(w, w).ravel()
    if cell == "tri":
        s, ws = (a + 1.0) / 2.0, w / 2.0
        S, T = np.meshgrid(s, s, indexing="ij")
        W = np.outer(ws, ws) * (1.0 - S)
        return np.column_stack([S.ravel(), (T * (1.0 - S)).ravel()]), W.ravel()
    raise ValueError(f"Celda de referencia {cell!r} no admitida (quad, tri).")


def error_norms(domain, U: np.ndarray, solution: ManufacturedElasticity2D,
                n_points: int = 6) -> dict:
    """Normas del error ``u* − u_h``: ``{"L2": ‖e‖_L2, "energy": ‖e‖_E}``,
    por unidad de espesor, integradas con ``n_points`` puntos de Gauss por
    dirección en cada elemento (más que la regla del elemento: el error no es
    polinómico)."""
    e_l2 = 0.0
    e_energy = 0.0
    for el in domain.elements.values():
        cls = type(el)
        pts, wts = _reference_rule(cls.REFERENCE_CELL, n_points)
        coords = np.ascontiguousarray(el.get_coordinate_matrix(ndim=2))
        u_e = np.asarray(el.get_local_displacements(U), dtype=np.float64)
        N = el.batch_shape_functions(pts)
        x = N @ coords
        u_h = N @ u_e.reshape(-1, 2)
        B = np.empty((3, u_e.size))
        for g in range(len(wts)):
            dA = cls.BATCH_KINEMATICS(pts[g], coords, B) * wts[g]
            du = solution.displacement(x[g]) - u_h[g]
            de = solution.strain(x[g]) - B @ u_e
            e_l2 += du @ du * dA
            e_energy += de @ solution.C @ de * dA
    return {"L2": float(np.sqrt(e_l2)), "energy": float(np.sqrt(e_energy))}


def interpolation_degree(element, max_degree: int = 6) -> int:
    """Grado ``p`` del polinomio completo que reproducen exactamente las
    funciones de forma en la celda de referencia (``Σ N_i(ξ)·m(ξ_i) = m(ξ)``
    para todo monomio ``m`` de grado ≤ p). Quad4 y Tri3: 1; Quad8, Quad9 y
    Tri6: 2. Es el que fija los órdenes de convergencia teóricos."""
    nodes = np.asarray(type(element).REFERENCE_NODES, dtype=np.float64)
    rng = np.random.default_rng(0)
    pts = rng.uniform(0.0, 1.0, (12, 2))
    if type(element).REFERENCE_CELL == "quad":
        pts = 2.0 * pts - 1.0
    else:
        fold = pts.sum(axis=1) > 1.0
        pts[fold] = 1.0 - pts[fold]
    N = element.batch_shape_functions(pts)
    degree = 0
    for k in range(1, max_degree + 1):
        for a in range(k + 1):
            exact = pts[:, 0] ** a * pts[:, 1] ** (k - a)
            interp = N @ (nodes[:, 0] ** a * nodes[:, 1] ** (k - a))
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


def convergence_study(element, solution: ManufacturedElasticity2D, divisions, *,
                      mapping=None, size=(1.0, 1.0), origin=(0.0, 0.0),
                      dirichlet=("xmin", "ymin"), n_points: int = 6,
                      linear_algebra=None, **element_kwargs) -> ConvergenceStudy:
    """Estudio de convergencia en mallas estructuradas de ``n × n`` celdas,
    ``n`` en ``divisions``, sobre el bloque transformado por ``mapping``.

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
        mesh = structured_mesh(element, (n, n), size=size, origin=origin,
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
