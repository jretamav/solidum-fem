"""Seguimiento global de la trayectoria de grieta (spec ``GlobalCrackTracking``).

Decide, al inicio de cada paso y con el estado convergido, qué elementos
``CST_Embedded2D`` se agrietan y con qué partición (el nodo solitario). La
partición la dicta la trayectoria (decisión del usuario del 2026-11-10): cada
grieta ``k`` tiene un campo escalar ``θ_k``, lineal por elemento, que minimiza

    Π(θ) = ½∫(t·∇θ)² dΩ + ½ε∫|∇θ|² dΩ ,

con ``t`` el campo de direcciones (tangente fijada en los elementos
agrietados; perpendicular a ``σ_I`` en los intactos) y condiciones de
Dirichlet en los nodos de la grieta ya formada (distancia con signo a su
recta). La trayectoria es la isolínea ``θ_k = 0`` y el lado de cada nodo es
el signo de ``θ_k``, el mismo para todos los elementos que lo comparten, de
modo que la partición es coherente por construcción.

Regla del vértice (decisión del usuario): un elemento cortado cuyo nodo
aislado dé ``cos ψ ≤ 0`` no se activa así; ese nodo cambia de lado y la
trayectoria pasa por el vértice. Si el nodo ya pertenece a la grieta formada
(no puede cambiar), cambia de lado otro nodo libre del elemento, el que deje
``cos ψ > 0``. Cada nodo cambia de lado como mucho una vez por paso.

Referencias: Oliver et al. (2004, 2008) para el problema auxiliar; Riccardi
et al. (2017) para la regularización y las condiciones de Dirichlet en la
grieta formada; Parvaneh y Foster (2016) para la actualización al inicio del
paso. Formulación completa en
``docs/user/discontinuities/specs/GlobalCrackTracking.md``.
"""
from __future__ import annotations

import weakref
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

from solidum.logging import get_logger
from solidum.user.discontinuities.constants import (
    TRACKING_EPSILON,
    TRACKING_MIN_DISTANCE,
)
from solidum.user.discontinuities.registry import CrackTrackerRegistry

_log = get_logger("user.discontinuities")

_EDGES = ((0, 1), (1, 2), (2, 0))


# =========================================================================
# Geometría del triángulo y problema auxiliar (funciones puras)
# =========================================================================


def triangle_gradients(X: np.ndarray, conn: np.ndarray):
    """Gradientes de las funciones de forma lineales y áreas.

    Returns
    -------
    grads : ndarray (ne, 3, 2)
        ``∇N_i`` de cada nodo de cada triángulo. ``∇N_i`` es normal al lado
        opuesto al nodo ``i`` y apunta hacia él, con módulo ``1/h_i``.
    area : ndarray (ne,)
    """
    x = X[conn]                                   # (ne, 3, 2)
    x0, x1, x2 = x[:, 0], x[:, 1], x[:, 2]
    two_a = ((x1[:, 0] - x0[:, 0]) * (x2[:, 1] - x0[:, 1])
             - (x2[:, 0] - x0[:, 0]) * (x1[:, 1] - x0[:, 1]))
    if np.any(np.abs(two_a) <= 0.0):
        raise ValueError("triángulo degenerado (área nula) en la región de seguimiento")
    dx = np.stack([x1[:, 1] - x2[:, 1], x2[:, 1] - x0[:, 1], x0[:, 1] - x1[:, 1]], axis=1)
    dy = np.stack([x2[:, 0] - x1[:, 0], x0[:, 0] - x2[:, 0], x1[:, 0] - x0[:, 0]], axis=1)
    grads = np.stack([dx, dy], axis=2) / two_a[:, None, None]
    return grads, 0.5 * np.abs(two_a)


def tracking_matrix(grads: np.ndarray, area: np.ndarray, tangents: np.ndarray,
                    n_nodes: int, conn: np.ndarray, epsilon: float) -> sp.csr_matrix:
    """Matriz ``K = Σ_e A_e ∇Nᵀ(t_e t_eᵀ + ε I)∇N`` del problema auxiliar (§6).

    ``tangents`` (ne, 2) unitarias; su signo no importa (campo de líneas).
    """
    t = tangents / np.linalg.norm(tangents, axis=1, keepdims=True)
    D = np.einsum("ei,ej->eij", t, t) + epsilon * np.eye(2)[None]
    ke = np.einsum("e,eai,eij,ebj->eab", area, grads, D, grads)
    rows = np.repeat(conn, 3, axis=1).ravel()
    cols = np.tile(conn, (1, 3)).ravel()
    return sp.csr_matrix((ke.ravel(), (rows, cols)), shape=(n_nodes, n_nodes))


def solve_tracking_field(K: sp.csr_matrix, nodes: np.ndarray, fixed: dict) -> np.ndarray:
    """Resuelve ``θ`` en ``nodes`` (índices) con Dirichlet ``fixed`` {índice: valor}.

    Los nodos fuera de ``nodes`` quedan en 0. Los de ``nodes`` deben formar
    una región conexa con al menos un valor prescrito (si no, ``K`` restringida
    es singular: el regularizador sólo fija ``θ`` salvo una constante).
    """
    theta = np.zeros(K.shape[0])
    if not fixed:
        raise ValueError("el campo de seguimiento necesita al menos un valor prescrito")
    d = np.fromiter(fixed.keys(), dtype=np.int64)
    vals = np.fromiter(fixed.values(), dtype=float)
    theta[d] = vals
    free = np.setdiff1d(nodes, d, assume_unique=False)
    if free.size:
        Kf = K[free]
        theta[free] = spsolve(Kf[:, free].tocsc(), -(Kf[:, d] @ vals))
    return theta


def partition(xs: np.ndarray, th: np.ndarray, normal: np.ndarray, inward: np.ndarray):
    """Partición de un triángulo por el signo de ``θ`` en sus nodos.

    Parameters
    ----------
    xs : (3, 2) coordenadas; th : (3,) valores de ``θ``;
    normal : (2,) normal de Rankine (cualquier signo);
    inward : (3, 2) normal unitaria de cada lado opuesto, hacia su nodo.

    Returns
    -------
    ``None`` si el elemento no está cortado; si no, ``(lone, n, x_m, cos ψ)``:
    el nodo solitario, la normal orientada hacia él (respecto del punto medio
    ``x_m`` de los dos cruces de la isolínea con sus lados) y ``cos ψ = n·n_e``.
    """
    s = th >= 0.0
    count = int(s.sum())
    if count in (0, 3):
        return None
    lone = int(np.argmax(s)) if count == 1 else int(np.argmin(s))
    crossings = [xs[lone] + th[lone] / (th[lone] - th[j]) * (xs[j] - xs[lone])
                 for j in range(3) if j != lone]
    x_m = 0.5 * (crossings[0] + crossings[1])
    n = normal / np.linalg.norm(normal)
    if n @ (xs[lone] - x_m) < 0.0:
        n = -n
    return lone, n, x_m, float(n @ inward[lone])


def root_partition(xs: np.ndarray, normal: np.ndarray):
    """Partición de la raíz de una grieta (§4.2 y §5).

    La recta de normal ``normal`` por el centroide aísla un nodo. Si ese nodo
    da ``cos ψ ≤ 0``, la recta pasa al otro lado del nodo intermedio y aísla
    el del extremo opuesto: de las dos particiones posibles con esa normal,
    al menos una tiene ``cos ψ > 0``.

    Returns
    -------
    ``(lone, n, point)``: nodo solitario, normal unitaria orientada hacia él y
    un punto de la recta.
    """
    n = np.asarray(normal, dtype=float)
    n = n / np.linalg.norm(n)
    grads, _ = triangle_gradients(xs, np.array([[0, 1, 2]]))
    inward = grads[0] / np.linalg.norm(grads[0], axis=1, keepdims=True)
    s = xs @ n
    centroid = xs.mean(axis=0)
    i_max, i_min = int(np.argmax(s)), int(np.argmin(s))
    i_mid = 3 - i_max - i_min
    at_centroid = (i_max, n) if int((s > s.mean()).sum()) == 1 else (i_min, -n)
    shifted = (i_min, -n) if at_centroid[0] == i_max else (i_max, n)
    for option in (at_centroid, shifted):
        lone, nn = option
        if nn @ inward[lone] > 0.0:
            if option is at_centroid:
                return lone, nn, centroid
            s_line = 0.5 * (nn @ xs[lone] + nn @ xs[i_mid])
            return lone, nn, centroid + nn * (s_line - nn @ centroid)
    raise RuntimeError("partición de raíz sin cos ψ > 0 (triángulo degenerado)")


def isoline_segments(X: np.ndarray, conn: np.ndarray, theta: np.ndarray,
                     elements=None) -> list:
    """Segmentos de la isolínea ``θ = 0`` por elemento: ``[(e, p, q)]``."""
    out = []
    for e in (range(len(conn)) if elements is None else elements):
        idx = conn[e]
        th = theta[idx]
        s = th >= 0.0
        if s.all() or not s.any():
            continue
        pts = []
        for a, b in _EDGES:
            if s[a] != s[b]:
                w = th[a] / (th[a] - th[b])
                pts.append(X[idx[a]] + w * (X[idx[b]] - X[idx[a]]))
        out.append((e, pts[0], pts[1]))
    return out


# =========================================================================
# Grieta y algoritmo de seguimiento
# =========================================================================


@dataclass
class Crack:
    """Una grieta del seguimiento global.

    Attributes
    ----------
    index : int
        Índice de la grieta (orden de nacimiento).
    root : int or None
        Id del elemento raíz (``None`` hasta la primera activación de una
        trayectoria prescrita).
    elements : list[int]
        Ids de sus elementos agrietados, en orden de activación.
    fixed : dict[int, float]
        Valores de Dirichlet de su campo de seguimiento por id de nodo: la
        distancia con signo a la recta de la grieta en el elemento que fijó
        el nodo (mínimo ``TRACKING_MIN_DISTANCE·h_e``). El signo es el lado.
    """
    index: int
    root: int | None
    elements: list = field(default_factory=list)
    fixed: dict = field(default_factory=dict)

    def side(self, node_id: int) -> int:
        """Lado (+1/−1) del nodo respecto de la grieta; 0 si no está fijado."""
        v = self.fixed.get(node_id)
        return 0 if v is None else (1 if v >= 0.0 else -1)


@CrackTrackerRegistry.register
class GlobalCrackTracking:
    """Seguimiento global de grietas para ``CST_Embedded2D`` (spec
    ``GlobalCrackTracking``).

    Parameters
    ----------
    epsilon : float, default 1e-6
        Regularización adimensional del problema auxiliar (§2).

    Notes
    -----
    Los elementos se dan de alta al construirse (``crack_tracker=`` en el
    elemento o, sin él, un seguimiento por omisión por región conexa). El
    solver llama a ``prepare_step`` de cada elemento al inicio de cada paso;
    el primero de cada ronda dispara :meth:`update`, que es idempotente para
    el mismo campo convergido (se itera hasta que no cambia nada).
    """

    def __init__(self, epsilon: float = TRACKING_EPSILON):
        epsilon = float(epsilon)
        if not epsilon > 0.0:
            raise ValueError(f"GlobalCrackTracking: epsilon debe ser > 0 (recibido {epsilon}).")
        self.epsilon = epsilon
        self.cracks: list[Crack] = []
        self.unactivated_warnings = 0
        self.side_conflict_warnings = 0
        self._members: list = []
        self._member_keys: set = set()
        self._built = False
        self._called: set = set()
        self._last_U = None
        self._last_copy = None

    # ------------------------------------------------------------------
    # Alta de elementos y disparo por ronda
    # ------------------------------------------------------------------

    @property
    def elements(self) -> tuple:
        return tuple(self._members)

    def add(self, element) -> None:
        """Da de alta un elemento en la región de seguimiento."""
        if id(element) in self._member_keys:
            return
        self._members.append(element)
        self._member_keys.add(id(element))
        self._built = False

    def prepare_step(self, element, U_committed: np.ndarray) -> None:
        """Gancho de inicio de paso de un elemento de la región.

        La primera llamada de cada ronda (o una con otro campo) ejecuta
        :meth:`update`; las demás de la misma ronda no hacen nada.
        """
        U = np.asarray(U_committed)
        same = (U is self._last_U
                or (self._last_copy is not None and U.shape == self._last_copy.shape
                    and np.array_equal(U, self._last_copy)))
        key = id(element)
        if same and key not in self._called:
            self._called.add(key)
            return
        self.update(U)
        self._last_U = U
        self._last_copy = U.copy()
        self._called = {key}

    # ------------------------------------------------------------------
    # Topología de la región
    # ------------------------------------------------------------------

    def _build(self) -> None:
        els = sorted(self._members, key=lambda e: e.id)
        self._els = els
        self._pos = {e.id: p for p, e in enumerate(els)}
        node_ids = []
        coords = {}
        for e in els:
            for nd in e.nodes:
                if nd.id not in coords:
                    coords[nd.id] = np.asarray(nd.coordinates, dtype=float)[:2]
                    node_ids.append(nd.id)
        self._node_ids = np.array(node_ids, dtype=np.int64)
        self._nidx = {nid: i for i, nid in enumerate(node_ids)}
        self._X = np.array([coords[nid] for nid in node_ids]).reshape(-1, 2)
        self._conn = np.array([[self._nidx[nd.id] for nd in e.nodes] for e in els],
                              dtype=np.int64).reshape(-1, 3)
        self._grads, self._area = triangle_gradients(self._X, self._conn)
        self._inward = self._grads / np.linalg.norm(self._grads, axis=2, keepdims=True)
        self._delta = TRACKING_MIN_DISTANCE * np.sqrt(2.0 * self._area)
        # Matriz B del CST en Voigt [ε_xx, ε_yy, γ_xy] y umbral de cada elemento.
        B = np.zeros((len(els), 3, 6))
        B[:, 0, 0::2] = self._grads[:, :, 0]
        B[:, 1, 1::2] = self._grads[:, :, 1]
        B[:, 2, 0::2] = self._grads[:, :, 1]
        B[:, 2, 1::2] = self._grads[:, :, 0]
        self._B = B
        self._sigma_t0 = np.array([float(getattr(e.cohesive_material, "sigma_t0", None) or 0.0)
                                   for e in els])

        self._edge_elems: dict = {}
        for p, tri in enumerate(self._conn):
            for a, b in _EDGES:
                key = (min(tri[a], tri[b]), max(tri[a], tri[b]))
                self._edge_elems.setdefault(key, []).append(p)

        # Componentes conexas por nodos compartidos (el problema auxiliar de
        # cada grieta se resuelve en la suya).
        parent = list(range(len(node_ids)))

        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i
        for tri in self._conn:
            r0 = find(tri[0])
            for j in tri[1:]:
                rj = find(j)
                if rj != r0:
                    parent[rj] = r0
        roots = np.array([find(i) for i in range(len(node_ids))])
        self._comp_nodes = {r: np.nonzero(roots == r)[0] for r in np.unique(roots)}
        self._node_comp = roots

        # Estado derivado de las grietas existentes (los ids son la verdad).
        self._fixed_idx = [{self._nidx[nid]: v for nid, v in c.fixed.items() if nid in self._nidx}
                           for c in self.cracks]
        self._crack_pos = [{self._pos[i] for i in c.elements if i in self._pos} for c in self.cracks]
        self._cracked_count = np.zeros(len(node_ids), dtype=np.int64)
        for pos in self._crack_pos:
            for p in pos:
                self._cracked_count[self._conn[p]] += 1
        self._built = True

    def _adopt_external(self) -> None:
        """Incorpora a una grieta los elementos agrietados fuera del
        seguimiento (activación explícita en Python): se une a la grieta con
        la que comparte nodos, o nace una nueva."""
        known = set()
        for c in self.cracks:
            known.update(c.elements)
        for p, e in enumerate(self._els):
            ds = e.discontinuity_state
            if ds is None or e.id in known:
                continue
            xs = self._X[self._conn[p]]
            th = ds.normal @ (xs - ds.centroid).T
            sign = -np.ones(3)
            sign[ds.solitary_node] = 1.0
            th = sign * np.maximum(np.abs(th), self._delta[p])
            k, th = self._adoption_crack(p, th)
            ds.crack = k
            self._record(k, p, th)
            known.add(e.id)

    def _adoption_crack(self, p: int, th: np.ndarray):
        """Grieta a la que se une un elemento agrietado fuera del seguimiento:
        la primera con la que comparte nodos (con ``θ`` orientado como ella) o
        una nueva con raíz en él."""
        for k, fixed in enumerate(self._fixed_idx):
            shared = [j for j in range(3) if self._conn[p][j] in fixed]
            if shared:
                j = shared[0]
                if np.sign(fixed[self._conn[p][j]]) != np.sign(th[j]):
                    th = -th
                return k, th
        return self._new_crack(self._els[p].id), th

    def _new_crack(self, root_id: int | None) -> int:
        k = len(self.cracks)
        self.cracks.append(Crack(index=k, root=root_id))
        self._fixed_idx.append({})
        self._crack_pos.append(set())
        self._flips.append({})
        self._theta.append(None)
        self._chain.append(set())
        self._tips.append(set())
        return k

    def _record(self, k: int, p: int, th: np.ndarray) -> None:
        """Fija los nodos de ``p`` en la grieta ``k`` con los valores ``th``
        (sin sobreescribir los ya fijados) y lo añade a sus elementos."""
        crack = self.cracks[k]
        for j, li in enumerate(self._conn[p]):
            self._flips[k].pop(li, None)
            if li in self._fixed_idx[k]:
                continue
            self._fixed_idx[k][li] = float(th[j])
            crack.fixed[int(self._node_ids[li])] = float(th[j])
        if crack.root is None:
            crack.root = self._els[p].id
        crack.elements.append(self._els[p].id)
        self._crack_pos[k].add(p)
        self._cracked_count[self._conn[p]] += 1

    # ------------------------------------------------------------------
    # Actualización (§4)
    # ------------------------------------------------------------------

    def update(self, U_committed: np.ndarray) -> None:
        """Activa los elementos que corresponda con el estado convergido."""
        if not self._built:
            self._build()
        els = self._els
        ne = len(els)
        self._flips = [{} for _ in self.cracks]
        self._theta = [None for _ in self.cracks]
        self._chain = [set() for _ in self.cracks]
        self._tips = [set() for _ in self.cracks]
        self._adopt_external()

        intact = np.array([e.discontinuity_state is None for e in els], dtype=bool)
        sigma_I, self._nrank = self._principal_stresses(np.asarray(U_committed), intact)
        ratio = np.where(intact & (self._sigma_t0 > 0.0),
                         sigma_I / np.where(self._sigma_t0 > 0.0, self._sigma_t0, 1.0), -np.inf)
        if not np.any(ratio > 1.0):
            return                  # nada que activar: no hace falta el campo
        tangents = np.column_stack([-self._nrank[:, 1], self._nrank[:, 0]])
        for p in np.nonzero(~intact)[0]:
            tangents[p] = els[p].discontinuity_state.tangent
        self._assemble_field(tangents)
        self._prepare_cracks()

        for k in range(len(self.cracks)):
            self._solve(k)
            self._vertex_pass(k)

        # La propagación tiene prioridad: una grieta que avanza reclama los
        # elementos de su trayectoria antes de que nazca otra delante de su
        # punta. Sólo cuando ninguna punta puede avanzar nace una grieta (la
        # del elemento más solicitado que no comparte nodos con una grieta) y
        # se vuelve a propagar. Un elemento que comparte un nodo con un
        # elemento agrietado sin estar en una punta no se activa en este paso
        # (exclusión, §4.3).
        order = [int(p) for p in np.argsort(-ratio, kind="stable") if ratio[p] > 1.0]
        warned: set = set()
        changed = True
        while changed:
            changed = False
            for p in order:
                if els[p].discontinuity_state is not None:
                    continue
                k = next((kk for kk in range(len(self.cracks)) if p in self._tips[kk]), None)
                if k is not None and self._propagate(k, p, warned):
                    changed = True
            if changed:
                continue
            for p in order:
                if els[p].discontinuity_state is None and self._can_nucleate(p):
                    self._nucleate(p)
                    changed = True
                    break

    def _principal_stresses(self, U: np.ndarray, intact: np.ndarray):
        """``(σ_I, n)`` de los elementos intactos con el campo ``U`` (criterio
        de Rankine). Vectorizado para bulks lineales (``batch_matrix``);
        elemento a elemento para los demás."""
        els = self._els
        ne = len(els)
        sigma_I = np.full(ne, -np.inf)
        normals = np.tile([1.0, 0.0], (ne, 1))
        idx = np.nonzero(intact)[0]
        if idx.size == 0:
            return sigma_I, normals
        dofs = np.array([els[p].get_global_dof_indices() for p in idx])
        strain = np.einsum("eij,ej->ei", self._B[idx], U[dofs])
        stress = np.empty_like(strain)
        groups: dict = {}
        for row, p in enumerate(idx):
            groups.setdefault(id(els[p].material), []).append(row)
        for rows in groups.values():
            material = els[idx[rows[0]]].material
            try:
                C = np.asarray(material.batch_matrix(), dtype=float)
            except NotImplementedError:
                C = None
            if C is not None and C.shape == (3, 3):
                stress[rows] = strain[rows] @ C.T
            else:
                for row in rows:
                    p = idx[row]
                    stress[row] = els[p].material.compute_state(strain[row], els[p].state.vars[0])[0]
        sxx, syy, sxy = stress[:, 0], stress[:, 1], stress[:, 2]
        sigma_I[idx] = 0.5 * (sxx + syy) + np.hypot(0.5 * (sxx - syy), sxy)
        theta_p = 0.5 * np.arctan2(2.0 * sxy, sxx - syy)
        normals[idx] = np.column_stack([np.cos(theta_p), np.sin(theta_p)])
        return sigma_I, normals

    # Ganchos que redefine una trayectoria prescrita.

    def _assemble_field(self, tangents: np.ndarray) -> None:
        self._K = tracking_matrix(self._grads, self._area, tangents, len(self._node_ids),
                                  self._conn, self.epsilon)

    def _prepare_cracks(self) -> None:
        """Grietas que deben existir antes de activar (ninguna en el global)."""

    def _can_nucleate(self, p: int) -> bool:
        return not self._cracked_count[self._conn[p]].any()

    def _fixed_values(self, p: int, n: np.ndarray, x_m: np.ndarray, th: np.ndarray) -> np.ndarray:
        """Valores de Dirichlet de los nodos de un elemento recién activado:
        distancia con signo a la recta de la grieta en él (§2)."""
        xs = self._X[self._conn[p]]
        return np.where(th >= 0.0, 1.0, -1.0) * np.maximum(np.abs((xs - x_m) @ n), self._delta[p])

    def _solve(self, k: int) -> None:
        fixed = dict(self._fixed_idx[k])
        fixed.update(self._flips[k])
        any_node = next(iter(self._fixed_idx[k]))
        nodes = self._comp_nodes[self._node_comp[any_node]]
        self._theta[k] = solve_tracking_field(self._K, nodes, fixed)
        self._update_chain(k)

    def _update_chain(self, k: int) -> None:
        """Elementos intactos cortados por la trayectoria de ``k`` conectados a
        la grieta a través de lados cruzados, y las puntas entre ellos."""
        s = self._theta[k] >= 0.0
        cracked = self._crack_pos[k]
        chain: set = set()
        tips: set = set()
        queue = list(cracked)
        visited = set(queue)
        while queue:
            p = queue.pop()
            tri = self._conn[p]
            for a, b in _EDGES:
                if s[tri[a]] == s[tri[b]]:
                    continue
                key = (min(tri[a], tri[b]), max(tri[a], tri[b]))
                for q in self._edge_elems[key]:
                    if q in visited or self._els[q].discontinuity_state is not None:
                        continue
                    visited.add(q)
                    chain.add(q)
                    queue.append(q)
        for q in chain:
            tri = self._conn[q]
            for a, b in _EDGES:
                if s[tri[a]] == s[tri[b]]:
                    continue
                key = (min(tri[a], tri[b]), max(tri[a], tri[b]))
                if any(r in cracked for r in self._edge_elems[key]):
                    tips.add(q)
                    break
        self._chain[k] = chain
        self._tips[k] = tips

    def _partition(self, p: int, th: np.ndarray):
        return partition(self._X[self._conn[p]], th, self._nrank[p], self._inward[p])

    def _is_free(self, k: int, li: int) -> bool:
        return li not in self._fixed_idx[k] and li not in self._flips[k]

    def _choose_flip(self, k: int, p: int, th: np.ndarray, lone: int):
        """Nodo que cambia de lado por la regla del vértice (§5): el aislado si
        está libre; si no, el nodo libre que deje ``cos ψ > 0``."""
        tri = self._conn[p]
        delta = self._delta[p]
        if self._is_free(k, tri[lone]):
            return tri[lone], (-delta if th[lone] >= 0.0 else delta)
        best = None
        for j in range(3):
            if j == lone or not self._is_free(k, tri[j]):
                continue
            th2 = th.copy()
            th2[j] = -delta if th[j] >= 0.0 else delta
            part = self._partition(p, th2)
            if part is not None and part[3] > 0.0 and (best is None or part[3] > best[2]):
                best = (tri[j], th2[j], part[3])
        return None if best is None else best[:2]

    def _vertex_pass(self, k: int) -> None:
        """Aplica la regla del vértice a los elementos de la trayectoria
        prevista de ``k`` con ``cos ψ ≤ 0``, y vuelve a resolver."""
        for _ in range(len(self._node_ids) + 1):
            theta = self._theta[k].copy()
            flipped = False
            for p in sorted(self._chain[k]):
                th = theta[self._conn[p]]
                part = self._partition(p, th)
                if part is None or part[3] > 0.0:
                    continue
                choice = self._choose_flip(k, p, th, part[0])
                if choice is not None:
                    self._flips[k][choice[0]] = choice[1]
                    theta[choice[0]] = choice[1]
                    flipped = True
            if not flipped:
                return
            self._solve(k)

    def _propagate(self, k: int, p: int, warned: set) -> bool:
        th = self._theta[k][self._conn[p]]
        part = self._partition(p, th)
        if part is None:
            return False
        lone, n, x_m, cos_psi = part
        if cos_psi <= 0.0:
            choice = self._choose_flip(k, p, th, lone)
            if choice is not None:
                self._flips[k][choice[0]] = choice[1]
                self._solve(k)
                self._vertex_pass(k)
                return True
            if p not in warned:
                warned.add(p)
                self.unactivated_warnings += 1
                _log.warning(
                    f"GlobalCrackTracking: el elemento {self._els[p].id} está en la punta "
                    f"de la grieta {k} con σ_I > σ_t0, pero su partición da cos ψ = "
                    f"{cos_psi:.3f} ≤ 0 y ningún nodo libre puede cambiar de lado; no se "
                    f"activa en este paso."
                )
            return False
        self._activate(k, p, lone, n, x_m, th)
        return True

    def _activate(self, k: int, p: int, lone: int, n: np.ndarray, x_m: np.ndarray,
                  th: np.ndarray) -> None:
        self._els[p].activate(n, lone, x_m, crack=k)
        self._record(k, p, self._fixed_values(p, n, x_m, th))
        for q in self._edge_conflicts(p):
            self.side_conflict_warnings += 1
            _log.warning(
                f"GlobalCrackTracking: el elemento {self._els[p].id} (grieta {k}) y su vecino "
                f"{self._els[q].id} (grieta {self._els[q].discontinuity_state.crack}) reparten "
                f"de forma distinta los nodos de su lado común: las grietas se encontraron y su "
                f"unión está fuera de alcance."
            )
        self._solve(k)
        self._vertex_pass(k)

    def _edge_conflicts(self, p: int) -> list:
        """Vecinos agrietados de ``p`` que reparten de otro modo los nodos de un
        lado común (uno los deja del mismo lado y el otro no)."""
        def same_side(q, a, b):
            lone = self._conn[q][self._els[q].discontinuity_state.solitary_node]
            return a != lone and b != lone
        out = []
        tri = self._conn[p]
        for a, b in _EDGES:
            na, nb = tri[a], tri[b]
            for q in self._edge_elems[(min(na, nb), max(na, nb))]:
                if q != p and self._els[q].discontinuity_state is not None \
                        and same_side(q, na, nb) != same_side(p, na, nb):
                    out.append(q)
        return out

    def side_conflicts(self) -> list:
        """Pares ``(id, id)`` de elementos agrietados vecinos con repartos
        incoherentes de su lado común. Vacía por construcción dentro de una
        grieta; sólo aparece si dos grietas se encuentran (fuera de alcance)."""
        if not self._built:
            return []
        pairs = set()
        for p, e in enumerate(self._els):
            if e.discontinuity_state is not None:
                for q in self._edge_conflicts(p):
                    pairs.add(tuple(sorted((e.id, self._els[q].id))))
        return sorted(pairs)

    def _nucleate(self, p: int) -> None:
        xs = self._X[self._conn[p]]
        lone, n, point = root_partition(xs, self._nrank[p])
        k = self._new_crack(self._els[p].id)
        th = (xs - point) @ n
        self._activate(k, p, lone, n, point, th)

    # ------------------------------------------------------------------
    # Consulta (posproceso y pruebas)
    # ------------------------------------------------------------------

    def field(self, k: int) -> dict:
        """Último campo ``θ_k`` calculado, por id de nodo (vacío si aún no hay)."""
        thetas = getattr(self, "_theta", [])
        theta = thetas[k] if k < len(thetas) else None
        if theta is None:
            return {}
        return {int(nid): float(theta[i]) for i, nid in enumerate(self._node_ids)}


def polyline_signed_distance(points: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Distancia con signo de los puntos ``X`` a la poligonal ``points``,
    prolongada por sus extremos: positiva a la izquierda del sentido de
    recorrido. El signo es el del tramo más cercano."""
    best = np.full(len(X), np.inf)
    out = np.zeros(len(X))
    last = len(points) - 2
    for i, (a, b) in enumerate(zip(points[:-1], points[1:])):
        d = b - a
        w = ((X - a) @ d) / float(d @ d)
        w = np.clip(w, -np.inf if i == 0 else 0.0, np.inf if i == last else 1.0)
        r = X - (a + w[:, None] * d)
        dist = np.hypot(r[:, 0], r[:, 1])
        side = np.where(d[0] * (X[:, 1] - a[1]) - d[1] * (X[:, 0] - a[0]) >= 0.0, 1.0, -1.0)
        closer = dist < best
        best[closer] = dist[closer]
        out[closer] = side[closer] * dist[closer]
    return out


@CrackTrackerRegistry.register
class PrescribedCrackPath(GlobalCrackTracking):
    """Trayectoria de grieta prescrita (spec ``PrescribedCrackPath``).

    Variante de :class:`GlobalCrackTracking` para los modelos en que la
    trayectoria es un dato, como la banda de elementos embebidos de la grieta
    prescrita de la tesis (ensayo de Van Vliet, Retama 2010 §8.1): ``θ`` es la
    distancia con signo a la poligonal ``points`` en lugar de la solución del
    problema auxiliar. Una sola grieta; se activa todo elemento cortado por la
    trayectoria con ``σ_I > σ_t0`` (criterio y normal de Rankine), sin
    propagación ni nucleación. La regla del vértice es la misma.

    Parameters
    ----------
    points : array_like (m, 2), m ≥ 2
        Vértices de la poligonal, prolongada por sus extremos.
    """

    def __init__(self, points):
        super().__init__()
        P = np.asarray(points, dtype=float)
        if P.ndim != 2 or P.shape[1] != 2 or len(P) < 2:
            raise ValueError("PrescribedCrackPath: 'points' debe ser una lista de al menos "
                             "dos puntos [x, y].")
        if np.any(np.hypot(*np.diff(P, axis=0).T) <= 0.0):
            raise ValueError("PrescribedCrackPath: dos puntos consecutivos coinciden.")
        self.points = P

    def _build(self) -> None:
        super()._build()
        self._theta_path = polyline_signed_distance(self.points, self._X)

    def _assemble_field(self, tangents: np.ndarray) -> None:
        """``θ`` es dato: no hay problema auxiliar que ensamblar."""

    def _prepare_cracks(self) -> None:
        if not self.cracks:
            self._new_crack(None)

    def _adoption_crack(self, p: int, th: np.ndarray):
        self._prepare_cracks()
        path = self._theta_path[self._conn[p]]
        j = int(np.argmax(np.abs(path)))
        return 0, (th if np.sign(path[j]) == np.sign(th[j]) else -th)

    def _can_nucleate(self, p: int) -> bool:
        return False

    def _solve(self, k: int) -> None:
        theta = self._theta_path.copy()
        for li, v in self._fixed_idx[k].items():
            theta[li] = v
        for li, v in self._flips[k].items():
            theta[li] = v
        self._theta[k] = theta
        self._update_chain(k)

    def _update_chain(self, k: int) -> None:
        """Todo elemento intacto cortado por la trayectoria puede activarse."""
        s = self._theta[k] >= 0.0
        cut = {p for p, tri in enumerate(self._conn)
               if self._els[p].discontinuity_state is None and 0 < int(s[tri].sum()) < 3}
        self._chain[k] = cut
        self._tips[k] = set(cut)

    def _fixed_values(self, p: int, n: np.ndarray, x_m: np.ndarray, th: np.ndarray) -> np.ndarray:
        return np.where(th >= 0.0, 1.0, -1.0) * np.maximum(np.abs(th), self._delta[p])


# =========================================================================
# Seguimiento por omisión: uno por región conexa de elementos sin referencia
# =========================================================================

_PENDING: "weakref.WeakSet" = weakref.WeakSet()
_DEFAULTS: "weakref.WeakSet" = weakref.WeakSet()


def register_default(element) -> None:
    """Deja el elemento a la espera de su seguimiento por omisión."""
    _PENDING.add(element)


def default_tracker(element) -> GlobalCrackTracking:
    """Seguimiento por omisión del elemento.

    Agrupa los elementos pendientes por nodos compartidos (identidad del
    objeto nodo): cada región conexa recibe un ``GlobalCrackTracking`` con
    los valores por omisión, o se une al de una región ya creada con la que
    comparta nodos.
    """
    pending = sorted(_PENDING, key=lambda e: e.id)
    _PENDING.clear()
    if not any(e is element for e in pending):
        pending.append(element)
    owner: dict = {}
    for tracker in list(_DEFAULTS):
        for e in tracker.elements:
            for nd in e.nodes:
                owner[id(nd)] = tracker
    for e in pending:
        found = []
        for nd in e.nodes:
            tr = owner.get(id(nd))
            if tr is not None and all(tr is not f for f in found):
                found.append(tr)
        if not found:
            tracker = GlobalCrackTracking()
            _DEFAULTS.add(tracker)
        else:
            tracker = found[0]
            for other in found[1:]:
                if other.cracks:
                    continue                  # no se funden regiones ya agrietadas
                for moved in other.elements:
                    tracker.add(moved)
                    moved.crack_tracker = tracker
                    for nd in moved.nodes:
                        owner[id(nd)] = tracker
                _DEFAULTS.discard(other)
        tracker.add(e)
        e.crack_tracker = tracker
        for nd in e.nodes:
            owner[id(nd)] = tracker
    return element.crack_tracker
