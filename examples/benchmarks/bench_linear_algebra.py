"""Benchmark reproducible de la fase algebraica (ADR 0017, ADR 0018).

Mide el coste de **resolver** el sistema, que es lo que domina el tiempo de
un análisis una vez que el ensamblaje por lotes (ADR 0014) dejó esa fase en
1–4 µs por elemento. Sobre mallas ``Hex8 n³`` empotradas en ``x = 0``
(sistema reducido, con apoyos reales: nada de regularizar la diagonal)
compara los tres tipos de backend:

- ``SuperLU`` (SciPy, directo monohilo, reordenamiento COLAMD) — el suelo,
  siempre disponible. Se omite por encima de ``--superlu-max`` DOF porque su
  coste crece muy deprisa (110 s a 5·10⁴ DOF).
- ``Pardiso`` (Intel MKL, directo multihilo, disección anidada) —
  ``pip install solidum-fem[fast]``. Reporta además su **pico de memoria**
  (``iparm`` 15–17 de la MKL).
- ``iterativo`` (CG + AMG con modos de cuerpo rígido) —
  ``pip install solidum-fem[iterative]``. Reporta el setup de AMG, las
  iteraciones y la memoria de ``K``, que es lo que el iterativo guarda.

y el relleno ``nnz(L+U)/nnz(K)`` de SuperLU, la razón de fondo por la que un
directo escala mal en tiempo y memoria.

Uso::

    python examples/benchmarks/bench_linear_algebra.py [--sizes 10 20 30 40] [--repeat 1]

Los backends no instalados se reportan como "—".
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, __file__.rsplit("bench_linear_algebra.py", 1)[0])

from bench_assembly import mesh_hex8  # noqa: E402

from solidum import Assembler  # noqa: E402
from solidum.materials.von_mises_3d import VonMises3D  # noqa: E402
from solidum.math.linalg import IterativeSolver, StiffnessProperties  # noqa: E402
from solidum.math.linalg.iterative import HAS_PYAMG  # noqa: E402

try:
    import pypardiso  # type: ignore[import-not-found]

    from solidum.math.linalg import PardisoSolver  # type: ignore[attr-defined]
    HAS_PARDISO = True
except ImportError:
    HAS_PARDISO = False


def reduced_system(n: int):
    """``K`` reducida de ``Hex8 n³`` empotrada en ``x = 0`` y su ensamblador
    (que provee los modos de cuerpo rígido para AMG)."""
    domain = mesh_hex8(n, VonMises3D(E=200.0e9, nu=0.3, sigma_y=250.0e6, H=1.0e9))
    x_min = min(node.coordinates[0] for node in domain.nodes.values())
    for node in domain.nodes.values():
        if node.coordinates[0] - x_min < 0.4 / n:
            for dof in list(node.dofs):
                node.fix_dof(dof, 0.0)
    assembler = Assembler(domain)
    assembler.assemble_system()
    K, _, _, _ = assembler.reduce(assembler.K_global.tocsr(), np.zeros(domain.total_dofs))
    return K.tocsr(), assembler


def timed_min(fn, repeat: int):
    best, out = float("inf"), None
    for _ in range(repeat):
        t0 = time.perf_counter()
        out = fn()
        best = min(best, time.perf_counter() - t0)
    return best, out


def mb(nbytes: float) -> float:
    return nbytes / 2**20


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--sizes", type=int, nargs="+", default=[10, 20, 30],
                        help="lados de las mallas Hex8 n³ (default 10 20 30)")
    parser.add_argument("--repeat", type=int, default=1,
                        help="repeticiones; se reporta el mínimo")
    parser.add_argument("--superlu-max", type=int, default=60_000,
                        help="omitir SuperLU por encima de estos DOF")
    args = parser.parse_args()

    if HAS_PARDISO:
        # La primera llamada localiza la DLL de la MKL; no se cronometra.
        PardisoSolver().solve(sp.eye(10, format="csr") * 2.0, np.ones(10))

    rng = np.random.default_rng(0)
    rows = []
    for n in args.sizes:
        K, assembler = reduced_system(n)
        dof = K.shape[0]
        b = rng.standard_normal(dof)
        mem_K = mb(K.data.nbytes + K.indices.nbytes + K.indptr.nbytes)
        row = dict(n=n, dof=dof, nnz=K.nnz, mem_K=mem_K)

        x_ref = None
        if dof <= args.superlu_max:
            K_csc = K.tocsc()
            row["t_lu"], lu = timed_min(lambda: spla.splu(K_csc), args.repeat)
            row["fill"] = (lu.L.nnz + lu.U.nnz) / K.nnz
            x_ref = lu.solve(b)

        if HAS_PARDISO:
            row["t_pd"], x_pd = timed_min(lambda: PardisoSolver().solve(K, b), args.repeat)
            handle = pypardiso.PyPardisoSolver()
            handle.factorize(K)
            row["mem_pd"] = max(handle.get_iparm(15), handle.get_iparm(16) + handle.get_iparm(17)) / 1024
            x_ref = x_pd if x_ref is None else x_ref

        if HAS_PYAMG:
            props = StiffnessProperties(True, True, dof, near_nullspace=assembler.near_nullspace)
            solver = IterativeSolver(props, preconditioner="amg")
            row["t_setup"], fact = timed_min(lambda: solver.factorize(K), args.repeat)
            row["t_krylov"], x_it = timed_min(lambda: fact.solve(b), args.repeat)
            row["it"] = fact.last_iterations
            if x_ref is not None:
                row["err"] = np.linalg.norm(x_it - x_ref) / np.linalg.norm(x_ref)
        rows.append(row)

    def fmt(v, spec, suffix=""):
        return "—" if v is None else f"{v:{spec}}{suffix}"

    print()
    print(f"Pardiso: {HAS_PARDISO} · pyamg: {HAS_PYAMG}")
    print()
    print("| Malla | DOF | fill-in SuperLU | SuperLU | Pardiso | Pico Pardiso | "
          "Iterativo (setup + CG) | Iter. | Memoria K | Δ relativa |")
    print("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in rows:
        it_time = (f"{r['t_setup']:.2f} + {r['t_krylov']:.2f} = {r['t_setup'] + r['t_krylov']:.2f} s"
                   if "t_setup" in r else "—")
        print(f"| Hex8 {r['n']}³ | {r['dof']} | {fmt(r.get('fill'), '.1f')} "
              f"| {fmt(r.get('t_lu'), '.2f', ' s')} | {fmt(r.get('t_pd'), '.2f', ' s')} "
              f"| {fmt(r.get('mem_pd'), '.0f', ' MB')} | {it_time} | {fmt(r.get('it'), 'd')} "
              f"| {r['mem_K']:.0f} MB | {fmt(r.get('err'), '.1e')} |")


if __name__ == "__main__":
    main()
