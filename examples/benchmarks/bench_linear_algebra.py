"""Benchmark reproducible de la fase algebraica (ADR 0017).

Mide el coste de **resolver** el sistema, que es lo que domina el tiempo de
un análisis una vez que el ensamblaje por lotes (ADR 0014) dejó esa fase en
1–4 µs por elemento. Sobre mallas ``Hex8 n³`` crecientes compara:

- ``SuperLU`` (SciPy, monohilo, reordenamiento COLAMD) — el backend base;
- ``Pardiso`` (Intel MKL vía ``pypardiso``, multihilo, disección anidada) —
  backend opcional, ``pip install solidum-fem[fast]``.

y reporta además el **relleno** (*fill-in*) ``nnz(L+U) / nnz(K)`` de SuperLU,
que es la razón de fondo por la que un solver directo escala mal: crece con
la talla y se lleva por delante el tiempo y la memoria.

Uso::

    python examples/benchmarks/bench_linear_algebra.py [--sizes 10 15 20] [--repeat 3]

Sin ``pypardiso`` instalado el script mide sólo SuperLU y lo dice.
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

try:
    from solidum.math.linalg import PardisoSolver  # type: ignore[attr-defined]
    HAS_PARDISO = True
except ImportError:
    HAS_PARDISO = False


def stiffness(n: int) -> sp.csr_matrix:
    """``K`` de una malla ``Hex8 n³``, regularizada para ser no singular.

    El desplazamiento de cuerpo rígido se elimina con un término diagonal en
    vez de con apoyos: el objetivo es medir álgebra sobre un patrón de
    dispersión realista, no resolver un problema físico concreto.
    """
    domain = mesh_hex8(n, VonMises3D(E=200.0e9, nu=0.3, sigma_y=250.0e6, H=1.0e9))
    assembler = Assembler(domain)
    assembler.assemble_system()
    K = assembler.K_global.tocsr()
    return (K + sp.eye(K.shape[0], format="csr") * 1.0e6).tocsr()


def timed_min(fn, repeat: int) -> float:
    best = float("inf")
    for _ in range(repeat):
        t0 = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t0)
    return best


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--sizes", type=int, nargs="+", default=[10, 15, 20],
                        help="lados de las mallas Hex8 n³ (default 10 15 20)")
    parser.add_argument("--repeat", type=int, default=3,
                        help="repeticiones; se reporta el mínimo")
    args = parser.parse_args()

    if HAS_PARDISO:
        # La primera llamada localiza la DLL de la MKL; no se cronometra.
        PardisoSolver().solve(sp.eye(10, format="csr") * 2.0, np.ones(10))

    rng = np.random.default_rng(0)
    rows = []
    for n in args.sizes:
        K = stiffness(n)
        b = rng.standard_normal(K.shape[0])
        K_csc = K.tocsc()

        t_lu = timed_min(lambda: spla.splu(K_csc).solve(b), args.repeat)
        lu = spla.splu(K_csc)
        fill = (lu.L.nnz + lu.U.nnz) / K.nnz
        x_ref = lu.solve(b)

        if HAS_PARDISO:
            solver = PardisoSolver()
            t_pd = timed_min(lambda: solver.solve(K, b), args.repeat)
            err = np.linalg.norm(solver.solve(K, b) - x_ref) / np.linalg.norm(x_ref)
        else:
            t_pd, err = float("nan"), float("nan")

        rows.append(dict(n=n, dof=K.shape[0], nnz=K.nnz, fill=fill,
                         t_lu=t_lu, t_pd=t_pd, err=err))

    print()
    print(f"Pardiso disponible: {HAS_PARDISO}"
          f"{'' if HAS_PARDISO else '  (pip install solidum-fem[fast])'}")
    print()
    print("| Malla | DOF | nnz(K) | fill-in SuperLU | SuperLU | Pardiso | Aceleración | Δ relativa |")
    print("|---|---:|---:|---:|---:|---:|---:|---:|")
    for r in rows:
        if HAS_PARDISO:
            speed, err = f"×{r['t_lu'] / r['t_pd']:.1f}", f"{r['err']:.1e}"
            t_pd = f"{r['t_pd']:.3f} s"
        else:
            speed, err, t_pd = "—", "—", "—"
        print(f"| Hex8 {r['n']}³ | {r['dof']} | {r['nnz']} | ×{r['fill']:.1f} "
              f"| {r['t_lu']:.3f} s | {t_pd} | {speed} | {err} |")


if __name__ == "__main__":
    main()
