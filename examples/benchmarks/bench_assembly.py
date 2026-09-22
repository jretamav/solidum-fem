"""Benchmark reproducible del ensamblaje (ADR 0014, fase 0).

Mide, por separado, el ensamblaje tangente, el commit del estado y la
factorización del sistema reducido en las dos mallas de referencia de la
propuesta de vectorización por lotes:

- ``100 × 100`` Quad4 con ``VonMises2D`` (10 000 elementos, 40 000 puntos
  de Gauss);
- ``20 × 20 × 20`` Hex8 con ``VonMises3D`` (8 000 elementos, 64 000 puntos).

Cada malla se ensambla por el camino por elemento (``batch=False``) y por
el camino por lotes (``batch=True``), y se comprueba que ``K`` y ``F_int``
coinciden. Las cifras se imprimen como tabla Markdown para pegarlas en
``docs/STATUS.md``.

Uso::

    python examples/benchmarks/bench_assembly.py [--nx 100] [--n3 20] [--repeat 3]

La primera evaluación de cada camino calienta el JIT (compila o carga
del caché en disco) y no se cronometra.
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np

import solidum
from solidum import Assembler, Domain
from solidum.elements.solid_3d import Hex8
from solidum.materials.von_mises_3d import VonMises3D


def mesh_quad4(nx: int, ny: int, material) -> Domain:
    d = Domain()
    nid = 1
    for j in range(ny + 1):
        for i in range(nx + 1):
            # Ligera distorsión para que ningún elemento sea un cuadrado exacto.
            d.add_node(nid, [i / nx + 0.1 * np.sin(j) / nx, j / ny + 0.1 * np.cos(i) / ny])
            nid += 1
    eid = 1
    for j in range(ny):
        for i in range(nx):
            n0 = 1 + i + j * (nx + 1)
            nodes = [d.nodes[n0], d.nodes[n0 + 1], d.nodes[n0 + nx + 2], d.nodes[n0 + nx + 1]]
            d.add_element(solidum.Quad4(eid, nodes, material))
            eid += 1
    d.generate_equation_numbers()
    return d


def mesh_hex8(n: int, material) -> Domain:
    d = Domain()
    nid = 1
    ids = {}
    for k in range(n + 1):
        for j in range(n + 1):
            for i in range(n + 1):
                d.add_node(nid, [i / n, j / n, k / n])
                ids[(i, j, k)] = nid
                nid += 1
    eid = 1
    for k in range(n):
        for j in range(n):
            for i in range(n):
                nodes = [d.nodes[ids[c]] for c in (
                    (i, j, k), (i + 1, j, k), (i + 1, j + 1, k), (i, j + 1, k),
                    (i, j, k + 1), (i + 1, j, k + 1), (i + 1, j + 1, k + 1), (i, j + 1, k + 1))]
                d.add_element(Hex8(eid, nodes, material))
                eid += 1
    d.generate_equation_numbers()
    return d


def displacement_field(domain: Domain, amplitude: float) -> np.ndarray:
    """Campo suave con gradiente no nulo para que el material plastifique."""
    U = np.zeros(domain.total_dofs)
    for node in domain.nodes.values():
        x = node.coordinates
        for k, dof in enumerate(node.dofs.values()):
            U[dof] = amplitude * (x[k % len(x)] ** 2 + 0.3 * x[0] * x[-1])
    return U


def time_call(fn, repeat: int) -> float:
    best = float("inf")
    for _ in range(repeat):
        t0 = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t0)
    return best


def bench(domain: Domain, label: str, repeat: int, amplitude: float) -> dict:
    U = displacement_field(domain, amplitude)
    out = {"label": label, "n_elem": len(domain.elements), "n_dof": domain.total_dofs}
    results = {}
    # 1) Ensamblaje desde el mismo estado committed (el inicial): ambos
    #    caminos escriben sólo el trial, así que K y F_int son comparables.
    for batch in (False, True):
        asm = Assembler(domain, batch=batch)
        asm.assemble_non_linear_system(U)        # calienta JIT y topología
        t_asm = time_call(lambda: asm.assemble_non_linear_system(U), repeat)
        K, F = asm.assemble_non_linear_system(U)
        results[batch] = [t_asm, None, K, F]
        # Al invalidar, los estados por lotes vuelven a diccionarios y el
        # siguiente Assembler parte del mismo estado.
        asm.invalidate()
    # 2) Commit (trial → committed) por ambos caminos.
    for batch in (False, True):
        asm = Assembler(domain, batch=batch)
        asm.assemble_non_linear_system(U)
        results[batch][1] = time_call(asm.commit_all_states, repeat)
        asm.invalidate()
    K0, F0 = results[False][2], results[False][3]
    K1, F1 = results[True][2], results[True][3]
    dK = abs(K1 - K0).max() / abs(K0).max()
    dF = np.abs(F1 - F0).max() / max(np.abs(F0).max(), 1e-300)
    out.update(
        t_elem=results[False][0], t_batch=results[True][0],
        c_elem=results[False][1], c_batch=results[True][1],
        dK=dK, dF=dF,
    )
    return out


def main() -> None:
    # La consola de Windows puede estar en cp1252: la tabla lleva µ y Δ.
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--nx", type=int, default=100, help="Quad4 por lado (default 100)")
    parser.add_argument("--n3", type=int, default=20, help="Hex8 por lado (default 20)")
    parser.add_argument("--repeat", type=int, default=3, help="repeticiones; se reporta el mínimo")
    args = parser.parse_args()

    rows = []
    mat2d = solidum.VonMises2D(E=200.0e9, nu=0.3, sigma_y=250.0e6, H=1.0e9)
    rows.append(bench(mesh_quad4(args.nx, args.nx, mat2d), f"Quad4 {args.nx}×{args.nx} + VonMises2D",
                      args.repeat, amplitude=5.0e-3))
    mat3d = VonMises3D(E=200.0e9, nu=0.3, sigma_y=250.0e6, H=1.0e9)
    rows.append(bench(mesh_hex8(args.n3, mat3d), f"Hex8 {args.n3}³ + VonMises3D",
                      args.repeat, amplitude=5.0e-3))

    print()
    print("| Malla | Elementos | Ensamblaje por elemento | Ensamblaje por lotes | Aceleración | Commit por elemento | Commit por lotes | max Δ K / Δ F_int |")
    print("|---|---:|---:|---:|---:|---:|---:|---:|")
    for r in rows:
        print(
            f"| {r['label']} | {r['n_elem']} | {r['t_elem']*1e3:.0f} ms "
            f"({r['t_elem']/r['n_elem']*1e6:.0f} µs/elem) | {r['t_batch']*1e3:.1f} ms "
            f"({r['t_batch']/r['n_elem']*1e6:.1f} µs/elem) | ×{r['t_elem']/r['t_batch']:.0f} "
            f"| {r['c_elem']*1e3:.0f} ms | {r['c_batch']*1e3:.1f} ms | {r['dK']:.1e} / {r['dF']:.1e} |"
        )


if __name__ == "__main__":
    main()
