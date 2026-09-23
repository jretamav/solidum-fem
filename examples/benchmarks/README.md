# Benchmarks

Scripts de medición reproducible. No son ejemplos de modelización: no
generan figuras ni resultados físicos, sino **cifras de rendimiento** que
se citan en `docs/STATUS.md` y en los ADR que las motivan.

| Script | Qué mide | Referencia |
|---|---|---|
| [`bench_assembly.py`](bench_assembly.py) | Ensamblaje tangente y commit por el camino por elemento, por lotes con kernel serie y por lotes con kernel paralelo (`prange`), con comparación a precisión de máquina entre caminos y bit a bit entre serie y paralelo, en `100×100` Quad4 + J2 y `20³` Hex8 + J2. | ADR 0014 (vectorización por lotes), deuda #19 de `docs/STATUS.md` |
| [`bench_linear_algebra.py`](bench_linear_algebra.py) | Solución del sistema lineal en `Hex8 n³` empotradas (apoyos reales): SuperLU, Pardiso (con su pico de memoria) y solver iterativo CG + AMG con modos de cuerpo rígido (setup, iteraciones, memoria de `K`), más el relleno de SuperLU. Los backends no instalados se reportan como «—». | ADR 0017 (backend multihilo), ADR 0018 (solver iterativo) |

## Cómo ejecutar

```bash
python examples/benchmarks/bench_assembly.py            # mallas de referencia
python examples/benchmarks/bench_assembly.py --nx 50 --n3 10 --repeat 5
python examples/benchmarks/bench_linear_algebra.py --sizes 10 20 30   # fase algebraica
```

Imprime una tabla Markdown lista para pegar. La primera evaluación de cada
camino calienta el JIT (compila o carga el caché de Numba en `__pycache__`)
y no se cronometra; las cifras son el mínimo de `--repeat` repeticiones.

Las cifras dependen de la máquina: en `docs/STATUS.md` se registran junto
con la fecha, la versión de Python y de NumPy/Numba con que se obtuvieron.
La columna paralela depende además del número de hilos, que el script
imprime; para fijarlo, `NUMBA_NUM_THREADS=4 python examples/benchmarks/bench_assembly.py`.
