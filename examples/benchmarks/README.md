# Benchmarks

Scripts de medición reproducible. No son ejemplos de modelización: no
generan figuras ni resultados físicos, sino **cifras de rendimiento** que
se citan en `docs/STATUS.md` y en los ADR que las motivan.

| Script | Qué mide | Referencia |
|---|---|---|
| [`bench_assembly.py`](bench_assembly.py) | Ensamblaje tangente, commit y comparación bit a bit entre el camino por elemento y el camino por lotes, en `100×100` Quad4 + J2 y `20³` Hex8 + J2. | ADR 0014 (vectorización por lotes), deuda #19 de `docs/STATUS.md` |

## Cómo ejecutar

```bash
python examples/benchmarks/bench_assembly.py            # mallas de referencia
python examples/benchmarks/bench_assembly.py --nx 50 --n3 10 --repeat 5
```

Imprime una tabla Markdown lista para pegar. La primera evaluación de cada
camino calienta el JIT (compila o carga el caché de Numba en `__pycache__`)
y no se cronometra; las cifras son el mínimo de `--repeat` repeticiones.

Las cifras dependen de la máquina: en `docs/STATUS.md` se registran junto
con la fecha, la versión de Python y de NumPy/Numba con que se obtuvieron.
