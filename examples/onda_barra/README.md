# Onda en una barra: integración explícita e implícita

Ejemplo 8 del [manual de ejemplos](../../manuals/Example_manual.pdf) (capítulo
fuente: [`08_onda_barra.md`](../../manuals/sources/examples/08_onda_barra.md)).

**Qué demuestra.** Una barra de acero empotrada, con una fuerza en escalón en
el extremo libre, integrada con `CentralDifferenceSolver` (masa concentrada),
`NewmarkSolver` y `HHTSolver`, frente a la solución exacta de d'Alembert. Mide
el límite de estabilidad del explícito (Δt = h/c), muestra que ese paso es
exacto en los nodos y que un paso menor introduce dispersión, y que la
disipación de HHT-α apaga la cola de oscilaciones pero no el sobrepaso del
frente.

**Cómo ejecutarlo.**

```bash
python examples/onda_barra/run.py
```

**Qué esperar.** Diferencias centradas con Δt = h/c: error menor que 1e-10;
con Δt = 0.5·h/c, un sobrepaso del 26 % en el esfuerzo del empotramiento; con
Δt = 1.005·h/c, divergencia.

![Esfuerzo en el empotramiento](fig_onda.png)
