# Muro calentado bruscamente: conducción transitoria y el método θ

Ejemplo 9 del [manual de ejemplos](../../manuals/Example_manual.pdf) (capítulo
fuente: [`09_muro_termico.md`](../../manuals/sources/examples/09_muro_termico.md)).

**Qué demuestra.** Un muro de concreto cuya cara interior pasa bruscamente de
20 a 80 °C, resuelto con `ThetaMethodSolver` y `Quad4Thermal`, frente a la
serie de Carslaw y Jaeger. Euler implícito (θ = 1) nunca sale del rango físico
y es de orden 1; Crank-Nicolson (θ = 1/2) es de orden 2 pero, con pasos
grandes, oscila y viola el principio del máximo (104 °C con la cara a 80 °C).
La capacidad consistente con pasos muy pequeños enfría el muro antes de
calentarlo.

**Cómo ejecutarlo.**

```bash
python examples/muro_termico/run.py
```

Los avisos sobre la condición inicial que contradice la temperatura impuesta
son esperables: es el choque térmico del ejemplo.

**Qué esperar.** Órdenes temporales 1.00 y 2.00; error de la malla 5·10⁻⁵;
Crank-Nicolson dentro del rango con Δt ≤ 2.5 min.

![Perfiles, historia y convergencia](fig_muro.png)
