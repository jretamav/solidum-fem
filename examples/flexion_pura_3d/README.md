# Flexión pura de un prisma: el bloqueo por cortante del Hex8

Ejemplo 6 del [manual de ejemplos](../../manuals/Example_manual.pdf) (capítulo
fuente: [`06_flexion_pura_3d.md`](../../manuals/sources/examples/06_flexion_pura_3d.md)).

**Qué demuestra.** Un prisma de sección cuadrada en flexión pura, con la
solución exacta de Saint-Venant (curvatura anticlástica incluida), resuelto con
`Hex8` y `Hex20`. La curvatura se impone con el desplazamiento axial de los
extremos y se mide el momento de reacción. El Hex20 es exacto con un solo
elemento; el Hex8 se bloquea por cortante: aparece un esfuerzo cortante
parásito G·κ·a/(2√3) en sus puntos de Gauss y necesita hasta 2.6 veces el
momento exacto para la misma curvatura.

**Cómo ejecutarlo.**

```bash
python examples/flexion_pura_3d/run.py
```

**Qué esperar.** Hex20 con error menor que 1e-10; Hex8 entre 2.58 (5×2×2) y
1.016 (80×4×4) veces el momento exacto, con el cortante parásito sobre la
recta teórica.

![Bloqueo por cortante](fig_bloqueo.png)
