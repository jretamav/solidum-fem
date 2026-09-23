# Cilindro elastoplástico: colapso y bloqueo volumétrico (J2)

Ejemplo 4 del [manual de ejemplos](../../manuals/Example_manual.pdf) (capítulo
fuente: [`04_cilindro_elastoplastico.md`](../../manuals/sources/examples/04_cilindro_elastoplastico.md)).

**Qué demuestra.** El cilindro del ejemplo 3 con plasticidad perfecta de von
Mises (`VonMises2D`, deformación plana), cargado hasta el colapso con
`ArcLengthSolver`. Verifica la presión de colapso analítica,
p_lim = (2/√3) σy ln(b/a), y muestra el bloqueo volumétrico del Quad4 con
integración completa en régimen plástico.

**Cómo ejecutarlo.**

```bash
python examples/cilindro_elastoplastico/run.py
```

Reutiliza `modelo` y `carga_de_presion` de [`../cilindro_lame/run.py`](../cilindro_lame/run.py).
Los avisos "trazado detenido ... tras agotar max_steps" son esperables: en la
meseta de colapso la carga nunca alcanza `max_lambda`.

**Qué esperar.** Tri6 y Quad8 se aplanan a menos de 10⁻³ de p_lim
(192.09 MPa); el Quad4 no forma meseta y su esfuerzo medio sale del intervalo
físico.

![Curvas presión-desplazamiento](fig_curvas.png)
