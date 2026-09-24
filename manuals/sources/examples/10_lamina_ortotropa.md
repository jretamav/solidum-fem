# Lámina ortótropa fuera de ejes: acoplamiento tracción-cortante

<!-- ejemplo: lamina_ortotropa -->

Las fibras de un material compuesto dan mucha rigidez en su dirección y poca en las demás. Cuando la carga no sigue las fibras pasan dos cosas: la rigidez cae mucho más de lo que sugiere el ángulo, y una tracción pura distorsiona la pieza, un acoplamiento que ningún material isótropo tiene. Este ejemplo mide los dos efectos en una lámina de carbono-epoxi y los compara con la teoría clásica de láminas. Después muestra una consecuencia práctica: en el ensayo de tracción, las mordazas impiden esa distorsión y alteran el módulo que se mide.

**Qué se aprende**

- Definir un material ortótropo con `Orthotropic2D` y orientar sus fibras con `theta`.
- Calcular las constantes aparentes de una lámina con las fibras a un ángulo cualquiera.
- Reconocer el acoplamiento tracción-cortante.
- Por qué la probeta de un ensayo fuera de ejes debe ser larga.

## Planteamiento

Una lámina unidireccional de carbono-epoxi tiene $E_1 =$ {{r:E1_GPa:.0f}} GPa en la dirección de las fibras (eje 1), $E_2 =$ {{r:E2_GPa:.1f}} GPa en la transversal (eje 2), $G_{12} =$ {{r:G12_GPa:.2f}} GPa y $\nu_{12} =$ {{r:nu12:.2f}}, valores típicos de un T300/5208. Es delgada, así que trabaja en esfuerzo plano. De ella se corta una tira de ancho $w =$ {{r:W_mm:.0f}} mm con las fibras a un ángulo $\theta$ del eje $x$, que es el eje de carga.

El caso se resuelve dos veces:

1. **Tira libre.** Se tracciona con $\sigma_x =$ {{r:sigma_MPa:.0f}} MPa aplicada como tracción en el borde derecho. El borde izquierdo sólo tiene impedido $u_x$ (y $u_y$ en una esquina), así que la tira puede estrecharse y distorsionarse sin restricción.
2. **Tira con mordazas.** Los dos extremos están sujetos: $u_y = 0$, y $u_x$ uniforme, nulo en un extremo y $\delta$ en el otro. Es lo que hacen las mordazas rígidas de una máquina de ensayos.

## Solución analítica

En los ejes del material la flexibilidad de esfuerzo plano es diagonal por bloques, sin acoplamiento entre alargamientos y distorsión:

$$
\mathbf S = \begin{bmatrix} 1/E_1 & -\nu_{12}/E_1 & 0 \\ -\nu_{12}/E_1 & 1/E_2 & 0 \\ 0 & 0 & 1/G_{12} \end{bmatrix},
\qquad
\frac{\nu_{21}}{E_2} = \frac{\nu_{12}}{E_1} .
$$

La segunda igualdad es la reciprocidad de Maxwell-Betti ($\mathbf S$ es simétrica). Da $\nu_{21} = \nu_{12}E_2/E_1 =$ {{r:nu21:.4f}}: con la carga transversal a las fibras, la lámina casi no se estrecha.

En los ejes de carga, la flexibilidad es $\bar{\mathbf S} = \mathbf T^{-1}\,\mathbf S\,\mathbf T^{-\top}$, con $\mathbf T$ la matriz que pasa las deformaciones de ingeniería de los ejes $x, y$ a los ejes $1, 2$. Una tracción $\sigma_x$ sola produce $\varepsilon_x = \bar S_{11}\sigma_x$, $\varepsilon_y = \bar S_{12}\sigma_x$ y $\gamma_{xy} = \bar S_{16}\sigma_x$. De ahí salen las tres constantes aparentes de la lámina:

$$
E_x = \frac{1}{\bar S_{11}}, \qquad
\nu_{xy} = -\frac{\bar S_{12}}{\bar S_{11}}, \qquad
\eta_{xy,x} = \frac{\gamma_{xy}}{\varepsilon_x} = \frac{\bar S_{16}}{\bar S_{11}} .
$$

El **coeficiente de influencia mutua** $\eta_{xy,x}$ mide el acoplamiento: la distorsión que produce una tracción, por unidad de alargamiento. Es nulo con las fibras a 0° o a 90° y en cualquier material isótropo. Con la tira libre el estado es uniforme, así que el cálculo debe reproducir estas fórmulas sin error de discretización.

## Modelo en Solidum

La tira, con el material orientado por `theta` en grados:

{{py:run.py#tira}}

La tracción libre. `compute_edge_traction` convierte la tracción del borde en fuerzas nodales equivalentes:

{{py:run.py#traccion_libre}}

Las mordazas. Se impone el desplazamiento $\delta$ y la fuerza $P$ es la suma de las fuerzas internas en los nodos del extremo derecho. El módulo aparente es el que calcularía el laboratorio, $E_{\mathrm{ap}} = (P/w)/(\delta/L)$:

{{py:run.py#traccion_mordazas}}

La tira libre usa 4 × 2 elementos `Quad8`, porque el estado es uniforme. La tira con mordazas usa 16 elementos por ancho y se repite con 8 para comprobar la convergencia.

## Resultados

![Izquierda: módulo aparente y coeficiente de influencia mutua frente al ángulo de las fibras; curvas de la teoría de láminas y puntos de Solidum. Derecha: tira de longitud $4w$ con las fibras a 45°, libre y con mordazas; las dos deformadas amplificadas hasta el mismo alargamiento.](../../../examples/lamina_ortotropa/fig_lamina.png){#fig:lamina}

[TABLA: Constantes aparentes de la tira libre calculadas con Solidum. Todas coinciden con la teoría de láminas con una diferencia menor que $10^{-9}$.]
| $\theta$ | $E_x$ [GPa] | $\nu_{xy}$ | $\eta_{xy,x}$ |
|---|---|---|---|
| 0° | {{r:libre.0.E_GPa:.1f}} | {{r:libre.0.nu:.3f}} | 0 |
| 10° | {{r:libre.10.E_GPa:.1f}} | {{r:libre.10.nu:.3f}} | {{r:libre.10.eta:.3f}} |
| 20° | {{r:libre.20.E_GPa:.1f}} | {{r:libre.20.nu:.3f}} | {{r:libre.20.eta:.3f}} |
| 30° | {{r:libre.30.E_GPa:.1f}} | {{r:libre.30.nu:.3f}} | {{r:libre.30.eta:.3f}} |
| 45° | {{r:libre.45.E_GPa:.1f}} | {{r:libre.45.nu:.3f}} | {{r:libre.45.eta:.3f}} |
| 60° | {{r:libre.60.E_GPa:.1f}} | {{r:libre.60.nu:.3f}} | {{r:libre.60.eta:.3f}} |
| 75° | {{r:libre.75.E_GPa:.1f}} | {{r:libre.75.nu:.3f}} | {{r:libre.75.eta:.3f}} |
| 90° | {{r:libre.90.E_GPa:.1f}} | {{r:libre.90.nu:.3f}} | 0 |

**La tira libre es exacta.** Con el estado uniforme, los `Quad8` reproducen la teoría de láminas a precisión de máquina en los ocho ángulos (@fig:lamina, izquierda). A 90° el Poisson aparente es $\nu_{21}$, la reciprocidad en acción.

**La rigidez cae deprisa.** Basta girar las fibras 10° para que $E_x$ baje de {{r:libre.0.E_GPa:.0f}} a {{r:libre.10.E_GPa:.0f}} GPa. A 45° queda en {{r:libre.45.E_GPa:.1f}} GPa, el {{r:E45_pct:.1f}} % de $E_1$, cuando el promedio de $E_1$ y $E_2$ sugeriría cerca de la mitad. La causa es el cortante: fuera de ejes, buena parte del alargamiento viene de distorsionar la matriz de resina, y $G_{12}$ es el menor de los tres módulos.

**Una tracción produce distorsión.** A 45°, cada unidad de alargamiento viene acompañada de {{r:libre.45.eta:.2f}} unidades de distorsión angular, y a 10° el acoplamiento pasa de 2. El signo tiene una lectura sencilla. Con las fibras a +45°, la tira se alarga menos en la dirección de las fibras, que es rígida, que en la perpendicular; esa diferencia es precisamente $\gamma_{xy}$, que resulta negativa. La tira se convierte en un paralelogramo (@fig:lamina, derecha).

```callout Advertencia
En Solidum, `nu12` es $-\varepsilon_2/\varepsilon_1$ con la carga en la dirección 1, la de las fibras. Algunas fuentes llaman $\nu_{12}$ al cociente contrario, que aquí vale {{r:nu21:.3f}}. Si se copia ese valor, el cálculo no falla y el material resultante es otro: su contracción transversal es $E_1/E_2$ veces menor que la real.
```

**Las mordazas impiden la distorsión.** Con los extremos sujetos, la tira no puede convertirse en paralelogramo. Junto a las mordazas no puede distorsionarse; lejos de ellas sí, y para hacer compatibles las dos zonas se flexiona en su plano: se deforma en S (@fig:lamina, derecha). Impedir una deformación rigidiza, así que el módulo aparente crece.

[TABLA: Módulo aparente con mordazas, relativo al de la tira libre $E_x$. La última fila es la cota superior $\bar Q_{11}/E_x$.]
| $L/w$ | $\theta = 30°$ | $\theta = 45°$ |
|---|---|---|
| 2 | {{r:mordazas.30.2.E_sobre_Ex:.2f}} | {{r:mordazas.45.2.E_sobre_Ex:.2f}} |
| 4 | {{r:mordazas.30.4.E_sobre_Ex:.2f}} | {{r:mordazas.45.4.E_sobre_Ex:.2f}} |
| 8 | {{r:mordazas.30.8.E_sobre_Ex:.2f}} | {{r:mordazas.45.8.E_sobre_Ex:.2f}} |
| Cota superior | {{r:cota.30:.2f}} | {{r:cota.45:.2f}} |

Este resultado es numérico, sin solución cerrada. Pasar de 8 a 16 elementos por ancho cambia el módulo aparente como mucho en {{r:mordazas.30.2.cambio_malla:.1e}} (30°, $L/w = 2$), así que las cifras de la tabla son fiables al 1 %. Además, los teoremas de energía lo acotan por los dos lados:

- **Por abajo, $E_x$.** El esfuerzo uniforme de la tira libre está en equilibrio y es compatible con las mordazas. El teorema de la mínima energía complementaria da $E_{\mathrm{ap}} \ge E_x$.
- **Por arriba, $\bar Q_{11}$.** El desplazamiento uniforme $\varepsilon_x = \delta/L$, $\varepsilon_y = \gamma_{xy} = 0$ cumple las condiciones de las mordazas. El teorema de la mínima energía potencial da $E_{\mathrm{ap}} \le \bar Q_{11}$, con $\bar{\mathbf Q} = \bar{\mathbf S}^{-1}$ la rigidez de la lámina en ejes de carga. Es el módulo de una tira infinitamente corta, sujeta en toda su longitud.

Los seis casos caen dentro de las cotas. Cuanto más corta es la tira, más rígida resulta: con $L/w = 2$ el módulo medido a 45° excede a $E_x$ en un {{r:mordazas.45.2.exceso_pct:.0f}} %, y a 30° en un {{r:mordazas.30.2.exceso_pct:.0f}} %. A 30° el efecto es mayor porque el acoplamiento es mayor, $\eta_{xy,x} =$ {{r:libre.30.eta:.2f}} frente a {{r:libre.45.eta:.2f}}: hay más distorsión que impedir.

```callout Advertencia
Un ensayo de tracción fuera de ejes con una probeta corta sobreestima el módulo y no avisa: la curva fuerza-desplazamiento sale perfectamente lineal. Aun con $L/w = 8$, el error es del {{r:mordazas.45.8.exceso_pct:.1f}} % a 45° y del {{r:mordazas.30.8.exceso_pct:.1f}} % a 30°. Por eso las probetas fuera de ejes se hacen largas y esbeltas.
```

## Para seguir

- Repetir la tira con mordazas a 10°, cerca del ángulo de mayor acoplamiento. ¿Cuánto crece el módulo aparente con $L/w = 8$?
- Medir el alargamiento con un extensómetro virtual: la diferencia de $u_x$ entre dos puntos del centro de la tira, en vez de $\delta/L$. ¿Se recupera $E_x$?
- Cambiar el material por uno de vidrio-epoxi, con $E_1/E_2$ mucho menor. ¿Qué pasa con $\eta_{xy,x}$?

**Referencias.** R. M. Jones, *Mechanics of Composite Materials*, 2.ª ed., Taylor & Francis, 1999, cap. 2. N. J. Pagano y J. C. Halpin, "Influence of end constraint in the testing of anisotropic bodies", *Journal of Composite Materials* 2 (1968) 18-31.

**Archivos del ejemplo.** La carpeta `examples/lamina_ortotropa/` contiene el script `run.py`, los resultados en `resultados.json` y la figura.
