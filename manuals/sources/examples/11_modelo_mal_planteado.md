# Un modelo mal planteado: la red de seguridad del análisis estático

<!-- ejemplo: modelo_mal_planteado -->

A un modelo le puede faltar un apoyo o una barra. Entonces la estructura es un mecanismo: su solución no es única, o para ciertas cargas no existe. El solver algebraico no se da cuenta. Ante una matriz de rigidez singular, un solver directo devuelve números sin avisar: a veces desplazamientos absurdos, de miles de kilómetros o más, y a veces valores que parecen razonables. En el análisis estático lineal, Solidum no deja que esos números lleguen al usuario. Comprueba el modelo antes y después de resolver, y si es un mecanismo lo rechaza con un mensaje que explica qué pasa. Este ejemplo comete a propósito los errores de apoyo y de conectividad más habituales en una armadura, muestra qué responde Solidum y lo contrasta con la estática, que llega a la misma conclusión sin usar ninguna rigidez.

**Qué se aprende**

- Distinguir un mecanismo rígido (apoyos insuficientes o mal orientados) de uno interno (una barra de menos).
- Leer los mensajes de `MechanismError` e `IllPosedSystemError` y corregir el modelo.
- Encontrar un mecanismo interno con un análisis modal.
- Por qué contar barras y apoyos no basta, y qué no detectan estas comprobaciones.

## Planteamiento

Una armadura Pratt de cuatro paneles de {{r:panel_m:.0f}} m y {{r:alto_m:.0f}} m de altura, con {{r:n_nudos}} nudos y {{r:n_barras}} barras de acero ($E =$ {{r:E_GPa:.0f}} GPa, $A =$ {{r:A_cm2:.0f}} cm²) modeladas con elementos `Truss2D`. Soporta $P =$ {{r:P_kN:.0f}} kN en cada uno de los tres nudos inferiores interiores. Bien apoyada, tiene un pasador en el extremo izquierdo y un rodillo en el derecho que impide el desplazamiento vertical.

Sobre ella se prueban cuatro variantes erróneas:

1. **Dos rodillos**: los dos apoyos sólo impiden el desplazamiento vertical.
2. **Rodillo horizontal**: el rodillo impide el desplazamiento horizontal en vez del vertical.
3. **Sin diagonal**: falta la diagonal del segundo panel.
4. **Sin diagonal, con carga horizontal**: el mismo modelo, cargado sólo con {{r:P_kN:.0f}} kN horizontales en el rodillo.

## Solución analítica

Cada nudo aporta dos ecuaciones de equilibrio: $2j =$ {{r:casos.correcto.dos_j}} en total. Las incógnitas son las $m$ fuerzas axiales y las $r$ reacciones. Escritas juntas, forman el **método de los nudos** en forma matricial,

$$
\mathbf A\,\mathbf s = -\mathbf F ,
$$

con $\mathbf s$ las fuerzas axiales (tracción positiva) y las reacciones, y $\mathbf A$ la matriz de equilibrio, que sólo depende de la geometría. Su traspuesta $\mathbf A^\top$ es, salvo el signo, la matriz de compatibilidad: relaciona los desplazamientos de los nudos con el alargamiento de las barras y el movimiento de los apoyos. De ahí salen dos conclusiones:

- Si el rango de $\mathbf A$ es menor que $2j$, hay $2j - \operatorname{rango}(\mathbf A)$ **mecanismos**: movimientos que no alargan ninguna barra ni mueven ningún apoyo. Son el núcleo de $\mathbf A^\top$. Para ciertas cargas el equilibrio es imposible.
- Una condición necesaria para evitarlo es la **cuenta de Maxwell**: $m + r \ge 2j$. Pero no es suficiente: las barras o los apoyos pueden estar mal colocados.

La armadura bien apoyada tiene $m + r = 2j$ y $\mathbf A$ de rango completo: es **isostática**. Sus fuerzas axiales salen sólo del equilibrio, sin intervenir $E$ ni $A$. Su flecha en el centro sale del teorema de los trabajos virtuales, $\delta = \sum N\,n\,L/(EA)$, con $n$ las fuerzas axiales que produce una carga unidad en ese punto.

## Modelo en Solidum

La armadura de cada caso:

{{py:run.py#armadura}}

El análisis, con los dos errores que Solidum puede lanzar:

{{py:run.py#analizar}}

`LinearSolver` comprueba el modelo en dos momentos:

- **Antes de resolver**: qué movimientos de sólido rígido dejan libres los apoyos. No factoriza nada. Si queda alguno, lanza `MechanismError` y lo describe en términos del modelo: una traslación, o un giro alrededor de un punto.
- **Después de factorizar**: si la factorización encontró pivotes numéricamente nulos, la matriz es singular. Además, la solución tiene que cumplir el equilibrio, $\lVert\mathbf F - \mathbf K\mathbf u\rVert \le 10^{-8}\,\lVert\mathbf F\rVert$. Si falla cualquiera de las dos cosas, lanza `IllPosedSystemError`. Es lo que atrapa los mecanismos internos, que los apoyos no pueden delatar.

El usuario no elige el solver algebraico: por omisión Solidum usa un solver directo, SuperLU o Pardiso según lo que esté instalado, y con los dos el diagnóstico es el mismo en los cinco casos.

## Resultados

![(a) Armadura bien apoyada: el grosor de cada barra crece con su fuerza axial; en gris, las barras sin fuerza. (b) a (d): los mecanismos de los tres modelos mal apoyados o mal conectados, calculados con la estática (el núcleo de la matriz de compatibilidad); en gris la geometría inicial y en color el mecanismo, amplificado.](../../../examples/modelo_mal_planteado/fig_armadura.png){#fig:armadura}

[TABLA: Los cinco casos. $m + r$ frente a $2j =$ {{r:casos.correcto.dos_j}} es la cuenta de Maxwell; los mecanismos y la existencia de equilibrio salen del rango de la matriz de equilibrio de los nudos.]
| Caso | $m + r$ / mecanismos | ¿Equilibrio? | Respuesta de Solidum |
|---|---|---|---|
| Bien apoyada | {{r:casos.correcto.m_mas_r}} / {{r:casos.correcto.mecanismos}} | {{r:casos.correcto.equilibrio_txt}} | solución |
| Dos rodillos | {{r:casos.dos_rodillos.m_mas_r}} / {{r:casos.dos_rodillos.mecanismos}} | {{r:casos.dos_rodillos.equilibrio_txt}} | `{{r:casos.dos_rodillos.resultado}}` |
| Rodillo horizontal | {{r:casos.rodillo_girado.m_mas_r}} / {{r:casos.rodillo_girado.mecanismos}} | {{r:casos.rodillo_girado.equilibrio_txt}} | `{{r:casos.rodillo_girado.resultado}}` |
| Sin diagonal | {{r:casos.sin_diagonal.m_mas_r}} / {{r:casos.sin_diagonal.mecanismos}} | {{r:casos.sin_diagonal.equilibrio_txt}} | `{{r:casos.sin_diagonal.resultado}}` |
| Sin diagonal, carga horizontal | {{r:casos.sin_diagonal_H.m_mas_r}} / {{r:casos.sin_diagonal_H.mecanismos}} | {{r:casos.sin_diagonal_H.equilibrio_txt}} | `{{r:casos.sin_diagonal_H.resultado}}` |

**Bien apoyada.** Las {{r:n_barras}} fuerzas axiales coinciden con el método de los nudos con una diferencia relativa menor que $10^{-10}$, y la flecha en el centro, {{r:correcto.flecha_mm:.2f}} mm, con los trabajos virtuales (@fig:armadura a). La mayor tracción, {{r:correcto.N_max_kN:.1f}} kN, la lleva la diagonal de cada panel extremo, que transmite todo el cortante del panel, la reacción $1.5P$: $1.5\sqrt 2\,P$. La mayor compresión, {{r:correcto.N_min_kN:.0f}} kN, la lleva la cuerda superior de los paneles centrales: el momento flector en el centro dividido entre la altura, $-M/h = -2P$.

**Dos rodillos.** Faltan reacciones: $m + r =$ {{r:casos.dos_rodillos.m_mas_r}}, menos que $2j$, y la cuenta ya lo delata. Nada impide que la armadura se deslice en horizontal (@fig:armadura b). Con estas cargas, todas verticales, la estática sí tiene solución, porque nada empuja la armadura hacia un lado. Pero el viento o una imperfección del montaje la desplazarían sin límite. Solidum la rechaza antes de resolver. Su diagnóstico: *{{r:casos.dos_rodillos.movimiento}}*.

**Rodillo horizontal.** Aquí la cuenta se cumple, $m + r = 2j$, y aun así hay un mecanismo. Las tres reacciones pasan por el pasador: las dos del pasador, y la del rodillo, que es horizontal y actúa a lo largo de la cuerda inferior. Tres fuerzas que pasan por un mismo punto no pueden equilibrar un momento respecto de él, y la armadura gira alrededor del pasador (@fig:armadura c). Con las cargas verticales ni siquiera hay equilibrio. Solidum encuentra el giro y dice dónde está su centro:

```
{{r:casos.rodillo_girado.mensaje}}
```

**Sin diagonal.** Los apoyos son correctos y la armadura no puede moverse como sólido rígido, así que la comprobación previa pasa. Pero sin su diagonal el segundo panel es un cuadrilátero articulado, que se distorsiona sin alargar ninguna barra. Las partes de la armadura a cada lado del panel giran lo mismo y se desplazan en vertical una respecto de la otra (@fig:armadura d). Las cuerdas del panel siguen resistiendo el momento flector; lo que el panel no puede transmitir es el cortante, y con estas cargas lo tiene. Es un mecanismo interno. Solidum lo detecta al factorizar, con {{r:casos.sin_diagonal.pivotes_nulos}} pivote nulo, y el mensaje apunta a las causas habituales:

```
{{r:casos.sin_diagonal.mensaje}}
```

**Dónde está el mecanismo.** El mensaje dice qué pasa, pero no dónde. Para encontrar un mecanismo interno basta un análisis modal del mismo modelo, con densidad en el material. Cada mecanismo aparece como una frecuencia nula, porque se mueve sin deformar nada, y su forma modal muestra qué parte se mueve:

{{py:run.py#buscar_mecanismo}}

Aquí `ModalSolver` da {{r:modal.frecuencias_nulas}} frecuencia nula, seguida de una de {{r:modal.f1_hz:.1f}} Hz, y su forma modal coincide con el mecanismo de la estática (@fig:armadura d).

**Sin diagonal, con carga horizontal.** Esta carga no activa el mecanismo: la cuerda inferior la lleva hasta el pasador y la estática tiene solución. El residuo del equilibrio no delataría nada, porque existe un campo de desplazamientos que lo cumple. El pivote nulo sí lo delata, y Solidum rechaza el modelo igual que antes. Tiene razón: esa armadura sólo resiste las cargas que no activan el mecanismo, y cualquier carga vertical en un nudo del vano la derrumbaría.

```callout Advertencia
Estas comprobaciones protegen de los modelos mal planteados, no de los bien planteados con datos equivocados. Si el área se escribe en cm² en unas unidades que la esperan en m², $A = 20$ en vez de $20\cdot 10^{-4}$, la armadura es $10^4$ veces más rígida. Las fuerzas axiales no cambian, porque en una armadura isostática sólo dependen del equilibrio, pero la flecha pasa de {{r:correcto.flecha_mm:.2f}} mm a {{r:unidades.flecha_um:.2f}} µm. El modelo está bien planteado y Solidum no avisa. Comprobar el orden de magnitud de los desplazamientos sigue siendo trabajo del usuario.
```

```callout Advertencia
La comprobación después de factorizar sólo se hace en el análisis estático lineal con un solver directo. Con `linear_algebra: iterative`, o con `NonlinearSolver`, el último caso, el del mecanismo que la carga no activa, pasa sin aviso: el resultado cumple el equilibrio, pero sus desplazamientos incluyen una parte arbitraria del mecanismo. En un análisis no lineal no puede ser de otro modo, porque cerca de un punto límite la matriz tangente es casi singular por naturaleza y rechazarla impediría seguir la curva. Ante la duda, resolver primero el modelo con `LinearSolver`.
```

## Para seguir

- Quitar todos los apoyos y calcular los modos de vibración. ¿Cuántas frecuencias nulas aparecen, y a qué movimientos corresponden? En un análisis modal, un modelo libre es legítimo.
- Poner el rodillo horizontal en el nudo superior derecho en vez del inferior. ¿Sigue siendo un mecanismo?
- Añadir en un panel una segunda diagonal, cruzada con la primera. Ahora $m + r > 2j$. ¿Dependen las fuerzas axiales del área de esa barra?

**Referencias.** J. C. Maxwell, "On the calculation of the equilibrium and stiffness of frames", *Philosophical Magazine* 27 (1864) 294-299. S. Pellegrino y C. R. Calladine, "Matrix analysis of statically and kinematically indeterminate frameworks", *International Journal of Solids and Structures* 22 (1986) 409-428.

**Archivos del ejemplo.** La carpeta `examples/modelo_mal_planteado/` contiene el script `run.py`, los resultados en `resultados.json` y la figura.
