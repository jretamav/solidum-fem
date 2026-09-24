# Edificio de cinco niveles: análisis modal y espectro de respuesta

<!-- ejemplo: edificio_espectro -->

El análisis sísmico por espectro de respuesta es la herramienta de diseño habitual de la ingeniería de edificios: se calculan los modos de vibrar de la estructura, se lee en el espectro la respuesta máxima de cada uno y se combinan. Este ejemplo lo hace sobre un marco de cinco niveles cuya versión idealizada, el edificio de cortante, tiene solución cerrada, y repite la combinación a mano con los modos analíticos para comprobar cada paso.

**Qué se aprende**

- Modelar un marco con `Frame2DEuler`, con la masa concentrada en las losas.
- Calcular frecuencias y modos con `ModalSolver`.
- Hacer un análisis espectral con `ResponseSpectrumSolver`: masas efectivas, combinación SRSS y CQC.
- Combinar correctamente las respuestas: las derivas de entrepiso se combinan modo a modo.

## Planteamiento

Un marco de concreto de {{r:N}} niveles y un vano, con entrepisos de $h =$ {{r:H_m:.0f}} m y claro de {{r:claro_m:.0f}} m. Las columnas son de {{r:lado_cm:.0f}} × {{r:lado_cm:.0f}} cm, con $E =$ {{r:E_GPa:.0f}} GPa, y cada nivel tiene una masa $m =$ {{r:m_t:.0f}} t.

Se adoptan las hipótesis del **edificio de cortante**: losas infinitamente rígidas y columnas inextensibles. En el modelo, las losas son vigas con una rigidez a flexión {{r:rigidez:.0e}} veces la de una columna, las columnas tienen un área {{r:rigidez:.0e}} veces la real, y toda la masa está en las losas. Así cada columna trabaja biempotrada, y la rigidez lateral de cada entrepiso es

$$
k = 2\,\frac{12\,E\,I}{h^3} = \text{ {{r:k_MN_m:.1f}} MN/m} .
$$

El sismo se representa con un espectro de diseño simplificado de aceleración, para un amortiguamiento del 5 %: una rampa desde la aceleración del terreno $a_0 =$ {{r:a0_g:.2f}} $g$ hasta una meseta de $2.5\,a_0$ entre $T_B =$ {{r:T_B:.1f}} s y $T_C =$ {{r:T_C:.1f}} s, y una caída proporcional a $1/T$ después.

## Solución analítica

El edificio de cortante es una cadena de masas $m$ y resortes $k$ iguales, empotrada en la base. Sus frecuencias y modos tienen expresión cerrada:

$$
\omega_j = 2\sqrt{\frac{k}{m}}\,\sin\frac{(2j-1)\,\pi}{2(2n+1)},
\qquad
\phi_j(i) = \sin\frac{(2j-1)\,i\,\pi}{2n+1},
\qquad i, j = 1, \dots, n .
$$

Con ellos se calculan a mano, sin Solidum, el factor de participación y la masa efectiva de cada modo,

$$
\Gamma_j = \frac{\sum_i m\,\phi_j(i)}{\sum_i m\,\phi_j(i)^2},
\qquad
M_j^* = \frac{\left(\sum_i m\,\phi_j(i)\right)^2}{\sum_i m\,\phi_j(i)^2},
$$

el desplazamiento máximo de cada modo, $u_j(i) = \Gamma_j\,\phi_j(i)\,S_a(T_j)/\omega_j^2$, y el cortante basal modal, $V_j = M_j^*\,S_a(T_j)$. La combinación SRSS toma la raíz de la suma de los cuadrados de las respuestas modales.

## Modelo en Solidum

El marco se construye con dos columnas y una losa por nivel. Las columnas usan un material con densidad nula; la densidad de la losa se elige para que su masa sea la del piso:

{{py:run.py#marco}}

El análisis espectral calcula internamente los modos, los combina y devuelve la respuesta máxima:

{{py:run.py#espectral}}

La dirección del sismo se declara por el nombre del grado de libertad (`ux` en todos los nodos), y el espectro como una función de la frecuencia que devuelve el desplazamiento espectral, $S_d = S_a/\omega^2$. También se puede dar como tabla de periodos y valores, que es lo habitual desde el YAML.

## Resultados

![Izquierda: los tres primeros modos, Solidum (puntos) sobre los analíticos (bandas). Derecha: el espectro de diseño con los periodos de los cinco modos.](../../../examples/edificio_espectro/fig_modos.png){#fig:modos}

[TABLA: Periodos, masas efectivas y aceleración espectral de cada modo.]
| Modo | $T$ exacto [s] | $T$ Solidum [s] | Masa efectiva | Acumulada | $S_a/g$ |
|---|---|---|---|---|---|
| 1 | {{r:modo.1.T_exacto:.4f}} | {{r:modo.1.T_fe:.4f}} | {{r:modo.1.masa_pct:.1f}} % | {{r:modo.1.masa_acum_pct:.1f}} % | {{r:modo.1.Sa_g:.3f}} |
| 2 | {{r:modo.2.T_exacto:.4f}} | {{r:modo.2.T_fe:.4f}} | {{r:modo.2.masa_pct:.1f}} % | {{r:modo.2.masa_acum_pct:.1f}} % | {{r:modo.2.Sa_g:.3f}} |
| 3 | {{r:modo.3.T_exacto:.4f}} | {{r:modo.3.T_fe:.4f}} | {{r:modo.3.masa_pct:.1f}} % | {{r:modo.3.masa_acum_pct:.1f}} % | {{r:modo.3.Sa_g:.3f}} |
| 4 | {{r:modo.4.T_exacto:.4f}} | {{r:modo.4.T_fe:.4f}} | {{r:modo.4.masa_pct:.1f}} % | {{r:modo.4.masa_acum_pct:.1f}} % | {{r:modo.4.Sa_g:.3f}} |
| 5 | {{r:modo.5.T_exacto:.4f}} | {{r:modo.5.T_fe:.4f}} | {{r:modo.5.masa_pct:.1f}} % | {{r:modo.5.masa_acum_pct:.1f}} % | {{r:modo.5.Sa_g:.3f}} |

**Modos.** Los periodos coinciden con los de la solución cerrada con una diferencia relativa de {{r:error_w_max:.1e}}, y las formas modales con {{r:error_formas:.0e}}. Esa pequeña diferencia es la rigidez finita de losas y columnas: con las losas diez veces más rígidas baja a {{r:error_w_max_rigido:.1e}}. El edificio de cortante es el límite del marco cuando las losas se vuelven rígidas. El primer modo concentra el {{r:modo.1.masa_pct:.0f}} % de la masa; con tres modos se supera el 99 %. Los reglamentos piden incluir modos hasta reunir al menos el 90 %.

[TABLA: Desplazamiento máximo de cada nivel con SRSS y deriva de entrepiso combinada de dos maneras.]
| Nivel | $u$ Solidum [mm] | $u$ a mano [mm] | Deriva, modo a modo [mm] | Deriva, restando desplazamientos [mm] |
|---|---|---|---|---|
| 1 | {{r:piso.1.u_mm:.2f}} | {{r:piso.1.u_exacto_mm:.2f}} | {{r:piso.1.deriva_bien_mm:.2f}} | {{r:piso.1.deriva_mal_mm:.2f}} |
| 2 | {{r:piso.2.u_mm:.2f}} | {{r:piso.2.u_exacto_mm:.2f}} | {{r:piso.2.deriva_bien_mm:.2f}} | {{r:piso.2.deriva_mal_mm:.2f}} |
| 3 | {{r:piso.3.u_mm:.2f}} | {{r:piso.3.u_exacto_mm:.2f}} | {{r:piso.3.deriva_bien_mm:.2f}} | {{r:piso.3.deriva_mal_mm:.2f}} |
| 4 | {{r:piso.4.u_mm:.2f}} | {{r:piso.4.u_exacto_mm:.2f}} | {{r:piso.4.deriva_bien_mm:.2f}} | {{r:piso.4.deriva_mal_mm:.2f}} |
| 5 | {{r:piso.5.u_mm:.2f}} | {{r:piso.5.u_exacto_mm:.2f}} | {{r:piso.5.deriva_bien_mm:.2f}} | {{r:piso.5.deriva_mal_mm:.2f}} |

**Respuesta espectral.** Los desplazamientos SRSS de Solidum coinciden con los combinados a mano con los modos analíticos a {{r:error_u_srss:.1e}}. El cortante basal, combinando los cortantes modales, vale {{r:V_fe_kN:.1f}} kN ({{r:V_srss_kN:.1f}} kN a mano): {{r:V_srss_sobre_W:.3f}} veces el peso del edificio.

**SRSS y CQC.** La combinación cuadrática completa (CQC) tiene en cuenta la correlación entre modos de frecuencias cercanas. Aquí los modos están bien separados y los desplazamientos CQC quedan entre {{r:cqc_sobre_srss_min:.4f}} y {{r:cqc_sobre_srss_max:.4f}} veces los SRSS. La diferencia importa cuando hay modos próximos, por ejemplo en edificios con torsión.

```callout Advertencia
Las respuestas máximas de los modos no ocurren a la vez y la combinación pierde su signo, así que **cada magnitud de diseño se combina por separado**. La deriva de un entrepiso se calcula en cada modo y después se combina; restar los desplazamientos ya combinados de dos niveles consecutivos da otro valor. En este edificio subestima la deriva del último entrepiso: {{r:piso.5.deriva_mal_mm:.2f}} mm en lugar de {{r:piso.5.deriva_bien_mm:.2f}} mm. Lo mismo vale para los esfuerzos y las fuerzas internas de los elementos.
```

## Para seguir

- Dar a las losas la rigidez de una viga real, por ejemplo de 30 × 60 cm. Los nodos giran, las columnas dejan de estar biempotradas y los periodos crecen: el marco deja de ser un edificio de cortante.
- Calcular sólo dos modos (`n_modes=2`). ¿Qué fracción del cortante basal se pierde? ¿Y de la deriva del último entrepiso?
- Declarar el mismo espectro en un YAML como tabla de periodos y valores, y resolverlo con `solidum.run_yaml`.

**Referencias.** A. K. Chopra, *Dynamics of Structures*, 4.ª ed., Pearson, 2012. E. L. Wilson, A. Der Kiureghian y E. P. Bayo, "A replacement for the SRSS method in seismic analysis", *Earthquake Engineering and Structural Dynamics* 9 (1981) 187-194.

**Archivos del ejemplo.** La carpeta `examples/edificio_espectro/` contiene el script `run.py`, los resultados en `resultados.json` y la figura.
