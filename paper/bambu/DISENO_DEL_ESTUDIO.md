# Diseño del estudio numérico — validez de la idealización isótropa en bambú

> **Fase 0 del plan.** Este documento fija las decisiones metodológicas **antes** de
> calcular nada. Su función es que el estudio sea reproducible y que las métricas no
> se elijan a posteriori para favorecer una conclusión.
>
> Fecha: 2026-09-10. Estado: **diseño cerrado, sin resultados**.

---

## 1. Pregunta de investigación

> **¿Bajo qué condiciones es admisible modelar el bambú como material isótropo, y qué
> error se comete al hacerlo?**

La pregunta no admite una respuesta escalar. Se descompone en tres dependencias que el
estudio trata por separado:

1. **De la cantidad predicha.** Un desplazamiento global promediado tolera mejor la
   idealización que un esfuerzo local. El acoplamiento tracción–cortante es el caso
   extremo: un modelo isótropo lo predice **idénticamente nulo** para cualquier
   orientación, de modo que el error relativo no es grande sino *indefinido* — el
   fenómeno no existe en su formulación.
2. **De la orientación de la fibra** respecto a la carga. Con la fibra alineada y `E`
   calibrado contra el ensayo longitudinal, el error es nulo por construcción; crece al
   girar.
3. **De qué idealización isótropa se elija.** No hay una única: `E = E_L`, `E = E_T`, o
   un ajuste óptimo dan errores distintos y en direcciones distintas. Determinar el
   *mejor isótropo posible* convierte un resultado obvio ("el isótropo se equivoca") en
   uno informativo ("ni siquiera el mejor isótropo baja del X % bajo estas
   condiciones").

## 2. Hipótesis de trabajo

**H1.** Existe un umbral de anisotropía `(E_L/E_T)*` por debajo del cual la idealización
isótropa mantiene el error de desplazamientos bajo un margen ingenieril (a fijar, del
orden del 10 %), y por encima del cual deja de ser admisible.

**H2.** Ese umbral depende fuertemente de la orientación de fibra: para θ = 0° la
idealización es exacta por calibración, y el error crece monótonamente hasta un máximo
en la zona intermedia.

**H3.** El error en cantidades locales (esfuerzos, acoplamiento) es sustancialmente mayor
que en cantidades globales (desplazamientos), de modo que un modelo isótropo calibrado
para reproducir rigidez global puede seguir siendo inadmisible para predecir esfuerzos.

**H4.** La sensibilidad del error a `G_LT` es comparable o mayor que a `E_L/E_T`, lo que
—de confirmarse— justifica priorizar la medición de `G_LT` en la campaña experimental.
Es la hipótesis con mayor valor práctico, porque la revisión de Chen et al. (2024)
señala que las constantes transversales y de cortante del bambú están **poco medidas**.

## 3. Alcance

**Dentro:**
- Elasticidad lineal, pequeñas deformaciones, estado plano de esfuerzo.
- Orientación de fibra uniforme (constante en el dominio).
- Comparación ortótropo vs. isótropo en tres niveles: constitutivo, elemento, estructura.

**Fuera** (declarado, no omitido):
- Plasticidad, daño y falla. El estudio es elástico.
- Ortotropía cilíndrica (ejes que siguen la geometría del culmo). Requiere orientación
  por punto de Gauss — ver ADR 0013 §2.
- Material gradado (propiedades función de la posición radial en la pared).
- Comportamiento 3D. `plane_strain` ortótropo no está disponible y no cierra sin
  constantes 3D (ADR 0013 §3).
- Validación contra datos experimentales propios. **Este estudio es puramente numérico
  y su propósito explícito es informar el diseño de la campaña, no sustituirla.**

## 4. Material: transversalmente isótropo, no ortótropo general

La literatura (Chen et al. 2024) indica que el bambú se comporta como material
**transversalmente isótropo**: los haces vasculares tienen disposición unidireccional y
las direcciones radial y tangencial son mecánicamente similares
(`E_R ≈ E_T`; las relaciones medidas `E_L : E_T : E_R` son 6.05 : 1.08 : 1 en la zona
exterior, 5.12 : 1.17 : 1 en la central y 8.02 : 1.58 : 1 en la interior).

**Consecuencia para este estudio: ninguna operativa.** En estado plano de esfuerzo, un
material transversalmente isótropo con el eje de simetría en el plano se describe con
las mismas cuatro constantes independientes que un ortótropo general
(`E_1`, `E_2`, `G_12`, `ν_12`). La distinción importa para (a) cómo se describe el
material en el artículo y (b) una eventual extensión 3D, donde reduce de nueve a cinco
las constantes necesarias.

El plano de análisis se toma como **L–T** (longitudinal–tangencial), que es el
accesible en tiras extraídas de la pared del culmo.

## 5. Rangos de parámetros

Acotados por literatura, no elegidos por conveniencia. Ver §9 para las fuentes.

| Parámetro | Rango | Base |
|---|---|---|
| `E_L` | 10 – 30 GPa | Moso 16.6 GPa (tracción) y 17.0–35.9 GPa según posición; *D. sinicus* 16.5 GPa (carne) y 30.8 GPa (con capa verde) |
| `E_L/E_T` | **3 – 12** | Relaciones medidas 5.12–8.02, con margen a ambos lados |
| `G_LT/E_L` | 0.02 – 0.10 | Poco documentado — es el hueco que señala la revisión |
| `ν_LT` | 0.18 – 0.34 | Moso, variando con edad, altura y posición radial |
| θ | 0° – 90° | Barrido completo |

> **Nota de calibración.** Los valores por defecto usados en los tests de
> `Orthotropic2D` (`E_1/E_2 = 18.75`) proceden de literatura general de materiales
> ortótropos y son **más severos que el bambú real**. Este estudio usa los rangos de
> arriba, respaldados por medidas sobre bambú.

## 6. Las tres idealizaciones isótropas a comparar

El punto metodológico central. Se compara el ortótropo contra:

| Id | Definición | Representa |
|---|---|---|
| **I-L** | `E = E_L`, `ν = ν_LT` | La práctica habitual: se mide el módulo longitudinal —el dato más disponible— y se usa como si el material fuera isótropo. |
| **I-T** | `E = E_T`, `ν = ν_TL` | La idealización conservadora, tomando la dirección débil. |
| **I-opt** | `E`, `ν` que **minimizan** la métrica de error sobre el rango de θ considerado | La cota inferior: el mejor isótropo posible. Su error residual es el que **no** puede eliminarse con mejor calibración. |

`I-opt` se obtiene por minimización numérica escalar sobre `(E, ν)`, no analíticamente,
para no introducir hipótesis adicionales.

## 7. Métricas de error

Todas **adimensionales** y con signo definido, para poder comparar entre niveles y entre
parámetros con magnitudes distintas.

### 7.1 Nivel constitutivo

**M1 — Error en el módulo aparente:**

```
    e_E(θ) = | E_x^iso − E_x^orto(θ) | / E_x^orto(θ)
```

**M2 — Error en norma de la constitutiva:**

```
    e_C(θ) = ‖ C_iso − C_orto(θ) ‖_F / ‖ C_orto(θ) ‖_F
```

Norma de Frobenius. Captura el desajuste completo de la matriz, no sólo de una
componente.

**M3 — Acoplamiento omitido:**

```
    η(θ) = γ_xy / ε_xx   bajo σ_xx puro
```

No se expresa como error relativo: el isótropo da `η ≡ 0` para todo θ, de modo que el
error relativo sería infinito. Se reporta el **valor absoluto del acoplamiento que la
idealización descarta**, que es la magnitud físicamente interpretable.

### 7.2 Nivel de elemento y estructural

**M4 — Error en desplazamiento:**

```
    e_u = ‖ u_iso − u_orto ‖_2 / ‖ u_orto ‖_2
```

**M5 — Error en esfuerzos, en norma L²  sobre los puntos de Gauss:**

```
    e_σ = [ Σ_g w_g ‖σ_g^iso − σ_g^orto‖² / Σ_g w_g ‖σ_g^orto‖² ]^(1/2)
```

Ponderado por el peso de cuadratura, para que no dependa del número de elementos.

**M6 — Error en la cantidad de diseño**, específica del caso estructural de la Fase 3
(p. ej. esfuerzo máximo en el borde de una perforación). Es la métrica que un ingeniero
usaría para decidir, y por eso se reporta aparte de las globales.

### 7.3 Criterio de admisibilidad

Se declara **antes** de ver resultados, para que el umbral no se ajuste a la conclusión:

- `e < 5 %` — idealización **admisible**.
- `5 % ≤ e < 15 %` — admisible **con reservas**, según la aplicación.
- `e ≥ 15 %` — **inadmisible**.

Estos cortes son convencionales en ingeniería estructural y se declaran como tales, no
como resultado del estudio.

## 8. Fases

| Fase | Contenido | Entregable |
|---|---|---|
| **0** | Este documento + esqueleto de código | Diseño cerrado |
| **1** | Nivel constitutivo: M1, M2, M3 sobre el barrido. Sin FEM. | Mapas de error en (`E_L/E_T`, θ); determinación de `I-opt` |
| **2** | Nivel de elemento: M4, M5 bajo tracción, cortante y biaxial | Verificación cruzada con Fase 1 |
| **3** | Nivel estructural: caso con gradiente real de esfuerzo | Figura del artículo; M6 |
| **4** | Síntesis: umbrales de admisibilidad y prioridades de medición | Respuesta a la pregunta + guía para la campaña |

Cada fase se cierra y se verifica antes de la siguiente. Cada una produce algo utilizable
aunque el estudio se detenga ahí.

## 9. Fuentes de los rangos

- Chen, M. et al. (2024). *Comprehensive review of the elastic constants of bamboo*.
  European Journal of Wood and Wood Products. DOI: 10.1007/s00107-024-02143-6 —
  transversal isotropía; señala que las constantes radiales y transversales están poco
  medidas.
- Estudio de Poisson en Moso (*Phyllostachys heterocycla*), Frontiers in Materials 9
  (2022) 896756 — `ν_LT` = 0.180–0.334 con variación por edad, altura y posición radial.
- *Evaluation of orthotropic elasticity of gradient-structured bamboo*, Industrial Crops
  and Products (2023) — relaciones `E_L : E_T : E_R` por posición radial.
- Estudio de constantes elásticas de *Dendrocalamus sinicus*, Forests 15 (2024) 2017 —
  módulos longitudinales con y sin capa verde.

> **Pendiente de verificación por el autor.** Las referencias anteriores se identificaron
> por búsqueda bibliográfica y sus valores se tomaron de resúmenes y páginas de
> resultados; **varias no fueron accesibles a texto completo**. Antes de publicar, los
> valores citados deben verificarse contra los artículos originales y completarse la
> revisión (Fase 2 del plan general de trabajo, posterior a este estudio).

## 10. Reproducibilidad

- Todo el código del estudio vive en `paper/bambu/estudio/`.
- Cada figura del artículo se genera por script, sin edición manual.
- Los resultados numéricos se serializan a CSV para que las figuras se puedan regenerar
  sin recalcular.
- El estudio usa la versión de `Orthotropic2D` en el repositorio; el commit se registra
  en la salida.
