# Orthotropic2D — Elasticidad lineal ortótropa 2D (plane stress)

> Spec **pre-redactada por la IA** a partir de la sesión 2026-09-09, según `Reglas.md §4` y el criterio de [[feedback-filtrar-puntos-abiertos]]: el plumbing arquitectural lo decide la IA con razón documentada; los puntos abiertos que quedan son **físicos y de alcance**, y los valida el usuario en §Diálogo.
>
> Decisión arquitectural de fondo en [ADR 0013](../adr/0013-orientacion-material-y-ortotropia.md) — dónde vive la orientación material.

---

## Especificación física

### 0. Descripción general

Modelo elástico lineal **ortótropo** bidimensional bajo hipótesis **plane stress**, parametrizado por cuatro constantes independientes $(E_1, E_2, G_{12}, \nu_{12})$ y un ángulo de orientación $\theta$ que sitúa los ejes principales del material respecto a los ejes globales.

Un material ortótropo posee tres planos de simetría material mutuamente perpendiculares. Sus **ejes principales** son las direcciones privilegiadas respecto a las cuales el comportamiento se describe sin acoplamiento. En bambú:

| Eje | Dirección física | Qué gobierna |
|---|---|---|
| **1** ($L$) | Longitudinal, a lo largo del culmo | Haces vasculares — dirección rígida |
| **2** ($R$ o $T$) | Transversal en el plano de análisis | Sólo parénquima |

Independiente de la velocidad, sin variables históricas. Es el primer material anisótropo del catálogo.

**`plane_strain` está fuera de alcance** — ver §Fuera de alcance y ADR 0013 §3.

### 1. Descomposición de la deformación

No aplica — toda la deformación es elástica: $\boldsymbol\varepsilon = \boldsymbol\varepsilon^e$.

### 2. Ley elástica

$$\boldsymbol\sigma = \mathbf C_e\,\boldsymbol\varepsilon$$

con la notación Voigt 2D del proyecto (`Reglas.md §5`): $\boldsymbol\varepsilon = [\varepsilon_{xx},\,\varepsilon_{yy},\,\gamma_{xy}]^\top$ con $\gamma_{xy} = 2\varepsilon_{xy}$ *engineering*, y $\boldsymbol\sigma = [\sigma_{xx},\,\sigma_{yy},\,\sigma_{xy}]^\top$.

#### 2.1 Flexibilidad en ejes del material

La forma inversa es donde las constantes tienen significado físico directo y medible:

$$
\begin{Bmatrix} \varepsilon_{11} \\ \varepsilon_{22} \\ \gamma_{12} \end{Bmatrix}
=
\begin{bmatrix}
1/E_1 & -\nu_{21}/E_2 & 0 \\
-\nu_{12}/E_1 & 1/E_2 & 0 \\
0 & 0 & 1/G_{12}
\end{bmatrix}
\begin{Bmatrix} \sigma_{11} \\ \sigma_{22} \\ \sigma_{12} \end{Bmatrix}
$$

Los ceros del bloque cortante expresan que **en ejes del material no hay acoplamiento** entre esfuerzos normales y distorsión angular.

#### 2.2 Reciprocidad

La simetría de la matriz de flexibilidad —consecuencia de la existencia de una energía de deformación— impone:

$$\boxed{\dfrac{\nu_{12}}{E_1} = \dfrac{\nu_{21}}{E_2}} \qquad\Longrightarrow\qquad \nu_{21} = \nu_{12}\,\frac{E_2}{E_1}$$

Por eso hay **cuatro** constantes independientes y no cinco: $\nu_{21}$ queda determinada.

Consecuencia contraintuitiva pero correcta: como $E_1 \gg E_2$ en bambú, $\nu_{21} \ll \nu_{12}$. Con $E_1/E_2 = 20$ y $\nu_{12} = 0.35$ resulta $\nu_{21} = 0.0175$. **Valores de Poisson diminutos son normales** en materiales muy ortótropos, no error de medición.

#### 2.3 Rigidez en ejes del material

Invirtiendo la flexibilidad:

$$
\mathbf C_{\text{mat}} =
\begin{bmatrix}
\dfrac{E_1}{1-\nu_{12}\nu_{21}} & \dfrac{\nu_{12}E_2}{1-\nu_{12}\nu_{21}} & 0 \\[6pt]
\dfrac{\nu_{12}E_2}{1-\nu_{12}\nu_{21}} & \dfrac{E_2}{1-\nu_{12}\nu_{21}} & 0 \\[6pt]
0 & 0 & G_{12}
\end{bmatrix}
$$

El factor $1-\nu_{12}\nu_{21}$ juega el papel estructural que $1-\nu^2$ tiene en el isótropo.

#### 2.4 Rotación a ejes globales

Con $\theta$ el ángulo de los ejes del material respecto a los globales (positivo antihorario, coherente con `Reglas.md §5`), $c=\cos\theta$, $s=\sin\theta$, la transformación de **deformaciones** en Voigt *engineering* es:

$$
\mathbf T =
\begin{bmatrix}
c^2 & s^2 & cs \\
s^2 & c^2 & -cs \\
-2cs & 2cs & c^2-s^2
\end{bmatrix}
\qquad\Longrightarrow\qquad
\boxed{\mathbf C_{\text{glob}} = \mathbf T^\top\,\mathbf C_{\text{mat}}\,\mathbf T}
$$

Para $\theta \neq 0$ la matriz **se llena por completo**:

$$
\mathbf C_{\text{glob}} =
\begin{bmatrix}
C_{11} & C_{12} & C_{16} \\
C_{12} & C_{22} & C_{26} \\
C_{16} & C_{26} & C_{66}
\end{bmatrix}
$$

Los términos $C_{16}, C_{26}$ —nulos en ejes del material— son el **acoplamiento tracción–cortante**: un $\varepsilon_{xx}$ puro genera $\sigma_{xy}$. Físicamente, una tira cortada a 45° de la fibra no sólo se alarga al traccionarla: se tuerce. Este acoplamiento no existe en ningún material isótropo.

$\mathbf C_{\text{glob}}$ **sigue siendo simétrica** (lo garantiza la energía de deformación), luego `IS_SYMMETRIC = True` y el despachador algebraico del ADR 0003 la trata igual que a las isótropas.

### 3-6. Fluencia, flujo, endurecimiento, Kuhn-Tucker

No aplican — el material es indefinidamente elástico.

### 7. Variables internas

Ninguna. `PRIMARY_STATE_VAR = None`, como el resto de materiales puramente elásticos del catálogo.

### 8. Notación Voigt

La del proyecto (`Reglas.md §5`), sin componentes auxiliares fuera del plano: bajo plane stress $\sigma_{33}=0$ y $\varepsilon_{33}$ no interviene en la respuesta en el plano.

**Convención de subíndices de Poisson — fijada explícitamente.** Existen dos convenciones opuestas en la literatura y confundirlas **invierte los papeles de $E_1$ y $E_2$** sin fallo ruidoso. Se adopta la de Jones/Tsai:

$$\nu_{12} = -\frac{\varepsilon_{22}}{\varepsilon_{11}} \quad\text{bajo carga uniaxial en la dirección 1}$$

**Primer índice = dirección de carga. Segundo = dirección de la contracción medida.**

---

## Formulación numérica

### 9-10. Esquema temporal y actualización de estado

No aplican — material sin historia. $\mathbf C_{\text{glob}}$ se **precalcula una sola vez en el constructor** (igual que `Elastic2D` precalcula `self.C` según la hipótesis) y `compute_state` devuelve $(\mathbf C\boldsymbol\varepsilon,\ \mathbf C,\ \text{state\_vars})$ sin recomputar nada.

### 11. Tangente

Constante e igual a $\mathbf C_{\text{glob}}$ en todo régimen. Simétrica y definida positiva dentro del rango admisible (§12).

### 12. Admisibilidad y caveats

**No se puede reutilizar la validación de `Elastic2D`.** Ese material valida $-1<\nu<0.5$, criterio correcto para isótropo e **incorrecto aquí**. La condición ortótropa es que $\mathbf C$ sea definida positiva:

$$E_1,\,E_2,\,G_{12} > 0 \qquad\text{y}\qquad 1-\nu_{12}\nu_{21} > 0 \iff |\nu_{12}| < \sqrt{E_1/E_2}$$

Con $E_1/E_2=20$ el límite es $|\nu_{12}|<4.47$: **$\nu_{12}>0.5$ es físicamente legítimo** en un material muy ortótropo. Un validador heredado del isótropo rechazaría datos válidos de bambú.

**Caveat de implementación — el factor 2 del cortante.** Con $\gamma$ *engineering*, las matrices de transformación de esfuerzo y de deformación **no coinciden**: difieren en factores 2 en el bloque cortante. La forma $\mathbf T^\top\mathbf C\mathbf T$ es válida precisamente para la $\mathbf T$ de deformaciones de §2.4. Es el error clásico al programar esto y debe quedar blindado por el test de invariancia (§`acceptance`).

**$G_{12}$ es independiente**: no se deriva de $E$ y $\nu$ como en isótropo. Debe medirse.

---

## Contrato de implementación

```yaml
name: Orthotropic2D
kind: material
status: validated

interface:
  strain_dim: 3
  primary_state_var: null
  is_symmetric: true

parameters:
  - { name: E1,    type: float, required: true,  desc: "Módulo de Young en la dirección 1 (fibra). > 0" }
  - { name: E2,    type: float, required: true,  desc: "Módulo de Young en la dirección 2 (transversal). > 0" }
  - { name: G12,   type: float, required: true,  desc: "Módulo de cortante en el plano 1-2. Independiente de E y nu. > 0" }
  - { name: nu12,  type: float, required: true,  desc: "Poisson mayor: carga en 1, contracción medida en 2. Admisible |nu12| < sqrt(E1/E2)" }
  - { name: theta, type: float, required: false, default: 0.0, desc: "Ángulo en GRADOS de los ejes del material respecto a los globales, positivo antihorario. ADR 0013: atajo consciente, migrará al elemento." }
  - { name: density, type: float, required: false, default: null, desc: "Densidad (ADR 0008). Opcional; requerida por peso propio y matriz de masa." }

state_variables: []

conventions:
  voigt: "[eps_xx, eps_yy, gamma_xy] con gamma_xy = 2*eps_xy (engineering); [sig_xx, sig_yy, sig_xy]"
  units: "E1, E2, G12 en Pa; nu12 adimensional; theta en grados; density en kg/m^3"
  sign: "Tracción positiva. theta positivo antihorario (Reglas.md §5). nu_ij: primer índice = dirección de carga (Jones/Tsai)."

validity:
  - "Pequeñas deformaciones, elasticidad lineal indefinida"
  - "Hipótesis plane_stress únicamente"
  - "Ortotropía con ejes principales constantes en todo el material (theta uniforme)"

out_of_scope:
  - "plane_strain — requiere E3, nu13, nu23; el problema plano no cierra sin datos 3D (ADR 0013 §3)"
  - "Orthotropic3D — nueve constantes independientes, rotación 6x6; cuando lo pida un caso real"
  - "Criterios de falla ortótropos (Tsai-Wu, Hashin, Hoffman) — componente distinto y posterior"
  - "Material gradado (propiedades función de la posición)"
  - "Ortotropía cilíndrica (ejes que siguen la geometría del culmo) — requiere orientación por punto de Gauss, ADR 0013 §2"
  - "Orientación variable por elemento — migración a ADR 0013 opción B"

acceptance:
  verification:
    - name: degeneracion_al_isotropo
      setup: "Construir con E1 = E2 = E, nu12 = nu, G12 = E/(2(1+nu)), theta = 0.
             Comparar C contra Elastic2D(E, nu, 'plane_stress').C"
      expect: "Coincidencia exacta: el isótropo es el caso particular del ortótropo"
      tol_rel: 1.0e-14
    - name: respuesta_elastica_en_ejes
      setup: "theta = 0. Deformación uniaxial eps_xx impuesta, resto nulo"
      expect: "sigma = C_mat @ eps exacto; sin acoplamiento: sigma_xy = 0"
      tol_rel: 1.0e-12
    - name: tangente_constante
      setup: "Evaluar compute_state con varias deformaciones distintas"
      expect: "La tangente devuelta es idéntica en todas: material sin historia"
      tol_rel: 1.0e-14
  specific:
    - name: acoplamiento_traccion_cortante
      setup: "theta = 45 grados, deformación uniaxial pura eps_xx"
      expect: "sigma_xy != 0 (acoplamiento). Verificar contra C_glob = T^T C_mat T
              calculada independientemente en el test"
      tol_rel: 1.0e-12
    - name: sin_acoplamiento_en_ejes
      setup: "theta = 0 y theta = 90 grados"
      expect: "C[0,2] = C[1,2] = 0 en ambos casos (ejes principales alineados);
              en theta=90 los papeles de E1 y E2 quedan intercambiados"
      tol_rel: 1.0e-14
    - name: invariancia_rotacion_material_y_carga
      setup: "Rotar material y estado de deformación el mismo ángulo; comparar la
             energía de deformación 0.5*eps^T C eps y los invariantes de sigma"
      expect: "Invariantes idénticos: la respuesta física no depende del sistema
              de referencia. Blinda el factor 2 del cortante engineering (§12)"
      tol_rel: 1.0e-12
    - name: simetria_y_definicion_positiva
      setup: "Barrido de theta en [0, 180) grados con datos típicos de bambú"
      expect: "C simétrica a precisión máquina y autovalores > 0 para todo theta"
      tol_rel: 1.0e-14
    - name: reciprocidad
      setup: "Verificar nu21 derivado internamente"
      expect: "nu12/E1 == nu21/E2"
      tol_rel: 1.0e-14
    - name: admisibilidad_rechaza_invalidos
      setup: "E1<=0, E2<=0, G12<=0, |nu12| >= sqrt(E1/E2)"
      expect: "ValueError con mensaje que nombra el parámetro y el criterio violado"
      tol_rel: null
    - name: poisson_mayor_que_medio_es_admisible
      setup: "E1/E2 = 20, nu12 = 0.6 (excede el límite isótropo de 0.5)"
      expect: "Se construye sin error: |nu12| < sqrt(20) = 4.47. Blinda que NO se
              copió la validación de Elastic2D"
      tol_rel: null

references:
  - "Jones R.M. (1999). Mechanics of Composite Materials, 2a ed. Taylor & Francis. Cap. 2 (constitutiva ortótropa 2D, reciprocidad, transformación de ejes)."
  - "Lekhnitskii S.G. (1963). Theory of Elasticity of an Anisotropic Elastic Body. Holden-Day. (Restricciones de admisibilidad por definición positiva.)"
  - "Ting T.C.T. (1996). Anisotropic Elasticity: Theory and Applications. Oxford University Press."
```

---

## Diálogo

### Puntos abiertos para el usuario

Ninguno bloqueante — las tres decisiones de alcance (sólo `plane_stress`, sólo 2D, sólo elasticidad) se acordaron en sesión el 2026-09-09 y están registradas en el ADR 0013 §3.

Queda una cuestión **física** que la campaña experimental del PAPIIT resolverá y que conviene anotar aquí cuando haya datos:

1. **Valores de referencia de la especie estudiada.** Los órdenes de magnitud usados en los tests son de literatura general, no de una especie concreta:

   | Constante | Orden típico en bambú |
   |---|---|
   | $E_1$ | 10–20 GPa |
   | $E_2$ | 0.5–1 GPa |
   | $G_{12}$ | 0.5–1 GPa |
   | $\nu_{12}$ | 0.3–0.4 |

   La relación $E_1/E_2 \approx 20$ sitúa al bambú en el mismo régimen que un laminado unidireccional de carbono (15–20) o la madera (15–25). Cuando existan valores medidos, sustituirlos en el test de barrido y anotar la especie.

2. **¿El plano de análisis es $L$–$R$ o $L$–$T$?** Ambos son ortótropos 2D válidos y la spec los cubre igual, pero $E_2$ toma valores distintos según cuál sea. Es dato de la campaña, no decisión de diseño.

---

## Implementación

*Rellena la IA tras programar.*

- **Archivo**: [`solidum/materials/orthotropic_2d.py`](../../solidum/materials/orthotropic_2d.py)
- **Clase**: `Orthotropic2D` (auto-registrada vía `@MaterialRegistry.register`)
- **Tests (verificación)**: [`tests/test_orthotropic_2d.py`](../../tests/test_orthotropic_2d.py) — 26 tests + 57 subtests.
- **Validación externa**: [`tests/validation/test_off_axis_orthotropic.py`](../../tests/validation/test_off_axis_orthotropic.py) — 16 tests + 60 subtests contra las formas cerradas de **Jones (1999) §2.8**: módulo aparente `E_x(θ)` (ec. 2.85, 1e-12 rel.), Poisson aparente `ν_xy(θ)` (ec. 2.86, 1e-10) y coeficiente de influencia mutua `η_xy,x(θ)`. Incluye consistencia FEM sobre un `Quad4` real (1e-10), no sólo sobre la constitutiva. Cifra citable: a 45° la rigidez cae **por debajo del 10 % de E₁**, muy por debajo de lo que sugiere `E₁/E₂ = 18.75`, porque el término cruzado `(1/G₁₂ − 2ν₁₂/E₁)` domina la zona intermedia. **Mutation test documentado**: quitar el factor 2 del bloque cortante de `T` produce 39 fallos — el benchmark no pasa por construcción.
- Suite global 1176 → 1220.
- **Notas de implementación**:
  - `C_mat` y `C` (rotada) se **precalculan en el constructor**; `compute_state` no recomputa nada. Mismo patrón que `Elastic2D` con su hipótesis.
  - **Tolerancias adimensionales en todos los tests** (ADR 0006). Primera versión usaba `assertAlmostEqual(..., places=6)` absoluto y fallaba en `theta ∈ {90°, 180°}` con residuos de ~1e-6 sobre magnitudes de ~1e10 Pa — que son **1e-16 en relativo**, precisión máquina. El origen es que `sin(180°)` no da cero exacto en punto flotante. Comparar contra cero absoluto convertía el test en una medida de la magnitud de `E1`, no de la física. Verificado que el criterio es invariante ante cambio de unidades en 15 órdenes de magnitud (Pa, MPa, micro-unidades).
  - **Bring-up end-to-end** con `Quad4` + `LinearSolver`: `E_aparente` = 15.10 GPa (θ=0°), 4.82 GPa (θ=45°), 0.805 GPa (θ=90°), con `σ_xy ≠ 0` sólo fuera de ejes principales y `K` simétrica en los tres casos.
- **Notas de traducción**: ninguna — la formulación de Jones se implementa en la convención del proyecto sin cambio de signos.
