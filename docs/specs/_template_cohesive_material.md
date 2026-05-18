# <Nombre> — <descripción corta>

> Orden de trabajo. El usuario escribe la **especificación física**, la **formulación numérica** y el **contrato**. La IA rellena **implementación** y responde en **diálogo** durante el trabajo.
>
> Los materiales cohesivos forman una jerarquía paralela a la de materiales continuos (ADR 0010): operan sobre el salto `[[u]]` y devuelven tracciones `t` sobre `Γ_d`. No tienen `STRAIN_DIM` ni `density` (la inercia es del bulk). En su lugar declaran `JUMP_DIM`. Viven en `fenix/cohesive_materials/`, heredan de `fenix.core.cohesive_material.CohesiveMaterial`, se registran con `@CohesiveMaterialRegistry.register` y se declaran en YAML bajo `cohesive_materials:` (sección separada de `materials:`).

---

## Especificación física

### 0. Descripción general
*Tipo de modelo cohesivo (daño escalar, daño anisótropo, plasticidad cohesiva, mixed mode, etc.); modo de fractura cubierto (Modo-I puro, mixto I–II, …); regímen de aplicación (cuasi-frágil, dúctil, cíclico). `JUMP_DIM` ∈ {2 (2D), 3 (3D)}. Si reutiliza la familia introducida por ADR 0010 indicarlo.*

### 1. Descomposición del salto
*Cómo se proyecta `[[u]] ∈ ℝ^{JUMP_DIM}` sobre el frame local de `Γ_d`. Convención del proyecto: `n` normal a `Γ_d` (fijada al activarse el elemento), `s` (2D) o `(s1, s2)` (3D) tangentes. Si el modelo separa parte elástica y plástica del salto (`[[u]] = [[u]]^e + [[u]]^p`), documentarlo aquí.*

### 2. Ley elástica del salto (sin daño / plasticidad)
*$\mathbf t = \mathbf T^{el}\,\llbracket\mathbf u\rrbracket$ con $\mathbf T^{el}$ rigidez elástica del salto. En Modo-I puro: $\mathbf T^{el} = K_e\,(\mathbf n \otimes \mathbf n)$, rigidez `K_e` en normal y cero en tangencial. `K_e` actúa como penalty: balance entre `κ_0` despreciable (`K_e` grande) y condicionamiento del sistema (`K_e` moderado).*

### 3. Variable equivalente y criterio de activación / fluencia
*Magnitud escalar derivada del salto que controla la activación o disipación (ej. `⟨[[u_n]]⟩` para Modo-I Rankine, `‖[[u]]‖` para mixto isótropo). Función de activación/fluencia $f(\llbracket\mathbf u\rrbracket, q) \le 0$ y variables internas $q$.*

### 4. Variable histórica y evolución (Kuhn-Tucker)
*Si aplica: $\kappa_{n+1} = \max(\kappa_n, ...)$. Condiciones de complementariedad $\dot\omega \ge 0$, $f \le 0$, $\dot\omega\,f = 0$. Comportamiento en descarga, recarga, cierre/compresión.*

### 5. Ley de evolución del daño / plasticidad cohesiva
*Forma cerrada o algorítmica de $\omega(\kappa)$ (daño escalar) o de las variables internas plásticas cohesivas. Ley de softening (lineal, exponencial, bilineal). Energía de fractura $G_F$: cómo cierra la curva $t–[[u_n]]$ analíticamente ($\int t\,d\llbracket u_n\rrbracket = G_F$).*

### 6. Variables internas
*Lista con tipo, dimensión y significado físico. Cuál es `PRIMARY_STATE_VAR` (la variable exportada al post-procesado: típicamente `damage`).*

### 7. Frame local y notación
*Convención del proyecto: `[[u]]` y `t` se reciben/devuelven en el frame local `(n, s)` o `(n, s1, s2)` de `Γ_d`. La transformación entre globales y locales es responsabilidad del elemento consumidor (no del material). Documentar signos: `[[u_n]] > 0 ⇔ apertura (tracción)`.*

---

## Formulación numérica

### 8. Algoritmo (explícito vs implícito local)
*Si el modelo permite actualización cerrada (ej. daño escalar Modo-I): describir los pasos. Si requiere Newton local (plasticidad cohesiva, modo mixto con superficie no proyectable): describir el sistema y el predictor.*

### 9. Tangente algorítmica consistente
*$\mathbf T^\text{alg} = \partial\mathbf t_{n+1}/\partial\llbracket\mathbf u\rrbracket_{n+1}$ derivada del algoritmo. Indicar simetría: si `IS_SYMMETRIC = False`, el despachador algebraico ADR 0003 elige LU. Para Modo-I escalar puede salir simétrica (`α·(n⊗n)`); modo mixto típicamente no lo es. Documentar las ramas (sin daño, carga activa, descarga, saturación) y la forma cerrada de cada una.*

### 10. Simetría / asimetría
*Justificar `IS_SYMMETRIC`. Si la tangente es simétrica por construcción, indicar cuándo dejaría de serlo (ej. extender a modo mixto).*

### 11. Decisión de régimen durante el Newton global
*Cómo decide el material carga vs descarga comparando con el estado committed `κ_n`. Si la transición es discontinua, advertirlo y sugerir solver (Newton estándar con line-search, `ArcLengthSolver`).*

### 12. Caveats numéricos
*Elección del penalty `K_e` (típicamente `K_e ≈ 10·E_bulk/ℓ_c`). Comportamiento en compresión (¿hay contacto unilateral o se reduce la rigidez?). Cap por `DAMAGE_MAX` para evitar singularidad: dónde y cómo se aplica (recordatorio del patrón en `CohesiveDamageIsotropic` — sólo a la rigidez tangente, no a `ω` ni a `t`). Mesh-objectivity en softening: se aborda a nivel del elemento (longitud característica), no del material.*

---

## Contrato de implementación

```yaml
name: <Nombre>
kind: cohesive_material  # nueva familia ADR 0010
status: draft            # draft → implemented → validated

interface:
  jump_dim: 2            # 2 en 2D, 3 en 3D
  primary_state_var: ""  # nombre de la variable interna exportada al post
  is_symmetric: true     # tangente algorítmica simétrica? (afecta al despachador algebraico ADR 0003)

parameters:
  - { name: , type: , required: true, desc:  }

state_schema:
  # Variables internas que el material guarda y actualiza entre commits.

conventions:
  sign: ""               # convención de signos (apertura positiva ⇒ tracción, etc.)
  frame: ""              # frame local (n, s) en 2D; transformación a globales por el elemento
  units: ""              # unidades de cada parámetro (SI por defecto)

validity:
  - ""                   # regímenes de aplicación válidos (ej. σ_t0 > 0, K_e > σ_t0² / (2·G_F))

out_of_scope:
  - ""                   # hipótesis no soportadas, características diferidas
                          # (modo mixto, contacto unilateral, anisotropía, rate-dependence, fatiga…)

numerical_caveats:
  - ""                   # penalty K_e, cap por DAMAGE_MAX, mesh-objectivity, etc.

acceptance:
  # Bloque de verificación obligatorio. La disciplina V&V para cohesivos:
  # respuesta elástica bajo el umbral; inicio del daño en κ_0; integral
  # ∫t·d[[u_n]] = G_F (validación energética); ciclos carga/descarga/recarga;
  # simetría/asimetría de la tangente; consistencia por diferencias finitas;
  # comportamiento en tangencial y compresión; saturación.
  verification:
    - name: respuesta_elastica_bajo_umbral
      setup: ""
      expect: ""
      tol_rel: 1.0e-12
    - name: inicio_dano_en_umbral
      setup: ""
      expect: ""
      tol_rel: 1.0e-10
    - name: energia_fractura
      setup: "carga monótona [[u_n]] hasta apertura completa; cuadratura de la curva t–[[u_n]]"
      expect: "∫ t·d[[u_n]] = G_F"
      tol_rel: 1.0e-3
    - name: descarga_secante
      setup: ""
      expect: ""
      tol_rel: 1.0e-10

  specific:
    - name: tangente_simetria
      setup: ""
      expect: ""
      tol_rel: 1.0e-14
    - name: tangente_diferencia_finita_consistente
      setup: ""
      expect: ""
      tol_rel: 1.0e-5
    - name: compresion_no_dana
      setup: ""
      expect: ""
      tol_rel: 1.0e-12

  arch:
    - name: registry_kind_es_cohesive
      expect: "Registro vía CohesiveMaterialRegistry, no MaterialRegistry; kind='cohesive_material' en YAML"
    - name: is_symmetric_attribute
      expect: "<NombreClase>.IS_SYMMETRIC == <true/false>"

references:
  - ""                   # bibliografía primaria; tesis del usuario si aplica
  - "ADR 0010 — Discontinuidades interiores embebidas (familia CohesiveMaterial, hoja de ruta)."
```

---

## Implementación

*Rellena la IA tras programar.*

- Archivo: —
- Clase: —
- Tests:
  - —
- Notas de traducción: —

---

## Diálogo

*Preguntas, aclaraciones y hallazgos durante la implementación. Entradas fechadas.*

- *(vacío)*
