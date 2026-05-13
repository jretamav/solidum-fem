# ISOTROPIC DAMAGE 2D — daño escalar con ablandamiento exponencial

> Orden de trabajo. El usuario escribe **especificación física**, **formulación** y **contrato**; la IA rellena **implementación** y responde en **diálogo**.
>
> Spec **retroactiva con ampliación**: el modelo está implementado desde sesiones anteriores con **tangente secante**. Esta spec lo documenta y añade la **tangente algorítmica consistente** para recuperar convergencia cuadrática del Newton global.

---

## Especificación física

### 0. Descripción general

Modelo de **daño isótropo escalar** en pequeñas deformaciones para sólidos 2D. Captura degradación progresiva e isótropa de la rigidez en respuesta a la deformación máxima histórica, con ablandamiento exponencial. Pensado para hormigón en tracción uniforme, cerámicos, materiales cuasi-frágiles donde la microfisuración distribuida puede modelarse mediante un escalar de daño en lugar de microestructura explícita.

### 1. Cinemática

Pequeñas deformaciones, descomposición elástica pura: la deformación plástica no existe en este modelo (sin plasticidad). El daño se manifiesta como degradación de la rigidez efectiva.

### 2. Esfuerzo nominal vs esfuerzo efectivo

$$\boldsymbol\sigma = (1 - d)\,\boldsymbol\sigma^\text{eff}, \qquad \boldsymbol\sigma^\text{eff} = \mathbf C_e\,\boldsymbol\varepsilon$$

con $d \in [0, d_\text{max}]$ la variable de daño isótropa y $\mathbf C_e$ el tensor constitutivo elástico isótropo (2D, plane stress o plane strain según `hypothesis` heredada del `Elastic2D` interno).

El esfuerzo efectivo $\boldsymbol\sigma^\text{eff}$ es el que el material intacto soportaría con la deformación corriente. El esfuerzo nominal $\boldsymbol\sigma$ es lo que efectivamente transmite considerando la degradación.

### 3. Deformación equivalente

Cantidad escalar que cuantifica la "intensidad" de la deformación. Convención del proyecto, norma simétrica simple en notación Voigt $[\varepsilon_{xx}, \varepsilon_{yy}, \gamma_{xy}]$ con $\gamma_{xy} = 2\varepsilon_{xy}$:

$$\varepsilon_\text{eq}(\boldsymbol\varepsilon) = \sqrt{\varepsilon_{xx}^2 + \varepsilon_{yy}^2 + \tfrac{1}{2}\gamma_{xy}^2} = \sqrt{\boldsymbol\varepsilon^\top \mathbf M\,\boldsymbol\varepsilon}$$

con $\mathbf M = \operatorname{diag}(1, 1, 1/2)$. Equivalente a la norma de Frobenius del tensor deformación 2D plana ($\boldsymbol\varepsilon : \boldsymbol\varepsilon$ con $\varepsilon_{xy} = \gamma_{xy}/2$).

**Limitación**: no distingue tracción de compresión. Para hormigón a compresión se necesitaría un split tipo Mazars o de Vree; aquí se diferencia explícitamente como *out-of-scope*.

### 4. Variable histórica y condición de carga (Kuhn-Tucker)

$$\kappa_{n+1} = \max(\kappa_n, \varepsilon_\text{eq}(\boldsymbol\varepsilon_{n+1})), \qquad \kappa_0 \text{ inicial}$$

Esto codifica:
- **Carga activa** ($\dot\kappa > 0$): $\kappa$ crece con $\varepsilon_\text{eq}$; el daño puede aumentar.
- **Descarga / recarga** ($\dot\kappa = 0$): $\kappa$ queda congelado; el daño no cambia.
- **Irreversibilidad**: $\kappa$ es monótono no decreciente.

### 5. Ley de evolución del daño

$$d(\kappa) = \begin{cases}
0 & \kappa \leq \kappa_0 \\
1 - \dfrac{\kappa_0}{\kappa}\,e^{-\alpha(\kappa - \kappa_0)} & \kappa > \kappa_0
\end{cases}$$

con saturación numérica $d \leq d_\text{max} = $ `DAMAGE_MAX` (`fenix.constants.DAMAGE_MAX`, evita rigidez nula).

Propiedades:
- $d(\kappa_0) = 0$ (continuidad).
- $d'(\kappa) > 0$ para $\kappa > \kappa_0$ (irreversibilidad).
- $\lim_{\kappa \to \infty} d = 1$ (asintótico).
- $\alpha$ controla la velocidad de degradación (ablandamiento más abrupto si $\alpha$ grande).

---

## Formulación numérica

### 6. Algoritmo en un paso

Dado $\{\kappa_n\}$ y $\boldsymbol\varepsilon_{n+1}$:

1. $\varepsilon_\text{eq} = \sqrt{\varepsilon_{xx}^2 + \varepsilon_{yy}^2 + \tfrac{1}{2}\gamma_{xy}^2}$.
2. Si $\varepsilon_\text{eq} > \kappa_n$: **carga activa**, $\kappa_{n+1} = \varepsilon_\text{eq}$. Si no: **descarga**, $\kappa_{n+1} = \kappa_n$.
3. Aplicar ley de daño con $\kappa_{n+1}$; saturar a $d_\text{max}$.
4. $\boldsymbol\sigma_{n+1} = (1 - d)\,\mathbf C_e\,\boldsymbol\varepsilon_{n+1}$.
5. Tangente: **secante** en descarga / sin daño / saturación; **algorítmica consistente** en carga activa con daño no saturado (ver §7).

El algoritmo es **explícito**: no requiere Newton local. La actualización de $\kappa$ y $d$ es directa.

### 7. Tangente algorítmica consistente (ampliación)

$$\mathbf C^\text{alg} = \frac{\partial \boldsymbol\sigma_{n+1}}{\partial \boldsymbol\varepsilon_{n+1}}$$

Derivando $\boldsymbol\sigma = (1-d)\,\mathbf C_e\,\boldsymbol\varepsilon$:

$$\mathbf C^\text{alg} = (1-d)\,\mathbf C_e - \boldsymbol\sigma^\text{eff} \otimes \frac{\partial d}{\partial \boldsymbol\varepsilon}$$

con $\boldsymbol\sigma^\text{eff} = \mathbf C_e\,\boldsymbol\varepsilon$. La derivada del daño respecto a la deformación, por regla de la cadena:

$$\frac{\partial d}{\partial \boldsymbol\varepsilon} = \frac{\partial d}{\partial \kappa}\,\frac{\partial \kappa}{\partial \boldsymbol\varepsilon}$$

donde:

- $\dfrac{\partial d}{\partial \kappa} = (1 - d)\left(\dfrac{1}{\kappa} + \alpha\right)$ para la ley exponencial. (Demostración: derivando $d = 1 - (\kappa_0/\kappa)e^{-\alpha(\kappa-\kappa_0)}$ y reagrupando vía $1 - d = (\kappa_0/\kappa)e^{-\alpha(\kappa-\kappa_0)}$.)

- $\dfrac{\partial \kappa}{\partial \boldsymbol\varepsilon}$ depende del régimen:
  - **Carga activa** ($\varepsilon_\text{eq} > \kappa_n$, $\kappa = \varepsilon_\text{eq}$):
    $$\dfrac{\partial \kappa}{\partial \boldsymbol\varepsilon} = \dfrac{\partial \varepsilon_\text{eq}}{\partial \boldsymbol\varepsilon} = \dfrac{1}{\varepsilon_\text{eq}}\begin{bmatrix}\varepsilon_{xx} \\ \varepsilon_{yy} \\ \tfrac{1}{2}\gamma_{xy}\end{bmatrix}$$
  - **Descarga** ($\varepsilon_\text{eq} \leq \kappa_n$): $\partial \kappa / \partial \boldsymbol\varepsilon = 0$.

**Tangente final**:

$$\mathbf C^\text{alg} = \begin{cases}
\mathbf C_e & \kappa \leq \kappa_0 \;(\text{sin daño}) \\
(1-d)\,\mathbf C_e & \text{descarga} \\
(1-d)\,\mathbf C_e \;-\; \dfrac{(1-d)(1/\kappa + \alpha)}{\varepsilon_\text{eq}}\,(\mathbf C_e\boldsymbol\varepsilon)(\mathbf M\boldsymbol\varepsilon)^\top & \text{carga activa, } d < d_\text{max} \\
(1-d_\text{max})\,\mathbf C_e & d = d_\text{max} \;(\text{saturación})
\end{cases}$$

con $\mathbf M = \operatorname{diag}(1, 1, 1/2)$.

### 8. Pérdida de simetría

El producto externo $(\mathbf C_e\boldsymbol\varepsilon)(\mathbf M\boldsymbol\varepsilon)^\top$ es una matriz $3\times 3$ que **no es simétrica en general**. $\mathbf C_e\boldsymbol\varepsilon$ y $\mathbf M\boldsymbol\varepsilon$ no son proporcionales salvo en estados de deformación muy particulares (e.g. uniaxial puro alineado con ejes).

Esto rompe la simetría heredada de la tangente secante $(1-d)\mathbf C_e$:

- **Tangente secante**: simétrica. `IS_SYMMETRIC = True` sería válido.
- **Tangente consistente**: no simétrica en general. **`IS_SYMMETRIC = False` obligatorio.**

Consecuencia algebraica (ADR 0003): el despachador elige LU en lugar de Cholesky/LDLᵀ para sistemas globales que incluyan este material. El coste por paso aumenta ~30-50% respecto a Cholesky, pero la convergencia cuadrática del Newton global con tangente consistente recupera ese coste con creces (típicamente 2-3 iteraciones en lugar de 5-10 con tangente secante).

### 9. Decisión de régimen durante el Newton

Durante una iteración del Newton global, el material recibe un $\boldsymbol\varepsilon$ trial y debe decidir carga/descarga comparando con $\kappa_n$ (el estado committed del paso anterior). La transición entre ramas (loading ↔ unloading) es **discontinua en la tangente**: al cruzar $\varepsilon_\text{eq} = \kappa_n$ el segundo término aparece/desaparece bruscamente. En la práctica esto no genera problemas porque el Newton converge antes de oscilar entre ramas; en problemas patológicos puede requerir continuation o arc-length.

---

## Contrato de implementación

```yaml
name: IsotropicDamage2D
kind: material
status: validated

interface:
  strain_dim: 3
  primary_state_var: damage
  is_symmetric: false      # tangente consistente no simétrica (ver §8)

parameters:
  - { name: E,          type: float, required: true,  desc: "Módulo de Young intacto" }
  - { name: nu,         type: float, required: true,  desc: "Coeficiente de Poisson" }
  - { name: kappa_0,    type: float, required: true,  desc: "Umbral de deformación equivalente; daño nulo mientras κ ≤ κ_0" }
  - { name: alpha,      type: float, required: true,  desc: "Velocidad de degradación exponencial (>0)" }
  - { name: hypothesis, type: str,   required: false, default: "plane_stress", desc: "plane_stress | plane_strain (delegado a Elastic2D interno)" }
  - { name: density,    type: float, required: false, default: null, desc: "Densidad (ADR 0008)" }

signature:
  compute_state: "(ε: ndarray(3,), state_vars=None) -> (σ: ndarray(3,), C_tan: ndarray(3,3), state_vars')"
  strain_kind: "plano Voigt — [ε_xx, ε_yy, γ_xy] con γ_xy = 2·ε_xy"

state_schema:
  kappa:  "float ≥ κ_0 — máxima deformación equivalente histórica"
  damage: "float ∈ [0, DAMAGE_MAX] — variable de daño escalar isótropa"

conventions:
  sign: "σ > 0 ⇔ tracción; el modelo no distingue tracción de compresión en la activación del daño"
  voigt: "γ_xy = 2·ε_xy en deformación; ε_eq usa norma M=diag(1,1,1/2)"
  damage_irreversibility: "κ monótono no decreciente vía max(κ_old, ε_eq); d(κ) creciente"

validity:
  - "E > 0, κ_0 > 0, α > 0"
  - "ν ∈ (-1, 0.5) (delegación a Elastic2D)"
  - "hypothesis ∈ {'plane_stress', 'plane_strain'}"

out_of_scope:
  - "split tracción/compresión (Mazars, de Vree, Mises modificado); el modelo daña por igual en cualquier régimen"
  - "anisotropía del daño (tensor de daño D vs escalar)"
  - "regularización para evitar mesh-dependency en regime de ablandamiento (longitud característica, no-local, gradient)"
  - "acoplamiento con plasticidad"
  - "fatiga / acumulación cíclica del daño"
  - "tangente consistente para IsotropicDamage1D — fuera del alcance de esta spec; misma deuda gemela"

numerical_caveats:
  - "Tangente consistente NO simétrica → IS_SYMMETRIC=False obliga al despachador algebraico a usar LU (ADR 0003). Costo por paso ~30-50% mayor que Cholesky pero convergencia cuadrática del Newton global lo compensa."
  - "Transición loading↔unloading es discontinua en la tangente. Newton converge antes de oscilar entre ramas en casos típicos. Para problemas patológicos cerca del umbral, considerar arc-length."
  - "Saturación d=DAMAGE_MAX corta la corrección consistente (∂d/∂κ deja de ser efectivamente (1-d)(1/κ+α)); el algoritmo devuelve tangente secante en ese caso para evitar términos espurios."
  - "Sin regularización, la solución es mesh-dependent en régimen de ablandamiento (localización en una banda de elementos). Para benchmark cuantitativo se requiere mismo refinamiento de malla; para comparaciones cualitativas (presencia de banda, dirección) basta con regularización de gradiente o longitud característica (fuera de scope)."

acceptance:
  # Cubierto por tests/test_materials_unit.py · TestIsotropicDamage2D (preexistentes)
  # y tests/test_damage2d_consistent_tangent.py (nuevos)
  unit_existing_no_regression:
    - name: secant_response_below_threshold
      expect: "κ ≤ κ_0 ⇒ d=0, σ = C_e·ε, tangente = C_e (sin cambio respecto a versión secante anterior)"
      tol_rel: 1.0e-12

    - name: damage_evolution_exponential
      expect: "para κ > κ_0, d cumple ley exponencial dentro de la fórmula"
      tol_rel: 1.0e-10

  unit_consistent_tangent:
    - name: tangent_equals_secant_on_unloading
      setup: "carga hasta κ>κ_0, descarga con ε de menor norma"
      expect: "C_tan = (1-d)·C_e exacto (segundo término se anula con ∂κ/∂ε = 0)"
      tol_rel: 1.0e-12

    - name: tangent_equals_secant_below_threshold
      setup: "ε pequeño, κ_old = κ_0"
      expect: "C_tan = C_e (sin daño activo, ni secante con d=0 ni corrección)"
      tol_rel: 1.0e-12

    - name: tangent_equals_secant_at_saturation
      setup: "carga hasta κ enorme tal que d alcanza DAMAGE_MAX"
      expect: "C_tan = (1-DAMAGE_MAX)·C_e (el término consistente se corta porque ∂d/∂κ efectivamente 0 al saturar)"
      tol_rel: 1.0e-12

    - name: tangent_consistent_term_on_loading
      setup: "carga activa con κ > κ_0 y d < DAMAGE_MAX"
      expect: "C_tan ≠ (1-d)·C_e — incluye el término rank-1 -(1-d)(1/κ+α)/ε_eq · (Ce·ε) ⊗ (M·ε)"
      verification: "compara contra fórmula explícita evaluada en el mismo estado"
      tol_rel: 1.0e-9

    - name: tangent_not_symmetric_in_loading
      setup: "carga activa con estado de deformación no degenerado (ε_xx ≠ ε_yy, γ_xy ≠ 0)"
      expect: "C_tan − C_tan^T ≠ 0 (al menos un elemento off-diagonal con diferencia significativa)"
      verification: "‖C_tan - C_tan^T‖_F / ‖C_tan‖_F ≥ 0.01 en el caso de prueba"

    - name: tangent_finite_difference_consistency
      setup: "carga activa; calcular tangente por diferencia finita centrada (perturbación 1e-7) y comparar con cerrada"
      expect: "‖C_tan_FD − C_tan_closed_form‖ / ‖C_tan‖ < 1e-5"
      tol_rel: 1.0e-5

  unit_arch:
    - name: is_symmetric_attribute_is_false
      expect: "IsotropicDamage2D.IS_SYMMETRIC == False (atributo declarado en clase)"

  integration:
    - name: quad4_damage_newton_converges_in_few_iter
      setup: "Quad4 unitario, IsotropicDamage2D plane stress, carga uniaxial hasta régimen plenamente dañado"
      expect: "Newton converge en ≤6 iteraciones por paso (evidencia de convergencia cuadrática con tangente consistente)"

references:
  - "Lemaitre J., Chaboche J.-L. (1990). Mechanics of Solid Materials. Cambridge UP. Cap. 7 (continuum damage mechanics)."
  - "Oliver J. (1989). A consistent characteristic length for smeared cracking models. IJNME 28, 461-474."
  - "Simó J.C., Ju J.W. (1987). Strain- and stress-based continuum damage models — I. Formulation. Int. J. Solids Struct. 23, 821-840 (tangente consistente para daño isótropo)."
  - "de Souza Neto E.A., Perić D., Owen D.R.J. (2008). Computational Methods for Plasticity. Wiley. §12 (damage models)."
```

---

## Implementación

- **Archivo**: [fenix/materials/damage_2d.py](../../fenix/materials/damage_2d.py).
- **Clase**: `IsotropicDamage2D`, registrada vía `@MaterialRegistry.register`. Declara `IS_SYMMETRIC = False` como ClassVar (sobreescribe el default `True` de `Material`).
- **Estado**: dict `{'kappa': float, 'damage': float}`. `PRIMARY_STATE_VAR = 'damage'`.
- **Algoritmo `compute_state`**:
  1. Recuperar $\kappa_n$ y calcular $\varepsilon_\text{eq}(\boldsymbol\varepsilon)$.
  2. Si $\varepsilon_\text{eq} > \kappa_n$ y $\kappa_{n+1} > \kappa_0$: rama de carga activa.
  3. Calcular $d$ por la ley exponencial, saturar a `DAMAGE_MAX`.
  4. Calcular $\boldsymbol\sigma = (1-d)\mathbf C_e\boldsymbol\varepsilon$.
  5. Tangente: rama según §7 (secante o consistente). En la rama consistente, calcular el término rank-1 una sola vez y restar a $(1-d)\mathbf C_e$.
- **`admissibility_scale`**: sin cambios respecto a versión previa, `E·κ_0` (umbral inicial — la escala no evoluciona con el estado porque el criterio se mide contra $\kappa_0$, no contra $\kappa$).
- **Compatibilidad**: el cambio rompe `IS_SYMMETRIC` para este material. El despachador algebraico (ADR 0003) verifica `domain_is_symmetric` agregando `material.IS_SYMMETRIC` sobre todos los materiales del dominio; basta un material con `IS_SYMMETRIC=False` para que el solver elija LU.
- **Tests**:
  - [tests/test_materials_unit.py](../../tests/test_materials_unit.py) · `TestIsotropicDamage2D` (existente, no regresión) + nueva clase `TestIsotropicDamage2DConsistentTangent`.
  - [tests/test_solid_2d_damage.py](../../tests/test_solid_2d_damage.py) nuevo: integración Quad4 + IsotropicDamage2D + NonlinearSolver con tangente consistente.

---

## Diálogo

- **2026-05-14** · Spec creada retroactivamente sobre material existente (tangente secante) y extendida con tangente algorítmica consistente. Sin ADR: extiende un material existente añadiendo un término a la tangente; no rompe contratos arquitecturales transversales. Sí cambia `IS_SYMMETRIC` para este material, lo cual el despachador algebraico ya soporta vía agregación sobre el dominio (no requiere refactor).
- **2026-05-14** · Decisión: la transición loading↔unloading se decide con $\kappa_n$ committed (no con $\kappa$ trial dentro del Newton). Esto es lo estándar y evita oscilaciones del Newton entre ramas en pasos con ambigüedad. Si el paso converge, $\kappa$ se compromete y la decisión es consistente con la trayectoria.
- **2026-05-14** · `IsotropicDamage1D` queda con tangente secante. Misma deuda gemela conceptual, pero fuera del alcance pactado (A2 era solo 2D). Si el usuario la prioriza luego, replicar es trivial: misma fórmula reducida a escalar.
