# ISOTROPIC DAMAGE 1D — daño escalar 1D con ablandamiento exponencial

> Orden de trabajo. El usuario escribe **especificación física**, **formulación** y **contrato**; la IA rellena **implementación** y responde en **diálogo**.
>
> Spec **retroactiva con ampliación**: modelo implementado desde sesiones anteriores con tangente secante; esta spec lo documenta y añade la tangente algorítmica consistente, replicando la derivación de [[IsotropicDamage2D]] reducida a 1D escalar.

---

## Especificación física

Versión 1D de [`IsotropicDamage2D`](IsotropicDamage2D.md). Modelo de daño escalar isótropo en pequeñas deformaciones para barras, cables y elementos axiales. Captura degradación progresiva de la rigidez axial con ablandamiento exponencial.

### Ecuaciones

- **Constitutiva**: $\sigma = (1 - d)\,E\,\varepsilon$, con $\sigma^\text{eff} = E\varepsilon$.
- **Deformación equivalente**: $\varepsilon_\text{eq} = |\varepsilon|$.
- **Variable histórica**: $\kappa_{n+1} = \max(\kappa_n, |\varepsilon_{n+1}|)$, con $\kappa_0$ inicial.
- **Daño**: $d(\kappa) = \begin{cases}0 & \kappa \leq \kappa_0 \\ 1 - (\kappa_0/\kappa)\,e^{-\alpha(\kappa - \kappa_0)} & \kappa > \kappa_0\end{cases}$, saturado a `DAMAGE_MAX`.

Todas las propiedades físicas (irreversibilidad, asintótica $d \to 1$, papel de $\alpha$) son las del modelo 2D restringidas a un escalar. No distingue tracción de compresión.

---

## Formulación numérica

### Tangente algorítmica consistente

$$\frac{\partial \sigma}{\partial \varepsilon} = (1-d)\,E - E\varepsilon\,\frac{\partial d}{\partial \kappa}\,\frac{\partial \kappa}{\partial \varepsilon}$$

con:

- $\dfrac{\partial d}{\partial \kappa} = (1-d)\left(\dfrac{1}{\kappa} + \alpha\right)$ (idéntico a 2D).
- $\dfrac{\partial \kappa}{\partial \varepsilon} = \operatorname{sign}(\varepsilon)$ en **carga activa** ($|\varepsilon| > \kappa_n$); $0$ en descarga.

**Tangente final**:

$$E^\text{alg} = \begin{cases}
E & \kappa \leq \kappa_0 \;(\text{sin daño}) \\
(1-d)\,E & \text{descarga} \\
(1-d)\,E\big[1 - (1/\kappa + \alpha)\,\kappa\big] = -(1-d)\,E\,\alpha\,\kappa & \text{carga activa, } d < d_\text{max} \\
(1-d_\text{max})\,E & d = d_\text{max} \;(\text{saturación})
\end{cases}$$

donde se ha usado $(1/\kappa + \alpha)\kappa = 1 + \alpha\kappa$ y $E\varepsilon\,\operatorname{sign}(\varepsilon) = E|\varepsilon| = E\kappa$ en carga activa.

### Signo de la tangente consistente

En carga activa **siempre es negativa**: $E^\text{alg} = -(1-d)\,E\,\alpha\,\kappa < 0$ porque $1-d > 0$, $E > 0$, $\alpha > 0$, $\kappa > 0$. Refleja físicamente el régimen de **softening** (ablandamiento): la curva $\sigma$–$\varepsilon$ desciende tras el umbral.

**Implicaciones numéricas**:

- La rigidez tangente del elemento ($\mathbf K_e = \int B^\top E^\text{alg} B\, dV$) puede tener autovalores negativos. La matriz **sigue siendo simétrica** (escalar × forma cuadrática) — `IS_SYMMETRIC = True` se mantiene.
- La matriz global puede dejar de ser **positiva definida**. El despachador algebraico (ADR 0003) detecta el fallo de Cholesky y degrada a LU automáticamente — no requiere cambio declarativo.
- Control de carga puede no converger en problemas con suficiente daño para producir snap-back; en ese régimen se requiere control de longitud de arco (`ArcLengthSolver`).

### Por qué no tests de integración en esta spec

A diferencia de [[IsotropicDamage2D]], no se añade benchmark de integración Truss2D + IsotropicDamage1D + NonlinearSolver: con un solo elemento axial el comportamiento global está dominado por la inestabilidad de softening ($E_{\text{alg}} < 0$) y la única vía limpia de seguir la respuesta post-pico es arc-length. La validación de la tangente consistente se cubre con tests unitarios (incluyendo diferencia finita centrada como verificación de la derivada cerrada). Si en el futuro aparece un caso real de problema con daño 1D ramificado, se añadirá un benchmark dedicado con `ArcLengthSolver`.

---

## Contrato de implementación

```yaml
name: IsotropicDamage1D
kind: material
status: validated

interface:
  strain_dim: 1
  primary_state_var: damage
  is_symmetric: true       # tangente escalar ⇒ contribución elemento simétrica (puede ser indefinida pero simétrica)

parameters:
  - { name: E,        type: float, required: true,  desc: "Módulo de Young intacto" }
  - { name: kappa_0,  type: float, required: true,  desc: "Umbral de deformación equivalente |ε|" }
  - { name: alpha,    type: float, required: true,  desc: "Velocidad de degradación exponencial (>0)" }
  - { name: density,  type: float, required: false, default: null, desc: "Densidad (ADR 0008)" }

signature:
  compute_state: "(ε: float, state_vars=None) -> (σ: float, E_tan: float, state_vars')"
  strain_kind: axial escalar

state_schema:
  kappa:  "float ≥ κ_0 — máxima |ε| histórica"
  damage: "float ∈ [0, DAMAGE_MAX]"

conventions:
  sign: "σ > 0 ⇔ tracción; el modelo no distingue tracción de compresión en la activación del daño"
  damage_irreversibility: "κ monótono no decreciente; d(κ) no decreciente"

validity:
  - "E > 0, κ_0 > 0, α > 0"
  - "|ε| ≲ 1e-2 (régimen de pequeñas deformaciones)"

out_of_scope:
  - "split tracción/compresión"
  - "fatiga / acumulación cíclica"
  - "regularización para mesh-dependency en problemas con varios elementos en serie"
  - "test de integración con NonlinearSolver — requiere arc-length para superar el snap-back, fuera de scope hasta caso real"

numerical_caveats:
  - "Tangente consistente NEGATIVA en carga activa (softening). La rigidez global puede no ser PD; el despachador algebraico (ADR 0003) degrada a LU automáticamente al detectar fallo de Cholesky."
  - "Control de carga con un elemento aislado tras alcanzar el pico se vuelve inestable. Usar arc-length para seguir la rama post-pico."
  - "Transición carga↔descarga discontinua en la tangente (rama loading vs unloading): aceptable para Newton convergente, problemático cerca del umbral en problemas patológicos."

acceptance:
  unit_existing_no_regression:
    - name: secant_response_below_threshold
      expect: "ε pequeño, κ ≤ κ_0 ⇒ d=0, σ=E·ε, tangente = E"
      tol_rel: 1.0e-12

    - name: damage_evolution_exponential
      expect: "ley exponencial reproducida con la fórmula"
      tol_rel: 1.0e-10

  unit_consistent_tangent:
    - name: tangent_equals_E_below_threshold
      setup: "ε < κ_0"
      expect: "E_tan = E (sin daño)"
      tol_rel: 1.0e-12

    - name: tangent_equals_secant_on_unloading
      setup: "cargar a κ > κ_0; descargar a |ε| < κ"
      expect: "E_tan = (1-d)·E exacto"
      tol_rel: 1.0e-12

    - name: tangent_negative_on_loading
      setup: "carga activa con κ > κ_0, d ∈ (0, DAMAGE_MAX)"
      expect: "E_tan = -(1-d)·E·α·κ < 0"
      tol_rel: 1.0e-12

    - name: tangent_equals_secant_at_saturation
      setup: "ε enorme tal que d = DAMAGE_MAX"
      expect: "E_tan = (1-DAMAGE_MAX)·E (sin término consistente al saturar)"
      tol_rel: 1.0e-12

    - name: tangent_matches_finite_difference
      setup: "carga activa; ∂σ/∂ε por FD centrada (h=1e-7) vs fórmula cerrada"
      expect: "error relativo < 1e-5"
      tol_rel: 1.0e-5

    - name: tangent_consistent_with_2D_when_uniaxial
      setup: "evaluar IsotropicDamage1D(E, κ_0, α) en ε=eps y comparar contra IsotropicDamage2D(E, ν=0, κ_0, α) plane stress en ε=[eps, 0, 0]; la primera componente de σ y la entrada [0,0] de C_alg deben coincidir"
      expect: "σ_1D = σ_xx_2D; E_tan_1D = C_alg_2D[0,0]"
      tol_rel: 1.0e-10

references:
  - "Lemaitre J., Chaboche J.-L. (1990). Mechanics of Solid Materials. Cambridge UP. Cap. 7."
  - "Simó J.C., Ju J.W. (1987). Strain- and stress-based continuum damage models — I. Int. J. Solids Struct. 23, 821-840."
  - "Spec hermana: [docs/specs/IsotropicDamage2D.md](IsotropicDamage2D.md) — modelo 2D, derivación detallada de la tangente consistente."
```

---

## Implementación

- **Archivo**: [fenix/materials/damage_1d.py](../../fenix/materials/damage_1d.py).
- **Clase**: `IsotropicDamage1D`, registrada vía `@MaterialRegistry.register`. `IS_SYMMETRIC` queda `True` (default heredado de `Material`) — la contribución elemental es simétrica aunque pueda ser indefinida.
- **Algoritmo `compute_state`**:
  1. Recuperar $\kappa_n$ y calcular $\varepsilon_\text{eq} = |\varepsilon|$.
  2. Si $\varepsilon_\text{eq} > \kappa_n$: carga activa, $\kappa_{n+1} = \varepsilon_\text{eq}$. Si no: descarga, $\kappa_{n+1} = \kappa_n$.
  3. Calcular $d$ por la ley exponencial, saturar a `DAMAGE_MAX`.
  4. $\sigma = (1 - d)\,E\,\varepsilon$.
  5. Tangente: ramas (sin daño / descarga / saturación / carga activa) idénticas en estructura al modelo 2D, con la fórmula escalar resultante. En carga activa con daño efectivo no saturado, `E_tan = -(1-d)·E·α·κ`.
- **Tests**:
  - [tests/test_materials_unit.py](../../tests/test_materials_unit.py) · `TestIsotropicDamage1D` (existente, no regresión) + nueva clase `TestIsotropicDamage1DConsistentTangent`.
  - Sin tests de integración (justificado arriba).

---

## Diálogo

- **2026-05-14** · Spec creada retroactivamente sobre material existente (tangente secante) y extendida con tangente algorítmica consistente. Réplica directa de [[IsotropicDamage2D]] reducida a escalar. Sin ADR — extensión local al material, no afecta a contratos.
- **2026-05-14** · `IS_SYMMETRIC` permanece `True`: la rigidez axial es escalar y la contribución elemental `B^T·E_tan·B` es simétrica aunque pueda ser indefinida cuando `E_tan < 0`. El sistema global puede dejar de ser PD; el fallback automático de Cholesky a LU del despachador algebraico (ADR 0003) lo gestiona sin cambios declarativos.
- **2026-05-14** · No se añade test de integración. Justificación: con tangente consistente <0 en carga activa, un benchmark Truss2D bajo control de carga se vuelve inestable tras el pico; el seguimiento de la rama post-pico requiere arc-length que cae fuera del scope inmediato. Tests unitarios (incluyendo FD numérica y consistencia con la versión 2D reducida) cubren la corrección de la fórmula.
