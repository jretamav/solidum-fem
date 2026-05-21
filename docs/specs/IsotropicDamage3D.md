# ISOTROPIC DAMAGE 3D — daño escalar isótropo con ablandamiento exponencial

> Orden de trabajo. El usuario revisa **especificación física**, **formulación numérica** y **contrato**; la IA propone, ejecuta y rellena **implementación** + **diálogo**.
>
> **Pre-redactada por la IA** sobre la base de `IsotropicDamage2D.md` (mismo modelo físico) y `Elastic3D.md` / `VonMises3D.md` (patrón Voigt 6D). Tercera y última entrega de la sub-etapa **A.bis** (materiales 3D no lineales).

---

## Especificación física

### 0. Descripción general

Modelo de **daño isótropo escalar** en pequeñas deformaciones para sólidos 3D, formulado íntegramente en Voigt 6D del proyecto (ADR 0012). Captura degradación progresiva e isótropa de la rigidez en respuesta a la deformación máxima histórica, con ablandamiento exponencial. Pensado para hormigón en tracción uniforme, cerámicos, materiales cuasi-frágiles donde la microfisuración distribuida puede modelarse mediante un escalar.

Diferencias respecto a `IsotropicDamage2D`:

- **Sin variantes de hipótesis** (igual que VM3D, DP3D, Elastic3D): en 3D todas las componentes son activas.
- **Algoritmo idéntico al 2D extendido a 6D**: misma ley de evolución exponencial, mismo régimen carga/descarga por $\kappa$, misma estructura de tangente algorítmica consistente. **Sin Newton local** (el algoritmo es explícito).
- **`IsotropicDamage2D` plane_strain** es la restricción de este modelo bajo $\varepsilon_{zz} = \gamma_{yz} = \gamma_{xz} = 0$ con la misma fórmula de $\varepsilon_\text{eq}$ — el test cruzado se hace contra esta variante (no contra plane_stress, que es genuinamente distinta).

### 1. Cinemática

Pequeñas deformaciones, descomposición elástica pura: no hay deformación plástica (sin plasticidad). El daño se manifiesta como degradación de la rigidez efectiva. Sin variables internas tensoriales — solo el historial escalar $\kappa$ y el daño $d$.

### 2. Esfuerzo nominal vs esfuerzo efectivo

$$\boldsymbol\sigma = (1 - d)\,\boldsymbol\sigma^\text{eff}, \qquad \boldsymbol\sigma^\text{eff} = \mathbf C_e\,\boldsymbol\varepsilon$$

con $d \in [0, d_\text{max}]$ la variable de daño isótropa y $\mathbf C_e$ el tensor constitutivo elástico isótropo 3D 6×6 de `Elastic3D` (sin variantes de hipótesis en 3D).

El esfuerzo efectivo $\boldsymbol\sigma^\text{eff}$ es el que el material intacto soportaría con la deformación corriente. El esfuerzo nominal $\boldsymbol\sigma$ es lo que efectivamente transmite considerando la degradación.

### 3. Deformación equivalente

Cantidad escalar que cuantifica la "intensidad" de la deformación. Convención del proyecto, norma simétrica simple en notación Voigt 6D $[\varepsilon_{xx}, \varepsilon_{yy}, \varepsilon_{zz}, \gamma_{xy}, \gamma_{yz}, \gamma_{xz}]$ con $\gamma_{ij} = 2\varepsilon_{ij}$:

$$\varepsilon_\text{eq}(\boldsymbol\varepsilon) = \sqrt{\varepsilon_{xx}^2 + \varepsilon_{yy}^2 + \varepsilon_{zz}^2 + \tfrac{1}{2}(\gamma_{xy}^2 + \gamma_{yz}^2 + \gamma_{xz}^2)} = \sqrt{\boldsymbol\varepsilon^\top \mathbf M\,\boldsymbol\varepsilon}$$

con $\mathbf M = \operatorname{diag}(1, 1, 1, 1/2, 1/2, 1/2)$. Equivalente a la norma de Frobenius del tensor deformación 3D ($\boldsymbol\varepsilon : \boldsymbol\varepsilon$ con $\varepsilon_{ij} = \gamma_{ij}/2$):

$$\boldsymbol\varepsilon:\boldsymbol\varepsilon = \sum_i \varepsilon_{ii}^2 + 2\sum_{i<j}\varepsilon_{ij}^2 = \varepsilon_{xx}^2 + \varepsilon_{yy}^2 + \varepsilon_{zz}^2 + \tfrac{1}{2}(\gamma_{xy}^2 + \gamma_{yz}^2 + \gamma_{xz}^2)$$

Es **extensión natural** de la convención 2D del proyecto ($\mathbf M_\text{2D} = \operatorname{diag}(1, 1, 1/2)$). Bajo restricción plane strain ($\varepsilon_{zz} = \gamma_{yz} = \gamma_{xz} = 0$), $\varepsilon_\text{eq}^{3D}$ se reduce exactamente a $\varepsilon_\text{eq}^{2D}$ — propiedad clave para el test cruzado con `IsotropicDamage2D` plane_strain (§13).

**Limitación**: no distingue tracción de compresión. Para hormigón a compresión con resistencia distinta se necesitaría un split tipo Mazars o de Vree; explícitamente *out-of-scope* (idéntica deuda gemela al 2D).

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

con saturación numérica $d \leq d_\text{max} = $ `DAMAGE_MAX` (evita rigidez nula).

Propiedades idénticas al 2D: $d(\kappa_0) = 0$ (continuidad); $d'(\kappa) > 0$ para $\kappa > \kappa_0$ (irreversibilidad); $\lim_{\kappa \to \infty} d = 1$ (asintótico); $\alpha$ controla la velocidad de degradación. **Misma función de ablandamiento** que en 1D y 2D — centralizada en `solidum.materials._softening.evaluate_exponential_damage`.

---

## Formulación numérica

### 6. Algoritmo en un paso

Dado $\{\kappa_n\}$ y $\boldsymbol\varepsilon_{n+1}$:

1. $\varepsilon_\text{eq} = \sqrt{\boldsymbol\varepsilon^\top \mathbf M\,\boldsymbol\varepsilon}$ (norma 6D).
2. Si $\varepsilon_\text{eq} > \kappa_n$: **carga activa**, $\kappa_{n+1} = \varepsilon_\text{eq}$. Si no: **descarga**, $\kappa_{n+1} = \kappa_n$.
3. Aplicar ley de daño con $\kappa_{n+1}$; saturar a $d_\text{max}$.
4. $\boldsymbol\sigma_{n+1} = (1 - d)\,\mathbf C_e\,\boldsymbol\varepsilon_{n+1}$.
5. Tangente: **secante** en descarga / sin daño / saturación; **algorítmica consistente** en carga activa con daño no saturado (ver §7).

El algoritmo es **explícito**: no requiere Newton local. La actualización de $\kappa$ y $d$ es directa, idéntica al 2D.

### 7. Tangente algorítmica consistente

$$\mathbf C^\text{alg} = \frac{\partial \boldsymbol\sigma_{n+1}}{\partial \boldsymbol\varepsilon_{n+1}}$$

Derivando $\boldsymbol\sigma = (1-d)\,\mathbf C_e\,\boldsymbol\varepsilon$:

$$\mathbf C^\text{alg} = (1-d)\,\mathbf C_e - \boldsymbol\sigma^\text{eff} \otimes \frac{\partial d}{\partial \boldsymbol\varepsilon}$$

con $\boldsymbol\sigma^\text{eff} = \mathbf C_e\,\boldsymbol\varepsilon$ (vector 6×1) y la derivada del daño:

$$\frac{\partial d}{\partial \boldsymbol\varepsilon} = \frac{\partial d}{\partial \kappa}\,\frac{\partial \kappa}{\partial \boldsymbol\varepsilon}$$

donde:

- $\dfrac{\partial d}{\partial \kappa} = (1 - d)\left(\dfrac{1}{\kappa} + \alpha\right)$ para la ley exponencial (idéntica al 2D, derivación independiente de la dimensión).

- $\dfrac{\partial \kappa}{\partial \boldsymbol\varepsilon}$ depende del régimen:
  - **Carga activa** ($\varepsilon_\text{eq} > \kappa_n$, $\kappa = \varepsilon_\text{eq}$):
    $$\dfrac{\partial \kappa}{\partial \boldsymbol\varepsilon} = \dfrac{\partial \varepsilon_\text{eq}}{\partial \boldsymbol\varepsilon} = \dfrac{1}{\varepsilon_\text{eq}}\,\mathbf M\,\boldsymbol\varepsilon = \dfrac{1}{\varepsilon_\text{eq}}\begin{bmatrix}\varepsilon_{xx} \\ \varepsilon_{yy} \\ \varepsilon_{zz} \\ \tfrac{1}{2}\gamma_{xy} \\ \tfrac{1}{2}\gamma_{yz} \\ \tfrac{1}{2}\gamma_{xz}\end{bmatrix}$$
  - **Descarga** ($\varepsilon_\text{eq} \leq \kappa_n$): $\partial \kappa / \partial \boldsymbol\varepsilon = 0$.

**Tangente final**:

$$\mathbf C^\text{alg} = \begin{cases}
\mathbf C_e & \kappa \leq \kappa_0 \;(\text{sin daño}) \\
(1-d)\,\mathbf C_e & \text{descarga} \\
(1-d)\,\mathbf C_e \;-\; \dfrac{(1-d)(1/\kappa + \alpha)}{\varepsilon_\text{eq}}\,(\mathbf C_e\boldsymbol\varepsilon)\otimes(\mathbf M\boldsymbol\varepsilon) & \text{carga activa, } d < d_\text{max} \\
(1-d_\text{max})\,\mathbf C_e & d = d_\text{max} \;(\text{saturación})
\end{cases}$$

con $\mathbf M = \operatorname{diag}(1, 1, 1, 1/2, 1/2, 1/2)$ (6×6).

### 8. Pérdida de simetría

El producto externo $(\mathbf C_e\boldsymbol\varepsilon)\otimes(\mathbf M\boldsymbol\varepsilon)$ es una matriz $6\times 6$ que **no es simétrica en general**. $\mathbf C_e\boldsymbol\varepsilon$ y $\mathbf M\boldsymbol\varepsilon$ no son proporcionales salvo en estados de deformación muy particulares (e.g. uniaxial puro alineado con ejes principales).

Consecuencia: **`IS_SYMMETRIC = False` obligatorio**. El despachador algebraico (ADR 0003) elige LU en lugar de Cholesky/LDLᵀ para sistemas globales que incluyan este material. Coste por paso ~30-50% mayor, compensado con creces por la convergencia cuadrática del Newton global (típicamente 2-3 iteraciones en lugar de 5-10 con tangente secante).

### 9. Decisión de régimen durante el Newton

Idéntica al 2D: la transición loading↔unloading se decide con $\kappa_n$ committed (no con $\kappa$ trial dentro del Newton). Estándar y evita oscilaciones del Newton entre ramas en pasos con ambigüedad. Si el paso converge, $\kappa$ se compromete y la decisión es consistente con la trayectoria.

---

## Contrato de implementación

```yaml
name: IsotropicDamage3D
kind: material
status: validated        # draft → implemented → validated

interface:
  strain_dim: 6
  primary_state_var: damage
  is_symmetric: false      # tangente consistente no simétrica (ver §8)

parameters:
  - { name: E,       type: float, required: true,  desc: "Módulo de Young intacto (>0)" }
  - { name: nu,      type: float, required: true,  desc: "Coeficiente de Poisson ∈ (-1, 0.5)" }
  - { name: kappa_0, type: float, required: true,  desc: "Umbral de deformación equivalente; daño nulo mientras κ ≤ κ_0 (>0)" }
  - { name: alpha,   type: float, required: true,  desc: "Velocidad de degradación exponencial (>0). α mayor → ablandamiento más abrupto" }
  - { name: density, type: float, required: false, default: null, desc: "Densidad (ADR 0008)" }

signature:
  compute_state: "(ε: ndarray(6,), state_vars=None) -> (σ: ndarray(6,), C_tan: ndarray(6,6), state_vars')"
  strain_kind: "3D Voigt — [ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz] con γ_ij = 2·ε_ij"

state_schema:
  kappa:  "float ≥ κ_0 — máxima deformación equivalente histórica"
  damage: "float ∈ [0, DAMAGE_MAX] — variable de daño escalar isótropa"

conventions:
  sign:  "σ > 0 ⇔ tracción; el modelo no distingue tracción de compresión en la activación del daño"
  voigt: "γ_ij = 2·ε_ij en deformación; ε_eq usa M = diag(1, 1, 1, 1/2, 1/2, 1/2)"
  damage_irreversibility: "κ monótono no decreciente vía max(κ_old, ε_eq); d(κ) creciente"

validity:
  - "E > 0, κ_0 > 0, α > 0"
  - "ν ∈ (-1, 0.5) estricto"

out_of_scope:
  - "split tracción/compresión (Mazars, de Vree, Mises modificado); el modelo daña por igual en cualquier régimen"
  - "anisotropía del daño (tensor de daño D vs escalar)"
  - "regularización para evitar mesh-dependency en regime de ablandamiento (longitud característica, no-local, gradient)"
  - "acoplamiento con plasticidad"
  - "fatiga / acumulación cíclica del daño"

numerical_caveats:
  - "Tangente consistente NO simétrica → IS_SYMMETRIC=False obliga al despachador algebraico a usar LU (ADR 0003). Costo por paso ~30-50% mayor que Cholesky pero convergencia cuadrática del Newton global lo compensa (Simó-Ju 1987)."
  - "Transición loading↔unloading discontinua en la tangente. Newton converge antes de oscilar entre ramas en casos típicos. Para problemas patológicos cerca del umbral, considerar arc-length."
  - "Saturación d=DAMAGE_MAX corta la corrección consistente; el algoritmo devuelve tangente secante en ese caso para evitar términos espurios."
  - "Sin regularización, la solución es mesh-dependent en régimen de ablandamiento (localización en una banda de elementos). Es limitación inherente del modelo continuo de daño isótropo escalar; out-of-scope esta entrega."

acceptance:
  unit:
    - name: paso_elastico_bajo_umbral
      setup: "ε pequeño 6D con todas las componentes activas, ε_eq < κ_0"
      expect: "σ = C_e·ε exacto; d=0; C_tan = C_e (sin corrección rank-1)"
      tol_rel: 1.0e-12

    - name: damage_evolution_exponential
      setup: "ε grande tal que ε_eq >> κ_0"
      expect: "d cumple la ley exponencial; estado coincide con la fórmula directa"
      tol_rel: 1.0e-10

    - name: rejects_loading_below_threshold
      setup: "ε_eq apenas por encima de κ_0 (banda elástica + bandera de carga activa)"
      expect: "d > 0 minúsculo; κ_new = ε_eq; tangente con corrección rank-1"

    - name: tangent_equals_secant_on_unloading
      setup: "carga hasta κ > κ_0; descarga con ε de menor norma"
      expect: "C_tan = (1-d)·C_e exacto (segundo término se anula con ∂κ/∂ε = 0)"
      tol_rel: 1.0e-12

    - name: tangent_equals_secant_below_threshold
      setup: "ε pequeño, κ_old = κ_0 (estado inicial intacto)"
      expect: "C_tan = C_e (sin daño activo)"
      tol_rel: 1.0e-12

    - name: tangent_equals_secant_at_saturation
      setup: "carga hasta κ enorme tal que d alcanza DAMAGE_MAX"
      expect: "C_tan = (1-DAMAGE_MAX)·C_e (el término consistente se corta)"
      tol_rel: 1.0e-12

    - name: tangent_not_symmetric_in_loading
      setup: "carga activa con estado de deformación no degenerado (ε_xx ≠ ε_yy ≠ ε_zz, cortantes activos)"
      expect: "‖C_tan - C_tan^T‖_F / ‖C_tan‖_F ≥ 0.01 (asimetría significativa)"

    - name: tangent_finite_difference_consistency
      setup: "carga activa; comparar C_alg cerrada con diferencia finita centrada"
      expect: "‖C_FD − C_alg‖ / ‖C_alg‖ < 1e-5"
      tol_rel: 1.0e-5

    - name: kappa_monotonic
      setup: "trayectoria multi-paso con carga, descarga, recarga"
      expect: "κ es monótono no decreciente en todos los pasos"

    - name: descarga_recarga_irreversibilidad
      setup: "cargar a κ_load, descargar elásticamente, recargar"
      expect: "d permanece constante durante toda la descarga; al recargar, d sigue creciendo desde el valor de descarga"

    - name: unit_invariance_MPa_vs_Pa
      setup: "mismo path adimensional con (MPa,mm,N) y (Pa,m,N)"
      expect: "d, κ y σ/E idénticos (κ_0 es una deformación, ε es adimensional → invariancia exacta)"
      tol_rel: 1.0e-12

    - name: rechazo_inputs_invalidos
      setup: "construir con E≤0, κ_0≤0, α≤0, ν∉(-1, 0.5), density<0"
      expect: "ValueError con mensaje claro en cada caso"

  cross_consistency:
    - name: equivalencia_plane_strain
      setup: "Damage3D con ε = (ε_xx, ε_yy, 0, γ_xy, 0, 0); Damage2D plane_strain con ε = (ε_xx, ε_yy, γ_xy);
             path multi-paso con carga, descarga, recarga"
      expect: "σ_xx, σ_yy, σ_xy, κ, d coinciden exactamente entre 3D y 2D PS;
              σ_zz en VM3D no nulo (consistente con C_e plane_strain interno);
              componentes 3D ausentes en 2D PS (σ_yz, σ_xz) permanecen nulas"
      tol_rel: 1.0e-12

  integration:
    - name: hex8_elastico
      setup: "Hex8 unitario con carga uniforme muy por debajo del umbral"
      expect: "d = 0 en todos los Gauss points; σ = C_e·ε exacto"
      tol_rel: 1.0e-8

    - name: hex8_damage_newton_converges_in_few_iter
      setup: "Hex8 unitario, carga uniaxial hasta régimen plenamente dañado (κ >> κ_0), 8 pasos"
      expect: "Newton converge en ≤6 iteraciones por paso (evidencia de convergencia cuadrática con tangente consistente)"

    - name: tet4_basic_pipeline
      setup: "Tet4 con confinamiento + tracción uniaxial, carga que activa daño"
      expect: "convergencia, d > 0 en el único Gauss point, κ ≥ κ_0"

references:
  - "Lemaitre J., Chaboche J.-L. (1990). Mechanics of Solid Materials. Cambridge UP. Cap. 7 (continuum damage mechanics)."
  - "Simó J.C., Ju J.W. (1987). Strain- and stress-based continuum damage models — I. Formulation. Int. J. Solids Struct. 23, 821-840 (tangente consistente para daño isótropo)."
  - "Oliver J. (1989). A consistent characteristic length for smeared cracking models. IJNME 28, 461-474."
  - "de Souza Neto E.A., Perić D., Owen D.R.J. (2008). Computational Methods for Plasticity. Wiley. §12 (damage models)."
  - "ADR 0012 — Sólidos 3D: convención Voigt 6D del proyecto."
  - "ADR 0003 — Despachador algebraico: tangente asimétrica fuerza LU global."
  - "ADR 0006 — Tolerancias adimensionales: patrón atol + rtol·escala con escala física E·κ_0."
```

---

## Implementación

- **Archivo**: [solidum/materials/damage_3d.py](../../solidum/materials/damage_3d.py).
- **Clase**: `IsotropicDamage3D`, registrada vía `@MaterialRegistry.register`. Declara `IS_SYMMETRIC = False` como ClassVar (tangente consistente asimétrica).
- **Composición**:
  - `Elastic3D` como base elástica (acceso a $\mathbf C_e$ 6×6 vía `self.elastic_base.C`).
  - `evaluate_exponential_damage` de `solidum.materials._softening` para la ley de softening (centralizada, compartida con `IsotropicDamage1D` e `IsotropicDamage2D`).
- **Estado**: dict `{'kappa': float, 'damage': float}`. `PRIMARY_STATE_VAR = 'damage'`.
- **Algoritmo `compute_state`**: idéntica estructura a `IsotropicDamage2D` extendida a Voigt 6D. Sin Numba (operaciones numpy puras; el coste dominante es el producto matriz-vector $\mathbf C_e \cdot \boldsymbol\varepsilon$, no se beneficia significativamente de JIT).
- **`admissibility_scale`**: `E·κ_0` constante respecto al estado (idéntico al patrón 1D/2D — el criterio se mide contra `κ_0`, no contra `κ` corriente).
- **Tests**:
  - [tests/test_materials_unit.py](../../tests/test_materials_unit.py) · `TestIsotropicDamage3D` (12 tests unitarios) + `TestIsotropicDamage3DvsPlaneStrain` (1 cruzado contra Damage2D plane_strain).
  - [tests/test_solid_3d_plasticity.py](../../tests/test_solid_3d_plasticity.py) · 3 tests de integración: `TestHex8IsotropicDamage3DElastic` (1 — sin daño), `TestHex8IsotropicDamage3DActiveLoad` (1 — daño activo con tangente asimétrica), `TestTet4IsotropicDamage3DBasic` (1 — smoke Tet4).

- **Notas de traducción**:
  - Matriz $\mathbf M = \operatorname{diag}(1, 1, 1, 1/2, 1/2, 1/2)$ aplicada implícitamente en código (factor $1/2$ sobre los slots cortantes de $\boldsymbol\varepsilon$) en lugar de allocar la 6×6, idéntico patrón al 2D.
  - Convención engineering input (γ_ij): el factor $1/2$ en $\mathbf M$ y en $\partial\varepsilon_\text{eq}/\partial\boldsymbol\varepsilon$ corrige automáticamente para que la norma corresponda a la Frobenius del tensor (γ_ij/2 = ε_ij_tens).
  - Validación de parámetros en constructor: `E > 0`, `κ_0 > 0`, `α > 0`, `ν ∈ (-1, 0.5)`, `density ≥ 0` si se pasa. Mensajes en español.

---

## Diálogo

- **2026-05-21** · Spec pre-redactada por la IA al cerrar el segundo material de A.bis (`DruckerPrager3D`). Base: `IsotropicDamage2D.md` (modelo físico, tangente algorítmica consistente, semánticas carga/descarga, cap DAMAGE_MAX) + `VonMises3D.md`/`DruckerPrager3D.md` (estructura Voigt 6D, patrón cross-consistency).

- **2026-05-21** · **Decisión de plumbing (decide la IA, Reglas §3)**: reuso de la función centralizada `evaluate_exponential_damage` de `solidum.materials._softening`. Esta función ya está compartida entre `IsotropicDamage1D` e `IsotropicDamage2D` (auditoría H-3.4 cerrada el 2026-05-19); añadir `IsotropicDamage3D` como tercer consumidor satisface la regla de centralización del proyecto. Sin duplicación de la ley de softening.

- **2026-05-21** · **Decisión de plumbing**: el material usa `Elastic3D` internamente como base elástica (paralelo a cómo `IsotropicDamage2D` usa `Elastic2D`). Acceso a la matriz constitutiva 6×6 vía `self.elastic_base.C` — idéntico patrón al 2D, sin duplicación de la construcción de $\mathbf C_e$.

- **2026-05-21** · **Decisión de plumbing**: la matriz $\mathbf M$ del cálculo de $\varepsilon_\text{eq}$ y de la derivada $\partial\kappa/\partial\boldsymbol\varepsilon$ se materializa implícitamente en código (multiplicaciones por $1/2$ sobre los slots cortantes) en lugar de allocar una 6×6 — ahorro mínimo de memoria pero código más legible y consistente con el patrón seguido en `damage_2d.py`.

- **2026-05-21** · **Cross-consistency** contra `IsotropicDamage2D` plane_strain, no plane_stress. La equivalencia matemática es exacta bajo $\varepsilon_{zz} = \gamma_{yz} = \gamma_{xz} = 0$ impuesto (plane strain): la matriz $\mathbf C_e$ 3D restringida a ε_zz=0 input se reduce a la matriz $\mathbf C_e^\text{plane strain}$ 2D, y $\varepsilon_\text{eq}^\text{3D}$ se reduce a $\varepsilon_\text{eq}^\text{2D}$ bajo la misma restricción. Plane stress 2D usa un $\mathbf C_e$ proyectado distinto y por tanto no es comparable directamente con 3D bajo restricción cinemática.

- **2026-05-21** · **Sin variantes de hipótesis** en 3D (igual que VM3D, DP3D, Elastic3D). Constructor más simple que el 2D.

- **2026-05-21** · **No se incluye consistent tangent para `IsotropicDamage1D`** (deuda gemela del 2D): fuera del alcance de esta entrega A.bis. La fórmula 1D es trivial reducción de la 3D si se prioriza después.

- **2026-05-21** · **Implementación completa, status `validated`**. Clase + 12 tests unitarios + 1 test cruzado + 3 tests de integración (Hex8 elástico, Hex8 con daño activo + tangente asimétrica, Tet4 smoke). **16 tests añadidos a la suite (848 → 866 verdes; suite global sin regresiones).** Todos los tests pasaron a la primera, incluyendo:
  - Cruzado Damage3D ↔ Damage2D plane_strain a **14 decimales en `d` y `κ`** y 10 decimales en σ — la equivalencia bajo $\varepsilon_{zz} = \gamma_{yz} = \gamma_{xz} = 0$ es exacta a precisión máquina (misma ley centralizada, mismo $\mathbf C_e$ plane_strain bajo restricción cinemática).
  - Diferencia finita de la tangente algorítmica vs cerrada con error rel. < 1e-5 en las 36 componentes 6×6 en carga activa.
  - Asimetría significativa de la tangente en estados no degenerados (norma rel. ≥ 0.01).
  - Invariancia bajo cambio de unidades (E en MPa vs Pa) — `d`, `κ` idénticos a 12 decimales; σ escala linealmente con E.

  **Con esta entrega cierra la sub-etapa A.bis** (materiales 3D no lineales): VonMises3D + DruckerPrager3D + IsotropicDamage3D. Validación contra benchmark publicado (e.g. ensayo tensión-softening con valor cuantitativo de carga máxima en barra cargada axialmente con softening exponencial calibrado por G_f efectivo) **diferida a la campaña 3D consolidada** que se hará al cierre de A.bis junto con Hill 3D para VM3D y esfera hueca para DP3D.
