# Catálogo de Modelos Constitutivos

## Elasticidad Lineal

### `Elastic1D` — elasticidad lineal axial

Ley de Hooke escalar $\sigma = E \cdot \varepsilon$. Sin variables internas. Compatible con todos los elementos 1D (armaduras y marcos).

```yaml
materials:
  - id: 1
    type: Elastic1D
    E: 200.0e9
    density: 7850.0
```

### `Elastic2D` — elasticidad lineal isótropa 2D

Tensor constitutivo isótropo $\boldsymbol\sigma = \mathbf C \cdot \boldsymbol\varepsilon$ en notación Voigt $[xx, yy, xy]$. Sin variables internas. Compatible con `Quad4`, `Tri3`.

```yaml
materials:
  - id: 1
    type: Elastic2D
    E: 200.0e9
    nu: 0.3
    hypothesis: plane_stress
    density: 7850.0
```

La hipótesis `plane_stress` se usa para placas delgadas sin confinamiento normal; `plane_strain` para secciones de presa, túneles o problemas con extensión infinita en el eje longitudinal.

## Plasticidad

### `Elastoplastic1D` — plasticidad J2 axial con endurecimiento isótropo

Plasticidad asociativa con criterio $f = |\sigma_{\text{trial}}| - (\sigma_y + H \alpha) \le 0$. Algoritmo de *return mapping* clásico 1D; tangente algorítmica consistente $E_t = E\,H / (E + H)$.

- **Parámetros**: `E`, `sigma_y` (fluencia inicial), `H` (módulo de endurecimiento; `H=0` = perfectamente plástico).
- **Variables internas**: $\varepsilon_p$ (deformación plástica), $\alpha$ (acumulada equivalente).
- **Compatible con**: armaduras y marcos 1D (en marcos solo se aplica al esfuerzo axial; ver advertencia en el capítulo *Catálogo de Elementos Finitos*).

```yaml
materials:
  - id: 1
    type: Elastoplastic1D
    E: 200.0e9
    sigma_y: 250.0e6
    H: 2.0e9
    density: 7850.0
```

### `VonMises2D` — plasticidad J2 plana con endurecimiento isótropo

Plasticidad asociativa con criterio J2 y endurecimiento isótropo lineal. Soporta **dos hipótesis cinemáticas** mutuamente excluyentes, seleccionadas con el campo `hypothesis` al construir el material:

- `plane_strain` ($\varepsilon_{zz} = 0$): return mapping radial cerrado sobre la parte desviadora 3D extendida (Simó-Hughes §3.3). Corrector $\Delta\gamma = f_{\text{trial}} / (2G + \tfrac{2}{3} H)$ con $N = s_{\text{trial}} / \lVert s_{\text{trial}}\rVert$ y actualización $\alpha_{\text{new}} = \alpha + \sqrt{2/3}\,\Delta\gamma$. Default histórico.
- `plane_stress` ($\sigma_{zz} = 0$): *plane stress projected algorithm* (Simó-Hughes §3.4.1). Función de fluencia proyectada $\bar f = \tfrac{1}{2}\boldsymbol\sigma^\top \mathbf P\,\boldsymbol\sigma - R^2/3$ y Newton local escalar sobre $\Delta\gamma$. Converge en 3–6 iteraciones por evaluación. Actualización físicamente correcta $\alpha_{\text{new}} = \alpha + \Delta\gamma\sqrt{2\boldsymbol\sigma^\top\!\mathbf P\,\boldsymbol\sigma / 3}$ (en plane stress $\Delta\gamma$ tiene unidades $[1/\text{esfuerzo}]$, por lo que la heurística $\sqrt{2/3}\,\Delta\gamma$ válida en plane strain rompe la invariancia bajo cambio de unidades). Incompresibilidad plástica cierra $e^p_{zz} = -(e^p_{xx} + e^p_{yy})$.

**Parámetros comunes**: `E`, `nu`, `sigma_y`, `H` ($\ge 0$, default 0 = perfectamente plástico), `hypothesis` (default `plane_strain`), `density` (opcional). Variables internas: $\varepsilon_p$ tensorial 4 componentes $[xx, yy, zz, xy_{\text{tens}}]$ y $\alpha$ acumulada equivalente adimensional. **Tangente algorítmica consistente cerrada** en ambas hipótesis (con corrección por $\mathrm d\alpha/\mathrm d\Delta\gamma \neq \sqrt{2/3}$ en plane stress). Compatible con todos los elementos 2D del catálogo.

```yaml
materials:
  # Hipótesis por defecto: plane_strain
  - id: 2
    type: VonMises2D
    E: 210.0e9
    nu: 0.3
    sigma_y: 250.0e6
    H: 2.0e9
    hypothesis: plane_strain
    density: 7850.0

  # Mismo material en plane_stress (placas delgadas sin confinamiento normal)
  - id: 3
    type: VonMises2D
    E: 210.0e9
    nu: 0.3
    sigma_y: 250.0e6
    H: 2.0e9
    hypothesis: plane_stress
    density: 7850.0
```

### `DruckerPrager2D` — plasticidad friccional cohesivo-friccional 2D

Modelo de Drucker-Prager: cono circular suave de Mohr-Coulomb con cohesión y fricción interna. Criterio de fluencia $f = \sqrt{J_2} + \eta_f\,I_1 - k(\alpha) \le 0$ con $k(\alpha) = k_0 + H\alpha$ (endurecimiento isótropo lineal en cohesión). Plasticidad *no* asociada por defecto: el ángulo de dilatancia $\psi$ es parámetro independiente del ángulo de fricción $\phi$.

**Algoritmo — dos ramas con detección automática**:

- **Return regular** (superficie del cono): $\Delta\gamma = f_{\text{trial}} / (G + 9K\eta_f\eta_g + H)$, dirección desviadora preservada. Tangente algorítmica $\mathbf C_{\text{alg}} = K\,\mathbf v\otimes\mathbf v + 2G(1-\beta)\mathbf I_{\text{dev}} + 4G\beta\,\hat{\mathbf n}\otimes\hat{\mathbf n} - (1/A)\,\mathbf b_g\otimes\mathbf b_f$ con $\beta = G\Delta\gamma / \sqrt{J_2^{\text{trial}}}$.
- **Return al ápice** (vértice puramente hidrostático): si tras el return regular el estado caería fuera del cono, conmuta a $\Delta\gamma_{\text{apex}} = (I_1^{\text{trial}}\eta_f - k(\alpha)) / (9K\eta_f\eta_g + H)$ con $\boldsymbol\sigma_{n+1} = (k/3\eta_f)\,\mathbf I$ puramente hidrostático. Tangente reducida $K H / (9K\eta_f\eta_g + H)\,\mathbf v\otimes\mathbf v$.

**Calibración con Mohr-Coulomb** (campo `variant`):

- `plane_strain_matched` (default): coincide exactamente con MC en plane strain. $\eta_f = \tan\phi / \sqrt{9 + 12\tan^2\phi}$, $k_0 = 3c_0 / \sqrt{9 + 12\tan^2\phi}$.
- `outer_cone`: circunscribe MC. $\eta_f = 2\sin\phi / [\sqrt 3 (3 - \sin\phi)]$.
- `inner_cone`: inscribe MC. $\eta_f = 2\sin\phi / [\sqrt 3 (3 + \sin\phi)]$.

**Parámetros**: `E`, `nu`, `cohesion` (cohesión inicial $c_0$), `phi_deg` ($0 \le \phi < 90$), `psi_deg` (dilatancia, $0 \le \psi \le \phi$; default $\psi = \phi$ = asociada), `H` ($\ge 0$, default 0), `variant` (default `plane_strain_matched`), `density` (opcional). Variables internas: $\varepsilon_p$ tensorial 4 componentes y $\alpha$ acumulada. `IS_SYMMETRIC = False` declarativo (no asociada $\Rightarrow$ tangente asimétrica; el despachador algebraico ADR 0003 elige LU).

```callout Advertencia
**Restricción cinemática**: solo `plane_strain`. La proyección con $\sigma_{zz} = 0$ acoplada a flujo dilatante es notoriamente delicada y queda *out-of-scope* en esta entrega.
```

```yaml
materials:
  - id: 4
    type: DruckerPrager2D
    E: 30.0e9
    nu: 0.25
    cohesion: 2.0e5      # Pa
    phi_deg: 30.0
    psi_deg: 10.0        # no asociada (psi < phi)
    H: 0.0               # plasticidad perfecta
    variant: plane_strain_matched
    density: 2000.0      # arena densa
```

## Daño Mecánico

### `IsotropicDamage1D` — daño isótropo 1D con softening exponencial

Daño escalar $d \in [0, 1)$ con esfuerzo nominal $\sigma = (1 - d)\,E\,\varepsilon$. Evolución por norma de la deformación equivalente $\varepsilon_{eq} = |\varepsilon|$ y umbral $\kappa_0$:

$$d = 1 - \frac{\kappa_0}{\kappa}\,\exp\bigl(-\alpha(\kappa - \kappa_0)\bigr) \quad \text{para } \kappa > \kappa_0$$

- **Parámetros**: `E`, `kappa_0` (umbral elástico), `alpha` (velocidad de softening), `density` (opcional).
- **Variables internas**: $\kappa$ (máxima histórica), $d$ (daño).
- **Tangente**: **algorítmica consistente** en carga activa con daño no saturado: $E_{\text{tan}} = -(1 - d)\,E\,\alpha\,\kappa$ — negativa, refleja la pendiente descendente $\sigma$-$\varepsilon$ post-pico. En descarga, sin daño activo o al saturar, vuelve a la secante $(1-d)E$.

### `IsotropicDamage2D` — daño isótropo 2D

Extensión 2D del modelo anterior. Esfuerzo nominal $\boldsymbol\sigma = (1-d)\,\mathbf C_e\,\boldsymbol\varepsilon$. Deformación equivalente simétrica:

$$\varepsilon_{eq} = \sqrt{\varepsilon_{xx}^2 + \varepsilon_{yy}^2 + \tfrac{1}{2}\gamma_{xy}^2}$$

- **Parámetros**: `E`, `nu`, `kappa_0`, `alpha`, `hypothesis` (`plane_stress` default o `plane_strain`), `density` (opcional).
- **Tangente algorítmica consistente**: $\mathbf C_{\text{alg}} = (1-d)\,\mathbf C_e - \frac{(1-d)(1/\kappa + \alpha)}{\varepsilon_{eq}}\,(\mathbf C_e\boldsymbol\varepsilon) \otimes (\mathbf M\boldsymbol\varepsilon)$ con $\mathbf M = \mathrm{diag}(1, 1, 1/2)$ en carga activa. **No simétrica** ($\mathbf C_e\boldsymbol\varepsilon$ y $\mathbf M\boldsymbol\varepsilon$ no son proporcionales); el atributo `IS_SYMMETRIC = False` hace que el despachador algebraico ADR 0003 elija LU. Recupera convergencia cuadrática del Newton global frente a la secante.
- **Limitación**: $\varepsilon_{eq}$ simétrica no distingue daño en tensión vs compresión; para hormigón se requiere modelo de Mazars o split tensión/compresión (no implementado). Sin regularización por longitud característica $\to$ *mesh-dependency* en régimen de ablandamiento.

```yaml
materials:
  - id: 3
    type: IsotropicDamage2D
    E: 30.0e9
    nu: 0.2
    kappa_0: 1.5e-4
    alpha: 800.0
    hypothesis: plane_stress
    density: 2500.0   # hormigón
```

## Elasticidad Unilateral (Cable)

### `CableMaterial1D` — cable lineal en tensión, nulo en compresión

Material 1D *memoryless* que captura la propiedad esencial del cable: solo transmite esfuerzo cuando está tensado.

$$\sigma(\varepsilon) = E\,\langle\varepsilon\rangle^+ = \begin{cases} E\,\varepsilon & \varepsilon > 0 \\ 0 & \varepsilon \le 0 \end{cases}$$

con $\langle \cdot \rangle^+$ el operador de MacAuley (parte positiva). Tangente discontinua en $\varepsilon = 0$ (asignación al régimen compresivo por convención).

- **Parámetros**: `E` (módulo de Young en tensión).
- **Sin variables internas** (la respuesta depende sólo de $\varepsilon$ instantánea).
- **Caveats numéricos**: en transiciones tensado ↔ destensado Newton-Raphson puede oscilar. Mitigación: usar `ArcLengthSolver` o pasos finos con `adaptive`.

```yaml
materials:
  - id: 1
    type: CableMaterial1D
    E: 150.0e9
    density: 7850.0
```

## Materiales 3D

Los materiales 3D operan sobre la convención **Voigt 6D** del proyecto (ADR 0012):

$$\boldsymbol\varepsilon = [\varepsilon_{xx},\ \varepsilon_{yy},\ \varepsilon_{zz},\ \gamma_{xy},\ \gamma_{yz},\ \gamma_{xz}]^\top, \qquad \boldsymbol\sigma = [\sigma_{xx},\ \sigma_{yy},\ \sigma_{zz},\ \sigma_{xy},\ \sigma_{yz},\ \sigma_{xz}]^\top$$

con $\gamma_{ij} = 2\varepsilon_{ij}$ (deformación angular *engineering*). Sólo son compatibles con elementos 3D ($\texttt{STRAIN\_DIM} = 6$).

**Ninguno admite `hypothesis`**: las variantes `plane_stress` / `plane_strain` no tienen sentido en 3D, donde todas las componentes están activas. Declararla es un error.

**Caveat de interoperabilidad con ABAQUS**: ABAQUS ordena Voigt como $[11, 22, 33, 12, 13, 23]$, es decir con el bloque cortante permutado ($yz \leftrightarrow xz$) respecto al proyecto. Importar tensores de ABAQUS exige intercambiar las componentes 5 y 6 una sola vez en el preprocesador.

### `Elastic3D` — elástico lineal isótropo 3D

$\boldsymbol\sigma = \mathbf{C}\,\boldsymbol\varepsilon$ con $\mathbf{C}$ isótropa $6\times6$. Tangente constante, sin variables internas.

- **Parámetros**: `E` ($>0$), `nu` $\in (-1, 0.5)$, `density` (opcional).

### `VonMises3D` — plasticidad J2 3D con endurecimiento isótropo lineal

Extensión 3D de `VonMises2D`. Flujo asociado y *return mapping* radial cerrado sobre la parte desviadora (Simó-Hughes §3.3), con tangente algorítmica consistente.

- **Parámetros**: `E`, `nu`, `sigma_y` ($>0$), `H` ($\ge 0$, default 0), `density` (opcional).
- **Variables internas**: `eps_p` (6 componentes) y `alpha` (deformación plástica equivalente).
- **Relación con el 2D**: bajo la restricción $\varepsilon_{zz} = \gamma_{yz} = \gamma_{xz} = 0$ se reduce exactamente a `VonMises2D` en `plane_strain`, verificado a 10 decimales.

### `DruckerPrager3D` — plasticidad friccional 3D

Cono circular suave de Mohr-Coulomb, criterio $f = \sqrt{J_2} + \eta_f I_1 - k(\alpha) \le 0$ con endurecimiento isótropo lineal en la cohesión. Dos ramas cerradas con detección automática: retorno regular a la superficie y retorno al ápice hidrostático.

- **Parámetros**: `E`, `nu`, `cohesion`, `phi_deg` (ángulo de fricción en grados), `psi_deg` (dilatancia en grados, default $=\phi$), `H` (default 0), `variant` (default `outer_cone`), `density` (opcional).
- **Calibraciones**: `outer_cone` (default, circunscribe Mohr-Coulomb) e `inner_cone` (inscribe). La variante `plane_strain_matched` del modelo 2D **no existe en 3D** y el constructor la rechaza con mensaje explícito: es una calibración 2D por construcción, porque Mohr-Coulomb en 3D depende del ángulo de Lode.
- **Plasticidad no asociada**: con $\psi \neq \phi$ la tangente es asimétrica y el despachador algebraico selecciona LU automáticamente.

### `IsotropicDamage3D` — daño isótropo 3D con softening exponencial

$\boldsymbol\sigma = (1-d)\,\mathbf{C}_e\,\boldsymbol\varepsilon$ con daño escalar $d$ gobernado por la deformación equivalente

$$\varepsilon_{eq} = \sqrt{\varepsilon_{xx}^2 + \varepsilon_{yy}^2 + \varepsilon_{zz}^2 + \tfrac{1}{2}\left(\gamma_{xy}^2 + \gamma_{yz}^2 + \gamma_{xz}^2\right)}$$

que es la norma de Frobenius del tensor de deformación y extiende de forma natural la convención 2D. Historial $\kappa$ máximo con condiciones de Kuhn-Tucker y ley de softening exponencial compartida con las versiones 1D y 2D.

- **Parámetros**: `E`, `nu`, `kappa_0` (umbral de daño), `alpha` (parámetro de la ley de softening), `density` (opcional).
- **Tangente**: algorítmica consistente en carga activa (no simétrica en general), secante en descarga y al saturar.

```yaml
materials:
  - {id: 1, type: Elastic3D,   E: 210.0e9, nu: 0.3, density: 7850.0}
  - {id: 2, type: VonMises3D,  E: 210.0e9, nu: 0.3, sigma_y: 250.0e6, H: 1.0e9}
```

## Materiales Cohesivos (salto de desplazamientos)

Familia **paralela** al resto del catálogo: no relacionan $\boldsymbol\sigma$ con $\boldsymbol\varepsilon$ sino la **tracción** $\mathbf{t}$ con el **salto de desplazamientos** $\llbracket u \rrbracket$ a través de una discontinuidad. Viven en su propio registro y se declaran en el bloque `cohesive_materials` del YAML, no en `materials`. Sólo los consumen elementos con discontinuidad embebida como `CST_Embedded2D`.

### `CohesiveDamageIsotropic` — daño cohesivo Modo-I

Daño escalar $\omega \in [0, 1]$ sobre el salto normal: $t_n = (1-\omega)K_e \llbracket u_n \rrbracket$, con $t_s = 0$ (Modo-I puro). La energía de fractura $G_f$ cierra la curva analíticamente, en variante lineal (apertura crítica $w_c = 2G_f/\sigma_{t0}$) o exponencial (asintótica).

- **Parámetros**: `sigma_t0` (resistencia a tracción, Pa), `G_f` (energía de fractura, N/m), `K_e` (rigidez del salto, Pa/m), `softening` $\in$ {`linear`, `exponential`}.
- **Sobre `K_e`**: es un *penalty* que representa la cohesión antes de agrietar, **no** una propiedad del bulk. No tiene default automático; como guía, $K_e \approx 10\,E_{bulk}/\ell_c$. Valores muy altos endurecen el sistema y dificultan la convergencia; valores bajos introducen flexibilidad espuria antes de la activación.
- **Variables internas**: `kappa` (historial del salto) y `damage` ($\omega$).

```yaml
cohesive_materials:
  - {id: 1, type: CohesiveDamageIsotropic, sigma_t0: 2.5e6, G_f: 100.0,
     K_e: 1.0e13, softening: linear}
```
