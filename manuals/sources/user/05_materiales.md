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
