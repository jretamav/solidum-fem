# Elastoplastic1D — Plasticidad 1D con endurecimiento isótropo lineal

> Spec **retroactiva**: el material existe desde la fase inicial del proyecto y está validado por tests unitarios y de integración. Esta spec lo documenta sin cambiar comportamiento — H-5.3 de la auditoría 2026-05-18.

---

## Especificación física

### 0. Descripción general

Modelo elastoplástico unidimensional independiente de la velocidad ("rate-independent") con criterio de fluencia $|\sigma| \leq \sigma_y + H\,\alpha$ y **endurecimiento isótropo lineal**. Es el análogo 1D de `VonMises2D` — la versión axial sobre la que se valida y benchmarkean los pipelines plásticos. Se usa en barras (`Truss2D/3D`) y como componente axial conceptual de marcos plásticos (pendiente `FiberSection`).

### 1. Descomposición aditiva

$$\varepsilon = \varepsilon^e + \varepsilon^p$$

con $\varepsilon^p$ histórica (memoria del material).

### 2. Ley elástica

$$\sigma = E\,\varepsilon^e = E\,(\varepsilon - \varepsilon^p)$$

### 3. Criterio de fluencia

$$f(\sigma, \alpha) = |\sigma| - (\sigma_y + H\,\alpha) \;\leq\; 0$$

donde $\sigma_y > 0$ es la fluencia inicial y $H \geq 0$ el módulo de endurecimiento isótropo lineal.

### 4. Regla de flujo (asociada)

$$\dot\varepsilon^p = \dot\gamma\,\operatorname{sign}(\sigma)$$

### 5. Endurecimiento isótropo lineal

$$\dot\alpha = \dot\gamma$$

$\alpha$ es la deformación plástica acumulada equivalente, escalar monótona no decreciente.

### 6. Condiciones de Kuhn-Tucker

$$\dot\gamma \geq 0, \qquad f \leq 0, \qquad \dot\gamma\,f = 0$$

### 7. Variables internas

| Símbolo | Tipo | Significado |
|---|---|---|
| $\varepsilon^p$ | escalar | deformación plástica acumulada con signo |
| $\alpha$ | escalar | deformación plástica acumulada equivalente (≥0) |

`alpha` es la variable principal exportable (`PRIMARY_STATE_VAR = 'alpha'`).

### 8. Notación

Escalares en convención axial: $\sigma > 0 \Leftrightarrow$ tracción, $\varepsilon > 0 \Leftrightarrow$ alargamiento (Reglas §5). $\varepsilon^p$ con signo (puede ser negativa tras carga compresiva).

---

## Formulación numérica

### 9. Esquema temporal — Backward Euler implícito

Dado el estado convergido $\{\varepsilon^p_n, \alpha_n\}$ al paso $n$ y la deformación total $\varepsilon_{n+1}$, resolver para $\{\sigma_{n+1}, \varepsilon^p_{n+1}, \alpha_{n+1}\}$.

### 10. Return mapping — forma cerrada

**Predictor elástico**:

$$\sigma^{\text{trial}} = E\,(\varepsilon_{n+1} - \varepsilon^p_n)$$

**Función de fluencia trial**:

$$f^{\text{trial}} = |\sigma^{\text{trial}}| - (\sigma_y + H\,\alpha_n)$$

- Si $f^{\text{trial}} \leq \text{tol}$ ⇒ **paso elástico**: $\sigma_{n+1} = \sigma^{\text{trial}}$, $E_t = E$, estado intacto.
- Si $f^{\text{trial}} > \text{tol}$ ⇒ **paso plástico**, retorno radial cerrado:

$$\Delta\gamma = \frac{f^{\text{trial}}}{E + H}, \qquad \sigma_{n+1} = \sigma^{\text{trial}} - \Delta\gamma\,E\,\operatorname{sign}(\sigma^{\text{trial}})$$

$$\varepsilon^p_{n+1} = \varepsilon^p_n + \Delta\gamma\,\operatorname{sign}(\sigma^{\text{trial}}), \qquad \alpha_{n+1} = \alpha_n + \Delta\gamma$$

Sin Newton local — la forma cerrada es exacta en 1D.

### 11. Tangente algorítmica consistente

En régimen elástico: $E_t = E$. En régimen plástico:

$$E_t = \frac{E\,H}{E + H}$$

con $H = 0$ recuperando $E_t = 0$ (plasticidad perfecta, descenso de carga indefinido bajo control de desplazamiento).

`IS_SYMMETRIC = True` (trivial en 1D).

### 12. Caveats numéricos

- **Plasticidad perfecta** ($H = 0$): tangente $E_t = 0$ tras la fluencia. Newton-Raphson global puede oscilar o requerir line search en problemas multibarra (modos de mecanismo). Mitigación: usar incremento pequeño cerca del umbral o activar `line_search=True` en el solver.
- **Switching exacto trial elástico ↔ plástico**: la tolerancia `is_admissible` (ADR 0006) usa `admissibility_scale = σ_y + H·α` y patrón atol + rtol·escala — invariancia bajo cambio de unidades.

---

## Contrato de implementación

```yaml
name: Elastoplastic1D
kind: material
status: validated

interface:
  strain_dim: 1
  primary_state_var: alpha   # deformación plástica acumulada equivalente
  is_symmetric: true

parameters:
  - { name: E,       type: float, required: true,  desc: "Módulo de Young (>0)" }
  - { name: sigma_y, type: float, required: true,  desc: "Esfuerzo de fluencia inicial (>0)" }
  - { name: H,       type: float, required: false, default: 0.0,
      desc: "Módulo de endurecimiento isótropo lineal (≥0). H=0 ⇒ plasticidad perfecta" }
  - { name: density, type: float, required: false, default: null,
      desc: "Densidad (kg/m³). Opcional al construir; obligatoria solo si se ensambla peso propio o masa (ADR 0008)" }

signature:
  compute_state: "(ε: float, state_vars=None) -> (σ: float, E_t: float, state_vars')"
  strain_kind: "axial escalar"

state_schema:
  eps_p: "float — deformación plástica acumulada con signo"
  alpha: "float ≥ 0 — deformación plástica acumulada equivalente"

conventions:
  sign: "σ > 0 ⇔ tracción, ε^p con signo (puede ser negativa)"
  units: "[E]=[σ_y]=[H]=Pa; α adimensional"

validity:
  - "E > 0, σ_y > 0, H ≥ 0"
  - "carga proporcional o no proporcional (sin restricción de trayectoria)"
  - "pequeñas deformaciones (|ε| ≲ 1e-2)"

out_of_scope:
  - "endurecimiento cinemático (back-stress) ⇒ requiere variante 1D dedicada"
  - "endurecimiento isótropo no lineal (Voce, exponencial)"
  - "ablandamiento (H<0) ⇒ requiere regularización"
  - "viscoplasticidad / dependencia de la velocidad"
  - "grandes deformaciones (descomposición multiplicativa F = F^e · F^p)"

acceptance:
  # Cubierto por tests/test_materials_unit.py y tests/test_truss_plastic.py
  unit:
    - name: paso_elastico
      setup: "ε = 0.5·σ_y/E (mitad del umbral); state_vars=None"
      expect: "σ = E·ε exacto, α=0, ε^p=0, tangente=E"
      tol_rel: 1.0e-12
    - name: fluencia_uniaxial_traccion
      setup: "ε = 1.5·σ_y/E (1.5× el umbral elástico)"
      expect: "σ = σ_y (con H=0); α=0.5·σ_y/E; ε^p > 0; tangente=0"
      tol_rel: 1.0e-10
    - name: endurecimiento_lineal
      setup: "carga monotónica hasta ε = 2·σ_y/E con H = E/10"
      expect: "σ_final = σ_y + H·α; tangente = E·H/(E+H) en régimen plástico"
      tol_rel: 1.0e-10
    - name: ciclo_carga_descarga
      setup: "ε creciente hasta plástico, luego decreciente"
      expect: "descarga elástica con tangente E; α permanece constante durante descarga; histéresis correcta"
      tol_rel: 1.0e-10
    - name: cambio_de_signo
      setup: "tracción plástica luego compresión hasta fluencia inversa"
      expect: "fluencia inversa en σ = -(σ_y + H·α) (endurecimiento isótropo es simétrico)"
      tol_rel: 1.0e-10
    - name: monotonicidad_alpha
      setup: "trayectoria arbitraria"
      expect: "α no decrece nunca"

  # Cubierto por tests de integración con barras plásticas
  integration:
    - name: truss_uniaxial_plastico
      setup: "Truss2D con Elastoplastic1D, carga axial creciente"
      expect: "respuesta fuerza-desplazamiento bilineal exacta (E hasta fluencia, E_t después)"
      tol_rel: 1.0e-6

references:
  - "Simó J.C., Hughes T.J.R. (1998). Computational Inelasticity. Springer. §1 (J2 1D return mapping)."
  - "de Souza Neto E.A., Perić D., Owen D.R.J. (2008). Computational Methods for Plasticity. Wiley. §6.2 (uniaxial)."
  - "Lubliner J. (1990). Plasticity Theory. Macmillan. §3.2."
```

---

## Implementación

- **Archivo**: [fenix/materials/plastic_1d.py](../../fenix/materials/plastic_1d.py).
- **Clase**: `Elastoplastic1D`, registrada vía `@MaterialRegistry.register`.
- **`admissibility_scale`**: $\sigma_y + H\,\alpha$ (ADR 0006). El switching elástico/plástico usa el patrón atol + rtol·escala.
- **State schema**: `{'eps_p': float, 'alpha': float}`. `state_vars=None` se interpreta como estado virgen ($\varepsilon^p = 0$, $\alpha = 0$).
- **Tests**:
  - [tests/test_materials_unit.py](../../tests/test_materials_unit.py) · clase `TestElastoplastic1D` (carga, descarga, endurecimiento, plasticidad perfecta).
  - Integración: cualquier test de truss con material plástico.

---

## Diálogo

- **2026-05-19** · Spec creada retroactivamente para cerrar el hueco H-5.3 de la auditoría 2026-05-18. Material anterior a la convención de specs; no se modifica código. La aceptación se considera cumplida por los tests unitarios existentes.
