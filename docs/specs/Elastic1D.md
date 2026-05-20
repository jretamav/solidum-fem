# Elastic1D — Elasticidad lineal axial 1D

> Spec **retroactiva**: el material existe desde el origen del proyecto y está validado por tests preexistentes. Esta spec lo documenta sin cambiar comportamiento — H-5.3 de la auditoría 2026-05-18.

---

## Especificación física

### 0. Descripción general

Modelo elástico lineal isótropo unidimensional. Captura la relación $\sigma$–$\varepsilon$ axial de barras, cables traccionados (no unilaterales — para unilateralidad ver `CableMaterial1D`), y la respuesta elástica de cualquier elemento 1D del catálogo (truss, frame Euler en su componente axial). Independiente de la velocidad, sin variables históricas.

Es el modelo más sencillo del catálogo: **una sola constante elástica** (módulo de Young $E$) y una sola ecuación constitutiva.

### 1. Descomposición de la deformación

No aplica — no hay parte plástica. Toda la deformación es elástica:

$$\varepsilon = \varepsilon^e$$

### 2. Ley elástica (Hooke 1D)

$$\sigma = E\,\varepsilon$$

con $E > 0$ módulo de Young. Tangente constante $E_t = E$ en todo régimen.

### 3. Superficie de fluencia / criterio de daño

No aplica — el material es indefinidamente elástico.

### 4. Regla de flujo

No aplica.

### 5. Endurecimiento

No aplica.

### 6. Condiciones de Kuhn-Tucker

No aplica.

### 7. Variables internas

Ninguna — `state_vars = None` se acepta y se devuelve intacto. Sin `PRIMARY_STATE_VAR`.

### 8. Notación

Esfuerzo y deformación escalares en convención axial: $\sigma > 0 \Leftrightarrow$ tracción, $\varepsilon > 0 \Leftrightarrow$ alargamiento (Reglas §5).

---

## Formulación numérica

### 9. Esquema temporal

No aplica — el material no tiene historia. Cada evaluación de `compute_state` es independiente.

### 10. Return mapping / actualización

No aplica. La evaluación se reduce a una multiplicación escalar:

$$\sigma_{n+1} = E\,\varepsilon_{n+1}, \qquad E_t = E$$

### 11. Tangente algorítmica consistente

$$C^{\text{alg}} = E$$

Coincide con la tangente continua y es simétrica trivialmente. `IS_SYMMETRIC = True` (por default del contrato base).

### 12. Caveats numéricos

Ninguno — el modelo es perfectamente bien condicionado para $E > 0$.

---

## Contrato de implementación

```yaml
name: Elastic1D
kind: material
status: validated

interface:
  strain_dim: 1
  primary_state_var: null    # sin variables internas
  is_symmetric: true

parameters:
  - { name: E,       type: float, required: true,  desc: "Módulo de Young (>0)" }
  - { name: density, type: float, required: false, default: null,
      desc: "Densidad (kg/m³). Opcional al construir; obligatoria solo si se ensambla peso propio o masa (ADR 0008)" }

signature:
  compute_state: "(ε: float, state_vars=None, **kwargs) -> (σ: float, E_t: float, state_vars)"
  strain_kind: "axial escalar"

state_variables: []        # ninguna

conventions:
  sign: "σ > 0 ⇔ tracción, ε > 0 ⇔ alargamiento (Reglas §5)"
  units: "[E] = [σ] / [ε] = Pa (o consistente con el sistema)"

validity:
  - "E > 0"
  - "régimen elástico indefinido — sin umbral superior"
  - "pequeñas deformaciones (|ε| ≲ 1e-2 por hipótesis cinemática global)"

out_of_scope:
  - "endurecimiento / fluencia ⇒ usar Elastoplastic1D"
  - "softening por daño ⇒ usar IsotropicDamage1D"
  - "unilateralidad (sin compresión) ⇒ usar CableMaterial1D"
  - "anisotropía o dependencia de temperatura"
  - "viscoelasticidad"

acceptance:
  # Cubierto por tests existentes de elementos 1D (Truss, Frame, Cable, etc.)
  # que usan Elastic1D como material por default. La verificación se hace
  # implícita en cada test de elemento.
  verification:
    - name: ley_de_hooke_exacta
      setup: "evaluar compute_state con ε arbitrario; comparar σ = E·ε"
      expect: "σ = E·ε exacto; tangente = E"
      tol_rel: 1.0e-14
    - name: tangente_constante
      setup: "evaluar dos veces con strain distintos"
      expect: "E_t devuelto idéntico en ambas evaluaciones"
      tol_rel: 1.0e-14
    - name: ciclo_carga_descarga_sin_historia
      setup: "ε creciente luego decreciente; verificar reversibilidad σ↔ε"
      expect: "σ(ε)=E·ε para todo ε (sin histéresis, sin state_vars retornados modificados)"
      tol_rel: 1.0e-14
    - name: rechazo_E_no_positivo
      setup: "construir Elastic1D(E=0) y Elastic1D(E=-1)"
      expect: "ValueError con mensaje claro"
    - name: rechazo_density_negativa
      setup: "construir Elastic1D(E=210e9, density=-1)"
      expect: "ValueError con mensaje claro"

references:
  - "Cualquier texto elemental de mecánica de materiales (Beer, Hibbeler, Timoshenko §1)."
```

---

## Implementación

- **Archivo**: [solidum/materials/elastic.py](../../solidum/materials/elastic.py).
- **Clase**: `Elastic1D`, registrada vía `@MaterialRegistry.register`.
- **Estado**: `state_vars=None` permitido; se acepta y se devuelve sin tocar (compatibilidad con elementos no lineales que pasan diccionarios históricos).
- **`admissibility_scale`**: hereda el comportamiento por defecto del contrato `Material` (no usado — sin criterio de admisibilidad).
- **Tests**: ningún test unitario propio de `Elastic1D` (es trivial); su correctitud se ejercita implícitamente en todos los tests de elementos 1D que lo usan como material (decenas de tests en `tests/test_truss*.py`, `tests/test_frame*.py`, `tests/test_cable*.py`, `tests/test_dynamic*.py`).

---

## Diálogo

- **2026-05-19** · Spec creada retroactivamente para cerrar el hueco H-5.3 de la auditoría 2026-05-18. El material es anterior a la convención de specs validadas; no se modifica código. La aceptación se considera cumplida por cobertura indirecta en tests de elementos.
