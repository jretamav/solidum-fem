# MATERIAL DE CABLE 1D — elasticidad unilateral

> Orden de trabajo. El usuario escribe **especificación física**, **formulación** y **contrato**; la IA rellena **implementación** y responde en **diálogo**.

---

## Especificación física

### 0. Descripción general
Material 1D memoryless que modela la respuesta axial de un cable: **elástico lineal en tensión y nulo en compresión**. No tiene variables internas ni historia — la respuesta depende exclusivamente de la deformación instantánea. Captura la propiedad esencial de un cable: que solo transmite esfuerzo cuando está tensado.

### 1. Ecuación constitutiva
$$\sigma(\varepsilon) \;=\; E\,\langle\varepsilon\rangle^+ \;=\; \begin{cases} E\,\varepsilon & \varepsilon > 0 \\ 0 & \varepsilon \leq 0 \end{cases}$$

Rampa lineal de pendiente $E$ en tensión ($\varepsilon > 0$), cero en compresión o estado sin deformación ($\varepsilon \leq 0$). $\langle\cdot\rangle^+$ denota la parte positiva (operador de MacAuley).

### 2. Módulo tangente
$$E_t(\varepsilon) \;=\; \frac{d\sigma}{d\varepsilon} \;=\; \begin{cases} E & \varepsilon > 0 \\ 0 & \varepsilon \leq 0 \end{cases}$$

**Discontinuo en $\varepsilon = 0$**. Por convención el valor en el punto de corte se asigna al régimen compresivo ($E_t = 0$), lo que evita generar rigidez espuria cuando el cable está exactamente sin tensión.

### 3. Variables internas
No hay. La respuesta es puramente función de $\varepsilon$ actual; el estado del material es el vacío.

### 4. Interpretación física
- $\varepsilon > 0$: cable tensado — se comporta como barra elástica.
- $\varepsilon < 0$: cable destensado — no transmite nada, la carga pasa por otros caminos del modelo.
- $\varepsilon = 0$: frontera neutra — no hay tensión, no hay rigidez.

La no linealidad es **constitutiva** (en la respuesta del material), no geométrica. Un elemento con este material puede tener cualquier cinemática (lineal o corotacional); la unilateralidad solo vive aquí.

---

## Contrato de implementación

```yaml
name: CableMaterial1D
kind: material
status: validated

interface:
  strain_dim: 1
  primary_state_var: null      # memoryless, no exporta variable principal

parameters:
  - { name: E, type: float, required: true, desc: "Módulo de Young (en tensión)" }

signature:
  compute_state: "(ε: float, state_vars=None) -> (σ: float, E_tangent: float, state_vars)"
  strain_kind: axial escalar

conventions:
  sign: "σ > 0 ⇔ ε > 0 (elongación); ε ≤ 0 ⇒ σ = 0"
  state_passthrough: "state_vars se devuelve intacto (el material no almacena historia)"

validity:
  - "E > 0"
  - "|ε| ≲ 1e-2 en el régimen tensado (elasticidad lineal)"
  - "IS_UNILATERAL = True ⇒ solo admitido en elementos con ACCEPTS_UNILATERAL = True (Cable2DCorot, Cable3DCorot). La base Element rechaza la combinación con Truss* en construcción."

out_of_scope:
  - plasticidad, fluencia, fatiga
  - viscoelasticidad / dependencia temporal
  - pretensión inicial (el usuario la modela con una longitud de referencia ajustada)
  - regularización numérica en compresión (no hay rigidez residual en ε ≤ 0)
  - anisotropía direccional; la respuesta es puramente axial escalar

numerical_caveats:
  - "En transiciones tensado ↔ destensado (ε cruza 0), E_tangent salta de E a 0 y Newton-Raphson puede oscilar. Mitigación: usar arc-length o pasos de carga finos si el problema tiene destensiones reales."
  - "Un cable completamente destensado tiene K_M = 0 y K_G = 0 ⇒ matriz tangente singular en los DOFs del elemento. El usuario debe garantizar que otros elementos del modelo proporcionen rigidez suficiente, o pretensar el cable antes de aplicar cargas que puedan destensarlo."

acceptance:
  - name: respuesta en tensión pura
    setup: "ε = 1e-4"
    expect: "σ = E · 1e-4, E_tangent = E"
    tol_rel: 1.0e-12

  - name: respuesta en compresión pura
    setup: "ε = -1e-4"
    expect: "σ = 0, E_tangent = 0"
    tol_abs: 1.0e-15

  - name: punto de corte
    setup: "ε = 0"
    expect: "σ = 0, E_tangent = 0 (régimen compresivo por convención)"
    tol_abs: 1.0e-15

  - name: state_vars se devuelve intacto
    setup: "state_vars = {'dummy': 42}, ε arbitrario"
    expect: "el tercer retorno es idéntico al state_vars de entrada"

references:
  - "Irvine H.M., Cable Structures, MIT Press, §1.2 (modelo elástico del cable)"
  - "MacAuley M.A., on the deflection of beams, 1919 (operador parte positiva)"
```

---

## Implementación

- **Archivo**: [solidum/materials/cable_1d.py](../../solidum/materials/cable_1d.py) — archivo propio, sin dependencia de otros módulos de material salvo la clase base `Material`.
- **Clase**: `CableMaterial1D` — hereda directamente de `Material` (no de `Elastic1D` ni de ningún otro material), registrada vía `@MaterialRegistry.register`.
- **Tests**: [tests/test_cable_material.py](../../tests/test_cable_material.py) · `TestCableMaterial1D` — seis tests:
  - `test_acceptance_tension_pura` (criterio 1).
  - `test_acceptance_compresion_pura` (criterio 2).
  - `test_acceptance_punto_de_corte` (criterio 3).
  - `test_acceptance_state_vars_passthrough` (criterio 4).
  - `test_registro_en_registry` — verifica autodiscover.
  - `test_validacion_E_positivo` — rechaza $E \leq 0$.
- **Notas de traducción**:
  - La ramificación en `compute_state` es un único `if strain > 0.0`. El `>` estricto implementa la convención del punto de corte (el cero se asigna al régimen compresivo). El `else` cubre tanto $\varepsilon < 0$ como $\varepsilon = 0$ con las mismas salidas.
  - `state_vars` se devuelve como llega (identidad, no copia). El test `test_acceptance_state_vars_passthrough` usa `assertIs` para confirmar identidad.
  - Validación de $E > 0$ al construir, con mensaje claro. Atrapa typos del input YAML temprano.
  - El material acepta `**kwargs` en `compute_state` por compatibilidad con elementos que pudieran pasar argumentos adicionales (convención del proyecto).

---

## Diálogo

- **2026-04-21** · Material creado como pieza independiente para la futura implementación de elementos de cable (`Cable2DCorot`, `Cable3DCorot`). Deliberadamente **no hereda** de `Elastic1D` aunque su rama en tensión sea idéntica: mantener el material autocontenido simplifica su lectura y evita acoplamientos futuros si la rama elástica cambia en otro archivo.
- **2026-04-21** · Discontinuidad en $\varepsilon = 0$: documentada explícitamente en `out_of_scope: regularización numérica`. Si aparece un caso de uso con destensiones frecuentes donde Newton-Raphson oscile, se añadirá una variante `CableMaterial1DRegularized` con $E_c \ll E$ en compresión — **no** se modificará este material (contrato limpio).
- **2026-05-12** · Restricción de emparejamiento elevada a contrato: `IS_UNILATERAL = True` en el material; `ACCEPTS_UNILATERAL = True` en `Cable2DCorot`/`Cable3DCorot`. La base `Element._validate_material_compatibility` rechaza la combinación con `Truss*` (cuyo tangente colapsa la rigidez global sin K_G que lo compense). Antes el contrato lo permitía y solo se desaconsejaba en docstring; ahora falla limpio al construir.
