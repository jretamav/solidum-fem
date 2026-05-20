# ResponseSpectrumSolver — análisis sísmico por combinación modal (ADR 0009 fase 7)

> Spec de **solver nuevo**. Combina las respuestas máximas modales bajo un espectro de respuesta dado (SRSS o CQC). Cierra el ADR 0009 completo y aplica la regla D de refactor arquitectural (free_vibration sale de `results.py` al nuevo `solidum/math/modal_response.py`).

---

## Especificación física

### 0. Descripción general

Análisis sísmico clásico: el suelo se caracteriza por un **espectro de respuesta** (en aceleración `S_a(T)` o desplazamiento `S_d(T)` vs período), no por un acelerograma temporal. Para cada modo del sistema se calcula la respuesta máxima individual `u_n^max = Γ_n·φ_n·S_d(ω_n)` y se combinan las respuestas modales en una envolvente máxima. La fase temporal entre modos se pierde por construcción — sólo se conservan amplitudes envolventes.

Es **lineal** y **estadístico**: aproxima los máximos absolutos sin reconstruir la historia temporal. Útil para diseño normativo donde el espectro reglamentario es la entrada y los máximos modales son la salida solicitada (e.g., ASCE 7, Eurocódigo 8, NCh 433, NTC-DS México).

Hipótesis: lineales (K, M constantes). Para no-linealidad usar transitorio (Newton-Newmark con acelerograma).

### 1. Ecuación de equilibrio resuelta

Para cada modo `n`:

$$u_n^\max = \Gamma_n \cdot \phi_n \cdot S_d(\omega_n)$$

con $\Gamma_n = \phi_n^T M\, r$ el factor de participación modal en la dirección de excitación `r`.

**Combinación SRSS** (Square Root of Sum of Squares):

$$u^\max_k = \sqrt{\sum_n (u_{k,n}^\max)^2}$$

**Combinación CQC** (Complete Quadratic Combination, Der Kiureghian 1980):

$$u^\max_k = \sqrt{\sum_i \sum_j \rho_{ij} \cdot u_{k,i}^\max \cdot u_{k,j}^\max}$$

$$\rho_{ij} = \frac{8\xi^2(1+r) r^{3/2}}{(1-r^2)^2 + 4\xi^2 r (1+r)^2}, \quad r = \omega_i/\omega_j$$

Cuando $\xi \to 0$ o $r \to 0$: $\rho_{ij} \to \delta_{ij}$ y CQC degenera a SRSS.

### 2. Condiciones iniciales / de contorno

- Apoyos Dirichlet constantes (eliminación directa, ADR 0004).
- **Dirección de excitación rígida** `r` (vector unitario o por DOF name): vector `r` con valor 1 en los DOFs traslacionales en la dirección sísmica, 0 en el resto.
- **No** hay condiciones iniciales — la respuesta es envolvente máxima.

### 3. Salidas físicas

`ResponseSpectrumResult` (inmutable):

- `u_combined[ndof]` respuesta máxima envolvente combinada (positiva).
- `u_per_mode[ndof, n_modes]` contribución modal individual con signo.
- `frequencies_rad[n_modes]` frecuencias naturales empleadas.
- `participation_factors[n_modes]` ($\Gamma_n$).
- `effective_masses[n_modes]` ($\Gamma_n^2$). Suma converge a la masa total proyectada en la dirección al incluir suficientes modos.
- `combination` (`"SRSS"` o `"CQC"`).
- `damping` (ξ usado para CQC; 0.0 si SRSS).
- `direction[ndof]` vector de excitación.
- Método `.cumulative_effective_mass_ratio()` para verificar el porcentaje de masa modal capturado.

---

## Formulación numérica

### 4. Esquema operativo

```
1. Ensamble interno vía ModalSolver(n_modes, sigma, lumping):
   → modes, frequencies_rad
2. Resolver direction (np.ndarray o {"dof_name": "ux"}).
3. Resolver spectrum (callable o dict tabulado/constante).
4. Calcular Γ_n = φ_nᵀ M r (participation_factors).
5. Calcular S_d(ω_n) para cada modo.
6. u_per_mode[:, n] = Γ_n · φ_n · S_d(ω_n).
7. Combinar: SRSS (raíz cuadrada de suma de cuadrados elemento a elemento)
   o CQC (forma cuadrática con ρ_ij(ω, ξ)).
8. Devolver ResponseSpectrumResult.
```

### 5. Predictor / corrector

No aplica — un solve por análisis.

### 6. Criterio de convergencia

No aplica como Newton — pero la calidad del análisis depende del número de modos incluidos. Se verifica con
`.cumulative_effective_mass_ratio() ≥ 0.9` (objetivo normativo típico).

### 7. Imposición de Dirichlet y MPC

Eliminación directa heredada del `ModalSolver` interno. MPC fuera de alcance.

### 8. Backend algebraico

Delega en `ModalSolver` para el cálculo de modos (`scipy.sparse.linalg.eigsh` con shift-invert) y luego operaciones de álgebra densa para la combinación modal (los modos son densos por columna).

### 9. Adaptatividad / control de paso

No aplica. El usuario controla `n_modes`. Recomendación: iterativamente subir `n_modes` hasta que la masa efectiva acumulada supere el 90% en la dirección de excitación.

### 10. Caveats numéricos

- **Cuántos modos**: el método pierde precisión si los modos significativos en la dirección de excitación no están incluidos. Verifica con `cumulative_effective_mass_ratio`.
- **CQC vs SRSS**: usar CQC cuando los modos están cercanos en frecuencia (cociente < 1.3 según práctica habitual). Para modos bien separados, SRSS es la opción estándar y CQC se degenera a SRSS automáticamente.
- **Espectro fuera de rango**: `spectrum_tabulated` aplica clamp a los valores extremos (no extrapola). Si el barrido de frecuencias del modelo cae fuera del rango tabulado, el usuario debe extender la tabla.
- **Excitación multi-direccional**: este solver maneja una sola dirección por análisis. Para combinación direccional (CQC3, regla 100/30/30 de norma) se requieren dos o tres análisis independientes y una combinación adicional fuera del solver.
- **Cargas con fase**: no aplica al método espectral por construcción — la fase se pierde en el cómputo de máximos modales.

---

## Contrato de implementación

```yaml
name: ResponseSpectrumSolver
kind: solver
status: validated

interface:
  yaml_type: ResponseSpectrumSolver
  output: ResponseSpectrumResult

parameters:
  - { name: n_modes, type: int, required: true, desc: "Número de modos a calcular y combinar" }
  - { name: direction, type: "ndarray | dict", required: true, desc: "Vector r o {'dof_name': 'ux'}" }
  - { name: spectrum, type: "callable | dict", required: true, desc: "S_d(ω) o {'type': 'tabulated'|'constant_acceleration', ...}" }
  - { name: combination, type: str, required: false, desc: "'SRSS' (default) o 'CQC'" }
  - { name: damping, type: float, required: false, desc: "ξ modal (default 0.05); sólo usado por CQC" }
  - { name: sigma, type: float, required: false, desc: "Shift para ARPACK (default 0)" }
  - { name: lumping, type: str, required: false, desc: "'consistent' (default) o 'lumped'" }

requirements:
  - "Materiales con `density` declarada (ADR 0008)."
  - "Problema lineal — K y M constantes."

conventions:
  units:    "heredadas del modelo (ADR 0008)."
  stability: "No iterativo; estabilidad numérica gobernada por ModalSolver interno."

out_of_scope:
  - "Excitación multi-direccional simultánea (CQC3) — usuario combina afuera."
  - "Espectros con dependencia explícita de ξ por modo (ξ_n distinto) — `damping` es global."
  - "MPC y no linealidad."

acceptance:
  verification:
    - name: single_mode_recovers_gamma_phi_Sd
      setup: "Modelo 2-DOF, n_modes=1, S_d constante=1"
      expect: "u_combined = |u_per_mode[:, 0]| (SRSS trivial sobre un modo)"
      tol_rel: 1.0e-12
    - name: laplaciano_1d_3dof_analitico
      setup: |
        Cadena 4-truss n1-n2-n3-n4-n5 con n1, n5 fijos. E=A=L=1, ρ=1, lumped →
        K Toeplitz tridiagonal [[2,-1,0],[-1,2,-1],[0,-1,2]], M=I. Autovalores
        cerrados del Laplaciano 1D discreto: ω_n² = 2 − 2·cos(nπ/4); modos
        φ_n[i] = √(2/4)·sin(n·i·π/4).
      expect: |
        (1) ω₁ = √(2−√2), ω₂ = √2 a precisión doble.
        (2) Con d = (1, 1/√2, 0): γ₁ = 1, γ₂ = 1/√2, γ₃ = 0 (φ₃ᵀd = 0 por construcción).
        (3) Para Sd ≡ 1 y SRSS: u_combined = (1/√2, 1/√2, 1/√2) exacto.
        (4) Con d = e_{ux₃} (centro): sólo modo 1 contribuye, u_combined = (c/(2√2), c/2, c/(2√2)).
      tol_rel: 1.0e-10
      ref: "tests/test_response_spectrum.py::TestResponseSpectrumAnalytic2DOF"
    - name: srss_manual_combination
      setup: "Modelo 3-DOF, 2 modos, comparar contra raíz cuadrada de suma de cuadrados componente a componente"
      expect: "Coincidencia exacta a precisión de máquina"
      tol_rel: 1.0e-12
    - name: well_separated_modes_cqc_equals_srss
      setup: "Modelo 3-DOF con cociente ω₂/ω₁ > 2"
      expect: "|u_CQC − u_SRSS| / |u_SRSS| < 1%"
      tol_rel: 0.01
  specific:
    - name: close_modes_cqc_differs_from_srss
      setup: "Sistema sintético con ω₂/ω₁ = 1.05 y modos no-ortogonales-en-componentes"
      expect: "Diferencia relativa CQC vs SRSS > 5%"
      tol_rel: 0
    - name: dof_name_builds_unit_vector
      setup: "direction={'dof_name': 'ux'}"
      expect: "direction[node.dofs['ux']] = 1 en todos los nodos con ux"
      tol_rel: 0
    - name: cumulative_effective_mass_monotonic
      setup: "n_modes=2 sobre 3-DOF chain"
      expect: "cumulative crece monotonamente y termina en 1.0"
      tol_rel: 1.0e-10
    - name: free_vibration_wrapper_unchanged
      setup: "Regla D aplicada — ModalResult.free_vibration sigue funcionando"
      expect: "Resultado idéntico a llamada directa a modal_response.free_vibration"
      tol_rel: 1.0e-12

references:
  - "Chopra, A. K. (2017). *Dynamics of Structures*, cap. 13."
  - "Wilson, E. L. (2002). *Three-Dimensional Static and Dynamic Analysis of Structures*, cap. 15."
  - "Der Kiureghian, A. (1980). 'Structural Response to Stationary Excitation'. *J. Eng. Mech. Div. ASCE*, 106(EM6), 1195-1213."
```

---

## Implementación

- Archivo solver: [solidum/math/solvers/response_spectrum.py](../../solidum/math/solvers/response_spectrum.py)
- Clase: `ResponseSpectrumSolver` (registrada en `SolverRegistry`)
- Atributo de despacho: `PIPELINE_KIND = "spectrum"` → `solidum.entry.run_response_spectrum` → `ResponseSpectrumResult`.
- Algoritmos centralizados: [solidum/math/modal_response.py](../../solidum/math/modal_response.py) (`response_spectrum_srss`, `response_spectrum_cqc`, `participation_factors`, helpers `spectrum_from_sa`, `spectrum_tabulated`).
- Resultado: [`ResponseSpectrumResult`](../../solidum/results.py) — dataclass inmutable.
- Tests: [tests/test_response_spectrum.py](../../tests/test_response_spectrum.py) (19 tests verdes).

---

## Diálogo

- *2026-05-18*: regla D aplicada al introducir este solver — segundo método de cómputo sobre `ModalResult` (el primero fue `free_vibration`). El algoritmo de free_vibration se movió a `solidum/math/modal_response.py`; `ModalResult.free_vibration` queda como wrapper delgado de 3 líneas, preservando la API histórica. Los nuevos `response_spectrum_srss` y `response_spectrum_cqc` viven en el mismo módulo. `ResponseSpectrumSolver` orquesta `ModalSolver` interno + delegación a `modal_response`.
- *2026-05-18*: el cuarto valor de `PIPELINE_KIND` (`"spectrum"`) extiende el despacho declarativo de `run_yaml` introducido en D2. Cierra el ADR 0009 completo (fases 1-7 + variante HHT-α).
