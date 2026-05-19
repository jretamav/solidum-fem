# Catálogo de materiales

> Referencia rápida de los modelos constitutivos implementados. Una entrada por material. Para detalles del algoritmo → código fuente.
>
> **Convenciones**: `STRAIN_DIM` = dimensión Voigt de la deformación de entrada (1 = escalar, 3 = 2D `[ε_xx, ε_yy, γ_xy]`, 6 = 3D). `PRIMARY_STATE_VAR` = variable interna que el `VtkExporter` lleva al output.
>
> **Densidad (ADR 0008)**: todos los materiales aceptan un parámetro opcional `density` (kg/m³ en las unidades del problema). No requerido al construir — análisis estáticos sin peso propio funcionan sin declararla. Cuando un consumidor que requiere masa la pide (`Assembler.assemble_self_weight(g)` y futuras `assemble_mass_matrix`, etc.) y la encuentra como `None`, **falla con `ValueError`** listando los materiales sin densidad. Valor `0.0` declarado explícitamente cubre el caso legítimo de material sin masa física (penalty, restricción) y emite warning informativo.

---

## Elastic1D — elástico lineal 1D

- **Ley**: `σ = E · ε` (Hooke 1D).
- **STRAIN_DIM**: 1.
- **Parámetros**: `E` (módulo de Young).
- **Variables internas**: ninguna (sin historia).
- **Tangente**: constante = `E`.
- **Compatible con**: `Truss2D`, `Truss3D`, `Frame2DEuler`, `Frame2DTimoshenko`.
- **Spec**: [docs/specs/Elastic1D.md](specs/Elastic1D.md)
- **Archivo**: [fenix/materials/elastic.py](fenix/materials/elastic.py)

---

## Elastic2D — elástico lineal 2D

- **Ley**: `σ = C · ε`, con `C` el tensor constitutivo isótropo en notación Voigt `[xx, yy, xy]`.
- **STRAIN_DIM**: 3.
- **Parámetros**: `E`, `nu`, `hypothesis ∈ {'plane_stress', 'plane_strain'}`.
- **Variables internas**: ninguna.
- **Tangente**: constante = `C`.
- **Compatible con**: `Quad4`, `Tri3`, `Tri6`, `Quad8`, `Quad9` (todos los sólidos 2D con `STRAIN_DIM = 3`).
- **Spec**: [docs/specs/Elastic2D.md](specs/Elastic2D.md)
- **Archivo**: [fenix/materials/elastic_2d.py](fenix/materials/elastic_2d.py)

---

## Elastic3D — elástico lineal isótropo 3D

- **Ley**: `σ = C · ε`, con `C` el tensor constitutivo isótropo 6×6 en notación Voigt 3D del proyecto `[xx, yy, zz, xy, yz, xz]` (ADR 0012, `Reglas.md §5`).
- **STRAIN_DIM**: 6.
- **Parámetros**: `E` (>0), `nu ∈ (-1, 0.5)`, `density` (opcional, ADR 0008).
- **Variables internas**: ninguna.
- **Tangente**: constante = `C`.
- **Hipótesis**: no aplican `plane_stress`/`plane_strain` en 3D — toda la deformación es activa.
- **Compatible con**: `Hex8`, `Tet4` (todos los sólidos 3D con `STRAIN_DIM = 6`).
- **Caveat de compatibilidad**: ABAQUS usa orden Voigt `[11, 22, 33, 12, 13, 23]` (permutación `yz ↔ xz` respecto al proyecto). Importar datos de ABAQUS requiere intercambiar los componentes 5↔6 una vez en el preprocesador.
- **Spec**: [docs/specs/Elastic3D.md](specs/Elastic3D.md)
- **Archivo**: [fenix/materials/elastic_3d.py](../fenix/materials/elastic_3d.py)

---

## IsotropicDamage1D — daño isótropo 1D con softening exponencial

- **Modelo**: daño escalar `d ∈ [0, DAMAGE_MAX]`, esfuerzo `σ = (1 − d) · E · ε`. Versión 1D de [`IsotropicDamage2D`](#isotropicdamage2d--daño-isótropo-2d-con-softening-exponencial).
- **Evolución del daño** (Kuhn-Tucker):
  - Deformación equivalente: `ε_eq = |ε|`.
  - Variable histórica: `κ = max(κ_old, ε_eq)`.
  - Si `κ ≤ κ_0` → `d = 0`. Si no: `d = 1 − (κ_0/κ) · exp(−α·(κ − κ_0))`, saturado en `DAMAGE_MAX`.
- **Tangente**: **algorítmica consistente** en carga activa con daño no saturado; **secante** `(1-d)·E` en descarga, sin daño activo o al saturar.
  - Fórmula consistente: `E_tan = -(1-d)·E·α·κ` — **negativa** durante toda la fase de softening (refleja la pendiente descendente σ-ε post-pico).
  - Derivación: de `(1-d)·E - E·ε·(∂d/∂κ)·sign(ε)` con `∂d/∂κ = (1-d)(1/κ + α)` y `sign(ε)·ε = |ε| = κ` en carga.
- **STRAIN_DIM**: 1 · **PRIMARY_STATE_VAR**: `'damage'` · **IS_SYMMETRIC**: `True` (la tangente escalar produce contribución elemental simétrica aunque pueda ser indefinida; el despachador algebraico ADR 0003 detecta no-PD y degrada Cholesky→LU automáticamente).
- **Parámetros**: `E`, `kappa_0` (umbral elástico), `alpha` (velocidad de degradación), `density` (opcional, ADR 0008).
- **Variables internas**: `kappa`, `damage`.
- **Admisibilidad (ADR 0006)**: `admissibility_scale = E · κ_0` — esfuerzo umbral inicial. Garantiza invariancia bajo cambio de unidades del check `f ≤ atol + rtol · escala`.
- **Limitaciones**: no distingue tracción/compresión en la activación del daño; control de carga global puede ser inestable tras el pico (requiere `ArcLengthSolver` para seguir la rama de softening).
- **Compatible con**: `Truss2D`, `Truss3D`.
- **Referencia**: ver `docs/specs/IsotropicDamage1D.md` y la hermana 2D para derivación detallada. Lemaitre & Chaboche, *Mechanics of Solid Materials*. Simó & Ju (1987, IJSS 23) marco de tangente consistente.
- **Archivo**: [fenix/materials/damage_1d.py](fenix/materials/damage_1d.py)

---

## IsotropicDamage2D — daño isótropo 2D con softening exponencial

- **Modelo**: extensión 2D de `IsotropicDamage1D`. Esfuerzo nominal `σ = (1 − d) · C_e · ε`; daño escalar `d ∈ [0, DAMAGE_MAX]`.
- **Deformación equivalente**: `ε_eq = √(ε_xx² + ε_yy² + ½·γ_xy²)` (norma simple sin distinguir tensión/compresión).
- **Evolución**: κ máxima histórica (Kuhn-Tucker `κ_n+1 = max(κ_n, ε_eq)`), ley exponencial `d(κ) = 1 - (κ_0/κ)·exp(-α(κ-κ_0))` para `κ > κ_0`, saturado a `DAMAGE_MAX`.
- **Tangente**: **algorítmica consistente** en carga activa con daño no saturado; **secante** `(1-d)·C_e` en descarga, sin daño activo (`κ ≤ κ_0`) o al saturar.
  - Fórmula consistente: `C_alg = (1-d)·C_e - [(1-d)·(1/κ + α)/ε_eq] · (C_e·ε) ⊗ (M·ε)` con `M = diag(1, 1, 1/2)`.
  - Derivación: `∂d/∂κ = (1-d)·(1/κ + α)` por la ley exponencial; `∂κ/∂ε = M·ε/ε_eq` en carga activa, 0 en descarga.
  - **No simétrica en general** (`(C_e·ε)` y `(M·ε)` no son proporcionales). Recupera convergencia cuadrática del Newton global frente a la secante (que daba convergencia lineal).
- **STRAIN_DIM**: 3 · **PRIMARY_STATE_VAR**: `'damage'` · **IS_SYMMETRIC**: `False` (la tangente consistente no es simétrica; el despachador algebraico ADR 0003 elige LU en lugar de Cholesky para el sistema global).
- **Parámetros**: `E`, `nu`, `kappa_0`, `alpha`, `hypothesis ∈ {'plane_stress' (default), 'plane_strain'}`, `density` (opcional, ADR 0008).
- **Variables internas**: `kappa`, `damage`.
- **Admisibilidad (ADR 0006)**: `admissibility_scale = E · κ_0` — idéntica al caso 1D, el criterio de daño tiene el mismo umbral característico independientemente de la dimensión.
- **Limitación**: la `ε_eq` simétrica no distingue daño en tensión vs compresión; para hormigón usar modelo de Mazars o split tensión/compresión (out-of-scope). Sin regularización por longitud característica → mesh-dependency en régimen de ablandamiento (banda localizada en una fila de elementos).
- **Compatible con**: `Quad4`, `Tri3`, `Tri6`, `Quad8`, `Quad9`.
- **Referencia**: ver `docs/specs/IsotropicDamage2D.md`. Simó & Ju (1987, IJSS 23) tangente consistente para daño isótropo. Lemaitre & Chaboche (1990) marco de continuum damage mechanics.
- **Archivo**: [fenix/materials/damage_2d.py](fenix/materials/damage_2d.py)

---

## Elastoplastic1D — plasticidad J2 1D con endurecimiento isótropo lineal

- **Modelo**: plasticidad asociativa con criterio `f = |σ_trial| − (σ_y + H·α) ≤ 0`.
- **Algoritmo**: return mapping clásico 1D.
  - Predictor elástico: `σ_trial = E · (ε − ε_p_old)`.
  - Corrector: `Δγ = f / (E + H)`, `σ = σ_trial − Δγ·E·sign(σ_trial)`.
  - Tangente algorítmica consistente: `E_t = E·H / (E + H)`.
- **STRAIN_DIM**: 1 · **PRIMARY_STATE_VAR**: `'alpha'`.
- **Parámetros**: `E`, `sigma_y` (fluencia inicial), `H` (módulo de endurecimiento; 0 = perfectamente plástico).
- **Variables internas**: `eps_p` (deformación plástica), `alpha` (acumulada equivalente).
- **Admisibilidad (ADR 0006)**: `admissibility_scale = σ_y + H · α` — fluencia corriente (adaptativa al endurecimiento). El check `f ≤ atol + rtol · escala` se hace contra la superficie de fluencia *en el estado de entrada del paso*.
- **Compatible con**: `Truss2D`, `Truss3D`, `Frame2DEuler`, `Frame2DTimoshenko` (en estos últimos solo se aplica al esfuerzo axial).
- **Referencia**: Simo & Hughes, *Computational Inelasticity*, cap. 1.
- **Spec**: [docs/specs/Elastoplastic1D.md](specs/Elastoplastic1D.md)
- **Archivo**: [fenix/materials/plastic_1d.py](fenix/materials/plastic_1d.py)

---

## CableMaterial1D — cable axial 1D con elasticidad unilateral

- **Ley**: `σ = E·ε` si `ε > 0`; `σ = 0` si `ε ≤ 0` (sin rigidez en compresión).
- **STRAIN_DIM**: 1.
- **Parámetros**: `E` (módulo de Young en tensión, estrictamente positivo).
- **Variables internas**: ninguna (respuesta memoryless: la rigidez depende solo del signo instantáneo de `ε`).
- **Tangente**: `E_t = E` en tensión, `E_t = 0` en compresión o sin tensión.
- **Implicación numérica**: la rigidez tangente colapsa a cero cuando el cable se afloja, lo que requiere precondicionamiento o regularización en el solver para no degenerar la matriz global. El uso típico es a través del elemento `Cable2DCorot`/`Cable3DCorot`, que detecta la situación y la maneja correctamente.
- **Compatible con**: `Cable2DCorot`, `Cable3DCorot`.
- **Referencia**: ver `docs/specs/CableMaterial1D.md`.
- **Archivo**: [fenix/materials/cable_1d.py](fenix/materials/cable_1d.py)

---

## DruckerPrager2D — plasticidad friccional cohesivo-friccional 2D plane strain

- **Modelo**: cono circular suave de Mohr-Coulomb con cohesión y fricción interna. Criterio `f = √J₂ + η_f·I₁ − k(α) ≤ 0` con `k(α) = k_0 + H·α` (endurecimiento isótropo lineal en cohesión).
- **Hipótesis**: solo `plane_strain` en esta entrega. `plane_stress` declarado *out-of-scope* (proyección con σ_zz=0 acoplada a flujo dilatante es notoriamente delicada).
- **Plasticidad no asociada por defecto**: ángulo de dilatancia `ψ` parámetro independiente, default `ψ = φ` (asociada). Para `ψ ≠ φ` la tangente algorítmica es **asimétrica** → `IS_SYMMETRIC = False`.
- **Calibración con Mohr-Coulomb** (variant):
  - `plane_strain_matched` (default): η_f = tan(φ)/√(9 + 12tan²(φ)), k_0 = 3c_0/√(9 + 12tan²(φ)). Coincide exactamente con MC en plane strain.
  - `outer_cone`: η_f = 2·sin(φ)/[√3·(3 - sin(φ))]. Circunscribe MC.
  - `inner_cone`: η_f = 2·sin(φ)/[√3·(3 + sin(φ))]. Inscribe MC.
- **Algoritmo — dos ramas cerradas con detección automática**:
  - **Return regular** (cone surface): `Δγ = f_trial / (G + 9K·η_f·η_g + H)`, dirección desviadora preservada. Tangente algorítmica `C_alg = K·v⊗v + 2G·(1-β)·I_dev + 4G·β·n̂⊗n̂ - (1/A)·b_g⊗b_f` con `β = G·Δγ/√J₂_trial`, `b_g = 2G·n̂ + 3K·η_g·v`, `b_f = 2G·n̂ + 3K·η_f·v`.
  - **Return al ápice** (vértice del cono, zona puramente hidrostática): `Δγ_apex = (I₁_trial·η_f - k(α))/(9K·η_f·η_g + H)`, σ_n+1 = (k/3η_f)·I puramente hidrostático. Tangente reducida a `K·H/(9K·η_f·η_g + H) · v⊗v` (simétrica, sin rigidez desviadora).
  - Detección regular ↔ apex: si tras `Δγ_regular` resulta `√J₂_new < 0` → cae en apex; cambia a la fórmula apex.
- **STRAIN_DIM**: 3 · **PRIMARY_STATE_VAR**: `'alpha'` · **IS_SYMMETRIC**: `False` (declarativo, conservador; el despachador algebraico ADR 0003 elige LU).
- **Parámetros** (físicos):
  - `E`, `nu`: elásticos.
  - `cohesion`: cohesión inicial `c_0` (esfuerzo).
  - `phi_deg`: ángulo de fricción interna (grados, `0 ≤ φ < 90`).
  - `psi_deg`: ángulo de dilatancia (grados, `0 ≤ ψ ≤ φ`; default `φ` ⇒ asociada).
  - `H` (≥0, default 0): endurecimiento isótropo lineal en cohesión.
  - `variant`: calibración con MC (default `'plane_strain_matched'`).
  - `density` (opcional, ADR 0008).
- **Variables internas**: `eps_p` (tensorial 4 componentes `[xx, yy, zz, xy_tensorial]`), `alpha` (multiplicador plástico acumulado, adimensional).
- **Admisibilidad (ADR 0006)**: `admissibility_scale = k(α) = k_0 + H·α` (cohesión efectiva corriente, unidades de esfuerzo igual que `f`).
- **Limitaciones** (out-of-scope):
  - Sin `plane_stress`, endurecimiento por fricción/dilatancia, endurecimiento cinemático, cap (modelo cap-cone), ablandamiento (`H<0` requiere regularización), grandes deformaciones.
  - Validación inicial cubierta por tests unitarios (incluyendo FD numérica de la tangente) y un benchmark de pipeline Quad4 (rama elástica + rama apex + caso asociado). El régimen "rama regular pura" en geometría confinada es difícil de calibrar sin caer en apex; se cubre por los tests unitarios de cortante puro.
- **Compatible con**: `Quad4`, `Tri3`, `Tri6`, `Quad8`, `Quad9` (todos con material 2D plane strain).
- **Referencia**: ver `docs/specs/DruckerPrager2D.md`. Drucker & Prager (1952). de Souza Neto, Perić & Owen (2008) cap. 8 (Drucker-Prager, tangentes algorítmicas).
- **Archivo**: [fenix/materials/drucker_prager_2d.py](fenix/materials/drucker_prager_2d.py)

---

## VonMises2D — plasticidad J2 2D con endurecimiento isótropo lineal

- **Modelo**: J2 (Von Mises) con regla de flujo asociada y endurecimiento isótropo lineal. Soporta **dos hipótesis cinemáticas** mutuamente excluyentes (selección por kwarg `hypothesis` al construir): `plane_strain` (`ε_zz = 0`) y `plane_stress` (`σ_zz = 0`). Cada hipótesis usa un kernel Numba especializado.
- **Algoritmo plane strain** (Simó-Hughes §3.3): return mapping radial cerrado sobre la parte desviadora 3D extendida. Predictor elástico desviador `s_trial = 2G·(e_dev − e_p)` con presión `p = K·tr(ε)`; criterio `f = ‖s_trial‖ − √(2/3)·(σ_y + H·α) ≤ 0`; corrector radial `Δγ = f_trial / (2G + ⅔·H)` con `N = s_trial/‖s_trial‖`; actualización `α_new = α + √(2/3)·Δγ`. Tangente algorítmica consistente cerrada.
- **Algoritmo plane stress** (Simó-Hughes §3.4.1, "plane stress projected algorithm"): predictor `σ_trial = C_e_ps · (ε − e_p)` en subespacio plane stress; función de fluencia proyectada `f̄ = ½·σ·P·σ − R²/3` con operador `P = (1/3)·[[2,-1,0],[-1,2,0],[0,0,6]]`. Newton local escalar sobre `Δγ` en base de autovectores ortonormales `v_1=[1,1,0]/√2`, `v_2=[1,-1,0]/√2`, `v_3=[0,0,1]` de `C_e_ps·P` con autovalores `μ_1=E/(3(1-ν))`, `μ_2=μ_3=2G`. Converge en 3-6 iteraciones (tangente cerrada). Actualización de `α` con la regla físicamente correcta `α_new = α + Δγ·√(2σPσ/3)` (la heurística `√(2/3)·Δγ` válida en plane strain rompe la invariancia bajo cambio de unidades en plane stress, donde `Δγ` tiene unidades `[1/esfuerzo]`). Incompresibilidad plástica cierra `e_p_zz = -(e_p_xx + e_p_yy)`. Tangente algorítmica consistente cerrada con corrección por `dα/dΔγ ≠ √(2/3)`.
- **STRAIN_DIM**: 3 · **PRIMARY_STATE_VAR**: `'alpha'`.
- **Parámetros**: `E`, `nu`, `sigma_y`, `H` (≥0, default 0), `hypothesis ∈ {'plane_strain' (default), 'plane_stress'}`, `density` (opcional, ADR 0008).
- **Variables internas**: `eps_p` (tensorial 4 componentes `[xx, yy, zz, xy_tensorial]`), `alpha` (acumulada equivalente, adimensional en ambas hipótesis).
- **Admisibilidad (ADR 0006)**:
  - plane_strain: `admissibility_scale = √(2/3)·(σ_y + H·α)` — radio desviador de la superficie J2.
  - plane_stress: `admissibility_scale = (σ_y + H·α)²/3` — escala de `f̄` proyectada.
  - En ambos, la tolerancia se precomputa fuera del kernel `@njit` y se pasa como `float`.
- **Limitaciones declaradas** (ver `out_of_scope` en spec): no endurecimiento cinemático/Voce, no ablandamiento (`H<0` requiere regularización), no regla de flujo no asociada, no daño acoplado, no viscoplasticidad, no grandes deformaciones, no anisotropía. Locking volumétrico en plane strain con elementos de bajo orden cuando `ν → 0.5` o dominio plástico (mitigaciones B-bar/F-bar diferidas).
- **Implementación**: kernels `_compute_j2_plane_strain` y `_compute_j2_plane_stress` con `@njit`. Despacho por `hypothesis` en construcción (sin coste runtime).
- **Compatible con**: `Quad4`, `Tri3`, `Tri6`, `Quad8`, `Quad9` (configurados en cualquiera de las dos hipótesis).
- **Referencia**: ver `docs/specs/VonMises2D.md`. Simó & Hughes, *Computational Inelasticity* (1998), §3.3 (plane strain), §3.4 (plane stress projected, Box 3.1). de Souza Neto, Perić & Owen, *Computational Methods for Plasticity* (2008), §9.4.
- **Archivo**: [fenix/materials/von_mises_2d.py](fenix/materials/von_mises_2d.py)

---

---

# Materiales cohesivos (familia paralela, ADR 0010)

Los materiales cohesivos son una **jerarquía paralela e independiente** de los materiales continuos: operan sobre el salto de desplazamientos `[[u]]` sobre una superficie de discontinuidad `Γ_d` y devuelven tracciones `t`, no esfuerzos sobre `ε`. Viven en `fenix/cohesive_materials/`, heredan de `fenix.core.cohesive_material.CohesiveMaterial`, se registran vía `CohesiveMaterialRegistry` y se declaran en YAML bajo la sección `cohesive_materials` (separada de `materials`). El bulk del elemento sigue siendo un `Material` continuo; el cohesivo entra al sistema sólo cuando el elemento activa una discontinuidad embebida (ver ADR 0010).

> **Convenciones de la familia**: `JUMP_DIM` = dimensión del vector de salto (2 en 2D, 3 en 3D); `PRIMARY_STATE_VAR` = variable interna que se exporta al post-proceso; `IS_SYMMETRIC` = simetría de la contribución a la tangente del elemento. Parámetros físicos típicos: `sigma_t0` (Pa), `G_f` (N/m), `K_e` (Pa/m). No tienen densidad — la inercia es del bulk.

---

## CohesiveDamageIsotropic — daño cohesivo isótropo Modo-I con softening lineal/exponencial

- **Modelo**: daño escalar `ω ∈ [0, 1]` sobre el salto. Tracción `t_n = (1 − ω)·K_e·[[u_n]]` en dirección normal, `t_s = 0` en tangencial (Modo-I puro). Activación tipo Rankine (`f = ⟨[[u_n]]⟩ − κ`); historial monótono `κ_{n+1} = max(κ_n, ⟨[[u_n]]⟩)` con `κ_0 = σ_t0/K_e`. Energía de fractura `G_f` cierra la curva analíticamente:
  - **Lineal**: `T_soft(κ) = σ_t0·(w_c − κ)/(w_c − κ_0)` con apertura crítica `w_c = 2·G_f/σ_t0`. `ω(κ) = 1 − σ_t0·(w_c − κ)/[K_e·κ·(w_c − κ_0)]`. `ω = 1` en `κ ≥ w_c`.
  - **Exponencial**: `T_soft(κ) = σ_t0·exp[−σ_t0·(κ − κ_0)/H]` con `H = G_f − σ_t0·κ_0/2`. `ω(κ) = 1 − T_soft(κ)/(K_e·κ)`. Asintótica a 1.
- **Tangente**: simétrica por construcción (rank-1 sobre `n⊗n`). Algorítmica consistente en carga activa: `T_tan = K_e·[(1−ω) − [[u_n]]·dω/dκ]·(n⊗n)`. Secante reducida `(1−ω)·K_e·(n⊗n)` en descarga y sin daño. Cap por `DAMAGE_MAX` aplicado **sólo a la rigidez tangente** cuando `ω ≥ DAMAGE_MAX`; el valor de `ω` reportado y la tracción usan el valor físico (puede llegar a 1.0 exacto) — esto distingue al material de `IsotropicDamage2D`, donde `K_e ≡ E` permite truncar ω sin sesgo.
- **JUMP_DIM**: 2 · **PRIMARY_STATE_VAR**: `'damage'` · **IS_SYMMETRIC**: `True`.
- **Parámetros**: `sigma_t0` (resistencia a tracción, Pa), `G_f` (energía de fractura, N/m), `K_e` (rigidez del salto / penalty, Pa/m; sin default automático, ver §12 de la spec para guía `K_e ≈ 10·E_bulk/ℓ_c`), `softening ∈ {'linear', 'exponential'}`.
- **Variables internas**: `kappa` (historial del salto equivalente, [m]), `damage` (ω).
- **Validación energética**: por construcción `∫_0^{w_c} t·d[[u_n]] = G_f` (lineal) o `∫_0^∞ ≈ G_f` (exponencial); cubierto por tests con cuadratura trapezoidal.
- **Limitaciones declaradas** (`out_of_scope` en la spec): sólo Modo-I (mixto I–II diferido a fase G del ADR 0010), sin contacto unilateral en compresión (la grieta dañada transmite compresión con rigidez `(1−ω)·K_e`, no rígida), sin anisotropía del daño, sin acoplamiento viscoso/cíclico, sin regularización para mesh-objectivity (se aborda a nivel del elemento `CST_Embedded2D`, no del material).
- **Compatible con**: pendiente de fase 2 del ADR 0010 (elemento `CST_Embedded2D`). En aislamiento, testeable con historial prescrito de `[[u]]`.
- **Referencia**: ver `docs/specs/CohesiveDamageIsotropic.md`. Retama (2010) Cap. 3; Hillerborg, Modéer & Petersson (1976); Simó & Ju (1987).
- **Archivo**: [fenix/cohesive_materials/damage_isotropic.py](fenix/cohesive_materials/damage_isotropic.py)

---

## Cómo añadir un material nuevo

`/fenix-new material <Name>` — genera archivo en `fenix/materials/`, decorador `@MaterialRegistry.register`, esqueleto de test.
Declarar **`STRAIN_DIM`** y, si tiene historia, **`PRIMARY_STATE_VAR`** (la variable que aparecerá en el VTK).
Tras implementar el modelo constitutivo, **añadir una entrada a este catálogo** siguiendo el formato de arriba.

Para un material **cohesivo** nuevo: el patrón es análogo pero el archivo vive en `fenix/cohesive_materials/`, la clase hereda de `CohesiveMaterial` (no `Material`) y se registra con `@CohesiveMaterialRegistry.register`. Declarar `JUMP_DIM`, `PRIMARY_STATE_VAR` y `IS_SYMMETRIC`. Añadir la entrada en la sección "Materiales cohesivos" de este mismo catálogo.
