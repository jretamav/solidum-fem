# Matriz de validación — Fenix FEM

> Artefacto vivo de la Fase A del plan de validación y verificación (sesión 2026-05-19). Cruza cada componente con sus criterios `acceptance` declarados en la spec y el tipo de cobertura real en `tests/`.

**Leyenda de cobertura**

- **Analítica** — comparación contra solución cerrada (Euler, von Mises uniaxial, etc.).
- **Benchmark** — comparación contra resultado publicado o de otro código.
- **Sanidad** — verifica que no crashea, signos correctos, invariantes triviales (simetría, equilibrio numérico).
- **NO** — criterio sin cobertura.

---

## Tabla maestra

| Componente | Tipo | Criterios (n) | Analítica | Benchmark | Sanidad | NO | Tests asociados |
|---|---|---:|---:|---:|---:|---:|---|
| Truss2D | elem | 1 | 1 | — | — | — | test_truss, test_integration |
| Truss3D | elem | 1 | 1 | — | — | — | test_truss |
| Truss2DCorot | elem | 3 | 2 | — | 1 | — | test_truss |
| Truss3DCorot | elem | 3 | 2 | — | 1 | — | test_truss |
| Frame2DEuler | elem | 3 | 3 | — | — | — | test_frame |
| Frame2DTimoshenko | elem | 3 | 3 | — | — | — | test_frame |
| Frame2DEulerCorot | elem | 6 | 5 | — | 1 | — | test_frame (Bathe-Bolourchi cuarto + medio círculo añadidos 2026-05-19) |
| Frame3D | elem | 5 | 5 | — | — | — | test_frame3d |
| Cable2DCorot | elem | 4 | 3 | — | 1 | — | test_cable_elements |
| Cable3DCorot | elem | 4 | 3 | — | 1 | — | test_cable_elements |
| Quad4 | elem | 5 | 5 | — | — | — | test_solid_2d, test_patch_solid_2d |
| Quad8 | elem | 5 | 4 | 1 | — | — | test_higher_order_solid_2d, test_cooks_membrane (2026-05-19) |
| Quad9 | elem | 3 | 2 | 1 | — | — | test_higher_order_solid_2d, test_cooks_membrane (2026-05-19) |
| Tri3 | elem | 4 | 4 | — | — | — | test_solid_2d, test_patch_solid_2d |
| Tri6 | elem | 4 | 3 | 1 | — | — | test_higher_order_solid_2d, test_cooks_membrane (Bathe/Hughes 4×4 + refinamiento monótono, 2026-05-19) |
| CST_Embedded2D | elem | 7 | 5 | — | 2 | — | test_cst_embedded, test_cst_embedded_integration |
| Elastic1D | mat | 5 | 5 | — | — | — | test_truss, test_cable_material (implícito) |
| Elastic2D | mat | 6 | 6 | — | — | — | test_solid_2d (implícito) |
| CableMaterial1D | mat | 4 | 4 | — | — | — | test_cable_material |
| VonMises2D | mat | 16 | 8 | 4 | 4 | — | test_materials_unit, test_solid_2d_plasticity |
| Elastoplastic1D | mat | 6 | 6 | — | — | — | test_materials_unit, test_truss |
| IsotropicDamage1D | mat | 6 | 6 | — | — | — | test_materials_unit |
| IsotropicDamage2D | mat | 8 | 5 | — | 3 | — | test_materials_unit, test_solid_2d_damage |
| DruckerPrager2D | mat | 10 | 8 | — | 2 | — | test_solid_2d_drucker_prager (DP-1/DP-2/DP-2bis analíticos añadidos 2026-05-19) |
| CohesiveDamageIsotropic | mat | 6 | 4 | — | 2 | — | test_cohesive_damage_isotropic |
| LinearSolver | solv | 3 | 3 | — | — | — | test_integration |
| NonlinearSolver | solv | 6 | 5 | 1 | — | — | test_integration, test_solid_2d_plasticity, test_solver_robustness |
| ArcLengthSolver | solv | 3 | 2 | 1 | — | — | test_integration |
| NewmarkSolver | solv | 4 | 2 | — | 2 | — | test_newmark |
| HHTSolver | solv | 6 | 3 | 1 | 2 | — | test_hht (amortiguamiento numérico contra matriz de amplificación analítica + ρ_∞ asintótico, 2026-05-19) |
| NewtonNewmarkSolver | solv | 4 | 2 | — | 2 | — | test_newmark_nonlinear |
| NewtonHHTSolver | solv | 4 | 2 | — | 2 | — | test_hht |
| ModalSolver | solv | 3 | 1 | 1 | 1 | — | test_modal |
| HarmonicSolver | solv | 3 | 1 | 1 | 1 | — | test_harmonic |
| ResponseSpectrumSolver | solv | 3 | 1 | 1 | 1 | — | test_response_spectrum (Laplaciano 1D 3-DOF con autovalores cerrados, 2026-05-19) |
| CentralDifferenceSolver | solv | 5 | 3 | — | 2 | — | test_central_difference (CFL = 2/ω_max contrastado contra autovalores cerrados 1-DOF y 2-DOF, 2026-05-19) |

---

## Huecos prioritarios (Fase B)

Componentes con criterios físicamente delicados cubiertos sólo por sanidad, ordenados por criticidad:

1. ~~**DruckerPrager2D**~~ — **CERRADO 2026-05-19 (Tanda 1)**: DP-1 (onset confinado), DP-2 (invariante de flujo no asociado tr(ε_p)=3η_g·α), DP-2bis (apex biaxial) añadidos a `test_solid_2d_drucker_prager.py`. DP-3 (cilindro Hill) diferido a sesión propia.
2. **CohesiveDamageIsotropic** — sin pipeline completo de fractura con `CST_Embedded2D` ni medida de G_F.
3. **Corotacionales**: `Frame2DEulerCorot` **CERRADO 2026-05-19 (Tanda 1)** (Bathe-Bolourchi cuarto + medio círculo). Pendientes: Truss2D/3DCorot, Cable2D/3DCorot — benchmark de snap-through (von Mises truss) requiere arc-length, sesión propia.
4. **CST_Embedded2D** — integración multi-elemento y modo II sólo en sanidad.
5. ~~**ResponseSpectrumSolver**~~ — **CERRADO 2026-05-19 (Tanda 2)**: cadena 4-truss Laplaciano 1D con autovalores cerrados ω_n² = 2 − 2·cos(nπ/4), valida ω, γ y SRSS contra fórmula analítica.
6. ~~**HHTSolver / NewtonHHTSolver**~~ — **CERRADO 2026-05-19 (Tanda 2)** para HHTSolver: matriz de amplificación 3×3 analítica para la convención Fenix (Hughes 1987 §9.3) contrastada contra contracción medida a Ω moderado + asíntota ρ_∞ = (1+α)/(1−α). NewtonHHTSolver hereda la verificación lineal. **CentralDifferenceSolver CERRADO 2026-05-19 (Tanda 2)**: frontera ``Δt_crit = 2/ω_max`` validada en 1-DOF y en cadena 2-DOF con autovalores cerrados ω_max=√3.
7. ~~**Quad8 / Quad9 / Tri6**~~ — **CERRADO 2026-05-19**: Quad8/Quad9 en Tanda 1 (Cook's 4×4 + refinamiento monótono). Tri6 en Tanda 2 (mismo trapezoide triangulado, diagonal c1→c3 con center node compartido). Pendiente: test de locking volumétrico ν → 0.5 (limitación arquitectural conocida, no bug).

## Limitaciones arquitecturales no documentadas (Fase D)

- Locking volumétrico en sólidos 2D para ν → 0.5 (B-bar / mixed no implementados).
- Apex return de DP en zona muy traccionante: sin smoothing si oscilara entre ramas.
- Newmark / HHT sin verificación de estabilidad analítica del esquema (sólo benchmarks numéricos).
