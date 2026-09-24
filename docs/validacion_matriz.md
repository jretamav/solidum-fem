# Matriz de validación — Solidum FEM

> Artefacto vivo de la Fase A del plan de validación y verificación (sesión 2026-05-19). Cruza cada componente con sus criterios `acceptance` declarados en la spec y el tipo de cobertura real en `tests/`.
>
> ⚠️ **Alcance (revisado 2026-09-23).** Esta matriz cubre los **36 componentes** que existían en mayo de 2026. Faltan 15 especificados después, cuya cobertura no se ha clasificado aquí criterio por criterio: `Elastic3D`, `VonMises3D`, `DruckerPrager3D`, `IsotropicDamage3D`, `Hex8`, `Hex20`, `Hex27`, `Tet4`, `Tet10` (Etapa 7 y sub-etapas A.bis/A.ter), `ThermalConduction`, `Quad4Thermal`, `Hex8Thermal`, `ThetaMethodSolver` (Etapa 8), `Orthotropic2D` (ADR 0013) y `DissipationArcLengthSolver`. Para ellos, la referencia vigente son los criterios `acceptance` de cada spec en `docs/specs/`, el recuento de specs validadas de `docs/STATUS.md` y los benchmarks de `tests/validation/README.md`. Ponerla al día exige clasificar cada criterio (analítica, benchmark, sanidad): es una auditoría de validación pendiente, no una actualización de texto. *Añadido 2026-09-24*: falta también `IndirectDisplacementSolver` (especificado el 2026-09-23), así que son 16.
>
> **Módulo de usuario (ADR 0020, 2026-09-24).** Dos de los 36 componentes, `CST_Embedded2D` y `CohesiveDamageIsotropic`, ya no son del programa principal: forman el módulo de usuario `discontinuities`, con specs en [`docs/user/discontinuities/specs/`](user/discontinuities/specs/) y tests en [`tests/user/discontinuities/`](../tests/user/discontinuities/), renombrados con el prefijo `test_disc_`. Sus filas siguen abajo con los nombres de test actuales; la clasificación de criterios es la de mayo.

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
| Truss2DCorot | elem | 3 | 2 | 1 | — | — | test_truss, test_snap_through_corot (von Mises 2-bar vs solución cerrada, 2026-05-19) |
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
| CST_Embedded2D (módulo de usuario) | elem | 7 | 5 | — | 2 | — | test_disc_cst_embedded, test_disc_cst_embedded_integration |
| Elastic1D | mat | 5 | 5 | — | — | — | test_truss, test_cable_material (implícito) |
| Elastic2D | mat | 6 | 6 | — | — | — | test_solid_2d (implícito) |
| CableMaterial1D | mat | 4 | 4 | — | — | — | test_cable_material |
| VonMises2D | mat | 16 | 8 | 4 | 4 | — | test_materials_unit, test_solid_2d_plasticity |
| Elastoplastic1D | mat | 6 | 6 | — | — | — | test_materials_unit, test_truss |
| IsotropicDamage1D | mat | 6 | 6 | — | — | — | test_materials_unit |
| IsotropicDamage2D | mat | 8 | 5 | — | 3 | — | test_materials_unit, test_solid_2d_damage |
| DruckerPrager2D | mat | 10 | 8 | — | 2 | — | test_solid_2d_drucker_prager (DP-1/DP-2/DP-2bis analíticos añadidos 2026-05-19) |
| CohesiveDamageIsotropic (módulo de usuario) | mat | 6 | 4 | — | 2 | — | test_disc_cohesive_damage_isotropic |
| LinearSolver | solv | 3 | 3 | — | — | — | test_integration |
| NonlinearSolver | solv | 6 | 5 | 1 | — | — | test_integration, test_solid_2d_plasticity, test_solver_robustness |
| ArcLengthSolver | solv | 3 | 2 | 1 | — | — | test_integration, test_snap_through_corot (curva carga-desplazamiento von Mises 2-bar punto a punto, 2026-05-19) |
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
2. **CohesiveDamageIsotropic** — sin pipeline completo de fractura con `CST_Embedded2D` ni medida de G_F. *Actualización 2026-09-23*: cubierto en parte por [`test_disc_embedded_uniaxial_softening.py`](../tests/user/discontinuities/test_disc_embedded_uniaxial_softening.py) (entonces `tests/validation/test_embedded_uniaxial_softening.py`): tracción uniaxial con `NonlinearSolver` a lo largo de la rama de ablandamiento frente a la solución exacta (lineal y exponencial) y, en el lineal, hasta la separación completa con el trabajo externo igual a `G_F·W·t_h` al 1 %. La clasificación de criterios de la tabla no se ha rehecho.
3. **Corotacionales**: `Frame2DEulerCorot` **CERRADO 2026-05-19 (Tanda 1)** (Bathe-Bolourchi cuarto + medio círculo). `Truss2DCorot` **CERRADO 2026-05-19 (Tanda 3)**: snap-through clásico del shallow von Mises 2-bar (h/L=0.1) con `ArcLengthSolver`; cada punto convergido (u_y, λ·F_ref) trazado por el solver se compara contra `P(w) = 2·E·A·(L_0 − L_d)/L_0 · (h − w)/L_d` (engineering strain) con tolerancia absoluta 0.01% de F_ref. Cobertura adicional: la malla de pasos cubre la rama elástica, la rama inestable (entorno del pico) y la rama invertida con λ < 0. Pendientes: `Truss3DCorot` (la formulación 3D del proyector perpendicular es plumbing del 2D, validación analítica heredada); `Cable2D/3DCorot` (no aplicable al snap-through — son unilaterales y no soportan compresión).
4. **CST_Embedded2D** — integración multi-elemento y modo II sólo en sanidad.
5. ~~**ResponseSpectrumSolver**~~ — **CERRADO 2026-05-19 (Tanda 2)**: cadena 4-truss Laplaciano 1D con autovalores cerrados ω_n² = 2 − 2·cos(nπ/4), valida ω, γ y SRSS contra fórmula analítica.
6. ~~**HHTSolver / NewtonHHTSolver**~~ — **CERRADO 2026-05-19 (Tanda 2)** para HHTSolver: matriz de amplificación 3×3 analítica para la convención Solidum (Hughes 1987 §9.3) contrastada contra contracción medida a Ω moderado + asíntota ρ_∞ = (1+α)/(1−α). NewtonHHTSolver hereda la verificación lineal. **CentralDifferenceSolver CERRADO 2026-05-19 (Tanda 2)**: frontera ``Δt_crit = 2/ω_max`` validada en 1-DOF y en cadena 2-DOF con autovalores cerrados ω_max=√3.
7. ~~**Quad8 / Quad9 / Tri6**~~ — **CERRADO 2026-05-19**: Quad8/Quad9 en Tanda 1 (Cook's 4×4 + refinamiento monótono). Tri6 en Tanda 2 (mismo trapezoide triangulado, diagonal c1→c3 con center node compartido). Locking volumétrico ν → 0.5 **CERRADO 2026-05-19 (Tanda 4)**: limitación arquitectural ahora blindada por test que demuestra que `Quad4` plane strain pierde >50% de deflexión al pasar de ν=0.3 a ν=0.4999, mientras que `Quad8` mantiene un ratio >1.5× mejor (lockea, pero menos). Ver `test_volumetric_locking.py`.

## Limitaciones arquitecturales documentadas (Fase D)

- ~~Locking volumétrico en sólidos 2D para ν → 0.5~~ ✅ **Blindado 2026-05-19**: `test_volumetric_locking.py` documenta el colapso del Quad4 y la mitigación parcial del Quad8 sobre cantilever esbelto. Mitigaciones B-bar / mixed siguen diferidas hasta caso de uso real.
- Apex return de DP en zona muy traccionante: sin smoothing si oscilara entre ramas.
- Newmark / HHT sin verificación de estabilidad analítica del esquema (sólo benchmarks numéricos).
