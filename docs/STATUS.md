# Estado de Fenix FEM

> Foto de una página del estado actual. **Se regenera o actualiza tras cada PR significativa** (cierre de etapa, nuevo componente validado, nuevo ADR aceptado, deuda saldada).
>
> Para visión por etapas y bifurcaciones: ver [ROADMAP.md](ROADMAP.md). Para combinaciones validadas: ver [MATRIZ.md](MATRIZ.md). Para arranque en frío: ver [ONBOARDING.md](ONBOARDING.md).

---

## Métricas

| Indicador | Valor |
|---|---|
| **Tests** | 804 pasan / 5 skipped / 0 fallos (+54 vs 750 baseline, +20 subtests) — **Etapa 7 cerrada 2026-05-19** (ADR 0012, sólidos 3D acotados a Hex8 + Tet4 + Elastic3D): 31 unitarios (`tests/test_solid_3d.py`), 10 modos rígidos (`tests/test_rigid_body_modes.py`), 3 validación cubo Lamé 3D + 2 MacNeal 3D (`tests/validation/`), 2 locking volumétrico 3D blindado. La campaña de validación externa anterior (Tandas 1-14, 181 tests) sigue verde. |
| **Elementos** | 18 (10 estructurales 1D + 5 sólidos 2D + 1 sólido 2D con discontinuidad embebida + **2 sólidos 3D**) |
| **Materiales** | 10 (9 continuos + 1 cohesivo traction-jump) |
| **Solvers** | 12 (4 estáticos + 1 modal + 5 transitorios + 1 armónico + 1 espectral) |
| **ADRs aceptados** | 12 (0001–0012) |
| **Specs `validated`** | 32 |
| **Etapas cerradas** | 7 completas (Etapa 4 = ADR 0009 fases 1, 3, 4; Etapa 5 = ADR 0010 embedded; Etapa 6 = cierre completo del ADR 0009 con fases 2, 5, 6, 7 + HHT-α + reglas C y D; **Etapa 7 = sólidos 3D acotados con ADR 0012**) |
| **Auditoría global** | Conducida 2026-05-18 — 51 hallazgos. 41 cerrados, 6 diferidos con rationale, 4 diferidos sin acción. Informe + addendum en [docs/auditorias/auditoria_global_2026-05-18.md](auditorias/auditoria_global_2026-05-18.md). Saneamiento completo en dos sesiones (19-mayo). |

---

## Capacidad de análisis hoy

**Análisis estáticos**
- Lineal (`LinearSolver`).
- No lineal con control de carga (`NonlinearSolver`, Newton-Raphson incremental con paso adaptativo). Atraviesa puntos límite suaves con bisección + cinemática corotacional (hallazgo empírico de la auditoría 2026-05-18); no atraviesa snap-back con `du/dλ < 0`.
- Trazado de curvas de equilibrio con snap-through/snap-back (`ArcLengthSolver`, Crisfield cilíndrico). Atraviesa softening de daño bulk continuo (1D/2D); el límite del solver para softening severo aplica solo al penalty cohesivo stiff del embedded.
- Trazado por **disipación de energía** (`DissipationArcLengthSolver`, Gutiérrez 2004 + switching automático cilíndrico↔disipación a la Verhoosel et al. 2009). Atraviesa softening de daño continuo (1D/2D) igual de bien que el padre cilíndrico, **con la ventaja** de que en problemas con disipación severa la restricción lineal por τ es más robusta numéricamente que la cuadrática por dl. Incluye sign-of-pivot tracking aproximado (signo del determinante). El caso cohesivo+embedded con penalty `K_e` stiff sigue diferido — requiere LDLᵀ Bunch-Kaufman real, deuda técnica documentada en la spec.

**Análisis dinámicos**
- Modal — autovalores generalizados con shift-invert ARPACK (`ModalSolver`).
- Transitorio lineal Newmark con amortiguamiento Rayleigh, β-γ parametrizables (`NewmarkSolver`); variante HHT-α (`HHTSolver`) con disipación numérica controlada (`-1/3 ≤ α ≤ 0`).
- Transitorio no lineal Newton-Newmark con Newton dentro de cada paso (`NewtonNewmarkSolver`); variante HHT-α no lineal (`NewtonHHTSolver`).
- Transitorio explícito diferencias centradas (`CentralDifferenceSolver`, ADR 0009 fase 5) — leapfrog Belytschko-Liu-Moran con masa lumped diagonal, lineal y no lineal en una sola clase con parámetro `nonlinear`. Estabilidad condicional CFL (detección a posteriori).
- Respuesta forzada armónica en el dominio de la frecuencia (`HarmonicSolver`, ADR 0009 fase 6) — barrido sobre `ω` resolviendo `(−ω²M + iωC + K)·û = F̂` con aritmética compleja. Barrido lineal/logarítmico/explícito. Devuelve `HarmonicResult` con amplitud y fase complejas; métodos `.amplitude()` y `.phase()`.
- Análisis sísmico por combinación modal espectral (`ResponseSpectrumSolver`, ADR 0009 fase 7) — combina las respuestas máximas de los primeros `n_modes` modos bajo un espectro `S_d(ω)` o `S_a(ω)` con SRSS (modos separados) o CQC (modos cercanos, Der Kiureghian 1980). Espectro como callable, tabulado o aceleración constante. Devuelve `ResponseSpectrumResult` con respuesta envolvente, contribución modal individual, factores de participación, masas efectivas, y método `.cumulative_effective_mass_ratio()` para diagnóstico.
- **Masa**: consistente (default) o lumped (ADR 0009 fase 2). Lumped vía HRZ canónico para sólidos isoparamétricos y vía lumping nodal directo para vigas/marcos (`ρAL/2` traslacional + `ρIL/2` rotacional). Diagonal en globales para todos los elementos excepto Frame3D en orientación oblicua (bloque-diagonal por nodo, limitación documentada estándar — Cook-Malkus-Plesha §11.4).
- **Dispatch declarativo**: `run_yaml` despacha por atributo de clase `PIPELINE_KIND` ∈ `{"static","modal","transient","harmonic","spectrum"}` (regla C de auditoría arquitectural 2026-05-13, aplicada 2026-05-18). Solvers nuevos no clásicos no requieren tocar `entry.py`.
- **Cómputo sobre modos centralizado**: `fenix/math/modal_response.py` agrupa `free_vibration`, `participation_factors`, `response_spectrum_srss`, `response_spectrum_cqc`, `spectrum_from_sa`, `spectrum_tabulated` (regla D de auditoría arquitectural 2026-05-13, aplicada 2026-05-18 al introducir el segundo método sobre `ModalResult`). `ModalResult.free_vibration` queda como wrapper delgado, preservando la API histórica.

**Geometrías cubiertas**
- Estructuras 1D: barras (`Truss2D/3D` lineales y corotacionales), cables (`Cable2D/3D` corotacionales), vigas (`Frame2DEuler`, `Frame2DTimoshenko`, `Frame2DEulerCorot`, `Frame3D`).
- Sólidos 2D: lineales (`Quad4`, `Tri3`) y cuadráticos (`Quad8`, `Quad9`, `Tri6`).
- **Sólidos 2D con discontinuidad interior embebida** (Etapa 5, ADR 0010): `CST_Embedded2D` — CST con cinemática KOS enriquecida, DOFs del salto condensados a nivel elemental, activación Rankine en `prepare_step` con `l_d = (A_e/h)·cos(θ−α)` (Cap. 6 Retama 2010).
- **Sólidos 3D lineales** (Etapa 7, ADR 0012): `Hex8` (hexaedro trilineal isoparamétrico, 8 puntos Gauss 2×2×2 default, 6 caras numeradas con normal saliente) y `Tet4` (tetraedro lineal CST 3D, 1 punto baricéntrico). Convención Voigt 6D del proyecto `[ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz]`. Locking volumétrico Hex8 + shear locking Tet4 declarados como limitaciones blindadas con tests, sin mitigación (política idéntica a Quad4/Tri3 en 2D). Cuadráticos 3D (`Hex20`, `Hex27`, `Tet10`) diferidos a sub-etapa posterior.

**Materiales cubiertos**
- Elásticos: `Elastic1D`, `Elastic2D` (plane stress + plane strain), **`Elastic3D`** (isótropo, sin variantes de hipótesis — Etapa 7).
- Cable unilateral: `CableMaterial1D`.
- Plasticidad J2: `Elastoplastic1D`, `VonMises2D` (kernels Numba especializados para plane strain y plane stress).
- Plasticidad friccional: `DruckerPrager2D` (plane strain, dos ramas regular/apex, plasticidad no asociada).
- Daño isótropo con softening exponencial: `IsotropicDamage1D`, `IsotropicDamage2D` (tangente algorítmica consistente, asimétrica en 2D).
- **Material cohesivo traction-jump** (familia paralela `CohesiveMaterial`, ADR 0010): `CohesiveDamageIsotropic` — Mode-I, daño isótropo con softening lineal/exponencial parametrizado por `σ_t0` y `G_F`. Penalty `K_e` separado de la rigidez del bulk.

**Infraestructura común**
- Auto-registro vía decoradores en `ElementRegistry`, `MaterialRegistry`, `SolverRegistry`.
- Despachador algebraico automático Cholesky/LU según simetría/SPD (ADR 0003).
- Caché de topología COO compartida entre `K`, `M`, masa, peso propio.
- Tolerancias adimensionales con escala física por material (ADR 0006) y convergencia dual desplazamiento + residuo (ADR 0007).
- Parser YAML con despacho automático estático / modal / transitorio.
- Tres manuales `.tex/.pdf` generados desde markdown: Reference (FF-MR), User (FF-MU), Architecture (FF-MA).
- Diagnóstico tipado de divergencia en `NonlinearSolver`/`NewtonNewmarkSolver` (ADR 0011): subclases de `SolverDivergedError` (`OscillatingNewtonError`, `SingularTangentError`, `LoadExceedsCapacityError`, `UnknownDivergenceError`) con métricas y `hint`. Line search por descenso no monótono (Grippo-Lampariello-Lucidi) disponible vía `line_search=True` — opt-in, default `False` por la justificación bibliográfica del ADR.

---

## Limitaciones declaradas (out-of-scope hoy)

- **Sólidos 3D cuadráticos** (`Hex20`, `Hex27`, `Tet10`): pendientes — los lineales `Hex8` y `Tet4` entraron en la Etapa 7. Los cuadráticos abren la base interna `_HigherOrderSolid3D` cuando aparezca el segundo caso real (regla de centralización post-dos-casos).
- **Materiales 3D no lineales** (`VonMises3D`, `DruckerPrager3D`, `IsotropicDamage3D`): pendientes — son el siguiente eslabón natural de la Etapa 7 (Voigt 6D ya canonicalizado en ADR 0012).
- **Placas y láminas**: pendiente.
- **Análisis térmico**: pendiente (decoupled → coupled).
- **Mohr-Coulomb 2D**, **FiberSection para frames no-lineales**: pendientes.
- **Contacto mecánico**: horizonte largo.
- **Grandes deformaciones en sólidos** (lagrangiano total/actualizado): horizonte largo.
- **Hiperelasticidad**, **plasticidad anisótropa**, **daño con regularización** (gradient damage, phase-field): horizonte largo.
- **Locking volumétrico** en `Quad4`/`Quad8` con ν → 0.5: declarado limitación y **blindado por test** (`test_volumetric_locking.py`, 2026-05-19): cantilever esbelto plane strain documenta el colapso del Quad4 (pierde >50% de deflexión al pasar ν=0.3→0.4999) y la mitigación parcial del Quad8 (lockea, pero >1.5× mejor que Q4). Mitigaciones B-bar/F-bar diferidas hasta caso de uso real con materiales incompresibles.
- **Plasticidad por flexión en frames**: hoy sólo plastifican axialmente (espera `FiberSection`).
- **Trazado completo de la rama post-pico con embedded discontinuity**: `NonlinearSolver` y `ArcLengthSolver` cilíndrico no atraviesan la transición elástico→softening **con penalty cohesivo stiff** (`K_t` casi singular en `κ_0`). El límite es específico al embedded — para softening de daño bulk continuo (1D/2D) el arc-length sí traza la rama post-pico (verificado en `tests/test_solver_robustness.py::test_arc_length_traverses_damage_softening`). El modelo físico del embedded está verificado en aislamiento; el límite es solver. Retoma vía mini-ADR de "solvers para softening" (dissipation arc-length, control CMOD/CTOD, sign-of-pivot tracking) cuando se priorice.
- **Modo mixto I-II en cohesivo**, **KSON no-simétrico**, **embedded sin tracking**, **3D embedded**, **orden superior**: fases F-J del ADR 0010 diferidas.

---

## Deuda técnica priorizada

| # | Item | Criterio de retoma |
|---|---|---|
| ~~1~~ | ~~`internal_forces` devuelve `None` en sólidos 2D (ADR 0002 incompleto).~~ ✅ **Cerrado 2026-05-19 por ADR 0012**: cierre por dominio explícito. `internal_forces` aplica solo a elementos estructurales 1D (truss, cable, frame 2D/3D); sólidos 2D y 3D exponen `compute_gauss_state(U)` como API canónica. Docstrings de `Element.internal_forces` y `ElementForces` actualizados con el dominio acotado. Decisión arquitectural, no deuda pendiente. | — |
| 2 | `FiberSection` para plasticidad por flexión en frames. | Si la próxima etapa se decide por la opción E (Mohr-Coulomb + FiberSection). |
| ~~3~~ | ~~Reglas de disparo C y D~~ ✅ **Ambas aplicadas 2026-05-18**: C con `PIPELINE_KIND` declarativo (commit `e66473d`); D con `fenix/math/modal_response.py` agrupando `free_vibration` + funciones de combinación espectral (commit final D4). | — |
| 4 | Solver para softening severo con embedded discontinuity (fase 4 ADR 0010 + curva descarga Van Vliet). | Spec [`DissipationArcLengthSolver`](specs/DissipationArcLengthSolver.md) **implementada y validada parcialmente 2026-05-19** — funciona para daño continuo 1D/2D (cierra parcialmente este item) pero **no para cohesivo+embedded con K_e stiff**: la activación discreta del Rankine produce un salto en F_int que el `sign` por producto escalar no maneja. Resolverlo completamente requiere LDLᵀ Bunch-Kaufman real para sign-of-pivot tracking exacto (item #7 abajo) o control por CMOD/CTOD del salto (spec separada). |
| 7 | LDLᵀ Bunch-Kaufman como backend algebraico para sign-of-pivot tracking exacto. | Cuando se quiera cerrar el caso cohesivo+embedded del item #4. Hoy `DissipationArcLengthSolver._negative_pivots` usa signo del determinante como aproximación; ArcLengthSolver tiene `_negative_pivots()` como placeholder retornando `None`. Un LDLᵀ verdadero (pypardiso, scikit-sparse) los habilitaría a ambos. |
| 5 | Stress recovery superconvergente (SPR de Zienkiewicz) en post-proceso. | Cuando NAFEMS LE1 u otros benchmarks con valor canónico en nodo de esquina deban alcanzar el valor citable con malla coarse. Hoy `tests/validation/test_nafems_le1.py` documenta σ_yy(D) = 90.1 (Q4 32×32) vs 92.7 canónico — diferencia atribuible al recovery sin SPR. |
| 6 | B-bar o F-bar para Quad4/Quad8 en plasticidad incompresible. | Cuando aparezca caso de uso real con material casi-incompresible. Hoy `test_volumetric_locking.py` documenta el locking. |

Ninguno bloquea el avance. Todos están documentados con su contexto en la memoria del proyecto o en [MATRIZ.md](MATRIZ.md). El item #4 tiene **spec en redacción** desde 2026-05-19, motivada por la campaña de validación externa.

**Deuda saldada 2026-05-14**: cobertura asimétrica de sólidos 2D no-Quad4 con materiales no lineales — cerrada con [`test_solid_2d_nonlinear_higher_order.py`](../tests/test_solid_2d_nonlinear_higher_order.py) (16 tests verdes a la primera; no había bug latente).

**Etapa cerrada 2026-05-18**: discontinuidades interiores embebidas (Etapa 5, [ADR 0010](adr/0010-discontinuidades-interiores-embebidas.md)). Subsistema completo de fractura computacional vía embedded discontinuity discreta, fiel a Retama 2010: material cohesivo `CohesiveDamageIsotropic`, elemento `CST_Embedded2D` con condensación local, validación numérica del `l_d` del Cap. 6 y bring-up end-to-end. Bug arquitectural cazado durante el bring-up (factor `thickness` faltante en términos cohesivos, invisible con tests de espesor unitario) documentado en el ADR §"Caveats y lecciones aprendidas" y blindado con cuatro tests de regresión. Fases 4 (benchmark Van Vliet faithful) y 5 (manuales explícitos) diferidas — modelo físico verificado por curva analítica 1D, límite en solver no en formulación.

**Rama de trabajo cerrada 2026-05-18 (segunda)**: robustez de solvers no lineales ([ADR 0011](adr/0011-robustez-newton-line-search.md)). Auditoría de 11 tests ([fase A](auditorias/solvers_robustez_fase_A.md)) que mapeó el comportamiento real de `NonlinearSolver`, `ArcLengthSolver` y `NewtonNewmarkSolver` en regímenes problemáticos; identificó como único modo de fallo no atribuible a limitación física la oscilación del Newton en Drucker-Prager perfectamente plástico con sobrecarga. Implementó (fase C) la infraestructura completa de globalización: helper de line search por descenso no monótono y jerarquía de excepciones tipadas (`SolverDivergedError` + 4 subclases con métricas y `hint`). Enmienda al ADR durante la implementación: `line_search=False` por defecto, ya que con full Newton-Raphson + tangente consistente la globalización es contraproducente (respaldado por Tabla 3.1 de de Borst & Sluys 1999 §3.4, añadido al repo en `docs/referencias/`). Hallazgos colaterales documentados en el catálogo de solvers y reflejados en las "Capacidades" y "Limitaciones" arriba.

**Sesión de saneamiento cerrada 2026-05-19**: post-auditoría 2026-05-18 ([informe + addendum](auditorias/auditoria_global_2026-05-18.md)). 23 hallazgos cerrados de los 51 listados en la primera vuelta: 3 críticos (todos en solvers — `singular_tangent_seen`, cache de factorización, doble resta `F_dir`), 4 altos (autodiscover recursivo, `_resolve_rayleigh` centralizado, doble contrato `internal_forces` documentado, tests dimensionales fail-fast), 12 medios + 3 cubiertos documentalmente. Suite +27 tests (558 → 585). APIs nuevas todas aditivas.

**Segunda sesión de housekeeping cerrada 2026-05-19**: 18 hallazgos adicionales cerrados de los 27 que quedaban abiertos tras la primera vuelta. Cubre las 7 specs faltantes (H-5.3 — `Elastic1D`, `Elastic2D`, `Elastoplastic1D`, `LinearSolver`, `NonlinearSolver`, `ArcLengthSolver`, `NewtonHHTSolver`), MATRIZ/ONBOARDING ya existentes (H-5.5, H-5.6 — limpieza de menciones obsoletas), nomenclatura "tensión"→"esfuerzo"/"tracción" (H-3.3), contrato Material trial/commit documentado (H-3.8), cap DAMAGE_MAX continuo vs cohesivo (H-3.7), test switching exacto CableMaterial1D (H-3.6), test tr(ε^p)=0 multi-paso plane stress (H-3.9), `step_callback` warning en run_yaml no-estático (H-1.8), Result "shallow-frozen" documentado (H-1.10), `ElementForces` movido a `fenix.core` (H-1.11), `MaterialRegistry.names()`/`get()` en yaml_parser (H-1.14), doble rotación Frame2D documentada (H-2.8), Frame3D warning vertical "una vez por sesión" (H-2.9), `yield_tol(alpha_old)` documentado en J2 PS (H-3.10), doble ensamblaje line_search dinámico documentado (H-4.9), `Iventalla.dat` borrado (H-4.10). Suite +17 tests + 3 (585 → 602 + 3 nuevos en batch 4). Total **41 hallazgos cerrados** de 51. Los 10 restantes: 6 diferidos con rationale técnico explícito (H-1.3/Harmonic-Spectrum, H-1.4, H-2.3, H-4.4, H-4.8, H-5.3 ya cerrado), 4 sin acción por su propia recomendación (H-1.12 "no actuar ahora", H-1.13 "monitorizar", H-1.15 sentinel `-1` en 49 archivos sin retorno, H-2.10 "cuando aparezca segundo bulk").

---

## Próximo hito

**Elección de Etapa 8** tras cerrar la Etapa 7 (sólidos 3D acotados con ADR 0012 el 2026-05-19). Opciones restantes del catálogo original + ramificaciones naturales de la Etapa 7:

- B. Placas y láminas.
- C. Análisis térmico desacoplado.
- E. Mohr-Coulomb 2D + `FiberSection`.
- **A.bis** (rama natural de A): materiales 3D no lineales (`VonMises3D`, `DruckerPrager3D`, `IsotropicDamage3D`) — extender plasticidad/daño a Voigt 6D sobre Hex8/Tet4 ya disponibles.
- **A.ter** (rama natural de A): elementos 3D cuadráticos (`Hex20`, `Hex27`, `Tet10`) — abrir base interna `_HigherOrderSolid3D` (regla de dos casos reales) y desbloquear benchmarks NAFEMS 3D con valor citable (LE10, LE2, LE3).

El argumentario completo de cada opción está en [ROADMAP.md](ROADMAP.md). La decisión la toma el usuario con base en la dirección que quiera dar al proyecto tras esta etapa.

---

## Cómo se regenera este documento

1. Actualizar la tabla de métricas contando ADRs, specs, componentes registrados.
2. Si entró un componente nuevo: añadir a "Capacidad de análisis hoy" en la sección correspondiente.
3. Si saldó deuda: tachar el item de la tabla con `~~strikethrough~~` y mover la línea bajo el cuerpo o eliminarla en la siguiente actualización.
4. Si cambió el próximo hito: reescribir la sección "Próximo hito".
5. Actualizar la fecha al pie.

---

*Última actualización: 2026-05-19 — **Etapa 7 cerrada: sólidos 3D acotados (ADR 0012)**. Alcance: Hex8 (hexaedro trilineal, 8 puntos Gauss 2×2×2, 6 caras con normal saliente, body load + face traction consistentes), Tet4 (CST 3D, 1 punto baricéntrico, masa consistente analítica), Elastic3D (isótropo, sin variantes plane stress/plane strain). Convención Voigt 6D fijada en `Reglas.md §5` y blindada por ADR. Cuadraturas 3D registradas (`hex_1x1x1`, `hex_2x2x2`, `hex_3x3x3`, `tet_1`). Cierre adicional: deuda técnica #1 (`internal_forces` en sólidos) por dominio explícito en ADR 0012 — sólidos exponen `compute_gauss_state`, estructurales 1D exponen `internal_forces`. Tests: +54 (suite 750 → 804): 31 unitarios + 10 modos rígidos 3D + 3 cubo Lamé 3D (tracción uniaxial + hidrostática exactos a precisión máquina) + 2 MacNeal-Harder 3D (locking documentado + convergencia h monótona) + 2 locking volumétrico 3D blindado. Specs `Elastic3D`, `Hex8`, `Tet4` validadas; entradas en catálogos de elementos y materiales; MATRIZ ampliada con columna Elastic3D y filas Hex8/Tet4.*

*Anterior 2026-05-19: **DissipationArcLengthSolver implementado** (spec `docs/specs/DissipationArcLengthSolver.md` status `implemented`, archivo `fenix/math/solvers/dissipation_arclength.py`, tests `tests/test_dissipation_arclength.py` 7/7 verde). Cierra parcialmente la deuda técnica #4 (solver para softening) — funciona para daño continuo bulk 1D/2D, con switching automático cilíndrico↔disipación (Gutiérrez 2004 + Verhoosel et al. 2009). El caso cohesivo+embedded con `K_e=1e13` queda diferido a una segunda iteración que requiera LDLᵀ Bunch-Kaufman como backend (nuevo item #7 de deuda técnica). Suite 741 → 750 verdes + 5 skipped.*

*Anterior 2026-05-19: **Campaña de validación externa Tandas 12-14** cerrada. Tres tandas, **5 benchmarks publicados** en `tests/validation/`, **38 tests** sumados a la suite (703 → 741): Lamé plane strain (12 tests, Q4/Tri3/Q8/Q9/Tri6, convergencia O(h²) Q8 verificada), NAFEMS LE1 elliptic membrane (10 tests, σ_yy(D)=90.1 vs 92.7 canónico con Q4 32×32 — error 2.8% atribuible al recovery sin SPR), MacNeal-Harder slender beam (8 tests, Frame Euler/Timoshenko exactos con 1 elemento, shear locking Q4 cuantificado), Bathe wave propagation con CentralDifference (4 tests, c_num=1.0017 vs c=1 analítico con N=200), Hill cylinder J2 perfecta (4 tests, contra solución analítica cerrada de Hill 1950, delimitación zona plástica/elástica verificada). La campaña queda registrada en `[[project_validacion_externa_tandas_12_14]]`. Los huecos restantes (NAFEMS LE10/LE2/LE3 3D, FV1/FV5 placas, FV2/FV12/FV32 cáscaras, G_F end-to-end, recovery LE1 canónico, locking volumétrico) **no se cierran con más tests** sino con componentes nuevos. El primero priorizado es **dissipation arc-length solver** (spec en redacción), que desbloquea fase 4 del ADR 0010 + validación cuantitativa de cohesivos en pipeline real.*
