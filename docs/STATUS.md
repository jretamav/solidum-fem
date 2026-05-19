# Estado de Fenix FEM

> Foto de una página del estado actual. **Se regenera o actualiza tras cada PR significativa** (cierre de etapa, nuevo componente validado, nuevo ADR aceptado, deuda saldada).
>
> Para visión por etapas y bifurcaciones: ver [ROADMAP.md](ROADMAP.md). Para combinaciones validadas: ver [MATRIZ.md](MATRIZ.md). Para arranque en frío: ver [ONBOARDING.md](ONBOARDING.md).

---

## Métricas

| Indicador | Valor |
|---|---|
| **Tests** | 703 pasan / 5 skipped / 0 fallos (708 colectados, +20 subtests) — Fase B validación 2026-05-19: Tanda 1 (8 benchmarks), Tanda 2 (14 benchmarks dinámicos), Tanda 3 (2: snap-through), Tanda 4 (3: locking ν→0.5), Tanda 5 (6: internos NonlinearSolver), Tanda 6 (18: consistent tangent FD + energía Newmark + YAML e2e), Tanda 7 (34: verificaciones canónicas), Tanda 8 (6: plasticidad cíclica), Tanda 9 (8: profundización física materiales), Tanda 10 (1: momento lineal CentralDifference) y Tanda 11 (1: momento angular CentralDifference) |
| **Elementos** | 16 (10 estructurales 1D + 5 sólidos 2D + 1 sólido 2D con discontinuidad embebida) |
| **Materiales** | 9 (8 continuos + 1 cohesivo traction-jump) |
| **Solvers** | 11 (3 estáticos + 1 modal + 5 transitorios + 1 armónico + 1 espectral) |
| **ADRs aceptados** | 11 (0001–0011) |
| **Specs `validated`** | 29 |
| **Etapas cerradas** | 6 completas (Etapa 4 = ADR 0009 fases 1, 3, 4; Etapa 5 = ADR 0010 embedded; Etapa 6 = cierre completo del ADR 0009 con fases 2, 5, 6, 7 + HHT-α + reglas C y D) |
| **Auditoría global** | Conducida 2026-05-18 — 51 hallazgos. 41 cerrados, 6 diferidos con rationale, 4 diferidos sin acción. Informe + addendum en [docs/auditorias/auditoria_global_2026-05-18.md](auditorias/auditoria_global_2026-05-18.md). Saneamiento completo en dos sesiones (19-mayo). |

---

## Capacidad de análisis hoy

**Análisis estáticos**
- Lineal (`LinearSolver`).
- No lineal con control de carga (`NonlinearSolver`, Newton-Raphson incremental con paso adaptativo). Atraviesa puntos límite suaves con bisección + cinemática corotacional (hallazgo empírico de la auditoría 2026-05-18); no atraviesa snap-back con `du/dλ < 0`.
- Trazado de curvas de equilibrio con snap-through/snap-back (`ArcLengthSolver`, Crisfield cilíndrico). Atraviesa softening de daño bulk continuo (1D/2D); el límite del solver para softening severo aplica solo al penalty cohesivo stiff del embedded.

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

**Materiales cubiertos**
- Elásticos: `Elastic1D`, `Elastic2D` (plane stress + plane strain).
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

- **Sólidos 3D**: pendiente (`Hex8`, `Tet4`, etc.).
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
| 1 | `internal_forces` devuelve `None` en sólidos 2D (ADR 0002 incompleto). | Cuando entren sólidos 3D, post-proceso avanzado o consumidor externo que pida `ElementForces` para sólidos. |
| 2 | `FiberSection` para plasticidad por flexión en frames. | Si la próxima etapa se decide por la opción E (Mohr-Coulomb + FiberSection). |
| ~~3~~ | ~~Reglas de disparo C y D~~ ✅ **Ambas aplicadas 2026-05-18**: C con `PIPELINE_KIND` declarativo (commit `e66473d`); D con `fenix/math/modal_response.py` agrupando `free_vibration` + funciones de combinación espectral (commit final D4). | — |
| 4 | Solver para softening severo con embedded discontinuity (fase 4 ADR 0010 + curva descarga Van Vliet). | Cuando se priorice un mini-ADR de "solvers para softening": dissipation arc-length, control CMOD/CTOD del crack, sign-of-pivot tracking. |

Ninguno de los tres bloquea el avance. Todos están documentados con su contexto en la memoria del proyecto o en [MATRIZ.md](MATRIZ.md), y no requieren acción proactiva.

**Deuda saldada 2026-05-14**: cobertura asimétrica de sólidos 2D no-Quad4 con materiales no lineales — cerrada con [`test_solid_2d_nonlinear_higher_order.py`](../tests/test_solid_2d_nonlinear_higher_order.py) (16 tests verdes a la primera; no había bug latente).

**Etapa cerrada 2026-05-18**: discontinuidades interiores embebidas (Etapa 5, [ADR 0010](adr/0010-discontinuidades-interiores-embebidas.md)). Subsistema completo de fractura computacional vía embedded discontinuity discreta, fiel a Retama 2010: material cohesivo `CohesiveDamageIsotropic`, elemento `CST_Embedded2D` con condensación local, validación numérica del `l_d` del Cap. 6 y bring-up end-to-end. Bug arquitectural cazado durante el bring-up (factor `thickness` faltante en términos cohesivos, invisible con tests de espesor unitario) documentado en el ADR §"Caveats y lecciones aprendidas" y blindado con cuatro tests de regresión. Fases 4 (benchmark Van Vliet faithful) y 5 (manuales explícitos) diferidas — modelo físico verificado por curva analítica 1D, límite en solver no en formulación.

**Rama de trabajo cerrada 2026-05-18 (segunda)**: robustez de solvers no lineales ([ADR 0011](adr/0011-robustez-newton-line-search.md)). Auditoría de 11 tests ([fase A](auditorias/solvers_robustez_fase_A.md)) que mapeó el comportamiento real de `NonlinearSolver`, `ArcLengthSolver` y `NewtonNewmarkSolver` en regímenes problemáticos; identificó como único modo de fallo no atribuible a limitación física la oscilación del Newton en Drucker-Prager perfectamente plástico con sobrecarga. Implementó (fase C) la infraestructura completa de globalización: helper de line search por descenso no monótono y jerarquía de excepciones tipadas (`SolverDivergedError` + 4 subclases con métricas y `hint`). Enmienda al ADR durante la implementación: `line_search=False` por defecto, ya que con full Newton-Raphson + tangente consistente la globalización es contraproducente (respaldado por Tabla 3.1 de de Borst & Sluys 1999 §3.4, añadido al repo en `docs/referencias/`). Hallazgos colaterales documentados en el catálogo de solvers y reflejados en las "Capacidades" y "Limitaciones" arriba.

**Sesión de saneamiento cerrada 2026-05-19**: post-auditoría 2026-05-18 ([informe + addendum](auditorias/auditoria_global_2026-05-18.md)). 23 hallazgos cerrados de los 51 listados en la primera vuelta: 3 críticos (todos en solvers — `singular_tangent_seen`, cache de factorización, doble resta `F_dir`), 4 altos (autodiscover recursivo, `_resolve_rayleigh` centralizado, doble contrato `internal_forces` documentado, tests dimensionales fail-fast), 12 medios + 3 cubiertos documentalmente. Suite +27 tests (558 → 585). APIs nuevas todas aditivas.

**Segunda sesión de housekeeping cerrada 2026-05-19**: 18 hallazgos adicionales cerrados de los 27 que quedaban abiertos tras la primera vuelta. Cubre las 7 specs faltantes (H-5.3 — `Elastic1D`, `Elastic2D`, `Elastoplastic1D`, `LinearSolver`, `NonlinearSolver`, `ArcLengthSolver`, `NewtonHHTSolver`), MATRIZ/ONBOARDING ya existentes (H-5.5, H-5.6 — limpieza de menciones obsoletas), nomenclatura "tensión"→"esfuerzo"/"tracción" (H-3.3), contrato Material trial/commit documentado (H-3.8), cap DAMAGE_MAX continuo vs cohesivo (H-3.7), test switching exacto CableMaterial1D (H-3.6), test tr(ε^p)=0 multi-paso plane stress (H-3.9), `step_callback` warning en run_yaml no-estático (H-1.8), Result "shallow-frozen" documentado (H-1.10), `ElementForces` movido a `fenix.core` (H-1.11), `MaterialRegistry.names()`/`get()` en yaml_parser (H-1.14), doble rotación Frame2D documentada (H-2.8), Frame3D warning vertical "una vez por sesión" (H-2.9), `yield_tol(alpha_old)` documentado en J2 PS (H-3.10), doble ensamblaje line_search dinámico documentado (H-4.9), `Iventalla.dat` borrado (H-4.10). Suite +17 tests + 3 (585 → 602 + 3 nuevos en batch 4). Total **41 hallazgos cerrados** de 51. Los 10 restantes: 6 diferidos con rationale técnico explícito (H-1.3/Harmonic-Spectrum, H-1.4, H-2.3, H-4.4, H-4.8, H-5.3 ya cerrado), 4 sin acción por su propia recomendación (H-1.12 "no actuar ahora", H-1.13 "monitorizar", H-1.15 sentinel `-1` en 49 archivos sin retorno, H-2.10 "cuando aparezca segundo bulk").

---

## Próximo hito

**Elección de Etapa 7** entre las opciones del catálogo de bifurcaciones original (A, B, C, E — la D se ejecutó como Etapa 6 el 2026-05-18):

- A. Sólidos 3D (`Hex8`, `Tet4`, …).
- B. Placas y láminas.
- C. Análisis térmico desacoplado.
- ~~D. Completar ADR 0009~~ ✅ ejecutada como Etapa 6 (2026-05-18) — fases 1-7 + variante HHT-α + reglas C y D aplicadas.
- E. Mohr-Coulomb 2D + `FiberSection`.

El argumentario completo de cada opción está en [ROADMAP.md](ROADMAP.md). La decisión la toma el usuario con base en la dirección que quiera dar al proyecto tras esta etapa.

---

## Cómo se regenera este documento

1. Actualizar la tabla de métricas contando ADRs, specs, componentes registrados.
2. Si entró un componente nuevo: añadir a "Capacidad de análisis hoy" en la sección correspondiente.
3. Si saldó deuda: tachar el item de la tabla con `~~strikethrough~~` y mover la línea bajo el cuerpo o eliminarla en la siguiente actualización.
4. Si cambió el próximo hito: reescribir la sección "Próximo hito".
5. Actualizar la fecha al pie.

---

*Última actualización: 2026-05-19 — Tanda 11: 1 test de conservación de momento angular en `CentralDifferenceSolver` (`test_central_difference_angular_momentum.py`). Cierra la pareja con la Tanda 10 (momento lineal). Triángulo equilátero libre de 3 truss (todos los DOFs libres → 3 modos rígidos + 3 elásticos). Excitación: rotación rígida `ω·k×r` + perturbación radial de respiración (no aporta a L). El integrador BLM conserva `L(t) = Σ M_ii·(x_i·v_iy − y_i·v_ix)` con drift relativo < 1e-10 sobre 200 pasos. Las fuerzas axiales del truss son centrales (parejas acción-reacción colineales con el eje), garantizando que las internas no aportan momento. Suite 702 → 703. Queda abierto un único hueco: cohesivo+embedded G_F end-to-end, bloqueado por mini-ADR previo de "solvers para softening severo".*
