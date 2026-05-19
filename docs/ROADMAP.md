# Roadmap de Fenix FEM

> Documento navegacional de alto nivel. **El detalle vive en ADRs** (`docs/adr/`) y **specs** (`docs/specs/`); este archivo es el índice temporal que las articula.
>
> **Para el usuario** — describe en una página el camino recorrido y las decisiones de bifurcación pendientes.
> **Para la IA** — en cada sesión, este es el primer documento que sitúa "dónde estamos" antes de auditar arquitectura o proponer trabajo nuevo.
>
> Convención: una etapa = un grupo coherente de capacidades nuevas que dejan el código en estado entregable (tests verdes, manuales actualizados, ejemplos representativos). Cerrar una etapa no implica cerrar todos los ADRs que abrió — es estado del proyecto, no contabilidad de tickets.
>
> **Nivel de granularidad** — este documento usa "etapa" para los hitos del proyecto entero. Los ADRs internos pueden subdividir su entrega en **fases internas** (terminología del propio ADR). Cuando este ROADMAP las cita, conserva la palabra "fase" tal como aparece en el ADR para no introducir traducción.

---

## Estado a fecha del último commit

Fenix resuelve hoy **estática lineal y no lineal** (material y geométrica) sobre **estructuras 1D (truss/cable/frame 2D y 3D)** y **sólidos 2D (Quad4/Quad8/Quad9/Tri3/Tri6)**, con un catálogo de materiales que cubre **elasticidad, plasticidad J2 (plane strain + plane stress), Drucker-Prager y daño isótropo (1D y 2D)**. Sobre la misma maquinaria se ha abierto recientemente la línea **dinámica**: **análisis modal** y **transitorio Newmark (lineal y no lineal)**. **401 tests verdes + 5 skipped intencionales**.

> Para una foto más detallada del estado actual (métricas, deuda técnica, próximos hitos) ver `docs/STATUS.md` (pendiente de crear; ver §"Documentos complementarios" al final).

---

## Etapa 1 — Núcleo estático lineal · estructuras 1D · cerrada

**Capacidad añadida**: ensamblaje sparse, resolución `K·U = F` en un paso, condiciones de Dirichlet, cargas puntuales y peso propio sobre elementos 1D.

**Componentes**:
- Elementos: `Truss2D`, `Truss3D`, `Frame2DEuler`, `Frame2DTimoshenko`, `Frame3D`.
- Materiales: `Elastic1D`.
- Solver: `LinearSolver`.

**ADRs que la articulan**:
- [ADR 0001 — Mecanismos anti caja-negra](adr/0001-mecanismos-anti-caja-negra.md): contratos declarativos (`STRAIN_DIM`, `DOF_NAMES`) + auto-registro vía decoradores.
- [ADR 0002 — API de resultados para consumidores](adr/0002-api-de-resultados-para-consumidores.md): `SolveResult`, `ElementForces`, `internal_forces` como contrato canónico. **Deuda conocida**: incompleto para sólidos 2D (ver §"Deuda técnica" más abajo).
- [ADR 0004 — Imposición de Dirichlet](adr/0004-imposicion-de-condiciones-de-dirichlet.md): eliminación directa de DOFs prescritos.
- [ADR 0005 — Logging configurable](adr/0005-logging-configurable.md).

**Cierre**: tests analíticos por elemento (flecha de voladizo, axial puro, simetría) verdes; manuales describen los elementos.

---

## Etapa 2 — No linealidad geométrica y material 1D · sólidos 2D base · cerrada

**Capacidad añadida**: grandes desplazamientos y rotaciones en 1D (corotacional), comportamiento unilateral (cable), plasticidad y daño escalar 1D, sólidos 2D lineales, y los solvers no lineales que los explotan.

**Componentes**:
- Elementos 1D corotacionales: `Truss2DCorot`, `Truss3DCorot`, `Cable2DCorot`, `Cable3DCorot`, `Frame2DEulerCorot`.
- Sólidos 2D: `Quad4`, `Tri3`.
- Materiales: `Elastic2D`, `Elastoplastic1D`, `IsotropicDamage1D`, `CableMaterial1D`.
- Solvers: `NonlinearSolver` (Newton-Raphson incremental con paso adaptativo), `ArcLengthSolver` (Crisfield cilíndrico).

**ADRs que la articulan**:
- [ADR 0003 — Despachador de solver algebraico](adr/0003-despachador-de-solver-algebraico.md): elección automática Cholesky/LU según simetría/definición positiva, reaprovechamiento de factorización.

**Cierre**: tests de aceptación por elemento + benchmark de softening (daño 1D vía arc-length); manuales describen formulaciones corotacionales y solvers no lineales.

---

## Etapa 3 — Plasticidad y daño 2D · sólidos 2D cuadráticos · cerrada

**Capacidad añadida**: el catálogo de sólidos 2D se completa con cuadráticos y el catálogo de materiales 2D incorpora plasticidad J2 (las dos hipótesis cinemáticas usuales), plasticidad friccional Drucker-Prager y daño isótropo 2D con tangente consistente.

**Componentes**:
- Sólidos 2D cuadráticos: `Quad8`, `Quad9`, `Tri6` (los tres comparten la base interna `_HigherOrderSolid2D`).
- Materiales 2D no lineales: `VonMises2D` (plane strain + plane stress en kernels especializados), `DruckerPrager2D` (plane strain, con detección automática regular ↔ apex), `IsotropicDamage2D` (tangente algorítmica consistente, asimétrica en general).

**ADRs que la articulan**:
- [ADR 0006 — Tolerancias de admisibilidad](adr/0006-tolerancias-admisibilidad.md): patrón `f ≤ atol + rtol · escala` con escala física por material (`σ_y + H·α`, `E·κ_0`, `k(α)`).
- [ADR 0007 — Tolerancias de convergencia de solvers](adr/0007-tolerancias-convergencia-solvers.md): criterio dual desplazamiento + residuo con normalización adimensional.
- [ADR 0008 — `density` como propiedad del material](adr/0008-densidad-propiedad-del-material.md): habilita peso propio y prepara dinámica.

**Cierre**: tests unitarios + benchmarks pipeline (Quad4 + Drucker-Prager, viga elastoplástica, daño 2D con tangente consistente vs FD numérica) verdes; manuales actualizados con ecuaciones, return mappings y tangentes algorítmicas.

---

## Etapa 4 — Análisis modal y dinámica transitoria · cerrada (2026-05-18)

**Capacidad añadida**: subsistema modal/dinámico/espectral **completo**. Cubre los regímenes canónicos: modal por autovalores generalizados, transitorio implícito (Newmark/HHT lineal y no lineal), transitorio explícito (diferencias centradas) y análisis en frecuencia (respuesta armónica + espectro sísmico).

**Componentes** (11 solvers en total):

- **Modal**: `ModalSolver` (ARPACK shift-invert sobre `K·φ = ω²M·φ`).
- **Transitorio implícito lineal**: `NewmarkSolver` (β-γ parametrizables, Rayleigh con calibración modal); variante `HHTSolver` con disipación numérica controlada.
- **Transitorio implícito no lineal**: `NewtonNewmarkSolver`, `NewtonHHTSolver` (Newton-Raphson dentro de cada paso, tangente algorítmica consistente, Rayleigh con `K_0`).
- **Transitorio explícito**: `CentralDifferenceSolver` (leapfrog Belytschko-Liu-Moran; lineal y no lineal con parámetro `nonlinear`; estabilidad condicional CFL con detección a posteriori).
- **Frecuencia**: `HarmonicSolver` (aritmética compleja sobre `(−ω²M + iωC + K)·û = F̂` con barrido lineal/log/explícito).
- **Espectral / sísmico**: `ResponseSpectrumSolver` (combinación SRSS/CQC sobre espectro `S_d(ω)` o `S_a(ω)`).

**Infraestructura compartida**:

- **Mass lumping** (`fenix/math/mass_lumping.py::lump_hrz`): HRZ canónico para sólidos isoparamétricos + fórmula nodal directa para frames. Disponible en todos los elementos vía `compute_mass_matrix(lumping="lumped")`.
- **Cómputos sobre modos** (`fenix/math/modal_response.py`): `free_vibration`, `participation_factors`, `response_spectrum_srss`, `response_spectrum_cqc`, helpers de espectros (`spectrum_from_sa`, `spectrum_tabulated`). `ModalResult.free_vibration` queda como wrapper delgado (regla D de auditoría aplicada).
- **Dispatch declarativo** en `entry.py::run_yaml` por atributo de clase `PIPELINE_KIND` ∈ `{"static","modal","transient","harmonic","spectrum"}` (regla C de auditoría aplicada). Solvers no clásicos futuros no requieren tocar `entry.py`.

**Resultados** (todos inmutables): `SolveResult`, `ModalResult`, `TransientResult`, `HarmonicResult`, `ResponseSpectrumResult`.

**ADR que la articula**:
- [ADR 0009 — Análisis modal y dinámico](adr/0009-analisis-modal-y-dinamico.md). Hoja de ruta en 7 fases **todas entregadas** + variante HHT-α + reglas C y D de refactor arquitectural saldadas.

**Validación**: ~100 tests específicos del subsistema (modal, modal_catalog, newmark, newmark_nonlinear, hht, mass_lumping, central_difference, harmonic, response_spectrum). Benchmarks: oscilador 1 GDL contra solución analítica, viga Bernoulli-Euler biapoyada, recuperación exacta del caso lineal en solvers no lineales, función de transferencia armónica analítica, CQC ↔ SRSS con modos cercanos vs separados.

---

## Etapa 5 — Discontinuidades interiores embebidas · cerrada (2026-05-18)

**Capacidad entregada**: fractura computacional vía **embedded discontinuity en aproximación discreta** (Retama 2010). Primera línea de fractura del proyecto. Subsistema completo:

- Familia paralela `CohesiveMaterial` (traction-jump `t = T[[u]]`, parámetros `σ_t0`, `G_F`) con `CohesiveDamageIsotropic` validado (Mode-I, daño isótropo, softening lineal/exponencial).
- Elemento `CST_Embedded2D` con **DOFs enriquecidos a nivel elemental** y **condensación estática local** (formulación KOS, Caps. 2/5/7 de la tesis).
- Aportación original `l_d = (A/h)·cos(θ−α)` (Cap. 6) implementada y validada con experimento numérico de stress locking que la tesis no contiene.
- Bring-up de integración end-to-end con `Assembler` + `NonlinearSolver`/`ArcLengthSolver`.

**Por qué se eligió entre las 5 opciones A-E**: maximizó el conocimiento doctoral del usuario (autor de la formulación), construyó abstracciones reutilizables (futuros CZM clásicos, futuros incompatible-modes, eventual XFEM-style), e introdujo fractura computacional a Fenix — capacidad que ninguno de los caminos A-E aportaba.

**ADR que la articula**:
- [ADR 0010 — Discontinuidades interiores embebidas](adr/0010-discontinuidades-interiores-embebidas.md). Fases 1, 2, 3, 3b completadas; fases 4 (benchmark Van Vliet faithful contra la curva experimental) y 5 (regeneración explícita de manuales) **diferidas con justificación**.

**Lo entregado** (commits `c7e68de`, `2e85a70`, `737a0e1`, `e98aebb`, `e54b5ce`):
- Material `CohesiveDamageIsotropic` y elemento `CST_Embedded2D` registrados, con specs `validated` y entradas en catálogos.
- 56 + 4 tests específicos del subsistema; 420 tests verdes en la suite global.
- Bug arquitectural cazado durante el bring-up: factor `thickness` faltante en el balance bulk↔cohesivo del Newton local, invisible mientras todos los tests usasen `thickness = 1.0`. Documentado en el ADR §"Caveats y lecciones aprendidas" y blindado con cuatro tests de regresión `TestThicknessDimensionalConsistency`.

**Lo diferido y por qué**:
- **Fase 4 — benchmark Van Vliet faithful**: el modelo físico (elemento + cohesivo + condensación) está verificado en aislamiento por los tests del subsistema. Pero los solvers actuales (`NonlinearSolver` Newton-Raphson, `ArcLengthSolver` cilíndrico) no atraviesan la rama post-pico — la transición elástico→softening con penalty cohesivo stiff (`K_e ≈ 1e15 N/m³`) hace `K_t` casi singular en `κ_0`. Es dificultad numérica conocida del campo, no defecto del modelo. Retoma vía un **mini-ADR 0011 de "solvers para softening severo"** (dissipation arc-length de Gutiérrez/de Borst, control por CMOD/CTOD vía MPC con reacción correctamente leída en el cabezal, indicador de sign-of-pivot vía LDLᵀ con inercia de Sylvester) cuando se priorice. **Sin solver no se justifica reabrir el benchmark.**
- **Fase 5 — regeneración explícita de manuales**: catálogos ya tienen las entradas; los manuales se auto-regeneran del MD. La sección dedicada "Familia de fractura computacional" en el manual User/Architecture queda diferida con la fase 4 (sin la comparativa Van Vliet la sección no aporta).

### Opciones diferidas (las 5 originalmente identificadas como bifurcación)

Las opciones A-E identificadas anteriormente quedan **diferidas como etapas futuras**, sin orden cerrado. Se enumeran abreviadas; el ROADMAP previo a esta versión contenía la argumentación completa, recuperable por `git log`.

- **A. Sólidos 3D** (Hex8, Tet4, Tet10/Hex20/Hex27): extender materiales J2/daño/Drucker-Prager a Voigt 6D; abrir locking volumétrico y hourglassing. Pre-requisito de prácticamente todas las extensiones posteriores.
- **B. Placas y láminas** (Mindlin, Kirchhoff, MITC): formulación shell con cortante transversal y drilling DOFs.
- **C. Análisis térmico desacoplado** (Laplaciano estacionario + transitorio): salto a problema escalar; abre la puerta a termomecánica acoplada.
- ~~**D. Completar ADR 0009**~~ **Cerrada 2026-05-18 como Etapa 6 in extenso** (HHT-α, mass lumping fase 2, diferencias centradas, harmonic, response spectrum + reglas C y D de auditoría arquitectural aplicadas).
- **E. Mohr-Coulomb 2D + FiberSection**: cierra dos huecos puntuales del catálogo 2D (geotecnia + plasticidad por flexión en frames).

Cuando la Etapa 5 se cierre, la decisión sobre cuál de las opciones A-E entra a continuación se retoma con argumentos de la dirección que tome el proyecto entonces.

---

## Etapa 6 — Cierre del subsistema dinámico (ADR 0009 completo) · cerrada (2026-05-18)

**Capacidad entregada** (en una sesión, cinco commits): el subsistema modal/dinámico de la Etapa 4 quedó parcial al cierre original (fases 1, 3, 4); esta etapa completa el ADR 0009 con las cuatro fases pendientes (2, 5, 6, 7) y la variante HHT-α, además de aplicar las dos reglas arquitecturales C y D de la auditoría 2026-05-13 que esperaban evento.

**Lo entregado**:
- **HHT-α** (`HHTSolver` + `NewtonHHTSolver` como variantes de la familia Newmark) con disipación numérica controlada en altas frecuencias. Spec corta tipo extensión (Reglas §4).
- **Fase 2 ADR 0009 — mass lumping**: `compute_mass_matrix(lumping="lumped")` operativo en todos los elementos (truss/cable, frames 2D/3D, sólidos 2D Tri3/Quad4/Tri6/Quad8/Quad9). Helper centralizado `fenix/math/mass_lumping.py::lump_hrz` (HRZ canónico) + fórmula nodal directa para frames.
- **Regla C de auditoría aplicada** (dispatch declarativo en `entry.py::run_yaml` por atributo `PIPELINE_KIND`) + **`CentralDifferenceSolver`** (fase 5 ADR 0009) — leapfrog Belytschko-Liu-Moran con `M⁻¹` trivial, lineal y no lineal en una sola clase con parámetro `nonlinear`.
- **`HarmonicSolver`** (fase 6) — respuesta forzada armónica en el dominio de la frecuencia con aritmética compleja, barrido configurable, `HarmonicResult` con métodos `.amplitude()` y `.phase()`.
- **Regla D de auditoría aplicada** (`fenix/math/modal_response.py` agrupa `free_vibration` (movido) y los nuevos `participation_factors`, `response_spectrum_srss`, `response_spectrum_cqc`, `spectrum_from_sa`, `spectrum_tabulated`) + **`ResponseSpectrumSolver`** (fase 7) — combinación modal SRSS/CQC para análisis sísmico contra espectros normativos.

**Por qué se eligió como Etapa 6 en lugar de A/B/C/E**: ningún otro camino aprovechaba el alineamiento natural — `CentralDifference` requería la fase 2 (lumping), `HarmonicSolver` requería la regla C (3ª rama no clásica), `ResponseSpectrumSolver` requería la regla D (2º método sobre `ModalResult`). Las cinco piezas se encadenan en una secuencia obligada que cierra el ADR 0009 en su totalidad. Cualquier orden distinto hubiera dejado deudas internas más caras de pagar luego.

**Commits**: `73742f4` (HHT-α), `c92bf4e` (mass lumping), `e66473d` (regla C + central differences), `47a1333` (harmonic), `5151d7d` (regla D + response spectrum).

**Validación**: 558 tests verdes en suite global. 19 + 10 + 11 + 19 + 19 = ~80 tests nuevos del subsistema (mass_lumping, central_difference, harmonic, response_spectrum + reorganización menor en test_modal). Benchmarks: oscilador 1 GDL contra solución analítica, función de transferencia `H(ω) = 1/(K-ω²M+iωC)`, recuperación SRSS = CQC con modos separados, CQC ≠ SRSS con modos cercanos.

**Nada diferido de esta etapa**: el ADR 0009 queda cerrado en su totalidad. Las únicas extensiones futuras (excitación sísmica multi-directional simultánea CQC3, Δt adaptativo, MPC en frecuencia, generalized-α) requieren caso de uso real específico — no son piezas del ADR pendientes.

---

## Etapa 7+ — Horizonte largo

Lo que el proyecto **previsiblemente** abrirá una vez consolidada la etapa 5, sin orden cerrado:

- **Análisis termomecánico acoplado** (si C entró antes).
- **Contacto mecánico** (penalización / Lagrangiano aumentado / mortar).
- **Grandes deformaciones** (lagrangiano total / actualizado en sólidos; corotacional 3D para frames; viscoplasticidad).
- **Materiales avanzados**: hiperelasticidad (Neo-Hooke, Mooney-Rivlin), plasticidad anisótropa, daño con regularización (gradient damage, phase-field), modelos para hormigón (Mazars, concrete damaged plasticity).
- **Optimización topológica y de forma** (SIMP, level-set).
- **Análisis estocástico / fiabilidad** (FORM, Monte Carlo sobre el modelo determinista).

Todos comparten un patrón: **incrementales sobre lo existente, no refactor estructural**. La arquitectura actual (contratos declarativos + registries + capa algebraica + caché de topología) está diseñada para soportarlos sin reabrir decisiones.

---

## Deuda técnica conocida

Items que no bloquean el avance pero conviene tener visibles:

- **ADR 0002 incompleto para sólidos 2D** ([memoria de proyecto](../../../.claude/projects/g--Mi-unidad-Proyectos-IA-fenix-fem/memory/project_adr_0002_incompleto_solidos.md)): `internal_forces` devuelve `None` en sólidos; `compute_internal_forces`/`compute_gauss_state` cubren la necesidad. Decisión sobre cerrar o replantear el ADR diferida hasta caso de uso real (post-proceso avanzado, sólidos 3D, consumidor externo).
- **`FiberSection` para frames no-lineales** ([memoria](../../../.claude/projects/g--Mi-unidad-Proyectos-IA-fenix-fem/memory/project_pendiente_fiber_section.md)): frames 2D/3D plastifican sólo en axial. Plasticidad por flexión espera caso de uso (entraría como componente de la etapa 5 si la opción E se elige).
- **Reglas de disparo C y D** ([memoria](../../../.claude/projects/g--Mi-unidad-Proyectos-IA-fenix-fem/memory/project_reglas_disparo_pendientes.md)): dos refactors arquitecturales que esperan eventos antes de ejecutarse.

---

## Documentos complementarios (pendientes de crear)

El presente ROADMAP es uno de cuatro documentos previstos para que la traceability del proyecto escale con el catálogo:

1. **`ROADMAP.md` (este archivo)**: visión por etapas, decisiones de bifurcación. **Cambia cuando se cierra una etapa o se decide una bifurcación.**
2. **`STATUS.md`** (pendiente): foto de una página del estado actual — métricas (#tests, #elementos, #materiales, #solvers, #ADRs, #specs), capacidades habilitadas, deuda técnica con prioridades, próximo hito. **Se regenera o actualiza tras cada PR significativa.**
3. **`MATRIZ.md`** (pendiente): tabla cruzada elemento × material × solver con celdas validadas por test, en blanco (combinación posible no testeada) o vetadas (combinación prohibida con razón). Sirve para detectar huecos del catálogo y para usuarios externos que quieren saber qué combinación es segura.
4. **`ONBOARDING.md`** (pendiente): documento de entrada para agente o humano que arranca sesión sin contexto previo. Indica qué leer primero (este ROADMAP, `Reglas.md`, último ADR), cómo correr tests, dónde está la memoria, qué patrones de cambio merecen ADR.

---

## Cómo se cierra una etapa

1. Todos los componentes nuevos tienen spec en `docs/specs/` con `status: validated`.
2. Tests verdes (unitarios + acceptance + smoke YAML representativo).
3. Catálogos (`docs/catalogo_*.md`) y manuales (`manuals/sources/*/`) actualizados.
4. ADR(s) de la fase en estado `Aceptado`.
5. Entrada en este ROADMAP movida a "cerrada" con resumen ejecutivo.
6. `STATUS.md` regenerado con la nueva métrica.

---

*Última actualización: 2026-05-18 — Cierre de Etapa 5 (discontinuidades embebidas, ADR 0010) **y** Etapa 6 (ADR 0009 completo: HHT-α + mass lumping + central differences + harmonic + response spectrum + reglas C y D aplicadas). Próxima decisión: Etapa 7 entre opciones A, B, C, E originalmente diferidas.*
