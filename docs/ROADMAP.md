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

## Etapa 4 — Análisis modal y dinámica transitoria · parcialmente cerrada

**Capacidad añadida**: matriz de masa consistente en todos los elementos del catálogo, autovalor generalizado `K·φ = ω²M·φ`, integración temporal Newmark lineal y no lineal con amortiguamiento Rayleigh.

**Componentes**:
- Solvers: `ModalSolver` (ARPACK shift-invert), `NewmarkSolver` (transitorio lineal, β-γ parametrizables, Rayleigh con calibración modal), `NewtonNewmarkSolver` (Newton dentro de cada paso temporal, hereda predictores/correctores de `NewmarkSolver` — patrón de variante de Reglas §4).
- Resultados: `ModalResult`, `TransientResult` (con historiales `u/u̇/ü`).
- Despacho automático en `fenix.run_yaml` por `isinstance` sobre el tipo de solver.

**ADR que la articula**:
- [ADR 0009 — Análisis modal y dinámico](adr/0009-analisis-modal-y-dinamico.md). Hoja de ruta en 7 fases internas; entregadas las fases **1 (modal)**, **3 (Newmark lineal)** y **4 (Newmark no lineal)**.

**Diferido dentro de la etapa 4** (fases del ADR 0009 sin entregar):
- Fase 2: modal con masa lumped — se abre cuando entre la fase 5 (explícita).
- Fase 5: integración explícita por diferencias centradas.
- Fase 6: respuesta en frecuencia (steady-state harmonic).
- Fase 7: análisis espectral / sísmico (CQC, SRSS).

**Cierre parcial**: la línea dinámica está abierta y validada con benchmarks (barra empotrada-libre, viga Bernoulli-Euler, recuperación exacta del caso lineal por `NewtonNewmarkSolver`). Las fases diferidas no bloquean el avance: el resto del proyecto puede progresar mientras se decide cuándo retomarlas.

---

## Etapa 5 — Discontinuidades interiores embebidas · abierta (2026-05-13)

**Capacidad que se abre**: fractura computacional vía **embedded discontinuity en aproximación discreta** (Retama 2010). Primera línea de fractura del proyecto. Introduce dos abstracciones nuevas: una familia paralela `CohesiveMaterial` (traction-jump `t = T[[u]]`, parámetros `σ_t0`, `G_F`) y una subjerarquía de elementos con **DOFs enriquecidos a nivel elemental** y **condensación estática local**. La etapa replica fielmente la formulación KOS sobre CST 2D de la tesis del usuario (incluyendo la aportación original `l_d = (A/h)·cos(θ−α)` del Cap. 6) y la valida con los benchmarks Van Vliet (Cap. 8.1) y el experimento numérico de stress locking que la tesis no contiene.

**Por qué esta entre las 5 opciones inicialmente identificadas (A-E, ahora diferidas)**: maximiza el uso del conocimiento doctoral del usuario (autor de la formulación), construye abstracciones reutilizables (futuros CZM clásicos, futuros incompatible-modes, eventual XFEM-style), e introduce fractura computacional a Fenix — capacidad que ninguno de los caminos A-E aportaba.

**ADR que la articula**:
- [ADR 0010 — Discontinuidades interiores embebidas](adr/0010-discontinuidades-interiores-embebidas.md). 5 fases activas (material cohesivo, elemento CST embedded, validación Cap. 6, benchmark Van Vliet, catálogo + manuales) + 5 fases diferidas (sin-tracking, modo mixto, KSON, 3D, orden superior).

**Estado actual (2026-05-13)**:
- ADR 0010 aceptado.
- Spec `CohesiveDamageIsotropic` redactada como borrador en [docs/specs/CohesiveDamageIsotropic.md](specs/CohesiveDamageIsotropic.md); `status: draft`. Tres puntos abiertos pendientes de decisión del usuario antes de implementar (ver sección *Diálogo* de la spec).
- Aún sin código.

### Opciones diferidas (las 5 originalmente identificadas como bifurcación)

Las opciones A-E identificadas anteriormente quedan **diferidas como etapas futuras**, sin orden cerrado. Se enumeran abreviadas; el ROADMAP previo a esta versión contenía la argumentación completa, recuperable por `git log`.

- **A. Sólidos 3D** (Hex8, Tet4, Tet10/Hex20/Hex27): extender materiales J2/daño/Drucker-Prager a Voigt 6D; abrir locking volumétrico y hourglassing. Pre-requisito de prácticamente todas las extensiones posteriores.
- **B. Placas y láminas** (Mindlin, Kirchhoff, MITC): formulación shell con cortante transversal y drilling DOFs.
- **C. Análisis térmico desacoplado** (Laplaciano estacionario + transitorio): salto a problema escalar; abre la puerta a termomecánica acoplada.
- **D. Completar ADR 0009** (HHT-α, diferencias centradas, harmonic): cierra los huecos de la línea dinámica.
- **E. Mohr-Coulomb 2D + FiberSection**: cierra dos huecos puntuales del catálogo 2D (geotecnia + plasticidad por flexión en frames).

Cuando la Etapa 5 se cierre, la decisión sobre cuál de las opciones A-E entra a continuación se retoma con argumentos de la dirección que tome el proyecto entonces.

---

## Etapa 6+ — Horizonte largo

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

*Última actualización: 2026-05-13 — apertura de Etapa 5 (discontinuidades interiores embebidas, ADR 0010); las 5 opciones A-E originalmente identificadas quedan diferidas.*
