# Roadmap de Solidum FEM

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

Solidum resuelve hoy **estática lineal y no lineal** (material y geométrica) sobre **estructuras 1D (truss/cable/frame 2D y 3D)**, **sólidos 2D (Quad4/Quad8/Quad9/Tri3/Tri6)** y **sólidos 3D lineales + cuadráticos (Hex8/Hex20/Hex27/Tet4/Tet10)**, con un catálogo de materiales que cubre **elasticidad, plasticidad J2 (1D/2D/3D), Drucker-Prager (2D/3D), daño isótropo (1D/2D/3D) y cohesivo traction-jump**. Sobre la misma maquinaria está abierta la línea **dinámica** completa: modal, transitorio implícito (Newmark/HHT lineal y no lineal), transitorio explícito (diferencias centradas), armónico en frecuencia y espectro sísmico. Desde la **Etapa 8** (2026-08-25) resuelve además **conducción de calor** —estacionaria y transitoria, 2D y 3D— sobre la misma infraestructura, sin acoplamiento con el campo mecánico. **1176 tests verdes + 8 skipped intencionales**.

**Próximo hito: elección de la Etapa 9.** Cerrada la opción C en su núcleo mínimo, quedan abiertas **B** (placas y láminas) y **E** (Mohr-Coulomb 2D + `FiberSection`), más la continuación de la propia línea térmica —convección, acoplamiento termomecánico o no linealidad térmica—. Es una decisión del usuario; el argumentario está en §"Opciones diferidas", donde se registran también **G** (elementos isogeométricos) y **H** (apoyos elásticos tipo Winkler), evaluadas el 2026-08-25 y no candidatas a Etapa 9.

> Para una foto más detallada del estado actual (métricas, deuda técnica, próximos hitos) ver [`docs/STATUS.md`](STATUS.md). Para combinaciones validadas: [`docs/MATRIZ.md`](MATRIZ.md). Para arranque en frío: [`docs/ONBOARDING.md`](ONBOARDING.md). Ver §"Documentos complementarios" al final para una guía del sistema.

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

- **Mass lumping** (`solidum/math/mass_lumping.py::lump_hrz`): HRZ canónico para sólidos isoparamétricos + fórmula nodal directa para frames. Disponible en todos los elementos vía `compute_mass_matrix(lumping="lumped")`.
- **Cómputos sobre modos** (`solidum/math/modal_response.py`): `free_vibration`, `participation_factors`, `response_spectrum_srss`, `response_spectrum_cqc`, helpers de espectros (`spectrum_from_sa`, `spectrum_tabulated`). `ModalResult.free_vibration` queda como wrapper delgado (regla D de auditoría aplicada).
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

**Por qué se eligió entre las 5 opciones A-E**: maximizó el conocimiento doctoral del usuario (autor de la formulación), construyó abstracciones reutilizables (futuros CZM clásicos, futuros incompatible-modes, eventual XFEM-style), e introdujo fractura computacional a Solidum — capacidad que ninguno de los caminos A-E aportaba.

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

- ~~**A. Sólidos 3D**~~ **Cerrada**: Etapa 7 (Hex8, Tet4, Elastic3D — ADR 0012, 2026-05-19), sub-etapa **A.bis** (materiales 3D no lineales, 2026-05-21) y sub-etapa **A.ter** (Hex20, Hex27, Tet10 + validación NAFEMS LE10 y Lamé 3D, 2026-05-27). Era pre-requisito de casi todas las extensiones posteriores, y ya no bloquea ninguna.
- **B. Placas y láminas** (Mindlin, Kirchhoff, MITC): formulación shell con cortante transversal y drilling DOFs. Desbloquearía además los benchmarks NAFEMS de placas y cáscaras (LE3 hemisferio, FV1/FV5, FV2/FV12/FV32), hoy inalcanzables por falta del componente.
- ~~**C. Análisis térmico desacoplado**~~ **Cerrada en su núcleo mínimo**: Etapa 8 (2026-08-25) — conducción de Fourier pura, estacionaria y transitoria, 2D y 3D, con Dirichlet y Neumann. El salto a problema escalar confirmó que la infraestructura es agnóstica al campo (el estacionario no necesitó solver nuevo). **La termomecánica acoplada de `Reglas.md §0` sigue sin existir**: era el C2 explícitamente excluido del alcance, y es ahora la continuación natural de esta línea junto con la convección (Robin).
- ~~**D. Completar ADR 0009**~~ **Cerrada 2026-05-18 como Etapa 6 in extenso** (HHT-α, mass lumping fase 2, diferencias centradas, harmonic, response spectrum + reglas C y D de auditoría arquitectural aplicadas).
- **E. Mohr-Coulomb 2D + FiberSection**: cierra dos huecos puntuales del catálogo 2D (geotecnia + plasticidad por flexión en frames). La más acotada de las tres abiertas; salda de paso el item #2 de deuda técnica.

### Opciones identificadas después (fuera de la bifurcación original A-E)

- **G. Elementos finitos isogeométricos (IGA).** Evaluada conceptualmente el 2026-08-25 a petición del usuario; **diferida sin fecha**, se registra para no repetir el análisis. **Veredicto: arquitecturalmente compatible, con una pieza faltante bien localizada** — no exige refactor transversal ni ruptura de contratos.
  - *A favor*: el contrato de `Element` es agnóstico a la naturaleza de las funciones de forma — sólo obliga a `compute_element_state(u_e) → (K_e, F_int_e)`. Nada en la clase base presupone interpolación Lagrangiana. `Node` es una estructura mínima (`id`, `coordinates`, `dofs`), así que sus coordenadas pueden ser las de un **punto de control** fuera del dominio físico; `get_global_dof_indices()` mapea `nodos × DOF_NAMES` genéricamente; el número de nodos por elemento lo fija cada subclase, de modo que un parche NURBS entra como entró el `Quad9`. Materiales, cuadratura, caché COO y solvers no lineales se reutilizan tal cual. La agnosticidad no es teórica: la Etapa 8 introdujo un campo físico nuevo sin tocar ensamblador, Dirichlet ni despachador algebraico.
  - *El obstáculo real está en las condiciones de frontera*, no en el ensamblaje. Las bases B-spline/NURBS **no son interpolatorias** (`u(x_A) ≠ u_A`), así que `node.fix_dof(name, value)` impondría el valor del **coeficiente de control**, no el del campo. Para empotramiento homogéneo coincide por partición de la unidad; para **asentamiento prescrito daría el valor equivocado en silencio** — justo el caso que el [ADR 0004](adr/0004-imposicion-de-condiciones-de-dirichlet.md) resolvió con esmero al descartar la penalización por "sólo aproximada". La vía canónica en IGA (proyección L² del dato de frontera, o Nitsche) no cabe en la firma `fix_dof(name, value)`.
  - *La abstracción correcta ya existe sin activar*: una BC IGA proyectada **es una restricción afín** entre coeficientes de control, es decir exactamente la forma `u_s = g_s + Σ α_si·u_mi` que `ConstraintSet` captura y cuyos MPC el ADR 0004 dejó contemplados como "extensión natural en una fase posterior".
  - *Trabajo restante, en orden de coste*: (1) activar los MPC de `ConstraintSet`, o aceptar sólo BC homogéneas en una primera versión; (2) familia de elementos con extracción de Bézier — contenida dentro de las subclases, es trabajo pero no obstáculo arquitectural; (3) **vía de entrada nueva** (knot vectors, grados, pesos, malla de control), el trozo más laborioso: el parser YAML actual lee `nodes` + conectividad y gmsh no sirve.
  - *Alcance*: **no es "un elemento más"** — es una etapa comparable a B o mayor, con ADR propio. Si lo que se busca es geometría exacta sobre superficies curvas, el camino barato ya existe: los cuadráticos representan curvatura isoparamétricamente, cuantificado en `tests/validation/test_lame_cylinder_3d.py` para `Hex20`/`Hex27`.

- **H. Apoyos elásticos y cimentaciones (Winkler).** Evaluada el 2026-08-25; **descartada como etapa propia**, se registra el argumentario. Hoy sólo existe Dirichlet exacto: no hay forma de expresar `F = −k·u`. Se distinguen tres piezas de coste muy distinto: (1) **resorte nodal concentrado** — una contribución `+k` a la diagonal de `K`, coste bajísimo, cubre apoyos discretos y aisladores; (2) **Winkler distribuido** `K_f = ∫ k_s·NᵀN dΓ` — misma estructura que una matriz de masa, coste medio, con validación analítica limpia disponible (Hetényi 1946, viga infinita sobre fundación elástica, longitud característica `β = ⁴√(k_s/4EI)`); (3) **unilateralidad** del terreno — problema de contacto no lineal, que sería el reflejo especular de `CableMaterial1D` y reutilizaría el patrón `ACCEPTS_UNILATERAL` ya probado; diferido hasta tener caso de uso real (regla de los dos casos). **Objeción de fondo**: Winkler es un modelo pobre — un resorte ignora lo que hace su vecino, así que no aporta física que el proyecto no tenga ya (mallar el terreno como continuo con `Quad4`/`Hex8` + Drucker-Prager es hoy posible). Su valor es de **economía de modelado**, no de capacidad física, y ese argumento pesa menos en un programa de investigación que en uno de despacho de cálculo. Si se retoma, encaja como complemento de la opción **E** (Winkler + Mohr-Coulomb + `FiberSection` forman un paquete geotécnico coherente), no como etapa independiente.

**Estado de la bifurcación**: cerradas **A**, **C** (en su núcleo mínimo) y **D**; quedan **B** y **E** como candidatas a Etapa 9, más la continuación de la línea térmica que la propia Etapa 8 habilita. Fuera de esa bifurcación original quedan registradas **G** (IGA, diferida sin fecha) y **H** (apoyos elásticos, descartada como etapa propia). La decisión es del usuario y se toma con la dirección que quiera dar al proyecto; ninguna está bloqueada técnicamente.

---

## Etapa 6 — Cierre del subsistema dinámico (ADR 0009 completo) · cerrada (2026-05-18)

**Capacidad entregada** (en una sesión, cinco commits): el subsistema modal/dinámico de la Etapa 4 quedó parcial al cierre original (fases 1, 3, 4); esta etapa completa el ADR 0009 con las cuatro fases pendientes (2, 5, 6, 7) y la variante HHT-α, además de aplicar las dos reglas arquitecturales C y D de la auditoría 2026-05-13 que esperaban evento.

**Lo entregado**:
- **HHT-α** (`HHTSolver` + `NewtonHHTSolver` como variantes de la familia Newmark) con disipación numérica controlada en altas frecuencias. Spec corta tipo extensión (Reglas §4).
- **Fase 2 ADR 0009 — mass lumping**: `compute_mass_matrix(lumping="lumped")` operativo en todos los elementos (truss/cable, frames 2D/3D, sólidos 2D Tri3/Quad4/Tri6/Quad8/Quad9). Helper centralizado `solidum/math/mass_lumping.py::lump_hrz` (HRZ canónico) + fórmula nodal directa para frames.
- **Regla C de auditoría aplicada** (dispatch declarativo en `entry.py::run_yaml` por atributo `PIPELINE_KIND`) + **`CentralDifferenceSolver`** (fase 5 ADR 0009) — leapfrog Belytschko-Liu-Moran con `M⁻¹` trivial, lineal y no lineal en una sola clase con parámetro `nonlinear`.
- **`HarmonicSolver`** (fase 6) — respuesta forzada armónica en el dominio de la frecuencia con aritmética compleja, barrido configurable, `HarmonicResult` con métodos `.amplitude()` y `.phase()`.
- **Regla D de auditoría aplicada** (`solidum/math/modal_response.py` agrupa `free_vibration` (movido) y los nuevos `participation_factors`, `response_spectrum_srss`, `response_spectrum_cqc`, `spectrum_from_sa`, `spectrum_tabulated`) + **`ResponseSpectrumSolver`** (fase 7) — combinación modal SRSS/CQC para análisis sísmico contra espectros normativos.

**Por qué se eligió como Etapa 6 en lugar de A/B/C/E**: ningún otro camino aprovechaba el alineamiento natural — `CentralDifference` requería la fase 2 (lumping), `HarmonicSolver` requería la regla C (3ª rama no clásica), `ResponseSpectrumSolver` requería la regla D (2º método sobre `ModalResult`). Las cinco piezas se encadenan en una secuencia obligada que cierra el ADR 0009 en su totalidad. Cualquier orden distinto hubiera dejado deudas internas más caras de pagar luego.

**Commits**: `73742f4` (HHT-α), `c92bf4e` (mass lumping), `e66473d` (regla C + central differences), `47a1333` (harmonic), `5151d7d` (regla D + response spectrum).

**Validación**: 558 tests verdes en suite global. 19 + 10 + 11 + 19 + 19 = ~80 tests nuevos del subsistema (mass_lumping, central_difference, harmonic, response_spectrum + reorganización menor en test_modal). Benchmarks: oscilador 1 GDL contra solución analítica, función de transferencia `H(ω) = 1/(K-ω²M+iωC)`, recuperación SRSS = CQC con modos separados, CQC ≠ SRSS con modos cercanos.

**Nada diferido de esta etapa**: el ADR 0009 queda cerrado en su totalidad. Las únicas extensiones futuras (excitación sísmica multi-directional simultánea CQC3, Δt adaptativo, MPC en frecuencia, generalized-α) requieren caso de uso real específico — no son piezas del ADR pendientes.

---

## Etapa 7 — Sólidos 3D acotados · cerrada (2026-05-19)

**Capacidad entregada**: primera línea de sólidos 3D en Solidum. Alcance acotado deliberadamente (Reglas.md §1 — regla de dos casos reales para centralizar) a los espejos naturales del Quad4/Tri3/Elastic2D:

- **`Hex8`** — hexaedro trilineal isoparamétrico (8 nodos, orden VTK_HEXAHEDRON). Default Gauss 2×2×2 (8 puntos), configurable a 3×3×3 o 1×1×1. Body load por cuadratura del elemento; face traction sobre 6 caras numeradas con normal saliente (ADR 0012) con cuadratura 2D Gauss 2×2.
- **`Tet4`** — tetraedro lineal (CST 3D, 4 nodos). 1 punto baricéntrico con peso 1/6. Caras numeradas por nodo opuesto. Masa consistente analítica `ρ·V·(1+δ_ij)/20`.
- **`Elastic3D`** — material elástico isótropo 3D. Sin variantes (en 3D no aplican plane_stress/plane_strain). Voigt 6D `[xx, yy, zz, xy, yz, xz]` con `γ_ij = 2·ε_ij` *engineering*.

**Infraestructura compartida**:

- **Cuadraturas 3D** registradas en `QuadratureRegistry`: `hex_1x1x1`, `hex_2x2x2`, `hex_3x3x3`, `tet_1`.
- **Mass lumping HRZ 3D**: `lump_hrz(M, total_mass=ρV, n_translational_dirs=3)` — el helper genérico ya soportaba 3D, solo había que invocarlo con la dimensión correcta. `_expand_scalar_mass_3d` en `solidum/elements/solid_3d/_shared.py` expande matriz escalar a bloque 3D vía Kronecker con I₃.
- **Convención Voigt 3D** fijada en `Reglas.md §5` y blindada por test (simetría C, positividad definida, cortantes puros en los tres planos).

**ADR que la articula**:
- [ADR 0012 — Sólidos 3D: Voigt 6D y cierre del contrato `internal_forces`](adr/0012-solidos-3d-y-voigt-6d.md). Tres decisiones simultáneas: (1) orden Voigt 3D `[xx, yy, zz, xy, yz, xz]` (extensión natural del 2D), (2) cierre de la deuda técnica #1 por dominio explícito (sólidos exponen `compute_gauss_state`, no `internal_forces`), (3) API de tracción de superficie con caras orientadas con normal saliente.

**Validación** (54 tests nuevos, suite 750 → 804):

- **Unitarios** (`tests/test_solid_3d.py`, 31 tests): Elastic3D (simetría/positividad de C, tracción uniaxial, los tres cortantes puros, compresión hidrostática con bulk modulus, rechazos de inputs); Hex8 y Tet4 (dimensiones/DOFs, simetría K_e, patch tracción uniaxial, jacobiano degenerado abortado, body load + face traction con balances de fuerza, masa consistente + lumped HRZ, gauss state).
- **Modos rígidos** (`tests/test_rigid_body_modes.py`, +10 tests): 6 modos rígidos en 3D (3 traslaciones + 3 rotaciones independientes en ejes x/y/z) para Hex8 y Tet4.
- **Cubo Lamé 3D** (`tests/validation/test_cube_lame_3d.py`, 3 tests): tracción uniaxial Hex8 (1 elemento) y Tet4 (cubo dividido en 5 tetraedros) — reproducen solución analítica exacta a precisión máquina porque el campo de desplazamientos es lineal en (x, y, z) y ambos elementos lo representan exactamente. Compresión hidrostática Hex8 — `σ = -p·I`, `ε = (-p/3K)·I` exactos.
- **MacNeal-Harder 3D** (`tests/validation/test_macneal_beam_3d.py`, 2 tests): cantilever esbelto `L/h = 30`. Hex8 12×1×1 documenta shear locking severo (ratio `u/u_EB` ∈ [0.10, 0.55]); refinamiento h reduce el locking monotonamente y la malla 48×4×4 alcanza ratio > 0.85.
- **Locking volumétrico 3D** (`tests/test_volumetric_locking_3d.py`, 2 tests): `ratio = u(ν=0.4999)/u(ν=0.3) < 0.6` con Hex8 4×1×1 documenta y blinda la limitación arquitectural; con ν=0.3 la flecha es razonable (`u/u_EB > 0.30`). Mitigaciones B-bar/F-bar diferidas (política idéntica al Quad4 2D).

**Cierre adicional**: deuda técnica #1 saldada por ADR 0012 — `internal_forces` aplica solo a estructurales 1D, sólidos exponen `compute_gauss_state`. Docstrings de `Element.internal_forces` y `ElementForces` actualizados con dominio acotado. STATUS.md tacha el item.

**Lo diferido y por qué**:

- ~~**Elementos 3D cuadráticos** (`Hex20`, `Hex27`, `Tet10`)~~ ✅ **Cerrado 2026-05-27 (sub-etapa A.ter)**: los tres elementos sobre base centralizada `_HigherOrderSolid3D` (regla de los dos casos reales aplicada al entrar el `Hex27`). Cuadraturas nuevas `tet_4` (Stroud orden 2 para K) y `tet_15` (Keast orden 5 para masa del Tet10). Suite 877 → 959 (+82 tests).
- ~~**Materiales 3D no lineales** (`VonMises3D`, `DruckerPrager3D`, `IsotropicDamage3D`)~~ ✅ **Cerrado 2026-05-21 (sub-etapa A.bis)**: tres materiales sobre Voigt 6D con tangente algorítmica consistente, cross-consistency vs 2D plane_strain a 10-14 decimales, integración Hex8/Tet4, y campaña de validación 3D consolidada (3 benchmarks publicados). Suite 804 → 877.
- ~~**NAFEMS LE10/LE2**~~ ✅ **Cerrado 2026-05-27 (sub-fase 4 de A.ter)**: **LE10** thick plate pressure entregado con cuadrante elíptico × espesor y σ_yy(D) ≈ −5.38 MPa contra el canónico para Hex20/Hex27 ([`test_nafems_le10.py`](../tests/validation/test_nafems_le10.py), 5 tests). **LE2 estricto descartado** por incompatibilidad de dominio: es un benchmark de *shell*, no de sólidos; se reemplazó conceptualmente por el **Lamé thick cylinder 3D** ([`test_lame_cylinder_3d.py`](../tests/validation/test_lame_cylinder_3d.py), 6 tests), extensión 3D del benchmark 2D del proyecto con solución analítica cerrada de Timoshenko-Goodier §28, que además es la primera demostración cuantitativa de la capacidad isoparamétrica curva. El `Tet10` sobre superficie curva queda como deuda técnica #8 (limitación de la malla, no del elemento).
- **NAFEMS LE3** (hemisferio) y los benchmarks de placas y cáscaras (FV1/FV5, FV2/FV12/FV32): **siguen abiertos**, pero no se cierran con más tests — requieren el componente placas/láminas de la opción B.

---

## Sub-etapa A.ter — Sólidos 3D cuadráticos · cerrada (2026-05-27)

**Capacidad entregada**: tres elementos sólidos 3D cuadráticos (Hex20, Hex27, Tet10) sobre base centralizada `_HigherOrderSolid3D`, paritarios con `_HigherOrderSolid2D` (Quad8/Quad9/Tri6). Cubre la totalidad de la matriz elemento 3D × material 3D (5 × 4 celdas, todas ✓).

- **`Hex20`** — serendípito 20 nodos (60 DOFs). Default Gauss `hex_3x3x3` (27 puntos), cara Quad8 con cuadratura `3x3`. Reduce drásticamente el shear locking del `Hex8`: malla 6×1×1 alcanza 97% u_EB en MacNeal beam (Hex8 12×1×1 < 55%). 6 modos de hourglass por elemento aislado con cuadratura reducida `hex_2x2x2`.
- **`Hex27`** — Lagrangiano triquadrático completo 27 nodos (81 DOFs). Default Gauss `hex_3x3x3`, cara Quad9. Patch test triquadrático `u_x = c·x²y²z²` exacto a precisión máquina (capacidad distintiva vs `Hex20`). 27 modos de hourglass por elemento aislado con `hex_2x2x2` (combinación más problemática del catálogo 3D — uso default obligatorio).
- **`Tet10`** — cuadrático tetraédrico 10 nodos (30 DOFs). Default Stroud `tet_4` (4 puntos, orden 2 — suficiente para K con J constante); masa con `tet_15` Keast (15 puntos, orden 5) fija independiente de K para garantizar correctitud en análisis modal/transitorio. Cara Tri6 con cuadratura `tri_3`. Cubo Lamé exacto en malla 5-tet con mid-edges compartidos.

**Centralización**: la base `_HigherOrderSolid3D` se introduce con la entrada del `Hex27` (sub-fase 2) — regla de los dos casos reales aplicada en el momento canónico, igual que `_HigherOrderSolid2D` cuando entró el `Quad9` después del `Quad8`. El `Hex20` original (sub-fase 1 standalone) se refactoriza a subclase de la base sin modificar tests. El `Tet10` (sub-fase 3) entra directamente sobre la base con face `Tri6` + cuadraturas tetraédricas.

**Cross-check elemento × material no lineal 3D** (sub-fase 5): 9 smoke tests (`tests/test_solid_3d_higher_order_nonlinear.py`) cierran las 9 celdas ○ → ✓ de la matriz 3D (`Hex20`/`Hex27`/`Tet10` × `VonMises3D`/`DruckerPrager3D`/`IsotropicDamage3D`). Patrón paritario con los smoke tests existentes para `Tet4` en `test_solid_3d_plasticity.py`.

**Sin ADR nuevo**: la decisión arquitectural (Voigt 6D, contrato `compute_gauss_state`, caras 3D con normal saliente) está tomada por ADR 0012. La centralización en `_HigherOrderSolid3D` se documenta en commits y specs (refactor zonal, no cambio arquitectural — Reglas.md §4).

**Validación** (82 tests nuevos, suite 877 → 959):

- **Unitarios** (`tests/test_solid_3d_higher_order.py`): patch lineal y cuadrático exactos para los tres elementos; patch triquadrático completo exclusivo para `Hex27`; conteos exactos de modos de hourglass con reducción (Hex20: 6, Hex27: 27); body load y face traction con balance global exacto; masa consistente y lumped HRZ con todas las entradas positivas.
- **Modos rígidos 3D** (`tests/test_rigid_body_modes.py`): 6 modos cero exactos para `Hex20`, `Hex27`, `Tet10`.
- **Cubo Lamé 3D** (`tests/validation/test_cube_lame_3d.py`): tracción uniaxial e hidrostática exactas a precisión máquina para los tres elementos. `Tet10` validado con malla 5-tet del cubo (mid-edges compartidos) — confirma compatibilidad multi-elemento.
- **MacNeal-Harder 3D** (`tests/validation/test_macneal_beam_3d.py`): `Hex20` 6×1×1 → ratio 0.97 (Hex8 12×1×1 → < 0.55); `Hex27` vs `Hex20` < 3% en flexión simple (el espacio extra triquadrático del Hex27 rara vez gobierna en flexión).
- **Locking volumétrico 3D** (`tests/test_volumetric_locking_3d.py`): atenuación cuantificada para `Hex20`/`Hex27` (ratio 0.80 vs < 0.6 Hex8), sin eliminación completa — B-bar/F-bar diferidos.
- **Cross-check no lineal** (`tests/test_solid_3d_higher_order_nonlinear.py`): 9 combinaciones (3 elementos × 3 materiales) pasan smoke tests end-to-end con `NonlinearSolver` y entran al régimen no lineal (α > 0 plasticidad, ω > 0 daño).

**Sub-fase 4 — validación externa cerrada 2026-05-27**:

- ~~**NAFEMS LE10 (placa elíptica gruesa) + LE2 (sólido cilíndrico)**: diferida~~ ✅ **Cerrada 2026-05-27** con un giro de alcance: LE2 estricto descartado por incompatibilidad (es benchmark shell, implementarlo con brick 3D sería trampear el setup), reemplazado por **Lamé thick cylinder 3D** — extensión natural del benchmark 2D ya validado en el proyecto, solución analítica cerrada de Timoshenko-Goodier §28 contra cada Gauss del modelo. Dos archivos nuevos en `tests/validation/`:
  - **`test_nafems_le10.py`** (5 tests verdes): NAFEMS LE10 thick plate pressure cuadrante elíptico × espesor, apoyo soft mid-plane, σ_yy(D) ≈ −5.38 MPa contra canónico para Hex20 6×6×2 y Hex27 5×5×2. Hallazgo de la sesión documentado en el módulo: el "thickness = 0.6 m" canónico es half-thickness, espesor total = 1.2 m.
  - **`test_lame_cylinder_3d.py`** (6 tests verdes Hex20+Hex27, 2 tests skip Tet10): error L²-relativo < 4% en σ_rr/σ_θθ/σ_zz y < 2% en u_r con malla 4×4×1; convergencia h monótona verificada. **Primera demostración cuantitativa de la capacidad isoparamétrica curva** de los elementos cuadráticos del proyecto (mid-edges automáticos sobre superficie cilíndrica por la no-linealidad del mapeo angular).
- **Tet10 sobre superficie curva** (limitación encontrada): la descomposición ingenua de cada celda hex en 5 tets rompe la representación isoparamétrica porque las aristas diagonales internas del hex se vuelven aristas del tet con mid-edges rectos. La cara del tet que toca la superficie curva tiene 2 de 3 mid-edges rectos → error σ_rr crece con refinamiento h (0.38 → 0.96 al refinar 2×2 → 8×8). Los tests Tet10 quedan en el archivo con `@pytest.mark.skip` y motivo arquitectural; código de la malla preservado para retoma futura. Nueva **deuda técnica #8** en STATUS — requiere mesher tetraédrico nativo (gmsh API) o descomposición específica del anillo cilíndrico. La capacidad del Tet10 sobre geometría plana ya está validada en `test_cube_lame_3d.py` (5-tet del cubo Lamé exacto a precisión máquina), así que no es fallo del elemento — es limitación de la malla.

Suite 959 → 970 (+11 tests verdes; +2 skip Tet10).

---

## Etapa 8 — Análisis térmico, núcleo mínimo C1 · cerrada (2026-08-25)

**Capacidad entregada**: conducción de calor de Fourier, estacionaria y transitoria, en 2D y 3D, con condiciones Dirichlet (temperatura impuesta, constante o variable en el tiempo) y Neumann (flujo prescrito). Es la opción **C** de la bifurcación, acotada por el usuario a su núcleo mínimo: **térmico puro, sin acoplamiento con el campo mecánico**.

**Alcance decidido por el usuario**, punto por punto:

1. Etapa 8 = análisis térmico (sobre B placas/láminas y E Mohr-Coulomb + FiberSection).
2. Sólo **C1 térmico puro** — sin deformación térmica `α·ΔT`, sin acoplamiento.
3. Sólo **conducción** — convección (Robin) y radiación explícitamente fuera.
4. **Conductividad tensorial** en el contrato; el constructor acepta un escalar y lo expande a `k·I`, de modo que la anisotropía no requerirá un material nuevo.
5. Ambas formas de capacidad, con **`lumped` por defecto** — invertido respecto a la dinámica estructural.
6. Dos elementos: **uno 2D y uno 3D** (`Quad4Thermal` + `Hex8Thermal`), no dos 2D.

### Componentes

- **`ThermalMaterial`** — familia con registro propio (`ThermalMaterialRegistry`), paralela a `CohesiveMaterial` (ADR 0010). Contrato `compute_flux(∇T) → (q, k)` en lugar de `compute_stress(ε)`. `FLUX_DIM` como **propiedad de instancia**, no `ClassVar`: la dimensión la fija el tensor `k` con que se construyó el material, no la clase.
- **`ThermalConduction`** — ley de Fourier `q = −k·∇T`. Valida simetría de `k` con tolerancia **escalada** a `1e-12·max|k|` y definición positiva por autovalores (no por determinante, que no la garantiza en 3D). `ρ` y `c` opcionales: sólo el transitorio los exige.
- **`Quad4Thermal`**, **`Hex8Thermal`** — un DOF escalar `T` por nodo, sobre base común `_ThermalSolid`. Comparten los kernels de forma y jacobiano de sus gemelos mecánicos. La base se introdujo con el **primer** elemento y no con el segundo, contra la regla habitual de los dos casos reales, porque la ecuación discretizada `C·Ṫ + K·T = F` es idéntica en 2D y 3D — la dimensión sólo cambia el tamaño de `B`. La decisión se validó a posteriori: el `Hex8Thermal` no necesitó tocar la base.
- **`ThetaMethodSolver`** — integración θ del sistema de **primer orden**. Familia propia con spec completa, **no variante de `NewmarkSolver`**: la familia Newmark integra la ecuación de segundo orden mediante hipótesis sobre la aceleración, y en conducción no existe segunda derivada temporal sobre la que aplicarlas. Resultado `ThermalTransientResult` propio.

### El hallazgo arquitectural

La infraestructura resultó **genuinamente agnóstica al campo**. El problema de pared plana se resolvió de extremo a extremo con el `Assembler` y el `LinearSolver` existentes —error 1.137e-13 contra el perfil lineal analítico— **sin tocar el ensamblador ni la imposición de Dirichlet**. El despachador algebraico reconoció por su cuenta la matriz térmica como candidata SPD a Cholesky.

Es la validación empírica de ADR 0003 (despacho algebraico por propiedades, no por tipo de análisis) y ADR 0004 (Dirichlet por eliminación sobre DOFs genéricos) frente a un campo físico que no existía cuando se tomaron esas decisiones. El régimen estacionario, de hecho, **no necesitó solver nuevo**.

### Dos decisiones físicas con el mismo criterio

Ambas invierten el default respecto al subsistema mecánico, y por la misma razón:

- **Capacidad `lumped` por defecto** (eje espacial) — la consistente produce oscilaciones espurias ante un frente térmico abrupto.
- **`θ = 1`, Euler implícito, por defecto** (eje temporal) — Crank-Nicolson es de orden 2 y más preciso, pero A-estable y **no L-estable**: su factor de amplificación tiende a `−1` para los modos altos en vez de a `0`.

El fenómeno común es la violación del **principio del máximo** de la ecuación de difusión: la solución exacta nunca excede los extremos de los datos iniciales y de frontera. Medido con `T_pared = 100` sobre un cuerpo a `0`: `θ = 1/2` alcanza `T_max = 151.27`; `θ = 1` se mantiene en `[0, 100]` y es monótono.

El criterio que decide no es el orden de convergencia sino **qué debe hacer el programa cuando el usuario no elige**: un resultado impreciso se detecta refinando el paso y se corrige; uno *cualitativamente imposible* desconcierta a quien no conozca la teoría de A- frente a L-estabilidad. Contrapartida asumida: el solver **reporta el orden efectivo** del esquema para que el coste en precisión del default robusto sea visible y no una penalización silenciosa.

### Validación

| Criterio | Resultado |
|---|---|
| Orden temporal θ=1 / θ=1/2 / θ=2/3 | **0.9944 / 2.0001 / 0.9933** |
| Convergencia al estacionario del `LinearSolver` | error **7.1e-14** |
| Sólido semi-infinito (Carslaw-Jaeger §2.4) | error **3.0e-3** |
| Cross-check 2D↔3D estacionario | error **9.9e-14** |
| Pared plana vs perfil lineal analítico | error **1.137e-13** |
| Estabilidad incondicional θ=1 vs divergencia θ=0 | confirmada con el mismo `Δt` |

El orden temporal se mide contra la solución **exacta del sistema semidiscreto** `exp(−C⁻¹K·t)·T₀`, no contra la del continuo: el error espacial actuaría como suelo y enmascararía la tasa. Es el test que realmente distingue los dos esquemas — un θ-method con signos cruzados aún converge al estacionario correcto, porque en el límite `t → ∞` el término temporal desaparece.

Las predicciones escritas en las specs se cumplieron sin ajuste: rango 7/8 del `Hex8Thermal` con un modo nulo (el de temperatura uniforme, análogo térmico del sólido rígido), 4 modos de hourglass con cuadratura reducida, y reparto `−q̄A/4` del flujo en las 6 caras.

**Sin ADR nuevo.** La etapa no rompió contratos ni introdujo un subsistema con decisiones arquitecturales propias: reutilizó el patrón de familia paralela ya establecido por ADR 0010 y la infraestructura de ADR 0003/0004 sin modificarla. `Reglas.md §5` recibió la convención de signo del flujo (`q̄ > 0` ⇔ saliente) y la constancia de que el problema térmico **no usa notación Voigt** — `∇T` es un vector genuino, no un tensor simétrico comprimido.

**Hallazgo colateral corregido**: la validación de `Assembler.assemble_mass_matrix` aconsejaba *"usa `0.0` si el material es sin masa por diseño"*, correcto en mecánica (ADR 0008) pero **físicamente incorrecto en térmico**, donde `density = 0` da capacidad calorífica nula y matriz singular. El ensamblador ahora delega el mensaje al material cuando éste sabe explicarse; ningún material mecánico cambia de comportamiento.

Suite 973 → 1176 (+203 tests, sin regresiones).

### Cierre pendiente de la etapa

No bloquea elegir la siguiente, pero queda anotado:

1. **Campaña de validación diferida**: cilindro hueco con perfil logarítmico (2D y 3D) y balance energético global. Ambos requieren mallar una **corona circular** — la misma carencia que difirió NAFEMS LE10 en su momento, y que sigue abierta.
2. ~~Sección térmica en los manuales~~ ✅ **cerrada 2026-08-25**, junto con el cableado de la vía YAML que la hacía documentable.

### Fuera de alcance, con rationale

Convección (Robin), radiación, conductividad `k(T)`, cambio de fase, paso de tiempo adaptativo y acoplamiento termomecánico. Las tres primeras y el cambio de fase exigen **Newton dentro de cada paso**; el `ThetaMethodSolver` actual asume el problema lineal y sin historia, que es justo lo que le permite factorizar `A = C + θΔt·K` una sola vez.

---

## Etapa 8+ — Horizonte largo

Lo que el proyecto **previsiblemente** abrirá tras la Etapa 8, sin orden cerrado:

- **Análisis termomecánico acoplado** — **desbloqueado** por la Etapa 8, que entregó el campo térmico sin acoplar. Es el C2 que quedó fuera del alcance de aquélla; su decisión pendiente es acoplamiento débil (secuencial) vs fuerte (monolítico), y toca el contrato `Material`.
- **Contacto mecánico** (penalización / Lagrangiano aumentado / mortar).
- **Grandes deformaciones** (lagrangiano total / actualizado en sólidos; corotacional 3D para frames; viscoplasticidad).
- **Materiales avanzados**: hiperelasticidad (Neo-Hooke, Mooney-Rivlin), plasticidad anisótropa, daño con regularización (gradient damage, phase-field), modelos para hormigón (Mazars, concrete damaged plasticity).
- **Optimización topológica y de forma** (SIMP, level-set).
- **Análisis estocástico / fiabilidad** (FORM, Monte Carlo sobre el modelo determinista).

Todos comparten un patrón: **incrementales sobre lo existente, no refactor estructural**. La arquitectura actual (contratos declarativos + registries + capa algebraica + caché de topología) está diseñada para soportarlos sin reabrir decisiones.

---

## Deuda técnica conocida

> **La autoridad sobre la deuda técnica es [STATUS.md](STATUS.md) §"Deuda técnica priorizada"**, donde vive la tabla numerada con su criterio de retoma. Esta sección solo da el resumen navegacional; si ambas discrepan, manda STATUS.

Al cierre de la sub-etapa A.ter quedan **6 items abiertos**, ninguno bloqueante (numeración de STATUS):

- **#2 — `FiberSection` para frames no-lineales**: los frames 2D/3D plastifican sólo en axial; la fluencia por flexión espera este componente. Entraría con la opción E.
- **#4 — Solver para softening severo con embedded discontinuity**: parcialmente cerrado por `DissipationArcLengthSolver` (funciona para daño continuo 1D/2D); el caso cohesivo+embedded con penalty `K_e` rígido sigue abierto y depende del #7.
- **#5 — Stress recovery superconvergente (SPR)**: NAFEMS LE1 da σ_yy(D) = 90.1 frente a 92.7 canónico, diferencia atribuible al recovery sin SPR.
- **#6 — B-bar / F-bar para Quad4/Quad8**: locking volumétrico documentado y blindado por test; espera caso de uso con material casi-incompresible.
- **#7 — LDLᵀ Bunch-Kaufman como backend algebraico**: habilitaría el *sign-of-pivot tracking* exacto que hoy se aproxima por el signo del determinante.
- **#8 — Validación del `Tet10` sobre superficie curva**: limitación de la malla (descomposición hex→5tets deja mid-edges rectos sobre la frontera curva), no del elemento.

**Cerrados** (se listan para que no se reabran por inercia):

- ~~**ADR 0002 incompleto para sólidos**~~ ✅ **Cerrado 2026-05-19 por el [ADR 0012](adr/0012-solidos-3d-y-voigt-6d.md)**: cierre por dominio explícito. `internal_forces` aplica sólo a elementos estructurales 1D; los sólidos 2D y 3D exponen `compute_gauss_state(U)` como API canónica. Es una decisión arquitectural tomada, no deuda pendiente.
- ~~**Reglas de disparo C y D**~~ ✅ **Ambas aplicadas 2026-05-18**: C con el despacho declarativo `PIPELINE_KIND`; D con `solidum/math/modal_response.py` agrupando el cómputo sobre modos. No quedan reglas de disparo arquitecturales activas.

---

## Documentos complementarios

El presente ROADMAP es uno de cuatro documentos navegacionales que escalan con el catálogo. Los cuatro existen y se mantienen sincronizados:

1. **[`ROADMAP.md`](ROADMAP.md) (este archivo)**: visión por etapas, decisiones de bifurcación. **Cambia cuando se cierra una etapa o se decide una bifurcación.**
2. **[`STATUS.md`](STATUS.md)**: foto de una página del estado actual — métricas (#tests, #elementos, #materiales, #solvers, #ADRs, #specs), capacidades habilitadas, deuda técnica con prioridades, próximo hito. **Se regenera o actualiza tras cada PR significativa.**
3. **[`MATRIZ.md`](MATRIZ.md)**: tabla cruzada elemento × material × solver con celdas validadas por test, en blanco (combinación posible no testeada) o vetadas (combinación prohibida con razón). Sirve para detectar huecos del catálogo y para usuarios externos que quieren saber qué combinación es segura.
4. **[`ONBOARDING.md`](ONBOARDING.md)**: documento de entrada para agente o humano que arranca sesión sin contexto previo. Indica qué leer primero (este ROADMAP, `Reglas.md`, último ADR), cómo correr tests, dónde está la memoria, qué patrones de cambio merecen ADR. **Punto de entrada recomendado para sesiones cold-start.**

---

## Cómo se cierra una etapa

1. Todos los componentes nuevos tienen spec en `docs/specs/` con `status: validated`.
2. Tests verdes (unitarios + acceptance + smoke YAML representativo).
3. Catálogos (`docs/catalogo_*.md`) y manuales (`manuals/sources/*/`) actualizados.
4. ADR(s) de la fase en estado `Aceptado`.
5. Entrada en este ROADMAP movida a "cerrada" con resumen ejecutivo.
6. `STATUS.md` regenerado con la nueva métrica.

---

*Última actualización: 2026-08-25 — **Etapa 8 cerrada: análisis térmico (núcleo mínimo C1)**. Nueva sección de etapa con el alcance decidido por el usuario punto por punto, los cuatro componentes entregados, el hallazgo arquitectural (la infraestructura resultó agnóstica al campo — el estacionario no necesitó solver nuevo) y la tabla de validación. Actualizado el estado de cabecera (1148 tests), cerrada la **opción C** de la bifurcación en su núcleo mínimo dejando constancia de que la termomecánica acoplada de `Reglas.md §0` sigue sin existir, y desbloqueado el acoplamiento en el horizonte largo (deja de estar condicionado a "si C entró antes"). El próximo hito pasa a ser la elección de la Etapa 9 entre B, E y la continuación de la línea térmica.*

*Anterior 2026-08-25 — **Sincronización con STATUS** (sin cambios de código). Este documento había quedado desfasado en cuatro puntos y, por ser lectura de arranque en sesiones cold-start, dirigía a trabajo ya hecho: (1) la sección "Deuda técnica conocida" listaba como abiertos el **ADR 0002** (cerrado el 19-mayo por el ADR 0012) y las **reglas de disparo C y D** (aplicadas el 18-mayo), y omitía los 6 items realmente abiertos — ahora delega la autoridad en STATUS y resume los vigentes; (2) los enlaces de esa sección apuntaban a `~/.claude/projects/.../memory/`, ruta vacía que `CLAUDE.md` prohíbe explícitamente (la memoria vive en `.claude/memory/` del repo) — el mismo error estaba en ONBOARDING §4 y se corrigió allí también; (3) la **opción A** de la bifurcación seguía listada como futura pese a haberse ejecutado en la Etapa 7 + A.bis + A.ter, y la frase de cierre condicionaba la decisión al cierre de una Etapa 5 cerrada en mayo; (4) **NAFEMS LE10/LE2** figuraban como diferidos cuando LE10 cerró en la sub-fase 4 de A.ter y LE2 fue descartado con razón documentada (es benchmark de shell), sustituido por el Lamé thick cylinder 3D. Actualizado también el estado de cabecera (973 tests) y añadido el próximo hito explícito. Colateralmente, MATRIZ §2 omitía `DissipationArcLengthSolver` de la tabla de solvers pese a ser seleccionable en YAML: añadido, y su cobertura verificada contra los cuatro registros (46/46 componentes).*

*Anterior 2026-05-27 — **Sub-etapa A.ter cerrada por completo**: sólidos 3D cuadráticos (`Hex20`, `Hex27`, `Tet10`) sobre base centralizada `_HigherOrderSolid3D`, más la sub-fase 4 de validación externa (NAFEMS LE10 + Lamé thick cylinder 3D). Suite 877 → 970. La matriz 3D queda con paridad funcional completa frente al 2D: 5 elementos × 4 materiales, todas las celdas en ✓.*

*Anterior 2026-05-21 — **Sub-etapa A.bis cerrada: materiales 3D no lineales**. Los tres materiales (`VonMises3D`, `DruckerPrager3D`, `IsotropicDamage3D`) cierran la paridad funcional 3D↔2D para plasticidad y daño. Suite 804 → 877 (+73 tests: 56 directos + 11 validación 3D consolidada + 6 colaterales).*

*Anterior 2026-05-19 — **Etapa 7 cerrada: sólidos 3D acotados (ADR 0012)**. Alcance: `Hex8` (trilineal), `Tet4` (CST 3D), `Elastic3D` (isótropo). Convención Voigt 6D `[xx, yy, zz, xy, yz, xz]` fijada en `Reglas.md §5`. Suite 750 → 804 (+54 tests verdes a la primera). Cierre adicional: deuda técnica #1 saldada por ADR 0012 (cierre por dominio explícito — sólidos exponen `compute_gauss_state`, no `internal_forces`).*

*Anterior 2026-05-19: sin nuevas etapas. Sesión de saneamiento post-auditoría 2026-05-18 cerró 23 hallazgos (3 críticos + 4 altos + 16 medios/bajos) sin tocar el ROADMAP estructural. Suite 558 → 585.*
