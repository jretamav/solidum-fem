# Decisiones de arquitectura (ADR)

Los Architecture Decision Records (ADR) son el registro persistente de las decisiones que han modelado la arquitectura de Solidum FEM. Cada ADR responde, en una página, a tres cuestiones: contexto, decisión y consecuencias. En este capítulo se resumen los ADR vigentes; el texto completo reside en `docs/adr/`.

## ADR 0001 — Mecanismos de transparencia conceptual

**Fecha**: 19 de abril de 2026. **Estado**: aceptado.

**Contexto.** El documento `Reglas.md` §2 establece el modelo "calculadora con manual disponible": el usuario delega la escritura del código pero conserva su comprensión conceptual. Para que dicho principio se opere, y no permanezca como mera declaración de intenciones, el proyecto requiere mecanismos prácticos que prevengan la deriva hacia la opacidad a medida que el código crece.

**Decisión.** Se introducen tres mecanismos vivos, de mantenimiento económico y mutuamente reforzantes:

1. *Mapa de arquitectura.* El archivo `docs/arquitectura.md` consiste en una página con diagrama de bloques y flujo de datos. El asistente lo actualiza al modificarse la topología del sistema.
2. *Conceptos clave.* El archivo `docs/conceptos_clave.md` contiene un inventario explícito de los conceptos sostenedores del sistema, cada uno con dos o tres frases conceptuales. Define el "test de no opacidad": el usuario debe poder explicar cualquier entrada en términos generales en cualquier momento.
3. *Modo "explícame X" bajo demanda.* Protocolo conversacional, no artefacto en repositorio. Constituye la válvula de escape ante opacidad puntual sin necesidad de sesiones programadas.

**Consecuencia para este manual.** El presente `Architecture_manual.pdf` es el siguiente eslabón en la cadena de transparencia conceptual: amplía el mapa de arquitectura a una lectura cómoda en formato de libro, manteniendo la disciplina de no entrar al detalle de cada componente.

## ADR 0002 — Interfaz pública de resultados para consumidores externos

**Fecha**: 22 de abril de 2026. **Estado**: aceptado.

**Contexto.** Solidum FEM se consume desde una GUI externa (SolidumBAR, llamado FenixBAR hasta el 2026-09-30, programa de análisis estructural con elementos barra). Tras una solución, el consumidor requiere información suficiente para visualizar el post-proceso estándar de Computer-Aided Engineering (CAE): deformada escalada, reacciones en apoyos y diagramas de fuerzas internas sobre cada barra. Esta información no se exponía de forma uniforme: el componente `VtkExporter` cubría únicamente algunos elementos y las fuerzas internas carecían de método público homogéneo.

**Decisión.** Interfaz pública estructurada en tres niveles:

1. *Contrato por elemento.* Todos los elementos tipo armadura, cable o marco implementan `internal_forces(U) -> ElementForces`, que devuelve las fuerzas internas en ejes locales con claves homogéneas por familia (`{N}` para armadura y cable; `{N, V, M}` para marco 2D; `{N, Vy, Vz, T, My, Mz}` para marco 3D). La convención de signos es la del proyecto (capítulo anterior).
2. *Agregado en `Domain`.* El objeto `SolveResult` se calcula de forma anticipada al final de `solver.solve()` y es inmutable. Expone los desplazamientos globales, las cargas aplicadas, las reacciones y un diccionario de fuerzas internas por elemento.
3. *Puntos de entrada oficiales.* Las funciones `solidum.run(case)` y `solidum.run_yaml(path)` devuelven directamente un `SolveResult`. Cualquier consumidor externo opera contra esta interfaz.

**Consecuencia.** El consumidor externo no recalcula las fuerzas internas a partir de `U`; la lógica de MEF reside íntegramente en Solidum FEM. La incorporación de un elemento nuevo obliga a implementar `internal_forces` con las claves de su familia: el contrato es explícito.

## ADR 0003 — Despachador de algoritmo algebraico

**Fecha**: 26 de abril de 2026. **Estado**: aceptado; fases 1 y 2 implementadas.

**Contexto.** Dentro de cada iteración del método de Newton-Raphson se requiere la resolución del sistema lineal `K · δU = R`. El algoritmo utilizado hasta el momento era SuperLU (factorización LU general), válido para todo caso pero subóptimo cuando `K` es simétrica y definida positiva, situación habitual en estática lineal y en numerosas no linealidades suaves. La factorización de Cholesky es típicamente entre 1.5 y 2 veces más rápida y consume considerablemente menos memoria en dichos casos.

**Decisión.** Se introduce un subsistema algebraico independiente del solver no lineal, con varios algoritmos disponibles y un despachador que selecciona el adecuado de forma automática:

- `LUSolver` (SuperLU) — algoritmo universal.
- `CholeskySolver` (CHOLMOD) — para matrices simétricas y definidas positivas; opcional según disponibilidad de `scikit-sparse`.
- `LDLTSolver` — reservado para la fase 2 (matrices simétricas indefinidas).

El despachador `select_solver(props, override)` recibe un objeto `StiffnessProperties` (atributos `is_symmetric`, `is_positive_definite`) y selecciona la factorización de Cholesky cuando aplica, o la factorización LU en caso contrario. El usuario puede forzar un algoritmo desde el archivo YAML mediante el campo `linear_algebra` en calidad de herramienta de diagnóstico, no como decisión de modelado.

**Consecuencia.** Los solvers no lineales (`LinearSolver`, `NonlinearSolver`, `ArcLengthSolver`) ya no invocan SuperLU directamente: solicitan un solver algebraico al despachador y se desentienden del algoritmo numérico subyacente. La incorporación de un nuevo algoritmo (Pardiso, MUMPS, métodos algebraicos multimalla iterativos) afecta exclusivamente a `solidum/math/linalg/`. La separación entre estrategia (solver no lineal) y táctica (subsistema algebraico) queda explícita en el código.

## ADR 0004 — Imposición de condiciones de Dirichlet por eliminación directa

**Fecha**: 30 de abril de 2026. **Estado**: aceptado; fases 1 (Dirichlet nodal) y 2 (restricciones lineales afines) implementadas. Queda abierta la interpretación de las reacciones con restricciones multipunto (punto 10 de la fase 2).

**Contexto.** Las condiciones de contorno de Dirichlet, tanto los empotramientos nodales clásicos como las restricciones multipunto (MPC: asentamientos prescritos, apoyos en plano oblicuo, periodicidad de celda unitaria, unión rígida master-slave), se imponían previamente por método de penalidad: añadir términos diagonales muy grandes al sistema. Dicho método introduce condicionamiento numérico artificialmente malo, contamina la simetría y la positividad de la matriz de rigidez y no es exacto a redondeo.

**Decisión.** Se adopta la imposición por eliminación directa. Toda restricción se expresa en forma afín `u_s = g_s + Σ α_si · u_mi` y se acumula en un `ConstraintSet` del subpaquete `solidum.bc`. El ensamblador construye un operador disperso `T` y un vector `g` tales que `u = T · u_libre + g`. El sistema entregado al subsistema algebraico es `K_red = TᵀKT`, `F_red = Tᵀ(F − K · g)`. Las restricciones encadenadas (un esclavo cuyo maestro es a su vez esclavo) se resuelven por cierre transitivo al construir `T`; los ciclos y redeclaraciones inconsistentes se detectan en validación temprana.

**Consecuencia.** La imposición es exacta a redondeo, preserva simetría y positividad definida de `K`, y deja el sistema reducido apto para los algoritmos del subsistema algebraico, incluidos solvers iterativos como gradiente conjugado o GMRES. La misma maquinaria cubre empotramiento, asentamiento, apoyo oblicuo, periodicidad y unión rígida sin ramificación específica por tipo.

## ADR 0005 — Logging configurable en lugar de `print` directo

**Fecha**: 7 de mayo de 2026. **Estado**: aceptado.

**Contexto.** La trazabilidad del solver durante un análisis (incrementos de carga, iteraciones de Newton, ratios de convergencia, factorizaciones reusadas) se emitía mediante `print` directo en la salida estándar. Esto impedía silenciar el ruido en pipelines de batch testing, redirigir la traza a un archivo de log paralelo y diferenciar niveles de gravedad (información ordinaria, advertencias, errores recuperables).

**Decisión.** Se introduce un *logger* configurable basado en el módulo estándar `logging` de Python, expuesto al resto del proyecto vía `solidum.logging.get_logger(name)`. Todos los componentes consumen este *logger* en lugar de `print`. El nivel y el formato se configuran de forma centralizada y pueden sobreescribirse desde el archivo YAML del caso o programáticamente.

**Consecuencia.** Tests automatizados pueden silenciar la traza sin perder la información de errores. La integración con herramientas externas (GUI SolidumBAR, batería de tests, pipelines CI) se hace mediante la API estándar de `logging`, sin acoplamiento a `sys.stdout`. Las advertencias internas relevantes (degradación de Cholesky a LU, bisección de incremento) quedan marcadas con su nivel adecuado.

## ADR 0006 — Tolerancias del criterio de admisibilidad constitutiva

**Fecha**: 12 de mayo de 2026. **Estado**: aceptado.

**Contexto.** Cada modelo constitutivo no lineal tiene un criterio de admisibilidad — una función `f(estado)` que decide si el material está en régimen reversible (`f ≤ 0`) o disipativo (`f > 0`). En aritmética de doble precisión la comparación no puede ser exacta: tras un return mapping, el ruido de redondeo deja `f` cerca de cero pero no exactamente cero, por lo que la comparación se escribe `f ≤ tol`. La implementación anterior usaba una tolerancia absoluta hardcodeada (`PLASTIC_YIELD_TOL = 1e-9`), frágil bajo cambio de unidades y no adaptativa al estado en modelos con endurecimiento severo o sin parámetro de referencia inicial.

**Decisión.** Se adopta el esquema combinado absoluto + relativo, estilo solver de ecuaciones diferenciales ordinarias: `f ≤ ADMISSIBILITY_TOL_ABS + ADMISSIBILITY_TOL_REL · escala(estado)`. Cada material declara su escala característica vía el método `admissibility_scale(state_vars)`, devolviendo un esfuerzo positivo característico del criterio en el estado corriente: `σ_y + H·α` para plasticidad J2 1D, `E·κ_0` para daño isótropo, `√(2/3)·(σ_y + H·α)` para Von Mises 2D. La comparación se encapsula en `Material.is_admissible(f, state_vars)` y en `Material.admissibility_tol(state_vars)` para que los modelos no la repliquen.

**Consecuencia.** Misma física en MPa, Pa, GPa o adimensional produce idéntico estado interno hasta precisión de doble. En endurecimiento severo la tolerancia crece junto con la fluencia. Modelos energéticos futuros (phase-field, gradient damage) encajan declarando su escala (`√(2·E·g_c)`) sin tocar la política. Los modelos con kernel compilado (Numba) consumen `admissibility_tol` desde código Python antes de invocar el kernel, preservando la centralización de la fórmula.

## ADR 0007 — Tolerancias de convergencia en solvers no lineales

**Fecha**: 13 de mayo de 2026. **Estado**: aceptado.

**Contexto.** El criterio de convergencia del `NonlinearSolver` y del `ArcLengthSolver` era puramente relativo: `‖residuo‖ / ref_force ≤ tol` y `‖δU‖ / ‖U‖ ≤ tol`, con un único parámetro `tol` para los dos criterios y sin término absoluto real. Esta forma arrastraba tres limitaciones estructurales: no había tolerancia absoluta cuando la escala relativa colapsaba transitoriamente, una sola tolerancia mezclaba dos criterios físicamente distintos (incrementabilidad del iterado vs equilibrio), y la fórmula se replicaba en los dos solvers a mano.

**Decisión.** Se aplica el mismo patrón del ADR 0006 a la convergencia, separado por criterio. La política reside en una clase pequeña `ConvergenceCriterion` en `solidum/math/convergence.py`, que encapsula configuración (cuatro tolerancias: `rtol_force`, `rtol_disp`, `atol_force_factor`, `atol_disp_factor`) y estado calibrado (los `atol` efectivos derivados de las escalas del primer ensamblaje). Los solvers no lineales reciben una instancia en el constructor (parámetro `convergence`), calibran al inicio de `solve()` y delegan en `evaluate()` la comparación. La semántica es AND, no OR (ambos criterios deben cumplirse simultáneamente), coherente con la convención de Bathe, Crisfield y Owen-Hinton.

**Consecuencia.** Invariancia bit-paritaria bajo cambio de unidades (el mismo problema en N/m, kN/mm o MN/m converge en el mismo número de iteraciones). Régimen degenerado (carga inicial nula, control en desplazamiento puro) deja de oscilar espuriamente: el piso absoluto `atol_force` autoderivado define cuándo se da por convergido un residuo numéricamente cero. Problemas casi-rígidos y plásticos perfectos se afinan independientemente vía `rtol_disp` vs `rtol_force` sin compromisos. La política vive en un único punto del código; un solver no lineal futuro la consume sin replicarla.

## ADR 0008 — Densidad como propiedad del material

**Fecha**: 13 de mayo de 2026. **Estado**: aceptado.

**Contexto.** El cableado de fuerza de cuerpo introducido previamente aceptaba un vector `body_force: [bx, by, bz]` global, aplicable uniformemente a todos los elementos. Esto servía para estructuras monomaterial — el usuario calculaba `ρ · g` una vez y lo declaraba — pero fallaba en estructuras multimaterial: con un único vector global, no había forma de distinguir el peso propio del acero del hormigón o de un cable tensado. Adicionalmente, la densidad es la magnitud física que abre la puerta a la matriz de masa elemental y por tanto a todo el subsistema dinámico futuro.

**Decisión.** La densidad es propiedad del material, no del elemento. Se introduce como atributo `density` en la clase base `Material` con valor por defecto `None` (no declarado). La densidad es físicamente intrínseca al medio continuo; vive en el mismo nivel jerárquico que `E`, `ν` o `σ_y`. En el YAML se acepta el bloque atajo `gravity: [0, -9.81]` complementario a `body_force` (mutuamente excluyentes). El ensamblador incorpora `assemble_self_weight(g)` que itera sobre elementos y aplica `b_element = material.density · g`. Si algún material tiene `density=None`, el método falla con `ValueError` listando los materiales afectados — enforcement diferida al uso, sin posibilidad de propagar masa cero silenciosa a un análisis dinámico.

**Consecuencia.** Peso propio físicamente correcto en estructuras multimaterial. API del YAML más intuitiva. Tests existentes con `body_force` siguen pasando. La fase de análisis modal y dinámico hereda la propiedad sin rediseño: `Assembler.assemble_mass_matrix` reutiliza el mismo patrón de validación.

## ADR 0009 — Análisis modal y dinámico

**Fecha**: 13 de mayo de 2026. **Estado**: aceptado; **todas las fases (1-7) implementadas** al 18 de mayo de 2026, junto con la variante HHT-α y las reglas C y D de refactor arquitectural aplicadas.

**Contexto.** El ADR 0003 había identificado análisis modal y dinámica implícita entre los análisis futuros que la capa algebraica iba a habilitar; el ADR 0008 había puesto la densidad como propiedad del material disponible. Quedaba abrir el subsistema con dos análisis distintos pero relacionados: el análisis modal — problema algebraico de valores y vectores característicos generalizado `K · φ = ω² · M · φ` —, cuyas soluciones admiten interpretación dinámica pero no constituyen análisis dinámico per se; y el análisis dinámico propiamente dicho, con dependencia temporal explícita — integración de `M · ü + C · u̇ + K · u = F(t)`.

**Decisión.** Se acometen ambos análisis con una hoja de ruta común y siete decisiones arquitecturales fijadas hoy para no reabrir debates en fases posteriores:

1. Matriz de masa **consistente** en la fase 1 (con contrato extensible a `lumped` mediante un parámetro `lumping` ya presente en la firma de `compute_mass_matrix`, suficiente para la futura integración explícita).
2. Ensamblaje de M paralelo al de K, reutilizando la topología COO cacheada del `Assembler` y la enforcement de densidad del ADR 0008.
3. `ModalSolver` registrado en el mismo `SolverRegistry` que los solvers estáticos: no se crea un `IntegratorRegistry` paralelo, manteniendo el patrón "un solver = un análisis = una clase registrada".
4. Dirichlet en el problema modal por eliminación directa (mismo mecanismo del ADR 0004); el término independiente `g` de restricciones lineales no aplica en modal.
5. Estado dinámico transitorio `(u, u̇, ü)` fuera de `Node`: vive en un `TransientResult` específico, preservando la pureza topológica/geométrica de la clase nodal y evitando que análisis estáticos paguen el coste de campos que no usan.
6. Amortiguamiento Rayleigh `C = α · M + β · K` como entrada estándar para la fase de dinámica transitoria, con coeficientes calibrables modalmente a partir de dos pares `(ξ, ω)`. Amortiguamiento modal por modo y otros esquemas (Caughey, local) diferidos.
7. Convención de unidades heredada del modelo (ADR 0008): el usuario es responsable de la consistencia; las tolerancias adimensionales del ADR 0007 garantizan invariancia bajo cambio de unidades.

Las siete fases quedan implementadas y validadas con tests contra solución analítica y recuperación a paridad de bits del caso lineal en ausencia de plasticidad:

- **Fase 1** (`ModalSolver` + ARPACK Lanczos con shift-invert). Cerrada 2026-05-13.
- **Fase 2** (mass lumping HRZ canónico en sólidos + nodal directo en frames). Helper centralizado `solidum.math.mass_lumping.lump_hrz`. Cerrada 2026-05-18 (commit `c92bf4e`).
- **Fase 3** (`NewmarkSolver` con `(β, γ)` parametrizables, factorización única, Rayleigh, cargas por callback). Cerrada 2026-05-13.
- **Fase 3-bis** (`HHTSolver` como variante de Newmark con disipación numérica controlada por `α ∈ [−1/3, 0]`). Cerrada 2026-05-18 (commit `73742f4`).
- **Fase 4** (`NewtonNewmarkSolver` con Newton interno, jacobiano `J = M + γΔt·C + βΔt²·K_t`, convergencia dual ADR 0007, Rayleigh constante con `K_0`; `NewtonHHTSolver` análogo). Cerrada 2026-05-13/-18.
- **Fase 5** (`CentralDifferenceSolver` — leapfrog explícito Belytschko-Liu-Moran, lineal y no lineal en una clase con `nonlinear`, estabilidad condicional CFL con detección a posteriori). Cerrada 2026-05-18 (commit `e66473d`).
- **Fase 6** (`HarmonicSolver` — respuesta forzada armónica `(−ω²M+iωC+K)·û=F̂` con barrido configurable, aritmética compleja). Cerrada 2026-05-18 (commit `47a1333`).
- **Fase 7** (`ResponseSpectrumSolver` — combinación modal SRSS/CQC sobre espectro de respuesta). Cerrada 2026-05-18 (commit `5151d7d`).

**Reglas arquitecturales aplicadas durante el cierre del ADR**:

- **Regla C** (auditoría 2026-05-13). Disparada por `CentralDifferenceSolver` como primer solver fuera de la jerarquía de Newmark. `entry.py::run_yaml` despacha por atributo de clase `PIPELINE_KIND ∈ {"static", "modal", "transient", "harmonic", "spectrum"}` en lugar de cadena de `isinstance`. Solvers no clásicos futuros no requieren tocar `entry.py`.
- **Regla D** (auditoría 2026-05-13). Disparada por `ResponseSpectrumSolver` como segundo método de cómputo sobre `ModalResult`. El algoritmo de `free_vibration` migra a `solidum/math/modal_response.py` junto con los nuevos `participation_factors`, `response_spectrum_srss`, `response_spectrum_cqc`, `spectrum_from_sa`, `spectrum_tabulated`. `ModalResult.free_vibration` queda como wrapper delgado, preservando la API histórica.

**Consecuencia.** El subsistema modal/dinámico/espectral queda **completo en su totalidad**. Las extensiones futuras (excitación sísmica multi-direccional simultánea CQC3, multi-support seismic, generalized-α) requieren caso de uso específico — no son piezas del ADR pendientes. La clase `Node` conserva su semántica original.

**Ampliaciones (9 de noviembre de 2026).** (1) §8, decisión del usuario: la razón de masa efectiva acumulada se mide contra la masa que puede participar, $M_\text{part} = \mathbf b^\top\mathbf M_{LL}^{-1}\mathbf b$, la suma de las masas efectivas de todos los modos, de modo que con todos los modos vale 1 exactamente (con masas de nudo, grados de libertad sin masa y restricciones multipunto); `ModalResult.participating_mass` y `cumulative_effective_mass_ratio`. (2) Paso adaptativo en `NewtonNewmarkSolver` y `NewtonHHTSolver`, decisión del usuario: el paso cuyo Newton no converge se repite desde el último estado convergido en dos mitades, hasta `max_subdivisions` niveles (6 por omisión, Δt/64); el historial queda en la malla de tiempos del usuario. (3) Los transitorios implícitos admiten grados de libertad sin masa (ADR 0024, §6 bis).

## ADR 0010 — Discontinuidades interiores embebidas

**Fecha**: 13 de mayo de 2026. **Estado**: aceptado; fases 1, 2, 3 y 3b implementadas; fase 4 (curva de Van Vliet) validada con grieta prescrita el 29 de septiembre de 2026; fase E (seguimiento de grieta) completada el 10 de noviembre de 2026; las demás, diferidas con justificación.

**Contexto.** Abrir la línea de fractura computacional. Entre las alternativas (fisura difusa, elementos cohesivos de interfaz, XFEM, discontinuidades embebidas) se elige la que el autor domina a nivel doctoral: discontinuidades interiores embebidas en aproximación discreta (Retama, 2010).

**Decisión.** Cinco abstracciones transversales: (1) una familia `CohesiveMaterial` **paralela** a `Material`, no subclase, porque relaciona tracción con salto de desplazamiento y no esfuerzo con deformación; (2) elementos con **grados de libertad enriquecidos a nivel elemental**, invisibles para el ensamblador; (3) **condensación estática local** del salto dentro del elemento, de modo que el sistema global no cambia de tamaño; (4) un `DiscontinuityState` que guarda la orientación y la historia de la discontinuidad; (5) una semántica de **elemento que cambia de modo** (continuo a agrietado), cuya activación se decide entre pasos con el hook `prepare_step` y no dentro del Newton, para evitar oscilaciones sobre el umbral. Se fija además la longitud de discontinuidad `l_d = (A/h)·cos(θ−α)`, aportación del capítulo 6 de la tesis, validada numéricamente por primera vez.

**Consecuencia.** Fractura computacional sin tocar el ensamblador, el despachador algebraico, el Newton ni el parser. Las cuatro abstracciones son reutilizables para modelos cohesivos clásicos, modos incompatibles o EAS. El límite conocido es el solver para ablandamiento severo con penalización cohesiva rígida (deuda técnica). *(Nota 2026-09-29: ese límite se diagnosticó antes de corregir la tangente cohesiva el 2026-09-23; con la tangente corregida, `NonlinearSolver` con control por desplazamiento recorre el ablandamiento hasta la separación, y el fallo de `ArcLengthSolver` no se ha vuelto a comprobar.)*

**Ubicación.** Desde el ADR 0020 la formulación vive fuera del programa principal, en el módulo de usuario `discontinuities` (`solidum/user/discontinuities/`), y un modelo la carga con `user_modules: [discontinuities]`. El ADR se conserva como registro de la formulación.

**Actualizaciones (10 de noviembre de 2026).** *Seguimiento de grieta* (fase E, decisiones del usuario): el nodo solitario se elige por la trayectoria de la grieta y no elemento a elemento, con un seguimiento global por región conexa (`GlobalCrackTracking`) o una trayectoria prescrita (`PrescribedCrackPath`), y `l_d` lleva signo. *Familia cohesiva*: pasó al programa principal como familia de materiales de interfaz (ADR 0033); lo que aquí se dice de `CohesiveMaterial` vale para `InterfaceMaterial`. El elemento acepta también la ley de modo mixto `CohesiveDamageMixedMode`, que transmite cortante por la grieta; la cinemática sigue siendo KOS. El módulo conserva el elemento, el estado de la discontinuidad y el seguimiento.

## ADR 0011 — Robustez del Newton: line search y diagnóstico de divergencia

**Fecha**: 18 de mayo de 2026. **Estado**: aceptado con una enmienda.

**Contexto.** Los solvers no lineales divergían con un `RuntimeError` genérico que no distinguía entre causas cualitativamente distintas, y no tenían globalización.

**Decisión.** Dos mejoras coordinadas: un **line search** de backtracking con condición de descenso sobre el residuo, y una **familia de excepciones tipadas** que clasifica la divergencia (oscilación del Newton, tangente singular, carga por encima de la capacidad, modo desconocido) con métricas y una pista textual. **Enmienda**: el line search, previsto activo por defecto, quedó **desactivado por defecto** al comprobarse que rompía casos de daño que convergían sin él; con Newton completo y tangente consistente el backtracking apenas aporta (de Borst y Sluys, tabla 3.1).

**Consecuencia.** Un fallo de convergencia se reporta con su causa probable y la acción típica, no con un mensaje genérico.

## ADR 0012 — Sólidos 3D y notación de Voigt en seis componentes

**Fecha**: 19 de mayo de 2026. **Estado**: aceptado.

**Contexto.** Abrir la familia de elementos sólidos 3D obligaba a fijar la notación de Voigt tridimensional y a decidir qué hace el contrato de fuerzas internas en un sólido.

**Decisión.** (1) Voigt 3D con orden `[xx, yy, zz, xy, yz, xz]` y deformaciones angulares *engineering*, extensión natural del orden 2D (difiere de ABAQUS en la permutación `yz ↔ xz`); (2) el contrato `internal_forces()` **no se aplica a sólidos**: un continuo no tiene resultantes seccionales, su salida canónica es el campo de esfuerzos por punto de Gauss; (3) una API de tracciones y caras 3D paralela a la de bordes 2D.

**Consecuencia.** Especificaciones de sólidos 3D sin ambigüedad de notación, y el contrato de resultados del ADR 0002 cerrado por dominio explícito en vez de quedar como deuda.

## ADR 0013 — Orientación material y ortotropía

**Fecha**: 9 de septiembre de 2026. **Estado**: aceptado.

**Contexto.** El primer material anisótropo (`Orthotropic2D`, para la investigación con bambú) plantea dónde vive el ángulo de fibra: en el material, en el elemento, o en un campo sobre la malla.

**Decisión.** El ángulo va **en el material**, como atajo consciente y desechable: es una entrega puramente aditiva que no toca los contratos `Material` ni `Element`, y resuelve el caso real (probetas a ángulo constante). Se declara explícitamente que la orientación es conceptualmente una propiedad del elemento, y se escribe la condición de migración: cuando aparezca un caso real con orientación variable en la malla. Se fijan también la validación de admisibilidad propia de la ortotropía (`|ν₁₂| < √(E₁/E₂)`, no la de `Elastic2D`) y la convención de subíndices de Poisson.

**Consecuencia.** Calibración experimental desbloqueada sin riesgo de regresión. El isótropo queda como caso particular verificable.

## ADR 0014 — Ensamblaje por lotes

**Fecha**: 22 de septiembre de 2026. **Estado**: aceptado; fases 1-3 implementadas el mismo día.

**Contexto.** El ensamblaje recorría los elementos uno a uno en Python y, dentro de cada uno, los puntos de Gauss uno a uno: más del 90 % del tiempo era sobrecarga de despacho, no aritmética (120–170 µs por Quad4 y ensamblaje; 0,4 s por commit copiando diccionarios), y el estado interno como diccionarios ocupaba más que la matriz global. La solución conocida desde los supercomputadores vectoriales es invertir el bucle y evaluar por lotes los elementos que comparten receta, pero el proyecto exige que la agrupación no fije hipótesis que un componente futuro tenga que deshacer.

**Decisión.** (1) La *familia* se deriva de los contratos declarativos —clase de elemento, instancia de material, cuadratura, número de nodos—, nunca de listas. (2) El camino por elemento sigue siendo el contrato obligatorio y la referencia física; el camino por lotes es una capacidad que el componente declara (`BATCH_KINEMATICS` en el elemento, `STATE_SCHEMA` + `BATCH_KERNEL` en el material) con fallback automático dentro del mismo ensamblador. (3) Un único kernel de familia, compilado una vez y cacheado en disco, recibe la cinemática y la constitutiva como funciones tipadas por su firma. (4) Ambos caminos ejecutan las mismas funciones compiladas por punto de Gauss, y el barrido de contratos exige equivalencia a precisión de máquina. (5) El estado interno se declara (`STATE_SCHEMA`, obligatorio) y se almacena por arreglos con vistas compatibles con `ElementState`. (6) Nada derivado de la configuración deformada se cachea. (7) Trozos con presupuesto de memoria y mapa COO → CSR cacheado. (8) Cierre del mismo día: kernel paralelo con `prange` bit a bit igual al serie, señalización del jacobiano degenerado por bandera (una excepción dentro de una región paralela se pierde) y post-proceso por familia con el mismo kernel de cinemática y constitutiva.

**Consecuencia.** Quad4 + J2 ×27 y Hex8 + J2 ×13 en el ensamblaje con el kernel serie, ×110 y ×84 con el paralelo en 16 hilos; la exportación VTK ×6–8; el commit pasa de cientos de milisegundos a una copia de arreglo; el estado interno ocupa 80 B por punto de Gauss en vez de 462. Un material nuevo declara `STATE_SCHEMA` (una línea) y, si quiere acelerarse, un adaptador de diez líneas sobre su propio return mapping; un elemento nuevo escribe su cinemática compilada con las mismas funciones que su `compute_element_state`. Diferidos: los elementos estructurales 1D y los adaptadores sin asignaciones de Drucker-Prager y daño.

## ADR 0015 — Corrector de Newton compartido

**Fecha**: 22 de septiembre de 2026. **Estado**: aceptado; los cinco solvers migrados el mismo día.

**Contexto.** Los cinco solvers iterativos escribían el mismo bucle de corrección cinco veces (ensamblaje, residuo, criterio dual, backend con degradación a LU, Newton modificado, line search, historial de divergencia, commit), con dos variantes del line search y tres gestiones del backend. La auditoría de esa mañana tuvo que aplicar la misma corrección semántica en los cinco, y todo lo previsto (cargas con historia, acoplamiento, backend multihilo) los tocaría de nuevo.

**Decisión.** Un solo `NewtonCorrector` posee el bucle y lo que persiste entre pasos (backend, degradación a LU, factor congelado); cada solver construye por paso un `NewtonProblem` —ensamblar, residuo y su norma, escalas del criterio, corrección con el `solve` del corrector, actualización, commit— y conserva su control de paso. El iterado es opaco (`U`; `(U, λ, ΔU)`; `(u, u̇, ü, ü_libre)`). La restricción del arc-length es un método del problema que la variante por disipación sobreescribe (después también `IndirectDisplacementSolver`; desde el 2026-09-29 la variante por disipación vive en el módulo de usuario `dissipation_control` y el problema se exporta como `solidum.math.solvers.ArcProblem`). `CorrectionAborted` permite abandonar un paso sin marcarlo como tangente singular.

**Consecuencia.** Ningún resultado cambia (las 1 447 pruebas de entonces, sin tocar un valor esperado). Un solver iterativo nuevo escribe su problema y su control de paso; un cambio del bucle, del backend o de la globalización se hace una vez. Retirados los ayudantes privados duplicados; `solver.corrector` es la superficie de instrumentación.

## ADR 0016 — Acoplamiento termomecánico (estudio de viabilidad)

**Fecha**: 22 de septiembre de 2026. **Estado**: aceptado; sin código.

**Contexto.** El acoplamiento termomecánico era la decisión pendiente que el problema térmico (Etapa 8) había dejado abierta: acoplamiento débil o monolítico, y cómo entra en el contrato de material.

**Decisión.** Una verificación previa sobre el código mostró que la infraestructura ya admite nodos con desplazamiento y temperatura a la vez: numeración, ensamblaje en una sola matriz, condiciones mixtas, masa y capacidad por el mismo operador y resolución estática completa. Lo que falta es sólo la física de los bloques fuera de la diagonal. Se decide: **acoplamiento débil unidireccional primero**, con el monolítico como extensión aditiva; una **deformación propia genérica** `ε₀` con `σ = C:(ε − ε₀)`, de la que la térmica es un caso particular (la higroscópica, formalmente idéntica, domina en madera y bambú); `ε₀` almacenada **por punto de Gauss**; y sin familia de elementos acoplados.

**Consecuencia.** El acoplamiento deja de ser una decisión pendiente y pasa a ser una etapa de trabajo acotada que no necesita solver nuevo.

## ADR 0017 — Backend algebraico multihilo

**Fecha**: 22 de septiembre de 2026. **Estado**: aceptado.

**Contexto.** Medido con un perfilado de extremo a extremo, cerrado el ADR 0014 la factorización del sistema se llevaba el 97 % del tiempo de un análisis y el ensamblaje el 1 %. La causa es el relleno de los factores, que crece con la talla.

**Decisión.** Un backend opcional sobre Intel MKL Pardiso (factorización multihilo con reordenamiento por disección anidada), con el patrón de Cholesky: si la dependencia falta, el sistema funciona igual con SuperLU. Ocupa el lugar de LU en la regla automática (Cholesky, luego Pardiso, luego LU). Un único manejador de la MKL por proceso, y una factorización invalidada falla de forma ruidosa en lugar de devolver la solución de otro sistema.

**Consecuencia.** Aceleraciones de ×5 a ×56 según la talla, crecientes con ella, con la misma solución a precisión de máquina.

## ADR 0018 — Solver iterativo

**Fecha**: 23 de septiembre de 2026. **Estado**: aceptado; la tolerancia (10⁻¹⁰ sobre el residuo verdadero, fija y no configurable desde el YAML), validada por el usuario.

**Contexto.** Los programas comerciales ofrecen un solver directo por defecto y uno iterativo para modelos grandes. La ventaja del iterativo no es la velocidad sino la memoria: la del directo crece más deprisa que el modelo.

**Decisión.** Un backend iterativo **sólo a petición** (`linear_algebra: iterative`): gradiente conjugado para matrices simétricas definidas positivas, MINRES para simétricas indefinidas (el cambio es automático al detectar curvatura negativa) y solver directo para las no simétricas. El precondicionador es un multimalla algebraico que recibe los **modos de cuerpo rígido** del modelo, derivados del nombre de los grados de libertad; sin ellos rinde peor que no precondicionar. Tolerancia sobre el residuo verdadero cinco órdenes por debajo de la del Newton, de modo que el método de Newton no cambia; si no se alcanza, el análisis se detiene con la causa.

**Consecuencia.** Número de iteraciones casi independiente de la talla y memoria proporcional al modelo; por encima de unos 10⁵ grados de libertad gana también en tiempo.

## ADR 0019 — Red de seguridad del análisis estático

**Fecha**: 23 de septiembre de 2026. **Estado**: aceptado.

**Contexto.** El usuario de Solidum no elige el solver algebraico ni tiene por qué conocer sus fallos. Medido: con un modelo mal apoyado, los solvers directos devolvían en silencio desplazamientos de miles de kilómetros.

**Decisión.** Tres capas en el análisis estático: antes de resolver, los **mecanismos rígidos** se detectan a partir de los modos de cuerpo rígido que las restricciones dejan libres y se describen en términos del modelo (traslación, giro alrededor de un eje por un punto, campo escalar sin valor prescrito); después de factorizar y de resolver, en el análisis lineal, se verifican la ausencia de **mecanismos internos** (desde el 2026-10-20, por la energía del modo de menor energía y no por el tamaño de un pivote, con los elementos que se mueven nombrados) y que el solver **resolvió el sistema** (desde el 2026-10-20, por el error hacia atrás normado, que no depende de las unidades ni del contraste de rigideces). Dentro de un Newton no se rechaza nada —cerca de un punto límite resolver un sistema casi singular es legítimo—, sólo se mejora el diagnóstico.

**Consecuencia.** Un modelo mal planteado se detiene con una explicación en vez de devolver resultados absurdos, sin falsos positivos en toda la batería de pruebas.

**Capa 4 (27 de octubre de 2026, decisión del usuario).** Un sistema bien resuelto y no singular puede haber perdido cifras si está mal condicionado, y el residuo no lo delata: un voladizo de 10 000 elementos daba la flecha errónea en un 51 % sin aviso. El análisis lineal estima ahora el número de condición de la rigidez escalada por la diagonal con la factorización ya hecha (Higham y Tisseur; unas seis resoluciones) y avisa de las cifras perdidas por encima de $10^{12}$, cuando la solución baja a unas 5 cifras; el valor queda en el resultado.

## ADR 0020 — Módulos de usuario: las formulaciones no estándar fuera del programa principal

**Fecha**: 24 de septiembre de 2026. **Estado**: aceptado.

**Contexto.** Al decidir completar su formulación de discontinuidades interiores, el usuario fijó que un elemento no estándar no debe contaminar el programa principal, como en FEAP, que ofrece formulaciones estándar y deja al usuario añadir las suyas sin tocarlo. La discontinuidad embebida (ADR 0010) sí lo contaminaba: su familia cohesiva y su estado vivían en `core/`, su registro en `registry.py`, el lector YAML la trataba como caso especial (una sección y una clave propias) y la nombraban también la raíz del paquete, el descubrimiento automático, la herramienta de specs y el manual de referencia. La familia térmica sufría la misma rigidez: el lector repetía tres veces el mismo bloque, y los elementos térmicos resolvían su material con una heurística y no podían salir de gmsh.

**Decisión.** Dos términos, los de FEAP: **programa principal** (Solidum estándar, todo `solidum/` salvo `solidum/user/`) y **módulo de usuario** (una formulación no estándar en `solidum/user/<nombre>/`, análoga a los elementos y materiales de usuario de FEAP). La embebida pasa íntegra al módulo `discontinuities`. La carga es explícita, nunca automática: `user_modules: [discontinuities]` en el YAML o `solidum.load_user_module("discontinuities")` en Python; `import solidum` no carga nada de `solidum/user/`. El programa principal ofrece a los módulos piezas genéricas, cada una justificada por casos reales ajenos a la embebida: un `Registry` público con metadatos de familia (`YAML_SECTION`, `YAML_LABEL`, `SPEC_KIND`) y `Registry.families()`; un lector YAML con un solo bucle sobre las familias, en el que una sección de primer nivel desconocida es un error; `Element.REFERENCE_KWARGS` para las referencias de un elemento a otras familias; la carga de módulos; el gancho `prepare_step`, documentado como contrato genérico e idempotente; la cinemática pública (`compute_kinematics_tri3`, `compute_integrands`); la receta de elemento por grupo físico de gmsh; y la exclusión del camino por lotes de los elementos con estado propio. Quedan fuera a propósito, sin un segundo caso real, el servicio de seguimiento de trayectoria, un bus de eventos y la carga automática por *entry points*. Se descartaron un paquete Python aparte (dos instalaciones en cada PC y un nombre que sugiere otro programa), un repositorio aparte y la carga automática al importar Solidum (el YAML dejaría de decir qué formulación necesita).

**Consecuencia.** El programa principal no conoce la embebida, y un test de pureza (`tests/test_main_program_purity.py`) lo vigila: ningún archivo del programa principal nombra al módulo, y tras `import solidum` no hay ningún módulo de usuario cargado ni ninguna de sus clases registrada. Un módulo de usuario nuevo no toca el lector YAML, gmsh ni la herramienta de specs; la familia térmica gana validación de parámetros y referencias y puede salir de gmsh. Ningún resultado cambia. En contra: `from solidum import CST_Embedded2D` deja de funcionar, sin transición (el import es `from solidum.user.discontinuities import CST_Embedded2D`), los YAML de la embebida necesitan `user_modules`, y las piezas genéricas pasan a ser contrato estable del programa principal.

**Revisión del 2026-09-29.** El criterio de admisión de Reglas.md §1 (textos de referencia y código general establecido) fija la frontera, y cada spec del programa principal lo declara en su bloque `admission`. Por él sale del programa principal un segundo módulo, `dissipation_control`, con `DissipationArcLengthSolver`; el programa principal gana una sola pieza genérica, el nombre público `solidum.math.solvers.ArcProblem` para el problema por paso del arco, que ya extendía `IndirectDisplacementSolver`, y el lector YAML da la pista de `user_modules` también ante un tipo de solver desconocido. La frontera se endurece en los dos sentidos: el test de pureza pasa a tres controles (estructural por `ast`, léxico con los nombres derivados de lo que registra cada módulo, sin lista a mano, y en ejecución), `load_user_module` es atómica y sólo admite añadir entradas a los registros, y `tests/test_user_module_public_imports.py` exige a los módulos el contrato público: la cinemática se importa de `solidum.elements.solid_2d` (`compute_kinematics_tri3`, `compute_integrands`, sin guion bajo) y `YAML_LABEL` es la etiqueta de todo registro, en lugar del atributo privado `_kind`, retirado.

**Actualización del 10 de noviembre de 2026.** La familia cohesiva del módulo `discontinuities` volvió al programa principal (ADR 0033), como preveía Reglas.md §4, al entrar el primer elemento de interfaz del programa principal (la adherencia barra–sólido). El módulo conserva el elemento, el estado de la discontinuidad y el seguimiento de grieta.

## ADR 0021 — Secciones de marco: familia de secciones y sección por fibras

**Fecha**: 29 de septiembre de 2026. **Estado**: aceptado.

**Contexto.** La auditoría de pureza verificó un defecto físico silencioso en los cuatro marcos: evaluaban el material con la deformación axial y construían toda la rigidez local con la $E_t$ resultante, de modo que con fluencia axial el flector se anulaba aunque la sección estuviera curvada y con ablandamiento cambiaba de signo. El usuario decidió la física, una cuestión por vez: sección por fibras, elemento en desplazamientos, etapa 1 con los tres marcos Euler-Bernoulli (Timoshenko en una etapa 2), Gauss-Legendre con 3 puntos por defecto y mínimo 2, torsión elástica desacoplada con `GJ` de la sección y sin regularización del ablandamiento en esta etapa.

**Decisión.** Una familia nueva, `SectionRegistry` (sección YAML `sections`, `kind: section`), con la base abstracta `Section`, **sin estado**: el estado de cada fibra en cada punto de Gauss vive en el elemento. `FiberSection` es la primera sección. Los marcos Euler-Bernoulli aceptan exactamente uno de dos caminos: `material` + propiedades geométricas (forma cerrada, sin cambio de resultados, sólo con materiales `IS_LINEAR_ELASTIC = True`) o `section` + `n_ip` (elemento en desplazamientos integrado en Gauss-Legendre); la guarda de material elástico lineal protege también a `Frame2DTimoshenko`. Piezas genéricas nuevas, cada una con casos reales: `Material.IS_LINEAR_ELASTIC`, `Element.materials()` para las comprobaciones agregadas de masa, peso propio y simetría, la referencia `section` en `REFERENCE_KWARGS` (la primera con nombre distinto de `material`) y el gancho `from_yaml(data, resolve)` del lector YAML para objetos de familia con referencias anidadas. Se descartaron elementos nuevos por cada marco (duplicarían cinemática, catálogo y tests para cambiar sólo la ley de sección), la sección como material (contratos distintos) y la sección integrada dentro de cada marco (no se podría declarar una vez ni compartir).

**Consecuencia.** Los marcos dejan de dar flectores incorrectos en silencio y ganan plasticidad distribuida en la sección con cualquier material 1D, secciones mixtas y flexión biaxial. En contra: un modelo con marcos y material inelástico que antes corría ahora falla al construir, y el ablandamiento en fibras depende de la malla hasta que exista la regularización genérica.

## ADR 0022 — Cargas variables, mallas estructuradas y verificación por soluciones manufacturadas

**Fecha**: 27 de octubre de 2026. **Estado**: aceptado; el plan en cuatro fases y su criterio de aceptación (orden observado ≥ teórico − 0.1), aprobados por el usuario.

**Contexto.** La validación de Solidum se apoyaba en soluciones analíticas y benchmarks publicados, que cubren casos simples y no detectan un error que sólo afecte al orden de convergencia. La carga de borde de los elementos cuadráticos se repartía con pesos fijos sobre la cuerda (deuda #13) y no había presión.

**Decisión.** Las cargas distribuidas de todos los elementos (sólidos, térmicos, armaduras y marcos) aceptan una constante o un campo de la posición (`b(x)`, `t(x, n)` con la normal exterior, `Q(x)`, `q̄(x, n)`), integrado sobre la geometría real. Un contrato declarativo, `REFERENCE_CELL` y `REFERENCE_NODES` (y, en los elementos de línea, su campo interior y las rigideces de la sección), permite un generador de mallas estructuradas (`solidum.mesh`) y la verificación por soluciones manufacturadas con SymPy (`solidum.verification`): se inventa una solución, se derivan las cargas que la hacen exacta y se mide el orden de convergencia al refinar.

**Consecuencia.** Los 17 elementos lineales del programa principal dan el orden teórico; en los de línea, además, los valores nodales son exactos. Los barridos recorren el registro, así que un elemento nuevo queda verificado sin escribir pruebas nuevas.

## ADR 0023 — Cargas distribuidas y expresiones en el YAML

**Fecha**: 27 de octubre de 2026. **Estado**: aceptado; el convenio de signo de la presión ($p > 0$ comprime), decisión del usuario.

**Contexto.** Las cargas variables del ADR 0022 sólo estaban en el API de Python. El YAML no admitía tracciones de borde ni presión (deudas #23 y #24), y los grupos de gmsh se guardaban como conjuntos de nodos, de los que no se deduce con seguridad qué bordes cargar.

**Decisión.** Bloques `tractions` (por elemento) y `tractions_by_group` (por grupo físico de gmsh, con sus celdas de contorno exactas), con tracción `t` o presión `pressure` ($\mathbf t = -p\,\mathbf n$); `boundary_flux_by_group` en el térmico. Todo valor de carga puede ser una expresión de la posición (y de la normal en las de contorno), compilada por un intérprete propio con lista blanca, sin ejecutar el texto como código. Las cargas se resuelven y validan al leer el archivo.

**Consecuencia.** El YAML alcanza al API de Python en cargas distribuidas; un archivo ajeno no puede ejecutar código, y las erratas se detectan antes de resolver.

## ADR 0024 — Masas concentradas en nudos y análisis modal con grados de libertad sin masa

**Fecha**: 9 de noviembre de 2026. **Estado**: aceptado; las masas rotacionales en los nudos, decisión del usuario.

**Contexto.** La única fuente de masa era la densidad del material (ADR 0008), pero en un edificio la masa de cada nivel viene de la losa, los acabados y la carga viva, y la práctica estándar la concentra en los nudos. Con masas sólo en algunos nudos, o elementos de densidad cero, la matriz de masa es semidefinida y el análisis modal fallaba: ARPACK lanzaba un error sin causa y el camino denso exigía masa definida.

**Decisión.** La masa de nudo es un dato del dominio (`Domain.add_nodal_mass`, bloques YAML `nodal_masses`, `nodal_masses_by_group` y `nodal_masses_by_coord`), por grado de libertad mecánico, traslacional o rotacional. Se suma a la diagonal de la matriz de masa con los dos esquemas, sobre el patrón ya cacheado, así que la ven todos los consumidores de la masa. Es sólo inercia: `gravity` no genera su peso (decisión del usuario, como el `mass` de OpenSees: la masa sísmica y la carga de gravedad salen de combinaciones distintas), y el peso de los pisos se declara como carga. En el análisis modal, los grados de libertad sin masa se condensan estáticamente de forma exacta ($\boldsymbol\phi_0 = -\mathbf K_{00}^{-1}\mathbf K_{0m}\boldsymbol\phi_m$), sin formar la matriz condensada en los modelos grandes: el operador de desplazamiento e inversión del problema condensado es el bloque de $(\mathbf K - \sigma\mathbf M)^{-1}$ del modelo completo. Si el problema no tiene solución, `MasslessDofsError` nombra los grados de libertad. Se descartaron un elemento de masa puntual (material, rigidez nula y catálogo para una entrada diagonal), la masa ficticia pequeña (frecuencias espurias) y el problema inverso con $\mathbf K$ definida (falla en estructuras libres).

**Consecuencia.** Análisis modal y espectral de edificios con masas de piso, con frecuencias exactas para el modelo discreto; los modelos con masa en todos los grados de libertad no cambian.

**Ampliación (9 de noviembre de 2026, §6 bis, decisión del usuario).** Los transitorios implícitos (Newmark, HHT-α y sus variantes de Newton) admiten grados de libertad sin masa: la matriz efectiva $\mathbf M + \gamma\Delta t\,\mathbf C + \beta\Delta t^2\mathbf K$ es regular, y la aceleración inicial se toma consistente, $\ddot{\mathbf u}_0 = -\mathbf K_{00}^{-1}\mathbf K_{0m}\ddot{\mathbf u}_m$ en los grados de libertad sin masa. Un mecanismo sin masa da `MasslessDofsError`.

## ADR 0025 — Terminación de los solvers estáticos incrementales y registro por omisión

**Fecha**: 9 de noviembre de 2026. **Estado**: aceptado; amplía los ADR 0005 y 0011.

**Contexto.** Quien usaba los solvers de trayectoria desde Python no sabía, sin leer el registro, si el análisis terminó ni por qué se detuvo. Un trazado que agotaba `max_steps` regresaba con normalidad. El fracaso del arco era un `RuntimeError` genérico con el diagnóstico sólo en el registro, y un solver reutilizado conservaba el estado de la corrida anterior. El control indirecto no podía detenerse en un desplazamiento. El registro escribía las iteraciones en `INFO`, registraba como `ERROR` eventos que el solver recuperaba y producía trazas con una consola en cp1252.

**Decisión.** Atributo `termination` (`max_lambda`, `control_target`, `max_steps` o `diverged`) con `completed`, `lambda_final`, `steps_done`, `n_bisections` y `U_final`. Se fijan también tras una excepción, se reinician en cada `solve` y llegan a `SolveResult.termination`. El fracaso es una `SolverDivergedError` tipada: `StepConstraintError` si la restricción del paso no tiene solución. `IndirectDisplacementSolver` admite `control_target`: el último paso se recorta y la llegada es exacta. Ampliado el mismo día al arco cilíndrico: `ArcLengthSolver(..., control=[...], control_target=c)`, donde `control` sólo mide; el paso que cruza el objetivo se repite con la restricción del control indirecto y la llegada también es exacta. La base ofrece los puntos de extensión `_target_reached`, `_limit_step` y `_step_summary`. Registro: una línea `INFO` por paso, iteraciones en `DEBUG`, lo recuperado en `WARNING` con su motivo y `ERROR` sólo al detenerse; un manejador de consola adapta los mensajes a la codificación del flujo.

**Consecuencia.** Cómo terminó un análisis se lee en un atributo y un fracaso llega con su causa. Un análisis de empuje lateral se detiene en el desplazamiento pedido. Ningún resultado numérico cambia.

**Ampliación (9 de noviembre de 2026).** Una auditoría de generalidad del mismo día extendió el contrato a los análisis en el tiempo: `termination = "t_end"` o `"diverged"`, `completed`, `steps_done`, `t_final`, `U_final` y `partial_result`; una línea `ERROR` y una excepción tipada al detenerse (`TimeIntegrationInstabilityError` para la integración inestable), con el historial parcial. El θ-método dejó de devolver un resultado con columnas en cero al divergir. `control_target` toma el sentido del objetivo, así que una etapa también descarga.

## ADR 0026 — Etapas de carga en los solvers estáticos

**Fecha**: 9 de noviembre de 2026. **Estado**: aceptado.

**Contexto.** Los solvers estáticos arrancaban siempre en `U = 0` y su residuo no tenía término constante. No se podía aplicar la gravedad y empujar después, y una segunda corrida sobre un modelo plastificado empezaba en `U = 0` con el estado plástico de la primera.

**Decisión.** `solve(F_ref, u0=…, constant_load=…)` en `NonlinearSolver` y la familia del arco, con residuo $\mathbf F_c + \lambda\mathbf F_\text{ref} - \mathbf F_\text{int}$ y λ = 0 al empezar. Los valores impuestos crecen desde los de partida, $\mathbf g(\lambda) = \mathbf g_0 + \lambda(\mathbf g - \mathbf g_0)$. Cada etapa comprueba su equilibrio inicial, y un modelo con historia que empieza en cero lo avisa. `solidum.run(..., previous=resultado)` continúa la etapa anterior con su carga aplicada como constante y sus cargas equivalentes por elemento; el resultado separa `F_constant` de `F_variable`. `DissipationArcLengthSolver` rechaza las etapas: su restricción de disipación es física por decidir.

**Consecuencia.** Gravedad y empuje, P-Δ con la carga vertical y continuación de un análisis, verificados contra la viga-columna con acortamiento, la interacción N–M y la superposición elástica. Las etapas desde el YAML se añadieron el mismo día: bloque `stages`, una etapa por entrada con sus cargas y su solver, y `run_yaml` devuelve un `StagedResult` (deuda #36 saldada). Las dos cargas distribuidas en una misma etapa (deuda #35) las resolvió el ADR 0027.

## ADR 0027 — Cargas que llevan sus cargas equivalentes

**Fecha**: 9 de noviembre de 2026. **Estado**: aceptado; sustituye el §5 del ADR 0026.

**Contexto.** Las fuerzas internas de extremo con carga distribuida son $\mathbf F_\text{int} - \mathbf f_\text{eq}$, y el ensamblador recordaba sólo las cargas equivalentes de la última carga de cuerpo ensamblada: dos cargas de cuerpo en una corrida, dos casos de carga con el mismo ensamblador y el historial de los transitorios daban fuerzas internas erróneas sin aviso (deuda #35).

**Decisión.** La carga lleva su propia información: `LoadVector`, un `np.ndarray` con las cargas equivalentes por elemento de su parte distribuida, que la aritmética de las cargas combina linealmente (`F_a ± F_b`, `a·F`, `F/a`; cualquier otra operación da un `np.ndarray` simple). Se guardan concatenadas con una disposición compartida por las cargas del mismo ensamblador, así que combinar es una operación vectorizada. Un `np.ndarray` simple que las perdió recurre al registro anterior y avisa si es ambiguo.

**Consecuencia.** La deuda #35 queda cerrada en estático, por etapas y en el tiempo, sin cambiar la API de las cargas.

## ADR 0028 — Longitud característica del punto de integración

**Fecha**: 2026-11-09 · **Estado**: aceptado.

Un material con ablandamiento regularizado por energía de fractura necesita la longitud `h` de cada punto de integración, que sólo conoce el elemento. `h` viaja en el estado del punto: el material declara si la usa (`uses_characteristic_length`) y construye el estado inicial con `initial_point_state(h)`, donde la valida; el elemento, al terminar de construirse, siembra una longitud por punto con `_seed_characteristic_lengths` (`L` en las barras, `w_g·L` en los marcos con sección, que la reparte a sus fibras con `initial_states(characteristic_length=h)`). Un `h` inadmisible se rechaza al construir el elemento. Se descartó pasar `h` en cada `compute_state` (más de treinta llamadas y validación tardía) y guardarlo en el material (lo comparten puntos de longitudes distintas). Primer consumidor: `ManderConcrete1D`; regularizar los materiales de daño de los sólidos sólo exigirá que sus elementos siembren las longitudes.

## ADR 0029 — Continuidad de los análisis en el tiempo

**Fecha**: 2026-11-09 · **Estado**: aceptado.

Extiende las etapas del ADR 0026 a los análisis en el tiempo. `continue_from(previous)`, común a los seis integradores (mezcla `TransientTerminationMixin`), parte de un estático (sus desplazamientos o temperaturas, en reposo, con su carga constante y sus cargas equivalentes; aviso si el estado no la equilibra) o de un transitorio del mismo tipo (su estado final en su instante final, en tiempo absoluto). Al continuar un transitorio se hereda la aceleración final, no la del equilibrio instantáneo, que con HHT-α no coincide: parar y continuar es igual que no parar. Los resultados guardan `F_constant`. Pendiente de una decisión física: la tangente con que se forma el amortiguamiento de Rayleigh al continuar.

## ADR 0030 — Reacciones de los apoyos en los análisis dinámicos

**Fecha**: 2026-11-09 · **Estado**: aceptado.

Extiende el ADR 0004 §3 a los análisis dinámicos. Por decisión del usuario, la reacción es la ecuación de movimiento en las filas de apoyo, `R = F_int + M·ü + C·u̇ − F`, y la parte de amortiguamiento `C·u̇` se reporta aparte (OpenSees deja elegir entre ambas con `-rayleigh`). `TransientTrace` la registra en cada paso guardado, común a los cinco integradores mecánicos; el armónico la da en amplitudes complejas y el espectral por modo (`(K − ωₙ²M)·uₙ`), combinada con `combine_modal`. `total_reaction(dof)` da el cortante basal; en el espectral suma antes de combinar. Se guardan sólo los grados de libertad de apoyo, con el mismo reparto entre apoyos y uniones que el estático.

## ADR 0031 — Análisis lineales sobre el estado de una etapa

**Fecha**: 2026-11-09 · **Estado**: aceptado.

Por decisión del usuario, el modal, el armónico y el espectral que continúan una etapa (`continue_from`, `previous=`, etapas `stages`) linealizan alrededor de su estado con la rigidez tangente que dan los elementos (P-Δ de los corrotacionales, descarga elástica de la plasticidad, rigidez dañada), sin cambiarlo. La mezcla `LinearizedStateMixin` lo comparte. Los autovalores negativos dejan de recortarse a cero: `ModalResult.eigenvalues` con signo y `unstable_modes`. El `LinearSolver` como etapa es una perturbación lineal estática que tampoco cambia el estado (`SolveResult.perturbation`): casos superponibles sobre la etapa, con reacciones y fuerzas internas lineales (`Element.internal_forces_from_nodal`).

## ADR 0032 — Restricciones generadas por geometría y refuerzo embebido

**Fecha**: 2026-11-10 · **Estado**: aceptado.

Un generador de restricciones (`ConstraintRegistry`, tipo de spec `constraint`, bloque YAML `constraints`) crea sobre un dominio ya mallado los nodos, los elementos y las MPC (ADR 0004) que su geometría pide. La localización de puntos (`solidum/mesh/point_location.py`) invierte el mapeo isoparamétrico por Newton con lo que todo sólido declara (`REFERENCE_CELL`, funciones de forma), sin código por tipo de elemento. El primer generador es el refuerzo embebido (`EmbeddedReinforcement`): la barra, con su propia malla cortada en los bordes de los anfitriones, queda atada al campo interpolado del sólido. Sus MPC son de unión, y su fuerza es la adherencia (`EmbeddedBar.bond_forces`). `Domain.generated_nodes` aparta esos nodos de los selectores por coordenada del YAML.

**Ampliaciones (10 de noviembre de 2026).** (1) La barra admite elementos de tres nodos (`Truss2D3`, `Truss3D3`): el número de nodos se lee de `REFERENCE_NODES` y cada tramo recibe su nodo intermedio, atado como los demás. (2) Adherencia-deslizamiento: con `bond` el refuerzo no ata los nodos de la barra, sino que genera, uno por tramo, elementos de interfaz `BondSlipInterface` que unen la barra con su anfitrión y llevan una ley de interfaz de línea (ADR 0033). Esos elementos no se registran (los crea el generador, como crea las MPC), no tienen masa ni celda VTK (`VTK_EXPORT = False`, pieza genérica del exportador), y el lector YAML admite en una familia las claves que traduce su `from_yaml` (`YAML_KEYS`).

## ADR 0033 — Materiales de interfaz en el programa principal

**Fecha**: 2026-11-10 · **Estado**: aceptado.

Una sola familia de leyes tracción–desplazamiento relativo, `InterfaceMaterial` (`InterfaceMaterialRegistry`, tipo de spec `interface_material`, sección YAML `interface_materials`), para las dos geometrías de interfaz: superficie (grieta, junta: `(n, s)` / `(n, s1, s2)`) y línea (barra embebida: `(n, s)` / `(n1, n2, s)`), normales primero. El paso a ejes globales y el ancho de cada componente son del elemento. La familia cohesiva del módulo `discontinuities` pasa a ella sin cambio de física, como preveía Reglas §4, porque el programa principal tiene ahora un elemento de interfaz propio (la adherencia barra–sólido). La sección YAML `cohesive_materials` pasa a llamarse `interface_materials`, sin alias. Leyes de la familia: `CohesiveDamageIsotropic` (modo I) y `CohesiveDamageMixedMode` (modo mixto con apertura efectiva de Camacho–Ortiz, energía independiente del modo, ablandamiento lineal, exponencial o bilineal del fib Model Code 2010), de superficie; `BondSlipModelCode` y `BondSlipMultilinear` (adherencia-deslizamiento), de línea. Las dos leyes cohesivas comparten la envolvente de ablandamiento.

## ADR 0034 — Pseudo-tiempo de los solvers incrementales

**Fecha**: 2026-11-11 · **Estado**: aceptado.

El ensamblador lleva el intento de paso y su pseudo-tiempo para los materiales integrados con IMPL-EX. Los solvers que avanzan por pasos declaran en cada intento el incremento de su parámetro de control ($\Delta\lambda$, longitud de arco, $\Delta t$), en `prepare_all_steps(U, pseudo_time_increment=Δτ)` o, si lo conocen tras el predictor, con `set_pseudo_time_increment(Δτ)`. El ensamblador entrega a los materiales con `uses_pseudo_time` la razón $r = \Delta\tau_{n+1}/\Delta\tau_n$ respecto del último paso confirmado (`set_pseudo_time_ratio`) y refresca sus `batch_params`; al confirmar, la referencia se actualiza y los materiales quedan «sin paso», con el daño congelado en el confirmado. Repetir un intento tras una bisección sólo cambia $r$. Un solver que abre un intento sin declarar su pseudo-tiempo se rechaza en el primer commit si el modelo tiene materiales que lo usan. Se descartaron la razón en el estado de cada punto (no es historia), un argumento más en `compute_state` y en la firma compilada (rompería el contrato de todos los materiales) y una variable global. Primer consumidor: `IsotropicDamageCrackBand2D/3D` con `integration: impl-ex`.

## ADR 0035 — Mallado desde la geometría con Gmsh, y lector de Gmsh para todas las celdas

**Fecha**: 2026-11-11 · **Estado**: aceptado.

El lector de Gmsh importa todas las celdas de los sólidos de Solidum en 2D y 3D (`triangle`, `triangle6`, `quad`, `quad8`, `quad9`, `tetra`, `tetra10`, `hexahedron`, `hexahedron20`, `hexahedron27`). El orden de los nodos no sale de tablas de permutación por tipo: una tabla única da las coordenadas de referencia de los nodos de cada celda (las de Gmsh con la permutación de meshio) y el lector las casa con las `REFERENCE_NODES` del elemento destino (ADR 0022); un elemento con jacobiano negativo se reorienta reflejando su celda de referencia. Las líneas (2D) y caras (3D) son contorno, salvo que su grupo físico nombre un elemento (barras incrustadas). La clave `mesh` del YAML admite además un `.geo` o un bloque con `geometry`, `parameters`, `order`, `element_shape`, `serendipity` y `output`: Solidum llama a Gmsh, fija los parámetros del `.geo` (`DefineConstant`), malla, eleva a orden 2 si se pide y lee el `.msh`. Gmsh es una dependencia opcional (grupo `mesh`). Se descartaron las tablas de permutación a mano (fuente clásica de errores en Hex20 y Tet10) y leer el `.msh` con la API de Gmsh (haría obligatoria la dependencia). Salda la deuda #18.

## ADR 0036 — Salida de resultados del YAML: historia de la respuesta y frecuencia de VTK

**Fecha**: 2026-11-11 · **Estado**: aceptado.

La escritura de resultados pasa del script `examples/ejecutar_yaml.py` al paquete: `OutputWriter` se engancha al callback de paso del despacho completo de `run_yaml` (etapas, peso propio y cargas térmicas incluidos) y escribe, sólo si se pide (`run_yaml(..., write_output=True)` o `python -m solidum modelo.yaml`), el modelo sin deformar, los VTK por paso (`last_step`, `all_steps`, `{every: N}` o `{at: [factores]}`) con su `.pvd`, la historia de la respuesta en CSV (desplazamientos y reacciones por nodo, grupo físico o coordenada, fila a fila) y la salida de texto estilo FEAP. La reacción de un paso es `F_int(U) − factor·F` en los grados de libertad prescritos. Se descartaron guardar la historia en `SolveResult` (depende de lo que el usuario observa), escribir siempre que haya `output` (las pruebas ensuciarían `examples/`) y mantener la escritura en el script (despacho duplicado y reducido).

## Evolución de esta lista

Cada decisión de arquitectura de gran calado — refactor transversal, subsistema nuevo, ruptura de contratos — produce un ADR adicional. Los siguientes ADR se prevén en las fases de diseño futuras:

- El problema térmico se incorporó en la Etapa 8 **sin ADR propio**, a diferencia de lo previsto aquí. La razón es que no rompió contrato alguno ni exigió la generalización del concepto de grado de libertad que se anticipaba: la infraestructura resultó ya agnóstica al campo físico, y la familia térmica reutilizó el patrón de familia paralela con registro propio que el ADR 0010 había establecido para los materiales cohesivos. Las convenciones específicas —signo del flujo de frontera saliente, ausencia de notación de Voigt en el problema escalar— se registraron en `Reglas.md` §5.
- El acoplamiento termomecánico, que sí obliga a decidir sobre el contrato de material, tiene ya su ADR (0016): estrategia decidida, implementación pendiente.
