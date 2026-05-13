# Decisiones de arquitectura (ADR)

Los Architecture Decision Records (ADR) son el registro persistente de las decisiones que han modelado la arquitectura de Fenix FEM. Cada ADR responde, en una página, a tres cuestiones: contexto, decisión y consecuencias. En este capítulo se resumen los ADR vigentes; el texto completo reside en `docs/adr/`.

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

**Contexto.** Fenix FEM se consume desde una GUI externa (FenixBAR, programa de análisis estructural con elementos barra). Tras una solución, el consumidor requiere información suficiente para visualizar el post-proceso estándar de Computer-Aided Engineering (CAE): deformada escalada, reacciones en apoyos y diagramas de esfuerzos internos sobre cada barra. Esta información no se exponía de forma uniforme: el componente `VtkExporter` cubría únicamente algunos elementos y los esfuerzos internos carecían de método público homogéneo.

**Decisión.** Interfaz pública estructurada en tres niveles:

1. *Contrato por elemento.* Todos los elementos tipo armadura, cable o marco implementan `internal_forces(U) -> ElementForces`, que devuelve los esfuerzos en ejes locales con claves homogéneas por familia (`{N}` para armadura y cable; `{N, V, M}` para marco 2D; `{N, Vy, Vz, T, My, Mz}` para marco 3D). La convención de signos es la del proyecto (capítulo anterior).
2. *Agregado en `Domain`.* El objeto `SolveResult` se calcula de forma anticipada al final de `solver.solve()` y es inmutable. Expone los desplazamientos globales, las cargas aplicadas, las reacciones y un diccionario de esfuerzos internos por elemento.
3. *Puntos de entrada oficiales.* Las funciones `fenix.run(case)` y `fenix.run_yaml(path)` devuelven directamente un `SolveResult`. Cualquier consumidor externo opera contra esta interfaz.

**Consecuencia.** El consumidor externo no recalcula esfuerzos a partir de `U`; la lógica de MEF reside íntegramente en Fenix FEM. La incorporación de un elemento nuevo obliga a implementar `internal_forces` con las claves de su familia: el contrato es explícito.

## ADR 0003 — Despachador de algoritmo algebraico

**Fecha**: 26 de abril de 2026. **Estado**: aceptado; fases 1 y 2 implementadas.

**Contexto.** Dentro de cada iteración del método de Newton-Raphson se requiere la resolución del sistema lineal `K · δU = R`. El algoritmo utilizado hasta el momento era SuperLU (factorización LU general), válido para todo caso pero subóptimo cuando `K` es simétrica y definida positiva, situación habitual en estática lineal y en numerosas no linealidades suaves. La factorización de Cholesky es típicamente entre 1.5 y 2 veces más rápida y consume considerablemente menos memoria en dichos casos.

**Decisión.** Se introduce un subsistema algebraico independiente del solver no lineal, con varios algoritmos disponibles y un despachador que selecciona el adecuado de forma automática:

- `LUSolver` (SuperLU) — algoritmo universal.
- `CholeskySolver` (CHOLMOD) — para matrices simétricas y definidas positivas; opcional según disponibilidad de `scikit-sparse`.
- `LDLTSolver` — reservado para la fase 2 (matrices simétricas indefinidas).

El despachador `select_solver(props, override)` recibe un objeto `StiffnessProperties` (atributos `is_symmetric`, `is_positive_definite`) y selecciona la factorización de Cholesky cuando aplica, o la factorización LU en caso contrario. El usuario puede forzar un algoritmo desde el archivo YAML mediante el campo `linear_algebra` en calidad de herramienta de diagnóstico, no como decisión de modelado.

**Consecuencia.** Los solvers no lineales (`LinearSolver`, `NonlinearSolver`, `ArcLengthSolver`) ya no invocan SuperLU directamente: solicitan un solver algebraico al despachador y se desentienden del algoritmo numérico subyacente. La incorporación de un nuevo algoritmo (Pardiso, MUMPS, métodos algebraicos multimalla iterativos) afecta exclusivamente a `fenix/math/linalg/`. La separación entre estrategia (solver no lineal) y táctica (subsistema algebraico) queda explícita en el código.

## ADR 0004 — Imposición de condiciones de Dirichlet por eliminación directa

**Fecha**: 4 de mayo de 2026. **Estado**: aceptado; fases 1 (Dirichlet nodal) y 2 (restricciones lineales afines) implementadas.

**Contexto.** Las condiciones de contorno de Dirichlet, tanto los empotramientos nodales clásicos como las restricciones multipunto (MPC: asentamientos prescritos, apoyos en plano oblicuo, periodicidad de celda unitaria, unión rígida master-slave), se imponían previamente por método de penalidad: añadir términos diagonales muy grandes al sistema. Dicho método introduce condicionamiento numérico artificialmente malo, contamina la simetría y la positividad de la matriz de rigidez y no es exacto a redondeo.

**Decisión.** Se adopta la imposición por eliminación directa. Toda restricción se expresa en forma afín `u_s = g_s + Σ α_si · u_mi` y se acumula en un `ConstraintSet` del subpaquete `fenix.bc`. El ensamblador construye un operador disperso `T` y un vector `g` tales que `u = T · u_libre + g`. El sistema entregado al subsistema algebraico es `K_red = TᵀKT`, `F_red = Tᵀ(F − K · g)`. Las restricciones encadenadas (un esclavo cuyo maestro es a su vez esclavo) se resuelven por cierre transitivo al construir `T`; los ciclos y redeclaraciones inconsistentes se detectan en validación temprana.

**Consecuencia.** La imposición es exacta a redondeo, preserva simetría y positividad definida de `K`, y deja el sistema reducido apto para los algoritmos del subsistema algebraico, incluidos solvers iterativos como gradiente conjugado o GMRES. La misma maquinaria cubre empotramiento, asentamiento, apoyo oblicuo, periodicidad y unión rígida sin ramificación específica por tipo.

## ADR 0005 — Logging configurable en lugar de `print` directo

**Fecha**: 7 de mayo de 2026. **Estado**: aceptado.

**Contexto.** La trazabilidad del solver durante un análisis (incrementos de carga, iteraciones de Newton, ratios de convergencia, factorizaciones reusadas) se emitía mediante `print` directo en la salida estándar. Esto impedía silenciar el ruido en pipelines de batch testing, redirigir la traza a un archivo de log paralelo y diferenciar niveles de gravedad (información ordinaria, advertencias, errores recuperables).

**Decisión.** Se introduce un *logger* configurable basado en el módulo estándar `logging` de Python, expuesto al resto del proyecto vía `fenix.logging.get_logger(name)`. Todos los componentes consumen este *logger* en lugar de `print`. El nivel y el formato se configuran de forma centralizada y pueden sobreescribirse desde el archivo YAML del caso o programáticamente.

**Consecuencia.** Tests automatizados pueden silenciar la traza sin perder la información de errores. La integración con herramientas externas (GUI FenixBAR, batería de tests, pipelines CI) se hace mediante la API estándar de `logging`, sin acoplamiento a `sys.stdout`. Las advertencias internas relevantes (degradación de Cholesky a LU, bisección de incremento) quedan marcadas con su nivel adecuado.

## ADR 0006 — Tolerancias del criterio de admisibilidad constitutiva

**Fecha**: 12 de mayo de 2026. **Estado**: aceptado.

**Contexto.** Cada modelo constitutivo no lineal tiene un criterio de admisibilidad — una función `f(estado)` que decide si el material está en régimen reversible (`f ≤ 0`) o disipativo (`f > 0`). En aritmética de doble precisión la comparación no puede ser exacta: tras un return mapping, el ruido de redondeo deja `f` cerca de cero pero no exactamente cero, por lo que la comparación se escribe `f ≤ tol`. La implementación anterior usaba una tolerancia absoluta hardcodeada (`PLASTIC_YIELD_TOL = 1e-9`), frágil bajo cambio de unidades y no adaptativa al estado en modelos con endurecimiento severo o sin parámetro de referencia inicial.

**Decisión.** Se adopta el esquema combinado absoluto + relativo, estilo solver de ecuaciones diferenciales ordinarias: `f ≤ ADMISSIBILITY_TOL_ABS + ADMISSIBILITY_TOL_REL · escala(estado)`. Cada material declara su escala característica vía el método `admissibility_scale(state_vars)`, devolviendo un esfuerzo positivo característico del criterio en el estado corriente: `σ_y + H·α` para plasticidad J2 1D, `E·κ_0` para daño isótropo, `√(2/3)·(σ_y + H·α)` para Von Mises 2D. La comparación se encapsula en `Material.is_admissible(f, state_vars)` y en `Material.admissibility_tol(state_vars)` para que los modelos no la repliquen.

**Consecuencia.** Misma física en MPa, Pa, GPa o adimensional produce idéntico estado interno hasta precisión de doble. En endurecimiento severo la tolerancia crece junto con la fluencia. Modelos energéticos futuros (phase-field, gradient damage) encajan declarando su escala (`√(2·E·g_c)`) sin tocar la política. Los modelos con kernel compilado (Numba) consumen `admissibility_tol` desde código Python antes de invocar el kernel, preservando la centralización de la fórmula.

## ADR 0007 — Tolerancias de convergencia en solvers no lineales

**Fecha**: 13 de mayo de 2026. **Estado**: aceptado.

**Contexto.** El criterio de convergencia del `NonlinearSolver` y del `ArcLengthSolver` era puramente relativo: `‖residuo‖ / ref_force ≤ tol` y `‖δU‖ / ‖U‖ ≤ tol`, con un único parámetro `tol` para los dos criterios y sin término absoluto real. Esta forma arrastraba tres limitaciones estructurales: no había tolerancia absoluta cuando la escala relativa colapsaba transitoriamente, una sola tolerancia mezclaba dos criterios físicamente distintos (incrementabilidad del iterado vs equilibrio), y la fórmula se replicaba en los dos solvers a mano.

**Decisión.** Se aplica el mismo patrón del ADR 0006 a la convergencia, separado por criterio. La política reside en una clase pequeña `ConvergenceCriterion` en `fenix/math/convergence.py`, que encapsula configuración (cuatro tolerancias: `rtol_force`, `rtol_disp`, `atol_force_factor`, `atol_disp_factor`) y estado calibrado (los `atol` efectivos derivados de las escalas del primer ensamblaje). Los solvers no lineales reciben una instancia en el constructor (parámetro `convergence`), calibran al inicio de `solve()` y delegan en `evaluate()` la comparación. La semántica es AND, no OR (ambos criterios deben cumplirse simultáneamente), coherente con la convención de Bathe, Crisfield y Owen-Hinton.

**Consecuencia.** Invariancia bit-paritaria bajo cambio de unidades (el mismo problema en N/m, kN/mm o MN/m converge en el mismo número de iteraciones). Régimen degenerado (carga inicial nula, control en desplazamiento puro) deja de oscilar espuriamente: el piso absoluto `atol_force` autoderivado define cuándo se da por convergido un residuo numéricamente cero. Problemas casi-rígidos y plásticos perfectos se afinan independientemente vía `rtol_disp` vs `rtol_force` sin compromisos. La política vive en un único punto del código; un solver no lineal futuro la consume sin replicarla.

## ADR 0008 — Densidad como propiedad del material

**Fecha**: 13 de mayo de 2026. **Estado**: aceptado.

**Contexto.** El cableado de fuerza de cuerpo introducido previamente aceptaba un vector `body_force: [bx, by, bz]` global, aplicable uniformemente a todos los elementos. Esto servía para estructuras monomaterial — el usuario calculaba `ρ · g` una vez y lo declaraba — pero fallaba en estructuras multimaterial: con un único vector global, no había forma de distinguir el peso propio del acero del hormigón o de un cable tensado. Adicionalmente, la densidad es la magnitud física que abre la puerta a la matriz de masa elemental y por tanto a todo el subsistema dinámico futuro.

**Decisión.** La densidad es propiedad del material, no del elemento. Se introduce como atributo `density` en la clase base `Material` con valor por defecto `None` (no declarado). La densidad es físicamente intrínseca al medio continuo; vive en el mismo nivel jerárquico que `E`, `ν` o `σ_y`. En el YAML se acepta el bloque atajo `gravity: [0, -9.81]` complementario a `body_force` (mutuamente excluyentes). El ensamblador incorpora `assemble_self_weight(g)` que itera sobre elementos y aplica `b_element = material.density · g`. Si algún material tiene `density=None`, el método falla con `ValueError` listando los materiales afectados — enforcement diferida al uso, sin posibilidad de propagar masa cero silenciosa a un análisis dinámico.

**Consecuencia.** Peso propio físicamente correcto en estructuras multimaterial. API del YAML más intuitiva. Tests existentes con `body_force` siguen pasando. La fase de análisis modal y dinámico hereda la propiedad sin rediseño: `Assembler.assemble_mass_matrix` reutiliza el mismo patrón de validación.

## ADR 0009 — Análisis modal y dinámico

**Fecha**: 13 de mayo de 2026. **Estado**: aceptado; fases 1 (modal), 3 (Newmark transitorio lineal) y 4 (Newton-Newmark transitorio no lineal) implementadas.

**Contexto.** El ADR 0003 había identificado análisis modal y dinámica implícita entre los análisis futuros que la capa algebraica iba a habilitar; el ADR 0008 había puesto la densidad como propiedad del material disponible. Quedaba abrir el subsistema con dos análisis distintos pero relacionados: el análisis modal — problema algebraico de valores y vectores característicos generalizado `K · φ = ω² · M · φ` —, cuyas soluciones admiten interpretación dinámica pero no constituyen análisis dinámico per se; y el análisis dinámico propiamente dicho, con dependencia temporal explícita — integración de `M · ü + C · u̇ + K · u = F(t)`.

**Decisión.** Se acometen ambos análisis con una hoja de ruta común y siete decisiones arquitecturales fijadas hoy para no reabrir debates en fases posteriores:

1. Matriz de masa **consistente** en la fase 1 (con contrato extensible a `lumped` mediante un parámetro `lumping` ya presente en la firma de `compute_mass_matrix`, suficiente para la futura integración explícita).
2. Ensamblaje de M paralelo al de K, reutilizando la topología COO cacheada del `Assembler` y la enforcement de densidad del ADR 0008.
3. `ModalSolver` registrado en el mismo `SolverRegistry` que los solvers estáticos: no se crea un `IntegratorRegistry` paralelo, manteniendo el patrón "un solver = un análisis = una clase registrada".
4. Dirichlet en el problema modal por eliminación directa (mismo mecanismo del ADR 0004); el término independiente `g` de restricciones lineales no aplica en modal.
5. Estado dinámico transitorio `(u, u̇, ü)` fuera de `Node`: vive en un `TransientResult` específico, preservando la pureza topológica/geométrica de la clase nodal y evitando que análisis estáticos paguen el coste de campos que no usan.
6. Amortiguamiento Rayleigh `C = α · M + β · K` como entrada estándar para la fase de dinámica transitoria, con coeficientes calibrables modalmente a partir de dos pares `(ξ, ω)`. Amortiguamiento modal por modo y otros esquemas (Caughey, local) diferidos.
7. Convención de unidades heredada del modelo (ADR 0008): el usuario es responsable de la consistencia; las tolerancias adimensionales del ADR 0007 garantizan invariancia bajo cambio de unidades.

La fase 1 (`ModalSolver` + ARPACK Lanczos con shift-invert), la fase 3 (`NewmarkSolver` con `(β, γ)` parametrizables, factorización única de la matriz efectiva, Rayleigh proporcional y cargas por callback Python) y la fase 4 (`NewtonNewmarkSolver` — subclase de `NewmarkSolver` que añade Newton-Raphson dentro de cada paso con jacobiano `J = M + γΔt·C + βΔt²·K_t`, convergencia dual del ADR 0007, Rayleigh constante calibrado con `K_0`) quedan implementadas y validadas con tests contra solución analítica y recuperación a paridad de bits del caso lineal en ausencia de plasticidad.

**Consecuencia.** El subsistema dinámico se cierra con el alcance descrito; las fases 2 (lumped), 5 (diferencias centradas), 6 (respuesta en frecuencia) y 7 (análisis espectral / sísmico modal) quedan tabuladas en la hoja de ruta del ADR sin reabrir las decisiones arquitecturales tomadas aquí. La clase `Node` conserva su semántica original.

## Evolución de esta lista

Cada decisión de arquitectura de gran calado — refactor transversal, subsistema nuevo, ruptura de contratos — produce un ADR adicional. Los siguientes ADR se prevén en las fases de diseño futuras:

- [PENDIENTE: ADR para la introducción del problema térmico, incluyendo la generalización del concepto de DOF a magnitudes escalares no mecánicas.]
- [PENDIENTE: ADR para el acoplamiento termo-mecánico, incluyendo la elección entre estrategia desacoplada y monolítica.]
