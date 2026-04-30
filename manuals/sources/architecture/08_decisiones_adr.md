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

## Evolución de esta lista

Cada decisión de arquitectura de gran calado — refactor transversal, subsistema nuevo, ruptura de contratos — produce un ADR adicional. Los siguientes ADR se prevén en las fases de diseño futuras:

- [PENDIENTE: ADR para la introducción del análisis dinámico transitorio (esquemas implícitos y explícitos), incluyendo el contrato de la matriz de masa por elemento y el tratamiento de amortiguamiento.]
- [PENDIENTE: ADR para la introducción del análisis modal (problema de autovalores generalizado), incluyendo la elección de algoritmo (Lanczos, Arnoldi).]
- [PENDIENTE: ADR para la introducción del problema térmico, incluyendo la generalización del concepto de DOF a magnitudes escalares no mecánicas.]
- [PENDIENTE: ADR para el acoplamiento termo-mecánico, incluyendo la elección entre estrategia desacoplada y monolítica.]
