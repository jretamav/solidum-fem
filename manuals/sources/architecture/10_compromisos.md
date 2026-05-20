# Compromisos de diseño

Los Architecture Decision Records (ADR) del capítulo 8 registran las decisiones de arquitectura que han sido **tomadas e implementadas**. Este capítulo registra el complemento simétrico: las decisiones que han sido **deliberadamente postergadas o rechazadas**, junto con la justificación de su no adopción. La función es preservar memoria sobre los caminos no recorridos, de modo que cuestiones ya resueltas no se reabran innecesariamente en sesiones futuras y que, si en algún momento las circunstancias del proyecto cambian, la revisión parta de un punto informado y no desde cero.

A diferencia de los ADR, las entradas de este capítulo no se materializan en archivos individuales en `docs/adr/`: viven exclusivamente aquí. Cada entrada se actualiza cuando las circunstancias que motivaron la decisión cambien sustancialmente.

## Paralelización con MPI

**Decisión**: postergada sin fecha.

**Contexto**. La paralelización mediante el estándar Message Passing Interface (MPI) sería el paso natural para llevar Solidum FEM a problemas de gran escala (mallas de millones de elementos, descomposición de dominio, ejecución en clúster).

**Razones de la postergación**. Un solo desarrollador trabaja en el proyecto y los problemas que motivan su uso (investigación en mecánica de sólidos, validación de modelos constitutivos, prototipado de algoritmos) se ejecutan satisfactoriamente en máquinas de un solo nodo. La introducción de MPI multiplicaría el coste de cualquier extensión por la necesidad de mantener consistencia de estado entre procesos, complicaría sustancialmente las pruebas (ejecución determinística reproducible bajo distintas particiones de malla) y exigiría rediseñar el ensamblaje y el subsistema algebraico. El compromiso entre coste de mantenimiento y beneficio inmediato no es favorable mientras el caso de uso del proyecto se mantenga en su perfil actual.

```callout Reapertura de la decisión
Procedería revisar este compromiso si Solidum FEM se aplicara de forma habitual a mallas mayores que 10⁶ elementos, o si se incorporasen colaboradores adicionales con experiencia en programación paralela.
```

## Adopción de PETSc o Trilinos como backend numérico

**Decisión**: rechazada para el alcance actual.

**Contexto**. PETSc y Trilinos son las dos bibliotecas estándar para álgebra lineal dispersa, solvers iterativos, descomposición de dominio y gestión de matrices distribuidas en programas de elementos finitos de gran escala. Su adopción daría acceso inmediato a multimalla algebraica, precondicionadores avanzados y soporte MPI maduro.

**Razones del rechazo**. Ambas son dependencias binarias pesadas con un proceso de instalación no trivial en Windows, plataforma de trabajo del desarrollador. Su modelo de programación es de matriz distribuida y orientación procedural, lo que choca con la arquitectura actual basada en `scipy.sparse` y orientación a objetos. La capacidad de Solidum FEM hoy se ajusta a lo que SuperLU y CHOLMOD ofrecen mediante `scipy.sparse` y `scikit-sparse`; añadir PETSc o Trilinos sin un caso de uso concreto que lo exija supondría una dependencia infrautilizada y un coste de mantenimiento desproporcionado.

```callout Reapertura de la decisión
Procedería si se adopta MPI (compromiso anterior) y, simultáneamente, los problemas tratados requirieran solvers iterativos con precondicionamiento avanzado. En ese escenario, PETSc se convertiría en la opción natural por encima de implementar dichos algoritmos en código propio.
```

## Arquitectura cliente-servidor o servicio web

**Decisión**: rechazada.

**Contexto**. Una arquitectura cliente-servidor expondría Solidum FEM como un servicio remoto al que la GUI externa (FenixBAR) y otros consumidores se conectarían vía protocolo HTTP, gRPC o similar. Es el modelo dominante en simulación industrial moderna como servicio.

**Razones del rechazo**. El consumidor actual (FenixBAR) es local, está bajo control del mismo desarrollador y no tiene requisitos de concurrencia ni de despliegue distribuido. La interfaz pública se ha estabilizado mediante el ADR 0002 con `SolveResult`, `solidum.run` y `solidum.run_yaml` en proceso, lo que cubre el caso de uso de consumidor externo sin la sobrecarga de un protocolo de red. La introducción de un servidor exigiría serialización de mallas y resultados, gestión de sesiones, autenticación si el servicio fuera multiusuario, y un ciclo de despliegue paralelo al del programa. Ninguna de estas necesidades existe hoy.

```callout Reapertura de la decisión
Procedería si Solidum FEM se ofreciera como servicio público a múltiples consumidores no relacionados con el desarrollo, o si la GUI se desacoplara hasta el punto de ejecutarse en máquinas distintas.
```

## Soporte de múltiples lenguajes en la interfaz pública

**Decisión**: rechazada para el alcance actual.

**Contexto**. La interfaz pública del programa (`solidum.run`, `solidum.run_yaml`, `SolveResult`) es accesible exclusivamente desde Python. Una alternativa sería exponer la misma interfaz en C, C++ o como biblioteca compartida llamable desde Fortran y MATLAB, lo que ampliaría el espectro de consumidores potenciales.

**Razones del rechazo**. El consumidor identificado (FenixBAR) está escrito en Python; el ecosistema de programas docentes y de investigación al que se dirige Solidum FEM utiliza Python como lenguaje dominante. Ofrecer un binding multilingüe duplicaría el esfuerzo de mantenimiento de la interfaz pública sin un caso de uso concreto que lo justifique. La integración con MATLAB, si fuera necesaria en algún momento, se realizaría mediante el mecanismo nativo de MATLAB para llamar a Python (`py.solidum.run_yaml(path)`), sin necesidad de modificar Solidum FEM.

```callout Reapertura de la decisión
Procedería si surgiera un consumidor estable escrito en un lenguaje distinto de Python y la solución *vía* el puente nativo de cada lenguaje no fuera viable.
```

## Sistema de plugins externo al repositorio

**Decisión**: rechazada.

**Contexto**. Algunos programas grandes de elementos finitos permiten que terceros publiquen materiales, elementos o solvers como paquetes externos que se descubren en tiempo de ejecución sin formar parte del repositorio principal. La alternativa es la actual: todo componente vive dentro del repositorio de Solidum FEM y se valida con sus pruebas.

**Razones del rechazo**. La filosofía del proyecto (ver capítulo 1) considera la transparencia y la auditabilidad del código como valor primordial. Un sistema de plugins externos abriría la puerta a componentes no auditados que podrían silenciosamente romper la convención de signos del proyecto, los contratos declarativos del catálogo o las hipótesis de validación contra soluciones analíticas. Mientras el desarrollo lo lleve un único equipo, la incorporación de cualquier componente nuevo pasa por el repositorio y por las pruebas correspondientes, lo que mantiene la trazabilidad completa.

```callout Reapertura de la decisión
Procedería si el proyecto creciera hasta tener una comunidad amplia de contribuyentes externos con necesidades especializadas no absorbibles en el catálogo principal. Hoy no es el caso ni se prevé que lo sea.
```

## Refactor a tipado estático estricto con `mypy --strict`

**Decisión**: postergada sin fecha.

**Contexto**. El código de Solidum FEM utiliza anotaciones de tipo en sus interfaces públicas pero no se valida con `mypy` en modo estricto. Adoptarlo proporcionaría garantías estáticas adicionales y detectaría algunas familias de errores antes de la ejecución.

**Razones de la postergación**. La librería `numpy` y, en menor medida, `scipy.sparse` tienen un soporte de tipado estático que ha mejorado mucho pero sigue dejando puntos donde `mypy --strict` produce ruido sin valor real. El esfuerzo de pulir todas esas fricciones no se compensa, en el alcance actual del proyecto, con el beneficio en detección de defectos. La política presente — anotaciones de tipo en interfaces públicas y validación de comportamiento mediante pruebas de física — se considera suficiente.

```callout Reapertura de la decisión
Procedería cuando el ecosistema NumPy/SciPy cierre las lagunas de tipado restantes y la activación del modo estricto deje de generar fricción artificial.
```

## Mantenimiento de este capítulo

Cada vez que se reabra una decisión y se modifique su estado (de "postergada" a "rechazada", o viceversa), o cuando una decisión postergada se materialice en un ADR formal, la entrada correspondiente se actualiza. Si la decisión pasa a aceptarse, su entrada se retira de este capítulo y se sustituye por el ADR correspondiente en el capítulo 8.
