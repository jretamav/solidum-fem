# Evolución de la arquitectura

Este capítulo presenta la trayectoria prevista de la arquitectura de Fenix FEM. Su función es preservar la imagen de la arquitectura del programa a medida que crece: cada fase incorpora un nuevo tipo de análisis con sus contratos propios, sus variables de estado y, eventualmente, capas o subsistemas adicionales. La diferencia con los pendientes del capítulo 4 es de granularidad: aquellos enumeran componentes individuales; este capítulo describe los **estados sucesivos de la arquitectura** del programa, cada uno con su mapa de capas correspondiente.

Las fases posteriores a la actual están marcadas como pendientes y se materializarán mediante los Architecture Decision Records (ADR) anunciados en el capítulo 8. La estructura de cada fase se redefinirá en su ADR de apertura; lo aquí descrito constituye una hipótesis razonada de partida.

## Fase 1 (actual) — Estática lineal y no lineal

**Estado**: implementado.

**Capas activas**: las seis capas descritas en el capítulo 2 (entrada, inicialización, interpretación del caso, dominio, capa numérica, salida).

**Tipos de análisis**: estática lineal (`LinearSolver`) y estática no lineal incremental (`NonlinearSolver`, `ArcLengthSolver`).

**Variables de estado**: desplazamientos `U`. Variables internas del material (deformaciones plásticas, daño, historia) gestionadas mediante `ElementState` con semántica *trial* / comprometido.

**Matrices globales**: matriz de rigidez `K` (única).

**Contrato de elemento**: `K_local`, `f_internal`, `internal_forces(U)`. La matriz `B` y el residuo se calculan en cada iteración.

```
[Entrada] - [Inicializacion] - [Interprete] - [Dominio] - [Numerica: K, solvers] - [Salida]
```

## Fase 2 — Dinámica estructural transitoria

**Estado**: [PENDIENTE: Implementación de la fase 2 — dinámica estructural transitoria.]

**Cambios respecto a la fase 1**:

- Aparece la matriz de masa `M`, ensamblada con la misma infraestructura que `K`. Cada elemento añade un método `M_local` consistente con su esquema de integración.
- Aparece (opcionalmente) la matriz de amortiguamiento `C`, típicamente como combinación de Rayleigh `C = α·M + β·K`.
- Las variables de estado se amplían: además del desplazamiento `U`, se almacenan velocidad `Ú` y aceleración `Ü`.
- Se introduce una nueva familia de solvers para integración temporal: implícitos (Newmark, HHT-α, Bossak-α) y explícitos (diferencias centradas).
- El subsistema algebraico recibe matrices efectivas del tipo `K_eff = a·M + b·C + c·K`, donde los coeficientes dependen del integrador y del paso de tiempo. El despachador del ADR 0003 sigue aplicando: `K_eff` mantiene las propiedades de simetría y, a menudo, definición positiva de `K`.

**Topología prevista**:

```
[Entrada] - [Inicializacion] - [Interprete] - [Dominio] - [Numerica: K, M, C, solvers] - [Salida]
                                                                         |
                                                            [Integradores temporales]
```

```callout Riesgo de arquitectura
La generalización del concepto de DOF se mantiene (siguen siendo grados de libertad mecánicos), pero `ElementState` debe ampliarse para almacenar las derivadas temporales. La semántica *trial* / comprometido se generaliza al estado completo `(U, Ú, Ü, variables internas)`.
```

## Fase 3 — Análisis modal

**Estado**: [PENDIENTE: Implementación de la fase 3 — análisis modal.]

**Cambios respecto a la fase 2**:

- No hay paso de tiempo ni iteración no lineal; el problema es un autovalor generalizado `(K − ω²·M)·φ = 0`.
- Se introduce un nuevo subsistema en la capa numérica: el resolutor de autovalores. Algoritmos típicos: Lanczos para problemas grandes, Arnoldi para sistemas no simétricos, ARPACK como biblioteca de referencia.
- La salida cambia de naturaleza: en lugar de una historia temporal, se devuelve un conjunto de frecuencias propias `ω_i` y modos de vibración `φ_i`.
- El `SolveResult` se generaliza o se especializa: una alternativa es introducir un `ModalResult` separado; otra es mantener un único agregado polimórfico. La decisión queda diferida al ADR de apertura.

**Topología prevista**:

```
[Entrada] - [Inicializacion] - [Interprete] - [Dominio] - [Numerica: K, M, eigensolver] - [Salida: ModalResult]
```

```callout Riesgo de arquitectura
Bajo. El análisis modal reutiliza `K` y `M` ya disponibles desde la fase 2 y solo añade un nuevo tipo de algoritmo numérico.
```

## Fase 4 — Problema térmico estacionario y transitorio

**Estado**: [PENDIENTE: Implementación de la fase 4 — problema térmico.]

**Cambios respecto a las fases anteriores**:

- Aparece un campo escalar nuevo: la temperatura `T`. Esto obliga a generalizar el concepto de DOF: hasta ahora todos los DOF eran mecánicos (desplazamientos, rotaciones); ahora se incorporan DOF térmicos.
- Aparece la matriz de conductividad `K_t` (análogo térmico de la rigidez) y la matriz de capacidad `C_t` (análogo térmico de la masa).
- Se introduce una nueva familia de elementos térmicos, paralela a la mecánica. Cada elemento térmico declara `DOF_NAMES = ['T']` y consume materiales con propiedades térmicas (conductividad, capacidad, fuente).
- Aparece una nueva familia de materiales térmicos con interfaz propia: `compute_flux_and_tangent(grad_T, state)`.
- En transitorio, se reutilizan los integradores temporales de la fase 2, adaptados al sistema térmico (Crank-Nicolson, theta-método).

**Topología prevista**:

```
[Entrada] - [Inicializacion] - [Interprete] - [Dominio] - [Numerica: K_m, M_m, K_t, C_t, solvers] - [Salida]
                                                  |
                              [DOF mecanicos] + [DOF termicos]
                              [Materiales mec.] + [Materiales term.]
                              [Elementos mec.]  + [Elementos term.]
```

```callout Riesgo de arquitectura
Medio. Es la primera fase que rompe la suposición implícita "todos los DOF son mecánicos" presente hoy en numerosos puntos del código. La generalización requiere atributo `DOF_KIND` o similar en cada DOF (mecánico, térmico, futuro presión, futuro químico); posible separación entre `MechanicalElement` y `ThermalElement` con base común `Element`, o un único `Element` polimórfico (decisión a tomar en el ADR de apertura); adaptación del `VtkExporter` para escribir campos escalares además de vectoriales.
```

## Fase 5 — Acoplamiento termo-mecánico

**Estado**: [PENDIENTE: Implementación de la fase 5 — acoplamiento termo-mecánico.]

**Cambios respecto a la fase 4**:

- Aparece la interacción entre los campos térmico y mecánico: dilatación térmica que induce deformaciones, disipación mecánica que genera calor.
- Dos estrategias distintas en términos de arquitectura, a decidir en el ADR de apertura:
  - *Estrategia desacoplada (staggered)*: en cada paso, se resuelve primero el problema térmico con los desplazamientos de la última iteración, después el problema mecánico con las temperaturas recién calculadas. Cada subsistema mantiene su propia matriz; el acoplamiento se realiza mediante intercambio de campos. Implementación de bajo riesgo de arquitectura pero convergencia limitada en problemas de fuerte interacción.
  - *Estrategia monolítica*: se ensambla un sistema único `K_acoplada` que combina los DOF mecánicos y térmicos. Convergencia robusta pero exige una nueva infraestructura de ensamblaje multifísico.
- Aparece una nueva familia de materiales acoplados con interfaz extendida: `compute_response(strain, T, grad_T, state)`.

**Topología prevista (estrategia monolítica)**:

```
[Entrada] - [Inicializacion] - [Interprete] - [Dominio] - [Numerica: K_acoplada, M_acoplada, solvers] - [Salida]
                                                  |
                                       [Materiales termo-mec.]
                                       [Elementos termo-mec.]
                                       [Solver acoplado]
```

```callout Riesgo de arquitectura
Alto. Es la primera fase con auténtica multifísica. Requiere generalización del subsistema algebraico para manejar matrices con bloques no homogéneos (un bloque mecánico simétrico positivo definido y un bloque térmico no simétrico, en la formulación habitual); generalización del concepto de prueba de validación, ya que no basta una solución analítica monodisciplinaria y los problemas de referencia deben tener acoplamiento conocido; posible introducción de precondicionadores específicos para sistemas acoplados (precondicionamiento por bloques, complemento de Schur).
```

## Mantenimiento de este capítulo

Cada vez que se cierre la implementación de una fase y se acepte su ADR de apertura, esta sección se actualiza:

1. La fase pasa de "[PENDIENTE: ...]" a "implementado".
2. La hipótesis inicial sobre la arquitectura se sustituye por la descripción real, ajustada a las decisiones que se hayan tomado durante la implementación.
3. El diagrama de capas se redibuja conforme a la topología efectiva.
4. Las consecuencias no previstas (ajustes en capas anteriores, contratos modificados) se registran en una subsección "Lecciones de la fase N", de cara a fases posteriores.

De este modo, el capítulo se convierte en una crónica ordenada del crecimiento del programa: en cualquier momento futuro, abrir el manual permite reconstruir no solo el estado actual de Fenix FEM sino la trayectoria completa que lo ha llevado hasta él.
