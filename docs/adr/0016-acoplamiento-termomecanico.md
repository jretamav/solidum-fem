# ADR 0016 — Acoplamiento termomecánico: deformación propia genérica sobre un dominio de campos mixtos

- **Estado**: aceptado (estudio de viabilidad; **sin código**)
- **Fecha**: 2026-09-22
- **Alcance**: decide la **estrategia** del acoplamiento termomecánico que la Etapa 8 dejó explícitamente fuera (`docs/ROADMAP.md` §C2). No implementa nada. Fija: (a) qué grado de acoplamiento se adopta y en qué orden, (b) dónde vive la deformación propia `ε₀` y cómo entra en el contrato `Material`, (c) qué piezas de la infraestructura sirven ya tal cual y cuáles hay que tocar. Las decisiones de formulación concretas de cada componente irán en sus specs.

## Contexto

`Reglas.md §0` define Solidum FEM como un programa para "simulación de problemas mecánico y térmico, **acoplados** o desacoplados". El desacoplado se entregó en la Etapa 8 (conducción de Fourier pura, 2D y 3D, estacionaria y transitoria). El acoplado nunca se abordó: era el C2 excluido del alcance de aquella etapa, y desde entonces figura como continuación natural de la línea térmica y como decisión pendiente ("débil vs fuerte, y toca el contrato `Material`").

Dos hechos condicionan la decisión y no estaban disponibles cuando se difirió:

1. **El hallazgo de la Etapa 8**: la infraestructura resultó ser **agnóstica al campo**. El problema estacionario de conducción no necesitó solver nuevo — el `LinearSolver` mecánico lo resolvió sin cambios, porque el ensamblaje se articula sobre nombres de DOF y no sobre una semántica cableada de "desplazamiento".
2. **El ADR 0015**: los cinco solvers iterativos comparten ya un único corrector de Newton, y su cierre anotaba que "un residuo acoplado (u, T) es otro `NewtonProblem`".

Falta comprobar si esos dos hechos se sostienen con nodos que lleven `u` y `T` **a la vez**, que es lo que el acoplamiento exige y lo que nunca se ha ejercitado: hoy los modelos térmicos son dominios aparte, con sus elementos (`FLUX_DIM`), sus materiales (`ThermalMaterial`) y su DOF `T`.

## Verificación previa: qué admite hoy la infraestructura

Antes de decidir se construyó un dominio con los **mismos cuatro nodos** compartidos por un `Quad4` mecánico y un `Quad4Thermal` superpuestos, es decir, nodos con `ux, uy, T` simultáneos. Resultados medidos (sondas desechables, no incorporadas al repositorio):

| Pieza | Resultado | Lectura |
|---|---|---|
| `Domain.generate_equation_numbers` | 12 DOF: `{ux:0, uy:1, T:2}`, `{ux:3, …}`, … | Numera por **nombre de DOF** iterando `node.dofs`; un nodo multicampo no es un caso especial. **Sirve tal cual.** |
| `Assembler.assemble_system` | `K` de 12×12, bloques `K_uu` y `K_TT` poblados, bloque `K_uT` **nulo** | Ensambla ambos campos en **una sola matriz global** sin distinguirlos. El bloque de acoplamiento está vacío porque nadie lo escribe todavía — no porque la estructura lo impida. **Sirve tal cual.** |
| `Assembler.reduce` | BCs mixtas `ux`+`T` sobre el mismo nodo ⇒ `K_red` 10×10 | El `ConstraintSet` (ADR 0004) opera sobre índices globales, ajeno al significado físico del DOF. **Sirve tal cual.** |
| `assemble_mass_matrix` | Bloque mecánico traza = 15 700 (= 7850 kg × 2 DOF); bloque térmico traza = 3 611 000 = `ρ·c·V` exacto | **El mismo operador ensambla masa `M` y capacidad calorífica `C`**: `compute_mass_matrix` del elemento térmico devuelve la capacidad y el ensamblador no necesita saberlo. **Sirve tal cual.** |
| `solidum.run` + `LinearSolver` + `build_solve_result` | Análisis estático mixto resuelto de extremo a extremo; reacciones mecánicas correctas y reacciones térmicas `±2500 W` con suma nula (balance energético exacto) | El pipeline estático **ya resuelve hoy** un problema de dos campos, sin acoplar. `reactions_by_node` reporta ambos sin cableado. **Sirve tal cual.** |
| `VtkExporter` | Exporta `Displacements`, `Temperature` y `Flux` del mismo modelo | Ya cubre el campo mixto (`_TEMPERATURE_DOF`). *Corrige una nota de la memoria de la Etapa 8 que lo daba por pendiente: se saldó en la auditoría del 2026-09-22.* |

**Conclusión de la verificación**: el hallazgo de la Etapa 8 se sostiene con campos mixtos. Lo que falta para acoplar **no es infraestructura de grados de libertad, ensamblaje, restricciones, solución ni post-proceso** — todo eso funciona hoy. Lo que falta es exclusivamente **la física de los bloques fuera de la diagonal** y el contrato que la transporta.

## El problema, en sus términos

El sistema termomecánico completo, en forma semidiscreta, es:

```
[ 0  0 ] [ ü ]   [ 0  0 ] [ u̇ ]   [ K_uu  K_uT ] [ u ]   [ F_ext ]
[ 0  0 ] [ T̈ ] + [ 0  C ] [ Ṫ ] + [ K_Tu  K_TT ] [ T ] = [ Q     ]
```

con los cuatro bloques:

- `K_uu = ∫ Bᵀ C B dΩ` — rigidez mecánica. **Existe.**
- `K_TT = ∫ B_θᵀ k B_θ dΩ` — conductividad. **Existe** (Etapa 8).
- `C = ∫ ρc NᵀN dΩ` — capacidad calorífica. **Existe** (Etapa 8, vía `compute_mass_matrix`).
- `K_uT = −∫ Bᵀ C β N dΩ` — **dilatación térmica**: el campo `T` genera fuerzas mecánicas. *No existe.*
- `K_Tu = −T₀ ∫ Nᵀ βᵀ C B dΩ` (actuando sobre `u̇`) — **calor de deformación** (efecto termoelástico de Gough-Joule): la velocidad de deformación genera calor. *No existe.*

donde `β = C:α` es el tensor de esfuerzos térmicos y `α` el de dilatación.

La asimetría física decide la arquitectura. `K_uT` es un efecto de **primer orden**: una barra restringida con `ΔT = 100 K` desarrolla `σ = E·α·ΔT ≈ 250 MPa` en acero — el orden de la fluencia. `K_Tu` es de **segundo orden**: el calentamiento por deformación elástica es de décimas de kelvin salvo en régimen adiabático rápido (impacto, conformado a alta velocidad) o con disipación plástica sostenida. **El grueso del valor de ingeniería del acoplamiento está en `K_uT` solo.**

## Decisión

### 1. Acoplamiento **débil (unidireccional) primero**, monolítico como extensión aditiva

Se adopta el acoplamiento **secuencial unidireccional**: se resuelve el problema térmico, su campo `T` se entrega al problema mecánico como deformación propia, y **no** se retroalimenta. No es una simplificación de conveniencia: es la física correcta mientras `K_Tu` sea despreciable, que es el caso de todo problema cuasi-estático de termoelasticidad estructural — el dominio real del proyecto y del PAPIIT que lo motiva.

Se descarta escribir el monolítico ahora, pero **no se cierra la puerta**: la verificación anterior demuestra que el sistema mixto `(u, T)` ya se ensambla y se resuelve en una sola matriz global. Cuando aparezca el caso real que exija `K_Tu` (termoplasticidad adiabática, conformado), el paso al monolítico consiste en **escribir los dos bloques fuera de la diagonal y un `NewtonProblem` con residuo acoplado** (ADR 0015), sobre la misma infraestructura de DOF, ensamblaje y reducción que ya existe. No es un refactor: es un componente N+1.

**Consecuencia inmediata**: el acoplamiento débil no requiere ningún solver nuevo. Es el mismo hallazgo de la Etapa 8, aplicado una vez más.

### 2. La deformación propia es **genérica `ε₀`**, no térmica

Se confirma y se eleva a decisión arquitectural el criterio fijado por anticipado el 2026-09-10 (`ROADMAP.md`, `STATUS.md` §"Próximo hito"). La constitutiva pasa de `σ = C : ε` a:

```
σ = C : (ε − ε₀)
```

con `ε₀` una **deformación propia genérica** (*eigenstrain*): deformación impuesta que no genera esfuerzo por sí misma y que se resta de la total antes de evaluar la constitutiva. La térmica es el **caso particular** `ε₀ = α·ΔT`.

La razón está documentada y es un caso de uso real fuera de este repositorio: la **deformación higroscópica** (hinchamiento y contracción por humedad) es formalmente idéntica — deformación impuesta por un campo escalar difusivo — y en madera y bambú **domina sobre la térmica en unos dos órdenes de magnitud**. Un acoplamiento atado a la temperatura resolvería el efecto secundario, dejaría fuera el principal, y obligaría a refactorizar el contrato `Material` una segunda vez. Es el mismo argumento que llevó a `ThermalConduction` a aceptar conductividad **tensorial** desde el contrato: *la generalidad que no cuesta al diseñarla, cuesta mucho al añadirla después.*

Con `ε₀` genérica, el acoplamiento higroscópico no es un componente nuevo: es otro proveedor del mismo campo, y **no vuelve a tocar el contrato**.

### 3. `ε₀` viaja por el **estado del punto de Gauss**, no por el material ni por el elemento

Es la decisión de diseño no evidente, y la que condiciona el coste de todo lo demás. Tres ubicaciones posibles:

| Ubicación | Por qué se descarta / se adopta |
|---|---|
| **Atributo del material** | Descartada. El material es **compartido por instancia** entre todos los elementos de una familia (`family_key` agrupa por `id(elem.material)`, ADR 0014). Un `ε₀` en el material sería el mismo para todo el modelo — justo lo contrario de un campo. Rompería además el agrupamiento por lotes. Es el mismo error que la orientación material del ADR 0013, ya reconocido allí como atajo desechable. |
| **Atributo del elemento** | Descartada. `ε₀` varía **dentro** del elemento (el campo `T` se interpola con `N`); un valor único por elemento degrada el orden de convergencia y produce esfuerzos térmicos espurios a trozos constantes. |
| **Fila del estado por punto de Gauss** | **Adoptada.** Es donde ya viven las variables internas con historia (`STATE_SCHEMA`, ADR 0014). `ε₀` se evalúa por punto de Gauss a partir del campo nodal `T` interpolado, se almacena por filas como cualquier otra variable de estado, y el kernel por lotes la lee sin cambiar de forma. |

El contrato `Material` gana un canal para `ε₀` en la evaluación constitutiva. La firma concreta (parámetro explícito de `compute_stress` frente a variable del `STATE_SCHEMA`) se fija en la spec del primer material que lo implemente, con el caso real delante; lo que este ADR fija es que **el dato es por punto de Gauss** y que **su semántica es genérica**, no térmica.

**El default es `ε₀ = 0`**, y con él todo material existente conserva su comportamiento bit a bit. Los quince materiales del catálogo no se tocan: la deformación propia es aditiva sobre el contrato, no una ruptura.

### 4. La interpolación de `T` la hace el **elemento mecánico**, no un elemento acoplado nuevo

Se descarta crear una familia de elementos `Quad4ThermoMechanical` con `DOF_NAMES = ["ux","uy","T"]`. Duplicaría cada elemento del catálogo (diez familias), multiplicaría la matriz de combinaciones válidas y forzaría a un usuario que ya tiene malla mecánica y malla térmica a construir una tercera.

En su lugar, en el acoplamiento débil el elemento mecánico **recibe el campo `T` ya resuelto** y lo interpola con sus propias funciones de forma para evaluar `ε₀` en cada punto de Gauss. Requiere que ambas mallas compartan nodos —restricción que este ADR acepta y declara— y aprovecha que las funciones de forma del `Quad4Thermal` son literalmente las del `Quad4` (comparten `_shared.py`). Mallas no conformes (interpolación entre discretizaciones distintas) quedan **fuera de alcance**: es un subsistema propio y no hay caso real que lo pida.

### 5. Orden de implementación propuesto

1. **`ε₀` en el contrato `Material` con default cero** — aditivo, sin consumidor todavía; validado con un test de barra restringida contra `σ = −E·α·ΔT`.
2. **Proveedor térmico de `ε₀`** — el elemento mecánico interpola `T` nodal y evalúa `ε₀ = α·ΔT`. Primer acoplamiento real, unidireccional.
3. **Validación** contra solución analítica: barra restringida, bimetal, y disco/cilindro con gradiente radial (este último bloqueado por el mismo mallador curvo que difirió la validación de la Etapa 8 y NAFEMS LE10).
4. *(Abierto, sin fecha)* **Proveedor higroscópico** — mismo canal, campo distinto. No toca el contrato.
5. *(Abierto, sin fecha)* **Monolítico** — `K_uT`, `K_Tu` y `NewtonProblem` acoplado, cuando exista el caso real.

Los pasos 1–3 son una etapa de trabajo acotada. Cada componente lleva su spec según `Reglas.md §4`.

## Consecuencias

**Lo que NO hay que tocar** (resultado principal de la verificación): `Domain`, `Node`, la numeración de ecuaciones, el `Assembler`, el `ConstraintSet` y las BCs, la matriz de masa/capacidad, los solvers, `build_solve_result`, `SolveResult` y el `VtkExporter`. Todos operan ya sobre dominios de campos mixtos. **El acoplamiento débil no necesita ningún solver nuevo.**

**Lo que sí se toca**: el contrato `Material` (canal aditivo para `ε₀`, default cero) y los elementos mecánicos (interpolación del campo nodal auxiliar y evaluación de `ε₀` por punto de Gauss). Ambos cambios son aditivos.

**Coste del componente N+1**: un segundo fenómeno de deformación impuesta (higroscópica, de retracción, de crecimiento) es un proveedor de `ε₀` y **no vuelve a tocar el contrato**. Ése es exactamente el criterio de `Reglas.md §1`.

**Restricciones declaradas**: (a) mallas mecánica y térmica **conformes**, compartiendo nodos; (b) sin retroacción mecánica sobre el campo térmico mientras el acoplamiento sea débil; (c) propiedades independientes de la temperatura (`C(T)`, `k(T)` y `α(T)` quedan fuera — exigen Newton en el problema térmico, ya registrado como línea aparte en `ROADMAP.md`).

**Sobre `Reglas.md §0`**: cerrados los pasos 1–3, el "acoplados o desacoplados" de la identidad del proyecto deja de ser una promesa pendiente en su caso de ingeniería dominante.

## Alternativas descartadas

- **Monolítico desde el principio.** Resuelve `K_Tu`, que es de segundo orden y sin caso real en el proyecto; a cambio obliga hoy a un residuo acoplado, a una tangente no simétrica (el sistema `(u,T)` pierde simetría salvo escalado) y a revisar el despachador algebraico (ADR 0003). Se pospone hasta tener el caso que lo exija — con la infraestructura ya verificada como capaz de soportarlo.
- **Acoplamiento débil bidireccional iterado** (resolver térmico y mecánico alternadamente hasta convergencia del par). Es el punto intermedio clásico, pero sólo aporta sobre el unidireccional cuando `K_Tu` importa; y cuando importa, el monolítico converge mejor (el iterado de punto fijo degrada precisamente con el acoplamiento fuerte). No hay ventana donde sea la mejor opción.
- **`ε₀` como deformación térmica específica.** Descartada por el criterio del §2: el caso higroscópico es real, domina dos órdenes de magnitud en el material del usuario, y obligaría a refactorizar el contrato una segunda vez.
- **Familia de elementos acoplados** (`Quad4ThermoMechanical`, …). Descartada por el §4: duplica el catálogo de elementos para resolver un problema de transporte de datos.
- **`ε₀` en el material.** Descartada por el §3: el material es compartido entre elementos de una familia; un campo no cabe en él.
