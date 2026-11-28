# ADR 0044 — Cargas de línea en barras y marcos y fuerzas internas a lo largo del elemento

- **Estado**: aceptado (2026-11-27). Plan revisado por `arquitecto-solidum` (`informes/arquitecto/revision_plan_cargas_linea_2026-11-27.md`, «cambios antes de programar», adoptados) y spec por `teorico-solidum` (`informes/teorico/revision_spec_cargas_linea_2026-11-27.md`, «sí con cambios menores», aplicados); código verificado por el teórico (`informes/teorico/verificacion_teorica_2026-11-27_cargas_linea.md`, «sí con cambios menores», atendidos) y revisado por el arquitecto (`informes/arquitecto/revision_diff_cargas_linea_2026-11-27.md`, «cambios antes del commit», atendidos salvo h12 y h14 en parte).
- **Fecha**: 2026-11-27
- **Alcance**:
  - `solidum/elements/_line.py` (`LineElementMixin`): carga de línea, muestras de la carga, ejes de las fuerzas internas;
  - `solidum/core/element_forces.py`: fuerzas internas a lo largo del elemento a partir de las de extremo;
  - `solidum/core/load_vector.py`: la carga a lo largo de cada elemento viaja con la carga;
  - `Assembler`: `assemble_line_loads` y las muestras de la fuerza de cuerpo y del peso propio;
  - `YamlParser`: bloque `line_loads` (nivel raíz y etapas), cargas de las etapas resueltas al leer;
  - salida: fuerzas internas una vez por instante (VTK, historia, diagramas), historia con `at`, archivo de diagramas, figura `diagram`.
- **Física**: [Frame2DEuler](../specs/Frame2DEuler.md) §14–§15 (formulación común), [Frame3D](../specs/Frame3D.md) §13–§14, [Frame2DTimoshenko](../specs/Frame2DTimoshenko.md) §14, [Frame2DEulerCorot](../specs/Frame2DEulerCorot.md) §11, [Truss2D](../specs/Truss2D.md) §12, [Truss3D](../specs/Truss3D.md) §12. Plan aprobado por el usuario el 2026-11-27 («apruebo el plan A y B»).

## Contexto

Los marcos y las barras reparten bien una carga a lo largo del elemento (ADR 0022: `compute_body_load` con un campo `b(x)`, 4 puntos de Gauss y las funciones de la rigidez), pero sólo la reciben como fuerza de cuerpo por unidad de volumen: sirve para el peso propio. El YAML no tiene una carga por unidad de longitud propia de la viga (en kN/m, en ejes locales, trapecial, proyectada), y `tractions` sólo actúa en el contorno de los sólidos.

Las fuerzas internas salen sólo en los extremos de cada elemento, descontada la carga equivalente (ADR 0002 y 0027). Con una carga repartida el flector máximo de un claro cae entre nodos (`qL²/8` en la viga apoyada) y Solidum no lo da; tampoco hay diagramas. Para calcularlas dentro del elemento por estática del tramo hace falta, además de las fuerzas de extremo, **la carga misma a lo largo del elemento, con el factor con que actúa en el paso**; `LoadVector` (ADR 0027) lleva sólo los vectores nodales equivalentes.

La revisión del arquitecto midió además dos defectos de la salida del YAML que este bloque no puede heredar:
- **el VTK no descuenta la carga equivalente**: escribe las fuerzas de extremo de los marcos como fuerzas nodales $K u$ (voladizo con fuerza de cuerpo: $M_i$ = −6018.3 en el VTK frente a −6280 = −qL²/2 en `SolveResult` y en la historia);
- **la historia de un análisis en el tiempo sin etapas no rechaza** las fuerzas de marco con carga distribuida, como dice el ADR 0036: las escribe sin descontarla (Newmark: −286.48 frente a −351.90), porque el escritor se crea con las cargas nodales solas.

## Decisión

### 1. Contrato de los elementos de línea

En `LineElementMixin`, que ya comparten armaduras y marcos. Capacidades **declaradas**, que leen el lector, la historia y el archivo de diagramas (no se sondean con `try/except` ni por nombre de clase):

| Atributo de clase | Significado |
|---|---|
| `LINE_LOADS: bool` | recibe `line_loads` |
| `LINE_LOAD_LOCAL_AXES: bool` | admite la carga en ejes locales (`False` en `Truss3D`: sin ejes de sección) |
| `INTERNAL_FORCES_ALONG: bool` | da fuerzas internas a lo largo del elemento (`False` en `Truss2D3`/`Truss3D3`: su nodo intermedio aplica una fuerza dentro del tramo, pendiente) |
| `LINE_LOAD_POINTS`, `LINE_LOAD_WEIGHTS: tuple` | abscisas $\xi_g \in [0, 1]$ y pesos de la regla con que se muestrea e integra la carga (hoy, Gauss-Legendre de 4 puntos en todas; nadie fuera de la clase supone el «4») |
| `LINE_LOAD_BY_COMPONENT: bool` | `True` en las armaduras: reparto por componente global con sus funciones axiales |

Métodos:
- **Una sola declaración** del área y de la densidad lineal: `_load_area()` y `_linear_density()` (con sección, $\sum A_i$ y $\sum\rho_i A_i$). Las muestras de la fuerza de cuerpo ($A\,\mathbf b$) y del peso propio ($\rho A\,\mathbf g$), `body_load_samples(b)` y `self_weight_samples(g)`, se derivan de ellas en el mixin y ningún elemento las redefine. Los vectores nodales siguen con sus formas cerradas (no cambian los resultados); una prueba de barrido exige `compute_body_load(b) = line_load_vector(body_load_samples(b))` y lo mismo con el peso propio, en todo elemento con la capacidad.
- `line_load_samples(q, *, local=False, projected=False) -> (n, d)`: la carga por unidad de longitud en los `LINE_LOAD_POINTS`, en **ejes globales** (local → global con los ejes de referencia; proyectada, con $\sqrt{1 - t_k^2}$).
- `line_load_vector(samples) -> f`: el vector nodal consistente. Los marcos proyectan a sus ejes de carga y reparten con `_line_interpolation`; las armaduras, componente a componente con las funciones axiales de `_line_interpolation` — **un solo reparto por componente** en el mixin para las armaduras de 2 y de 3 nodos. `_field_line_load` pasa a ser `line_load_vector(area·b(x_g))` con las mismas operaciones (bit a bit en los marcos); en las armaduras cambia el redondeo (medido: 19 050 de 20 000 casos difieren en algún bit) y se compara con tolerancia declarada ($4\varepsilon$ relativa a la escala).
- `compute_line_load(q, ...) = line_load_vector(line_load_samples(q, ...))`.
- `internal_force_axes(U) -> (R, L0, l)`: los ejes locales (filas, en globales) y las longitudes de referencia y corriente con que se expresan las fuerzas internas. Los lineales, los de referencia y $l = L_0$; los corotacionales, los de la cuerda corriente. Lo usan `internal_forces_from_nodal` y los diagramas: los ejes de cada elemento viven en un lugar.
- `Frame3D._line_load_vector` (forma cerrada de la carga uniforme) pasa a `_uniform_line_load_vector`, para no confundirse con el del mixin.

Los cables no son `LineElementMixin`: admiten la fuerza de cuerpo constante y el peso propio con su reparto propio, aportan vector nodal pero no muestras, y no reciben `line_loads` en este bloque. El mixin y sus ganchos son **capacidad del programa principal**: hoy ningún módulo de usuario tiene un elemento de línea; cuando uno lo necesite, el contrato pasa a ser una pieza pública del ADR 0020 (ruta y nombres públicos, con su prueba), como las demás piezas (Reglas §1: segundo caso real).

### 2. Fuerzas internas a lo largo del elemento

Se construyen con **las fuerzas de extremo ya calculadas** (`ElementForces`, convención pública), no reevaluando el elemento: así coinciden con `element_forces` también en la perturbación lineal (ADR 0031) y en una etapa ya superada, no tocan el estado de prueba y no cuestan una evaluación constitutiva. Con los valores del nodo i, las muestras proyectadas a los ejes de `internal_force_axes` y las longitudes, una **tabla por familia** (`ElementKind`) en `solidum/core/element_forces.py`, junto a `_VALID_COMPONENTS`, da los términos de carga (specs §15 y §14):

| Familia | Fuerzas internas en $\xi$ ($s = \xi l$; $Q$ con $L_0$, brazos con $l$) |
|---|---|
| `frame2d` | $N = N_i - Q_{x'}$, $V = V_i + Q_{y'}$, $M = M_i + s V_i + R_{y'}$ |
| `frame3d` | $N = N_i - Q_{x'}$, $V_y = V_{y,i} - Q_{y'}$, $V_z = V_{z,i} - Q_{z'}$, $T = T_i$, $M_y = M_{y,i} + s V_{z,i} - R_{z'}$, $M_z = M_{z,i} - s V_{y,i} + R_{y'}$ |
| `truss` | $N = N_i - Q_{x'}$ |

$\tilde{\mathbf q}$ es el polinomio que interpola las muestras, así que las fuerzas internas son polinomios exactos. Resultado tipado `ElementForcesAlong(kind, xi, s, components)`. Para la salida, una evaluación **por lotes** de las estaciones uniformes (coeficientes con la inversa de Vandermonde de las abscisas, en caché; $Q$ y $R$ con bases de integración); los ceros del cortante, elemento por elemento con `polyroots`. Medido por el arquitecto con 2·10⁴ marcos: 82 µs por elemento y 583 B por elemento en el archivo (1.63 s por instante; las fuerzas de extremo, 12.3 µs por elemento); por eso el archivo tiene su propia frecuencia, sólo el final por omisión. Hacer por lotes los ceros (matrices compañeras apiladas) queda como mejora si hace falta.

### 3. La carga a lo largo del elemento viaja con la carga (ADR 0027)

`EquivalentLoads` pasa a tener **canales** con una sola álgebra (`scaled`, `plus`, serialización): el vector nodal (como hoy) y las muestras de los elementos de línea (`LoadSamples`), por elemento. `F = F_peso + λ·F_viento` lleva $\tilde{\mathbf q} = \tilde{\mathbf q}_\text{peso} + \lambda\,\tilde{\mathbf q}_\text{viento}$. Las muestras se guardan en globales porque no dependen del giro: cada elemento las proyecta a sus ejes. `samples_of(id)` devuelve el arreglo `(n, d)`.

La combinación de cargas con **disposiciones distintas** (una carga de línea en una parte de los elementos y el peso propio en todos) deja de pasar por diccionarios: un mapa de índices hacia la unión, en caché por disposición (`WeakKeyDictionary`: se suelta cuando la otra disposición muere). Medido por el arquitecto: con 1e5 elementos, 156 ms antes, 42 ms la primera vez y 4.3 ms en caché; con 1e4, 11.5, 3.6 y 0.06 ms.

`Assembler.assemble_line_loads(entradas)` devuelve un `LoadVector` con los dos canales y actualiza el registro heredado del ensamblador, como las cargas de cuerpo; `assemble_body_load` y `assemble_self_weight` añaden las muestras (coste medido: `assemble_self_weight` de 100 a 247 ms con 2·10⁴ marcos). `build_solve_result` las toma siempre de `F_applied`; `SolveResult.line_load_samples` es `{id: LineLoadSamples}` (abscisas `points` y valores, con `at(xi)`, ya con el factor aplicado) y `SolveResult.internal_forces_along(element, xi)` las usa. `as_load_vector` acepta las muestras. Si hay vector nodal de un elemento con `INTERNAL_FORCES_ALONG` pero no sus muestras (una carga que pasó por `np.asarray`), pedir sus fuerzas a lo largo del elemento es un error claro, nunca un diagrama sin la carga.

### 4. Bloque YAML `line_loads`

```yaml
line_loads:
  - {elements: {group: Viga}, q: [0.0, -20.0e3]}                    # globales [N/m]
  - {elements: {group: Faldon}, q: [0.0, -1.5e3], projected: true}  # por proyección
  - {elements: {ids: [4, 5]}, q_local: [0.0, "-0.6e3*(1 + 0.1*y)"]} # locales
```

- `elements` (obligatorio, sin «todos» implícito): `id`, `ids`, `where` (como la historia) o `group`, que aquí significa **las celdas del grupo físico** con coincidencia exacta por nodos, como las tracciones (ADR 0023 §2): error si ninguna celda coincide con un elemento de línea o si una coincide con dos. El selector se extrae de `output_writer` a `solidum/utils/element_selection.py`; la historia conserva su `group` por nodos.
- Una de `q` (globales) o `q_local`; `projected: true` sólo con `q`. Valores numéricos o expresiones de `x`, `y`, `z` (ADR 0023 §3), sin la normal.
- Errores al leer, todos juntos: claves, `q` y `q_local` a la vez o ninguna, `projected` con `q_local`, componentes ≠ dimensión, elemento sin `LINE_LOADS` o sin `LINE_LOAD_LOCAL_AXES`. Mensajes factuales (qué no se admite y por qué), sin recomendar alternativas (Reglas §3).
- **Las cargas de todas las etapas se resuelven al leer** (`parse()`), con el prefijo `stages[i].loads.`; hoy una tracción inválida de la etapa 2 se descubre tras resolver la etapa 1 (medido). Una sola declaración de las claves de carga distribuida y un solo `get_distributed_loads(assembler)` (fuerza de cuerpo o peso propio más cargas de línea) para el estático, las etapas y los transitorios.
- Se aplica donde la fuerza de cuerpo, con la misma regla: estático lineal y no lineal (crece con λ), etapas, transitorios (escalón desde $t = 0$); el armónico, el modal y el espectral no toman cargas distribuidas, tampoco éstas. Hoy no avisan de las que el archivo declara (la fuerza de cuerpo se ignora en silencio en el armónico); avisar cambia lo que Solidum le dice a quien lo usa, así que queda como pregunta para el usuario (Reglas §3), con la misma regla para todas las cargas distribuidas.

### 5. Salida

- **Una sola fuente por instante.** El escritor calcula una vez las `ElementForces` del paso con su carga equivalente (constante de las etapas anteriores más λ veces la de referencia) y con ellas escribe el VTK, la historia (`end`, `at`) y el archivo de diagramas. Corrige el defecto del VTK.
- **Análisis en el tiempo.** El escritor sabe si el archivo declara carga distribuida en todos los pipelines (bandera `distributed_loads`, de `YamlParser.has_distributed_loads`). Las fuerzas internas de marco con carga distribuida (`end`, `at`) y el archivo de diagramas se rechazan al leer en los análisis en el tiempo: la estática del tramo no incluye la inercia ni el amortiguamiento repartidos (spec Frame2DEuler §15). Corrige el defecto de la historia. El VTK decide por instante: en una etapa estática (también dentro de un análisis con etapas en el tiempo), las fuerzas del instante; en un instante en el tiempo con carga distribuida, `NaN` («no aplica», como la historia); sin carga distribuida, `internal_forces(U)`.
- **Historia**: `element: {id: 3, var: M, at: 0.5}`, `at` $\in [0, 1]$ desde el nodo i, excluyente con `end`; sólo en elementos con `INTERNAL_FORCES_ALONG`.
- **Archivo de diagramas**: `output.results.diagrams: {at: [factores] | final, stations: n}` (11 estaciones por omisión, sólo el instante final por omisión; frecuencia propia, no la del VTK) escribe `<nombre>_diagramas_<NNNN>.csv` y `<nombre>_diagramas_final.csv`: una fila por estación con `element, s, xi, x, y, z` y las componentes (vacías donde no aplican). Estaciones: `n` uniformes más los ceros del cortante dentro del elemento, donde el flector es extremo exacto. Tamaño: unas 13 filas por elemento e instante (≈ 1.3 kB por elemento); por eso no sigue la frecuencia del VTK.
- **Figura** `{type: diagram, var: M, at: [final], scale: auto, zoom, title, file}` (2D, como las demás): la estructura sin deformar y, sobre cada elemento, el diagrama perpendicular a su eje, con los extremos rotulados. Se dibuja del archivo de diagramas (con su índice `<nombre>_diagramas_indice.csv`, factor → archivo), que el escritor produce en los instantes que pide la figura aunque `results.diagrams` no lo declare. **El flector, del lado de la fibra en tracción** (decisión del usuario del 2026-11-27, spec Frame2DEuler, Diálogo, Q1), rotulado sin signo; `N` y `V`, con su signo y los positivos hacia unos ejes de dibujo fijos (arriba en las horizontales, izquierda en las verticales). Nada depende de la numeración de los nodos.

## Consecuencias

- Las cargas repartidas de vigas y marcos se declaran como se piensan (kN/m, ejes locales, proyectadas), y Solidum da el flector máximo dentro del claro y los diagramas, exactos por estática con cualquier material en los análisis estáticos.
- El VTK de los marcos con carga distribuida pasa a dar las fuerzas internas de extremo (antes, fuerzas nodales), y la historia de un transitorio deja de escribir fuerzas de marco sin descontar la carga.
- Un `LoadVector` combinado (casos de carga, etapas, función de carga de un transitorio) sigue sabiendo qué contiene, ahora también a lo largo de cada elemento, y combinar cargas de disposiciones distintas es vectorizado.
- **Coste del componente N+1**: un elemento con `LineElementMixin` declara sus capacidades, `_load_area`, `_linear_density`, `internal_force_axes` y sus ganchos de interpolación; nada en el lector ni en la salida. La skill `/solidum-new` gana la plantilla de elemento de línea.
- **Límites**: sin cargas puntuales dentro de un elemento ni momentos repartidos (un tramo cargado en parte o una carga concentrada van a un nodo); carga muerta, sin cargas de línea seguidoras; sin diagramas en las armaduras cuadráticas ni en los análisis en el tiempo; figura sólo 2D; cables sin cargas de línea. En una perturbación lineal (ADR 0031) sobre un elemento corrotacional, las fuerzas de extremo se proyectaron con los ejes del estado base y el diagrama usa los del resultado ($u_b + \Delta U$): el cierre de $V$ tiene un error de orden $|\Delta U|$ (medido por el teórico: 0.49·|ΔU| relativo a qL), y el de $M$ no se anula ni con los ejes del estado base, porque la rigidez tangente linealizada incluye términos geométricos que la estática del tramo no ve (medido por el arquitecto: M_j = 0 frente a 7.96 con |ΔU| = 0.35); sólo con corrotacional, perturbación y carga repartida, pendiente de decidir qué debe dar el diagrama en ese caso. `diagrams.at` con etapas: el índice guarda el factor sin la etapa, ambiguo si dos etapas pasan por el mismo factor.
- **Pruebas** (`tests/test_line_loads.py`, 69 más 2 omitidas): formas cerradas y casos de aceptación de las specs (los 13 del teórico); barrido sobre las 13 clases de línea registradas (con una guarda que exige incluir toda clase nueva): vector ↔ muestras (constante y campo, material y sección) y cierre en los extremos y derivadas con desplazamientos al azar; VTK = historia = diagrama; VTK de una etapa estática en un análisis mixto y `NaN` en el tiempo; rechazo en el tiempo (también con fuerza de cuerpo); álgebra de los canales, caché de la unión que se suelta y `np.asarray` que da error; errores de las etapas al leer (también con tracciones); `group` por celdas y `where`; bit a bit en los marcos; figura (lado del flector, ejes de N y V, invariancia a la numeración). Mutaciones: 30/30 (`informes/prototipos_2026-11-27/cargas_linea/mutar_cargas_linea.py`). La pureza del programa principal, sin cambios.

## Alternativas consideradas

- **Guardar la carga como función** (el campo `q(x)` y su factor). No se combina linealmente sin una lista creciente de términos, no se serializa con las expresiones compiladas y no es vectorizable en un transitorio de muchos pasos.
- **Derivar el vector nodal de las muestras** (una sola fuente). Cambia en el redondeo los resultados de las formas cerradas que hoy usan el peso propio y la fuerza de cuerpo constante; se prefiere conservarlos y proteger la igualdad con una prueba de barrido.
- **Coeficientes sobre «casos» ya ensamblados** (la carga como combinación de casos de carga con sus factores). Combinar costaría $O(\text{casos})$ en vez de $O(\text{elementos})$, pero es el modelo de cargas que el ADR 0027 dejó para el futuro; la fusión vectorizada basta hoy.
- **Diagramas reevaluando el elemento en `U`.** Discrepa de `element_forces` en la perturbación lineal y en etapas superadas, toca el estado de prueba y cuesta una evaluación constitutiva por consulta.
- **Diagramas con la interpolación del elemento** (la curvatura de Hermite; la fase 2 del ADR 0002, `sample_internal`). Es lineal en el elemento y no ve la carga; con fibras, las fuerzas de sección no cumplen el equilibrio punto a punto. La estática del tramo es exacta y sustituye esa fase.
- **Elementos subdivididos sólo para dibujar.** Cambia la malla y los resultados; no da el extremo exacto.

## Referencias

- ADR 0002 (fuerzas internas de extremo; su fase 2 queda sustituida por §2), 0020 (piezas para módulos de usuario), 0022 y 0023 (cargas variables y del YAML), 0027 (cargas con sus cargas equivalentes), 0026, 0029 y 0031 (etapas y perturbación), 0036, 0038 y 0039 (salida y figuras).
