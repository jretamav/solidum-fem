# Archivo de Entrada YAML

La interacción exclusiva con el motor de Solidum FEM se realiza mediante el archivo de entrada `.yaml`. Un archivo válido contiene los bloques siguientes (algunos opcionales según el caso de uso):

- `nodes` *o* `mesh` — geometría.
- `materials` — catálogo de leyes constitutivas instanciadas.
- `sections` — secciones de marco (sección por fibras), sólo si algún marco usa `section`; sus fibras referencian materiales de `materials`.
- `elements` — conectividad (omisible si se importa malla Gmsh).
- `boundary_conditions_by_node` / `_by_coord` / `_by_group` — restricciones de Dirichlet.
- `point_loads_by_node` / `_by_coord` / `_by_group` — cargas de Neumann puntuales.
- `solver` — tipo de solver y sus parámetros.
- `output` — frecuencia y formato de exportación de resultados.
- `user_modules` — módulos de usuario que el modelo necesita (sólo si usa una formulación no estándar).

Una sección de primer nivel que el lector no conoce es un error, no se ignora: una errata en el nombre de un bloque detiene la lectura en vez de dejar el modelo sin ese bloque.

## Geometría: `nodes` y `elements`

Para modelos pequeños o estructuras reticulares (armaduras, marcos, cables), la geometría se escribe directamente en el YAML:

```yaml
nodes:
  - {id: 1, coords: [0.0, 0.0]}
  - {id: 2, coords: [2.0, 0.0]}
  - {id: 3, coords: [4.0, 0.0]}

elements:
  - id: 1
    type: Frame2DEuler
    material: 1
    nodes: [1, 2]
    A: 0.01
    I: 8.33e-6
```

Las coordenadas pueden ser de 2 o 3 componentes; los elementos 1D 3D (`Truss3D`, `Cable3DCorot`, `Frame3D`) requieren las tres.

## Importación de Malla Gmsh: `mesh`

Para problemas continuos 2D, se sustituye el bloque `nodes`/`elements` por una referencia al archivo `.msh`:

```yaml
mesh: "geometria_placa.msh"
mesh_material: 1
mesh_thickness: 0.1
mesh_quadrature: "2x2"
```

| Parámetro | Descripción |
|---|---|
| `mesh` | Ruta al archivo `.msh` de Gmsh, relativa al YAML. |
| `mesh_material` | ID del material por defecto que se asocia a toda la malla. |
| `mesh_thickness` | Espesor por defecto (estado plano 2D). |
| `mesh_quadrature` | Regla de Gauss: `"2x2"` (completa, recomendada) o `"1x1"` (reducida; produce *hourglassing*). La masa consistente y la capacidad térmica se integran siempre con la regla completa, con independencia de la de `K`. Los elementos sin cuadratura (`Tri3`) la ignoran. |
| `mesh_element` | Opcional. Tipo de elemento de toda la malla, p. ej. `Quad4Thermal` para un modelo térmico. Sin él, cada celda gmsh da el elemento estándar: `quad` → `Quad4`, `triangle` → `Tri3`. El número de nodos del elemento debe coincidir con el de la celda. |

`mesh_material` se busca en la familia de material del elemento: en `materials` para los mecánicos, en `thermal_materials` para los térmicos.

### Mapeo Avanzado: `mesh_physical_groups`

Si la geometría tiene zonas con materiales, espesores o elementos distintos definidos como *Physical Groups* en Gmsh, se mapean explícitamente:

```yaml
mesh_physical_groups:
  "Zona_Hormigon":
    material: 1
    thickness: 0.20
    quadrature: "2x2"
  "Zona_Acero":
    material: 2
    thickness: 0.05
```

Un grupo también puede elegir su propio `element`, con cualquier parámetro que acepte ese elemento. Sin `element`, un grupo sólo admite `material`, `thickness` y `quadrature`: cualquier otra clave es un error, no se ignora.

## Materiales

Cada material se declara con un `id` único y un `type` registrado en el catálogo:

```yaml
materials:
  - id: 1
    type: Elastic1D
    E: 200.0e9
    density: 7850.0
  - id: 2
    type: VonMises2D
    E: 210.0e9
    nu: 0.3
    sigma_y: 250.0e6
    H: 2.0e9
    hypothesis: plane_strain
    density: 7850.0
```

Los parámetros aceptados dependen del material; consulte el capítulo *Catálogo de Modelos Constitutivos* para el catálogo completo.

## Secciones de marco: `sections`

Una **sección de marco** es la ley constitutiva de un marco a nivel de sección transversal: recibe la deformación axial y las curvaturas del eje en cada punto de Gauss y devuelve las fuerzas internas de la sección ($N$, $M_z$ en 2D; $N$, $M_y$, $M_z$ en 3D). Hoy existe un tipo, `FiberSection` (ADR 0021): la sección se discretiza en **fibras**, cada una con su propio material 1D y su propia historia, de modo que la fluencia o el daño avanzan fibra a fibra. Es la forma de modelar **plasticidad o daño en un marco**: con `material` + `A`, `I`, los marcos sólo admiten un material elástico lineal.

Una sección se declara una vez en el bloque `sections` y la comparten los marcos que la citan con `section: <id>` (`Frame2DEuler`, `Frame2DEulerCorot` y `Frame3D`; ver el capítulo *Catálogo de Elementos Finitos*). Las fibras se generan con **parches** y **fibras sueltas**:

[TABLA: Primitivas de una FiberSection.]
| Primitiva | Claves | Fibras que genera |
|---|---|---|
| `shape: rectangle` | `material`, `y: [y_min, y_max]`, `z: [z_min, z_max]`, `ny`, `nz` | Rejilla de `ny × nz` celdas iguales, una fibra en el centro de cada celda con el área de la celda. |
| `shape: circle` | `material`, `r_outer`, `r_inner` (default `0.0`, maciza), `center: [y, z]` (default `[0.0, 0.0]`), `n_rad`, `n_circ` | Corona dividida en `n_rad` anillos y `n_circ` sectores; una fibra en el centroide de cada sector, con su área exacta. Ángulos desde `+y` hacia `+z`; el primer sector empieza en 0. |
| `shape: layer` | `material`, `start: [y, z]`, `end: [y, z]`, `n`, `area`, `displaces` (opcional) | `n` barras de área `area` equiespaciadas de `start` a `end`, extremos incluidos (`n: 1`: una barra en `start`). |
| fibra suelta (en `fibers`) | `y`, `z`, `area`, `material`, `displaces` (opcional) | Una fibra. |

Además, `GJ` (opcional): rigidez torsional elástica $G J$ de la sección, **obligatoria si la usa un `Frame3D`**. La torsión es elástica y no interactúa con la flexión ni se degrada cuando las fibras fluyen.

**Coordenadas y signos.** `y`, `z` son coordenadas en los **ejes locales del elemento**, medidas desde el **eje del elemento** (la línea que une sus nodos), no desde el centroide: una sección excéntrica acopla axial y flexión. En 2D, $y$ es la dirección transversal del plano ($z$ sólo da el ancho); en 3D los ejes locales los fija `ref_vector`. $M_z > 0$ tracciona las fibras con $y < 0$ y $M_y > 0$ las fibras con $z > 0$. Materiales de fibra admitidos: cualquiera de `STRAIN_DIM = 1` (`Elastic1D`, `Elastoplastic1D`, `IsotropicDamage1D`…); una misma sección puede mezclarlos (concreto y acero).

**Ejemplo 2D mínimo** — voladizo de acero elastoplástico, sección rectangular de 0.3 m × 0.5 m en 20 capas, llevado más allá de la primera fluencia ($M(0) = 3.9$ MN·m entre $M_y = 3.125$ MN·m y $M_p = 4.6875$ MN·m):

```yaml
nodes:
  - {id: 1, coords: [0.0, 0.0]}
  - {id: 2, coords: [1.0, 0.0]}
  - {id: 3, coords: [2.0, 0.0]}
  - {id: 4, coords: [3.0, 0.0]}

materials:
  - {id: 1, type: Elastoplastic1D, E: 2.0e11, sigma_y: 2.5e8, H: 2.0e9, density: 7850.0}

sections:
  - id: 1
    type: FiberSection
    patches:
      - {shape: rectangle, material: 1, y: [-0.25, 0.25], z: [-0.15, 0.15], ny: 20, nz: 1}

elements:
  - {id: 1, type: Frame2DEuler, nodes: [1, 2], section: 1}
  - {id: 2, type: Frame2DEuler, nodes: [2, 3], section: 1}
  - {id: 3, type: Frame2DEuler, nodes: [3, 4], section: 1, n_ip: 4}

boundary_conditions_by_node:
  - {node_id: 1, ux: 0.0, uy: 0.0, rz: 0.0}

point_loads_by_node:
  - {node_id: 4, uy: -1.3e6}

solver:
  type: NonlinearSolver
  num_steps: 10
```

Con `nz: 1` el parche es una pila de capas: basta para la flexión en el plano. `n_ip` (puntos de Gauss a lo largo del elemento) vale 3 si no se declara.

**Ejemplo 3D con `GJ`** — voladizo de un tubo de acero (radios 0.15 m y 0.13 m) con carga transversal en $z$ y un momento torsor en la punta; `ref_vector: [0, 1, 0]` hace coincidir los ejes locales con los globales:

```yaml
nodes:
  - {id: 1, coords: [0.0, 0.0, 0.0]}
  - {id: 2, coords: [1.5, 0.0, 0.0]}
  - {id: 3, coords: [3.0, 0.0, 0.0]}

materials:
  - {id: 1, type: Elastoplastic1D, E: 2.0e11, sigma_y: 2.5e8, H: 2.0e9}

sections:
  - id: 1
    type: FiberSection
    GJ: 2.67e7                  # obligatoria con Frame3D
    patches:
      - {shape: circle, material: 1, r_outer: 0.15, r_inner: 0.13, n_rad: 2, n_circ: 32}

elements:
  - {id: 1, type: Frame3D, nodes: [1, 2], section: 1, ref_vector: [0.0, 1.0, 0.0]}
  - {id: 2, type: Frame3D, nodes: [2, 3], section: 1, ref_vector: [0.0, 1.0, 0.0]}

boundary_conditions_by_node:
  - {node_id: 1, ux: 0.0, uy: 0.0, uz: 0.0, rx: 0.0, ry: 0.0, rz: 0.0}

point_loads_by_node:
  - {node_id: 3, uz: -1.0e5, rx: 1.0e4}

solver:
  type: NonlinearSolver
  num_steps: 10
```

En el empotramiento resultan $M_y = +P L = 3.0\times 10^5$ N·m (tracción en las fibras con $z > 0$, por encima del momento de primera fluencia del tubo, $2.89\times 10^5$ N·m) y $T = 10^4$ N·m.

**Concreto armado: barras que desplazan el concreto.** Una sección de concreto armado combina un parche de concreto y capas de barras de acero dentro de él. Las fibras del parche cubren toda la sección, también el sitio que ocupan las barras, así que, si no se indica otra cosa, en la posición de cada barra el área cuenta dos veces: la fuerza axial de la barra es $(\sigma_s + \sigma_c)\,A_b$ en lugar de $\sigma_s A_b$, y lo mismo la rigidez, la masa y el peso propio. Es la práctica de OpenSees y el comportamiento por defecto. El error suele ser pequeño (cuantía de acero del 1 al 3 %), pero es sistemático. Para descontarlo, cada capa de barras o fibra suelta declara con `displaces` el material que desplaza, el del parche que la contiene:

```yaml
materials:
  - {id: 1, type: IsotropicDamage1D, E: 3.0e10, kappa_0: 1.0e-4, alpha: 0.99, density: 2400.0}   # concreto
  - {id: 2, type: Elastoplastic1D, E: 2.0e11, sigma_y: 4.2e8, H: 2.0e9, density: 7850.0}        # acero

sections:
  - id: 2
    type: FiberSection
    patches:
      - {shape: rectangle, material: 1, y: [-0.25, 0.25], z: [-0.15, 0.15], ny: 40, nz: 1}
      - {shape: layer, material: 2, start: [-0.21, -0.11], end: [-0.21, 0.11], n: 3, area: 5.1e-4, displaces: 1}
      - {shape: layer, material: 2, start: [0.21, -0.11], end: [0.21, 0.11], n: 2, area: 3.1e-4, displaces: 1}
```

En la posición de cada barra con `displaces`, la sección añade una **fibra de compensación** de área $-A_b$ con el material desplazado. Esta fibra se evalúa con la ley de ese material y con su propia historia (daño, fluencia), así que el concreto queda con su área neta $A_c - \sum A_b$ en las fuerzas internas, en la rigidez, en la masa y en el área bruta, también en régimen inelástico. `displaces` es opcional en cada primitiva (sólo `layer` y las fibras sueltas lo admiten) y es un id de `materials`: el del material de alguna fibra declarada de la sección, distinto del de la propia barra. Sin `displaces`, las barras se suman al concreto.

**Validación al leer el archivo.** Son errores, con el contexto de la entrada (p. ej. `sections[0] (id=7), patches[0].material: referencia a material inexistente (9)`): una clave desconocida en un parche, una `shape` inexistente, un material de fibra inexistente o de dimensión distinta de 1, y una sección sin fibras. También lo son un `displaces` que nombre un material inexistente, de dimensión distinta de 1, igual al de la barra o que no sea el de ninguna fibra de la sección, y un `displaces: true`: no es una bandera sino el id del material desplazado. Un marco que declara a la vez `section` y `material` (o `A`, `I`…) también es un error.

**Masa y peso propio.** Salen de las fibras: $\rho A = \sum \rho_i A_i$ y los momentos de inercia másicos $\sum \rho_i y_i^2 A_i$, $\sum \rho_i z_i^2 A_i$, con las fibras de compensación incluidas (con `displaces`, $\rho A = \rho_c (A_c - \sum A_b) + \rho_s \sum A_b$). Los marcos Euler-Bernoulli desprecian la inercia rotatoria de flexión, así que su masa sólo usa $\rho A$; los momentos de inercia másicos entran únicamente en la inercia torsional de `Frame3D` ($\sum \rho_i (y_i^2 + z_i^2) A_i$). Cada material de fibra necesita `density` si el análisis usa `gravity` o la matriz de masa.

```callout Advertencia
**Límites de la sección por fibras.** (1) La curvatura es lineal dentro del elemento: una rótula plástica se reparte sobre el elemento y la carga de colapso **converge por arriba** al refinar (en una viga biempotrada, +35 % con 4 elementos en el vano y +3.5 % con 32). Use varios elementos por miembro cerca de las rótulas. (2) Con fibras en **ablandamiento** (`IsotropicDamage1D`) el daño se localiza en uno o dos puntos de Gauss y la respuesta posterior al pico depende de la malla y de `n_ip`: no hay regularización. (3) Sin deformación por cortante: `Frame2DTimoshenko` no admite sección todavía.
```

## Módulos de usuario: `user_modules`

El **programa principal** de Solidum contiene los elementos, materiales y solvers de elementos finitos estándar. Las formulaciones no estándar viven aparte, en **módulos de usuario** (`solidum/user/<nombre>/`), a la manera de los elementos, materiales y macros de usuario de FEAP (ADR 0020). Un módulo de usuario puede aportar tipos de elemento, materiales, solvers y secciones de primer nivel propias, pero el programa principal no lo carga por sí solo: el modelo que lo necesita lo declara.

```yaml
user_modules: [discontinuities]
```

Los módulos declarados se cargan al principio de la lectura, antes de validar el resto del archivo, de modo que sus tipos y secciones se validan como los del programa principal. La carga es explícita a propósito: el archivo dice qué formulación necesita, y el mismo modelo no funciona en una sesión y falla en otra según lo que se haya importado antes. Si el YAML usa un tipo o una sección de un módulo sin declararlo, el error lo dice y lista los módulos disponibles; un nombre de módulo inexistente también es un error.

Hoy existen dos módulos de usuario:

- `discontinuities`: el elemento con discontinuidad embebida `CST_Embedded2D` y su familia de materiales cohesivos, con la sección `cohesive_materials` (capítulos *Catálogo de Elementos Finitos* y *Catálogo de Modelos Constitutivos*).
- `dissipation_control`: el solver de arco por disipación de energía `DissipationArcLengthSolver` (capítulo *Esquemas de Solución*).

Un modelo puede declarar varios: `user_modules: [discontinuities, dissipation_control]`.

## Condiciones de Frontera

Solidum FEM ofrece tres mecanismos para imponer restricciones cinemáticas (Dirichlet), elegibles según la facilidad operativa:

### Por nodo: `boundary_conditions_by_node`

```yaml
boundary_conditions_by_node:
  - {node_id: 1, ux: 0.0, uy: 0.0, rz: 0.0}
  - {node_id: 4, ux: 0.0, uy: 0.0}
```

### Por coordenada: `boundary_conditions_by_coord`

Escanea todos los nodos y aplica la restricción a aquellos cuya coordenada coincide (dentro de `tol`):

```yaml
boundary_conditions_by_coord:
  x_min: {ux: 0.0, uy: 0.0, tol: 1.0e-5}      # empotrar lado izquierdo
  rodillos: {coord: 'y', val: 10.0, uy: 0.0, tol: 1.0e-4}
```

Las claves predefinidas `x_min`, `x_max`, `y_min`, `y_max`, `z_min`, `z_max` usan los extremos detectados automáticamente del dominio; `coord` admite `'x'`, `'y'` o `'z'`. Los selectores funcionan igual con nodos 2D y 3D. Un selector que no casa con **ningún** nodo es un error de modelo y detiene la lectura con `ValueError` (no se ignora en silencio).

### Por grupo físico: `boundary_conditions_by_group`

Aplica restricciones a todos los nodos de un *Physical Group* de Gmsh:

```yaml
boundary_conditions_by_group:
  "Borde_Empotrado": {ux: 0.0, uy: 0.0}
```

### Restricciones multipunto (MPC): `linear_constraints`

Para apoyos en plano oblicuo, periodicidad de celda unitaria, uniones rígidas master-slave y simetrías no alineadas con los ejes globales, Solidum FEM admite restricciones afines lineales de la forma $u_s = g_s + \sum_i \alpha_{si} u_{m_i}$ (ADR 0004 fase 2). Se declaran en el bloque `linear_constraints`:

```yaml
linear_constraints:
  # Apoyo a 45 grados en el nodo 4: u_y = u_x.
  - slave: {node: 4, dof: uy}
    masters:
      - {node: 4, dof: ux}
    coefficients: [1.0]

  # Periodicidad: u_x del nodo 7 igual al del nodo 3.
  - slave: {node: 7, dof: ux}
    masters:
      - {node: 3, dof: ux}
    coefficients: [1.0]

  # Union rigida: u_x del satelite N9 sigue al maestro N5 con offset r_y=0.5.
  - slave: {node: 9, dof: ux}
    masters:
      - {node: 5, dof: ux}
      - {node: 5, dof: rz}
    coefficients: [1.0, -0.5]
    g: 0.0
```

Las restricciones se imponen por eliminación directa (no por penalización), de modo que la imposición es exacta a redondeo y preserva la simetría y la positividad definida del sistema. Las cadenas master-slave se resuelven automáticamente por cierre transitivo: un maestro puede ser a su vez esclavo de otra restricción, o tener un apoyo, y entonces el esclavo hereda su valor prescrito.

Cada entrada admite sólo las claves `slave`, `masters`, `coefficients` y `g`, y cada referencia sólo `node` y `dof`; una clave con errata (`G: 5` por `g: 5`) es un error, no se ignora. Los coeficientes y `g` deben ser números. Al leer el archivo, antes de ensamblar nada y también con malla gmsh, se comprueba además que:

- cada nodo existe y tiene el grado de libertad nombrado (`ux` no existe en un nodo térmico);
- el esclavo no tiene además un apoyo: un grado de libertad no puede ser esclavo de dos restricciones. Si el apoyo es coherente con la restricción (por ejemplo, los dos nodos coincidentes de una interfaz atada tomados por un mismo apoyo por coordenada o por grupo), basta con quitarlo del esclavo, porque la restricción le transmite el valor del maestro;
- el esclavo no figura entre sus propios maestros ni se declara dos veces con restricciones distintas (declarar dos veces la misma es inofensivo);
- las cadenas master-slave no forman ciclos ($a \leftarrow b$, $b \leftarrow a$).

Todos los problemas se informan juntos, cada uno con su entrada (`linear_constraints[i]`) y el par (nodo, dof) afectado. Con la API Python (`Domain.add_linear_constraint`) el ensamblador detecta los mismos errores al construir las restricciones, de uno en uno, y su mensaje también nombra el esclavo y los maestros por (nodo, dof).

En los solvers estáticos no lineales (`NonlinearSolver`, `ArcLengthSolver`, `IndirectDisplacementSolver`) el término independiente `g` se aplica proporcionalmente al factor de carga, como $\lambda\,g$, igual que un desplazamiento prescrito: en un trazado de arco sigue a $\lambda$ también cuando la carga baja.

**Reacciones y fuerzas de restricción.** Una restricción multipunto ejerce una fuerza sobre la estructura en los grados de libertad de su esclavo y de sus maestros. Solidum distingue dos papeles:

- **Apoyo**: la restricción impide un movimiento de sólido rígido. Es el caso del apoyo a 45 grados del ejemplo ($u_y = u_x$ en un mismo nodo): su fuerza es una **reacción**, normal al plano de rodadura, y aparece en `R` y en `reactions_by_node` con las componentes `ux` y `uy` del nodo. Una relación con coeficiente distinto de uno entre nodos distintos ($u_s = 2\,u_m$) también cuenta como apoyo.
- **Unión**: una traslación rígida de toda la estructura (o una temperatura uniforme, en el problema térmico) la cumple. Son uniones el atado de mallas, la periodicidad (también con `g` distinto de cero) y la unión rígida del ejemplo. Su fuerza es una **fuerza de restricción**, que transmite fuerza entre partes de la estructura y no la toma del exterior: aparece en `constraint_forces` y `constraint_forces_by_node` del resultado, con valores iguales y opuestos en los nodos atados, y no en las reacciones.

Un apoyo con valor prescrito reporta su reacción completa, incluida la fuerza que le llega por los nodos atados a él: atar dos nodos coincidentes con $u_s = u_m$ en todos sus grados de libertad da las mismas reacciones que fundirlos en uno. Con esta convención se cumple el equilibrio global $\sum R + \sum F = 0$ por dirección.

## Cargas Puntuales

Misma estructura que las condiciones de frontera, en bloques `point_loads_by_node`, `point_loads_by_coord` y `point_loads_by_group`. Los valores son fuerzas nodales (en las unidades del modelo, típicamente Newtons):

```yaml
point_loads_by_node:
  - {node_id: 3, uy: -15000.0}

point_loads_by_group:
  "Borde_Carga": {ux: 1500000.0}
```

Las cargas son **estrictas**: un `node_id` inexistente se rechaza en la validación del archivo, y un nombre de DOF que el nodo no tiene (`fy` en vez de `uy`, `rz` sobre una armadura) o un grupo físico inexistente detienen la lectura con `ValueError`. Antes se ignoraban en silencio y el análisis corría sin carga.

En los **solvers transitorios mecánicos** (`NewmarkSolver`, `HHTSolver`, `NewtonNewmarkSolver`, `NewtonHHTSolver`, `CentralDifferenceSolver`) las cargas del archivo —puntuales más peso propio o fuerza de cuerpo— se aplican como **escalón constante en el tiempo** desde $t = 0$ (`F_func = cte`). Una historia $\mathbf F(t)$ arbitraria requiere construir el solver desde Python con su propio `F_func` (§10).

## Fuerzas de Cuerpo

Una *fuerza de cuerpo* es una fuerza por unidad de volumen $\mathbf b$ que actúa sobre todo el dominio material del elemento, no solo sobre su frontera. Cada elemento del catálogo (armaduras, cables, marcos 2D/3D, sólidos 2D) integra la contribución consistente $\int \mathbf{N}^\top \mathbf{b}\,A\,dx$ y la acumula al vector global de fuerzas. Para vigas Hermite cúbicas el reparto incluye fuerza nodal $qL/2$ y momento $\pm qL^2/12$; para barras axiales es mitad-mitad por nodo.

Casos típicos en mecánica estructural: peso propio (un caso particular con $\mathbf b = \rho\,\mathbf g$), fuerzas centrífugas en cuerpos rotatorios ($\mathbf b = \rho\,\omega^2\,\mathbf r$), fuerzas electromagnéticas sobre materiales conductores, expansión térmica modelada como pseudocarga. Solidum expone dos formas declarativas en YAML que cubren todos estos casos.

### Forma general: `body_force`

Aplica un vector $\mathbf b$ uniforme (fuerza por unidad de volumen, en ejes globales) a todos los elementos del modelo:

```yaml
body_force: [0.0, -77008.5]          # 2D: bx, by en N/m^3
```

o, para problemas 3D:

```yaml
body_force: [0.0, 0.0, -77008.5]     # 3D: bx, by, bz
```

Es la forma genérica y la única necesaria para fuerzas de cuerpo que no son peso propio (campos electromagnéticos precalculados, centrífuga estática, etc.). También sirve para peso propio cuando todo el modelo es monomaterial y el usuario prefiere precalcular $\rho\cdot g$.

Cada componente puede ser también una expresión de la posición (ver «Valores Variables: Expresiones»), por ejemplo una fuerza centrífuga $\mathbf b = \rho\,\omega^2\,\mathbf r$:

```yaml
body_force: ["7850*100*x", "7850*100*y"]   # ρ·ω²·r con ω = 10 rad/s
```

La admiten los sólidos, las armaduras y los marcos; los cables, no (reparten la carga con su longitud sin tensión).

### Caso particular — peso propio: `gravity` + `density`

El peso propio es el caso especial $\mathbf b = \rho\,\mathbf g$. Para problemas multimaterial — donde cada elemento debe ver *su propio* $\rho$ — Solidum expone un atajo declarativo: la densidad como propiedad del material, y un vector global de gravedad.

```yaml
materials:
  - {id: 1, type: Elastic1D, E: 210.0e9, density: 7850.0}     # acero
  - {id: 2, type: Elastic1D, E:  30.0e9, density: 2500.0}     # hormigón

gravity: [0.0, -9.81]                # m/s^2, +y hacia arriba
```

**Qué hace Solidum**. Para cada elemento computa $\mathbf b_e = \rho_e\cdot\mathbf g$, donde $\rho_e$ es la densidad del material asignado, e integra la contribución consistente con la misma maquinaria de `body_force`. Resultado: estructuras multimaterial dan peso propio físicamente correcto, cada elemento con su propio $\rho$.

**Densidad solo cuando se usa**. El campo `density` es **opcional** al declarar un material (ADR 0008): los análisis estáticos puros que no invocan peso propio ni matriz de masa no necesitan declararla. Pero al invocar `gravity` (o, en su día, análisis modal/dinámico), Solidum exige la densidad de cada material involucrado: si algún material no la trae declarada, el cálculo falla con `ValueError` listando los materiales afectados por nombre. Imposible propagar masa cero silenciosa. Si un material legítimamente carece de masa física (penalty, restricción), declararlo explícitamente como `density: 0.0` — caso en el cual `gravity` emite un `WARNING` informativo y su contribución al peso propio es nula.

### Detalles transversales

**Exclusividad**. `body_force` y `gravity` son dos formas declarativas para expresar fuerzas de cuerpo; declarar ambas en el mismo archivo lanza error de validación. La forma general (`body_force`) es la subyacente; `gravity` es atajo para su caso particular más común.

**Padding 2D↔3D**. Si el modelo mezcla elementos 2D y 3D, declarar el vector con tres componentes; el ensamblador pasa el slice apropiado a cada elemento según la dimensionalidad declarada por su `DOF_NAMES`.

**Reservado para análisis dinámico**. El atributo `density` introducido aquí es la misma magnitud que consume la matriz de masa $\mathbf M_e = \int \rho\,\mathbf N^\top\mathbf N\,dV$ en análisis modales y dinámicos. Declararlo prepara el modelo para esas extensiones sin coste adicional (ADR 0008).

## Cargas Distribuidas de Contorno: `tractions` y `tractions_by_group`

Una *tracción* $\mathbf t$ es una fuerza por unidad de superficie sobre la frontera del sólido (en 2D, por unidad de longitud de borde y de espesor). Se integra con la geometría real del borde o de la cara, $\int_\Gamma \mathbf N^\top\,\mathbf t\;d\Gamma$, también sobre bordes curvos de elementos cuadráticos. La admiten los sólidos 2D y 3D.

**Por grupo físico** (mallas gmsh): las celdas de contorno del grupo (las líneas de gmsh en 2D) son los bordes exactos que se cargan.

```yaml
tractions_by_group:
  Borde_Carga: {t: [1.0e6, 0.0]}        # tracción en ejes globales [Pa]
  Agujero:     {pressure: 2.0e6}         # presión normal [Pa]
```

**Por elemento** (cualquier malla, también 3D): el borde (`edge`, 2D) o la cara (`face`, 3D) con la numeración local del elemento.

```yaml
tractions:
  - {element: 12, edge: 1, t: [1.0e6, 0.0]}
  - {element: 40, face: 3, pressure: 2.5e3}
```

Cada entrada lleva **una** de dos formas:

- `t`: el vector de tracción en ejes globales, con 2 ó 3 componentes según la dimensión del elemento.
- `pressure`: la presión $p$, normal al contorno. **Convención de signo** (Reglas §5): $p > 0$ es **compresión**, la presión empuja hacia el interior del cuerpo, $\mathbf t = -p\,\mathbf n$ con $\mathbf n$ la normal exterior. Una succión se declara con $p < 0$. Es el convenio de Abaqus y ANSYS.

La carga es **muerta**: actúa sobre la geometría de referencia y no sigue la deformación. Un grupo que no sea de contorno (por ejemplo, una superficie de la malla) o un borde compartido por dos elementos se rechazan al leer el archivo.

## Valores Variables: Expresiones

Cualquier valor de `body_force`, `tractions`, `tractions_by_group` y de las cargas térmicas (`body_source`, `boundary_flux`, `boundary_flux_by_group`) puede ser, en vez de un número, una expresión entre comillas de la posición `x`, `y`, `z`; en las cargas de contorno, también de la normal exterior `nx`, `ny`, `nz`:

```yaml
body_force: [0.0, "-78.5e3*(1 + 0.1*y)"]          # campo b(x)
tractions_by_group:
  Agujero: {pressure: "1.0e6*(1 - y/0.2)"}         # presión variable
  Borde:   {t: ["1.0e6*nx", 0.0]}                  # con la normal
```

Se admiten números, `+ - * / **`, paréntesis, la constante `pi` y las funciones `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `atan2`, `sinh`, `cosh`, `tanh`, `exp`, `log`, `log10`, `sqrt`, `abs`, `min` y `max`. Nada más: la expresión no se ejecuta como código Python, así que un archivo ajeno no puede hacer otra cosa que calcular una carga, y una errata (`w` en lugar de `y`, una función que no existe) se detecta al leer el archivo, con la entrada y lo que no se admite. Una expresión que no se puede evaluar en algún punto (`log` de un negativo) se informa con el punto. Las cargas que varían en el tiempo no se declaran en el YAML (ADR 0023).

## Solver

```yaml
solver:
  type: NonlinearSolver
  num_steps: 20
  max_iter: 15
  adaptive: true
  convergence:
    rtol_force: 1.0e-5
    rtol_disp: 1.0e-5
```

Los solvers disponibles y sus parámetros se documentan en los capítulos *Esquemas de Solución* y *Análisis Dinámico*.

El solver del sistema lineal que cada uno resuelve internamente **lo elige Solidum** según las propiedades de la matriz; no hace falta declararlo. El campo opcional `linear_algebra` sirve para diagnóstico o para pedir el solver iterativo en modelos muy grandes:

```yaml
solver:
  type: LinearSolver
  linear_algebra: iterative    # auto (default), lu, pardiso, cholesky, iterative[:amg|jacobi|none]
```

Un valor desconocido, o `iterative:amg` sin `pyamg` instalado, se rechaza al leer el archivo.

### Bloque `convergence`

Los solvers no lineales (`NonlinearSolver`, `ArcLengthSolver`) aceptan un bloque `convergence` que configura el criterio de parada de las iteraciones de Newton. La política sigue el patrón estándar de tolerancias mixtas absoluta + relativa, separado por dos criterios físicos: *fuerza* (residuo de equilibrio) y *desplazamiento* (incremento del iterado). La forma exacta de la comparación es:

```
||R_red|| <= atol_force + rtol_force * max(||F_ext||, ||F_int||)
||dU||    <= atol_disp  + rtol_disp  * ||U||
```

`R_red` es el residuo de equilibrio en los grados de libertad libres: las reacciones de los apoyos no cuentan y, con restricciones multipunto, el residuo de cada grado de libertad esclavo se suma a sus maestros multiplicado por su coeficiente ($\mathbf T^\top\mathbf R$, ADR 0004).

Ambos criterios deben cumplirse **simultáneamente** (semántica AND) para dar el paso por convergido.

**Parámetros disponibles** (todos opcionales, con defaults razonables):

- `rtol_force` (default `1.0e-5`): tolerancia relativa del residuo de fuerza. Es el parámetro que normalmente ajustarás: bajarlo a `1.0e-7` para análisis de precisión, subirlo a `1.0e-3` para exploración rápida.
- `rtol_disp` (default `1.0e-5`): tolerancia relativa de la corrección de desplazamiento. En problemas casi-rígidos puede afinarse independientemente del de fuerza.
- `atol_force_factor` (default `1.0e-9`): factor adimensional que define el piso absoluto de fuerza. La tolerancia absoluta efectiva se autoderiva como `atol_force_factor · escala_inicial`, donde la escala se calcula en el primer ensamblaje. Solo conviene ajustarlo si el análisis tiene tramos donde la carga característica colapsa transitoriamente y aparece falsa no-convergencia.
- `atol_disp_factor` (default `1.0e-9`): análogo para desplazamiento.

**Por qué los `atol` se autoderivan.** El término absoluto vive en las unidades del problema (newtons para fuerza, metros para desplazamiento). Codificarlo como constante global obligaría a ajustarlo al cambiar de unidades (N/m, kN/mm, MPa/mm…). En su lugar, Solidum calcula la escala característica del problema en el primer ensamblaje y multiplica por el factor adimensional. El resultado: **el mismo análisis planteado en distintos sistemas de unidades converge en el mismo número de iteraciones hasta paridad de bits**, sin tocar las tolerancias. La justificación arquitectural está en el ADR 0007.

**Cuándo afinar**. Para la inmensa mayoría de análisis los defaults son suficientes. Considera ajustar:

- `rtol_force` y `rtol_disp` más estrictos cuando publiques resultados o compares contra solución analítica (`1.0e-7` a `1.0e-9`).
- `rtol_disp` más laxo que `rtol_force` en problemas casi-rígidos donde el incremento de desplazamiento es naturalmente pequeño y baja muy rápido, dejando que el criterio de fuerza domine.
- Los factores `atol_*` prácticamente nunca; solo si encuentras un análisis específico donde la escala colapsa y los logs muestran iteraciones oscilando cerca del umbral.

**Nota sobre el solver algebraico interno.** El sistema lineal $\mathbf{K}\,\delta\mathbf{U}=\mathbf{R}$ que aparece dentro de cada iteración se resuelve con un backend algebraico (Cholesky, LU…) seleccionado *automáticamente* según las propiedades de $\mathbf{K}$. Esta elección **no** es decisión de modelado y no requiere configuración. Si por motivos de diagnóstico necesitas forzar un backend específico, consulta el capítulo *Anexos técnicos — Capa algebraica* del Manual de Referencia.

## Salida: `output`

Configura qué se exporta y cuándo:

```yaml
output:
  file_name: "resultados_placa"
  pre_process: {export: true}     # exporta el modelo sin deformar (paso 0)
  results:
    frequency: "all_steps"        # "last_step" o "all_steps"
    nodal_results: ["Displacements"]
    element_results: ["Von_Mises", "Internal_State", "Sigma_XX"]
  text_export:                     # opcional: salida estilo FEAP en .txt
    export: true
    nodes: all
    elements: all
```

```callout Nota
Cuando `frequency: "all_steps"` y el solver es no lineal, se genera un archivo `.vtu` por paso convergido más un archivo maestro `.pvd` que ParaView interpreta como animación.
```
