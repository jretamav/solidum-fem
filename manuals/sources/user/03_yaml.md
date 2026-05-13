# Archivo de Entrada YAML

La interacción exclusiva con el motor de Fenix FEM se realiza mediante el archivo de entrada `.yaml`. Un archivo válido contiene los bloques siguientes (algunos opcionales según el caso de uso):

- `nodes` *o* `mesh` — geometría.
- `materials` — catálogo de leyes constitutivas instanciadas.
- `elements` — conectividad (omisible si se importa malla Gmsh).
- `boundary_conditions_by_node` / `_by_coord` / `_by_group` — restricciones de Dirichlet.
- `point_loads_by_node` / `_by_coord` / `_by_group` — cargas de Neumann puntuales.
- `solver` — tipo de solver y sus parámetros.
- `output` — frecuencia y formato de exportación de resultados.

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
| `mesh_quadrature` | Regla de Gauss: `"2x2"` (completa, recomendada) o `"1x1"` (reducida; produce *hourglassing*). |

### Mapeo Avanzado: `mesh_physical_groups`

Si la geometría tiene zonas con materiales o espesores distintos definidos como *Physical Groups* en Gmsh, se mapean explícitamente:

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

## Condiciones de Frontera

Fenix FEM ofrece tres mecanismos para imponer restricciones cinemáticas (Dirichlet), elegibles según la facilidad operativa:

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

Las claves predefinidas `x_min`, `x_max`, `y_min`, `y_max` usan los extremos detectados automáticamente del dominio.

### Por grupo físico: `boundary_conditions_by_group`

Aplica restricciones a todos los nodos de un *Physical Group* de Gmsh:

```yaml
boundary_conditions_by_group:
  "Borde_Empotrado": {ux: 0.0, uy: 0.0}
```

### Restricciones multipunto (MPC): `linear_constraints`

Para apoyos en plano oblicuo, periodicidad de celda unitaria, uniones rígidas master-slave y simetrías no alineadas con los ejes globales, Fenix FEM admite restricciones afines lineales de la forma $u_s = g_s + \sum_i \alpha_{si} u_{m_i}$ (ADR 0004 fase 2). Se declaran en el bloque `linear_constraints`:

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

Las restricciones se imponen por eliminación directa (no por penalización), de modo que la imposición es exacta a redondeo y preserva la simetría y la positividad definida del sistema. Las cadenas master-slave se resuelven automáticamente por cierre transitivo; los ciclos y las redeclaraciones inconsistentes se detectan en construcción y abortan el análisis con un mensaje específico.

## Cargas Puntuales

Misma estructura que las condiciones de frontera, en bloques `point_loads_by_node`, `point_loads_by_coord` y `point_loads_by_group`. Los valores son fuerzas nodales (en las unidades del modelo, típicamente Newtons):

```yaml
point_loads_by_node:
  - {node_id: 3, uy: -15000.0}

point_loads_by_group:
  "Borde_Carga": {ux: 1500000.0}
```

## Fuerzas de Cuerpo

Una *fuerza de cuerpo* es una fuerza por unidad de volumen $\mathbf b$ que actúa sobre todo el dominio material del elemento, no solo sobre su frontera. Cada elemento del catálogo (armaduras, cables, marcos 2D/3D, sólidos 2D) integra la contribución consistente $\int \mathbf{N}^\top \mathbf{b}\,A\,dx$ y la acumula al vector global de fuerzas. Para vigas Hermite cúbicas el reparto incluye fuerza nodal $qL/2$ y momento $\pm qL^2/12$; para barras axiales es mitad-mitad por nodo.

Casos típicos en mecánica estructural: peso propio (un caso particular con $\mathbf b = \rho\,\mathbf g$), fuerzas centrífugas en cuerpos rotatorios ($\mathbf b = \rho\,\omega^2\,\mathbf r$), fuerzas electromagnéticas sobre materiales conductores, expansión térmica modelada como pseudocarga. Fenix expone dos formas declarativas en YAML que cubren todos estos casos.

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

### Caso particular — peso propio: `gravity` + `density`

El peso propio es el caso especial $\mathbf b = \rho\,\mathbf g$. Para problemas multimaterial — donde cada elemento debe ver *su propio* $\rho$ — Fenix expone un atajo declarativo: la densidad como propiedad del material, y un vector global de gravedad.

```yaml
materials:
  - {id: 1, type: Elastic1D, E: 210.0e9, density: 7850.0}     # acero
  - {id: 2, type: Elastic1D, E:  30.0e9, density: 2500.0}     # hormigón

gravity: [0.0, -9.81]                # m/s^2, +y hacia arriba
```

**Qué hace Fenix**. Para cada elemento computa $\mathbf b_e = \rho_e\cdot\mathbf g$, donde $\rho_e$ es la densidad del material asignado, e integra la contribución consistente con la misma maquinaria de `body_force`. Resultado: estructuras multimaterial dan peso propio físicamente correcto, cada elemento con su propio $\rho$.

**Densidad solo cuando se usa**. El campo `density` es **opcional** al declarar un material (ADR 0008): los análisis estáticos puros que no invocan peso propio ni matriz de masa no necesitan declararla. Pero al invocar `gravity` (o, en su día, análisis modal/dinámico), Fenix exige la densidad de cada material involucrado: si algún material no la trae declarada, el cálculo falla con `ValueError` listando los materiales afectados por nombre. Imposible propagar masa cero silenciosa. Si un material legítimamente carece de masa física (penalty, restricción), declararlo explícitamente como `density: 0.0` — caso en el cual `gravity` emite un `WARNING` informativo y su contribución al peso propio es nula.

### Detalles transversales

**Exclusividad**. `body_force` y `gravity` son dos formas declarativas para expresar fuerzas de cuerpo; declarar ambas en el mismo archivo lanza error de validación. La forma general (`body_force`) es la subyacente; `gravity` es atajo para su caso particular más común.

**Padding 2D↔3D**. Si el modelo mezcla elementos 2D y 3D, declarar el vector con tres componentes; el ensamblador pasa el slice apropiado a cada elemento según la dimensionalidad declarada por su `DOF_NAMES`.

**Reservado para análisis dinámico**. El atributo `density` introducido aquí es la misma magnitud que consume la matriz de masa $\mathbf M_e = \int \rho\,\mathbf N^\top\mathbf N\,dV$ en análisis modales y dinámicos. Declararlo prepara el modelo para esas extensiones sin coste adicional (ADR 0008).

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

### Bloque `convergence`

Los solvers no lineales (`NonlinearSolver`, `ArcLengthSolver`) aceptan un bloque `convergence` que configura el criterio de parada de las iteraciones de Newton. La política sigue el patrón estándar de tolerancias mixtas absoluta + relativa, separado por dos criterios físicos: *fuerza* (residuo de equilibrio) y *desplazamiento* (incremento del iterado). La forma exacta de la comparación es:

```
||R[libres]|| <= atol_force + rtol_force * max(||F_ext||, ||F_int||)
||dU||       <= atol_disp  + rtol_disp  * ||U||
```

Ambos criterios deben cumplirse **simultáneamente** (semántica AND) para dar el paso por convergido.

**Parámetros disponibles** (todos opcionales, con defaults razonables):

- `rtol_force` (default `1.0e-5`): tolerancia relativa del residuo de fuerza. Es el parámetro que normalmente ajustarás: bajarlo a `1.0e-7` para análisis de precisión, subirlo a `1.0e-3` para exploración rápida.
- `rtol_disp` (default `1.0e-5`): tolerancia relativa de la corrección de desplazamiento. En problemas casi-rígidos puede afinarse independientemente del de fuerza.
- `atol_force_factor` (default `1.0e-9`): factor adimensional que define el piso absoluto de fuerza. La tolerancia absoluta efectiva se autoderiva como `atol_force_factor · escala_inicial`, donde la escala se calcula en el primer ensamblaje. Solo conviene ajustarlo si el análisis tiene tramos donde la carga característica colapsa transitoriamente y aparece falsa no-convergencia.
- `atol_disp_factor` (default `1.0e-9`): análogo para desplazamiento.

**Por qué los `atol` se autoderivan.** El término absoluto vive en las unidades del problema (newtons para fuerza, metros para desplazamiento). Codificarlo como constante global obligaría a ajustarlo al cambiar de unidades (N/m, kN/mm, MPa/mm…). En su lugar, Fenix calcula la escala característica del problema en el primer ensamblaje y multiplica por el factor adimensional. El resultado: **el mismo análisis planteado en distintos sistemas de unidades converge en el mismo número de iteraciones hasta paridad de bits**, sin tocar las tolerancias. La justificación arquitectural está en el ADR 0007.

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
