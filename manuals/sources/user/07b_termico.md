# Análisis Térmico

Solidum FEM resuelve **conducción de calor** por la ley de Fourier, en régimen estacionario y transitorio, sobre dominios 2D y 3D. El campo incógnita es un escalar $T$ por nodo.

**Alcance actual**: conducción pura con condiciones de Dirichlet (temperatura impuesta) y Neumann (flujo prescrito). **No** incluye convección (Robin), radiación, conductividad dependiente de la temperatura $k(T)$, cambio de fase ni **acoplamiento termomecánico** — no existe deformación térmica $\alpha\,\Delta T$: el análisis térmico y el mecánico son hoy independientes.

## Formulación

La ecuación de balance de energía, con $\mathbf q = -\mathbf k\cdot\nabla T$ (Fourier):

$$\rho c\,\dot T = \nabla\cdot(\mathbf k\,\nabla T) + Q$$

Discretizada por elementos finitos:

$$\mathbf C\,\dot{\mathbf T} + \mathbf K\,\mathbf T = \mathbf F$$

donde $\mathbf K = \int_\Omega \mathbf B^\top\mathbf k\,\mathbf B\,d\Omega$ es la matriz de **conductividad**, $\mathbf C = \int_\Omega \rho c\,\mathbf N^\top\mathbf N\,d\Omega$ la de **capacidad calorífica** y $\mathbf F$ recoge la fuente volumétrica y el flujo de frontera.

Nótese que el sistema es de **primer orden** en el tiempo, a diferencia del mecánico $\mathbf M\ddot{\mathbf u} + \mathbf C\dot{\mathbf u} + \mathbf K\mathbf u = \mathbf F$. No hay aceleración ni inercia: no existe un análogo térmico de la velocidad con significado físico en este modelo.

**Sin notación Voigt.** El gradiente $\nabla T$ es un vector genuino de 2 ó 3 componentes, no un tensor simétrico comprimido, así que la matriz $\mathbf B$ térmica es el gradiente crudo $\partial N/\partial x$ y no lleva los factores $\gamma = 2\varepsilon$ de la notación Voigt mecánica. Los elementos térmicos declaran `FLUX_DIM` (2 ó 3) en lugar de `STRAIN_DIM`.

## Material: `ThermalConduction`

Se declara en un bloque **`thermal_materials`**, paralelo a `materials`. Los dos bloques tienen espacios de nombres de `type:` independientes, y un elemento resuelve su material contra el bloque donde se declaró el `id`.

[TABLA: Parámetros de `ThermalConduction`.]
| Parámetro | Unidades | Obligatorio | Descripción |
|---|---|---|---|
| `k` | W/(m·K) | sí | Conductividad. Escalar (isótropo) o tensor completo |
| `dim` | — | no (2) | Dimensión del problema cuando `k` es escalar |
| `c` | J/(kg·K) | sólo transitorio | Calor específico |
| `density` | kg/m³ | sólo transitorio | Densidad |

```yaml
thermal_materials:
  - id: 1
    type: ThermalConduction
    k: 45.0            # acero estructural
    c: 460.0           # sólo hace falta en transitorio
    density: 7850.0    # sólo hace falta en transitorio
```

**`c` y `density` son opcionales**: el régimen estacionario no los necesita, porque el término $\mathbf C\dot{\mathbf T}$ desaparece. Si faltan y se pide un transitorio, el error lo dice explícitamente y sugiere el solver estacionario.

**Conductividad anisótropa**: `k` admite un tensor completo, sin necesidad de un material distinto. Debe ser simétrico y definido positivo — ambas cosas se validan al construir.

```yaml
thermal_materials:
  - id: 1
    type: ThermalConduction
    k: [[45.0, 0.0], [0.0, 12.0]]   # 3.75× más conductor en x que en y
```

## Elementos: `Quad4Thermal` y `Hex8Thermal`

Un DOF escalar `T` por nodo. Comparten geometría y numeración de nodos con sus gemelos mecánicos `Quad4` y `Hex8`.

[TABLA: Elementos térmicos.]
| Elemento | `FLUX_DIM` | Nodos | Cuadratura default | Parámetros |
|---|---|---|---|---|
| `Quad4Thermal` | 2 | 4 | `2x2` | `thickness` (default 1.0), `quadrature` |
| `Hex8Thermal` | 3 | 8 | `hex_2x2x2` | `quadrature` |

```yaml
elements:
  - id: 1
    type: Quad4Thermal
    material: 1          # id del bloque thermal_materials
    nodes: [1, 2, 3, 4]
    thickness: 0.3
```

Pasar un material mecánico a un elemento térmico —o al revés— se rechaza al construir, con un mensaje que nombra las dos familias. La dimensión también se comprueba: un material construido con `dim: 2` no entra en un `Hex8Thermal`.

## Condiciones de frontera

### Dirichlet — temperatura impuesta

Se declaran como cualquier otra condición de frontera, usando el nombre del DOF:

```yaml
boundary_conditions:
  - {node_id: 1, T: 100.0}
  - {node_id: 2, T: 20.0}
```

**Todo modelo térmico necesita al menos un Dirichlet.** Sin ninguno, la matriz de conductividad es singular: un campo de temperatura uniforme no produce gradiente ni flujo, y el nivel absoluto de temperatura queda indeterminado. Es el análogo térmico de un modo de sólido rígido, y el error lo nombra así en vez de reportar un fallo algebraico genérico.

### Neumann — flujo prescrito

```yaml
thermal_loads:
  boundary_flux:
    - {element: 1, edge: 2, q: 500.0}      # 2D: borde del elemento
    - {element: 7, face: 4, q: -200.0}     # 3D: cara del elemento
```

**Convención de signo**: $\bar q > 0$ es flujo **saliente** del dominio, es decir **enfriamiento**. Un flujo entrante (calentamiento) se declara con $\bar q < 0$. El signo sale de la forma débil, donde el término de frontera aparece como $-\int_\Gamma w\,\bar q\,d\Gamma$.

**Un borde o cara sin declarar es adiabático** ($\bar q = 0$), igual que un borde sin carga es libre de tracción en mecánica. No hay que declarar el flujo nulo.

Los índices de borde y cara son los mismos que en los elementos mecánicos: `edge: 0` conecta los nodos locales 0-1, `edge: 1` los 1-2, etc.

### Fuente volumétrica

```yaml
thermal_loads:
  body_source:
    - {Q: 1.0e5}                       # W/m³, a todo el dominio
    - {Q: 5.0e4, elements: [3, 4, 5]}  # sólo a esos elementos
```

`Q` es generación de calor por unidad de volumen (calor de hidratación, efecto Joule, reacción química). Es **carga del elemento y no propiedad del material**: varía en el espacio y con el tiempo, y como propiedad obligaría a duplicar materiales.

También pueden aplicarse fuentes **nodales concentradas** [W] por la vía estándar:

```yaml
point_loads:
  - {node_id: 12, T: 250.0}
```

## Régimen estacionario

$\mathbf K\,\mathbf T = \mathbf F$ lo resuelve el **`LinearSolver` sin ninguna configuración especial**:

```yaml
solver:
  type: LinearSolver
```

El resultado es un `SolveResult` como el de un análisis mecánico; el campo `U` contiene las temperaturas nodales, y `reactions_by_node` los flujos de calor en los nodos con temperatura impuesta.

## Régimen transitorio: `ThetaMethodSolver`

Integra $\mathbf C\dot{\mathbf T} + \mathbf K\mathbf T = \mathbf F(t)$ ponderando el estado dentro del paso con el parámetro $\theta$:

$$\left(\mathbf C + \theta\Delta t\,\mathbf K\right)\mathbf T_{n+1} = \left[\mathbf C - (1-\theta)\Delta t\,\mathbf K\right]\mathbf T_n + \Delta t\left[\theta\mathbf F_{n+1} + (1-\theta)\mathbf F_n\right]$$

```yaml
solver:
  type: ThetaMethodSolver
  dt: 60.0
  n_steps: 500
  T_initial: 20.0
```

[TABLA: Parámetros de `ThetaMethodSolver`.]
| Parámetro | Default | Descripción |
|---|---|---|
| `dt` | — | Paso temporal [s]. Obligatorio |
| `n_steps` | — | Número de pasos. Obligatorio |
| `T_initial` | — | Condición inicial: escalar (campo uniforme) o vector. Obligatorio |
| `theta` | `1.0` | Peso del esquema $\in [0,1]$ |
| `lumping` | `"lumped"` | Forma de la matriz de capacidad |
| `output_every` | `1` | Almacenar el campo cada N pasos |

**No hay condición inicial de velocidad**: la ecuación es de primer orden y $\mathbf T_0$ determina la evolución por completo. Tampoco existe el parámetro `rayleigh`: el amortiguamiento de Rayleigh modela disipación en la ecuación de segundo orden, y aquí la disipación es la propia conducción, ya contenida en $\mathbf K$.

### Elección de $\theta$

[TABLA: Casos particulares del esquema $\theta$.]
| $\theta$ | Nombre | Orden | Estabilidad |
|---|---|---|---|
| `0` | Euler explícito | 1 | Condicional ($\Delta t \le 2/\lambda_{\max}$) |
| `0.5` | Crank-Nicolson | **2** | Incondicional (A-estable) |
| `2/3` | Galerkin | 1 | Incondicional |
| `1` | Euler implícito | 1 | Incondicional, **L-estable** (default) |

**Estabilidad no significa ausencia de oscilaciones**, y esta distinción decide el valor por defecto.

Crank-Nicolson ($\theta = 1/2$) es de segundo orden y por tanto el más preciso, pero su factor de amplificación tiende a $-1$ para los modos de frecuencia alta en lugar de a $0$. Ante un **cambio brusco** —un escalón de temperatura impuesto en la frontera, típico en el arranque— produce oscilaciones amortiguadas lentamente, con temperaturas que pueden **salirse del rango de los datos**. Eso viola el *principio del máximo* de la ecuación de difusión, según el cual la solución exacta nunca excede los extremos de los datos iniciales y de frontera: no es un resultado impreciso, es un resultado imposible.

Euler implícito ($\theta = 1$) es sólo de primer orden, pero **L-estable**: amortigua los modos altos por completo en un paso y nunca oscila. Por eso es el default.

En un ensayo con pared a $100$ y cuerpo inicialmente a $0$, con paso grande, $\theta = 1$ se mantiene en $[0, 100]$ y es monótono, mientras $\theta = 1/2$ alcanza $151.27$.

Si su problema es suave y busca la precisión de segundo orden, use `theta: 0.5` conscientemente. El resultado expone el campo `order` (2 con Crank-Nicolson, 1 en el resto) para que el coste en precisión del default robusto sea visible.

### Elección de $\Delta t$

Con $\theta \ge 0.5$ **ningún paso diverge**, pero eso no significa que cualquier paso sirva. La escala física la da la difusividad $\alpha = k/(\rho c)$ y el tamaño de elemento $h$:

$$\Delta t_{\text{car}} \sim \frac{h^2}{\alpha}$$

Es el tiempo que tarda el frente térmico en atravesar un elemento. Un $\Delta t$ mucho mayor es estable pero **se salta el transitorio** que se quiere observar. El solver reporta este valor al arrancar y avisa si el paso elegido lo supera ampliamente. Es información, no restricción.

En el otro extremo, el transitorio no ha **terminado** hasta $t \gg L^2/\alpha$ con $L$ la dimensión global del dominio. Integrar por debajo de esa escala y comparar con el estacionario da diferencias que no son error del esquema, sino un análisis inacabado.

### Capacidad `lumped` por defecto

A diferencia del análisis dinámico estructural, donde el default es `consistent`, aquí la matriz de capacidad se agrupa (`lumped`) por defecto. El motivo es el mismo principio del máximo: la capacidad consistente produce oscilaciones espurias ante un frente térmico abrupto. Los dos defaults invertidos —`lumped` en el eje espacial y $\theta = 1$ en el temporal— atacan el mismo fenómeno desde lados distintos.

### Resultado

`ThetaMethodSolver` devuelve un **`ThermalTransientResult`**, no el `TransientResult` del análisis dinámico: no hay velocidad, aceleración ni coeficientes de Rayleigh que reportar.

[TABLA: Campos de `ThermalTransientResult`.]
| Campo / método | Descripción |
|---|---|
| `t_history` | Instantes almacenados, shape `(n_stored,)` |
| `T_history` | Campo global, shape `(n_dof, n_stored)` |
| `T_final` | Campo en el último instante |
| `temperature_at(dof)` | Historia temporal de un DOF concreto |
| `extremes()` | `(T_min, T_max)` sobre todos los nodos e instantes |
| `n_steps`, `theta`, `dt`, `order` | Parámetros efectivos del análisis |

`extremes()` es la herramienta de diagnóstico del principio del máximo: en un problema sin fuente, cualquier valor fuera del rango de los datos señala oscilación espuria del esquema, no un fenómeno físico.

## Ejemplo completo: pared plana en régimen transitorio

Pared de acero de 2 m, inicialmente a 20 °C, con una cara llevada súbitamente a 100 °C.

```yaml
nodes:
  - {id: 1, coords: [0.0, 0.0]}
  - {id: 2, coords: [1.0, 0.0]}
  - {id: 3, coords: [2.0, 0.0]}
  - {id: 4, coords: [0.0, 0.5]}
  - {id: 5, coords: [1.0, 0.5]}
  - {id: 6, coords: [2.0, 0.5]}

thermal_materials:
  - {id: 1, type: ThermalConduction, k: 45.0, c: 460.0, density: 7850.0}

elements:
  - {id: 1, type: Quad4Thermal, nodes: [1, 2, 5, 4], material: 1, thickness: 0.3}
  - {id: 2, type: Quad4Thermal, nodes: [2, 3, 6, 5], material: 1, thickness: 0.3}

boundary_conditions:
  - {node_id: 1, T: 100.0}
  - {node_id: 4, T: 100.0}
  - {node_id: 3, T: 20.0}
  - {node_id: 6, T: 20.0}

solver:
  type: ThetaMethodSolver
  dt: 2000.0
  n_steps: 200
  T_initial: 20.0
```

Con estos valores el nodo central llega a $59.998$, frente al $60$ exacto del perfil lineal estacionario. La diferencia **no es error del esquema**: con $\Delta t = 2000$ s y 200 pasos el análisis cubre $4\times10^5$ s, mientras el tiempo de difusión global de esta pared es $L^2/\alpha \approx 3.2\times10^5$ s — el transitorio está terminando pero aún no ha terminado. Es exactamente el efecto descrito en «Elección de $\Delta t$»: alargando el análisis, el campo converge al perfil lineal $100 \to 20$ que da el `LinearSolver` sobre el mismo modelo, a precisión de máquina.

Como la condición inicial (20 °C) contradice el Dirichlet impuesto (100 °C) en la cara izquierda, el solver **avisa** de la discontinuidad en $t=0$ sin abortar: un choque térmico es modelización legítima, y es precisamente el caso que motiva el default $\theta = 1$.

## Limitaciones actuales

- **Sin convección ni radiación**. Una frontera con convección $q = h(T - T_\infty)$ debe modelarse hoy como flujo prescrito equivalente, lo que exige conocer $T$ de antemano — es decir, no es un sustituto general.
- **Sin acoplamiento termomecánico**: el campo de temperatura no produce deformación. Un análisis termoelástico requiere hoy resolver el térmico y trasladar el resultado a mano.
- **Material lineal**: $k$ no depende de la temperatura, y no hay cambio de fase. Cualquiera de las dos cosas exigiría iteraciones de Newton dentro de cada paso temporal.
- **Paso de tiempo constante**: no hay adaptación automática.
- **Sólo elementos lineales**: `Quad4Thermal` y `Hex8Thermal`. No existen versiones cuadráticas ni triangulares/tetraédricas térmicas.
