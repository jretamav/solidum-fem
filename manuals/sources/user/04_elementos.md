# Catálogo de Elementos Finitos

El motor expone dos familias de elementos: **1D** (estructuras reticulares: armaduras, cables, marcos) y **2D** (continuo plano: cuadrilátero y triángulo). Todos se referencian desde el YAML por su nombre exacto en `type:`.

## Elementos 1D

### Armaduras: `Truss2D` / `Truss2DCorot` / `Truss3D` / `Truss3DCorot`

Barra biarticulada que transmite exclusivamente fuerza axial. Cuatro variantes según dimensión y régimen geométrico:

[TABLA: Familia de armaduras.]
| Elemento | DOFs/nodo | Dimensión | Régimen geom. | Uso |
|---|---|---|---|---|
| `Truss2D` | `ux,uy` | 2D | lineal | cargas pequeñas, rotaciones pequeñas |
| `Truss2DCorot` | `ux,uy` | 2D | corotacional | grandes desplazamientos y rotaciones planas |
| `Truss3D` | `ux,uy,uz` | 3D | lineal | armaduras espaciales en régimen lineal |
| `Truss3DCorot` | `ux,uy,uz` | 3D | corotacional | grandes rotaciones en el espacio |

**Parámetros**: `A` (área de la sección transversal). Convención de signos: $\sigma > 0 \Leftrightarrow \varepsilon > 0$ (elongación).

**Régimen de validez**: $|\varepsilon| \lesssim 10^{-2}$. Las variantes `Corot` aceptan rotaciones de cualquier magnitud entre commits; las versiones lineales asumen rotaciones pequeñas.

```yaml
elements:
  - id: 1
    type: Truss2DCorot
    material: 1
    nodes: [1, 2]
    A: 0.0005
```

### Cables: `Cable2DCorot` / `Cable3DCorot`

Elemento 1D que transmite *únicamente* fuerza axial de tensión (no resiste compresión). La unilateralidad la aporta el material (típicamente `CableMaterial1D`, ver capítulo *Catálogo de Modelos Constitutivos*); el elemento implementa la cinemática corotacional y la transferencia de la deformación al material.

- **DOFs/nodo**: `ux, uy` (2D); `ux, uy, uz` (3D).
- **Parámetros**: `A`.
- **Cinemática**: corotacional (longitud y cosenos directores se recalculan en cada evaluación).
- **Régimen de validez**: $|\varepsilon| \lesssim 10^{-2}$ tensado. En estado destensado el elemento aporta rigidez nula al sistema global; otros elementos deben garantizar la estabilidad numérica.

```callout Advertencia
Un cable completamente destensado ($\sigma = 0$) tiene $\mathbf K_T = \mathbf 0$ en sus DOFs. Si todos los caminos de carga del modelo dependen del cable y este se destensa, la matriz tangente global se vuelve singular. Pretensar el cable o garantizar redundancia estructural.
```

```yaml
materials:
  - id: 1
    type: CableMaterial1D
    E: 150.0e9
    density: 7850.0

elements:
  - id: 1
    type: Cable2DCorot
    material: 1
    nodes: [1, 2]
    A: 5.0e-5
```

### Marcos 2D: `Frame2DEuler` / `Frame2DTimoshenko` / `Frame2DEulerCorot`

Elementos viga 2D que transmiten axial, cortante y momento flector. Tres variantes:

- **`Frame2DEuler`**: vigas esbeltas ($L/h \gtrsim 10$), Euler-Bernoulli, régimen geométricamente lineal. Parámetros: `A`, `I`.
- **`Frame2DTimoshenko`**: vigas peraltadas o cortas ($L/h \lesssim 10$). Incluye deformación por cortante; corrige automáticamente *shear locking* en el límite esbelto. Parámetros: `A`, `I`, `As` (área efectiva de cortante), `nu` (opcional si el material expone `nu`).
- **`Frame2DEulerCorot`**: variante corotacional de Euler-Bernoulli. Captura grandes desplazamientos y grandes rotaciones rígidas del elemento (rotaciones deformacionales nodales moderadas, $|\bar\theta| \lesssim 30°$). Parámetros: `A`, `I`.

DOFs/nodo: `ux, uy, rz`. Convención: $\sigma > 0 \Leftrightarrow \varepsilon > 0$ (elongación) en el eje axial; $r_z$ positivo antihorario.

```callout Advertencia
**Plasticidad por flexión.** Los elementos marco pasan al material únicamente la deformación axial centroidal $(u_j - u_i)/L$; el módulo tangente $E_t$ devuelto escala toda la matriz local (axial y flexión por igual). Esto significa que los marcos reproducen correctamente *plasticidad axial pura* pero *no* forman rótulas plásticas por flexión: la fluencia por momento requiere integración de $\sigma(y)$ sobre la sección, característica que se incorporará en una clase `FiberSection` futura.
```

```yaml
elements:
  - id: 1
    type: Frame2DEulerCorot
    material: 1
    nodes: [1, 2]
    A: 0.01
    I: 8.33e-6

  - id: 2
    type: Frame2DTimoshenko
    material: 1
    nodes: [2, 3]
    A: 0.01
    I: 8.33e-6
    As: 0.00833
```

### Marco 3D: `Frame3D`

Viga 3D Euler-Bernoulli con 6 DOFs/nodo (`ux, uy, uz, rx, ry, rz`). Transmite axial, cortantes en dos planos, flexiones en dos planos y torsión de Saint-Venant pura. Régimen geométricamente lineal (la variante corotacional 3D queda como pendiente).

**Parámetros**:

- `A` — área de la sección.
- `Iy`, `Iz` — momentos de inercia respecto a los ejes locales $y$, $z$.
- `J` — constante torsional de Saint-Venant.
- `nu` — coeficiente de Poisson (opcional; se toma del material si lo expone).
- `ref_vector` — vector de referencia (opcional, default `[0,0,1]`) que fija la orientación de los ejes locales $y$, $z$ de la sección. Para barras casi-verticales se sugiere fijarlo explícitamente.

```yaml
elements:
  - id: 1
    type: Frame3D
    material: 1
    nodes: [1, 2]
    A: 0.01
    Iy: 8.33e-6
    Iz: 8.33e-6
    J: 1.66e-5
    ref_vector: [0.0, 0.0, 1.0]
```

## Elementos 2D

La importación de mallas desde Gmsh instancia automáticamente los elementos 2D según la topología detectada (cuadriláteros → `Quad4`; triángulos → `Tri3`). También pueden declararse manualmente desde el bloque `elements`.

### Cuadrilátero Bilineal: `Quad4`

Elemento isoparamétrico de 4 nodos, integración Gauss configurable.

- **DOFs/nodo**: `ux, uy` ($\texttt{STRAIN\_DIM} = 3$, Voigt $[\varepsilon_{xx}, \varepsilon_{yy}, \gamma_{xy}]$).
- **Nodos**: 4, ordenados en sentido antihorario.
- **Cuadratura**: `"2x2"` por defecto (4 puntos). `"1x1"` disponible pero produce modos espurios (*hourglassing*).
- **Parámetros**: `thickness`, `quadrature` (opcional).
- **Implementación**: kernels críticos compilados con `@njit` (Numba).
- **Limitación**: bloqueo volumétrico con materiales casi-incompresibles ($\nu \to 0.5$); en ese régimen requeriría formulación mixta (no implementada).

### Triángulo de Deformación Constante: `Tri3`

Triángulo CST de 3 nodos, un único punto de integración.

- **DOFs/nodo**: `ux, uy` ($\texttt{STRAIN\_DIM} = 3$).
- **Cuadratura**: 1 punto (deformación uniforme).
- **Parámetros**: `thickness`.
- **Limitación**: *shear locking* severo, convergencia lenta. **Preferir `Quad4`** salvo en transiciones donde Quad4 no encaja geométricamente.

### Cuadrilátero Serendípito de Orden 2: `Quad8`

Cuadrilátero isoparamétrico cuadrático de 8 nodos (4 vértices antihorarios + 4 nodos medios de borde). Funciones de forma serendípitas (sin el término $\xi^2\eta^2$). Reproduce campos cuadráticos exactamente y mitiga el bloqueo por cortante en problemas de flexión.

- **DOFs/nodo**: `ux, uy` ($\texttt{STRAIN\_DIM} = 3$).
- **Cuadratura**: Gauss 3×3 (9 puntos) por defecto; 2×2 disponible vía `quadrature` pero con riesgo de modos espurios.
- **Parámetros**: `thickness`, `quadrature` (opcional).
- **Tracción de borde**: reparte $1/6$, $4/6$, $1/6$ entre los nodos del borde (vértice, medio, vértice).

### Cuadrilátero Lagrangiano de Orden 2: `Quad9`

Cuadrilátero Lagrangiano cuadrático de 9 nodos (los 8 de `Quad8` más un noveno nodo central interior). Las funciones de forma son producto tensorial Lagrange 1D-1D, espacio polinómico completo $Q_2$.

- **DOFs/nodo**: `ux, uy` ($\texttt{STRAIN\_DIM} = 3$).
- **Cuadratura**: Gauss 3×3.
- **Parámetros**: `thickness`, `quadrature` (opcional).
- **Diferencia con `Quad8`**: el nodo central interior añade el término $\xi^2\eta^2$ y mejora el comportamiento en problemas con campos no separables.

### Triángulo Cuadrático Completo P₂: `Tri6`

Triángulo isoparamétrico cuadrático de 6 nodos (3 vértices + 3 medios de borde). Resuelve el *shear locking* severo de `Tri3` y reproduce campos cuadráticos exactamente.

- **DOFs/nodo**: `ux, uy` ($\texttt{STRAIN\_DIM} = 3$).
- **Cuadratura**: 3 puntos en los puntos medios (regla `tri_3`).
- **Parámetros**: `thickness`.
- **Tracción de borde**: reparte $1/6$, $4/6$, $1/6$ entre los nodos del borde (vértice, medio, vértice).
- **Cuándo usarlo**: transiciones de malla cuadráticas, geometrías curvadas donde Quad8/Quad9 no encajan.

```yaml
elements:
  - {id: 1, type: Quad8, material: 1, thickness: 0.1, nodes: [1, 2, 3, 4, 5, 6, 7, 8]}
  - {id: 2, type: Quad9, material: 1, thickness: 0.1, nodes: [1, 2, 3, 4, 5, 6, 7, 8, 9]}
  - {id: 3, type: Tri6,  material: 1, thickness: 0.1, nodes: [1, 2, 3, 4, 5, 6]}
```
