# Catálogo de Elementos Finitos

El programa principal expone varias familias de elementos: **1D** (estructuras reticulares: armaduras, cables, marcos), **2D** (continuo plano: cuadrilátero y triángulo), **3D** (continuo tridimensional: hexaedros y tetraedros, lineales y cuadráticos) y **térmicos** (conducción de calor, con un grado de libertad escalar por nodo; se documentan en el capítulo «Análisis Térmico»). Todos se referencian desde el YAML por su nombre exacto en `type:`.

Las formulaciones no estándar viven fuera del programa principal, en **módulos de usuario** que el modelo carga de forma explícita (ADR 0020). El elemento con discontinuidad embebida (fractura con salto interior) pertenece al módulo de usuario `discontinuities` y se documenta en la última sección de este capítulo.

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

Elemento 1D que transmite *únicamente* fuerza axial de tracción (no resiste compresión). La unilateralidad la aporta el material (típicamente `CableMaterial1D`, ver capítulo *Catálogo de Modelos Constitutivos*); el elemento implementa la cinemática corotacional y la transferencia de la deformación al material.

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

- **`Frame2DEuler`**: vigas esbeltas ($L/h \gtrsim 10$), Euler-Bernoulli, régimen geométricamente lineal. Parámetros: `A`, `I` (con `material`) o `section`, `n_ip`.
- **`Frame2DTimoshenko`**: vigas peraltadas o cortas ($L/h \lesssim 10$). Incluye deformación por cortante; corrige automáticamente *shear locking* en el límite esbelto. Parámetros: `A`, `I`, `As` (área efectiva de cortante), `nu` (opcional si el material expone `nu`). **Sólo material elástico lineal**; todavía no admite `section`.
- **`Frame2DEulerCorot`**: variante corotacional de Euler-Bernoulli. Captura grandes desplazamientos y grandes rotaciones rígidas del elemento (rotaciones deformacionales nodales moderadas, $|\bar\theta| \lesssim 30°$). Parámetros: `A`, `I` (con `material`) o `section`, `n_ip`; con sección, ésta trabaja en el sistema corrotado del elemento.

DOFs/nodo: `ux, uy, rz`. Convención: $\sigma > 0 \Leftrightarrow \varepsilon > 0$ (elongación) en el eje axial; $r_z$ positivo antihorario.

**Dos formas de dar la constitutiva** (ADR 0021), excluyentes entre sí:

- **`material` + propiedades geométricas** (`A`, `I`): rigidez de forma cerrada. Sólo admite un material **elástico lineal** (`Elastic1D`); un material con plasticidad o daño se rechaza al construir el modelo, con un mensaje que remite a la sección por fibras.
- **`section`** (id de una sección del bloque `sections`, capítulo *Archivo de Entrada YAML*) **+ `n_ip`** (puntos de Gauss a lo largo del elemento, 3 por defecto, mínimo 2): el elemento integra la sección por fibras en cada punto de Gauss. Cada fibra lleva su material 1D y su historia, de modo que se forman **rótulas plásticas por flexión** y el daño avanza fibra a fibra. Con fibras elásticas y `n_ip ≥ 2` se obtiene exactamente la misma rigidez que con `A` e `I`. Con `section` no se declaran `material`, `A` ni `I`: rigidez, masa y peso propio salen de las fibras.

```callout Advertencia
**Cambio desde el 2026-09-29.** Hasta esa fecha los marcos aceptaban un material inelástico con `A` e `I`, le pasaban sólo la deformación axial y escalaban toda la matriz local con su módulo tangente: con fluencia axial el momento flector se anulaba y con ablandamiento cambiaba de signo, sin aviso. Un modelo así ahora falla al construirse; su equivalente físico es el mismo marco con `section`. Con la sección por fibras, la carga de colapso converge por arriba al refinar la malla: use varios elementos por miembro cerca de las rótulas.
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

  - id: 3                     # con sección por fibras (bloque sections)
    type: Frame2DEuler
    nodes: [3, 4]
    section: 1
    n_ip: 3
```

### Marco 3D: `Frame3D`

Viga 3D Euler-Bernoulli con 6 DOFs/nodo (`ux, uy, uz, rx, ry, rz`). Transmite axial, cortantes en dos planos, flexiones en dos planos y torsión de Saint-Venant pura. Régimen geométricamente lineal (la variante corotacional 3D queda como pendiente).

**Parámetros**, con las mismas dos formas excluyentes de dar la constitutiva que los marcos 2D:

- Con `material` (sólo elástico lineal): `A` — área de la sección; `Iy`, `Iz` — momentos de inercia respecto a los ejes locales $y$, $z$; `J` — constante torsional de Saint-Venant; `nu` — coeficiente de Poisson (opcional; se toma del material si lo expone).
- Con `section` (sección por fibras con `GJ` declarado) y `n_ip` (default 3, mínimo 2): flexión en los dos planos y fuerza axial salen de las fibras, con su interacción; la torsión es elástica, $T = GJ\,\theta_x'$, con el `GJ` de la sección. Sin `nu`.
- `ref_vector` — vector de referencia (opcional, default `[0,0,1]`) que fija la orientación de los ejes locales $y$, $z$ de la sección, que son también los ejes en que se declaran las fibras. Para barras casi-verticales se sugiere fijarlo explícitamente.

Convención de fuerzas internas (Reglas.md §5): $M_y > 0$ tracciona las fibras con $z > 0$; $M_z > 0$, las fibras con $y < 0$.

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

  - {id: 2, type: Frame3D, nodes: [2, 3], section: 1, n_ip: 3, ref_vector: [0.0, 0.0, 1.0]}
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

## Elementos 3D

Sólidos tridimensionales isoparamétricos (ADR 0012). Todos comparten DOFs `ux, uy, uz` y $\texttt{STRAIN\_DIM} = 6$ sobre la convención Voigt 3D del proyecto, $[\varepsilon_{xx}, \varepsilon_{yy}, \varepsilon_{zz}, \gamma_{xy}, \gamma_{yz}, \gamma_{xz}]$, y exigen un material 3D (`Elastic3D`, `VonMises3D`, `DruckerPrager3D`, `IsotropicDamage3D`).

A diferencia de los 2D, **no llevan `thickness`**: el volumen sale de la geometría. Las cargas de superficie se aplican por **cara numerada con normal saliente** mediante `compute_face_traction(face, t̄)`, con `t̄` expresada en ejes globales.

### Hexaedro Trilineal: `Hex8`

Hexaedro isoparamétrico de 8 nodos, espejo natural del `Quad4`. Orden de nodos VTK_HEXAHEDRON.

- **Nodos**: 8 · **Cuadratura**: Gauss 2×2×2 (8 puntos) por defecto.
- **Parámetros**: `quadrature` (opcional; `hex_3x3x3` para no lineales severos, `hex_1x1x1` reducida con riesgo de *hourglass*).
- **Caras**: 6, numeradas 0 (−ζ), 1 (+ζ), 2 (−η), 3 (+ξ), 4 (+η), 5 (−ξ).
- **Limitaciones**: bloqueo volumétrico con $\nu \to 0.5$ y *shear locking* con una sola capa de elementos en la sección. Sin mitigación, por política idéntica a la del `Quad4`.

### Tetraedro Lineal: `Tet4`

Tetraedro de 4 nodos, CST 3D — espejo del `Tri3`. Deformación uniforme en el elemento.

- **Nodos**: 4 · **Cuadratura**: 1 punto baricéntrico.
- **Caras**: 4 triangulares; la cara $i$ es la opuesta al nodo $i$.
- **Limitación**: *shear locking* severo, peor que el `Hex8` por la pobreza del espacio lineal, y bloqueo volumétrico aún más acusado. **Preferir `Hex8`** en mallas hexaédricas; reservar `Tet4` para transiciones o geometrías que no admitan hexaedros.

### Hexaedro Serendípito de Orden 2: `Hex20`

Hexaedro cuadrático de 20 nodos (8 vértices + 12 medios de arista, sin nodos de cara ni centroide). Análogo 3D del `Quad8`.

- **Nodos**: 20 · **Cuadratura**: Gauss 3×3×3 (27 puntos) por defecto.
- **Caras**: 6, de 8 nodos cada una (cara `Quad8`).
- **Tracción de cara**: reparto serendípito — cada vértice recibe $-A\bar{t}/12$ y cada medio $4A\bar{t}/12$. Los signos negativos en los vértices son el fenómeno serendípito conocido del `Quad8` (Cook-Malkus-Plesha-Witt §6.5), no un error.
- **Ventaja medida**: en la viga esbelta de MacNeal-Harder, una malla 6×1×1 de `Hex20` alcanza el 97 % de la deflexión de Euler-Bernoulli, frente a menos del 55 % de un `Hex8` 12×1×1.
- **Precaución**: con `hex_2x2x2` aparecen 6 modos de *hourglass* por elemento aislado, sin estabilización. Usar la cuadratura por defecto.

### Hexaedro Lagrangiano Triquadrático: `Hex27`

Hexaedro Lagrangiano completo de 27 nodos (20 del `Hex20` + 6 centros de cara + 1 centro de cuerpo). Análogo 3D del `Quad9`.

- **Nodos**: 27 · **Cuadratura**: Gauss 3×3×3.
- **Caras**: 6, de 9 nodos cada una (cara `Quad9`).
- **Diferencia con `Hex20`**: reproduce **todos** los polinomios triquadráticos, incluidos los términos $\xi^2\eta^2\zeta^2$ que faltan en el serendípito. En flexión simple con geometría rectilínea ambos dan resultados prácticamente idénticos, con mayor coste para el `Hex27`; la diferencia se paga en geometrías de curvatura severa (Bathe FEP §5.3.2).
- **Precaución**: con `hex_2x2x2` aparecen 27 modos de *hourglass* por elemento aislado — la combinación más problemática del catálogo 3D. Usar 3×3×3 siempre.

### Tetraedro Cuadrático: `Tet10`

Tetraedro de 10 nodos (4 vértices + 6 medios de arista). Análogo 3D del `Tri6`.

- **Nodos**: 10 · **Cuadratura**: Stroud `tet_4` (4 puntos) por defecto; `tet_15` (Keast, 15 puntos) para elementos distorsionados.
- **Caras**: 4 triangulares de 6 nodos (cara `Tri6`).
- **Masa**: la matriz consistente usa `tet_15` fijo, independientemente de la cuadratura elegida para $K$, para integrar exactamente el producto cuadrático×cuadrático en análisis modal y transitorio.
- **Cuándo usarlo**: mallas no estructuradas y geometrías complejas donde un mapeo hexaédrico no es viable. Mitiga drásticamente ambos bloqueos del `Tet4`.
- **Limitación conocida**: sobre **superficies curvas** malladas por descomposición de hexaedros en tetraedros, la representación isoparamétrica se degrada porque parte de los nodos medios quedan sobre aristas rectas. Sobre geometría plana el elemento es exacto. Para geometría curva con tetraedros conviene un mallador nativo que sitúe los nodos medios sobre la superficie.

```yaml
elements:
  - {id: 1, type: Hex8,  material: 1, nodes: [1, 2, 3, 4, 5, 6, 7, 8]}
  - {id: 2, type: Tet4,  material: 1, nodes: [1, 2, 3, 4]}
  - {id: 3, type: Hex20, material: 1, nodes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                                              11, 12, 13, 14, 15, 16, 17, 18, 19, 20]}
  - {id: 4, type: Tet10, material: 1, nodes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]}
```

## Módulo de usuario `discontinuities`: elemento con discontinuidad embebida

Formulación de **fractura computacional** de Retama (2010): el elemento materializa una discontinuidad interna $\Gamma_d$ cuando se cumple un criterio de activación, y enriquece su cinemática con un salto de desplazamientos $\llbracket u \rrbracket$ gobernado por un material cohesivo. Los grados de libertad del salto son **elementales**: se condensan dentro del elemento y nunca llegan al ensamblador, de modo que el tamaño del sistema global no cambia.

No es un elemento finito estándar, así que no forma parte del programa principal: vive en el **módulo de usuario** `discontinuities` (`solidum/user/discontinuities/`), junto con su familia de materiales cohesivos (ADR 0020). `import solidum` no lo carga; el YAML que use `CST_Embedded2D` o la sección `cohesive_materials` debe declararlo:

```yaml
user_modules: [discontinuities]
```

Sin esa línea, la lectura se detiene con un error que nombra el tipo o la sección desconocidos y lista los módulos de usuario disponibles.

### Triángulo CST con Discontinuidad Interior: `CST_Embedded2D`

CST de 3 nodos con cinemática KOS enriquecida (Retama 2010, Caps. 2, 5, 6 y 7).

- **DOFs/nodo**: `ux, uy` ($\texttt{STRAIN\_DIM} = 3$). Los 2 DOFs del salto son elementales, no globales.
- **Cuadratura**: 1 punto (hereda del `Tri3`).
- **Parámetros**: `material` (bulk), `cohesive_material` y `thickness` (opcional, default `1.0`).
- **Dos materiales**: a diferencia del resto del catálogo, requiere **un material de bulk** que gobierna el continuo y **uno cohesivo** que gobierna el salto en $\Gamma_d$. En YAML son dos campos distintos: `material` se busca en `materials` y `cohesive_material` en `cohesive_materials`.
- **Estado intacto**: idéntico al `Tri3` hasta que la discontinuidad se activa.
- **Activación**: criterio de Rankine ($\sigma_I > \sigma_{t0}$ del cohesivo) evaluado en el centroide con el estado convergido del paso anterior; es el único criterio y no se declara en el YAML. Es **irreversible**: una vez activada persiste aunque el paso siguiente descargue.
- **Bulk aceptado**: sólo `Elastic2D` en la fase actual — la *discrete approach* presupone bulk elástico (ADR 0010).
- **Post-proceso**: `compute_gauss_state(U)` añade la clave `'discontinuity'` con normal, tangente, centroide, $l_d$, salto, tracción y daño cuando el elemento está agrietado.
- **Solver para la rama post-pico**: con desplazamiento impuesto, `NonlinearSolver` recorre el ablandamiento hasta la separación completa con penalización cohesiva rígida (tracción uniaxial frente a la solución exacta y ensayo de Van Vliet con grieta prescrita, en `tests/user/discontinuities/`). Las limitaciones registradas de `ArcLengthSolver` y de `DissipationArcLengthSolver` (módulo de usuario `dissipation_control`, al final del capítulo «Esquemas de Solución») con esa penalización son anteriores a la corrección de la tangente cohesiva del 2026-09-23 y no se han vuelto a comprobar.

```yaml
user_modules: [discontinuities]
cohesive_materials:
  - {id: 1, type: CohesiveDamageIsotropic, sigma_t0: 2.5e6, G_f: 100.0,
     K_e: 1.0e13, softening: linear}
elements:
  - {id: 1, type: CST_Embedded2D, nodes: [1, 2, 3],
     material: 1, cohesive_material: 1}
```

La formulación completa (cinemática KOS, condensación, longitud $l_d$) está en la spec del elemento, en el apéndice «Módulos de usuario» del Manual de Referencia.
