# ADR 0022 — Cargas distribuidas variables, mallas estructuradas y verificación por soluciones manufacturadas

- **Estado**: aceptado (2026-10-27)
- **Fecha**: 2026-10-27
- **Alcance**: `solidum/math/distributed_loads.py`; base `_PlaneSolidLoads` de los cinco sólidos 2D (`Quad4`, `Quad8`, `Quad9`, `Tri3`, `Tri6`); `Assembler.assemble_body_load` con campos; contrato `REFERENCE_CELL` / `REFERENCE_NODES` de `Element`; paquete `solidum.mesh` (`structured_mesh`); paquete `solidum.verification` (soluciones manufacturadas, SymPy opcional); extra `verification` de `pyproject.toml`. Fase 1 de 4 (sólidos 2D). Fase 2 (sólidos 3D, mismo día): base `_SolidLoads3D` de `Hex8`, `Hex20`, `Hex27`, `Tet4` y `Tet10`, celdas `hex` y `tet` del generador, verificación 3D.

## Contexto

El 2026-10-27 el usuario aprobó verificar Solidum por el **método de soluciones manufacturadas** (MMS). Se inventa una solución suave `u*`, se deducen de ella la fuerza de cuerpo `b* = −∇·σ(u*)` y las tracciones `t* = σ(u*)·n`, y se mide el orden con que converge la solución de elementos finitos. Una sola prueba ejerce a la vez las piezas que las soluciones analíticas comprueban por separado: malla distorsionada, cargas variables, Dirichlet no homogéneo, contorno curvo y material anisótropo.

El plan tiene cuatro fases: (1) sólidos 2D, (2) sólidos 3D, (3) térmico estacionario y transitorio, (4) opcional, elementos 1D. Quedan fuera los materiales no lineales, los corrotacionales y la discontinuidad embebida.

**Criterio de aceptación** (decisión física del usuario): el orden observado entre los dos últimos niveles de refinamiento es ≥ teórico − 0.1, en L² y en energía, y el error baja en cada nivel.

Faltaban tres piezas:

1. **Cargas variables en el espacio.** `compute_body_load(b)` y `compute_edge_traction(edge, t)` sólo admitían valores constantes por elemento. Además, la carga de borde de los elementos cuadráticos repartía la cuerda con pesos fijos 1-4-1, exactos sólo con el nodo medio centrado sobre un borde recto (deuda #13).
2. **Mallas para todos los elementos.** Había unos doce generadores improvisados en pruebas y ejemplos, cada uno para uno o dos elementos. Un barrido sobre el registro necesita mallar cualquier elemento, con su propio orden de nodos.
3. **Medir el error y deducir el orden teórico** de cualquier elemento.

## Decisión

### 1. Una carga distribuida es una constante o un campo

Las firmas no cambian; cambia lo que aceptan. La fuerza de cuerpo es una constante o un campo `b(x)`. La tracción de borde es una constante o un campo `t(x, n)`, donde `n` es la **normal exterior unitaria del borde discreto** en el punto de Gauss.

Que el campo de borde reciba la normal sirve para dos cosas:

- Una presión es `t(x, n) = −p·n`, sin un método aparte; cubre en el API de Python la deuda #24.
- Una solución manufacturada impone `t = σ(u*)·n` con la normal de la malla. Así `u*` resuelve exactamente el problema en el dominio discreto, también con contornos curvos. El error que se mide es sólo de discretización, sin un error de consistencia geométrica que lo contamine.

`Assembler.assemble_body_load` pasa el campo sin cambios a cada elemento, que lo evalúa con la posición en su dimensión.

### 2. Integración de las cargas

Las cargas se integran en `solidum/math/distributed_loads.py`, una implementación compartida por los cinco elementos (`_PlaneSolidLoads`).

- **Fuerza de cuerpo.** Se integra con `_LOAD_QUADRATURE`, que es la regla de la masa consistente del elemento:

  | Elemento | Regla |
  |---|---|
  | Quad4 | 2×2 |
  | Tri3 | `tri_3` |
  | Quad8, Quad9 | 3×3 |
  | Tri6 | `tri_6` |

  Esa regla integra exactamente `N_i·N_j`, y por tanto `N_i·b` con `b` del grado de la interpolación, sea cual sea la regla de la rigidez. Antes se usaba la regla de la rigidez, que puede ser reducida.

  Con `b` constante y lados rectos el resultado no cambia. Cambia, para mejor, en tres casos:
  - con `b` variable;
  - en el Quad4 con 1×1 sobre elementos distorsionados;
  - en el Tri6 con lados curvos.
- **Carga de borde.** El borde se interpola con sus propios nodos sobre `s ∈ [−1, 1]`: lineal con 2 nodos, cuadrática con 3. Se integra `N_k·t·|dx/ds|` con Gauss de 3 puntos sobre la geometría real del borde. Es exacta para tracción constante sobre un borde recto, con el nodo medio centrado o no. La resultante de una presión uniforme es exacta sobre cualquier borde cuadrático, porque `∫ n ds` es la cuerda girada −90°.

  Esto sustituye al reparto 1-4-1 sobre la cuerda, y la **deuda #13 queda saldada**. En los bordes de dos nodos con tracción constante da lo mismo que el `L/2` anterior.

  Evidencia: con el 1-4-1 antiguo, el MMS de Quad8 y Tri6 cae de orden 3 a orden 2 en L², tanto en el cuadrado distorsionado como en la corona. El orden en energía no cambia, como predice la teoría para un error de consistencia `O(h^p)` en las cargas.

### 3. Contrato `REFERENCE_CELL` / `REFERENCE_NODES`

Son dos datos declarativos de `Element` (Reglas §1), `None` por omisión:

- `REFERENCE_CELL`: la celda de referencia del elemento (`"quad"` `[-1, 1]²`, `"tri"` de coordenadas de área, `"hex"` o `"tet"`).
- `REFERENCE_NODES`: las coordenadas naturales de los nodos en su orden local.

Los consumen el generador de mallas y la deducción del orden teórico. `tests/test_reference_nodes_contract.py` comprueba tres cosas:

- `N_i(ξ_j) = δ_ij` con las funciones de forma del propio elemento;
- los bordes de `EDGE_NODES` son lados de la celda, recorridos en sentido antihorario y con el nodo medio en el centro;
- ningún sólido plano deja de declararlos, porque quedaría fuera del barrido sin aviso.

### 4. `solidum.mesh.structured_mesh`

Malla el bloque `[0, 1]^d` en celdas y coloca los nodos de cada elemento a partir de `REFERENCE_NODES`. Funde los nodos compartidos por sus coordenadas de referencia y lleva cada nodo a `x = mapping(origin + ξ·size)`. Con una transformación curva, los nodos intermedios caen sobre la curva.

Devuelve:
- el `Domain`, sin numeración de ecuaciones;
- los bordes de contorno por lado (`xmin`, `xmax`, `ymin`, `ymax`), como pares (elemento, borde), listos para `compute_edge_traction`. Desde la fase 2 están en `boundary`, que en 3D guarda las caras;
- los nodos de cada lado, con `nodes_on(side)`.

Celdas admitidas:
- `"quad"`;
- `"tri"`: dos por celda, separadas por la diagonal.

Los hexaedros y los tetraedros, en la fase 2 (abajo).

No es un mallador general: las geometrías reales siguen yendo por Gmsh. Los generadores improvisados que ya existen no se migran, porque solo se movería código; las pruebas nuevas usan este.

### 5. `solidum.verification`

- **`ManufacturedElasticity(u*, C)`** (`ManufacturedElasticity2D` hasta la fase 2): deriva con SymPy, de forma exacta, `ε*` (Voigt del proyecto), `σ* = C·ε*`, `b*` y `t*(x, n)`. Con `from_material` toma la matriz `C` del material y rechaza un material que no sea elástico lineal.
- **`solve_manufactured`**: impone `u*` en los nodos de los lados Dirichlet, `t*` en los demás lados y `b*` en el dominio.
- **`error_norms`**: integra ‖e‖_L2 y ‖e‖_E = (∫(ε*−ε_h)ᵀC(ε*−ε_h) dA)^½. Usa Gauss de 6×6 puntos en el cuadrado o una regla colapsada de Duffy en el triángulo. El campo de elementos finitos se evalúa con el contrato del camino por lotes (`batch_shape_functions`, `BATCH_KINEMATICS`), así que el módulo no conoce ningún elemento concreto.
- **`interpolation_degree`**: deduce el grado `p` del polinomio completo que reproducen las funciones de forma. De él salen los órdenes teóricos `p + 1` (L²) y `p` (energía), sin un dato nuevo en el elemento.
- **`convergence_study`**: encadena lo anterior y devuelve los errores, los órdenes y una tabla.

Ni `import solidum` ni el descubrimiento automático lo cargan. SymPy no es dependencia del programa principal: está en el extra `verification` y en `dev`, porque la suite la necesita. El módulo funciona para cualquier elemento plano que cumpla el contrato, incluido uno de un módulo de usuario (ADR 0020).

### 6. El barrido

`tests/validation/test_mms_solid_2d.py` recorre el registro. Cubre:
- todo sólido plano de desplazamientos con celda declarada;
- además, el Quad8 con 2×2, que es la mitigación del bloqueo de la deuda #6.

Lo ejecuta con dos geometrías:
- **Cuadrado distorsionado.** Lados rectos con los nodos deslizados a lo largo de ellos, de modo que los nodos medios quedan descentrados en el contorno. Elementos interiores curvos. Se ejecuta con deformación plana, esfuerzo plano y ortótropo girado 30°.
- **Cuarto de corona.** Tracción sobre los arcos; se ejecuta con deformación plana.

En las dos, la transformación es suave y las mallas son asintóticamente afines. Por eso los serendípitos deben dar el orden completo.

Una prueba de la parcela con `u*` lineal exige error de redondeo. Con ella se distingue un fallo de cargas o de contorno de uno de aproximación.

Duración: unos 15 s.

## Fase 2 (2026-10-27): sólidos 3D

Mismas decisiones, extendidas a los cinco sólidos 3D (`Hex8`, `Hex20`, `Hex27`, `Tet4`, `Tet10`).

- **Cargas.** La base `_SolidLoads3D` reúne `compute_body_load(b)` y `compute_face_traction(face, t)` de los cinco elementos:
  - la carga de cara es una constante o un campo `t(x, n)`, con `n` la normal exterior de la cara discreta, `∂x/∂s × ∂x/∂r` normalizada;
  - la orientación sale de los `FACE_NODES`, antihorarios vistos desde fuera (ADR 0012);
  - se integra sobre la geometría real con las funciones de forma de la cara (`face_load` en `distributed_loads.py`).

  La fuerza de cuerpo usa la regla de la masa:

  | Elemento | Fuerza de cuerpo | Caras |
  |---|---|---|
  | Hex8 | 2×2×2 | Gauss 2×2 |
  | Hex20, Hex27 | 3×3×3 | Gauss 3×3 |
  | Tet4 | `tet_4` | `tri_3` |
  | Tet10 | `tet_15`, antes `tet_4` | `tri_6`, antes `tri_3` |

  Las reglas antiguas del Tet10 también daban el orden óptimo, porque una regla exacta de grado 2p − 2 basta. El cambio mejora la exactitud, no el orden, y se hace por coherencia con 2D.
- **Contrato.** Los cinco declaran `REFERENCE_CELL` (`"hex"`, `"tet"`) y `REFERENCE_NODES`. La prueba de contrato comprueba además las caras:
  - están sobre la frontera de la celda;
  - con los nodos intermedios en su sitio, la transformación de la cara es la lineal o bilineal de sus vértices;
  - `∂x/∂s × ∂x/∂r` sale de la celda.
- **Generador.** Celdas `"hex"` (una por celda) y `"tet"` (seis tetraedros de Kuhn alrededor de la diagonal principal, todos con volumen positivo). Lados `zmin` y `zmax`. Los pares (elemento, cara) del contorno quedan en `StructuredMesh.boundary`, que en 2D tiene los bordes; antes el atributo se llamaba `boundary_edges`.
- **Verificación.** `ManufacturedElasticity`, antes `ManufacturedElasticity2D`, toma la dimensión del número de componentes de `u*`. Para elasticidad 3D usa el Voigt del proyecto y comprueba el orden de los cortantes (`γ_yz`, `γ_xz`) con una prueba que fallaría si estuvieran permutados. Las normas del error evalúan la solución exacta de una vez en todos los puntos de cada elemento y reutilizan las funciones de forma por clase. Para las cargas, que se evalúan punto a punto, la solución usa una versión compilada con `math`.
- **Barrido** (`tests/validation/test_mms_solid_3d.py`). Cubre los cinco elementos y el Hex20 con 2×2×2, con `Elastic3D`, en dos geometrías:
  - un cubo distorsionado: caras planas con nodos deslizados y elementos interiores curvos;
  - un sector de cilindro hueco, con tracción sobre las caras cilíndricas.

  Todos dan el orden teórico en L² y en energía. Los cuadráticos usan n = 2, 4, 8. Los lineales usan n = 4, 8, 16, porque el Tet4 converge desde abajo y despacio: 1.50, 1.75 y 1.91 en L².
- **Hallazgo: la prueba de la parcela en lados curvos.**
  - Con `u*` lineal, todas las combinaciones son exactas en una malla afín.
  - En las mallas curvas, el Tet10 con su regla por omisión (`tet_4`) y el Hex20 con 2×2×2 dejan un error de consistencia que baja como h^3.3.
  - La causa: sobre un elemento de lados curvos, el integrando `Bᵀ·detJ` es de grado mayor que el que esas reglas integran exactamente. Con la regla completa (`tet_15`, 3×3×3) la parcela es exacta.
  - No afecta al orden de convergencia, que el barrido comprueba con esas mismas reglas.
  - La spec del Tet10 ya recomendaba `tet_15` en mallas muy curvas; ahora el efecto está medido y blindado por la prueba (`SOLO_AFIN`).

## Consecuencias

- **Deuda #13 saldada.** La #24 queda cubierta en el API de Python mediante `t(x, n) = −p·n`, en los bordes 2D y en las caras 3D. La declaración en YAML sigue pendiente con la #23.
- **Coste del componente N+1.** Un sólido nuevo, plano o 3D, declara `REFERENCE_CELL` y `REFERENCE_NODES` y obtiene sin más código la malla estructurada y la verificación por MMS.
- **Limitaciones:**
  - la MMS de este módulo es de elasticidad lineal;
  - la pérdida de orden de los serendípitos en mallas que siguen trapezoidales al refinar es teoría conocida y no se prueba.
- **Pendiente:**
  - fase 3: `Q(x)` y `q̄(x, n)` en el térmico, transitorio con `T(x, t)`;
  - fase 4, opcional: elementos 1D.

## Alternativas descartadas

- **Gmsh para las mallas de prueba.** Es una dependencia externa en el CI, obliga a traducir el orden de los nodos y da menos control sobre la subdivisión.
- **Un método `compute_edge_pressure`.** El campo `t(x, n)` lo cubre con un solo contrato.
- **Un campo de borde `t(x)` sin la normal.** La MMS sobre un contorno curvo arrastraría un error de consistencia geométrica, y una presión no se podría expresar.
- **Integrar la fuerza de cuerpo con la regla de la rigidez.** El centroide del Tri3 no integra `N_i·b` con `b` lineal, y el 1×1 del Quad4 tampoco sobre elementos distorsionados.
- **Declarar el orden teórico en cada elemento.** Se deduce de las funciones de forma, así que un dato declarado podría contradecirlas.

## Referencias

Verificadas por el bibliotecario el 2026-10-27 contra Crossref, OpenAlex y OSTI.

- P. J. Roache, "Code verification by the method of manufactured solutions", *Journal of Fluids Engineering* 124(1) (2002) 4–10. DOI [10.1115/1.1436090](https://doi.org/10.1115/1.1436090).
- K. Salari, P. Knupp, *Code verification by the method of manufactured solutions*, informe SAND2000-1444, Sandia National Laboratories, 2000. DOI [10.2172/759450](https://doi.org/10.2172/759450). Muestra que el MMS detecta todo error que afecte al orden de precisión.
- D. N. Arnold, D. Boffi, R. S. Falk, "Approximation by quadrilateral finite elements", *Mathematics of Computation* 71(239) (2002) 909–922. DOI [10.1090/S0025-5718-02-01439-4](https://doi.org/10.1090/S0025-5718-02-01439-4). Sobre cuadriláteros de lados rectos:
  - los serendípitos de grado 2 bajan a orden 2 en L² y 1 en H¹ en mallas que no son asintóticamente paralelogramos;
  - en mallas asintóticamente paralelogramo recuperan el orden óptimo;
  - el bicuadrático de 9 nodos lo conserva siempre.

  Las mallas de este barrido tienen lados curvos, un caso que el artículo no trata. Para ellas la evidencia es la medición: Quad8 da orden 3 y 2.
