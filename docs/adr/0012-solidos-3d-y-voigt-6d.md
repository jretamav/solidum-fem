# ADR 0012 — Sólidos 3D: convención Voigt 6D y cierre del contrato `internal_forces`

- **Estado**: aceptado
- **Fecha**: 2026-05-19
- **Alcance**: apertura de la familia de elementos sólidos 3D (Etapa 7). Toca: clase base `Material` (nueva dimensión Voigt), clase base `Element` (contratos declarativos `STRAIN_DIM = 6` y nuevas caras), `Reglas.md §5` (orden Voigt 3D), `ElementForces` y el contrato de `internal_forces` (cierre de deuda técnica #1 heredada del ADR 0002), catálogos de elementos y materiales.

## Contexto

El proyecto entra en **Etapa 7 — Sólidos 3D** con un alcance acotado decidido en sesión: `Hex8` (hexaedro trilineal isoparamétrico), `Tet4` (tetraedro lineal) y `Elastic3D` (elástico isótropo, sin variantes plane stress/plane strain — en 3D no aplican). Los materiales no lineales 3D (`VonMises3D`, `DruckerPrager3D`, `IsotropicDamage3D`) quedan para una sub-etapa posterior; los elementos cuadráticos 3D (`Hex20`, `Hex27`, `Tet10`) también. La justificación del alcance está en el chat: arrancar con el espejo natural de `Quad4 + Tri3 + Elastic2D`, validar arquitectura, y dejar las extensiones para iteraciones incrementales sobre maquinaria ya probada (regla de dos casos reales para centralizar — ver `feedback_consistencia_arquitectural` en memoria).

Tres decisiones arquitecturales que no son evidentes de la mera "réplica del 2D" exigen registro persistente:

1. **Convención Voigt 3D**. El proyecto fija en `Reglas.md §5` los signos y la convención stress-resultant/RHR para esfuerzos seccionales, y deja explícita la convención Voigt 2D `[ε_xx, ε_yy, γ_xy]` con `γ_xy = 2·ε_xy` (deformación angular *engineering*). Para 3D existen dos órdenes que conviven en la literatura: `[xx, yy, zz, xy, yz, xz]` (extensión natural del 2D, usada por Crisfield, Bathe, Cook-Malkus-Plesha) y `[xx, yy, zz, yz, xz, xy]` (Voigt clásico de cristalografía, usado por algunos textos de mecánica del continuo y por ABAQUS internamente). Hay que elegir una y blindarla.

2. **Cierre de la deuda técnica #1** (`internal_forces` devuelve `None` en sólidos). El ADR 0002 introdujo `ElementForces` para consumidores de la API pública (diagramas, post-proceso) con `internal_forces() → ElementForces` como contrato canónico. Los sólidos 2D devolvieron `None` desde el origen porque la primitiva natural de salida en sólidos es el campo $\boldsymbol\sigma$ por punto de Gauss, no un escalar/vector "esfuerzo seccional". `compute_gauss_state(U)` cubre esa necesidad. Pero la entrada de sólidos 3D fuerza el momento de decisión: o el contrato `internal_forces` deja de aplicar a sólidos (decisión documentada), o se le da semántica útil para sólidos (e.g. promedio o valor representativo). Procrastinarlo otra vez deja la asimetría heredada del 2D extendida al 3D, multiplicando el coste futuro.

3. **API de tracción de superficie en 3D**. Los sólidos 2D exponen `compute_edge_traction(edge, t̄)` con índice entero de borde y conectividad nodal local. La extensión natural en 3D es `compute_face_traction(face, t̄)`, pero "cara" es un objeto bidimensional (3 o 4 nodos según el elemento) y la cuadratura para integrar la tracción consistente sobre la cara ya no es trivial (1D vs 2D superficial). El patrón debe quedar uniforme entre `Hex8` y `Tet4`, con identificación clara de caras y reparto consistente.

Resolver las tres juntas en este ADR evita ADRs sucesivos cuando entren `Hex20`, materiales 3D no lineales, etc.

## Decisión

### 1. Convención Voigt 3D

Se adopta el orden **`[xx, yy, zz, xy, yz, xz]`** con factores *engineering* en los cortantes. Es decir:

$$\boldsymbol\varepsilon = [\varepsilon_{xx},\ \varepsilon_{yy},\ \varepsilon_{zz},\ \gamma_{xy},\ \gamma_{yz},\ \gamma_{xz}]^\top, \qquad \gamma_{ij} = 2\,\varepsilon_{ij}\ (i\neq j)$$

$$\boldsymbol\sigma = [\sigma_{xx},\ \sigma_{yy},\ \sigma_{zz},\ \sigma_{xy},\ \sigma_{yz},\ \sigma_{xz}]^\top$$

Los esfuerzos son los componentes tensoriales reales sin factor; solo las deformaciones angulares llevan el factor 2.

**Por qué este orden y no `[xx, yy, zz, yz, xz, xy]`**:

- **Es la extensión natural del Voigt 2D del proyecto** (`[xx, yy, xy]`). La columna del cortante en el plano sigue siendo la primera del bloque cortante. Tests, materiales y elementos 2D leen el último componente como cortante "principal" del plano. La extensión a 3D simplemente añade `yz` y `xz` detrás. Materiales mixtos 2D↔3D (códigos que comparten subrutinas) no rompen invariantes.

- **Es la convención de los libros de referencia activos del proyecto**: Bathe §6.6, Crisfield §12, Cook-Malkus-Plesha §6.8. Las ecuaciones del catálogo y del manual referenciarán estos textos sin necesidad de traducción.

- **Es el orden de los códigos de FEM académicos más documentados** (CalculiX, FEniCS, FEAP por defecto). ABAQUS expone `[11, 22, 33, 12, 13, 23]` que coincide salvo por `13↔23` — diferencia menor que se anota en el catálogo de materiales por si alguien transfiere datos.

**Matriz constitutiva elástica isótropa 6×6** (`Elastic3D`) en esta convención:

$$\mathbf C_e = \frac{E}{(1+\nu)(1-2\nu)}\begin{bmatrix}
1-\nu & \nu & \nu & 0 & 0 & 0 \\
\nu & 1-\nu & \nu & 0 & 0 & 0 \\
\nu & \nu & 1-\nu & 0 & 0 & 0 \\
0 & 0 & 0 & \tfrac{1-2\nu}{2} & 0 & 0 \\
0 & 0 & 0 & 0 & \tfrac{1-2\nu}{2} & 0 \\
0 & 0 & 0 & 0 & 0 & \tfrac{1-2\nu}{2}
\end{bmatrix}$$

El factor $\tfrac{1-2\nu}{2}$ en los cortantes (no $1-2\nu$) es consecuencia directa de usar $\gamma_{ij}$ *engineering* en la deformación: $\sigma_{ij} = G\,\gamma_{ij} = \tfrac{E}{2(1+\nu)}\,\gamma_{ij}$, y $\tfrac{E}{2(1+\nu)} = \tfrac{E}{(1+\nu)(1-2\nu)}\cdot\tfrac{1-2\nu}{2}$.

**Tests blindan la convención**: el primer test analítico de `Elastic3D` verifica simetría de $\mathbf C_e$, positividad definida en $\nu \in \{0, 0.3, 0.49\}$, y la traducción de tracción uniaxial $\varepsilon = (\varepsilon_0, -\nu\varepsilon_0, -\nu\varepsilon_0, 0, 0, 0) \Rightarrow \sigma = (E\varepsilon_0, 0, 0, 0, 0, 0)$ exacta. Si alguien en el futuro permuta sin querer las columnas del bloque cortante, el test salta.

**Documentación**: `Reglas.md §5` recibe una subsección nueva "Voigt 3D" con el orden y el factor *engineering*, en línea con la subsección Voigt 2D ya existente implícitamente. La regla queda obligatoria para todo el proyecto a partir de este ADR.

### 2. Cierre del contrato `internal_forces` para sólidos

**Decisión**: el contrato `internal_forces() → ElementForces` **no se aplica a sólidos** ni en 2D ni en 3D. Es semánticamente vacío: un sólido continuo no tiene "esfuerzo seccional" — el equivalente en sólidos es el campo $\boldsymbol\sigma(\mathbf x)$, expuesto por `compute_gauss_state(U) → dict[gauss_point → {ε, σ, coords}]`. Tratar de coercer un sólido a `ElementForces` (e.g. devolviendo el promedio de $\boldsymbol\sigma$, o el valor en el centroide) inventa un dato que ningún consumidor pide y que enmascara la primitiva real.

**Concretamente**:

- `Element` (clase base) explicita en docstring que `internal_forces` es **opcional**, no abstracta. Default: `return None`. Sólidos 2D y 3D heredan el default sin sobreescribir.
- Estructuras 1D (truss, cable, frame en 2D y 3D) siguen sobreescribiendo `internal_forces` con su semántica seccional (`N`, `V`, `M`, `T`) — sin cambios para ellos.
- `compute_gauss_state(U)` queda **garantizada** para todo elemento sólido (contrato abstracto en `Element` solo para subclases con `STRAIN_DIM ∈ {3, 6}`; para 1D estructurales es opcional y por convención no se implementa). Consumidores externos que necesiten salida de sólidos (post-proceso VTK, suavizado nodal, integradores de error) usan `compute_gauss_state`.
- `ElementForces` (dataclass en `solidum.core`) recibe una nota en docstring explicitando que **aplica solo a elementos estructurales 1D**. No se mueve ni se renombra — el ADR 0002 sigue válido para su dominio.

**Por qué cierre por exclusión y no por extensión**:

- La asimetría conceptual no es defecto de la implementación; es propiedad del dominio. Un puente continuo no tiene momento flector en cada punto en el sentido de viga; tiene un campo tensorial completo.
- Inventar un `ElementForces` para sólidos (promedio, centroide, valor representativo) introduce ambigüedad: tres consumidores podrían querer tres definiciones distintas. Mejor que cada consumidor procese `compute_gauss_state` con la semántica que necesite.
- El coste de mantener dos APIs (`internal_forces` para estructurales, `compute_gauss_state` para sólidos) es bajo porque ya están ambas implementadas y los consumidores actuales no se confunden — el catálogo distingue las dos familias.

**Migración**: el código actual no cambia. La diferencia es documental — el ADR 0002 deja de ser "incompleto"; pasa a ser "completo con dominio explícitamente acotado a elementos estructurales 1D". El item #1 de la deuda técnica de `STATUS.md` se tacha con la nota "cerrado por ADR 0012 con dominio explícito; sólidos 2D/3D exponen `compute_gauss_state`".

### 3. API de tracción de superficie y caras 3D

Cada elemento sólido 3D expone:

```python
def compute_face_traction(self, face: int, t_bar: np.ndarray) -> np.ndarray:
    """Tracción uniforme t̄ (vector global 3D) sobre la cara identificada
    por índice entero. Devuelve el vector de fuerzas nodales equivalente
    consistente (longitud = n_nodes * 3 = ndof_e). t̄ se especifica
    siempre en coordenadas globales."""
```

**Identificación de caras**:

- `Hex8` — 6 caras cuadrilaterales numeradas según la convención estándar de hexaedro (consistente con VTK_HEXAHEDRON y con la mayoría de pre-procesadores):
  - Cara 0: `(0, 3, 2, 1)` — −ζ (cara inferior).
  - Cara 1: `(4, 5, 6, 7)` — +ζ (cara superior).
  - Cara 2: `(0, 1, 5, 4)` — −η (cara frontal).
  - Cara 3: `(1, 2, 6, 5)` — +ξ (cara derecha).
  - Cara 4: `(2, 3, 7, 6)` — +η (cara trasera).
  - Cara 5: `(3, 0, 4, 7)` — −ξ (cara izquierda).

  El orden de los 4 nodos en cada cara es tal que la normal sale del elemento (regla de la mano derecha). Crítico para definir presión normal en una iteración posterior.

- `Tet4` — 4 caras triangulares numeradas por el nodo opuesto:
  - Cara 0: `(1, 2, 3)` (opuesta al nodo 0).
  - Cara 1: `(0, 3, 2)` (opuesta al nodo 1).
  - Cara 2: `(0, 1, 3)` (opuesta al nodo 2).
  - Cara 3: `(0, 2, 1)` (opuesta al nodo 3).

  Orden de nodos con normal saliente.

**Cuadratura para la integración de la tracción**:

- Hex8 — cuadratura 2D Gauss `2×2` sobre la cara (4 puntos por cara). Exacta para tracción constante sobre cara plana; degrada con tracción variable o cara alabeada, pero esos casos están fuera de alcance hoy.
- Tet4 — cuadratura 1 punto centroide de la cara (peso `A_cara`). Exacta para tracción constante por linealidad de las funciones de forma.

**Tracción variable y presión normal**: fuera de alcance en este ADR — patrón idéntico al 2D actual (`compute_edge_traction` solo acepta tracción uniforme en globales). Cuando entre presión hidrostática u otro caso real, se decide la API en su propio ADR.

## Consecuencias

**Inmediatas (esta etapa)**:

- Tres specs nuevas (`Elastic3D`, `Hex8`, `Tet4`) implementables sin ambigüedad sobre la convención Voigt o el contrato de salida.
- `STATUS.md` cierra deuda técnica #1 con nota explícita.
- `Reglas.md §5` ampliada con la subsección Voigt 3D — referencia obligatoria para tests, manuales, formulaciones futuras.

**Hacia adelante**:

- Materiales 3D no lineales (`VonMises3D`, `DruckerPrager3D`, `IsotropicDamage3D`) heredan la convención sin discusión adicional. Sus tangentes algorítmicas son 6×6 con el mismo orden.
- Elementos cuadráticos 3D (`Hex20`, `Hex27`, `Tet10`) reutilizan las caras y la cuadratura por cara (escalando a 2D 3×3 para hexaedros de orden 2). El patrón base ya está definido.
- Cuando entre embedded discontinuity 3D (`Tet4_Embedded`), el plano de discontinuidad es un triángulo cuyo área y normal se construyen con la misma maquinaria de caras tetraédricas.
- Modelos plásticos 3D pueden derivar su tensor desviador como $\mathbf s = \boldsymbol\sigma - \tfrac{1}{3}\text{tr}(\boldsymbol\sigma)\mathbf I$ con los tres componentes diagonales primeros del Voigt sin más permutaciones. El return mapping J2 3D queda más limpio.

**Costes**:

- Una pequeña incompatibilidad si alguien transfiere datos en formato ABAQUS (orden `[11, 22, 33, 12, 13, 23]` en lugar del nuestro `[11, 22, 33, 12, 23, 13]`): se documenta en el catálogo de materiales una vez y queda como referencia. Coste de un script de permutación si llega el caso.
- La convención Voigt 3D pasa a ser **invariante del proyecto**; cambiarla retroactivamente costaría reescribir todos los materiales 3D que lleguen después. Es decisión deliberada — el invariante temprano paga rendimientos compuestos.

**Deuda diferida (no creada, ya existente)**:

- Locking volumétrico en Hex8 con $\nu \to 0.5$ — patrón idéntico a Quad4. Se documenta en la spec de Hex8 como limitación, se blinda con un test (`test_volumetric_locking_3d.py`), B-bar/F-bar diferidas. Misma política que el 2D.
- Hourglass en Hex8 con integración reducida 1 punto — esquema disponible vía `quadrature` pero con warning explícito; estabilización Flanagan-Belytschko no implementada. Misma política que Quad4.

## Alternativas consideradas

- **Voigt 3D con orden `[xx, yy, zz, yz, xz, xy]`** (Voigt cristalográfico clásico). Rechazada: rompe la extensión natural del Voigt 2D del proyecto; obliga a documentar la permutación en cada test y en cada material mixto 2D↔3D.

- **`internal_forces` con semántica "promedio de σ por elemento"** para sólidos. Rechazada: ya hay tres definiciones razonables (promedio aritmético, valor en centroide, promedio ponderado por volumen del punto de Gauss); cualquier elección satisface algún consumidor y traiciona a otros dos. `compute_gauss_state` deja la elección al consumidor con coste de una línea de Python.

- **Diferir el cierre de la deuda #1 a una etapa posterior**. Rechazada: el criterio explícito de retoma documentado en STATUS era "cuando entren sólidos 3D"; saltarse el cierre cuando se cumple el criterio convierte la deuda en deuda perpetua, exactamente lo que el patrón de "fechas-gatillo" de la memoria del proyecto pretende evitar.

- **Caras Hex8 con la convención GMSH `(0,1,2,3)` para la base y `(4,5,6,7)` para la tapa, con orientación independiente de la normal**. Rechazada: pierde el invariante de normal saliente; complica la futura introducción de presión normal y la depuración geométrica.

## Referencias

- Bathe K.-J. (2014). *Finite Element Procedures*. §6.6 (sólidos 3D isoparamétricos).
- Crisfield M.A. (1991). *Non-linear Finite Element Analysis of Solids and Structures*, vol. 1, §12 (sólidos 3D).
- Cook R.D., Malkus D.S., Plesha M.E., Witt R.J. (2002). *Concepts and Applications of FEA*. §6.8, §11.3 (hexaedros y tetraedros).
- Zienkiewicz O.C., Taylor R.L. (2005). *The Finite Element Method: Its Basis and Fundamentals*. §6.7 (3D elements).
- ADR 0002 — API de resultados para consumidores (este ADR cierra explícitamente el dominio que aquel dejaba parcial).
- ADR 0008 — Densidad como propiedad del material (patrón heredado: propiedad intrínseca del medio, no de la discretización).
- `Reglas.md §5` — convenciones de signos del proyecto, ampliada con la subsección Voigt 3D en este commit.
