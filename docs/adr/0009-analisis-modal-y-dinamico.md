# ADR 0009 — Análisis modal y dinámico

- **Estado**: Aceptado
- **Fecha**: 2026-05-13
- **Alcance**: clase base `Element` y todas las subclases; `Assembler`; `SolverRegistry`; parser YAML; catálogo de solvers; manual de usuario.

## Terminología — distinción rigurosa

Este ADR cubre dos análisis distintos que conviene separar terminológicamente:

- **Análisis modal**: problema algebraico de valores y vectores característicos generalizado `K·φ = ω²·M·φ`. Las soluciones `(ω_n, φ_n)` admiten **interpretación** como frecuencias propias y formas de vibración libre no amortiguada, pero el análisis en sí es algebraico, no temporal.
- **Análisis dinámico**: análisis con dependencia temporal explícita — integración de `M·ü + C·u̇ + K·u = F(t)` (transitorio) o resolución de `(−ω²M + iωC + K)·û = F̂` para barridos en `ω` (respuesta en frecuencia). Requiere integrador temporal o aritmética compleja en el dominio de la frecuencia.

La fase 1 implementa **solo análisis modal**. Las fases 3-7 abren propiamente el análisis dinámico. El ADR las cubre conjuntamente porque comparten infraestructura: matriz de masa, amortiguamiento Rayleigh derivado de frecuencias modales, estado dinámico fuera de `Node`.

## Resumen ejecutivo

Solidum abre la línea de trabajo con un primer hito acotado — **análisis modal** — que introduce la matriz de masa **M** consistente y un solver de autovalor generalizado `K·φ = ω²·M·φ`. Este ADR documenta también la **hoja de ruta completa** hacia análisis dinámico transitorio (lineal y no lineal) y respuesta en frecuencia, fijando hoy las decisiones arquitecturales que condicionarán esas fases: M consistente con contrato extensible a lumped, todo el catálogo (modal y dinámico) en `SolverRegistry` (sin fragmentación en `IntegratorRegistry`), amortiguamiento Rayleigh como entrada estándar para transitorio, estado transitorio (`u̇`, `ü`) en el resultado del solver sin contaminar `Node`. Las fases posteriores reutilizan piezas ya existentes — capa algebraica (ADR 0003), `Material.density` (ADR 0008), convergencia calibrada (ADR 0007) — sin reabrir debates.

## Contexto

Tres ADRs previos ya prepararon el terreno:

- **ADR 0003** identificó análisis modal y dinámica implícita entre los análisis que la capa algebraica iba a habilitar (filas "Análisis modal" y "Dinámica implícita" de su tabla de hoja de ruta) y reservó `EigenSolver` como capacidad fase 3.
- **ADR 0008** introdujo `Material.density` como propiedad de primera clase, con enforcement diferido (`None` al construir → `ValueError` al pedir peso propio o matriz de masa). La densidad es el único dato físico nuevo que tanto modal como dinámico necesitan; ya está disponible.
- **ADR 0007** dejó tolerancias y convergencia calibradas para cualquier solver iterativo: aplicable inmediatamente al Newton-Raphson dentro de cada paso de Newmark.

Hoy Solidum no tiene ni matriz de masa, ni solver de autovalores, ni integrador temporal. Pero todas las piezas que faltan son **incrementales sobre lo que ya existe**, no refactor estructural. Ese es el supuesto operativo que este ADR formaliza.

## Decisión

### 1. Matriz de masa consistente, contrato extensible a lumped

Se añade al contrato de `Element` un método:

```python
def compute_mass_matrix(self, lumping: str = "consistent") -> np.ndarray:
    """Matriz de masa elemental en ejes globales, ordenada según DOF_NAMES × nodos.

    Parameters
    ----------
    lumping : {"consistent"}, default "consistent"
        Estrategia de discretización de la inercia. Solo "consistent" está
        implementado en esta fase. Llamadas con otro valor lanzan NotImplementedError.
    """
```

**Por qué consistente como única implementación inicial**:

- Para análisis modal y para integración temporal **implícita** (Newmark, HHT-α), la masa consistente da convergencia óptima en autovalores y es la referencia analítica. Es la elección por defecto de SAP2000, OpenSees, Abaqus, Code_Aster.
- En **frames y vigas**, la masa lumped obliga a decidir qué hacer con la inercia rotacional nodal. Si se asigna cero, M es singular sobre los DOF rotacionales y el problema generalizado `K·φ = ω²M·φ` deja de tener forma estándar (eigsh requiere shift-invert con cuidado). Si se asigna por HRZ, introduce una heurística que ensucia la validación frente a fórmulas analíticas. La consistente sale directamente del Galerkin sobre las funciones de forma Hermitianas — sin heurística. **La fase 1 incluye además la inercia rotacional propia de sección** (`ρI·L/6·[[2,1],[1,2]]` aplicada a los DOFs rotacionales): es término real del medio continuo, indispensable en vigas peraltadas (consistente con Timoshenko) y mejora pequeña pero correcta en Bernoulli (el efecto relativo es `O(r²/L²)`). Su omisión sería inconsistente con `Frame2DTimoshenko`.
- Lumped solo paga claramente en **integración explícita** (`ü = M⁻¹(F − Cu̇ − Ku)` con M⁻¹ trivial). Explícita en barras es marginal: la frecuencia más alta de un frame es enorme (modos rotacionales), el paso condicional es minúsculo. El dominio natural de Solidum (sólidos con plasticidad, no propagación de ondas) tira hacia implícito. Lumped tendrá sentido cuando entre explícita en sólidos 2D/3D.

**Por qué el parámetro `lumping` está en la firma desde el día uno**:

El día que llegue diferencias centradas no es necesario tocar la firma del método ni el llamador (`Assembler.assemble_mass_matrix`). Solo se añade la rama `"lumped"` (y eventualmente `"hrz"`) en cada elemento. La extensión es aditiva, no rompedora. Coherente con Reglas §1 (extensión barata) sin caer en código muerto presente.

### 2. Ensamblaje de M en `Assembler`

```python
def assemble_mass_matrix(self, lumping: str = "consistent") -> sp.csr_matrix:
    """Ensambla M global iterando sobre elementos. Reutiliza la topología COO
    cacheada (misma sparsity que K) y la enforcement de density del ADR 0008
    (ValueError si algún material no la declara)."""
```

Tres observaciones de diseño:

- **Misma topología que K**. M y K comparten patrón de sparsity porque ambos se ensamblan sobre los mismos pares (DOF_i, DOF_j) elemento a elemento. El caché de `_coo_rows`, `_coo_cols`, `_elem_dof_indices` que ya existe para K se reutiliza tal cual. Coste de ensamblaje de M ≈ coste de ensamblaje de K.
- **Enforcement de density delegada**. Si algún material tiene `density=None`, `assemble_mass_matrix` falla con `ValueError` listando los materiales afectados — exactamente el mismo patrón que `assemble_self_weight` (ADR 0008). No se replica la lógica; se centraliza en un helper si surge tercera consumidor de density (centrífuga, amortiguamiento másico, etc.).
- **M lineal y constante**. En análisis modal lineal y en transitorio con formulación lagrangeana total, M depende solo de ρ (constante en el tiempo) y de la configuración de referencia (inicial, fija). Se ensambla **una vez por análisis**, no por paso. La caché vive en el `Assembler` análoga a `_topology_built`.

### 3. Modal como solver en `SolverRegistry` (no en un registry nuevo)

Se introduce `ModalSolver` como una clase más registrada en `SolverRegistry`, paralela a `LinearSolver`, `NonlinearSolver`, `ArcLengthSolver`. No se crea ningún registry paralelo `IntegratorRegistry`.

**Razón**: en Solidum los solvers actuales ya orquestan ambos roles que en otros códigos (OpenSees, FEAP) están separados — control de avance (pseudo-temporal, factor de carga, longitud de arco) y resolución algebraica. La distinción "solver vs integrator" es semántica útil al hablar pero no aporta valor arquitectural cuando se materializa en dos registries distintos. Fragmentar abriría dos catálogos paralelos que el usuario tendría que conocer; mantenerlos unificados conserva el patrón "un solver = un análisis = una clase registrada".

`ModalSolver` internamente delega en `EigenSolver` de la capa algebraica (`solidum.math.linalg.eigen`), del mismo modo que `LinearSolver` delega en `LUSolver` / `CholeskySolver` (ADR 0003 §1). La separación de capas se mantiene:

```
ModalSolver (capa de orquestación: ensambla K, M, aplica Dirichlet, gestiona modos)
    │
    └── EigenSolver (capa algebraica: scipy.sparse.linalg.eigsh, shift-invert)
            │
            └── FactorizedSolver para (K − σ·M), del LUSolver/CholeskySolver existente
```

### 4. Tratamiento de Dirichlet en problema modal

Las condiciones de contorno de Dirichlet en modal se imponen por **eliminación directa** (mismo mecanismo que estática, ADR 0004 fase 1):

```
T.T @ K @ T · φ_red = ω² · T.T @ M @ T · φ_red
```

donde `T` es el operador que selecciona DOFs libres. El término `g_indep` de restricciones lineales no homogéneas no aplica en modal (modos son solución del problema homogéneo asociado; los modos del problema con apoyos prescritos a desplazamiento no nulo son los mismos que con apoyos nulos — los apoyos solo definen qué DOFs son libres). El `Assembler.reduce(K, F=0)` actual sirve directamente; se añade un caso `reduce_pair(K, M)` que aplica `T` simultáneamente a ambas matrices.

Tras resolver el problema reducido, los modos completos se obtienen con `Assembler.expand(φ_red, T, g=0)`. En DOFs prescritos el modo vale cero (consistente con que esos DOFs no oscilan).

### 5. Estado transitorio fuera de `Node`

Para análisis modal el estado `(u̇, ü)` no se necesita — modal devuelve solo (`ω`, `φ`). Pero la decisión arquitectural se fija hoy para que la fase de transitorio no requiera tocar `Node`:

- `Node` permanece con `dofs`, `boundary_conditions`, `point_loads`. **No se le añaden** atributos `velocity`, `acceleration`.
- El estado dinámico vive en `SolveResult` (transitorio) como arrays `u_history`, `udot_history`, `uddot_history` indexados por paso temporal. Mismo patrón que `lambda_history` en arc-length.
- El solver transitorio mantiene `u_curr`, `udot_curr`, `uddot_curr` como estado interno mientras avanza; al cierre del análisis se exporta como `SolveResult`.

Esto preserva la pureza de `Node` como datos topológicos+geométricos+condiciones de contorno y evita que análisis estáticos paguen el coste de campos que no usan.

### 6. Amortiguamiento — Rayleigh primero, modal después

Para la fase de transitorio (no para modal, que no usa C) se adopta el siguiente camino:

- **Fase 3 (Rayleigh)**: amortiguamiento `C = α·M + β·K` definido por dos pares (ξ₁, ω₁) y (ξ₂, ω₂). Los coeficientes (α, β) se derivan automáticamente a partir de las primeras dos frecuencias propias obtenidas de un análisis modal previo. Es la entrada estándar de la industria y cierra el lazo con la fase modal.
- **Fase 4 (amortiguamiento modal)**: ξₙ por modo individual, requiere descomposición previa y proyección. Diferido a cuando aparezca un caso concreto que lo requiera (sismicidad con espectro de respuesta, vibración con ξ no uniforme entre modos).
- **Otras formas (Caughey series, amortiguamiento local elemento a elemento)**: diferidas. Solo se introducen si surge demanda concreta — Reglas §1.

### 7. Convención de unidades — heredada

Como con `Material.density` (ADR 0008), todas las magnitudes dinámicas heredan la convención: el usuario es responsable de la consistencia de unidades. Si trabaja en kg/N/m, las frecuencias salen en rad/s; si trabaja en t/N/mm, también — el sistema no convierte. Las tolerancias adimensionales del ADR 0007 garantizan invariancia bajo cambio de unidades.

## Hoja de ruta — fases

| Fase | Análisis | Pieza algebraica clave | Solver Solidum | Estado |
|---|---|---|---|---|
| 1 | **Análisis modal** (frecuencias y modos del problema generalizado) | `eigsh(K, M, sigma, which="LM")` con shift-invert | `ModalSolver` | **Implementada** |
| 2 | Modal con concentración de masa (lumped) | M diagonal | `ModalSolver(lumping="lumped")` | **Implementada** (2026-05-18) |
| 3 | Transitorio lineal Newmark | `factorize(M + γΔt·C + βΔt²·K)` reutilizable | `NewmarkSolver` | **Implementada** |
| 3-bis | HHT-α (variante de Newmark con disipación numérica) | mismo, con `(1+α)γΔt`, `(1+α)βΔt²` | `HHTSolver` | **Implementada** (variante, spec corta) |
| 4 | Transitorio no lineal Newton-Newmark | Newton dentro de cada paso, residuo `R = F_ext − F_int − Mü − Cu̇` | `NewtonNewmarkSolver` (+ `NewtonHHTSolver`) | **Implementada** |
| 5 | Transitorio explícito (diferencias centradas) | `M⁻¹` trivial con M lumped | `CentralDifferenceSolver` | **Implementada** (2026-05-18) |
| 6 | Respuesta en frecuencia (steady-state harmonic) | `(−ω²M + iωC + K)·û = F̂` complejo, barrido en ω | `HarmonicSolver` | **Implementada** (2026-05-18) |
| 7 | Análisis espectral / sísmico (combinación modal) | CQC, SRSS sobre espectro de respuesta | `ResponseSpectrumSolver` | **Implementada** (2026-05-18) — **ADR completo** |

Cada fase es un commit cerrado con tests, no bloquea las siguientes y deja Solidum en estado funcional.

## Fase 1 — Análisis modal (implementada)

Piezas concretas entregadas:

1. **Spec** — `docs/specs/ModalSolver.md` desde `_template_solver.md`. Define interface YAML, parámetros (`n_modes`, `sigma`, `which`, tolerancia), acceptance (barra empotrada-libre, viga Bernoulli-Euler) y referencias bibliográficas.
2. **Contrato `compute_mass_matrix`** en `Element` base con signatura `(self, lumping: str = "consistent")` y default `NotImplementedError` en la base abstracta.
3. **Implementación elemento por elemento**:
   - `Truss` (axial: `(ρAL/6)·[[2,1],[1,2]]`).
   - `Cable` (mismo patrón que truss).
   - `Frame2D` (axial consistente + flexión Hermitiana consistente).
   - `Frame3D` (axial + flexión en dos planos + torsional).
   - `Solid2D` (`∫NᵀρN·t dΩ` con la regla de cuadratura ya declarada en `N_INTEGRATION_POINTS` y la cinemática `STRAIN_DIM`).
4. **`Assembler.assemble_mass_matrix(lumping="consistent")`** — patrón análogo a `assemble_system`, reusa la topología COO cacheada, valida `material.density` con el mismo error de ADR 0008.
5. **`EigenSolver`** en `solidum/math/linalg/eigen.py` — envuelve `scipy.sparse.linalg.eigsh(K, k=n_modes, M=M, sigma=sigma, which="LM")`. Devuelve `(eigenvalues, eigenvectors)`.
6. **`ModalSolver`** en `solidum/math/solvers.py` registrado en `SolverRegistry`. Orquesta: `Assembler.assemble_system()` + `assemble_mass_matrix()` + `reduce_pair(K, M)` + `EigenSolver.solve` + `expand(φ_red)`. Devuelve `ModalResult(frequencies_hz, frequencies_rad, periods, modes, effective_masses)`.
7. **Cableado YAML** — `entry.py` / `yaml_parser.py` aceptan:
   ```yaml
   solver:
     type: modal
     n_modes: 10
     sigma: 0.0          # shift para shift-invert (0 = frecuencias más bajas)
     tolerance: 1.0e-9   # tolerancia ARPACK
   ```
8. **Tests** — `tests/test_modal_truss.py` y `tests/test_modal_frame.py`. Ver sección "Criterios de aceptación" en la spec del `ModalSolver`.

## Fase 2 — Mass lumping (implementada 2026-05-18)

La fase 1 dejó el parámetro `lumping` en la firma del método y rechazaba `"lumped"` con `NotImplementedError`. La fase 2 lo habilita siguiendo dos estrategias distintas según la familia de elemento.

### Decisiones de esquema

**Sólidos isoparamétricos (Tri3, Quad4, Tri6, Quad8, Quad9)** — **HRZ canónico** (Hinton-Rock-Zienkiewicz 1976). Aplica:

```
α = (D · m_total) / Σ_{i traslacional} M_ii^consistent
M_lumped[i,i] = α · M_consistent[i,i]
off-diagonals → 0
```

Centralizado en `solidum/math/mass_lumping.py::lump_hrz`. Para Tri3 y Quad4 coincide con row-sum por simetría; para Tri6/Quad8/Quad9 HRZ evita las masas negativas que daría row-sum puro en vértices con funciones de forma cuadráticas (Bathe FEP §9.2.4).

**Vigas y marcos (Truss2D/3D, Cable2D/3D, Frame2D Euler/Timoshenko/EulerCorot, Frame3D)** — **lumping nodal directo**:

- Masa traslacional: `ρAL/2` por DOF y nodo (la misma en cada dirección local — condición necesaria para que el bloque traslacional `m_t·I_d` sea invariante a SO(d) y la matriz global resulte diagonal tras `Tᵀ·M·T`).
- Inercia rotacional flexional: `ρI·L/2` por DOF rotacional y nodo (Frame2D); `ρIy·L/2`, `ρIz·L/2` (Frame3D).
- Inercia rotacional torsional (Frame3D): `ρ·Jp·L/2` con `Jp = Iy + Iz` (momento polar geométrico).

Implementado en helpers cerrados (`_frame2d_lumped_mass_local`, `Frame3D._build_local_mass_lumped`). El HRZ canónico aplicado a la diagonal consistente local daría masas distintas en `ux_local` (lineal) vs `uy_local` (Hermitiana) — diagonal en local pero **no en global tras rotar**. La diagonalidad en global es la propiedad operacional relevante para diferencias centradas (fase 5); por eso se prefiere el esquema nodal directo, estándar en frames (SAP2000, OpenSees, Cook-Malkus-Plesha §11.4).

### Limitación documentada (Frame3D oblicuo)

Para Frame3D con eje del elemento desalineado de los ejes globales, el bloque rotacional 3×3 por nodo no es escalar (`m_rx ≠ m_ry ≠ m_rz` en general, incluso con `Iy = Iz`: la torsional usa `Jp = Iy + Iz`, las flexionales `Iy` o `Iz`). La rotación general mezcla los tres DOFs rotacionales locales con los globales y el bloque rotacional queda lleno tras `Tᵀ·M·T`. M_global resulta **bloque-diagonal por nodo** con bloque traslacional diagonal + bloque rotacional 3×3 lleno; los nodos siguen sin acoplarse entre sí. Es la limitación estándar del lumping en frames 3D — aceptable para Newmark/HHT (factorización requerida igualmente) y para diferencias centradas con inversión local de bloques 6×6 baratos por nodo.

### Tests

- `tests/test_mass_lumping.py` (19 tests): contrato del helper HRZ, diagonalidad estricta por elemento, masa total preservada, entradas positivas, invariancia a la orientación del elemento (frames 2D), bloque-diagonalidad en Frame3D oblicuo, recuperación de la fundamental modal con masa lumped dentro de tolerancia.
- `tests/test_modal.py::TestModalSolverContract::test_lumped_runs_after_phase2`: cableado modal con `lumping="lumped"`.

## Consecuencias

**Positivas**

- Solidum gana análisis modal sin tocar arquitectura estática. Toda la maquinaria de ensamblaje, Dirichlet, expansión, y capa algebraica se reutiliza.
- M consistente con contrato extensible deja el camino libre a lumped sin breaking changes.
- La hoja de ruta queda registrada: las fases 3–7 saben dónde colgar sin reabrir debates de registry, nomenclatura o ubicación del estado.
- Cierra el lazo con ADR 0003 (capacidad eigen prevista), ADR 0007 (convergencia para Newton dentro de Newmark), ADR 0008 (density disponible).

**Negativas / costes**

- Cada elemento existente necesita `compute_mass_matrix` implementado. Esfuerzo localizado, ~30–60 líneas por elemento, mecánico para los lineales (truss, frame) y rutinario para Solid2D (la cuadratura ya existe).
- Cuando se introduzca el estado transitorio en fase 3, habrá que añadir `SolveResult` específico para transitorio. Acotado por la decisión 5 de este ADR (no toca `Node`).

**Alternativas consideradas**

- *Implementar lumped y consistente a la vez en fase 1*. Rechazada entonces: lumped no aportaba valor en modal de barras y bloqueaba validación contra fórmulas analíticas. La fase 2 (2026-05-18) materializó el lumped cuando se convirtió en precondición para diferencias centradas — sin reabrir debates de firma porque el ADR ya lo había previsto en `compute_mass_matrix(self, lumping="consistent")`.
- *`IntegratorRegistry` separado de `SolverRegistry`*. Rechazada: fragmenta el catálogo sin beneficio funcional. La distinción solver/integrator es semántica conversacional, no estructural. Los solvers actuales ya orquestan ambos roles.
- *Almacenar `velocity` y `acceleration` en `Node`*. Rechazada: contamina la clase con campos que solo usan análisis transitorios (~5–10% del uso esperado del catálogo). El patrón "estado dinámico en el resultado del solver" preserva la pureza de `Node` y replica el modo en que arc-length almacena `lambda_history` sin tocar `Node`.
- *Densidad como propiedad del elemento (no del material)*. Ya rechazada en ADR 0008. Se mantiene por completitud.
- *Empezar por transitorio lineal Newmark*. Rechazada como punto de entrada: introduciría M, C e integrador simultáneamente, sin validación analítica intermedia que aísle la fuente de un eventual fallo. Modal primero es la red de seguridad.

## Referencias

- Bathe, K.-J. (2014), *Finite Element Procedures*, cap. 9 (matriz de masa consistente y lumped), cap. 11 (autovalores generalizados, métodos de Lanczos y subspace iteration).
- Cook, Malkus & Plesha, *Concepts and Applications of Finite Element Analysis*, §11.3–§11.7.
- Hughes, T. J. R., *The Finite Element Method*, cap. 7 (algoritmos de integración temporal: Newmark, HHT-α, generalized-α).
- Crisfield, M. A., *Non-Linear Finite Element Analysis of Solids and Structures*, vol. 2 cap. 24 (dinámica no lineal).
- Wilson, E. L., *Three-Dimensional Static and Dynamic Analysis of Structures*, cap. 10–14.
- ADR 0003 — Capa algebraica (eigen como capacidad fase 3 prevista).
- ADR 0007 — Convergencia y tolerancias (aplicables a Newton dentro de Newmark).
- ADR 0008 — `Material.density` (única magnitud física nueva consumida).
