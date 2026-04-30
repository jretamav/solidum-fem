# ADR 0004 — Imposición de condiciones de Dirichlet por eliminación directa con `ConstraintSet` afín

**Fecha**: 2026-04-30
**Estado**: Aceptado (2026-04-30)

## Resumen ejecutivo

Fenix sustituye la imposición de condiciones de Dirichlet por penalización (hoy: `K_ii ← K_ii + 10^{15}`, hard-coded) por **eliminación directa** mediante una abstracción única, `ConstraintSet`, que captura toda restricción afín de la forma `u_s = g_s + Σ α_si · u_mi`. La penalización desaparece del código. La nueva capa cubre Dirichlet homogéneo, asentamientos prescritos, simetrías inclinadas y, por extensión natural en una fase posterior, los multipoint constraints (MPC: uniones rígidas, periodicidad). El solver no lineal y los solvers algebraicos no se enteran del cambio: reciben un sistema reducido `K_red · u_libre = F_red` y la solución completa se reconstruye por post-proceso. Las reacciones se calculan a partir del residuo `K · u − F`, sin depender de un parámetro de penalización.

## Contexto

### Estado actual

[fenix/math/assembly.py:115-132](../../fenix/math/assembly.py#L115-L132) implementa el método de penalización:

```
K  ←  K + P · I_bc           (P = 10^{15} en filas/columnas de DOFs prescritos)
R[bc]  ←  P · (u_pre − u_cur)
```

Funciona porque (a) los tres solvers actuales son directos (LU/Cholesky vía la capa algebraica del ADR 0003), (b) el cálculo de reacciones en [fenix/results.py:155-161](../../fenix/results.py#L155-L161) usa la vía limpia `R = F_int − F_applied` y no `P·(u_pre−u_cur)`, y (c) los tests existentes son casi todos empotramientos puros (`u_pre = 0`), donde el sesgo del método no se manifiesta.

### Problemas con la penalización tal como está

1. **Penalización no escalada con la rigidez**. El valor `P = 10^{15}` es independiente de `max(diag(K))`. Para `EA/L` típico de 10⁶–10¹⁰, `κ(K)` resultante ronda 10⁹–10¹³. Sobrevive con LU/Cholesky directos; rompe cualquier solver iterativo (CG, GMRES) que llegará en la fase 5 del ADR 0003.
2. **Imposición sólo aproximada**. Para BC no homogéneas (`u_pre ≠ 0`) la solución cumple la restricción con error relativo ≈ ‖f‖/P. El sesgo no se ve en los tests actuales porque casi no hay casos de asentamiento prescrito.
3. **Conceptualmente errónea**. Una BC de Dirichlet es una restricción exacta, no un muelle de rigidez muy grande.
4. **No extiende a MPC limpiamente**. Penalizar acoplamientos (`u_s − Σ α_i u_mi = 0`) compone los problemas de condicionamiento.
5. **Acoplamiento con el solver**. El parámetro `penalty` viaja por la firma de `apply_bcs_to_system` y por los tres solvers de [fenix/math/solvers.py](../../fenix/math/solvers.py) — plumbing innecesario.

### Estado del arte y filosofía elegida

Los códigos FEM de referencia (ANSYS, Abaqus, Code_Aster, OpenSees, FEniCS, deal.II) usan **eliminación directa** como mecanismo por defecto para Dirichlet y MPC lineales, y reservan multiplicadores de Lagrange (o Lagrangiano aumentado) para restricciones no eliminables como contacto unilateral. El modelo arquitectural más limpio es el de deal.II con su clase `AffineConstraints`: un único objeto que captura Dirichlet, MPC, periodicidad y hanging nodes y los aplica al sistema con la misma maquinaria.

Para el dominio de Fenix — mecánica de sólidos con mallas conformes, sin contacto inminente, sin métodos inmersos — la trayectoria natural es:

- **Eliminación directa por defecto**, cubriendo Dirichlet y MPC lineales bajo una sola abstracción.
- **Lagrange / Lagrange aumentado diferido** al ADR que abra contacto.
- **Nitsche** queda fuera de alcance: sólo aporta cuando entran mallas no conformes o métodos inmersos, ninguno previsto.

## Decisión

### 1. Una abstracción única — `ConstraintSet` afín

Toda restricción lineal sobre los DOFs se expresa como una relación afín:

```
u_s  =  g_s  +  Σ α_si · u_mi
```

donde `s` es el DOF *esclavo* (eliminado del sistema reducido), `m_i` los DOFs *maestro* (sobreviven), `α_si` los coeficientes y `g_s` el término independiente.

Casos cubiertos por esta forma:

| Tipo de restricción | DOFs maestro | `g_s` |
|---|---|---|
| Empotramiento, apoyo en eje | (ninguno) | 0 |
| Asentamiento prescrito | (ninguno) | δ |
| Periodicidad `u_izq = u_der` | un maestro, α=1 | 0 |
| Simetría inclinada (apoyo en plano oblicuo) | combinación lineal de DOFs globales | 0 |
| Unión rígida master-slave (`u_s = u_m + r×θ_m`) | varios maestros, coeficientes geométricos | 0 |

**No incluido en `ConstraintSet`**: releases de momento en marcos (es ausencia de DOF, no relación entre DOFs — se gestiona en el elemento), contacto unilateral (es desigualdad, no igualdad — entrará por Lagrange aumentado en su ADR).

### 2. Forma matricial unificada

Con `n` DOFs totales y `n_f` libres tras eliminar esclavos, se construye una matriz sparse `T ∈ ℝ^{n × n_f}` (una columna por DOF libre) y un vector `g ∈ ℝ^n` (no nulo solo en filas esclavas con término independiente) tales que:

```
u  =  T · u_libre  +  g
```

El sistema reducido que ve el solver es:

```
K_red  =  Tᵀ · K · T
F_red  =  Tᵀ · ( F − K · g )
```

`K_red` es simétrica si `K` lo es (porque `T` es real), y positiva definida si `K` lo es y las restricciones son no redundantes.

Para Dirichlet pura, `T` es la matriz identidad rectangular que selecciona DOFs libres, y la operación `Tᵀ K T` se implementa como *fancy indexing* sparse `K[free, :][:, free]` — sin coste algebraico real.

Para MPC, `T` lleva los coeficientes `α_si` en las filas esclavas. Sigue siendo sparse y la composición es un par de productos sparse–sparse de coste despreciable comparado con la factorización.

### 3. Reacciones por residuo

Las reacciones se calculan post-hoc:

```
R_full  =  K_full · u_full  −  F_full
```

evaluado solo en los DOFs prescritos. No depende de ningún parámetro de penalización. Sustituye la rama actual de [fenix/results.py:155-161](../../fenix/results.py#L155-L161). Es exacta a redondeo, simétrica para sistemas simétricos, y transparente para MPC futuros (la "reacción" en un nodo maestro de un MPC es la suma ponderada de las contribuciones esclavas — caso que se documentará cuando entren MPC).

### 4. Estructura del módulo

Módulo nuevo `fenix.bc.constraints`:

```
fenix/bc/
    __init__.py
    constraints.py     # ConstraintSet, AffineConstraint, build_T_and_g(...)
```

Interfaz mínima:

```python
@dataclass(frozen=True)
class AffineConstraint:
    slave_dof: int
    masters: tuple[int, ...]      # vacío para Dirichlet pura
    coefficients: tuple[float, ...]
    g: float                      # término independiente

class ConstraintSet:
    def add_dirichlet(self, slave_dof: int, value: float = 0.0): ...
    def add_linear(self,
                   slave_dof: int,
                   masters: Sequence[int],
                   coefficients: Sequence[float],
                   g: float = 0.0): ...
    def build(self, ndof: int) -> tuple[sp.csr_matrix, np.ndarray]:
        """Devuelve (T, g) tras cierre transitivo y validación."""
```

`Assembler` expone:

```python
def reduce(self, K: sp.csr_matrix, F: np.ndarray
           ) -> tuple[sp.csr_matrix, np.ndarray, sp.csr_matrix, np.ndarray]:
    """Devuelve (K_red, F_red, T, g) listos para el solver."""
```

### 5. Validación temprana en `ConstraintSet.build`

Coherente con Reglas.md §1 (validación en construcción, no en runtime):

- **DOF prescrito dos veces** con valores distintos → error.
- **Ciclo master-slave** (`s ← m`, `m ← s`) → error.
- **Cadena master-slave** (`s ← m`, `m ← m'`) → cierre transitivo automático para que el sistema reducido contenga solo DOFs verdaderamente libres, replicando el algoritmo `closure` de `AffineConstraints` de deal.II.
- **Restricción redundante** (dos restricciones sobre el mismo esclavo, consistentes) → warning + deduplicación.
- **Restricción inconsistente** (igualdades incompatibles) → error.

Para Dirichlet pura los tres últimos casos son triviales o no aparecen; la validación entra en juego cuando se activen MPC.

### 6. Quién aplica `T` al pasar de `u_libre` a `u_full`

**Decisión**: el solver devuelve `u_libre`; un componente externo (el `Assembler`, vía un método `expand(u_libre)`) compone `u = T · u_libre + g` antes de poblar `SolveResult`.

Razón: mantiene los solvers en su rol de pura álgebra lineal, sin conocer `ConstraintSet`. La capa algebraica del ADR 0003 sigue recibiendo un sistema más pequeño y mejor condicionado, sin saber por qué.

### 7. Lo que se elimina del código existente

- `Assembler._build_bc_arrays`, `apply_bcs_to_system`, `apply_dirichlet_bcs`, `_bc_dofs`, `_bc_vals`, `_bc_built` en [fenix/math/assembly.py](../../fenix/math/assembly.py).
- Parámetro `penalty` en `LinearSolver`, `NonlinearSolver`, `ArcLengthSolver` ([fenix/math/solvers.py](../../fenix/math/solvers.py)).
- Rama de cálculo de reacciones penalty-aware en [fenix/results.py:144-161](../../fenix/results.py#L144-L161).

### 8. Lo que se preserva

- `Node.boundary_conditions` como API de entrada del usuario. Es la fuente declarativa; cambia el mecanismo, no la forma de declararlas.
- API pública de resultados (`SolveResult`, `internal_forces`).
- Capa algebraica `fenix.math.linalg`. Se beneficia gratis: recibe sistemas más pequeños, simétricos y mejor condicionados.

## Hoja de ruta — qué tipos de restricción habilita esta arquitectura

| Restricción | Cubierta por | Fase |
|---|---|---|
| Empotramiento, apoyo en eje (hoy) | `add_dirichlet(dof, 0)` | 1 |
| Asentamiento prescrito (hoy parcial, blindado por test nuevo) | `add_dirichlet(dof, δ)` | 1 |
| Simetría en plano alineado con ejes globales (hoy) | `add_dirichlet(dof_normal, 0)` | 1 |
| Simetría en plano oblicuo | `add_linear(slave, [m1, m2], [c1, c2], 0)` | 2 |
| Periodicidad celda unitaria | `add_linear` con un maestro y α=1 | 2 |
| Unión rígida master-slave (rigid link) | `add_linear` con coeficientes geométricos | 2 |
| Release de momento en marco | (no es `ConstraintSet` — es ausencia de DOF en el elemento) | independiente |
| Contacto unilateral | Lagrange aumentado, **ADR futuro** | — |
| Mallas no conformes (mortar) | Lagrange dualizado, **fuera de alcance previsible** | — |

## Fases de implementación

### Fase 1 — Eliminación directa para Dirichlet (esta es la orden de trabajo si apruebas el ADR)

1. Crear `fenix/bc/constraints.py` con `AffineConstraint`, `ConstraintSet`, `build(ndof) → (T, g)`. En fase 1 solo `add_dirichlet` está expuesto.
2. Añadir `Assembler.reduce(K, F) → (K_red, F_red, T, g)` y `Assembler.expand(u_libre) → u_full`. Construye internamente el `ConstraintSet` a partir de `Node.boundary_conditions`.
3. Migrar los tres solvers de [fenix/math/solvers.py](../../fenix/math/solvers.py) a llamar a `assembler.reduce(...)` en lugar de `apply_bcs_to_system(...)`. Los solvers retornan `u_libre`; el `Assembler` (o el wrapper del solver) llama `expand` antes de poblar `SolveResult`.
4. Sustituir cálculo de reacciones en [fenix/results.py](../../fenix/results.py) por `R = K_full · u_full − F_full` evaluado en `bc_dofs`.
5. Borrar todo el código de penalización: `_build_bc_arrays`, `apply_bcs_to_system`, `apply_dirichlet_bcs`, parámetros `penalty` en los solvers, rama penalty-aware en `results.py`.
6. Tests:
   - Todos los benchmarks existentes deben pasar sin modificación con tolerancias **más estrictas** que las anteriores (la solución es exacta a redondeo, antes era exacta a `‖f‖/P`).
   - Test nuevo: viga simplemente apoyada con asentamiento prescrito en un apoyo. Validación contra solución analítica de Bernoulli-Euler.
   - Test nuevo: comprobar simetría numérica de `K_red` (`‖K_red − K_redᵀ‖∞ ≤ tol·‖K_red‖∞`).
   - Test nuevo: reacciones suman a la fuerza externa en cada eje (equilibrio global) a tolerancia de redondeo.

### Fase 2 — MPC lineales

7. Exponer `ConstraintSet.add_linear(...)` y el cierre transitivo de cadenas master-slave.
8. Validaciones de la sección 5 activas para MPC (ciclos, redundancia, inconsistencia).
9. Tests:
   - Periodicidad en celda unitaria 2D bajo carga uniaxial: comparar con solución de homogeneización analítica.
   - Apoyo en plano oblicuo (45°) en armadura: comparar con la misma armadura rotada y empotrada en eje global.
   - Unión rígida master-slave entre dos nodos: comprobar que la cinemática se preserva exactamente.
10. Documentar la interpretación de reacciones en nodos maestro de un MPC.

### Fase 3 — Lagrange aumentado para contacto

ADR independiente cuando entre el primer caso de contacto. `ConstraintSet` se mantiene como primera capa (eliminación); las restricciones unilaterales viven en una segunda capa que se compone con la primera al ensamblar el sistema final.

## Consecuencias

**Positivas**

- Imposición exacta a redondeo, sin parámetro de penalización.
- κ(K_red) ≤ κ(K) — habilita iterativos cuando llegue la fase 5 del ADR 0003 sin tocar nada.
- Reglas.md §1 cumplido: añadir un nuevo *tipo* de restricción es un método nuevo en `ConstraintSet`. No toca solver, no toca capa algebraica, no toca API pública.
- API pública (`Node.boundary_conditions`, `SolveResult`, `internal_forces`) inalterada — refactor transparente para el usuario.
- Reacciones limpias por construcción.
- Reduce código: se borra más de lo que se añade en fase 1.

**Negativas / costes**

- Trabajo de migración localizado pero no trivial: tres solvers, módulo nuevo, cuatro o cinco tests nuevos.
- Recuperación post-solver (`expand`) es un paso adicional, aunque despreciable en coste.
- La detección y cierre de cadenas master-slave (fase 2) requiere implementar el algoritmo de cierre — no es complejo pero hay que testearlo bien.

**Alternativas consideradas**

- *Mantener penalización con escalado automático* `P = α·max(diag(K))`: corrige el problema de magnitud pero deja la imposición aproximada y no escala a iterativos. Rechazado: pequeña mejora a coste arquitectural similar a la solución correcta.
- *Multiplicadores de Lagrange como mecanismo único*: cubre todo (Dirichlet, MPC, contacto) bajo una sola maquinaria pero pierde positividad y simetría del sistema; obligaría a usar LDLᵀ con pivoteo o GMRES siempre, incluso para problemas SPD donde Cholesky es 2× más rápido. Rechazado para BC nodales puras y MPC lineales; reservado para contacto.
- *Eliminación directa + penalización opcional como fallback*: dos caminos paralelos sin justificación física. Rechazado por coste arquitectural inútil.
- *Nitsche*: irrelevante mientras Fenix use mallas conformes. Diferido sin fecha.

**Implicaciones para la documentación**

- `arquitectura.md`: nuevo subsistema `fenix.bc` y diagrama de flujo `Domain → Assembler.reduce → Solver → Assembler.expand → SolveResult`.
- `conceptos_clave.md`: entrada nueva *"Imposición de condiciones de frontera por eliminación directa"* explicando el modelo afín `u = T·u_libre + g` en lenguaje conceptual.
- Manual de arquitectura (FF-MA): capítulo o sección dedicada a la capa de restricciones, con la tabla de tipos cubiertos por fase.
- Manual de usuario: la declaración de BC en YAML/script no cambia. Mención breve de que la imposición es exacta y sin parámetros de tuning.

## Criterio de cierre del ADR

Se cierra como `Aceptado` cuando el usuario valida el resumen ejecutivo y la hoja de ruta. La fase 1 se ejecuta inmediatamente tras la aceptación bajo este ADR. Las fases 2 y 3 se ejecutan cuando entren los casos de uso que las consumen (primer MPC, primer contacto).
