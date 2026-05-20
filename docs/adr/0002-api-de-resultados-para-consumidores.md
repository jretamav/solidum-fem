# ADR 0002 — API pública de resultados para consumidores externos

**Fecha**: 2026-04-22
**Estado**: Aceptado (2026-04-22)

## Contexto

Solidum FEM se empieza a consumir desde una GUI externa (FenixBAR, análisis estructural con elementos barra). FenixBAR necesita, tras una solución, información suficiente para visualizar posproceso estándar de CAE: deformada escalada, reacciones en apoyos y diagramas de esfuerzos internos (N, V, M, y en 3D también T, My, Mz) sobre cada barra.

Hoy esa información no está expuesta de forma uniforme:

- El único "runner" es el script de ejemplo `examples/ejecutar_yaml.py`; no hay entrypoint oficial.
- `VtkExporter` solo maneja `Quad4`, `Tri3`, `Truss2D`, `Truss3D` (estos últimos como celdas `line` sin esfuerzos). `Frame2DEuler`, `Frame2DEulerCorot`, `Frame2DTimoshenko`, `Frame3D`, `Cable2DCorot`, `Cable3DCorot` no se exportan. `displacements` solo escribe `ux`/`uy` (sin `uz` ni rotaciones). _(Cerrado 2026-05-04: el exportador acepta cualquier elemento de 2 nodos como celda `line`, escribe desplazamientos 3D y rotaciones nodales cuando hay DOFs rotacionales; los esfuerzos internos N/V/M de barras siguen consumiéndose vía `SolveResult.element_forces`, no por VTK.)_
- Cada tipo de elemento conoce su formulación cinemática y constitutiva (Euler vs Timoshenko, corotacional vs lineal, ejes locales en 3D, tracción-sólo en cables), pero no ofrece método público homogéneo para devolver esfuerzos internos dado un `U`.

Alternativa descartada: que el consumidor calcule N/V/M desde `U`. Duplicaría la formulación fuera del repo que la define, con riesgo de divergencia en cada elemento nuevo. Contrario al principio de que la lógica FEM vive en Solidum FEM.

## Decisión

API pública de resultados en tres niveles.

### 1. Contrato por elemento — `internal_forces`

Todos los elementos barra/viga exponen:

```python
def internal_forces(self, U: np.ndarray) -> ElementForces:
    """Esfuerzos internos en ejes locales, en los dos nodos i, j.

    Respeta las convenciones de signos de Reglas.md §5.
    Para corotacionales, usa el state ya comprometido tras solve().
    """
```

Retorno homogéneo por familia, con claves fijas y valor `np.ndarray` de shape `(2,)` (índice 0 = nodo i, índice 1 = nodo j). Componentes no aplicables se omiten del dict (no se rellenan con ceros — su ausencia es información).

```python
# Definido en solidum.results
@dataclass(frozen=True)
class ElementForces:
    kind: Literal["truss", "cable", "frame2d", "frame3d"]
    components: dict[str, np.ndarray]  # shape (2,) por clave

# Claves por familia:
#   truss, cable   -> {"N"}
#   frame2d        -> {"N", "V", "M"}           # V en convención viga (§5)
#   frame3d        -> {"N", "Vy", "Vz", "T", "My", "Mz"}  # RHR (§5)
```

**Por qué homogéneo con `ndarray(2,)`**: el consumidor itera siempre igual (`for key, arr in ef.components.items(): arr[0], arr[1]`), sin ramificar por tipo más allá de qué claves existen. Evita el mix `float`/tupla/dict de la versión anterior.

**Convenciones de signo**: son contrato público. Definidas en Reglas.md §5, no se re-definen aquí ni por elemento. Resumen operativo:
- 2D: convención de viga estructural (`N` tracción+, `V` rota diferencial en horario, `M` sagging+).
- 3D: stress-resultant / RHR pura.
- `M_2D` ≡ `Mz_3D`; `V_2D` y `Vy_3D` difieren en signo en el API público (deliberado, ver §5).

**Muestreo interno para diagramas (fase 2, no bloqueante)**: los elementos frame admitirán opcionalmente `sample_internal(U, n_points) -> dict[str, np.ndarray]` que evalúa `N(x), V(x), M(x)` usando las funciones de forma del propio elemento (cúbico en Euler, lineal a trozos en Timoshenko reducido, etc.). Evita que la GUI asuma la forma del diagrama. Se deja fuera del alcance inicial; la API se diseña para admitirlo sin ruptura.

### 2. Agregado en `Domain` — `SolveResult`

Calculado **eager** al final de `solver.solve()` (no lazy). Inmutable. Accesible como `domain.last_result`.

```python
@dataclass(frozen=True)
class SolveResult:
    U: np.ndarray                               # desplazamientos globales, shape (n_dof,)
    F_applied: np.ndarray                       # cargas aplicadas por el usuario (nodales + equivalentes de elemento ya ensambladas), shape (n_dof,)
    R: np.ndarray                               # reacciones globales, shape (n_dof,); ceros en DOF libres
    reactions_by_node: dict[int, dict[str, float]]  # vista conveniente: node_id -> {"ux": .., "uy": ..}
    element_forces: dict[int, ElementForces]    # elem_id -> resultado de internal_forces()
    converged: bool
    num_steps: int                              # 1 para solver lineal; N para Newton-Raphson/arc-length
```

Decisiones:
- `F_applied` **incluye** cargas equivalentes de elemento (distribuidas, térmicas) ya ensambladas a nodos — es el vector que entra al sistema. No incluye reacciones. Esto corresponde a "lo que el sistema lineal ve como F" en `K·U = F_applied`. Si el usuario quiere "lo que puse en el YAML", lo tiene en `domain` (entrada).
- `R` se da como vector global con ceros en libres, y `reactions_by_node` como vista conveniente filtrada a nodos restringidos. Ambos redundantes a propósito: barato de dar ambos, GUIs grandes prefieren vector, scripts prefieren dict.
- `frozen=True` para evitar que un consumidor modifique resultados y sorprenda a otro consumidor posterior.

### 3. Entrypoints públicos — `solidum.run` y `solidum.run_yaml`

Dos funciones, no una con unión de tipos:

```python
def run(domain: Domain, *, solver_config: SolverConfig | None = None) -> SolveResult: ...
def run_yaml(path: str | Path, *, solver_config: SolverConfig | None = None) -> SolveResult: ...
```

Encapsulan `Assembler → solver.solve → recolección de resultados`. `run_yaml` delega en `run` tras parsear.

El script `examples/ejecutar_yaml.py` pasa a ser una demo delgada (≤ 20 líneas) que llama a `solidum.run_yaml`. La vía oficial para consumidores es `import solidum; result = solidum.run(domain)`.

## Consecuencias

**Positivas**
- FenixBAR (y futuros consumidores) tienen contrato estable que no depende de `.vtu` ni de scripts de ejemplo.
- La lógica por tipo de elemento queda encapsulada donde vive la formulación; añadir un elemento nuevo implica implementar su `internal_forces` y se integra sin tocar consumidores.
- `VtkExporter` sigue existiendo como capa de exportación a ParaView, pero deja de ser el cuello de botella para GUIs que trabajan en memoria.
- Convenciones de signo centralizadas en Reglas.md §5; auditables y versionables como parte del contrato.

**Negativas / costes**
- Implementar `internal_forces` en 10 tipos de elemento (8 barra + 2 sólidos si se cubren). Acotado pero no trivial en los corotacionales, que deben usar `state` comprometido.
- El elemento 2D tiene que aplicar traducción B→A (solo cambio de signo en `V`) dentro de `internal_forces()`. Un punto por tipo 2D, fácil de testear.
- Redundancia `R` + `reactions_by_node` es deliberada; documentarla para que no se perciba como inconsistencia.

**Alternativas consideradas**
- Extender `VtkExporter` a todos los tipos de barra: no resuelve el problema — GUIs interactivas no quieren pasar por disco ni por un formato de ParaView.
- Dejar al consumidor calcular N/V/M desde `U`: duplicaría formulación fuera del repo. Descartado.
- Convención stress-resultant pura también en API 2D: más limpio internamente pero rompe la intuición clásica de diagramas. Rechazado a favor de A interno + B expuesto en 2D.
- `solidum.run(domain_or_yaml)` con unión de tipos: rechazado por complicar el tipado y la documentación; dos funciones son más claras.

**Implicaciones para `arquitectura.md` y `conceptos_clave.md`**
- Añadir `SolveResult`, `ElementForces` y el contrato `Element.internal_forces` al mapa arquitectural.
- Entrada nueva en `conceptos_clave.md` sobre "API de resultados" que remita a este ADR y a Reglas.md §5.

## Plan de implementación (orden sugerido)

1. Definir `ElementForces` y `SolveResult` en `solidum/results.py`.
2. Añadir método abstracto `internal_forces` a `Element` base con default `NotImplementedError`.
3. Implementar en truss (trivial), cable (trivial con rama slack), frame2D Euler/Timoshenko, frame2D corot, frame3D.
4. Ensamblar `SolveResult` en `Domain` al final de `solver.solve`; cablear `last_result`.
5. Añadir `solidum.run` y `solidum.run_yaml`; refactorizar `ejecutar_yaml.py` como demo.
6. (Fase 2, cuando FenixBAR lo pida) `sample_internal` en frames para diagramas suaves.

Tests por punto (física no obvia, ver Reglas.md §6): ménsula con carga puntual para `Frame2DEuler` (M lineal, V constante), viga biempotrada con carga uniforme, torsión pura en Frame3D, arco con cable pretensado. Cada test compara signos contra Reglas.md §5, no contra la referencia original.
