# Arquitectura — mapa navegable

> Documento vivo. Función: índice mental del sistema. Una página, sin profundizar en implementación. Para el detalle, ir al código.

## 1. Capas

```
┌───────────────────────────────────────────────────────────────┐
│  Entrada del usuario                                          │
│  • case.yaml         (descripción declarativa del problema)   │
│  • mallas .msh       (geometría desde Gmsh)                   │
└───────────────────────────────┬───────────────────────────────┘
                                │
┌───────────────────────────────▼───────────────────────────────┐
│  CAPA DE INICIALIZACIÓN  —  fenix/__init__.py                 │
│  Al hacer `import fenix` se dispara, una sola vez:            │
│    autodiscover.initialize()                                  │
│      └─ recorre fenix/materials, fenix/elements, fenix/math   │
│         e importa cada módulo; los decoradores @*Registry     │
│         registran las clases en MaterialRegistry,             │
│         ElementRegistry, SolverRegistry.                      │
└───────────────────────────────┬───────────────────────────────┘
                                │
┌───────────────────────────────▼───────────────────────────────┐
│  CAPA DE PARSING  —  fenix/utils/yaml_parser.py               │
│  Lee el YAML y construye objetos consultando los Registry.    │
│  Generic: introspecciona kwargs del constructor → no necesita │
│  saber qué materiales/elementos/solvers existen.              │
└───────────────────────────────┬───────────────────────────────┘
                                │
┌───────────────────────────────▼───────────────────────────────┐
│  CAPA DE DOMINIO  —  fenix/core/                              │
│   Domain ── Node ── DOF (numeración global)                   │
│      │                                                         │
│      └── Element  (base abstracta, contratos declarativos:    │
│              DOF_NAMES, STRAIN_DIM, N_INTEGRATION_POINTS)     │
│              │                                                 │
│              └── ElementState (trial/commit, stresses,        │
│                                state_vars, history)           │
│      └── Material (base abstracta, STRAIN_DIM,                │
│              PRIMARY_STATE_VAR)                               │
└───────────────────────────────┬───────────────────────────────┘
                                │
┌───────────────────────────────▼───────────────────────────────┐
│  CAPA DE MATEMÁTICA  —  fenix/math/                           │
│   assembly.py     →  ensamblaje sparse con cache COO          │
│   integration.py  →  cuadraturas de Gauss                     │
│   solvers.py      →  LinearSolver / NonlinearSolver /         │
│                      ArcLengthSolver  (registrados)           │
│   linalg/         →  capa algebraica K·x = b (ADR 0003)       │
│       base.py        →  Protocol + StiffnessProperties        │
│       lu.py          →  LUSolver (SuperLU, fallback universal)│
│       cholesky.py    →  CholeskySolver (CHOLMOD, opcional)    │
│       ldlt.py        →  LDLTSolver (placeholder fase 2)       │
│       dispatcher.py  →  select_solver(props, override)        │
└───────────────────────────────┬───────────────────────────────┘
                                │
┌───────────────────────────────▼───────────────────────────────┐
│  CAPA DE SALIDA  —  fenix/utils/vtk_exporter.py               │
│  Exporta U, σ, y la PRIMARY_STATE_VAR del material a VTK.     │
└───────────────────────────────────────────────────────────────┘
```

## 2. Flujo de datos típico (caso no lineal)

1. **Usuario** lanza `python ejecutar_yaml.py case.yaml`.
2. **`import fenix`** dispara `autodiscover` → todos los registries quedan poblados.
3. **`YamlParser`** lee el YAML, instancia materiales, elementos, BCs, solver.
4. **`Domain`** numera DOFs globalmente recorriendo nodos y elementos.
5. **`Solver`** (en cada iteración Newton-Raphson o paso de arc-length):
   - Pide a cada elemento su `K_local` y `f_internal` (vía `ElementState` trial).
   - **`Assembler`** ensambla `K_global` sparse usando topología COO cacheada.
   - Aplica BCs Dirichlet por método de penalidad (vectorizado).
   - **`linalg.select_solver(props)`** elige el backend algebraico adecuado
     (Cholesky para SPD, LU general en el resto) y resuelve `K·δU = R`.
   - Calcula residuo, evalúa convergencia (criterio dual: desplazamientos + fuerza).
   - Si converge el paso → `commit_state()` en cada elemento (trial → committed).
6. **`VtkExporter`** escribe el resultado del paso.

## 3. Puntos de extensión (dónde añadir piezas nuevas)

| Quiero añadir… | Dónde vive | Cómo se registra | Hace falta tocar `__init__.py`, parser, etc. |
|---|---|---|---|
| Material nuevo | `fenix/materials/<snake>.py` | `@MaterialRegistry.register` | No |
| Elemento nuevo | `fenix/elements/<snake>.py` | `@ElementRegistry.register` | No |
| Solver nuevo | `fenix/math/solver_<snake>.py` | `@SolverRegistry.register` | No |

El skill `/fenix-new <kind> <Name>` (`.claude/skills/fenix-new/`) genera el esqueleto completo (archivo + decorador + test).

## 4. Mantenimiento de este documento

La IA actualiza este mapa cuando cambia la **topología** (nueva capa, nueva carpeta canónica, nueva responsabilidad transversal), no cuando cambia la implementación interna de un módulo. Si tras un refactor mediano este diagrama ya no refleja la realidad, es un bug del refactor.
