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
│  CAPA DE INICIALIZACIÓN  —  solidum/__init__.py                 │
│  Al hacer `import solidum` se dispara, una sola vez:            │
│    autodiscover.initialize()                                  │
│      └─ recorre solidum/materials, solidum/elements, solidum/math   │
│         e importa cada módulo; los decoradores @*Registry     │
│         registran las clases en MaterialRegistry,             │
│         ElementRegistry, SolverRegistry.                      │
└───────────────────────────────┬───────────────────────────────┘
                                │
┌───────────────────────────────▼───────────────────────────────┐
│  CAPA DE PARSING  —  solidum/utils/yaml_parser.py               │
│  Lee el YAML y construye objetos consultando los Registry.    │
│  Generic: introspecciona kwargs del constructor → no necesita │
│  saber qué materiales/elementos/solvers existen.              │
└───────────────────────────────┬───────────────────────────────┘
                                │
┌───────────────────────────────▼───────────────────────────────┐
│  CAPA DE DOMINIO  —  solidum/core/                              │
│   Domain ── Node ── DOF (numeración global)                   │
│      │                                                         │
│      └── Element  (base abstracta, contratos declarativos:    │
│              DOF_NAMES, STRAIN_DIM, N_INTEGRATION_POINTS)     │
│              │                                                 │
│              └── ElementState (trial/commit, stresses,        │
│                                state_vars, history)           │
│      └── Material (base abstracta, STRAIN_DIM,                │
│              PRIMARY_STATE_VAR, STATE_SCHEMA)                 │
│      └── familias paralelas: CohesiveMaterial (ADR 0010),     │
│              ThermalMaterial (FLUX_DIM, Etapa 8)              │
└───────────────────────────────┬───────────────────────────────┘
                                │
┌───────────────────────────────▼───────────────────────────────┐
│  CAPA DE MATEMÁTICA  —  solidum/math/                         │
│   assembly.py     →  ensamblaje sparse con cache COO; modos   │
│                      de cuerpo rígido (rigid_basis)           │
│   batch/          →  ensamblaje por lotes: familias + un      │
│                      kernel Numba (serie / prange) ADR 0014   │
│   integration.py  →  cuadraturas de Gauss                     │
│   solvers/        →  un módulo por solver (13 registrados)    │
│       corrector.py   →  NewtonCorrector compartido (ADR 0015) │
│       model_checks.py→  red de seguridad estática (ADR 0019)  │
│       diagnostics.py →  errores tipados (ADR 0011, 0019)      │
│   linalg/         →  capa algebraica K·x = b (ADR 0003)       │
│       base.py        →  Protocol + StiffnessProperties        │
│       lu.py          →  LUSolver (SuperLU, siempre)           │
│       cholesky.py    →  CholeskySolver (CHOLMOD, opcional)    │
│       pardiso.py     →  PardisoSolver (MKL, opcional; 0017)   │
│       iterative.py   →  CG/MINRES + AMG, a petición (0018)    │
│       nullspace.py   →  modos de cuerpo rígido por DOF_NAMES  │
│       ldlt.py        →  LDLTSolver (placeholder, deuda #7)    │
│       dispatcher.py  →  select_solver(props, override)        │
└───────────────────────────────┬───────────────────────────────┘
                                │
┌───────────────────────────────▼───────────────────────────────┐
│  CAPA DE SALIDA  —  solidum/utils/vtk_exporter.py             │
│  Exporta U, σ, T, flujo y la PRIMARY_STATE_VAR a VTK.         │
└───────────────────────────────────────────────────────────────┘
```

## 2. Flujo de datos típico (caso no lineal)

1. **Usuario** lanza `python ejecutar_yaml.py case.yaml`.
2. **`import solidum`** dispara `autodiscover` → todos los registries quedan poblados.
3. **`YamlParser`** lee el YAML, instancia materiales, elementos, BCs, solver.
4. **`Domain`** numera DOFs globalmente recorriendo nodos y elementos.
5. **`Solver`** (en cada iteración Newton-Raphson o paso de arc-length):
   - Pide a cada elemento su `K_local` y `f_internal` (vía `ElementState` trial).
   - **`Assembler`** ensambla `K_global` sparse usando topología COO cacheada.
   - **`Assembler.reduce(K, F, ...)`** elimina los DOFs prescritos del sistema
     (ADR 0004): construye `T, g` tal que `u = T·u_libre + g` y devuelve
     `K_red = TᵀKT`, `F_red = Tᵀ(F − K·g)`. Imposición exacta, sin penalización.
   - **`linalg.select_solver(props)`** elige el backend algebraico adecuado
     (Cholesky para SPD, LU general en el resto) y resuelve `K_red·δU_red = R_red`.
   - **`Assembler.expand(δU_red, …)`** reconstruye el incremento completo.
   - Calcula residuo, evalúa convergencia (criterio dual: desplazamientos + fuerza).
   - Si converge el paso → `commit_state()` en cada elemento (trial → committed).
6. **`VtkExporter`** escribe el resultado del paso.

## 3. Puntos de extensión (dónde añadir piezas nuevas)

| Quiero añadir… | Dónde vive | Cómo se registra | Hace falta tocar `__init__.py`, parser, etc. |
|---|---|---|---|
| Material nuevo | `solidum/materials/<snake>.py` | `@MaterialRegistry.register` | No |
| Elemento nuevo | `solidum/elements/<snake>.py` | `@ElementRegistry.register` | No |
| Backend algebraico nuevo | `solidum/math/linalg/<snake>.py` | entrada en `_REGISTRY` del despachador (import opcional absorbido si falta la dependencia) | Sólo `dispatcher.py` y `linalg/__init__.py` |
| Solver nuevo | `solidum/math/solvers/<snake>.py` | `@SolverRegistry.register` | No |

El skill `/solidum-new <kind> <Name>` (`.claude/skills/solidum-new/`) genera el esqueleto completo (archivo + decorador + test).

## 4. Mantenimiento de este documento

La IA actualiza este mapa cuando cambia la **topología** (nueva capa, nueva carpeta canónica, nueva responsabilidad transversal), no cuando cambia la implementación interna de un módulo. Si tras un refactor mediano este diagrama ya no refleja la realidad, es un bug del refactor.
