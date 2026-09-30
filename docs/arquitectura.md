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
│  CAPA DE INICIALIZACIÓN  —  solidum/__init__.py               │
│  Al hacer `import solidum` se dispara, una sola vez:          │
│    autodiscover.initialize()                                  │
│      └─ recorre solidum/materials, solidum/elements y         │
│         solidum/math e importa cada módulo; los decoradores   │
│         @*Registry registran las clases en MaterialRegistry,  │
│         ThermalMaterialRegistry, ElementRegistry,             │
│         SolverRegistry.                                       │
│  NO recorre solidum/user/ (módulos de usuario, ADR 0020):     │
│  cada uno se carga sólo si el modelo lo pide —                │
│  `user_modules:` en el YAML o solidum.load_user_module().     │
└───────────────────────────────┬───────────────────────────────┘
                                │
┌───────────────────────────────▼───────────────────────────────┐
│  CAPA DE PARSING  —  solidum/utils/yaml_parser.py             │
│  Lee el YAML y construye objetos consultando los Registry.    │
│  Generic: introspecciona kwargs del constructor → no necesita │
│  saber qué materiales/elementos/solvers existen. Las familias │
│  de material (Registry.families(), cada una con su sección    │
│  YAML) se recorren en un solo bucle; las referencias de un    │
│  elemento a otra familia las declara Element.REFERENCE_KWARGS.│
│  Una sección de primer nivel desconocida es error.            │
└───────────────────────────────┬───────────────────────────────┘
                                │
┌───────────────────────────────▼───────────────────────────────┐
│  CAPA DE DOMINIO  —  solidum/core/                            │
│   Domain ── Node ── DOF (numeración global)                   │
│      │                                                        │
│      └── Element  (base abstracta, contratos declarativos:    │
│              DOF_NAMES, STRAIN_DIM, N_INTEGRATION_POINTS)     │
│              │                                                │
│              └── ElementState (trial/commit, stresses,        │
│                                state_vars, history)           │
│      └── Material (base abstracta, STRAIN_DIM,                │
│              PRIMARY_STATE_VAR, STATE_SCHEMA)                 │
│      └── familia paralela: ThermalMaterial (FLUX_DIM, Etapa 8)│
│              (un módulo de usuario puede declarar la suya,    │
│               como la cohesiva de discontinuities)            │
└───────────────────────────────┬───────────────────────────────┘
                                │
┌───────────────────────────────▼───────────────────────────────┐
│  CAPA DE MATEMÁTICA  —  solidum/math/                         │
│   assembly.py     →  ensamblaje sparse con cache COO; modos   │
│                      de cuerpo rígido (rigid_basis)           │
│   batch/          →  ensamblaje por lotes: familias + un      │
│                      kernel Numba (serie / prange) ADR 0014   │
│   integration.py  →  cuadraturas de Gauss                     │
│   solvers/        →  un módulo por solver (14 registrados)    │
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

**Programa principal y módulos de usuario** ([ADR 0020](adr/0020-modulos-de-usuario.md)). Todo `solidum/` salvo `solidum/user/` es el **programa principal**: Solidum estándar, elementos finitos clásicos. Una formulación no estándar —la que no cumple el criterio de admisión de [Reglas.md](../Reglas.md) §1— vive en un **módulo de usuario**, `solidum/user/<nombre>/`, como los elementos, materiales y macros de usuario de FEAP. Hoy hay dos: `discontinuities`, la discontinuidad interior embebida del ADR 0010, y `dissipation_control`, el arco por disipación de Gutiérrez (`DissipationArcLengthSolver`, desde el 2026-09-29). El programa principal no lo importa ni lo nombra; lo vigilan los tres controles de `tests/test_main_program_purity.py`: estructural (ninguna forma de importar un módulo, salvo el cargador genérico `load_user_module` / `available_user_modules`), léxico (ninguna mención de lo que el módulo registra, con nombres derivados al cargarlo) y en ejecución (`import solidum` no carga ninguno). La carga es atómica y sólo añade entradas a los registros: si un módulo reasigna una existente, `load_user_module` lanza `ValueError` y deshace la carga. El módulo se conecta al programa principal sólo por piezas públicas (lo comprueba `tests/test_user_module_public_imports.py`): los registros (y, si la necesita, una familia de material propia con su sección YAML), `Element.REFERENCE_KWARGS`, el gancho `Element.prepare_step`, la cinemática pública de `solidum.elements.solid_2d`, el problema por paso del arco `solidum.math.solvers.ArcProblem` y la ruta gmsh con elemento por grupo físico. Sus specs viven en `docs/user/<nombre>/specs/` y sus tests en `tests/user/<nombre>/`.

## 2. Flujo de datos típico (caso no lineal)

1. **Usuario** lanza `python ejecutar_yaml.py case.yaml`.
2. **`import solidum`** dispara `autodiscover` → los registries del programa principal quedan poblados. Ningún módulo de usuario se carga todavía.
3. **`YamlParser`** lee el YAML: primero carga los módulos de usuario que declare `user_modules`, después valida e instancia materiales (cada familia desde su sección), elementos, BCs, solver.
4. **`Domain`** numera DOFs globalmente recorriendo nodos y elementos.
5. **`Solver`** (en cada iteración Newton-Raphson o paso de arc-length):
   - Al empezar cada paso, los solvers que avanzan por pasos llaman a `Assembler.prepare_all_steps(U_committed)`, que invoca el gancho `Element.prepare_step` en cada elemento (no-op salvo en los que toman decisiones discretas entre pasos, como activar una grieta; ADR 0020, P5).
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

## 3. Dónde añadir piezas nuevas

| Quiero añadir… | Dónde vive | Cómo se registra | Hace falta tocar `__init__.py`, parser, etc. |
|---|---|---|---|
| Material nuevo | `solidum/materials/<snake>.py` | `@MaterialRegistry.register` (térmico: `@ThermalMaterialRegistry.register`) | No |
| Elemento nuevo | `solidum/elements/<snake>.py` | `@ElementRegistry.register` | No |
| Backend algebraico nuevo | `solidum/math/linalg/<snake>.py` | entrada en `_REGISTRY` del despachador (import opcional absorbido si falta la dependencia) | Sólo `dispatcher.py` y `linalg/__init__.py` |
| Solver nuevo | `solidum/math/solvers/<snake>.py` | `@SolverRegistry.register` | No |
| Componente no estándar (formulación de investigación) | módulo de usuario `solidum/user/<nombre>/<snake>.py` (ADR 0020) | los mismos registros; una familia de material propia, con `Registry` y su `YAML_SECTION` en el `registry.py` del módulo | No: ni el programa principal ni el lector YAML. El modelo lo carga con `user_modules` o `solidum.load_user_module` |

El skill `/solidum-new <kind> <Name>` (`.claude/skills/solidum-new/`) genera el esqueleto completo (archivo + decorador + test); con `--user <nombre>`, dentro de un módulo de usuario.

## 4. Mantenimiento de este documento

La IA actualiza este mapa cuando cambia la **topología** (nueva capa, nueva carpeta canónica, nueva responsabilidad transversal), no cuando cambia la implementación interna de un módulo. Si tras un refactor mediano este diagrama ya no refleja la realidad, es un bug del refactor.
