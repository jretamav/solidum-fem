---
name: solidum-new
description: Genera scaffolding para nuevos materiales, elementos o solvers de Solidum FEM (archivo + decorador de registro + test esqueleto), en el programa principal o dentro de un módulo de usuario (solidum/user/<nombre>/). Úsalo cuando el usuario pida añadir/crear un material, elemento o solver derivado de su investigación.
---

# solidum-new

Skill de scaffolding para Solidum FEM. Reduce el coste de incorporar nuevos componentes a la arquitectura: la física la fija el usuario (en la spec), esta skill se encarga del boilerplate, del registro y del test esqueleto.

**Reparto de decisiones** (Reglas.md §2 y §3, desde el 2026-09-29): la arquitectura y la programación las decide la IA por completo. Al usuario sólo se le pregunta lo que cambia una ecuación o un resultado físico. No le presentes decisiones de software (carpeta, registro, estado propio, piezas genéricas): tómalas, documéntalas y repórtalas en una frase.

## Cuándo invocar

Activar siempre que el usuario pida añadir un nuevo:
- **Material** (constitutivo): elasticidad alternativa, plasticidad, daño, viscoelasticidad, hiperelasticidad, etc.
- **Elemento**: nueva formulación 1D/2D/3D (sólido, viga, cáscara, contacto…).
- **Solver**: lineal, no-lineal, path-following, dinámico, etc.

Disparadores típicos: "añade un material …", "crea un elemento …", "implementa un solver …", "/solidum-new …".

## Destino: programa principal o módulo de usuario

Lo decide la IA aplicando el **criterio de admisión** de Reglas.md §1 (ADR 0020). Una formulación es **estándar** si cumple las dos condiciones:

1. aparece en los textos de referencia del área: Zienkiewicz-Taylor, Hughes, Bathe, Belytschko-Liu-Moran, de Souza Neto-Perić-Owen, de Borst-Crisfield-Remmers-Verhoosel, Chopra;
2. está implementada en al menos un código general establecido: Abaqus, ANSYS, DIANA, Code_Aster, OpenSees.

Estándar → **programa principal** (sin `--user`). Si no cumple las dos → **módulo de usuario** (`--user <modulo>`, ver *Dentro de un módulo de usuario*). Sólo si el criterio no decide, pregunta al usuario antes de generar nada, con una pregunta y tu recomendación justificada: es un juicio sobre la literatura del área (Reglas.md §3). La clasificación queda anotada con su justificación en la entrada del catálogo (Reglas.md §1) y en la spec: en un componente estándar, en el bloque `admission` del contrato (`textbooks` y `general_codes`, al menos una entrada en cada uno; Reglas.md §4); en uno de módulo, con la condición que no cumple.

Las entradas de `admission` se escriben con los **identificadores canónicos** de `solidum/tools/spec.py`, no con los nombres de arriba: `textbooks` ∈ `ADMISSION_TEXTBOOKS` (`zienkiewicz_taylor`, `hughes`, `bathe`, `belytschko_liu_moran`, `de_souza_neto`, `de_borst_crisfield`, `chopra`) y `general_codes` ∈ `ADMISSION_GENERAL_CODES` (`abaqus`, `ansys`, `diana`, `code_aster`, `opensees`). Cada entrada es el identificador solo o `identificador: "nota breve"` (capítulo o apartado del texto; comando, elemento u opción en el código). La plantilla de spec lo recuerda y el validador rechaza un nombre que no esté en la lista, enumerando los válidos.

Ejemplos según el criterio:

- **Estándar**: XFEM (Belytschko-Liu-Moran, 2.ª ed.; Abaqus); leyes cohesivas tracción-salto genéricas con elementos de interfaz (de Borst-Crisfield-Remmers-Verhoosel; elementos cohesivos de Abaqus).
- **No estándar**: la discontinuidad embebida con cinemática KOS de Retama (2010) (módulo `discontinuities`); el seguimiento de trayectoria de grieta; el control del paso por liberación de energía de Gutiérrez (2004) (módulo `dissipation_control`: ningún código general lo implementa de forma nativa).

Una pieza estándar que hoy sólo consume un módulo sigue en él. Es el caso de la familia de materiales cohesivos de `discontinuities`, cuyo único consumidor es el elemento embebido: mientras el programa principal no tenga un elemento de interfaz, una ley cohesiva nueva va a ese módulo. Cuando entre un elemento de interfaz estándar, la familia pasa con él al programa principal, sin física nueva, y el módulo la reutiliza.

## Argumentos

`/solidum-new <kind> <Name> [--user <modulo>]`

- `kind` ∈ {`material`, `element`, `solver`}
- `Name`: PascalCase (`HyperelasticOgden`, `Shell4`, `DynamicNewmark`).
- `--user <modulo>`: genera el componente dentro del módulo de usuario `solidum/user/<modulo>/` (lo crea si no existe). `<modulo>` en snake_case y sin `_` inicial (`available_user_modules()` ignora esos subpaquetes). Lo decides tú con el criterio de admisión (sección *Destino*), aunque el usuario no lo escriba.

Si faltan `kind` o `Name`, dedúcelos de la petición o de la spec; pregunta sólo si no se deducen.

## Flujo

0. **Clasifica** el componente con el criterio de admisión (sección *Destino*) y localiza su spec, que es la orden de trabajo (Reglas.md §4): `docs/specs/<Name>.md` en el programa principal, `docs/user/<modulo>/specs/<Name>.md` en un módulo. El esqueleto puede generarse antes, pero la física no se escribe hasta que la spec está completa.
1. **Recoge los datos físicos** de la spec y pregunta al usuario sólo lo que falte, y sólo física (no satures):
   - **material**: dimensión de la deformación (1, 3 ó 6 componentes en Voigt) · parámetros físicos del constructor, con unidades · variables internas con significado físico (κ, ω, εᵖ…), y cuál gobierna el estado si hay una principal.
   - **element**: grados de libertad por nodo · número de nodos · dimensión de la deformación · regla de integración · parámetros geométricos (`A`, `I`, `thickness`…).
   - **solver**: algoritmo (lineal, Newton, control de carga o de desplazamiento, longitud de arco, dinámico…) · parámetros con significado físico y sus valores por omisión (tolerancias, incrementos).

   Los atributos declarativos los deduce la IA de esos datos, sin preguntarlos: `STRAIN_DIM`, `DOF_NAMES`, `N_INTEGRATION_POINTS`, `GEOMETRICALLY_NONLINEAR` (elementos: `True` sólo si la cinemática es corrotacional o de grandes desplazamientos; la red de seguridad del ADR 0019 lo usa y el barrido de contrato lo contrasta con la rigidez tras un giro rígido), `PRIMARY_STATE_VAR` y `STATE_SCHEMA` (`{nombre: forma}` de las variables internas, `{}` sin historia — obligatorio, ADR 0014). Lo mismo con el carácter del elemento: si la formulación guarda historia fuera de los puntos de Gauss (un salto, una grieta, un contacto) es un **elemento con estado propio**, y si activa algo entre pasos toma **decisiones discretas** (ver Convenciones).
2. **Genera el archivo** en su carpeta canónica (ver Templates; con `--user`, en `solidum/user/<modulo>/`) usando `snake_case(Name)` como nombre de archivo.
3. **Genera el test esqueleto** en `tests/test_<snake>.py` (con `--user`, en `tests/user/<modulo>/test_<pref>_<snake>.py`) con un test mínimo de "se construye + compute_state/compute_element_state/solve devuelve shapes correctas".
4. **Verifica**: ejecuta `pytest tests/test_<snake>.py -x` (con `--user`: `pytest tests/user/<modulo>/ tests/test_main_program_purity.py tests/test_user_module_public_imports.py -x`) y confirma que el esqueleto pasa.
5. **Comunica al usuario**: la clasificación (estándar / módulo de usuario) con su justificación en una línea, los archivos creados, dónde queda el `# TODO: implementar la física aquí`, y el ejemplo YAML mínimo para usarlo (si aplica a material/element; con `--user`, incluye `user_modules: [<modulo>]`).

**No hace falta** tocar `__init__.py`, `registry.py` ni `yaml_parser.py` del programa principal: el decorador y el autodescubrimiento (o `load_user_module`, en un módulo de usuario) hacen el resto.

## Convenciones del proyecto

- Material: hereda de `solidum.core.material.Material`. Decora con `@MaterialRegistry.register`. Implementa `compute_state(strain, state_vars=None) -> (stress, tangent, new_state_vars)` y declara `STATE_SCHEMA` (el barrido de contratos comprueba que coincide con el estado devuelto). El kernel por lotes (`BATCH_KERNEL`, `batch_params`, `batch_matrix`) es opcional; ver ADR 0014.
- Element: hereda de `solidum.core.element.Element`. Decora con `@ElementRegistry.register`. Implementa `compute_element_state(u_e) -> (K_e, F_int_e)`. La base se encarga de `state_vars`, validación STRAIN_DIM, registro de DOFs, creación de `ElementState` y del `commit_state` **de lo que vive en `ElementState`** (variables internas y esfuerzos por punto de Gauss).
  - **Elemento con estado propio** (historia fuera de `ElementState`: un salto, una grieta, un contacto; ADR 0020, P8): la base no lo conoce. El elemento sobrescribe `commit_state` (llama a `super().commit_state()` y confirma su propio estado) y `compute_global_stiffness` (la evaluación auxiliar en `u = 0` no debe dejar huella en su trial propio: lo guarda y lo restaura al salir). No va por el camino por lotes: `element_is_batchable` lo excluye si la clase redefine `compute_element_state` o `commit_state` por debajo de la que declaró `BATCH_KINEMATICS`, así que **no declares `BATCH_KINEMATICS` en la propia clase con estado propio** (la protección dejaría de verla). Modelo: `solidum/user/discontinuities/embedded_cst.py`.
  - **Decisiones discretas entre pasos** (activar una grieta o un contacto): sobrescribe `prepare_step(U_committed)` en vez de decidir dentro del Newton; debe ser idempotente para el mismo `U_committed` (tras una bisección, los solvers estáticos que bisecan repiten la llamada con el mismo campo). ADR 0020, P5.
  - **Referencias a otras familias**: `REFERENCE_KWARGS` declara qué parámetros del constructor son el `id` de un objeto de una familia de material (por omisión `{"material": MaterialRegistry}`; un elemento térmico declara `ThermalMaterialRegistry`; uno con dos materiales declara las dos). El lector YAML valida y resuelve la referencia. ADR 0020, P3.
- Solver: clase normal (no abstracta). Decora con `@SolverRegistry.register`. Constructor recibe `assembler` como primer arg. Implementa `solve(F_ext_global, step_callback=None) -> U_global`.
  - Si avanza por pasos, invoca `self.assembler.prepare_all_steps(U_committed)` al empezar cada intento de paso, con el campo convergido del paso anterior y antes del primer ensamblaje, también al repetir un paso tras una bisección (ADR 0020, P5).
  - Si es una variante de arco que sólo cambia la restricción del paso, hereda de `ArcLengthSolver` y extiende `ArcProblem` redefiniendo `constraint` (`from solidum.math.solvers import ArcLengthSolver, ArcProblem`), como `IndirectDisplacementSolver`.
  - Si construye `StiffnessProperties` a mano (en vez de usar el corrector), pasa `near_nullspace=assembler.near_nullspace` (ADR 0018). Si resuelve un estático lineal sin corrector, verifica la solución con `check_linear_solution` de `solidum/math/solvers/model_checks.py` (ADR 0019).
- Comentarios: solo cuando el "por qué" no sea obvio (sigue las reglas del proyecto). No documentar el "qué" — los nombres ya lo dicen.
- Idioma: docstrings y mensajes al usuario en español; código en inglés.

## Templates

Con `--user <modulo>`, cada plantilla va en `solidum/user/<modulo>/<snake>.py` con dos diferencias:

- Un material de una familia propia del módulo hereda de la clase base de esa familia (no de `Material`) y se decora con su registro (ver *Dentro de un módulo de usuario*).
- Del programa principal sólo se importan rutas y nombres públicos, sin guion bajo; `tests/test_user_module_public_imports.py` lo comprueba. La plantilla de solver importa de `solidum.math.solvers._shared`, que sólo vale dentro del programa principal: en un módulo, el logger sale de `solidum.logging.get_logger("user.<modulo>")`, y una utilidad que sólo exista en un módulo privado (hoy `domain_is_symmetric`) se hace pública en el programa principal antes de usarla (reexport en el `__init__.py` del paquete y en su `__all__`), como pieza genérica.

### material — `solidum/materials/<snake>.py`

```python
import numpy as np
from solidum.core.material import Material
from solidum.registry import MaterialRegistry


@MaterialRegistry.register
class {{Name}}(Material):
    """{{Descripción de una línea de la ley constitutiva}}.

    Parameters
    ----------
    {{param}} : {{tipo}}
        {{descripción}}
    """
    STRAIN_DIM = {{1|3|6}}
    PRIMARY_STATE_VAR = {{'kappa' o None}}
    # Variables internas como {nombre: forma}, en el orden en que se guardan por
    # filas en el ensamblaje por lotes (ADR 0014). {} si no hay historia.
    STATE_SCHEMA = {{{'kappa': (), 'damage': ()} o {}}}
    # Opcional: kernel puntual @njit con firma MAT_SIG (solidum/math/batch/signatures.py)
    # que llame al return mapping; sin él el material sigue el camino por elemento.
    # BATCH_KERNEL = _{{snake}}_batch

    def __init__(self, {{params}}, density: float | None = None):
        # density es opcional (ADR 0008): el material la declara solo si el
        # análisis la requiere. Consumidores que la necesitan (peso propio,
        # matriz de masa) fallan con ValueError si la encuentran como None.
        # TODO: almacenar parámetros físicos como atributos.
        self.density = density

    def compute_state(self, strain, state_vars=None):
        # TODO: implementar la física
        # 1) Recuperar variables internas del state_vars (con defaults si es None).
        # 2) Calcular trial elástico.
        # 3) Aplicar criterio (fluencia, daño, etc.) y corregir si toca.
        # 4) Devolver (stress, tangent, new_state_vars).
        raise NotImplementedError
```

### element — `solidum/elements/<familia>/<snake>.py`

`<familia>` es la subcarpeta existente que le corresponda (`solid_2d`, `solid_3d`, `frame`, `thermal`); los elementos sueltos (`truss.py`, `cable.py`, `frame3d.py`) viven en `solidum/elements/`.

```python
import numpy as np
from typing import List
from solidum.core.element import Element
from solidum.core.node import Node
from solidum.core.material import Material
from solidum.registry import ElementRegistry


@ElementRegistry.register
class {{Name}}(Element):
    """{{Descripción de una línea}}."""
    DOF_NAMES = {{['ux', 'uy', ...]}}
    STRAIN_DIM = {{1|3|6}}
    N_INTEGRATION_POINTS = {{int}}
    # GEOMETRICALLY_NONLINEAR = True   # sólo si la cinemática es corrotacional o de grandes
    #                                  # desplazamientos (ADR 0019): su rigidez puede crecer al deformarse
    # Opcional (ADR 0014): cinemática compilada `detJ = kin(pt, coords, B)` con firma
    # KIN_SIG, la misma que use compute_element_state, para entrar en el camino por lotes.
    # BATCH_KINEMATICS = _batch_kin_{{snake}}      # firma KIN_SIG; NO lanza: devuelve det J <= 0 si degenera
    # BATCH_SHAPE_FUNCTIONS = staticmethod(_N_{{snake}})   # N(xi, eta[, zeta]) para points_global del post-proceso

    def __init__(self, element_id: int, nodes: List[Node], material: Material,
                 {{params_extra}}):
        if len(nodes) != {{n_nodes}}:
            raise ValueError(f"El elemento {{Name}} requiere exactamente {{n_nodes}} nodos.")
        # Asignar atributos extra ANTES de super() para que estén disponibles
        # si la subclase necesita usarlos durante init.
        # TODO: self.{{param}} = {{param}}
        super().__init__(element_id, nodes, material)
        # TODO: precalcular geometría inicial (L0, c, s, T, etc.)

    def compute_element_state(self, u_e: np.ndarray):
        # TODO: implementar la cinemática y la rigidez tangente.
        # 1) Calcular strain a partir de u_e (matriz B, transformación T, etc.).
        # 2) Llamar self.material.compute_state(strain, self.state.vars[ip]).
        # 3) Guardar en self.state.vars_trial[ip] y self.state.stresses_trial[ip].
        # 4) Ensamblar K_e y F_int_e con dV = detJ * weight * thickness.
        raise NotImplementedError

    # Sólo si el elemento guarda estado propio fuera de ElementState (ADR 0020, P8):
    # def commit_state(self) -> None:
    #     super().commit_state()
    #     # TODO: confirmar el estado propio (trial -> convergido).
    #
    # def compute_global_stiffness(self) -> np.ndarray:
    #     # TODO: guardar el trial propio, llamar a super() y restaurarlo al salir.
    #
    # Sólo si toma decisiones discretas entre pasos (ADR 0020, P5):
    # def prepare_step(self, U_committed: np.ndarray) -> None:
    #     # TODO: decidir con el estado convergido; idempotente para el mismo U_committed.
```

### solver — `solidum/math/solvers/<snake>.py`

```python
from typing import Optional

import numpy as np
import scipy.sparse.linalg as spla

from solidum.constants import ZERO_TOL
from solidum.math.convergence import ConvergenceCriterion
from solidum.math.solvers._shared import _log, domain_is_symmetric
from solidum.math.solvers.corrector import NewtonCorrector, default_calibration_scales
from solidum.math.solvers.model_checks import ensure_statically_restrained
from solidum.registry import SolverRegistry


@SolverRegistry.register
class {{Name}}:
    """{{Descripción de una línea}}."""

    def __init__(self, assembler, convergence: Optional[ConvergenceCriterion] = None,
                 {{otros_params}}, linear_algebra: str = "auto"):
        self.assembler = assembler
        # Política de convergencia (ADR 0007). Si es None, defaults del proyecto.
        self.convergence = convergence if convergence is not None else ConvergenceCriterion()
        self.linear_algebra = linear_algebra
        # Corrector de Newton compartido (ADR 0015): bucle, backend con
        # degradacion a LU, Newton modificado, line search y telemetria.
        self.corrector = NewtonCorrector(
            self.convergence, max_iter=20,
            is_symmetric=domain_is_symmetric(assembler.domain),
            is_positive_definite=True, linear_algebra=linear_algebra,
            # Modos de cuerpo rígido para el AMG del solver iterativo
            # (ADR 0018); perezoso: sólo se calcula si se pide 'iterative'.
            near_nullspace=assembler.near_nullspace,
        )
        # TODO: parámetros adicionales. Imposición de Dirichlet: usa
        # assembler.reduce(K, F, U_current=..., load_factor=...) y
        # assembler.expand(u_red, T, g) — eliminación directa (ADR 0004).

    def solve(self, F_ext_global: np.ndarray, step_callback=None) -> np.ndarray:
        domain = self.assembler.domain
        ndof = domain.total_dofs
        U = np.zeros(ndof)
        # Red de seguridad (ADR 0019) — SÓLO si el solver es ESTÁTICO: un
        # mecanismo rígido no tiene equilibrio y se rechaza aquí con el
        # movimiento libre descrito. En un solver dinámico o modal NO se
        # llama: un modelo libre es legítimo cuando hay masa o capacidad.
        ensure_statically_restrained(self.assembler, type(self).__name__)
        # Por paso: construir un objeto con el protocolo NewtonProblem
        # (assemble, residual, residual_norm, calibration_scales,
        # reference_force, x_norm, correction, apply, on_converged) y
        #   res = self.corrector.run(problem, U, check_initial=False)
        #   if not res.converged: raise res.divergence_error(last_load_factor=...)
        #   U = res.x
        # Ver solidum/math/solvers/nonlinear.py (_IncrementalProblem) como
        # caso minimo; el corrector calibra el criterio (ADR 0007) con
        # calibration_scales -> default_calibration_scales(F_ext, F_int, K).
        # TODO: implementar el control de paso
        raise NotImplementedError
```

### test esqueleto — `tests/test_<snake>.py`

Para **material**:
```python
import numpy as np
from solidum.materials.{{snake}} import {{Name}}


def test_{{snake}}_construye_y_responde():
    mat = {{Name}}({{params_minimos}})
    assert mat.STRAIN_DIM in (1, 3, 6)
    strain = 0.0 if mat.STRAIN_DIM == 1 else np.zeros(mat.STRAIN_DIM)
    stress, tangent, new_state = mat.compute_state(strain)
    # TODO: aserción de física esperada en el caso trivial
```

Para **element**:
```python
import numpy as np
from solidum.core.domain import Domain
from solidum.elements.{{familia}}.{{snake}} import {{Name}}  # sin {{familia}} si vive en solidum/elements/
from solidum.materials.elastic import Elastic1D  # ajustar al material correcto


def test_{{snake}}_construye_y_calcula_K():
    d = Domain()
    nodes = [d.add_node(i+1, [float(i), 0.0]) for i in range({{n_nodes}})]
    mat = Elastic1D(E=1.0)  # o el material adecuado a STRAIN_DIM
    el = {{Name}}(1, nodes, mat, {{params_minimos}})
    d.add_element(el)
    d.generate_equation_numbers()

    n_dof = len(el.DOF_NAMES) * len(nodes)
    K, F = el.compute_element_state(np.zeros(n_dof))
    assert K.shape == (n_dof, n_dof)
    assert F.shape == (n_dof,)
```

Para **solver**:
```python
import numpy as np
from solidum.core.domain import Domain
from solidum.math.assembly import Assembler
from solidum.math.solvers.{{snake}} import {{Name}}


def test_{{snake}}_resuelve_caso_trivial():
    # TODO: construir un dominio mínimo (1 elemento, condiciones de frontera)
    # y verificar que solve() devuelve un vector de la dimensión correcta.
    pass
```

## Dentro de un módulo de usuario (`--user <modulo>`)

Un módulo de usuario es una formulación no estándar (no cumple el criterio de admisión, sección *Destino*) en su propia carpeta, que el programa principal no importa ni nombra, como los elementos, materiales y macros de usuario de FEAP (ADR 0020). Modelos: `discontinuities` (`solidum/user/discontinuities/`: elemento con estado propio y familia de material propia) y `dissipation_control` (`solidum/user/dissipation_control/`: un solver, variante de arco).

**Rutas**

| | Programa principal | Módulo de usuario `<modulo>` |
|---|---|---|
| Código | `solidum/materials/`, `solidum/elements/<familia>/`, `solidum/math/solvers/` | `solidum/user/<modulo>/<snake>.py` |
| Test | `tests/test_<snake>.py` | `tests/user/<modulo>/test_<pref>_<snake>.py` |
| Spec | `docs/specs/<Name>.md` | `docs/user/<modulo>/specs/<Name>.md` |
| Catálogo | `docs/catalogo_<elementos\|materiales\|solvers>.md` | sección con H1 exacto «# Módulos de usuario — <modulo>» al final del mismo catálogo (el build del manual de referencia la lleva al apéndice del módulo; sin el nombre del módulo, el build aborta) |

`<pref>` es un prefijo propio del módulo (`disc` en `discontinuities`): los nombres de archivo de test deben ser únicos en toda la suite, porque pytest los importa por nombre de archivo y los tests se importan entre sí por ese nombre.

**Carga explícita, nunca automática.** `import solidum` no carga nada de `solidum/user/`. El modelo que lo necesita lo pide:
- YAML: `user_modules: [<modulo>]` en el primer nivel. Sin esa clave, un tipo de elemento, un tipo de solver o una sección YAML del módulo son un error que sugiere declararla.
- Python y tests: `solidum.load_user_module("<modulo>")` (importa recursivamente el módulo para que sus decoradores registren) o importar desde `solidum.user.<modulo>`.

**Pasos propios del módulo**
1. Si el módulo no existe, crea `solidum/user/<modulo>/__init__.py` (docstring con la formulación, la referencia bibliográfica, el contenido y el uso; reexports con `__all__`), `docs/user/<modulo>/specs/` y `tests/user/<modulo>/`.
2. **Registro.** Los registros del programa principal (`MaterialRegistry`, `ElementRegistry`, `SolverRegistry`) valen tal cual. Si el módulo necesita una **familia de material propia** (una ley que no relaciona σ con ε, como la cohesiva `t`–`[[u]]`), la declara en `solidum/user/<modulo>/registry.py` como subclase de `solidum.registry.Registry` con sólo atributos públicos: `SPEC_KIND`, `YAML_SECTION` y `YAML_LABEL` (ADR 0020, P1). El almacén `_items` lo crea la base; no declares atributos con guion bajo. El lector YAML recorre su sección sin conocerla, `tools/spec.py` acepta su `kind` y los mensajes de error la nombran por su `YAML_LABEL`. Sus materiales se decoran con ese registro, y su plantilla de spec va en `docs/user/<modulo>/specs/_template_<SPEC_KIND>.md`. Modelo: `solidum/user/discontinuities/registry.py`.
3. **Constantes** del módulo en `solidum/user/<modulo>/constants.py`, no en `solidum/constants.py`.
4. **Del programa principal, sólo su contrato público**: `Element`, `Material` y `Node` (`solidum.core`), `Registry` y los registros (`solidum.registry`), `prepare_step`, `REFERENCE_KWARGS`, la cinemática `compute_kinematics_tri3` / `compute_integrands` (`from solidum.elements.solid_2d import …`), `ArcLengthSolver` y su problema por paso `ArcProblem` para una variante de arco que sólo cambia la restricción (`from solidum.math.solvers import …`; modelo: `solidum/user/dissipation_control/dissipation_arclength.py`), `get_logger` (`solidum.logging`)… Nunca una ruta ni un nombre con guion bajo (`_shared`, `_kind`): `tests/test_user_module_public_imports.py` lo rechaza. **No edites el programa principal para nombrar el módulo.** Si falta una pieza genérica, la IA la decide y la añade al programa principal, justificada por casos reales ajenos al módulo (Reglas.md §1 y §4), y la documenta (en el ADR 0020, o en un ADR nuevo si el cambio es grande). No se le presenta al usuario como decisión. Nunca un caso especial.
5. **Tests.** Empiezan con `solidum.load_user_module("<modulo>")`. Los barridos de contrato del programa principal sólo recorren sus propias clases: el módulo trae el suyo y valida sus specs (modelo: `test_disc_specs.py`). Para los elementos, reutiliza la batería de `tests/test_element_contract_sweep.py` con un arnés propio, `arnes(nombres, fabricas, …)`, puesto delante en la herencia de cada clase `Test*` (modelo: `tests/user/discontinuities/test_disc_contratos.py`). No parchees variables globales de aquel archivo ni uses sus nombres con guion bajo.
6. **Pureza.** `tests/test_main_program_purity.py` aplica tres controles: (1) estructural, por análisis sintáctico: el programa principal no importa el módulo de ninguna forma (absoluta, relativa, perezosa, por cadena ni con `load_user_module` y nombre literal); de `solidum.user` sólo usa el cargador genérico; (2) léxico: no menciona lo que aporta el módulo, con los nombres derivados al cargarlo (lo que añade a los registros, sus familias con su `YAML_SECTION` y `SPEC_KIND`, sus referencias nuevas, su `__all__` y su ruta); (3) en ejecución: `import solidum` no lo carga y toda clase registrada es del programa principal. No hay que dar de alta el módulo. Un nombre característico que la derivación no vea (p. ej. una constante que no reexporta) va, si se quiere, a `NOMBRES_ADICIONALES` de ese archivo.
7. **Sólo añadir.** Un módulo añade entradas a los registros; nunca reasigna una existente (p. ej. otra clase como `Elastic2D` o un solver con el nombre de uno del programa principal). Si lo hace, `load_user_module` lanza `ValueError` y deshace la carga entera. Da a sus clases nombres propios.

Test esqueleto en un módulo de usuario (cabecera; el cuerpo, como en las plantillas de arriba):
```python
import numpy as np
import solidum

solidum.load_user_module("{{modulo}}")
from solidum.user.{{modulo}}.{{snake}} import {{Name}}  # noqa: E402
```

## Recordatorios al cerrar

Tras generar y validar el esqueleto, comunica al usuario:
- La clasificación (estándar / módulo de usuario) y su justificación en una línea.
- Lista de archivos creados con paths clickables.
- Que el componente ya está **automáticamente registrado y disponible en YAML** sin tocar otros archivos (en un módulo de usuario, en cuanto el YAML declara `user_modules: [<modulo>]`).
- Dónde está el `TODO: implementar la física`.
- (Opcional) Snippet YAML mínimo de cómo invocar el nuevo tipo desde `examples/*.yaml`.

## Actualizar el catálogo

Una vez implementada la física de la spec y en verde los tests de `acceptance`, **añadir una entrada** al catálogo correspondiente (un componente de módulo de usuario, en su sección «Módulos de usuario — <modulo>»), con la clasificación del criterio de admisión y su justificación (Reglas.md §1):
- material → [docs/catalogo_materiales.md](../../../docs/catalogo_materiales.md)
- element → [docs/catalogo_elementos.md](../../../docs/catalogo_elementos.md)
- solver → [docs/catalogo_solvers.md](../../../docs/catalogo_solvers.md)

Seguir el formato de las entradas existentes (terso: propósito, DOFs/STRAIN_DIM, integración, parámetros, hipótesis/limitaciones, referencia bibliográfica si aplica, link al archivo). Esto cierra el ciclo extensión → documentación y mantiene el catálogo como índice navegable de las primitivas físicas del programa.
