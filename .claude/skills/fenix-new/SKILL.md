---
name: fenix-new
description: Genera scaffolding para nuevos materiales, elementos o solvers de Fenix FEM (archivo + decorador de registro + test esqueleto). Úsalo cuando el usuario pida añadir/crear un material, elemento o solver derivado de su investigación.
---

# fenix-new

Skill de scaffolding para Fenix FEM. Reduce el coste de incorporar nuevos componentes a la arquitectura: tú escribes la física, esta skill se encarga del boilerplate, del registro y del test esqueleto.

## Cuándo invocar

Activar siempre que el usuario pida añadir un nuevo:
- **Material** (constitutivo): elasticidad alternativa, plasticidad, daño, viscoelasticidad, hiperelasticidad, etc.
- **Elemento**: nueva formulación 1D/2D/3D (sólido, viga, cáscara, contacto…).
- **Solver**: lineal, no-lineal, path-following, dinámico, etc.

Disparadores típicos: "añade un material …", "crea un elemento …", "implementa un solver …", "/fenix-new …".

## Argumentos

`/fenix-new <kind> <Name>`

- `kind` ∈ {`material`, `element`, `solver`}
- `Name`: PascalCase (`HyperelasticOgden`, `Shell4`, `DynamicNewmark`).

Si faltan argumentos, pídelos al usuario.

## Flujo

1. **Recoge metadatos mínimos** preguntando al usuario lo estrictamente necesario (no satures):
   - **material**: `STRAIN_DIM` (1, 3 ó 6) · `PRIMARY_STATE_VAR` (str o `None`) · parámetros del constructor (lista `nombre: tipo`).
   - **element**: `DOF_NAMES` (lista) · `STRAIN_DIM` (1, 3 ó 6) · `N_INTEGRATION_POINTS` · número exacto de nodos requeridos · parámetros adicionales (`A`, `I`, `thickness`, etc.).
   - **solver**: tipo (`lineal` / `no-lineal` / `path-following` / `otro`) · parámetros del constructor con sus defaults.
2. **Genera el archivo** en su carpeta canónica (ver Templates) usando `snake_case(Name)` como nombre de archivo.
3. **Genera el test esqueleto** en `tests/test_<snake>.py` con un test mínimo de "se construye + compute_state/compute_element_state/solve devuelve shapes correctas".
4. **Verifica**: ejecuta `pytest tests/test_<snake>.py -x` y confirma que el esqueleto pasa.
5. **Comunica al usuario**: lista los archivos creados, marca dónde queda el `# TODO: implementar la física aquí`, y dale el ejemplo YAML mínimo para usarlo (si aplica a material/element).

**No hace falta** tocar `__init__.py`, `registry.py`, ni `yaml_parser.py`: el decorador + autodiscover hacen el resto.

## Convenciones del proyecto

- Material: hereda de `fenix.core.material.Material`. Decora con `@MaterialRegistry.register`. Implementa `compute_state(strain, state_vars=None) -> (stress, tangent, new_state_vars)`.
- Element: hereda de `fenix.core.element.Element`. Decora con `@ElementRegistry.register`. Implementa `compute_element_state(u_e) -> (K_e, F_int_e)`. La base se encarga de `commit_state`, `state_vars`, validación STRAIN_DIM, registro de DOFs y creación de `ElementState`.
- Solver: clase normal (no abstracta). Decora con `@SolverRegistry.register`. Constructor recibe `assembler` como primer arg. Implementa `solve(F_ext_global, step_callback=None) -> U_global`.
- Comentarios: solo cuando el "por qué" no sea obvio (sigue las reglas del proyecto). No documentar el "qué" — los nombres ya lo dicen.
- Idioma: docstrings y mensajes al usuario en español; código en inglés.

## Templates

### material — `fenix/materials/<snake>.py`

```python
import numpy as np
from fenix.core.material import Material
from fenix.registry import MaterialRegistry


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

    def __init__(self, {{params}}):
        # TODO: almacenar parámetros como atributos
        ...

    def compute_state(self, strain, state_vars=None):
        # TODO: implementar la física
        # 1) Recuperar variables internas del state_vars (con defaults si es None).
        # 2) Calcular trial elástico.
        # 3) Aplicar criterio (fluencia, daño, etc.) y corregir si toca.
        # 4) Devolver (stress, tangent, new_state_vars).
        raise NotImplementedError
```

### element — `fenix/elements/<snake>.py`

```python
import numpy as np
from typing import List
from fenix.core.element import Element
from fenix.core.node import Node
from fenix.core.material import Material
from fenix.registry import ElementRegistry


@ElementRegistry.register
class {{Name}}(Element):
    """{{Descripción de una línea}}."""
    DOF_NAMES = {{['ux', 'uy', ...]}}
    STRAIN_DIM = {{1|3|6}}
    N_INTEGRATION_POINTS = {{int}}

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
```

### solver — `fenix/math/solver_<snake>.py`

```python
import numpy as np
import scipy.sparse.linalg as spla
from fenix.constants import CONVERGENCE_TOL, ZERO_TOL
from fenix.registry import SolverRegistry


@SolverRegistry.register
class {{Name}}:
    """{{Descripción de una línea}}."""

    def __init__(self, assembler, tol=CONVERGENCE_TOL, {{otros_params}},
                 linear_algebra: str = "auto"):
        self.assembler = assembler
        self.tol = tol
        self.linear_algebra = linear_algebra
        # TODO: parámetros adicionales. Imposición de Dirichlet: usa
        # assembler.reduce(K, F, U_current=..., load_factor=...) y
        # assembler.expand(u_red, free_dofs, g) — eliminación directa (ADR 0004).

    def solve(self, F_ext_global: np.ndarray, step_callback=None) -> np.ndarray:
        domain = self.assembler.domain
        ndof = domain.total_dofs
        U = np.zeros(ndof)
        # TODO: implementar el algoritmo
        raise NotImplementedError
```

### test esqueleto — `tests/test_<snake>.py`

Para **material**:
```python
import numpy as np
from fenix.materials.{{snake}} import {{Name}}


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
from fenix.core.domain import Domain
from fenix.elements.{{snake}} import {{Name}}
from fenix.materials.elastic import Elastic1D  # ajustar al material correcto


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
from fenix.core.domain import Domain
from fenix.math.assembly import Assembler
from fenix.math.solver_{{snake}} import {{Name}}


def test_{{snake}}_resuelve_caso_trivial():
    # TODO: construir un dominio mínimo (1 elemento, condiciones de frontera)
    # y verificar que solve() devuelve un vector de la dimensión correcta.
    pass
```

## Recordatorios al cerrar

Tras generar y validar el esqueleto, comunica al usuario:
- Lista de archivos creados con paths clickables.
- Que el componente ya está **automáticamente registrado y disponible en YAML** sin tocar otros archivos.
- Dónde está el `TODO: implementar la física`.
- (Opcional) Snippet YAML mínimo de cómo invocar el nuevo tipo desde `examples/*.yaml`.

## Actualizar el catálogo

Una vez el usuario implemente la física y los tests pasen, **añadir una entrada** al catálogo correspondiente:
- material → [docs/catalogo_materiales.md](../../../docs/catalogo_materiales.md)
- element → [docs/catalogo_elementos.md](../../../docs/catalogo_elementos.md)
- solver → [docs/catalogo_solvers.md](../../../docs/catalogo_solvers.md)

Seguir el formato de las entradas existentes (terso: propósito, DOFs/STRAIN_DIM, integración, parámetros, hipótesis/limitaciones, referencia bibliográfica si aplica, link al archivo). Esto cierra el ciclo extensión → documentación y mantiene el catálogo como índice navegable de las primitivas físicas del programa.
