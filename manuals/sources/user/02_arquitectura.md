# Arquitectura del Programa

## Separación Elemento ↔ Material

Cada elemento finito declara dos contratos sobre el material que acepta:

- **`STRAIN_DIM`**: dimensión Voigt de la deformación que entrega al material. `1` para elementos axiales (truss, cable, marcos — la fibra centroidal); `3` para 2D (`[εxx, εyy, γxy]`); `6` para 3D continuo (no implementado todavía).
- **`DOF_NAMES`**: lista de grados de libertad por nodo. Por ejemplo, `['ux', 'uy', 'rz']` para un marco 2D.

El protocolo del material es uniforme:

```python
sigma, E_tangent, new_state = material.compute_state(epsilon, state_vars)
```

donde `epsilon` tiene dimensión `STRAIN_DIM`, `sigma` es el esfuerzo en la misma representación, `E_tangent` es la matriz tangente consistente y `state_vars` encapsula las variables internas (deformación plástica, daño acumulado, etc.).

Un material declara igualmente su `STRAIN_DIM`; el dominio rechaza al construir cualquier emparejamiento incompatible.

## Auto-Registro de Componentes

Elementos, materiales y solvers se registran automáticamente al arrancar el programa mediante decoradores:

```python
@ElementRegistry.register
class Truss2D(Element):
    DOF_NAMES = ['ux', 'uy']
    STRAIN_DIM = 1
    ...
```

Esto permite que el archivo YAML los referencie por nombre (`type: Truss2D`) sin que el usuario tenga que importar nada manualmente. El módulo `fenix.autodiscover` recorre los paquetes `elements/`, `materials/` y `math/solvers/` en la primera importación y puebla los registries.

## Validación Temprana del Input

El parser YAML acumula todos los errores detectados en una sola pasada (no aborta al primero) y los presenta juntos. Comprueba:

- Ausencia de bloques obligatorios (`nodes`, `materials`, `elements`, `solver`).
- IDs duplicados de nodo, material o elemento.
- Referencias a materiales o nodos inexistentes desde `elements` y `boundary_conditions`.
- Tipos de elemento, material y solver no registrados.
- **Parámetros no aceptados por el constructor del elemento.** Si un elemento define `A` e `I` como parámetros, escribir `large_strains: true` en su entrada YAML produce un error explícito que indica los parámetros válidos para ese tipo.
