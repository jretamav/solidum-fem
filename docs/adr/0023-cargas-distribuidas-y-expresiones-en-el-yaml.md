# ADR 0023 — Cargas distribuidas y expresiones en el YAML

- **Estado**: aceptado (2026-10-27)
- **Fecha**: 2026-10-27
- **Alcance**: `solidum/utils/expressions.py` (expresiones seguras), `solidum/utils/yaml_loads.py` (resolución de las cargas distribuidas), `YamlParser` (claves `tractions`, `tractions_by_group`, `thermal_loads.boundary_flux_by_group`; `body_force`, `body_source` y `boundary_flux` con expresiones), `GmshParser` (celdas de cada grupo físico). Salda las deudas #23 y #24.

## Contexto

El ADR 0022 dio a los elementos cargas variables en el API de Python: fuerza de cuerpo `b(x)`, tracción `t(x, n)` sobre la geometría real del borde o de la cara (con ella, la presión `−p·n`), fuente `Q(x)` y flujo `q̄(x, n)`. El YAML seguía como antes:

- sólo cargas nodales, fuerza de cuerpo uniforme y peso propio;
- ninguna tracción de borde (deuda #23): el ejemplo de Kirsch del manual tuvo que cargar por desplazamiento impuesto;
- ninguna presión (deuda #24);
- en el térmico, fuente y flujo constantes, éste sólo por elemento.

Quien trabaja con archivos YAML no tenía nada de lo anterior.

Los grupos físicos de gmsh se guardaban como conjuntos de nodos. Con los nodos de un grupo de contorno no se sabe con seguridad qué bordes cargar: en un elemento con dos nodos del grupo que no forman un borde, o en una franja de un elemento de espesor, se cargarían bordes que no son del contorno.

## Decisión

### 1. Las cargas distribuidas del YAML

```yaml
body_force: [0.0, "-78.5e3*(1 + 0.1*y)"]
tractions_by_group:
  Borde_Carga: {t: [1.0e6, 0.0]}
  Agujero:     {pressure: "1e6*(1 - y)"}
tractions:
  - {element: 12, edge: 1, t: ["1e3*nx", 0.0]}
  - {element: 40, face: 3, pressure: 2.5e3}
thermal_loads:
  body_source: {Q: "5e3*exp(-x)"}
  boundary_flux_by_group:
    Cara_caliente: {q: -200.0}
```

- `tractions` (por elemento, `edge` en 2D y `face` en 3D) y `tractions_by_group` (por grupo físico de gmsh) llaman a `compute_edge_traction` / `compute_face_traction`. Cada entrada lleva una de dos formas: `t`, en ejes globales, o `pressure`.
- **Convención de signo de la presión** (decisión del usuario, Reglas §5): `p > 0` es compresión, `t = −p·n` con `n` la normal exterior; es la de Abaqus y ANSYS y la que ya usaba el API. La carga es muerta: actúa sobre la geometría de referencia.
- En el térmico, `boundary_flux_by_group` completa `boundary_flux`, con la convención de siempre (`q̄ > 0` saliente).
- Las tracciones entran en `get_external_forces`, junto a las cargas nodales. Por eso se escalan y se aplican igual que ellas en todos los análisis: estático, no lineal, dinámico y armónico.

### 2. Grupos de contorno exactos

`GmshParser` guarda, además de los nodos, las celdas de cada grupo físico (`domain.physical_group_cells`). Una carga por grupo se aplica sobre los bordes o caras de elemento que **coinciden por sus nodos** con una celda del grupo y pertenecen a un único elemento. El índice se construye sólo con los elementos que admiten la carga, así que un modelo con elementos mecánicos y térmicos sobre los mismos nodos no confunde un borde de contorno con uno compartido.

Se rechazan al leer el archivo dos casos:
- las celdas de un grupo de superficie, que no son bordes;
- las celdas compartidas por dos elementos, que no son contorno.

### 3. Expresiones seguras

Un valor de carga es un número o una cadena con una expresión de `x`, `y`, `z` y, en las cargas de contorno, de la normal exterior `nx`, `ny`, `nz`.

**Cómo se evalúa.** La expresión se analiza con `ast`, y se compila sólo si contiene números, esas variables, `+ − * / **`, el signo, `pi` y una lista cerrada de funciones (`sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `atan2`, `sinh`, `cosh`, `tanh`, `exp`, `log`, `log10`, `sqrt`, `abs`, `min`, `max`). Se compila a una función de Python a partir del árbol ya validado, sin `__builtins__`; nunca se hace `eval` sobre el texto. Las constantes enteras pasan a reales, para que `9**9**9` desborde de inmediato en vez de colgar el programa.

**Lo que se rechaza al leer el archivo**, con la entrada y lo no admitido:
- atributos, índices, llamadas a otras funciones, `lambda`, condicionales y comparaciones;
- nombres desconocidos;
- la normal en una carga de volumen.

Una expresión que no se puede evaluar en un punto se informa con el punto.

### 4. Validación al leer, no al resolver

Las cargas se resuelven en `YamlParser.parse`, tras construir la malla (`_parse_distributed_loads`), y todos los errores se informan juntos (`YamlValidationError`). Se detectan así:
- elemento inexistente, índice fuera de rango, `t` y `pressure` a la vez o ninguno;
- número de componentes distinto de la dimensión, elemento que no admite la carga;
- grupo inexistente o que no es de contorno, claves no admitidas y expresiones inválidas.

Antes, un `boundary_flux` sin `edge` ni `face` sólo fallaba al ensamblar.

## Consecuencias

- **Deudas #23 y #24 saldadas.**
- **Coste del componente N+1.** Un elemento con `compute_edge_traction` / `compute_face_traction` (o `compute_edge_flux` / `compute_face_flux`) y `EDGE_NODES` / `FACE_NODES` recibe estas cargas desde el YAML sin tocar el lector.
- **Pruebas** (`tests/test_yaml_distributed_loads.py`):
  - lista blanca de las expresiones;
  - igualdad exacta con el API de Python para tracción, presión, fuerza de cuerpo, fuente y flujo;
  - resultantes analíticas con la malla de Kirsch: tracción uniforme en `Borde_Carga`, `ΣFx = σ·H·e`; presión en el `Agujero`, `ΣF = p·a·e·(1, 1)`, que confirma a la vez el signo y los bordes del arco;
  - presión en una cara alabeada de un hexaedro;
  - placa a tracción uniforme resuelta con `run_yaml`;
  - todos los errores.
- **Limitaciones:**
  - sin cargas que varían en el tiempo (variable `t`) en el YAML;
  - sin presión seguidora;
  - `tractions_by_group` sólo en 2D, porque `GmshParser` sólo importa mallas 2D lineales; en 3D se declara cara por cara;
  - el ejemplo de Kirsch del manual sigue cargado por desplazamiento impuesto: pasarlo a tracción cambiaría el capítulo, y queda para después.

## Alternativas descartadas

- **`eval` sobre el texto o `sympy.sympify`.** El primero ejecuta código arbitrario de un archivo ajeno. El segundo arrastra SymPy, que no es dependencia del programa principal, y no impide por sí solo llamadas peligrosas.
- **Deducir los bordes de un grupo de sus nodos.** Falla en los casos del contexto; las celdas de gmsh dicen exactamente qué bordes son.
- **Un bloque `pressure` separado de `tractions`.** La presión es una forma de la tracción y comparte entidades, grupos y validación.
