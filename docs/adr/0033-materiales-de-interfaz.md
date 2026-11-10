# ADR 0033 — Materiales de interfaz en el programa principal

- **Estado**: aceptado (2026-11-10)
- **Fecha**: 2026-11-10
- **Alcance**:
  - `solidum/core/interface_material.py`: clase base `InterfaceMaterial`;
  - `InterfaceMaterialRegistry` (`solidum/registry.py`), tipo de spec `interface_material`, sección YAML `interface_materials`;
  - `solidum/materials/interface/`: las leyes de la familia;
  - el módulo de usuario `discontinuities`, que deja de declarar su familia cohesiva.
- **Extiende**: el ADR 0010 (familia cohesiva paralela a `Material`) y el ADR 0020 (piezas genéricas para los módulos de usuario).

## Contexto

La ley cohesiva tracción-salto nació en el ADR 0010 como familia paralela a `Material` y en el ADR 0020 se fue al módulo de usuario `discontinuities`, porque su único consumidor era el elemento con discontinuidad embebida. Reglas.md §4 dejó escrito cuándo volvería: «cuando un componente del programa principal la necesite (un elemento de interfaz estándar), pasa al programa principal sin física nueva y el módulo la reutiliza».

Ese componente es la interfaz de adherencia-deslizamiento entre una barra embebida y el sólido (decisiones físicas del usuario del 2026-11-10). Su ley relaciona la tracción con el desplazamiento relativo, como la cohesiva: es el mismo tipo de objeto. Las dos difieren sólo en la geometría de la interfaz. Una grieta es una superficie, con una normal y d − 1 tangentes. La interfaz de una barra es una línea, con d − 1 normales y una tangente.

## Decisión

1. **Una familia, no dos.** `InterfaceMaterial` relaciona la tracción `t` con el desplazamiento relativo `Δu`, ambos en el sistema local de la interfaz: `compute_traction(Δu, estado) → (t, ∂t/∂Δu, estado')`, con el estado de prueba y confirmado como los materiales. Se registra en `InterfaceMaterialRegistry` y se declara en la sección YAML `interface_materials`.
2. **Geometría declarada.** `GEOMETRY ∈ {"surface", "line"}` fija qué componentes son normales y cuáles tangenciales, siempre en el orden normales primero y tangenciales después: superficie `(n, s)` / `(n, s1, s2)`, línea `(n, s)` / `(n1, n2, s)`. En 2D las dos geometrías tienen el mismo orden. `JUMP_DIM` es el número de componentes que admite la ley; `None` significa que admite las de 2D y las de 3D. El consumidor valida la ley con `accepts(geometry, n_components)` al construirse.
3. **Lo que es del elemento.** El paso entre ejes globales y locales, y el ancho sobre el que actúa cada componente, son del elemento. En una barra, la tangencial actúa sobre el perímetro y la normal sobre el diámetro. Así la ley no conoce la geometría de quien la usa.
4. **La ley cohesiva pasa con la familia, sin cambio de física.** `CohesiveDamageIsotropic` es estándar según el criterio de admisión: de Borst et al. 2012, cap. 6, y los elementos cohesivos con daño de Abaqus. Ahora está en `solidum/materials/interface/`, con su spec en `docs/specs/` y sus pruebas en `tests/`. El elemento del módulo conserva su parámetro `cohesive_material`, que ahora resuelve en `interface_materials`.

## Consecuencias

- **El YAML cambia de nombre de sección**: `cohesive_materials` pasa a ser `interface_materials`. Un modelo que use la discontinuidad embebida sigue declarando `user_modules: [discontinuities]` para el elemento, pero su ley cohesiva ya no la necesita.
- **El módulo queda con lo suyo**: el elemento, el estado de la discontinuidad y el seguimiento de grieta. Las pruebas de pureza (`tests/test_main_program_purity.py`) siguen verdes: el programa principal no nombra nada que aporte el módulo, y la lista de familias del programa principal gana `interface_materials`.
- **El manual de referencia** tiene un capítulo nuevo, «Modelos Constitutivos — Interfaz».
- **La ley de adherencia** y, más adelante, la de modo mixto de las grietas son materiales de esta familia. La segunda servirá también al elemento del módulo, que usa la geometría de superficie.

## Alternativas descartadas

- **Una familia nueva para la adherencia, junto a la cohesiva del módulo.** Habría dos familias para el mismo tipo de ley, y el programa principal dependería de un módulo de usuario o duplicaría su contrato (Reglas §1: centralizar cuando aparece el segundo caso real).
- **Mantener `cohesive_materials` como alias de la sección.** Un alias es una segunda forma de decir lo mismo que hay que documentar y probar. El único consumidor de la sección vieja era el módulo, cuyos modelos y pruebas se actualizan en este mismo cambio.
