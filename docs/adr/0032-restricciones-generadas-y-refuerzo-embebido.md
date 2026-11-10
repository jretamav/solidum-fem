# ADR 0032 — Restricciones generadas por geometría y refuerzo embebido

- **Estado**: aceptado (2026-11-10)
- **Fecha**: 2026-11-10
- **Alcance**:
  - `ConstraintRegistry` (`solidum/registry.py`), un tipo de spec nuevo, `constraint`;
  - `solidum/mesh/point_location.py`: inversa del mapeo isoparamétrico y `PointLocator`;
  - `solidum/bc/embedded.py`: `EmbeddedReinforcement`, `EmbeddedBar` y `embed_reinforcement`;
  - bloque YAML `constraints`;
  - `Domain.generated_nodes`.
- **Spec**: [`EmbeddedReinforcement`](../specs/EmbeddedReinforcement.md).
- **Extiende**: el ADR 0004 (restricciones afines, MPC).

## Contexto

Las MPC del ADR 0004 se declaraban una a una: el usuario escribía cada esclavo, sus maestros y sus coeficientes. El refuerzo embebido necesita miles: una restricción por traslación en cada nodo de la barra, con coeficientes que son las funciones de forma del elemento de sólido que contiene el nodo. Esos coeficientes salen de la geometría (en qué elemento cae el punto y con qué coordenadas naturales), no del usuario. Los acoplamientos de dimensión mixta (viga-sólido) tienen la misma forma: restricciones afines calculadas a partir de la geometría de una cara. Son dos casos reales con la misma estructura (Reglas §1).

## Decisión

1. **Generador de restricciones.** Un objeto que, aplicado a un dominio ya mallado con `apply(domain)`, crea los nodos, los elementos y las MPC que su geometría pide, y devuelve lo que creó. Se registra en `ConstraintRegistry` (`SPEC_KIND = "constraint"`) y se declara en el YAML en el bloque `constraints` (lista de `{type, …}`). El lector aplica el bloque **después** de crear los elementos y **antes** de los apoyos. Igual que las secciones (ADR 0021 §2), una clase que referencia objetos de las familias declara `from_yaml(data, resolve)`. El lector no conoce ningún tipo concreto.
2. **Localización de puntos genérica.** La inversa de $\mathbf x = \sum_i N_i(\boldsymbol\xi)\mathbf X_i$ se resuelve por Newton con lo que todo sólido ya declara: `REFERENCE_CELL` y `batch_shape_functions`. El jacobiano se aproxima por diferencias centradas. Es exacto en los elementos de Solidum, cuyas funciones de forma son a lo sumo cuadráticas en cada coordenada natural. Así sirve a los diez sólidos y a los de los módulos de usuario sin código por tipo. `PointLocator` filtra los candidatos por la caja envolvente, ampliada un 5 % para los lados curvos, y asigna el punto de un borde común al elemento de menor id: el resultado es determinista.
3. **Refuerzo embebido** (decisiones físicas del usuario en la spec: formulación A, área sin descontar). La barra es una cadena de elementos 1D ordinarios, con nodos propios. La polilínea se corta en los bordes de los anfitriones; el punto de salida se busca por bisección sobre la recta, con la pertenencia a la celda. Cada nodo queda atado con $u_d(\mathbf X_b) = \sum_i N_i(\boldsymbol\xi_b)\,u_{d,i}$ por traslación. Son MPC de **unión**, porque una traslación rígida las cumple ($\sum N_i = 1$). Por eso su fuerza sale en `constraint_forces` (ADR 0004 §3) y es la fuerza de adherencia concentrada en cada nodo de la barra: `EmbeddedBar.bond_forces(result)`. No hace falta nada nuevo en el ensamblador ni en los solvers: la eliminación $\mathbf T^\top\mathbf K\mathbf T$, $\mathbf T^\top\mathbf M\mathbf T$ ya existe en la estática, en el modal y en los transitorios.
4. **Nodos generados.** `Domain.generated_nodes` guarda los nodos que crea Solidum. Los selectores por coordenada del YAML (`*_by_coord`) no los toman: un apoyo o una carga por coordenada es del sólido. Un apoyo en `x_min` fijaría además el nodo extremo de la barra, que ya es esclavo, y daría un conflicto.

## Alternativas descartadas

- **La barra sin nodos propios, integrada en el anfitrión** (variante B de ANSYS y DIANA). Es idéntica en los anfitriones lineales y difiere ≤ 0.18 % en los demás casos medidos. Pero exige código dentro de cada familia de anfitrión y no deja nodos a los que atar una ley de adherencia-deslizamiento. Decisión del usuario.
- **Penalización** (como `ASDEmbeddedNodeElement` de OpenSees). Introduce un parámetro sin significado físico y empeora el condicionamiento. La eliminación del ADR 0004 es exacta.
- **Una familia YAML (`YAML_SECTION`)**. Las familias se construyen antes que los elementos y se referencian por id; un generador necesita los elementos ya creados.

## Medido

`tests/test_embedded_reinforcement.py` (29 pruebas):
- inversa del mapeo a 1e-9 en los diez sólidos, con mallas curvas;
- deformación uniforme con la barra a un ángulo cualquiera: $\varepsilon_b = \mathbf t^\top\boldsymbol\varepsilon\,\mathbf t$ y energía $\tfrac12 E_sA_s(\mathbf t^\top\boldsymbol\varepsilon\mathbf t)^2L_b$ a 1e-10 (Tri3, Quad4 distorsionado, Hex8, Tet4);
- 3 y 6 modos rígidos;
- malla de la barra independiente (1e-10);
- placa en tracción con $E_cA_g + E_sA_s$ a 1e-10;
- fluencia en $A_s\sigma_y$ a 1e-8;
- viga contra la sección transformada sin descuento: Q8 20×4 a +0.06 %; Quad4 converge de −10 % a −0.65 %;
- fuerzas de adherencia en equilibrio con la barra;
- masa;
- YAML igual a Python en 2D y 3D.

## Fuera de este ADR

- La adherencia-deslizamiento, que es otro componente: sustituye la MPC de cada nodo por una ley τ–s.
- El acoplamiento viga-sólido, que usará `ConstraintRegistry` y `PointLocator`.
- `hosts` por grupo físico en el YAML: hoy son ids de elementos.
