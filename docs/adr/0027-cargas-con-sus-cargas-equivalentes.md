# ADR 0027 — Cargas que llevan sus cargas equivalentes por elemento

- **Estado**: aceptado (2026-11-09)
- **Fecha**: 2026-11-09
- **Alcance**:
  - `solidum/core/load_vector.py`: `LoadVector`, `EquivalentLoads` y `EquivalentLoadLayout`;
  - `Assembler.assemble_self_weight` y `Assembler.assemble_body_load`, que devuelven un `LoadVector`;
  - `build_solve_result` y `solidum.run` (también por etapas, ADR 0026);
  - el historial de fuerzas internas de los transitorios (`ElementForcesRecorder`).
- **Cierra**: la deuda #35. **Sustituye**: el §5 del ADR 0026 (composición por la versión del registro del ensamblador).

## Contexto

Las fuerzas internas de extremo de un elemento estructural con carga distribuida son $\mathbf F_\text{int} - \mathbf f_\text{eq}$, donde $\mathbf f_\text{eq}$ es la carga nodal equivalente de esa carga en el elemento (ADR 0002). Para descontarla hay que saber **qué cargas distribuidas contiene el vector aplicado y con qué factor**.

Hasta ahora el ensamblador guardaba sólo las cargas equivalentes de la **última** carga de cuerpo ensamblada. Medido:
- con el peso propio y una fuerza de cuerpo lateral en una misma corrida, la fuerza axial de una columna salía constante (−18 193 N) cuando debía variar de −20 792 N a −15 594 N, sin aviso;
- dos casos de carga ensamblados con el mismo ensamblador descontaban las cargas del último ensamblado en los dos;
- el historial de fuerzas internas de un transitorio no descontaba ninguna, y no había forma de saber con qué factor actuaba cada carga en cada instante de `F(t)`.

El ADR 0026 lo resolvió a medias para las etapas, con un contador de versión. El problema es general: lo tiene cualquier combinación de cargas, en cualquier análisis.

## Decisión

La carga lleva su propia información. `LoadVector` es un `np.ndarray` (se usa igual que antes) con un atributo `equivalent_loads`: las cargas equivalentes por elemento de su parte distribuida. La aritmética de las cargas las combina linealmente:

| Operación | Cargas equivalentes del resultado |
|---|---|
| `F1 + F2`, `F1 − F2`, `−F` | suma, resta, signo |
| `a·F`, `F/a` con `a` escalar | escaladas |
| `F + P` con `P` un `np.ndarray` simple (cargas puntuales) | las de `F` |
| cualquier otra (producto elemento a elemento, normas, comparaciones, indexado, `np.asarray`) | ninguna: el resultado es un `np.ndarray` simple, porque ya no es una carga |

Copia, `deepcopy` y serialización las conservan.

Así, `F = 1.2·F_peso + 1.6·F_viento` lleva $1.2\,\mathbf f_\text{eq,peso} + 1.6\,\mathbf f_\text{eq,viento}$ sin que el usuario haga nada. Lo mismo vale para la carga total de una etapa ($\mathbf F_c + \lambda\,\mathbf F_\text{ref}$) y para la carga de un transitorio en cada instante ($\mathbf F_\text{peso} + a(t)\,\mathbf F_\text{viento}$).

**Representación.** Los vectores elementales se guardan concatenados en un arreglo, con una disposición (ids de elemento y offsets) compartida por todas las cargas del mismo ensamblador. Combinar dos cargas con la misma disposición es una operación vectorizada, también en un transitorio de muchos pasos. Con disposiciones distintas, se fusionan por elemento.

**Consumidores.**
- `build_solve_result` toma las cargas equivalentes de la carga aplicada. `SolveResult.F_applied` y `F_constant` son `LoadVector`, así que una etapa siguiente hereda las suyas (ADR 0026).
- El historial de fuerzas internas de los transitorios (`record_internal_forces=True`) descuenta en cada paso las de `F(t)`.

**Compatibilidad.** Un `np.ndarray` simple perdió sus cargas equivalentes: lo convirtieron con `np.asarray`/`np.array`, o se sumó con `+=` sobre un arreglo de ceros.
- En un análisis estático se recurre al registro anterior del ensamblador, la última carga de cuerpo, como antes. Si se ensambló más de una, se **avisa**, porque la suposición puede ser falsa.
- En un transitorio no se puede saber el factor en cada instante: se avisa una vez y no se descuentan.

## Consecuencias

- La deuda #35 queda cerrada, en estático, por etapas y en el tiempo, sin un modelo de cargas nuevo ni cambios en la API: las funciones de carga devuelven el mismo vector, que ahora sabe lo que contiene.
- Un resultado falso sin aviso deja de ser posible: o la carga lleva sus cargas equivalentes, o se avisa.
- Una carga distribuida nueva (por ejemplo, una carga de línea en un grupo de elementos) sólo tiene que devolver un `LoadVector` con las suyas.
- Pruebas: `tests/test_load_vector.py` (álgebra, dos cargas de cuerpo contra la suma a mano a 1e-10, dos casos con el mismo ensamblador, avisos, transitorio paso a paso a 1e-12).

## Alternativas consideradas

- **Registro acumulativo en el ensamblador.** Suma bien dos cargas en una corrida, pero falla con dos casos de carga ensamblados por separado: el ensamblador no sabe qué vectores entran en `F_applied`.
- **Un modelo de cargas con objetos «caso de carga» y combinaciones declaradas** (mejora 4 de la lista del 2026-09-22). Es la solución de los códigos comerciales, pero cambia la API de todas las cargas. `LoadVector` da el mismo resultado sobre la API actual, y un modelo de cargas futuro puede construirse encima.
- **Pasar las cargas equivalentes aparte en cada llamada**. Traslada al usuario la contabilidad que causaba el error.

## Referencias

- ADR 0002 (fuerzas internas, `F_int − f_eq`), ADR 0008 (peso propio), ADR 0022 y 0023 (cargas variables y distribuidas), ADR 0026 (etapas de carga).
