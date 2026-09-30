# ADR 0021 — Secciones de marco: familia de secciones y sección por fibras

- **Estado**: aceptado (2026-09-29)
- **Fecha**: 2026-09-29
- **Alcance**: registro nuevo `SectionRegistry` (familia YAML `sections`) y clase `FiberSection`; `Frame2DEuler`, `Frame2DEulerCorot` y `Frame3D` ganan el camino con sección; guarda de material elástico lineal en los cuatro marcos; `Material.IS_LINEAR_ELASTIC`; `Element.materials()`; gancho `from_yaml` para objetos de familia con referencias. Formulación física en [`docs/specs/FiberSection.md`](../specs/FiberSection.md).

## Contexto

La auditoría de pureza del 2026-09-29 verificó un defecto físico silencioso en los cuatro marcos: evalúan el material con la deformación axial `(u_j − u_i)/L` y construyen **toda** la rigidez local con la `E_t` resultante, y la fuerza interna de flexión es `K_local(E_t)·u`. Con fluencia axial (`E_t = 0`) el flector se anula aunque la sección esté curvada; con daño en ablandamiento (`E_t < 0`) cambia de signo; `Frame3D` pierde además el torsor (`G = E_t/(2(1+ν))`). Ningún texto ni código de referencia usa esa ley.

El usuario decidió el mismo día, una cuestión por vez (Reglas §3, física):

1. Implementar la **sección por fibras** en lugar de restringir los marcos a material elástico.
2. **Elemento en desplazamientos** (coherente con Reglas §0), no en fuerzas.
3. **Etapa 1**: `Frame2DEuler`, `Frame2DEulerCorot`, `Frame3D`. **Etapa 2**: `Frame2DTimoshenko`, con su propia decisión de cinemática; mientras tanto sólo material elástico.
4. **Gauss-Legendre**, 3 puntos por defecto, mínimo 2.
5. **Torsión elástica desacoplada**, `T = GJ·θ′` con `GJ` declarado en la sección.
6. **Sin regularización** del ablandamiento en esta etapa; la regularización genérica de todos los materiales con ablandamiento queda como deuda.
7. Ocho casos de validación (spec).

Lo que sigue es la arquitectura que aloja esa física. La decide la IA (Reglas §3).

## Decisión

### 1. La sección es un objeto de familia, sin estado

`SectionRegistry` es una familia (ADR 0020, P1) con `YAML_SECTION = "sections"`, `YAML_LABEL = "sección"` y `SPEC_KIND = "section"`. Una sección se declara una vez y la comparten muchos elementos, como un material. **No guarda estado**: el estado de cada fibra en cada punto de Gauss vive en el elemento, igual que el de un material en un sólido. La base abstracta `Section` (`solidum/core/section.py`) fija el contrato que consumen los marcos:

- `response_2d(e, states) → (s, k_s, states')` con `e = [ε₀, κ_z]`, `s = [N, M_z]`;
- `response_3d(e, states) → (s, k_s, states')` con `e = [ε₀, κ_y, κ_z]`, `s = [N, M_y, M_z]`;
- `materials()`, `mass_properties()` (`ρA`, `ρI_y`, `ρI_z`) y el atributo `GJ` (opcional; lo exige `Frame3D`).

`FiberSection` es la primera implementación: fibras `(y, z, A, material)` referidas al **eje del elemento** (línea de nodos), generadas por parches rectangulares, parches circulares (macizos o huecos), capas de barras o fibras sueltas. Cualquier material de `STRAIN_DIM = 1` del catálogo sirve como material de fibra, y una sección puede mezclarlos.

### 2. Referencias entre familias: `from_yaml(data, resolve)`

Una sección referencia materiales dentro de sus parches. Los objetos de familia se construyen hoy con `cls(**datos)`; para las referencias anidadas, una clase de familia puede declarar el método de clase `from_yaml(data, resolve)`. El lector YAML lo llama si existe y le pasa `resolve(registro, id, contexto)`, que valida que el `id` exista en una familia ya construida y devuelve el objeto, con un mensaje que nombra el contexto si no. Las familias se construyen en orden de definición y `SectionRegistry` se define después de `MaterialRegistry`. El lector no conoce las secciones por nombre.

### 3. Los marcos: dos caminos, uno por construcción

Los tres marcos Euler-Bernoulli aceptan **exactamente uno** de:

- `material` + propiedades geométricas (`A, I` o `A, Iy, Iz, J`): la forma cerrada actual, **sin cambios de resultado** con material elástico. Desde este ADR sólo admite materiales con `IS_LINEAR_ELASTIC = True`; cualquier otro se rechaza al construir con un mensaje que remite a la sección por fibras. La misma guarda se aplica a `Frame2DTimoshenko` hasta la etapa 2.
- `section` (+ `n_ip`, 3 por defecto, mínimo 2): elemento en desplazamientos con integración de Gauss-Legendre (spec §5).

`REFERENCE_KWARGS = {"material": MaterialRegistry, "section": SectionRegistry}`, las dos opcionales (ADR 0020, P3): es la primera referencia del programa principal con nombre distinto de `material`. En el camino con sección, `N_INTEGRATION_POINTS` es atributo de instancia (`= n_ip`) antes de construir el `ElementState`; `vars[g]` guarda la lista de estados de fibra del punto `g` y `stresses[g]` el vector de fuerzas internas de la sección `s_g`.

### 4. Contratos genéricos nuevos

- `Material.IS_LINEAR_ELASTIC: ClassVar[bool] = False`; `True` en `Elastic1D`, `Elastic2D`, `Elastic3D` y `Orthotropic2D`. Declarativo (Reglas §1), lo recorre el barrido de contrato de materiales.
- `Element.materials()`: los materiales que usa el elemento (por defecto `[self.material]`). La prevalidación de masa y de peso propio del ensamblador y la simetría agregada de los solvers (ADR 0003) lo recorren en lugar de leer `element.material`. Caso real: los marcos con sección (un material por fibra). *Corrección del 2026-09-29, tras la implementación*: el texto original citaba como segundo caso el elemento embebido del módulo `discontinuities` (material del sólido y cohesivo). No se aplicó: el material cohesivo no tiene densidad propia (la inercia es del sólido), así que exigírsela en la prevalidación de masa sería incorrecto; ese elemento sigue devolviendo sólo el material del sólido. `Element.compute_self_weight(g)` (por defecto `compute_body_load(ρ·g)`) completa la generalización del peso propio sin casos especiales para secciones.

### 5. Masa y peso propio con sección

`ρA = Σρ_i A_i`, `ρI_z = Σρ_i y_i² A_i`, `ρI_y = Σρ_i z_i² A_i`, con las mismas fórmulas de masa consistente y concentrada de los marcos elásticos (la cuestión física abierta sobre su inercia rotatoria, hallazgo clas-4 de la auditoría, se resuelve aparte y aplica igual a ambos caminos). Una fibra sin `density` produce el error agregado de siempre.

### 6. Rendimiento

La sección evalúa sus fibras agrupadas por material. Los marcos no van por el camino por lotes del ADR 0014 (los elementos 1D quedaron diferidos allí). Se mide antes de vectorizar más (memoria: optimizar midiendo primero).

## Alternativas descartadas

- **Elementos nuevos** (`Frame2DEulerFiber`…): duplican cinemática, catálogo y tests de tres elementos para cambiar sólo la ley de sección.
- **La sección como material**: un material recibe una deformación escalar y devuelve un esfuerzo; una sección recibe `e = [ε₀, κ]` y devuelve fuerzas internas. Mezclar los contratos obligaría a casos especiales en el barrido de materiales y en el exportador.
- **Integrar la sección dentro de cada marco**: la geometría de fibras se repetiría por elemento y no podría compartirse ni declararse una vez en el YAML.
- **Elemento en fuerzas** (Spacone et al. 1996): decisión física del usuario; queda como variante futura si hace falta un elemento por miembro.

## Consecuencias

- **A favor**: los marcos dejan de dar flectores incorrectos en silencio; plasticidad distribuida en la sección con cualquier material 1D; secciones mixtas (concreto y acero); primera referencia con nombre distinto de `material` en el programa principal, que prueba el contrato P3 sin depender de un módulo de usuario.
- **En contra**: un modelo con marcos y material inelástico que hoy corre pasa a fallar al construir (antes daba resultados incorrectos); el ablandamiento en fibras depende de la malla hasta que exista la regularización genérica.
