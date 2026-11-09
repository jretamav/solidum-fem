# ADR 0028 — Longitud característica del punto de integración

- **Estado**: aceptado (2026-11-09)
- **Fecha**: 2026-11-09
- **Alcance**:
  - `Material.uses_characteristic_length` y `Material.initial_point_state(h)`;
  - `Element._seed_characteristic_lengths(lengths)` y el atributo `characteristic_lengths`;
  - `Section.initial_states(characteristic_length=None)` (`FiberSection`);
  - las barras (`Truss2D`, `Truss3D` y sus corrotacionales) y los marcos con sección (`Frame2DEuler`, `Frame2DEulerCorot`, `Frame3D`).
- **Primer consumidor**: `ManderConcrete1D` con `G_f` o `G_fc` (spec `ManderConcrete1D` §5).

## Contexto

Un material con ablandamiento concentra la deformación en un punto de integración, y la energía que disipa es el área bajo su curva por la longitud de ese punto: sin corrección, el resultado depende de la malla. La regularización por energía de fractura (banda de fisuración, Bažant y Oh 1983; en marcos, Coleman y Spacone 2001) ajusta la rama descendente para que el punto disipe `G_f/h` por unidad de volumen. Para eso el material necesita la **longitud característica `h` de cada punto**, que sólo conoce el elemento (su longitud, el peso del punto de Gauss).

El contrato `Material.compute_state(strain, state_vars)` no la transporta, y un mismo objeto material lo comparten todos los puntos de todos los elementos que lo usan.

## Decisión

`h` viaja **en el estado del punto**, sembrado una sola vez por el elemento:

1. El material declara si la usa (`uses_characteristic_length`, propiedad de la instancia: el concreto sólo la usa si se dio la energía de fractura) y construye el estado inicial de un punto con `initial_point_state(h)`, donde valida `h` (p. ej. el límite de retroceso `h < 2·G_f·E/f_t²`) y lo guarda.
2. El elemento, al terminar de construirse (cuando ya conoce su geometría), llama a `_seed_characteristic_lengths(lengths)` con una longitud por punto de integración: `L` en una barra, `w_g·L` en un marco con sección. Si ningún material la usa, no cambia nada.
3. Con sección, `FiberSection.initial_states(characteristic_length=h)` reparte `h` a las fibras cuyo material la usa; las demás siguen con estado `None`.
4. El material arrastra `h` en cada `compute_state`. Si un elemento no la sembró y la rama regularizada se alcanza, el material lanza un error que nombra a los elementos que sí la proporcionan.

Un `h` inadmisible se rechaza **al construir el elemento**, con el nombre del elemento.

## Alternativas descartadas

- **Pasar `h` en cada `compute_state`** (argumento nuevo): obliga a tocar las más de treinta llamadas de todos los elementos y secciones, y la validación de `h` llegaría tarde, en la primera iteración.
- **Guardar `h` en el material**: un material lo comparten puntos de longitudes distintas.
- **Un registro global elemento → longitud**: el material no sabe a qué punto pertenece una llamada.

## Consecuencias

- El estado del punto lleva una constante más (`h`), declarada en el `STATE_SCHEMA` del material que la usa.
- Añadir la regularización a otro material (los de daño 2D/3D, con `h` del área o del volumen del elemento sólido) sólo exige que su elemento siembre las longitudes: la pieza ya es general.
- En los elementos de marco en desplazamientos, `h = w_g·L` es la elección de Coleman y Spacone para elementos en fuerzas. La prueba de objetividad de malla de la spec `ManderConcrete1D` mide su efecto; usar `h = L` da prácticamente la misma respuesta.
