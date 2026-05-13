# Diagnóstico de Problemas

- **`YamlValidationError: parámetro 'X' no aceptado por 'Y'.`**
  El parser detecta un parámetro que el constructor del elemento no admite. La salida lista los parámetros aceptados; revisar el catálogo o la spec del componente.

- **`Matriz singular detectada.`**
  Restricciones cinemáticas insuficientes: el dominio puede trasladarse o rotar como cuerpo rígido. Verificar que al menos un nodo tiene fijos los suficientes DOFs para anclar el dominio. En modelos con cables o materiales con softening pronunciado, comprobar que existen caminos de carga alternativos al elemento que se destensa o degrada.

- **`Newton-Raphson no convergió / bisecando incremento.`**
  El paso de carga es demasiado grande. Aumentar `num_steps`, asegurar `adaptive: true`. Verificar tipeo de los parámetros constitutivos ($\sigma_y$, $\kappa_0$, etc.). Si el problema tiene snap-through o softening, cambiar a `ArcLengthSolver`.

- **`No se importó ningún elemento (Solid2D).`**
  El parser de Gmsh no encontró superficies bidimensionales. En Gmsh, crear una *Plane Surface* y generar la malla 2D antes de exportar.

- **`Cable totalmente destensado — KT = 0.`**
  Pretensar el cable o garantizar que la estructura tiene rigidez residual por otros elementos. Considerar también pasos de carga finos para capturar la transición tensado → destensado.
