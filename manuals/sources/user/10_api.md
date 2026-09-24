# Uso Avanzado: API Programática

Aunque el flujo principal es *data-driven* (YAML), el motor se puede usar como librería Python para escenarios donde el archivo no es suficiente: optimización paramétrica, generación procedural de mallas, acoplamiento con bibliotecas externas, scripting de batería de tests.

## Patrón Fachada

Los componentes principales se exponen desde la raíz del paquete:

```python
import numpy as np
from solidum import Domain, VonMises2D, Quad4, Assembler, NonlinearSolver, VtkExporter

# 1. Dominio y nodos
domain = Domain()
n1 = domain.add_node(1, [0.0, 0.0])
n2 = domain.add_node(2, [1.0, 0.0])
n3 = domain.add_node(3, [1.0, 1.0])
n4 = domain.add_node(4, [0.0, 1.0])

# 2. Material y elemento
acero = VonMises2D(E=200e9, nu=0.3, sigma_y=250e6, H=2e9, hypothesis='plane_strain')
placa = Quad4(element_id=1, nodes=[n1, n2, n3, n4], material=acero, thickness=0.1)
domain.add_element(placa)

# 3. Restricciones (Dirichlet)
n1.fix_dof('ux', 0.0); n1.fix_dof('uy', 0.0)
n4.fix_dof('ux', 0.0); n4.fix_dof('uy', 0.0)

domain.generate_equation_numbers()

# 4. Cargas (Neumann)
F_ext = np.zeros(domain.total_dofs)
F_ext[n2.dofs['ux']] = 150000.0
F_ext[n3.dofs['ux']] = 150000.0

# 5. Resolución
assembler = Assembler(domain)
solver = NonlinearSolver(assembler, num_steps=10, adaptive=True)
U = solver.solve(F_ext)

# 6. Exportación
VtkExporter(domain).export("resultados.vtu", U=U, F_ext=F_ext)
```

`Assembler(domain)` ensambla **por lotes** (ADR 0014): agrupa los elementos que comparten clase, material y cuadratura en familias y los evalúa en un kernel compilado, con el mismo resultado que el bucle elemento a elemento a precisión de máquina. `Assembler(domain, batch=False)` fuerza el camino por elemento en todo el modelo (útil para comparar o diagnosticar); `parallel=False` usa el kernel serie en vez del paralelo (el resultado es bit a bit el mismo; el número de hilos lo fija Numba con la variable de entorno `NUMBA_NUM_THREADS` o `numba.set_num_threads`); `batch_memory_budget=` acota, en bytes, los temporales por trozo. Los elementos sin kernel por lotes (1D) y los que guardan estado propio fuera de `ElementState` (como el `CST_Embedded2D` del módulo de usuario `discontinuities`) se ensamblan como siempre dentro del mismo `Assembler`. Si se sustituye el material o el espesor de un elemento después del primer ensamblaje, llamar a `assembler.invalidate()`. Un elemento con nodos en orden incorrecto (jacobiano negativo) detiene el ensamblaje con `ValueError` que nombra su `id`.

Para obtener deformaciones y esfuerzos en los puntos de Gauss de todo el modelo, `solidum.math.batch.gauss_states(domain, U)` devuelve `{id: compute_gauss_state(U)}` evaluando cada familia de una vez (las entradas son vistas sobre arreglos) y llamando a `compute_gauss_state` sólo en los elementos sin familia; el exportador VTK hace lo mismo por dentro.

El exportador evalúa los esfuerzos en el `U` que recibe (`compute_gauss_state`), así que funciona igual tras `LinearSolver.solve` directo, tras `solidum.run` o en un paso intermedio del `step_callback`. Si el modelo lleva peso propio o fuerza de cuerpo, usar `solidum.run(...)` en vez de `solver.solve(...)` para que `SolveResult.element_forces` reste la carga nodal equivalente y devuelva fuerzas internas de extremo coherentes con las reacciones.

## Módulos de usuario

La raíz del paquete sólo reexporta el programa principal. Las formulaciones no estándar viven en módulos de usuario (`solidum/user/<nombre>/`, ADR 0020) que `import solidum` no carga; en Python se cargan importándolos directamente o con `solidum.load_user_module`, el equivalente de `user_modules` en el YAML:

```python
import solidum
from solidum import Domain, Elastic2D
from solidum.user.discontinuities import CST_Embedded2D, CohesiveDamageIsotropic

solidum.available_user_modules()                # ['discontinuities']
solidum.load_user_module("discontinuities")     # equivalente; cargarlo dos veces no hace nada

bulk = Elastic2D(E=30e9, nu=0.2, hypothesis="plane_strain")
cohesivo = CohesiveDamageIsotropic(sigma_t0=2.5e6, G_f=100.0, K_e=1.0e13, softening="linear")
grieta = CST_Embedded2D(element_id=1, nodes=[n1, n2, n3], material=bulk,
                        cohesive_material=cohesivo, thickness=1.0)
```

`from solidum import CST_Embedded2D` ya no funciona: el elemento, la ley cohesiva (`CohesiveDamageIsotropic`), su clase base (`CohesiveMaterial`), su registro (`CohesiveMaterialRegistry`) y `DiscontinuityState` se importan de `solidum.user.discontinuities`.

## Errores tipados

Desde un script conviene capturar los errores por su tipo, porque cada uno indica una causa y una salida distintas:

```python
from solidum.math.solvers.diagnostics import (
    MechanismError,          # el modelo no está suficientemente apoyado (antes de resolver)
    IllPosedSystemError,     # solución no fiable: mecanismo interno o matriz singular
    SolverDivergedError,     # familia del Newton: oscilación, tangente singular, capacidad...
)
from solidum.math.linalg import IterativeNotConvergedError   # sólo con linear_algebra: iterative

try:
    U = solver.solve(F_ext)
except MechanismError as err:
    print(err.motions)       # p. ej. ['traslación en la dirección x']
```

`MechanismError.motions` lista los movimientos libres ya descritos en palabras. Los mensajes completos y cómo corregir cada caso están en el capítulo *Diagnóstico de Problemas*.

El solver del sistema lineal se puede fijar también desde Python con el mismo parámetro que en el YAML, p. ej. `LinearSolver(assembler, linear_algebra="iterative")`.
