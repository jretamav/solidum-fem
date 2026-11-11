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

## Un modelo YAML desde Python: `load_yaml`

`solidum.run_yaml(path)` lee el archivo, arma el modelo y lo resuelve, y devuelve sólo el resultado. Para conservar el modelo y consultarlo antes o después de resolver —los grados de libertad de cada nodo, la matriz de masa, la exportación VTK, el estado de los elementos— se usa `solidum.load_yaml(path)`:

```python
import solidum

modelo = solidum.load_yaml("modelo.yaml")   # YamlModel: domain, assembler, parser
res = modelo.solve()                          # lo mismo que run_yaml("modelo.yaml")
tope = modelo.domain.nodes[12]
print(res.U[tope.dofs["ux"]])
M = modelo.assembler.assemble_mass_matrix()
```

`modelo.solve(step_callback=f)` acepta el mismo `step_callback` que `run_yaml`: `f(step, U, λ)` en los estáticos y `f(step, U, t)` en los análisis en el tiempo. `modelo.domain` tiene los grados de libertad ya numerados y `modelo.assembler` es el ensamblador con el que `solve()` resuelve el análisis declarado (estático, modal, transitorio, armónico, espectral, térmico o por etapas). `modelo.parser` da las cargas del archivo (`get_external_forces()`, `get_body_load(assembler)`, `get_thermal_loads()`) y el solver declarado (`get_solver(assembler)`).

## Resultado de un análisis estático: `SolveResult`

`solidum.run(...)` (y `solidum.run_yaml(...)` con un solver estático) devuelve un `SolveResult` inmutable, que queda también en `domain.last_result`:

| Campo | Contenido |
|---|---|
| `U` | Solución global (desplazamientos y giros, o temperaturas), shape `(n_dof,)`. |
| `F_applied` | Cargas que ve el sistema (nodales más equivalentes de elemento), escaladas por el factor de carga final. |
| `R` | Reacciones de los **apoyos**, shape `(n_dof,)`: los grados de libertad con valor prescrito y las restricciones multipunto que actúan como apoyo (un rodillo oblicuo). Ceros en el resto. |
| `reactions_by_node` | `{nodo: {dof: valor}}` de `R`: los nodos con algún valor prescrito y los grados de libertad de las restricciones multipunto de apoyo. |
| `constraint_forces` | Fuerzas de restricción de las restricciones multipunto de **unión** (atado de mallas, periodicidad, unión rígida), shape `(n_dof,)`; ceros si no hay. |
| `constraint_forces_by_node` | `{nodo: {dof: valor}}` de `constraint_forces` en los grados de libertad de esas uniones. |
| `element_forces` | `{elemento: ElementForces}` con las fuerzas internas (N, V, M, T) de los elementos estructurales. |
| `converged`, `num_steps` | Si el análisis llegó al final y cuántos pasos convergidos dio. |
| `termination` | Por qué se detuvo un solver incremental (ADR 0025): `"max_lambda"`, `"control_target"` o `"max_steps"`. `None` en el `LinearSolver`. |
| `load_factor` | Factor de carga final λ de la carga que crece en el análisis o en la etapa. |
| `perturbation` | `True` en un `LinearSolver` sobre una etapa (ADR 0031): perturbación lineal que no cambió el estado; no se puede usar como `previous`. |
| `F_constant`, `F_variable` | En un análisis por etapas (`solidum.run(..., previous=…)`, ADR 0026), la carga ya aplicada que se mantuvo constante y la que creció (`F_applied = F_constant + F_variable`). Fuera de una etapa, `F_constant` es `None`. |
| `equivalent_loads` | Cargas equivalentes de la carga distribuida restadas, por elemento, para dar `element_forces`: las que lleva la carga aplicada (`LoadVector`, ADR 0027); una etapa siguiente las hereda con su carga constante. |
| `condition_estimate` | Número de condición estimado de la rigidez (escalada por la diagonal) en el análisis lineal: $\log_{10}$ de este valor acota las cifras que la solución puede haber perdido de unas 16; por encima de $10^{12}$ Solidum avisa (capítulo «Diagnóstico de Problemas»). `None` en los análisis no lineales y con el solver iterativo. |

Las reacciones y las fuerzas de restricción cumplen `F_int − F_applied = R + constraint_forces`: el equilibrio global es `ΣR + ΣF_applied = 0` por dirección, y las fuerzas de restricción suman cero porque sólo transmiten fuerza entre partes de la estructura. El criterio que separa una restricción multipunto de apoyo de una de unión se explica en la sección `linear_constraints` del capítulo «Archivo de Entrada YAML». En el problema térmico `R` es el flujo de calor nodal en los nodos con temperatura impuesta.

## Cargas distribuidas y sus cargas equivalentes

`assemble_self_weight` y `assemble_body_load` devuelven un `LoadVector` (ADR 0027): un vector de cargas que se usa como cualquier `np.ndarray`, pero que lleva las cargas nodales equivalentes de su parte distribuida en cada elemento. Las fuerzas internas de extremo descuentan las de la carga aplicada, y la aritmética de las cargas las combina:

```python
F_peso = asm.assemble_self_weight([0.0, -9.81])
F_viento = asm.assemble_body_load([2.0e3, 0.0])
F = 1.2 * F_peso + 1.6 * F_viento + F_puntual     # lleva 1.2·f_eq,peso + 1.6·f_eq,viento
res = solidum.run(dom, assembler=asm, F_applied=F)
```

Vale igual para dos casos de carga con el mismo ensamblador, para la carga de una etapa y para la función de carga de un transitorio (`F_func=lambda t: F_peso + a(t) * F_viento`, con `record_internal_forces=True`). Las cargas puntuales en un `np.ndarray` simple no aportan cargas equivalentes. Lo que **pierde** la información es convertir la carga con `np.asarray`/`np.array` o sumarla con `+=` sobre un arreglo de ceros: entonces Solidum supone, como antes, que contiene la última carga de cuerpo ensamblada, y avisa si se ensambló más de una.

## Secciones por fibras

El bloque `sections` del YAML tiene su equivalente directo en Python (ADR 0021): las primitivas generan las fibras y la sección se pasa al marco con `section`, en lugar de `material` y las propiedades geométricas.

```python
from solidum import Domain, Elastic1D, Elastoplastic1D, Frame2DEuler
from solidum.sections import FiberSection, RectangularPatch, BarLayer

concreto = Elastic1D(E=3.0e10, density=2400.0)
acero = Elastoplastic1D(E=2.0e11, sigma_y=4.2e8, H=2.0e9, density=7850.0)
sec = FiberSection(
    patches=[RectangularPatch(concreto, y=(-0.25, 0.25), z=(-0.15, 0.15), ny=40, nz=1),
             BarLayer(acero, start=(-0.21, -0.11), end=(-0.21, 0.11), n=3, area=5.1e-4,
                      displaces=concreto)],   # opcional: descuenta el concreto que ocupan
    GJ=None,                        # obligatorio sólo si la usa un Frame3D
)

domain = Domain()
n1 = domain.add_node(1, [0.0, 0.0])
n2 = domain.add_node(2, [3.0, 0.0])
domain.add_element(Frame2DEuler(1, [n1, n2], section=sec, n_ip=3))
```

`CircularPatch(material, r_outer=..., r_inner=0.0, center=(0.0, 0.0), n_rad=..., n_circ=...)` y `Fiber(y, z, area, material, displaces=None)` completan las primitivas. Con `displaces` (el mismo objeto material del parche que contiene la barra), la sección añade al final de sus arreglos una fibra de compensación de área negativa por barra, marcada en `sec.is_compensation` (ver `displaces` en el bloque `sections` del capítulo del YAML). La sección expone `sec.y`, `sec.z`, `sec.area` (arreglos por fibra), `sec.n_fibers`, `sec.materials()` y `sec.mass_properties()`; su respuesta a nivel de sección, `sec.response_2d([eps0, kappa_z], estados)` y `sec.response_3d([eps0, kappa_y, kappa_z], estados)`, devuelve las fuerzas internas, su tangente y los estados nuevos de las fibras (con `sec.initial_states()` como estado virgen), útil para trazar diagramas momento–curvatura sin montar un modelo. Tras resolver, `elem.compute_internal_forces(U)` añade `section_deformations` y `section_forces` en los puntos de Gauss.

## Módulos de usuario

La raíz del paquete sólo reexporta el programa principal. Las formulaciones no estándar viven en módulos de usuario (`solidum/user/<nombre>/`, ADR 0020) que `import solidum` no carga; en Python se cargan importándolos directamente o con `solidum.load_user_module`, el equivalente de `user_modules` en el YAML:

```python
import solidum
from solidum import Domain, Elastic2D
from solidum.materials.interface import CohesiveDamageIsotropic   # programa principal (ADR 0033)
from solidum.user.discontinuities import CST_Embedded2D

solidum.available_user_modules()                # ['discontinuities', 'dissipation_control']
solidum.load_user_module("discontinuities")     # equivalente; cargarlo dos veces no hace nada

bulk = Elastic2D(E=30e9, nu=0.2, hypothesis="plane_strain")
cohesivo = CohesiveDamageIsotropic(sigma_t0=2.5e6, G_f=100.0, K_e=1.0e13, softening="linear")
grieta = CST_Embedded2D(element_id=1, nodes=[n1, n2, n3], material=bulk,
                        cohesive_material=cohesivo, thickness=1.0)
```

`from solidum import CST_Embedded2D` ya no funciona: el elemento y `DiscontinuityState` se importan de `solidum.user.discontinuities`. La ley cohesiva es un material de interfaz del programa principal (`solidum.materials.interface`; clase base `solidum.core.interface_material.InterfaceMaterial`, registro `solidum.registry.InterfaceMaterialRegistry`).

## Masas concentradas en nudos

```python
domain.add_nodal_mass(5, ux=12.0e3, uy=12.0e3)   # masa [kg]
domain.add_nodal_mass(5, rz=4.0e3)               # inercia rotacional [kg·m²]
```

Equivale al bloque YAML `nodal_masses` (ADR 0024): un grado de libertad mecánico por argumento, valores finitos no negativos, y varias llamadas sobre el mismo grado de libertad se suman. Las masas quedan en `domain.nodal_masses` (`{(nodo, dof): valor}`) y el `Assembler` las suma a la diagonal de la matriz de masa. Son sólo inercia: no generan peso con la gravedad.

## Refuerzo embebido

`solidum.embed_reinforcement(domain, points, area, material, *, element=None, hosts=None, subdivisions=1, element_kwargs=None, bond=None)` embebe una barra en los sólidos del dominio (ADR 0032, spec `EmbeddedReinforcement`), con adherencia perfecta si `bond` es `None`, y devuelve un `EmbeddedBar`:

```python
bar = solidum.embed_reinforcement(domain, [[0.0, 0.05], [4.0, 0.05]], 1.6e-3,
                                  solidum.Elastic1D(E=200e9))
r = solidum.run(domain, F_applied=F)
r.element_forces[bar.element_ids[0]].components["N"]   # fuerza axial del primer tramo
bar.bond_forces(r)    # (n_nodos, ndim): fuerza de adherencia en cada nodo de la barra
```

**Adherencia-deslizamiento** (spec `EmbeddedReinforcement` §12): `bond=solidum.Bond(law, diameter, perimeter=None, integration="nodal")` (o el diccionario equivalente) une cada tramo de la barra a su anfitrión con una interfaz distribuida gobernada por una ley de interfaz de línea (`BondSlipModelCode`, `BondSlipMultilinear`); `perimeter` es π·d_b por omisión e `integration` es `"nodal"` (recomendada) o `"gauss"`. Las interfaces generadas quedan en `bar.interface_ids`, y `bar.bond_state(domain, r.U)` devuelve, por tramo, los puntos de integración con el deslizamiento `slip`, el esfuerzo de adherencia `bond_stress`, la apertura y tracción normales y la fuerza de adherencia por unidad de longitud `bond_force`. Con `bond`, `bond_forces` no aplica (no hay restricciones de embebido) y lanza `ValueError`:

```python
ley = solidum.BondSlipModelCode.model_code_2010(30e6, clear_rib_spacing=0.0105,
                                                 diameter=0.016, E_c=30e9)
bar = solidum.embed_reinforcement(domain, [[0.0, 0.05], [4.0, 0.05]], 2.01e-4,
                                  solidum.Elastic1D(E=200e9), element="Truss2D3",
                                  bond=solidum.Bond(ley, diameter=0.016))
r = solidum.run(domain, F_applied=F)
puntos = bar.bond_state(domain, r.U)          # [[InterfacePoint, ...] por tramo]
s_max = max(abs(p.slip) for tramo in puntos for p in tramo)
```

`EmbeddedBar` guarda los nodos y elementos de la barra, el anfitrión y las coordenadas naturales de cada nodo, y la longitud. La llamada va después de crear los elementos de sólido; si el dominio ya estaba numerado, se renumera. Los generadores se registran en `ConstraintRegistry`; la localización de puntos en elementos (`solidum.mesh.point_location`: `natural_coordinates`, `PointLocator`) sirve a cualquier sólido que declare `REFERENCE_CELL` y sus funciones de forma.

## Errores tipados

Desde un script conviene capturar los errores por su tipo, porque cada uno indica una causa y una salida distintas:

```python
from solidum.math.solvers.diagnostics import (
    MechanismError,          # el modelo no está suficientemente apoyado (antes de resolver)
    IllPosedSystemError,     # solución no fiable: mecanismo interno o matriz singular
    SolverDivergedError,     # Newton, arco y transitorios: el análisis se detuvo
    StepConstraintError,     # arco: la restricción del paso no tiene solución
    TimeIntegrationInstabilityError,  # diferencia central o θ < 1/2: paso inestable
    MasslessDofsError,       # modal o transitorio: grados de libertad sin masa
    Termination,             # motivos de parada: solver.termination == Termination.MAX_STEPS
)
from solidum.math.linalg import IterativeNotConvergedError   # sólo con linear_algebra: iterative

try:
    U = solver.solve(F_ext)
except MechanismError as err:
    print(err.motions)       # p. ej. ['traslación en la dirección x']
```

`MechanismError.motions` lista los movimientos libres ya descritos en palabras. Para el fracaso de un trazado (Newton o arco) se captura `SolverDivergedError`, no `RuntimeError`: `IllPosedSystemError`, el veredicto de la red de seguridad, también es un `RuntimeError`. Tras la excepción, `solver.termination == "diverged"` y `solver.U_final`, `lambda_final` y `steps_done` dicen hasta dónde llegó (ADR 0025). En un análisis en el tiempo, `exc.partial_result` (también `solver.partial_result`) es el resultado hasta el último paso convergido, y `solver.t_final` su instante; un transitorio que llega al final deja `solver.termination == "t_end"` y `result.termination == "t_end"`. `MasslessDofsError` (ADR 0024) lleva la causa en `cause` (en el modal, `too_many_modes`, `no_mass`, `massless_mechanism`, `singular_mass`, `indefinite_mass`; en los transitorios implícitos, las mismas salvo `too_many_modes`; en la diferencia central, `explicit_integration`), los grados de libertad sin masa en `massless_dofs` y cuántos tienen masa en `n_mass`. Los mensajes completos y cómo corregir cada caso están en el capítulo *Diagnóstico de Problemas*.

El solver del sistema lineal se puede fijar también desde Python con el mismo parámetro que en el YAML, p. ej. `LinearSolver(assembler, linear_algebra="iterative")`.
