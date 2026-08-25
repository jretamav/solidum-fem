# Quad4Thermal — cuadrilátero bilineal de conducción de calor 2D

> Spec pre-redactada por la IA (2026-08-25). **El usuario revisa la física; el plumbing está resuelto.** Los puntos abiertos de §Diálogo son físicos o de alcance, nunca arquitecturales.
>
> Segundo componente de la Etapa 8 — C1 núcleo mínimo. Primer elemento de la familia térmica. Consume [`ThermalConduction`](ThermalConduction.md).

---

## Especificación física

### 0. Descripción general

Cuadrilátero isoparamétrico de 4 nodos que resuelve **conducción de calor en 2D**. Misma geometría, mismas funciones de forma bilineales y misma cuadratura que el [`Quad4`](Quad4.md) mecánico; lo que cambia es el campo que interpola y la ecuación que discretiza.

**Un grado de libertad por nodo**: la temperatura $T$. `DOF_NAMES = ['T']`.

Es una **clase distinta** del `Quad4` mecánico, no un modo de operación suyo. `DOF_NAMES` es atributo de clase y el registro, la validación temprana y el cruce spec↔código dependen de poder leerlo sin construir el objeto; un elemento que decidiera sus DOFs según el material recibido rompería ese contrato declarativo. Es también la separación que adoptan los códigos de referencia (Abaqus `CPE4` vs `DC2D4`; ANSYS `PLANE182` vs `PLANE55`).

### 1. Interpolación del campo

$$T(\xi,\eta) = \sum_{i=0}^{3} N_i(\xi,\eta)\,T_i, \qquad N_i = \tfrac14(1+\xi_i\xi)(1+\eta_i\eta)$$

Las $N_i$ son **las mismas** del `Quad4` mecánico. Isoparamétrico: geometría y temperatura comparten interpolación.

### 2. Gradiente de temperatura

$$\nabla T = \mathbf B\,\mathbf T_e, \qquad \mathbf B = \begin{bmatrix} \partial N_0/\partial x & \cdots & \partial N_3/\partial x \\[2pt] \partial N_0/\partial y & \cdots & \partial N_3/\partial y \end{bmatrix} \in \mathbb R^{2\times 4}$$

**Diferencia esencial con el mecánico**: aquí $\mathbf B$ es literalmente el gradiente de las funciones de forma, sin reordenamiento. En el `Quad4` mecánico las mismas derivadas $\partial N_i/\partial x_j$ se reordenan en una matriz $3\times 8$ que produce la deformación en notación Voigt. La materia prima es idéntica —el kernel `_compute_kinematics` ya la calcula—; el ensamblaje difiere.

Consecuencia: el gradiente $\nabla T$ es **bilineal por elemento** y su valor en cada punto de Gauss no es constante, a diferencia del CST.

### 3. Ecuación discretizada — forma débil

Partiendo del balance $\rho c\,\dot T = \nabla\cdot(\mathbf k\nabla T) + Q$, multiplicando por una función de peso $w$ e integrando por partes:

$$\int_\Omega \rho c\,w\,\dot T\,d\Omega + \int_\Omega (\nabla w)^\top \mathbf k \nabla T\,d\Omega = \int_\Omega w\,Q\,d\Omega - \int_{\partial\Omega} w\,\bar q\,d\Gamma$$

que en forma matricial es el sistema semidiscreto

$$\mathbf C\,\dot{\mathbf T} + \mathbf K\,\mathbf T = \mathbf F$$

**de primer orden en el tiempo** — a diferencia del mecánico $\mathbf M\ddot{\mathbf u} + \mathbf C\dot{\mathbf u} + \mathbf K\mathbf u = \mathbf F$, que es de segundo. Por eso Newmark no aplica y el transitorio requiere solver propio (spec separada).

En **régimen estacionario** el término $\mathbf C\dot{\mathbf T}$ se anula y queda $\mathbf K\mathbf T = \mathbf F$, que resuelve el `LinearSolver` existente sin modificación.

### 4. Matriz de conductividad

$$\mathbf K_e = \int_\Omega \mathbf B^\top\,\mathbf k\,\mathbf B\;t\;d\Omega = \sum_{g} \mathbf B_g^\top\,\mathbf k\,\mathbf B_g\;t\;|\mathbf J_g|\,w_g$$

$4\times4$, **simétrica** (por serlo $\mathbf k$) y **semidefinida positiva**: tiene un modo nulo, el campo de temperatura uniforme $\mathbf T = c\,[1,1,1,1]^\top$, que no produce flujo. Es el análogo térmico exacto de los modos de sólido rígido en mecánica, y como allí, el sistema global sólo se vuelve resoluble al imponer al menos una condición de Dirichlet.

El **espesor** $t$ interviene igual que en el mecánico: el problema plano representa una rebanada de profundidad $t$, y el flujo total escala con ella.

### 5. Matriz de capacidad

$$\mathbf C_e = \int_\Omega \rho c\;\mathbf N^\top\mathbf N\;t\;d\Omega$$

$4\times4$, simétrica y **definida positiva**. Estructuralmente idéntica a la matriz de masa consistente con $\rho c$ en lugar de $\rho$; de hecho $\mathbf C_e = \rho c\,\mathbf M_e^{escalar}/\rho$ sobre el mismo integrando.

**Por decisión de alcance, el default es `lumped`, al revés que en dinámica estructural.** La razón es física, no de conveniencia:

- La capacidad **consistente** produce **oscilaciones espurias** en los primeros pasos cuando existe un frente térmico abrupto (un salto de temperatura impuesto en la frontera). El campo calculado oscila alrededor de la solución y puede dar temperaturas **fuera del rango de los datos** —por debajo del mínimo o por encima del máximo impuestos—, lo que viola el **principio del máximo** de la ecuación de difusión: una propiedad matemática de la solución exacta, no un artefacto tolerable.
- La capacidad **lumped** amortigua esas oscilaciones y preserva la monotonía, a cambio de algo de precisión en régimen suave.

Ambas quedan disponibles con el mismo parámetro `lumping` del resto del catálogo. La coherencia arquitectural se preserva en el **mecanismo**; el default lo fija el criterio físico del dominio (Lewis-Nithiarasu-Seetharamu §6; Zienkiewicz-Taylor Vol.1 §7).

### 6. Condiciones de frontera y cargas

**Dirichlet — temperatura impuesta.** $T = \bar T$ en $\Gamma_T$. Se impone por eliminación directa, con la maquinaria existente del ADR 0004, sin cambios: el mecanismo es agnóstico al significado del DOF.

**Neumann — flujo impuesto.** $\bar q$ [W/m²] normal a un borde. El vector consistente sobre el borde $\Gamma_q$:

$$\mathbf f^{\bar q}_e = -\int_{\Gamma_q} \bar q\,\mathbf N^\top\,t\;d\Gamma$$

Para $\bar q$ uniforme sobre un borde recto de longitud $L$, la integral da $\bar q\,L\,t/2$ a cada uno de los dos nodos del borde, y cero a los otros dos — el mismo reparto mitad-mitad que `compute_edge_traction` produce en el mecánico, con la misma numeración de bordes `EDGE_NODES = ((0,1),(1,2),(2,3),(3,0))`.

**Convención de signo**: $\bar q > 0$ significa flujo **saliente** del dominio (enfriamiento), coherente con que la normal $\mathbf n$ es saliente y con la deducción de la forma débil, donde el término de frontera aparece con signo negativo. Un flujo entrante se declara con $\bar q < 0$.

Un borde sin condición declarada queda con $\bar q = 0$ — **frontera adiabática**, aislada. Es el default natural de Neumann homogéneo, exactamente como un borde mecánico sin tracción queda libre.

**Fuente volumétrica.** $Q$ [W/m³] se aplica como carga del elemento (decisión de alcance; ver §Diálogo de [`ThermalConduction`](ThermalConduction.md)):

$$\mathbf f^{Q}_e = \int_\Omega Q\;\mathbf N^\top\,t\;d\Omega$$

Para $Q$ uniforme la suma de las componentes nodales es exactamente $Q\,A_e\,t$, invariante ante distorsión de la geometría — el mismo criterio de verificación que se usa en `compute_body_load` del mecánico.

### 7. Salida por punto de Gauss

`compute_gauss_state(T_global)` devuelve, en paralelo al contrato de los sólidos (ADR 0012):

- `points_natural`, `points_global` — posición de los puntos de Gauss.
- `grad_T` — gradiente $(n_g, 2)$.
- `flux` — vector de flujo $\mathbf q = -\mathbf k\nabla T$, $(n_g, 2)$.

No expone `internal_forces`: por el cierre por dominio del ADR 0012, ese método corresponde a elementos estructurales 1D. La familia térmica sigue el mismo criterio que los sólidos.

---

## Formulación numérica

### 8. Cuadratura

Gauss $2\times2$ (4 puntos) por defecto, la misma del `Quad4` mecánico, mediante `QuadratureRegistry`. Es **exacta** para $\mathbf K_e$ y $\mathbf C_e$ con jacobiano constante (elemento paralelogramo).

La reducción $1\times1$ queda disponible por el mismo parámetro, pero **no se recomienda**: al igual que en el mecánico introduce modos espurios. En el caso térmico el rango de $\mathbf K_e$ con un punto es 2, frente a los 3 necesarios (4 DOFs menos el modo de temperatura uniforme), de modo que aparece **un modo hourglass** de temperatura no detectado por la matriz. Se documenta y se blinda con test de conteo de modos nulos.

### 9. Reutilización de kernels

`_shape_functions_quad4`, `_det_jacobian_quad4` y las derivadas de `_compute_kinematics` ya existen compiladas con Numba en `solidum/elements/solid_2d/_shared.py` y se reutilizan sin duplicar. Lo específico térmico es el ensamblaje $\mathbf B^\top\mathbf k\,\mathbf B$, que es más simple que el mecánico.

### 10. Caveats numéricos

- **Modo nulo de temperatura uniforme**: $\mathbf K$ global es singular sin Dirichlet. El sistema debe tener al menos un nodo con temperatura impuesta, o el solver falla con matriz singular. Es el análogo exacto de un modelo mecánico sin apoyos, y el diagnóstico debe nombrarlo como tal.
- **Sin locking**: el problema escalar no tiene el análogo de la incompresibilidad ni del bloqueo por cortante. No aplica ninguna de las limitaciones que arrastra el `Quad4` mecánico.
- **Número de Péclet**: sin término convectivo en la ecuación (fuera de alcance C1), no hay riesgo de oscilaciones por advección dominante. Se anota porque será la primera consideración a revisar si en el futuro entra convección forzada del medio.
- **Distorsión de la malla**: como en el mecánico, un jacobiano casi degenerado degrada la solución; se valida $|\mathbf J| > 0$ en todos los puntos de Gauss.

---

## Contrato de implementación

```yaml
name: Quad4Thermal
kind: element
status: validated        # draft → implemented → validated

interface:
  field: temperature            # campo primario (Etapa 8) — no displacement
  dof_names: [T]
  n_nodes: 4
  flux_dim: 2                   # dimensión de ∇T y de q en 2D (no hay strain_dim)
  n_integration_points: 4       # default Gauss 2×2

parameters:
  - { name: thickness,  type: float, required: false, default: 1.0,
      desc: "Espesor de la rebanada plana [m]; el flujo total escala con él" }
  - { name: quadrature, type: str,   required: false, default: "2x2",
      desc: "Regla desde QuadratureRegistry. 1x1 disponible pero introduce
             un modo hourglass de temperatura — sin estabilización" }

material_contract:
  signature: "compute_flux(∇T) -> (q, k)"
  gradient_kind: "gradiente físico [∂T/∂x, ∂T/∂y] — sin notación Voigt"
  accepted: [ThermalConduction]

conventions:
  sign_flux: "q̄ > 0 ⇔ flujo SALIENTE del dominio (enfriamiento); normal exterior"
  adiabatic_default: "borde sin condición declarada ⇒ q̄ = 0 (aislado)"
  node_orientation: "antihorario (asegura det J > 0) — idéntico al Quad4 mecánico"
  edge_numbering: "EDGE_NODES = ((0,1),(1,2),(2,3),(3,0)) — paritario con Quad4"
  capacity_default: "lumping='lumped' por default, al contrario que en dinámica
                     estructural; razón física en §5"

validity:
  - "conducción lineal; k independiente de T"
  - "material con field='temperature' (familia ThermalMaterial)"
  - "det J > 0 en todos los puntos de Gauss"
  - "el modelo global debe imponer al menos un Dirichlet de temperatura"

out_of_scope:
  - "convección (Robin) y radiación en la frontera — Etapa 8 posterior"
  - "acoplamiento termomecánico y deformación térmica"
  - "conducción con advección (término convectivo del medio)"
  - "cambio de fase / calor latente"
  - "estabilización del modo hourglass con cuadratura reducida"

acceptance:
  verification:
    - name: patch_test_termico_lineal
      setup: "5 Quad4Thermal distorsionados (4 exteriores + 1 interior);
             imponer T = a + b·x + c·y en los nodos del contorno"
      expect: "los nodos interiores adoptan el campo lineal exacto y ∇T =
              (b, c) constante en TODOS los puntos de Gauss; q = -k·(b,c)"
      tol_rel: 1.0e-10
    - name: modo_nulo_temperatura_uniforme
      setup: "K_e de un elemento aislado; campo T = [1,1,1,1]"
      expect: "K_e·T = 0 exacto; rango(K_e) = 3; exactamente un autovalor nulo"
      tol_abs: 1.0e-12
    - name: simetria_y_semidefinicion
      setup: "K_e y C_e con k isótropo y ortótropo, elemento regular y distorsionado"
      expect: "K_e == K_e.T y C_e == C_e.T exactas; autovalores de K_e ≥ 0 con un
               único cero; autovalores de C_e > 0 estrictos"
      tol_abs: 1.0e-12
  specific:
    - name: conduccion_1d_pared_plana
      setup: "franja de Quad4Thermal entre x=0 y x=L; T(0)=T_1, T(L)=T_2;
             bordes superior e inferior adiabáticos; estacionario"
      expect: "perfil lineal exacto T(x) = T_1 + (T_2-T_1)·x/L en todos los nodos;
              flujo uniforme q_x = k·(T_1-T_2)/L en todos los Gauss"
      tol_rel: 1.0e-12
    - name: cilindro_hueco_logaritmico
      setup: "corona circular r ∈ [a, b] mallada con Quad4Thermal;
             T(a)=T_a, T(b)=T_b; estacionario"
      expect: "perfil logarítmico analítico T(r) = T_a + (T_b-T_a)·ln(r/a)/ln(b/a);
              convergencia con refinamiento h. Es el análogo térmico exacto del
              benchmark de Lamé que el proyecto ya usa en mecánica"
      tol_rel: 5.0e-3
    - name: fuente_volumetrica_uniforme
      setup: "Quad4Thermal regular y distorsionado con Q uniforme; sumar
             componentes nodales de f^Q"
      expect: "Σf = Q·A_e·t exacto, invariante ante distorsión de la geometría"
      tol_rel: 1.0e-12
    - name: flujo_neumann_en_un_borde
      setup: "Quad4Thermal con q̄ constante sobre un borde"
      expect: "Σf = -q̄·L·t con reparto L/2 a cada nodo del borde y 0 en los demás"
      tol_rel: 1.0e-12
    - name: balance_energetico_estacionario
      setup: "dominio con fuente Q y Dirichlet en parte de la frontera; estacionario"
      expect: "el calor generado iguala al evacuado por las reacciones nodales:
              Σ Q·A·t + Σ reacciones = 0"
      tol_rel: 1.0e-10
    - name: capacidad_lumped_vs_consistente
      setup: "C_e en ambos modos para elemento regular"
      expect: "ambas suman ρ·c·A_e·t (capacidad total conservada); la lumped es
              diagonal y estrictamente positiva"
      tol_rel: 1.0e-12
    - name: anisotropia_desvia_el_flujo
      setup: "k = diag(k1, k2) con k1 ≠ k2; campo lineal a 45° de los ejes"
      expect: "q no es antiparalelo a ∇T; coincide con -k·∇T calculado a mano"
      tol_rel: 1.0e-12
    - name: jacobiano_degenerado_abortado
      setup: "elemento con nodos colineales o en orden inverso"
      expect: "ValueError con mensaje claro"
    - name: rechazo_material_mecanico
      setup: "construir Quad4Thermal con un material de la familia mecánica (Elastic2D)"
      expect: "TypeError/ValueError identificando la incompatibilidad de familia;
               nunca aceptación silenciosa"

references:
  - "Lewis R.W., Nithiarasu P., Seetharamu K.N. (2004). Fundamentals of the FEM for Heat and Fluid Flow. Wiley. §3 y §6 (capacidad lumped vs consistente, oscilaciones)."
  - "Zienkiewicz O.C., Taylor R.L. (2000). The Finite Element Method, Vol. 1. §7 (problemas de campo escalar)."
  - "Reddy J.N., Gartling D.K. (2010). The Finite Element Method in Heat Transfer and Fluid Dynamics. CRC Press. §4."
  - "Carslaw H.S., Jaeger J.C. (1959). Conduction of Heat in Solids. Oxford. §7.2 (cilindro hueco, solución logarítmica)."
  - "Bathe K.-J. (2014). Finite Element Procedures. §7 (formulación de campo escalar; analogía con elasticidad)."
```

---

## Implementación

- Archivo: [`solidum/elements/thermal/quad4_thermal.py`](../../solidum/elements/thermal/quad4_thermal.py)
- Clase: `Quad4Thermal`, registrada en `ElementRegistry`
- Base compartida: [`solidum/elements/thermal/_shared.py`](../../solidum/elements/thermal/_shared.py) · `_ThermalSolid`
- Tests: [`tests/test_quad4_thermal.py`](../../tests/test_quad4_thermal.py) — 31 verdes

**Notas de implementación:**

- **Base `_ThermalSolid` introducida con el primer elemento**, no con el segundo. No contradice la regla de los dos casos reales: el segundo caso —`Hex8Thermal`, de esta misma etapa— está especificado y difiere **únicamente en geometría**. La ecuación discretizada $\mathbf C\dot{\mathbf T} + \mathbf K\mathbf T = \mathbf F$, la validación del material, la capacidad, la fuente y el post-proceso son idénticos en 2D y 3D; la dimensión sólo cambia el tamaño de $\mathbf B$. La subclase concreta aporta funciones de forma, gradiente, cuadratura por defecto y topología de frontera.

- **`STRAIN_DIM = None` y validación propia.** `Element._validate_material_compatibility` compara `STRAIN_DIM` del elemento contra el del material; un elemento térmico no tiene deformación. Se sobreescribe el método para comprobar **familia** (`isinstance(mat, ThermalMaterial)`, con mensaje que explica la diferencia de contrato) y **dimensión** (`FLUX_DIM`). Un material mecánico pasado por error se detecta al construir, no más tarde con un `AttributeError` sobre `compute_flux`.

- **`_init_state` devuelve `None`.** La conducción de Fourier no tiene variables internas, así que no se crea `ElementState`: no hay historia que promover trial → committed. Cuando entre un material térmico con memoria, se reimplementa con el caso real delante.

- **Kernel `_compute_gradient_kinematics_quad4` añadido a `solid_2d/_shared.py`**, no a la familia térmica. `_compute_kinematics` ya calculaba `dN_dx` pero sólo devolvía la `B` mecánica ensamblada en Voigt; el nuevo kernel expone el gradiente crudo, que **es** la `B` térmica. Vive junto a su gemelo mecánico porque comparte el jacobiano y su inversión: separarlos habría duplicado esa aritmética.

- **`lump_hrz` con `total_mass = ρc·V_e`.** El campo es escalar (un DOF por nodo), así que `n_translational_dirs=1`. El esquema escala la diagonal para conservar la capacidad total, que es la propiedad que importa: la energía almacenable no debe depender de cómo se reparta entre nodos. Verificado en test para ambos modos.

- **`compute_edge_flux` devuelve valores negativos para $\bar q > 0$.** El signo sale de la forma débil, donde el término de frontera aparece como $-\int w\,\bar q\,d\Gamma$: un flujo saliente positivo **extrae** energía. Blindado con test explícito de signo en ambos sentidos.

- **Aliases `compute_global_stiffness` y `compute_mass_matrix`** delegan en conductividad y capacidad, de modo que el elemento encaja en el ensamblador y el pipeline existentes sin ramas especiales.

**Validación end-to-end confirmada**: el problema de pared plana resuelto con `Assembler` + `LinearSolver` reproduce el perfil lineal analítico con error máximo $1.1\times10^{-13}$, sin ninguna modificación del ensamblador ni de la imposición de Dirichlet (ADR 0004). Confirma que la infraestructura era genuinamente agnóstica al tipo de campo.

**Pendiente de esta spec — benchmark del cilindro hueco.** El caso `cilindro_hueco_logaritmico` del bloque `acceptance` **no está implementado todavía**: requiere mallar una corona circular, y se agrupa con la campaña de validación del cierre de etapa junto al `balance_energetico_estacionario`. Los tests entregados cubren el patch test, el modo nulo, las cargas, la capacidad en ambos modos y la pared plana end-to-end. El estado `validated` refleja que la formulación está verificada contra solución analítica; los dos casos restantes añaden cobertura sobre geometría curva y conservación global, no verifican la formulación básica.

---

## Diálogo

*Puntos abiertos para el usuario. Sólo física y alcance.*

**2026-08-25 — IA → usuario. Dos puntos:**

1. ~~**Convención de signo del flujo de frontera.**~~ ✅ **Resuelto 2026-08-25 — $\bar q > 0$ es flujo SALIENTE del dominio** (decisión del usuario).

   Es la que resulta de la deducción de la forma débil con normal exterior $\mathbf n$, donde el término de frontera aparece con signo negativo, y coincide con Abaqus y con la mayoría de la literatura de MEF. Un flujo entrante (calentamiento) se declara con $\bar q < 0$.

   **Consecuencia**: es convención del proyecto y debe registrarse en `Reglas.md §5` junto a las de signos mecánicas, al cerrar la etapa. Coherente además con `q = -k·∇T` del material: ambos signos salen de la misma deducción, no son elecciones independientes.

2. ~~**Espesor `thickness`.**~~ ✅ **Resuelto 2026-08-25 — se mantiene, con `default = 1.0`** (decisión del usuario).

   Misma semántica que en el `Quad4` mecánico: el problema plano representa una rebanada de profundidad $t$ y el flujo total escala con ella. Tres razones:

   - **Consistencia** con los cinco sólidos planos del catálogo, que ya lo declaran.
   - **Balances energéticos en unidades absolutas**: el test `balance_energetico_estacionario` da watts reales, no "watts por metro de profundidad".
   - **Extensión futura**: si entra el acoplamiento termomecánico, el elemento térmico y el mecánico deben compartir profundidad; tenerlo desde el principio evita una asimetría posterior.

   El coste es nulo: quien no lo declare obtiene exactamente el comportamiento "por unidad de profundidad" que asume la literatura térmica. Omitirlo habría sido particularizar al caso de hoy.
