# Hex8Thermal — hexaedro trilineal de conducción de calor 3D

> Spec pre-redactada por la IA (2026-08-25). **Spec de extensión** (`Reglas.md §4`): documenta **sólo lo que cambia** respecto a [`Quad4Thermal`](Quad4Thermal.md), que es el hermano 2D y donde vive la formulación completa de la familia térmica. Sin ADR propio: reusa las decisiones de alcance de la Etapa 8.
>
> Tercer componente de la Etapa 8 — C1 núcleo mínimo. Consume [`ThermalConduction`](ThermalConduction.md).

---

## Qué comparte con `Quad4Thermal`

Todo lo esencial de la formulación, que **no se repite aquí**:

- Un DOF por nodo, la temperatura $T$. `DOF_NAMES = ['T']`.
- Sistema semidiscreto $\mathbf C\dot{\mathbf T} + \mathbf K\mathbf T = \mathbf F$, de primer orden en el tiempo. Estacionario $\mathbf K\mathbf T = \mathbf F$ con el `LinearSolver` existente.
- $\mathbf K_e = \int \mathbf B^\top\mathbf k\,\mathbf B\,d\Omega$ simétrica y **semidefinida positiva**, con el modo nulo de temperatura uniforme.
- $\mathbf C_e = \int \rho c\,\mathbf N^\top\mathbf N\,d\Omega$, con **`lumped` por defecto** y la justificación física del principio del máximo.
- Convención de signo: $\bar q > 0$ es flujo **saliente** del dominio; frontera sin declarar es **adiabática**.
- Fuente $Q$ como carga del elemento.
- Sin locking de ningún tipo; sin `internal_forces` (cierre por dominio del ADR 0012).
- Clase separada del `Hex8` mecánico, por el mismo argumento del `ClassVar DOF_NAMES`.

---

## Qué cambia respecto a `Quad4Thermal`

### 1. Dimensión del gradiente y de la matriz B

$$\nabla T = \mathbf B\,\mathbf T_e, \qquad \mathbf B = \begin{bmatrix} \partial N_0/\partial x & \cdots & \partial N_7/\partial x \\[2pt] \partial N_0/\partial y & \cdots & \partial N_7/\partial y \\[2pt] \partial N_0/\partial z & \cdots & \partial N_7/\partial z \end{bmatrix} \in \mathbb R^{3\times 8}$$

`flux_dim = 3`. La conductividad $\mathbf k$ es $3\times3$; el material `ThermalConduction` la expande a $k\mathbf I_3$ cuando recibe un escalar.

$\mathbf K_e$ y $\mathbf C_e$ son $8\times8$. El modo nulo sigue siendo único (temperatura uniforme), luego $\text{rango}(\mathbf K_e) = 7$.

### 2. Sin espesor

**`Hex8Thermal` no lleva `thickness`.** El volumen sale de la geometría, igual que en los sólidos 3D mecánicos. Es la diferencia visible más inmediata para el usuario del YAML respecto al hermano 2D.

### 3. Interpolación trilineal

$$N_i(\xi,\eta,\zeta) = \tfrac18(1+\xi_i\xi)(1+\eta_i\eta)(1+\zeta_i\zeta)$$

Las mismas del `Hex8` mecánico, orden de nodos VTK_HEXAHEDRON. El gradiente resulta **trilineal** por elemento.

### 4. Cuadratura

Gauss $2\times2\times2$ (8 puntos) por defecto — `hex_2x2x2` del `QuadratureRegistry`, idéntica al `Hex8` mecánico. Exacta para $\mathbf K_e$ y $\mathbf C_e$ con jacobiano constante.

Reducida `hex_1x1x1` disponible por el mismo parámetro pero **no recomendada**: con un solo punto el rango de $\mathbf K_e$ es 3 frente a los 7 necesarios, dejando **4 modos hourglass** de temperatura sin detectar (frente al único del `Quad4Thermal` reducido). Se blinda con test de conteo.

### 5. Flujo de frontera por cara, no por borde

El flujo Neumann se aplica sobre **caras** con normal saliente, no sobre bordes:

$$\mathbf f^{\bar q}_e = -\int_{\Gamma_q} \bar q\,\mathbf N^\top\,d\Gamma$$

Numeración **paritaria con el `Hex8` mecánico** — misma `FACE_NODES`, mismas 6 caras: 0 ($-\zeta$), 1 ($+\zeta$), 2 ($-\eta$), 3 ($+\xi$), 4 ($+\eta$), 5 ($-\xi$). Reutiliza la maquinaria de integración sobre caras que ya existe.

Para $\bar q$ uniforme sobre una cara plana de área $A$, el reparto es $\bar q A/4$ a cada uno de los 4 nodos de la cara y cero a los otros cuatro.

### 6. Validación específica 3D

Los benchmarks del hermano 2D se extienden, más uno propio de la dimensión:

- **Pared plana 3D**: perfil lineal exacto; verifica que la tercera dimensión no introduce error.
- **Cilindro hueco 3D**: perfil logarítmico, extensión del 2D — análogo térmico del Lamé 3D que el proyecto ya valida en mecánica.
- **Cross-check contra `Quad4Thermal`**: la misma pared plana resuelta en 2D y con una capa de `Hex8Thermal` con las caras $z$ adiabáticas debe dar campos idénticos nodo a nodo. Es el análogo del cross-check 3D↔2D `plane_strain` que en A.bis blindó los materiales a 10-14 decimales, y **cubre los errores en la tercera dimensión sin necesidad de mallador**.

La **esfera hueca** —único caso con flujo genuinamente tridimensional— queda **diferida** de esta etapa (decisión del usuario, 2026-08-25): exige mallar un octante de corona esférica, y la cobertura que aportaría ya la da el cross-check anterior. Ver §Diálogo.

---

## Contrato de implementación

```yaml
name: Hex8Thermal
kind: element
status: validated        # draft → implemented → validated

interface:
  field: temperature
  dof_names: [T]
  n_nodes: 8
  flux_dim: 3                   # ∇T y q en 3D
  n_integration_points: 8       # default Gauss 2×2×2

parameters:
  - { name: quadrature, type: str, required: false, default: "hex_2x2x2",
      desc: "Regla desde QuadratureRegistry. hex_1x1x1 disponible pero deja 4
             modos hourglass de temperatura — sin estabilización" }
  # Sin `thickness`: el volumen sale de la geometría (igual que los sólidos 3D).

material_contract:
  signature: "compute_flux(∇T) -> (q, k)"
  gradient_kind: "gradiente físico [∂T/∂x, ∂T/∂y, ∂T/∂z] — sin notación Voigt"
  accepted: [ThermalConduction]

conventions:
  sign_flux: "q̄ > 0 ⇔ flujo SALIENTE del dominio; normal exterior"
  adiabatic_default: "cara sin condición declarada ⇒ q̄ = 0 (aislada)"
  node_orientation: "orden VTK_HEXAHEDRON, volumen positivo — idéntico al Hex8"
  face_numbering: "FACE_NODES paritaria con Hex8: 0(−ζ) 1(+ζ) 2(−η) 3(+ξ) 4(+η) 5(−ξ)"
  capacity_default: "lumping='lumped' por default; razón física en Quad4Thermal §5"

validity:
  - "conducción lineal; k independiente de T"
  - "material con field='temperature' (familia ThermalMaterial)"
  - "det J > 0 en todos los puntos de Gauss"
  - "el modelo global debe imponer al menos un Dirichlet de temperatura"

out_of_scope:
  - "convección (Robin) y radiación en la frontera — Etapa 8 posterior"
  - "acoplamiento termomecánico"
  - "conducción con advección; cambio de fase"
  - "estabilización del hourglass con cuadratura reducida"

acceptance:
  verification:
    - name: patch_test_termico_lineal_3d
      setup: "malla de Hex8Thermal distorsionados con nodo(s) interior(es);
             imponer T = a + b·x + c·y + d·z en los nodos del contorno"
      expect: "los nodos interiores adoptan el campo lineal exacto y ∇T =
              (b, c, d) constante en TODOS los puntos de Gauss"
      tol_rel: 1.0e-10
    - name: modo_nulo_temperatura_uniforme
      setup: "K_e de un elemento aislado; campo T = ones(8)"
      expect: "K_e·T = 0 exacto; rango(K_e) = 7; exactamente un autovalor nulo"
      tol_abs: 1.0e-12
    - name: hourglass_con_integracion_reducida
      setup: "K_e con quadrature='hex_1x1x1' sobre elemento aislado"
      expect: "rango(K_e) = 3; 5 autovalores nulos (1 físico + 4 hourglass).
               Documenta la limitación, no la corrige"
      tol_abs: 1.0e-12
    - name: simetria_y_semidefinicion
      setup: "K_e y C_e con k isótropo y ortótropo, elemento regular y distorsionado"
      expect: "K_e == K_e.T y C_e == C_e.T exactas; autovalores de K_e ≥ 0 con un
               único cero; autovalores de C_e > 0 estrictos"
      tol_abs: 1.0e-12
  specific:
    - name: conduccion_1d_pared_plana_3d
      setup: "bloque de Hex8Thermal entre x=0 y x=L; T(0)=T_1, T(L)=T_2;
             las otras cuatro caras adiabáticas; estacionario"
      expect: "perfil lineal exacto en todos los nodos; q uniforme y sin
              componentes transversales (q_y = q_z = 0)"
      tol_rel: 1.0e-12
    - name: cilindro_hueco_logaritmico_3d
      setup: "sector de corona circular extruido; T(a)=T_a, T(b)=T_b;
             caras planas adiabáticas; estacionario"
      expect: "perfil logarítmico T(r) = T_a + (T_b−T_a)·ln(r/a)/ln(b/a);
              convergencia con refinamiento h"
      tol_rel: 5.0e-3
    # Esfera hueca (perfil T(r) = T_a + (T_b−T_a)·(1/a − 1/r)/(1/a − 1/b)):
    # DIFERIDA de esta etapa por decisión del usuario (2026-08-25). Es el único
    # caso con flujo genuinamente tridimensional, pero exige mallar un octante
    # de corona esférica —trabajo geométrico que ya difirió la campaña de
    # validación 3D de A.bis por el mismo motivo— y la cobertura que aportaría
    # sobre la tercera dimensión ya la da `consistencia_con_quad4thermal`, que
    # no necesita mallador. Retomar si entra un mallador 3D nativo (gmsh API)
    # o si aparece un caso de uso con flujo radial esférico.
    - name: fuente_volumetrica_uniforme
      setup: "Hex8Thermal regular y distorsionado con Q uniforme"
      expect: "Σf = Q·V_e exacto, invariante ante distorsión"
      tol_rel: 1.0e-12
    - name: flujo_neumann_en_una_cara
      setup: "Hex8Thermal con q̄ constante sobre cada una de las 6 caras"
      expect: "Σf = -q̄·A_cara con reparto A/4 a cada nodo de la cara y 0 en el
              resto; verificado en las 6 caras"
      tol_rel: 1.0e-12
    - name: balance_energetico_estacionario
      setup: "dominio con fuente Q y Dirichlet en parte de la frontera"
      expect: "Σ Q·V + Σ reacciones = 0"
      tol_rel: 1.0e-10
    - name: capacidad_lumped_vs_consistente
      setup: "C_e en ambos modos para elemento regular"
      expect: "ambas suman ρ·c·V_e; la lumped es diagonal y estrictamente positiva"
      tol_rel: 1.0e-12
    - name: anisotropia_desvia_el_flujo_3d
      setup: "k = diag(k1, k2, k3) distintos; campo lineal oblicuo a los ejes"
      expect: "q coincide con -k·∇T calculado a mano; no antiparalelo a ∇T"
      tol_rel: 1.0e-12
    - name: consistencia_con_quad4thermal
      setup: "problema de pared plana resuelto con Quad4Thermal (2D) y con una
             capa de Hex8Thermal con las caras z adiabáticas"
      expect: "campos de temperatura idénticos nodo a nodo; es el cross-check
              dimensional análogo al que el proyecto usa entre materiales 3D y
              2D plane_strain"
      tol_rel: 1.0e-12
    - name: jacobiano_degenerado_abortado
      setup: "elemento con nodos coplanares o numeración invertida"
      expect: "ValueError con mensaje claro"
    - name: rechazo_material_mecanico
      setup: "construir Hex8Thermal con un material de la familia mecánica"
      expect: "TypeError/ValueError identificando la incompatibilidad de familia"

references:
  - "Lewis R.W., Nithiarasu P., Seetharamu K.N. (2004). Fundamentals of the FEM for Heat and Fluid Flow. Wiley. §3, §6."
  - "Zienkiewicz O.C., Taylor R.L. (2000). The Finite Element Method, Vol. 1. §7."
  - "Carslaw H.S., Jaeger J.C. (1959). Conduction of Heat in Solids. Oxford. §7.2 (cilindro hueco). §9.3 cubre la esfera hueca, diferida de esta etapa."
  - "Reddy J.N., Gartling D.K. (2010). The Finite Element Method in Heat Transfer and Fluid Dynamics. CRC Press. §4."
```

---

## Implementación

- Archivo: [`solidum/elements/thermal/hex8_thermal.py`](../../solidum/elements/thermal/hex8_thermal.py)
- Clase: `Hex8Thermal`, registrada en `ElementRegistry`
- Base compartida: [`solidum/elements/thermal/_shared.py`](../../solidum/elements/thermal/_shared.py) · `_ThermalSolid`
- Tests: [`tests/test_hex8_thermal.py`](../../tests/test_hex8_thermal.py) — 26 verdes

**Notas de implementación:**

- **La spec de extensión se confirmó en el código**: la clase concreta ocupa poco más de cien líneas y sólo aporta funciones de forma, gradiente, `FACE_NODES` y el flujo por cara. Todo lo demás —validación del material, conductividad, capacidad en ambos modos, fuente volumétrica, post-proceso— se hereda intacto de `_ThermalSolid`. Es la evidencia de que centralizar con el primer elemento fue correcto: el segundo caso real no necesitó tocar la base.

- **Kernel `_compute_gradient_kinematics_hex8` en `solid_3d/_shared.py`**, junto a su gemelo mecánico, por la misma razón que en 2D: comparte el jacobiano 3×3 y su inversión por cofactores, que es la parte cara.

- **`thickness` heredado como 1.0 y no expuesto** en el constructor. En 3D el volumen es geométrico; permitir el parámetro habría abierto la puerta a escalar el volumen dos veces. Blindado con test: escalar la geometría ×2 multiplica la fuente por 8, no por otro factor.

**Predicciones de la spec confirmadas por test:**

| Predicción | Medido |
|---|---|
| `rango(K_e) = 7`, un único modo nulo | ✔ exacto |
| Reducida `hex_1x1x1` ⇒ rango 3, 4 hourglass | ✔ 5 autovalores nulos (1 físico + 4 espurios) |
| Flujo por cara reparte `-q̄·A/4` | ✔ en las 6 caras |
| Patch test lineal 3D exacto | ✔ también en elemento distorsionado |

**Cross-check 2D↔3D — el test central de esta spec.** La misma pared plana resuelta con `Quad4Thermal` y con una capa de `Hex8Thermal` con las caras $z$ adiabáticas da campos **idénticos nodo a nodo con error $9.9\times10^{-14}$**, y ambos coinciden con el perfil lineal analítico. Es el análogo térmico del cross-check 3D↔2D `plane_strain` que en A.bis blindó los materiales, y es la razón por la que el benchmark de la esfera hueca pudo diferirse sin dejar hueco de cobertura: detecta errores en la tercera dimensión sin necesitar mallador.

**Pendiente de esta spec — cilindro hueco 3D y balance energético.** Igual que en el hermano 2D, se agrupan con la campaña de validación del cierre de etapa por requerir mallado de corona circular. El estado `validated` refleja que la formulación está verificada contra solución analítica y contra el elemento 2D.

---

## Diálogo

*Puntos abiertos para el usuario. Sólo física y alcance.*

**2026-08-25 — IA → usuario. Un punto:**

1. ~~**Benchmark de la esfera hueca.**~~ ✅ **Resuelto 2026-08-25 — DIFERIDA de esta etapa** (decisión del usuario).

   La esfera hueca es el único caso propuesto con flujo genuinamente tridimensional: su perfil $\propto 1/r$ no lo producen ni la pared plana ni el cilindro, donde la temperatura varía en una sola dirección. Pero exige **mallar un octante de corona esférica**, trabajo geométrico que ya difirió la campaña de validación 3D de A.bis por el mismo motivo.

   Se difiere porque **la cobertura que aportaría está cubierta por otra vía**: `consistencia_con_quad4thermal` detecta errores en la tercera dimensión comparando el mismo problema resuelto en 2D y en 3D, y no necesita mallador. La esfera daría confirmación adicional, no una comprobación ausente.

   Nota técnica: cuando se retome, el `Hex8Thermal` **no** sufrirá el problema que bloquea al `Tet10` sobre superficie curva (deuda técnica #8) — al ser trilineal no tiene nodos intermedios que puedan quedar sobre aristas rectas. El bloqueante es sólo de mallado.

   **Criterio de retoma**: si entra un mallador 3D nativo (gmsh API) o aparece un caso de uso con flujo radial esférico.
