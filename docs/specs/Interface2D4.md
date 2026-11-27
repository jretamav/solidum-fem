# INTERFACE 2D4 — elemento de interfaz de espesor nulo (familia de interfaz, ADR 0043)

> Orden de trabajo pre-redactada por la IA el 2026-11-26 dentro del bloque de **contacto acotado** aprobado por el usuario el mismo día. Esta spec contiene la formulación de toda la familia (`Interface2D4`, `Interface3D6`, `Interface3D8`); las specs de los elementos 3D sólo documentan lo que cambia. La arquitectura (paquete, contrato, generación desde la malla, salida) está en el ADR 0043. El usuario revisa la física y aprueba antes de que la IA escriba código; las decisiones abiertas están en el **Diálogo**.

---

## Especificación física

### 0. Descripción general

Elemento de interfaz de **espesor nulo** (Goodman, Taylor y Brekke 1968; de Borst et al. 2012 §13.1, Zienkiewicz-Taylor-Fox §4.11.2, confirmados por el bibliotecario) que une dos caras A y B por pares de nodos conocidos y transmite entre ellas la tracción de un **material de interfaz de superficie** (ADR 0033): contacto (`ContactPenalty`), grieta cohesiva con trayectoria conocida (`CohesiveDamageIsotropic`, `CohesiveDamageMixedMode`) o cualquier otra ley de superficie registrada.

- `Interface2D4`: 2D (deformación plana o esfuerzo plano: la interfaz no distingue), caras lineales de 2 nodos, 4 nodos, DOF `(ux, uy)`.
- `Interface3D6` / `Interface3D8`: 3D, caras triangulares de 3 nodos o cuadriláteras de 4 nodos, DOF `(ux, uy, uz)`.

Desplazamientos pequeños: la geometría, la normal, los ejes tangentes y los pares son los de la configuración de referencia. El elemento no tiene masa ni rigidez propia: todo lo pone la ley.

### 1. Nodos y caras

Orden de los elementos cohesivos de Abaqus (ADR 0043 §1):

- `Interface2D4`: cara A = nodos (1, 2); cara B = nodos (3, 4), con **4 frente a 1 y 3 frente a 2** (los cuatro nodos recorren el contorno en sentido antihorario cuando B está del lado de `n`).
- `Interface3D6`: cara A = (1, 2, 3); cara B = (4, 5, 6), con 4 frente a 1, 5 frente a 2 y 6 frente a 3.
- `Interface3D8`: cara A = (1, 2, 3, 4); cara B = (5, 6, 7, 8), con 5 frente a 1, etc.

Con nodos coincidentes (interfaz de espesor nulo) las dos caras tienen las mismas coordenadas; en un emparejamiento con holgura (`pair`, ADR 0043 §6) cada nodo de B está frente a su par de A a lo largo de la normal de A.

### 2. Ejes locales

Sobre la **cara A** en la configuración de referencia, en cada punto de integración `ξ`:

- 2D: `s = (x_2 − x_1)/‖x_2 − x_1‖` y `n = e_z × s` (`s` girado +90°). Ejes `(n, s)`.
- 3D: `a_ξ = ∂x_A/∂ξ`, `a_η = ∂x_A/∂η`; `n = a_ξ × a_η / ‖·‖`; `s1 = a_ξ/‖a_ξ‖`; `s2 = n × s1`. Ejes `(n, s1, s2)`, dextrógiros. En la cara triangular (plana) son constantes; en la cuadrilátera, que puede estar alabeada, se evalúan en cada punto de Gauss y, con la regla nodal, los del nodo `i` en `(ξ_i/3, η_i/3)`, donde `a_ξ × a_η` es el vector de área tributaria consistente `∫N_i (a_ξ × a_η)` (spec `Interface3D8`, decisión del usuario del 2026-11-26).

**`n` apunta de A hacia B.** Con nodos coincidentes no hay geometría que lo diga: lo fija el orden de los nodos, y el generador (ADR 0043 §6) ordena la cara A de modo que `n` sea la normal exterior del sólido del lado A. En 2D eso significa recorrer la cara A **en sentido contrario** al de esa arista en el contorno antihorario del sólido A, que es el sentido de la misma arista en el sólido B (*calculado* por el teórico).

**Qué cara es A.** Con pares a distancia, A es la cara **plana o la más regular** (en Hertz, el plano: `n = e_y`). Si A fuera la cara curva facetada, el par no quedaría sobre su normal (desfase tangencial `g_0·tan θ`) y el generador lo rechazaría.

### 3. Desplazamiento relativo y holgura

Con `N_i(ξ)` las funciones de forma de la cara (lineales en 2D; de área en el triángulo; bilineales en el cuadrilátero), iguales en las dos caras:

$$\mathbf u_A(\xi) = \sum_i N_i(\xi)\,\mathbf u_{A,i}, \qquad \mathbf u_B(\xi) = \sum_i N_i(\xi)\,\mathbf u_{B,i}, \qquad \Delta\mathbf u = \mathbf u_B - \mathbf u_A$$

**Holgura inicial** en el punto de integración:

$$g_0(\xi) = \mathbf n(\xi)\cdot\bigl(\mathbf x_B(\xi) - \mathbf x_A(\xi)\bigr)$$

con `x_A`, `x_B` interpoladas con las mismas `N_i`. Es 0 con nodos coincidentes. El desfase tangencial inicial `(x_B − x_A)·s` debe ser nulo salvo redondeo (el generador lo comprueba); no cuenta como deslizamiento.

Vector que recibe la ley, en ejes locales (ADR 0043 §4):

$$\mathbf g = \begin{bmatrix} g_0 + \mathbf n\cdot\Delta\mathbf u \\ \mathbf s\cdot\Delta\mathbf u\end{bmatrix} \ \ \text{(2D)}, \qquad \mathbf g = \begin{bmatrix} g_0 + \mathbf n\cdot\Delta\mathbf u \\ \mathbf s_1\cdot\Delta\mathbf u \\ \mathbf s_2\cdot\Delta\mathbf u\end{bmatrix} \ \ \text{(3D)}$$

Una ley que no declara `ACCEPTS_INITIAL_GAP` (cohesivas, adherencia) sólo se admite si `g_0 = 0` en todos los puntos (error al construir).

### 4. Equilibrio

Trabajo virtual de la interfaz, sumado al de los sólidos:

$$\delta W_\Gamma = \int_{\Gamma} \delta\mathbf g\cdot\mathbf t(\mathbf g)\; d\Gamma$$

con `Γ` medida sobre la cara A en la referencia (en 2D, `dΓ = e·dl`, con `e` el espesor fuera del plano de los sólidos adyacentes). Con pares a distancia (cilindro sobre plano), la medida es la de la cara A.

---

## Formulación numérica (FEM)

### 5. Matriz de salto

$$\mathbf g = g_0\,\mathbf e_1 + \mathbf B(\xi)\,\mathbf u_e, \qquad \mathbf B(\xi) = \mathbf R(\xi)\,\bigl[\,-\mathbf N(\xi)\otimes\mathbf I_d \ \big|\ \mathbf N(\xi)\otimes\mathbf I_d\,\bigr]$$

con `R` la matriz de filas `n, s` (2D) o `n, s1, s2` (3D), `u_e` ordenado según los nodos de §1 (cara A, cara B), y la cara B con sus `N_i` asignadas a los nodos frente a cada nodo de A.

### 6. Fuerzas internas y rigidez

$$\mathbf f_e = \sum_g w_g\,J_g\,e\;\mathbf B_g^{\top}\,\mathbf t_g, \qquad \mathbf K_e = \sum_g w_g\,J_g\,e\;\mathbf B_g^{\top}\,\mathbf T_g\,\mathbf B_g$$

con `J_g` el jacobiano de la cara A (longitud/2 en 2D; ‖a_ξ × a_η‖ en 3D), `e` el espesor en 2D (1 en 3D) y `T_g = ∂t/∂g` la tangente de la ley. `K_e` es simétrica si la ley lo es.

**Espesor en 2D**: el de los sólidos adyacentes; el generador lo toma de ellos y rechaza una cara cuyos dos lados tengan espesores distintos.

### 7. Integración

- **Nodal** (Newton–Cotes / Lobatto, puntos en los nodos): 2D, `ξ = ±1`, pesos 1, 1; triángulo, los tres vértices con peso 1/6; cuadrilátero, las cuatro esquinas con el área tributaria consistente `∫N_i J` (abajo, decisión Q3). Cada punto ve un solo par de nodos: la matriz no acopla pares distintos, es el **contacto nodo a nodo**, y no hay oscilación espuria de la tracción con rigideces grandes (Schellekens y de Borst 1993, pp. 57–58, *leído* por el bibliotecario; `docs/referencias/Schellekens_deBorst_1993_integracion_interfaces.pdf`).
- **Gauss**: 2D, 2 puntos; triángulo, 3 puntos interiores; cuadrilátero, 2×2. Integra exactamente la rigidez elástica, pero con penalización grande la tracción oscila entre puntos.

El mismo número de puntos en las dos reglas (2, 3 y 4).

**Pesos nodales del cuadrilátero** (hallazgo 1 del teórico, *calculado*): con peso 1 en cada esquina y el `J` de la esquina, el área que recibe cada nodo no es la consistente `∫N_i J dξdη` si la cara no es un paralelogramo. En una cara trapezoidal plana el cociente va de 0.675 a 1.245 (el área total sí coincide); en una malla de 3×3 caras con los nodos perturbados, de 0.78 a 1.13. Bajo esfuerzo uniforme la presión de contacto de cada par se desviaría lo mismo. En la línea y en el triángulo la regla nodal ya da el área consistente (`L/2`, `A/3`). Opciones: (a) pesos nodales iguales al área tributaria consistente `w_i = ∫N_i J` (conserva el desacople por pares y pasa la parcela en cualquier cara plana); (b) Newton–Cotes, exacta sólo en paralelogramos. **Decidido (a)** por el usuario (Diálogo Q3).

Precedente: en la adherencia barra–sólido el usuario eligió integración **nodal por omisión con Gauss como opción** (2026-11-10, spec `EmbeddedReinforcement`). **Decidido igual aquí** (Diálogo Q1).

### 7b. Escalera de validación (la corre `validador-solidum`)

1. Punto material (pruebas de la ley).
2. **Parcela de contacto por pares**: dos bloques apilados en deformación plana, de alturas `H₁`, `H₂` y materiales distintos, presión `p` arriba. Apoyos que dejan libre el efecto de Poisson: `u_y = 0` en toda la base y `u_x = 0` en un solo nodo de cada bloque (sin fricción el bloque superior tiene un modo rígido en `x`). Exacto: `t_n = −p` en todos los pares, interpenetración `p/K_n`, desplazamiento de la cara superior `−p·[(1−ν₁²)H₁/E₁ + (1−ν₂²)H₂/E₂] − p/K_n − g_0` (corrección del signo de `g_0` por el teórico).
3. **Despegue**: un par en serie con una barra; `F = 0` hasta cerrar la holgura y `F = K_eff·(u − g_0)` después, `1/K_eff = L/(EA) + 1/(K_n·A)`.
4. **Bloque con fricción**: presión prescrita y rampa tangencial; umbral `μ·p` exacto, descarga elástica, deslizamiento residual; en 3D, trayectoria no colineal.
5. **Hertz, cilindro sobre plano rígido en deformación plana** (Wriggers 2006, §14.6.1): `E = 7000`, `ν = 0.3`, `R = 1`, `P = 100` por unidad de espesor; `a = √(4PR/(πE*)) = 0.128655`, `p₀ = 2P/(πa) = 494.83`, `p(x) = p₀√(1−(x/a)²)`. Se comparan `a`, `p₀`, `p(x)` y la reacción, no el acercamiento (indefinido en deformación plana). Pares con paso `h ≲ a/10` en la zona de contacto; la cara A es el plano; holguras desde la geometría del círculo (la parábola `x²/(2R)` difiere 0.4 % en `x = a`); con media simetría, `P/2`. Con `K_n = 10E/h` y `h = a/10`, la interpenetración `p₀/K_n` es el 1.1 % de `g_0(a)` (*calculado* por el teórico). Criterio: convergencia al refinar a pocos % de `a` y `p₀`.
6. Viga SEN con la placa apoyada, como observación.

**Resultado de la escalera.** El 2026-11-27, con `validador-solidum` (*medido*). Casos 2–5 en 2D con `Interface2D4` y `ContactPenalty` (`informes/reportes/reporte_contacto_2026-11-27.pdf`): parcela exacta (error ≤ 1.4e-16, `pair` y `split`, nodal y Gauss); despegue con error 1e-10 (la tolerancia de Newton); fricción con el umbral `μp` exacto, descarga elástica y deslizamiento residual (error 2e-15); Hertz con `p₀` a −1.37/−0.77/−0.52 % y `a` a +0.06/+0.40/+0.57 % en tres mallas. Con etapas (`informes/reportes/reporte_contacto_etapas_2026-11-27.pdf`): despegue con carga y descarga, fricción 2D de ida y vuelta (energía disipada exacta en todos los pasos) y fricción 3D con trayectoria en L con `Interface3D8` (retorno radial discreto a 4e-10·`μpA`; frente a la solución continua 1.2/0.30/0.08 % con 20/80/320 pasos, orden 1). Caso 6 (`informes/reportes/reporte_sen_placa_contacto_2026-11-27.pdf`): sin fricción funciona y coincide con la carga como presión a 1.1 % del pico; el contacto es coherente en todos los pasos (equilibrio 4.8e-10, `t_n ≤ 0`, `|t_s|/μ|t_n|` = 1.0000 al deslizar, interpenetración ≤ 1 µm = `p/K_n`). Con μ = 0.5 y la placa libre en horizontal, el análisis se detiene en la cola (tangente singular). *Hipótesis*: al cruzar la grieta, las dos partes de la viga sólo quedan unidas en horizontal por la fricción de la placa, y si sus dos bordes deslizan a la vez la posición horizontal de la placa queda indeterminada. La representación de la placa con fricción queda abierta (Diálogo, 2026-11-27). `Interface2D6` e `Interface3D6` no se corrieron en la escalera.

Fuentes de pago para ampliar (que el usuario puede bajar en la UNAM): Taylor y Papadopoulos 1991 (parcela de contacto) y NAFEMS R0081 (Prinja y Feng 2001; parcela, Hertz y cuña con fricción).

### 8. Energía almacenada

$$E_e = \sum_g w_g\,J_g\,e\;\psi\bigl(\mathbf g_g, \text{estado}_g\bigr)$$

con `ψ` la de la ley (`stored_energy`, ADR 0043 §4–§5). En la historia, `E = Σ_e E_e` y `D = W − E`.

### 9. Advertencias y límites

- **Giro rígido con holgura**: con `g_0 ≠ 0`, un giro rígido pequeño `ω` de todo el conjunto desplaza la cara B respecto de la A en `g_s = −ω·g_0` tangencialmente (*calculado*: lineal en `ω`, no de segundo orden). Su dual: con fricción, las fuerzas tangenciales `±t_s` de las dos caras forman un par espurio `t_s·g_0` por unidad de área en la geometría de referencia. Es inocuo con la interfaz abierta (el deslizamiento se olvida) o sin fricción; con nodos coincidentes no aparece. Propio de los desplazamientos pequeños.
- **Cara cuadrilátera alabeada**: la normal varía dentro de la cara. Con la regla nodal, cada nodo toma la normal y la medida de su vector de área tributaria consistente, y la cara transmite una presión uniforme con el error de la penalización (spec `Interface3D8`; con la normal de la esquina el error era de 3 % en los esfuerzos con 2.4° de alabeo, *medido*). Con torsión y curvatura a la vez, en un nodo se juntan caras de normales distintas y queda un residuo de segundo orden en la tracción nodal, que converge al refinar (1.49e-3, 3.75e-4, 9.4e-5 con 3, 6 y 12 divisiones; Gauss no lo tiene; *medido* por el teórico). **Con holgura, la cara A debe ser plana** (§2): sobre una cara alabeada el par deja de estar sobre la normal del nodo y Solidum lo rechaza en cuanto `g_0·sen θ_i > 10⁻⁶·h` (con 2.4° de alabeo, toda holgura mayor que 1.4e-4·h).
- **Interfaz abierta**: sin rigidez; un cuerpo sostenido sólo por interfaces abiertas es un mecanismo.
- **Análisis de rigidez fija** (lineal, modal, armónico, espectral, transitorios lineales) con una ley unilateral: error, «el contacto requiere un análisis no lineal» (decisión Q4; ADR 0043 «Consecuencias»). Con cohesivas, admitidos.

---

## Contrato de implementación

```yaml
name: Interface2D4
kind: element
status: validated

interface:
  field: interface              # campo de interfaz: salto en ejes locales, JUMP_DIM (ADR 0043 §3)
  dof_names: [ux, uy]
  n_nodes: 4
  jump_dim: 2                   # (n, s)
  n_integration_points: 2       # nodal o Gauss, los dos con 2 puntos

admission:                      # criterio de admisión al programa principal (Reglas.md §1); evidencia en informes/bibliotecario/busqueda_contacto_2026-11-26.md
  textbooks:
    - de_borst_crisfield: "§13.1 Interface Elements (pp. 428–436): cinemática, matriz B_d, rigidez, integración (p. 432)"
    - zienkiewicz_taylor: "vol. 2, 7.ª ed., §4.11.2 elementos de junta (pp. 127–128)"
  general_codes:
    - abaqus: "COH2D4 (cohesivo de 4 nodos, ley tracción-separación)"
    - ansys: "INTER202 (interfaz lineal de 4 nodos)"
    - diana: "cap. 33, Structural Interface Elements"

parameters:
  - { name: material,    type: int,   required: true,  desc: "Id en interface_materials; ley de superficie con 2 componentes" }
  - { name: integration, type: str,   required: false, desc: "'nodal' (por omisión, Diálogo Q1) | 'gauss'" }
  - { name: thickness,   type: float, required: false, desc: "Espesor fuera del plano; lo pone el generador desde los sólidos adyacentes" }

material_contract:
  signature: "compute_traction(g, estado) -> (t, T, estado')"
  strain_kind: "salto local (g_n con holgura inicial, g_s)"

conventions:
  sign: "g_n > 0 separación (n de A hacia B); t en ejes locales con la convención de la ley"
  node_orientation: "A = (1, 2), B = (3, 4) con 4 frente a 1 y 3 frente a 2 (antihorario si B está del lado de n)"

validity:
  - "desplazamientos y deslizamientos pequeños"
  - "pares alineados según la normal de A (desfase tangencial nulo salvo redondeo)"
  - "g_0 = 0 con leyes sin ACCEPTS_INITIAL_GAP"

out_of_scope:
  - "Caras cuadráticas 3D (Interface3D15/20): después (la 2D, Interface2D6, entró en el bloque de contacto)"
  - "Desplazamientos grandes, deslizamiento grande, búsqueda de pares"
  - "Interfaces térmicas (conductancia de contacto)"
  - "Inserción de interfaces durante el análisis"

acceptance:
  verification:
    - name: parcela_de_contacto
      setup: "dos bloques con mallas coincidentes en la interfaz, presión uniforme p sobre el superior, ContactPenalty sin fricción; nodal y Gauss; espesor ≠ 1"
      expect: "presión de contacto p en todos los puntos; esfuerzo σ_yy = −p uniforme en los dos bloques; interpenetración p/K_n"
      tol_rel: 1.0e-10
    - name: salto_uniforme
      setup: "una interfaz aislada; cara A fija, cara B con desplazamiento uniforme (traslación); ley lineal de prueba (t = K·g)"
      expect: "t uniforme; fuerzas nodales = t·J·w·e en cada nodo (nodal) y suma = t·área en las dos reglas"
      tol_rel: 1.0e-12
  specific:
    - name: cuerpo_rigido
      setup: "traslación y giro rígido pequeño de las dos caras juntas con nodos coincidentes"
      expect: "g = 0 y f_e = 0"
      tol_rel: 1.0e-12
    - name: holgura_inicial
      setup: "pares separados g_0 a lo largo de n; ContactPenalty; cerrar la holgura con desplazamiento impuesto"
      expect: "t = 0 mientras Δu_n > −g_0; contacto exactamente al cerrarse"
      tol_rel: 1.0e-12
    - name: tangente_del_elemento
      setup: "K_e frente a diferencias centradas de f_e en estados abiertos, cerrados, en adherencia y en deslizamiento"
      expect: "coinciden"
      tol_rel: 1.0e-6
    - name: orientacion_invertida
      setup: "nodos coincidentes; intercambiar las caras A y B (normal opuesta)"
      expect: "mismas fuerzas nodales (con holgura variable no vale: J_A ≠ J_B y n_A ≠ −n_B)"
      tol_rel: 1.0e-12
    - name: escalado_con_el_espesor
      setup: "espesor e ≠ 1"
      expect: "K_e, f_e y E_e proporcionales a e"
      tol_rel: 1.0e-12
    - name: holgura_variable_en_el_elemento
      setup: "cara B inclinada respecto de A; nodal y Gauss"
      expect: "g_0 en cada punto = n·(x_B − x_A) interpolado"
      tol_rel: 1.0e-12
    - name: simetria_de_la_rigidez
      setup: "ContactPenalty en deslizamiento con μ > 0, y con μ = 0"
      expect: "K_e no simétrica con fricción en deslizamiento; simétrica con μ = 0"
      tol_rel: 1.0e-14
    - name: giro_rigido_con_holgura
      setup: "giro rígido pequeño ω con g_0 ≠ 0"
      expect: "g_s = −ω·g_0 (límite documentado de §9)"
      tol_rel: 1.0e-10
    - name: energia_al_cerrar_la_holgura
      setup: "cerrar una holgura en 1 paso y en N pasos con ContactPenalty sin fricción; energía de la historia"
      expect: "E = Σ w·ψ exacta; D = W − E tiende a 0 al refinar, acotada por ½·|t_{n+1}|·g_n·área en el paso del cierre (trapecio del ADR 0038)"
      tol_rel: 1.0e-12
    - name: nodal_desacopla_pares
      setup: "integración nodal"
      expect: "K_e no acopla nodos de pares distintos"
      tol_rel: 0.0
    - name: cohesiva_con_holgura_rechazada
      setup: "CohesiveDamageIsotropic con g_0 ≠ 0"
      expect: "error al construir"

references:
  - "Goodman, R. E., Taylor, R. L., Brekke, T. L. (1968). A model for the mechanics of jointed rock. J. Soil Mech. Found. Div. 94(3), 637–659. DOI 10.1061/JSFEAQ.0001133"
  - "Schellekens, J. C. J., de Borst, R. (1993). On the numerical integration of interface elements. IJNME 36(1), 43–66. DOI 10.1002/nme.1620360104"
  - "de Borst, R., Crisfield, M. A., Remmers, J. J. C., Verhoosel, C. V. (2012). Non-Linear Finite Element Analysis of Solids and Structures, 2.ª ed., §13.1. Wiley. DOI 10.1002/9781118375938"
  - "Zienkiewicz, O. C., Taylor, R. L., Fox, D. D. (2014). The Finite Element Method for Solid and Structural Mechanics, 7.ª ed., §4.11.2 y cap. 8. DOI 10.1016/B978-1-85617-634-7.00008-9"
  - "Wriggers, P. (2006). Computational Contact Mechanics, 2.ª ed., §14.6.1 (Hertz en deformación plana: E = 7000, ν = 0.3, R = 1, P = 100, p_max = 494.83). Springer. DOI 10.1007/978-3-540-32609-0"
  - "ADR 0033 — Familia de materiales de interfaz."
  - "ADR 0043 — Elementos de interfaz de espesor nulo."
```

---

## Implementación

- Archivo: `solidum/elements/interface/interface_2d.py` (clase); lo común de la familia en `solidum/elements/interface/_shared.py` (`InterfaceElement` sobre `LocalJumpElement`, compartida con `BondSlipInterface`: ejes, holgura inicial, reglas `nodal`/`gauss`, áreas tributarias consistentes, `stored_energy`, `interface_points`). Generación en `solidum/mesh/interfaces.py`; salida en `solidum/utils/vtk_exporter.py` y `output_writer.py` (ADR 0043 §6–§7).
- Clase: `Interface2D4` (`FACE_CELL = "line"`, `A_NODES = (0, 1)`, `B_NODES = (3, 2)`, `VTK_CELL_TYPE = "quad"`).
- Tests:
  - `tests/test_interface_elements.py`: casos de `acceptance` de la familia y parcela de contacto 2D (dos materiales, `g_0 = 0` con presión y `g_0 > 0` con bajada impuesta, nodal y Gauss).
  - `tests/test_interface_generation.py` (`pair`/`split`, parcela desde el YAML con espesor 3), `test_interface_output.py`, `test_interface_energy.py`, `test_point_quantity.py`, barrido de contrato.
- Verificación: `informes/teorico/verificacion_teorica_2026-11-26_contacto.md` (medidas, `g_0`, `f_e`, `K_e`, `E_e` frente a una referencia independiente ≤ 1.2e-15; Newton con la tangente ensamblada en 2–4 iteraciones; *calculado*).
- Notas de traducción: orden de nodos de COH2D4 (Abaqus).

---

## Diálogo

*Preguntas, aclaraciones y hallazgos. Entradas fechadas.*

- **2026-11-26 — Decisiones abiertas para el usuario** (una a la vez, después de la revisión de `teorico-solidum`):
  - **Q1 — DECIDIDO por el usuario el 2026-11-26 («si»): integración nodal por omisión, Gauss como opción.** Pregunta: integración por omisión: ¿nodal con Gauss como opción, como en la adherencia? Datos (*leído*, informe del bibliotecario): de Borst et al. (p. 432) y Schellekens y de Borst (1993) recomiendan Newton–Cotes / Lobatto (en elementos lineales, la nodal) porque eliminan la oscilación de la tracción; ANSYS (INTER202) y DIANA (§33.4) usan Gauss por omisión, y DIANA advierte de esas oscilaciones con rigidez alta. Con contacto unilateral, la nodal da un estado abierto o cerrado por par de nodos (contacto nodo a nodo de Zienkiewicz §8.2). Recomendación de la IA: nodal por omisión, Gauss como opción. *Corrección posterior del bibliotecario* (2026-11-26, `informes/bibliotecario/busqueda_interfaz_cuadratica_2026-11-26.md`, *leído*): DIANA declara **Newton–Cotes por omisión en cada elemento de interfaz** (§33.7: L8IF 2 puntos, CL12I 3, Q24IF 2×2, CQ48I 3×3); el «Gauss por omisión» de §33.4 es una ambigüedad del manual. Refuerza la decisión; no la cambia.
  - **Q2 — DECIDIDO por el usuario el 2026-11-26 («1»): normal, ejes y medida de la cara A, la que se declara como `a` en el YAML, que debe ser la plana o la más regular; `g_0 = n·(x_B − x_A)` en cada punto.** Pregunta: definición de `g_0` y de la normal con pares a distancia: normal y medida de la **cara A**, `g_0 = n·(x_B − x_A)` en cada punto (§2–§3), con A la cara plana o la más regular. La alternativa sería la superficie media; el teórico: no aporta nada en desplazamientos pequeños.
  - **Q3 — DECIDIDO por el usuario el 2026-11-26 («si»): pesos nodales = área tributaria consistente `w_i = ∫N_i J dξdη`** (en línea y triángulo coincide con Newton–Cotes; cambia sólo el cuadrilátero). Pregunta: pesos nodales del cuadrilátero: área tributaria consistente o Newton–Cotes (§7).
  - **Q4 — DECIDIDO por el usuario el 2026-11-26 («tu recomendación»): en el análisis lineal y en el modal, un modelo con leyes de interfaz unilaterales es error («el contacto requiere un análisis no lineal»).** Razón: nunca un resultado falso en silencio; permitirlo después (estado congelado con aviso en el modal) no cambia modelos existentes. Pregunta y opinión del teórico: análisis lineal y modal con interfaces: el teórico recomienda error en el estático lineal con leyes unilaterales (con la rama cerrada en `g_0 = 0` sería un resorte bilateral que admite tracción sin avisar), o comprobar después `g_n ≤ 0` en todos los puntos; en el modal, estado congelado del paso base (cerrado: `K_n`, y `K_t` con fricción; abierto: 0) con aviso.
  - **Q5 — DECIDIDO por el usuario el 2026-11-26 («si»): tolerancia `10⁻⁶·h` (h, tamaño de la cara). Una holgura inicial `g_0` con `−10⁻⁶·h ≤ g_0 < 0` se lleva a 0; por debajo de `−10⁻⁶·h`, error «las caras se interpenetran al inicio».** Las interferencias intencionales (ajustes a presión) quedan fuera de este bloque. **Precisión DECIDIDA por el usuario el 2026-11-26 («ok 1», hallazgo 2 de la verificación de `teorico-solidum`): la tolerancia es simétrica, `|g_0| ≤ 10⁻⁶·h` se lleva a 0** (también las holguras positivas pequeñas): sin ello, dos caras coincidentes salvo redondeo se rechazarían con una ley cohesiva («las caras están separadas»); el contacto se cierra como mucho `10⁻⁶·h` antes (`10⁻⁵·E` de tracción con `K_n = 10E/h`, *calculado*). Es lo que ya hacía el código (`_shared.py`, `GAP_TOLERANCE`). Pregunta y datos: interpenetración inicial (`g_0 < 0` da tracción en `U = 0`: con `K_n = 10E/h`, `10⁻⁶·h` produce `10⁻⁵·E`, 0.3 MPa con E = 30 GPa, *calculado*): el teórico recomienda llevar a 0 lo que esté bajo una tolerancia geométrica relativa al tamaño de la cara y dar error por encima; el valor de la tolerancia es del usuario.
  - **Q6 — DECIDIDO por el usuario el 2026-11-26 («si»): nota en el manual de usuario** (sección de energías y de interfaces) que explique la disipación espuria del trapecio cuando un paso cruza un quiebre de la ley y que se reduce con pasos más pequeños. Pregunta: ¿el manual advierte que `D = W − E` da una disipación espuria cuando un paso cruza un quiebre de la ley (cierre de una holgura, paso de adherencia a deslizamiento)? Se reduce al refinar el paso (1 paso: `D/ψ = 0.50`; 10 pasos: 5e-3). *Precisión de la verificación* (`teorico-solidum`, 2026-11-26, *calculado*): el signo depende del quiebre: positiva al cerrar una holgura (el trapecio sobrestima `W`), negativa al pasar de adherencia a deslizamiento bajo presión constante (el trapecio subestima `W` en `½·μp·(g − g_y)`: la disipación de ese paso sale la mitad de la real). La nota del manual es neutra.
  - Decididos por la IA, con el acuerdo del teórico (rechazos de casos no soportados, sin recomendación física): una ley sin `ACCEPTS_INITIAL_GAP` con `g_0 ≠ 0` es error (si algún día se admite, lo correcto sería pasarle sólo `Δu`); `split` con elementos térmicos sobre los mismos nodos es error (conductancia perfecta y adiabático son los dos límites y ninguno es correcto en general).
- **2026-11-26 — Revisión de la spec por `teorico-solidum`** (`informes/teorico/revision_spec_contacto_2026-11-26.md`, «cambios antes de programar»): B, orientación, pesos de línea y triángulo, despegue y Hertz correctos (*calculado*). Aplicado: orientación de la cara A para el generador y cara A plana con pares a distancia (§2), signo de `g_0` y apoyos de la parcela (§7b.2), giro rígido de primer orden y par espurio con fricción (§9), datos de Hertz (§7b.5), pruebas nuevas; `orientacion_invertida` sólo con nodos coincidentes. Pendiente de comprobar al programar (*hipótesis* del teórico, de memoria): la celda `wedge` de VTK podría tener la base con la normal al revés que este orden (sólo visualización).
- **2026-11-26 — Revisión del plan por `arquitecto-solidum`** (`informes/arquitecto/revision_plan_contacto_2026-11-26.md`): orden de nodos de Abaqus, `JUMP_DIM` en lugar de `STRAIN_DIM`, celda VTK propia, generación con `pair`/`split`, energía almacenada por elemento: ADR 0043.
- **2026-11-27 — Validación (`validador-solidum`) y cierre del bloque de contacto** (decisión del usuario, «SI LA 1»): escalera §7b hecha (resultados en §7b); `status: validated` para `Interface2D4` y `ContactPenalty`. Viga SEN con la placa apoyada: el resultado de referencia es la placa sin fricción. **Queda abierto (decide el usuario cuando se retome)**: cómo representar la placa con fricción. Con la placa libre en horizontal, el análisis se detiene en la cola; con la rótula fija en horizontal, la rótula toma ~5 kN horizontales (11–12 % de 10P/11) y el pico sube hasta 17 %. Una rigidez horizontal del sistema de carga (resorte en la rótula) no está en los datos del ensayo. Sin explicar (de la IA): con la rótula fija y la malla de 2.5 mm, el Newton no converge en el paso 386; con la placa libre y la de 10 mm, la tangente es singular aunque dos puntos estén adheridos.
