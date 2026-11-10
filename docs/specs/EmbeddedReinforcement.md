# EmbeddedReinforcement — refuerzo embebido en sólidos con adherencia perfecta

> Orden de trabajo. El usuario escribe la **especificación física**, la **formulación numérica** y el **contrato**. La IA rellena **implementación** y responde en **diálogo** durante el trabajo.

---

## Especificación física

### 0. Descripción general
Barras de refuerzo rectas o poligonales dentro de un sólido 2D (esfuerzo o deformación plana) o 3D, con adherencia perfecta. La barra es una cadena de elementos 1D ordinarios de Solidum (`Truss2D`, `Truss3D`, sus variantes corrotacionales, los cables y las barras de tres nodos `Truss2D3` y `Truss3D3`), con cualquier material 1D registrado. Su malla es independiente de la del sólido: la define una polilínea, no los nodos del sólido. La adherencia perfecta se impone con restricciones afines (ADR 0004): cada nodo de la barra sigue el desplazamiento interpolado del elemento de sólido que lo contiene.

Es el concepto de refuerzo embebido de de Borst et al. (2012, §6.8.2): «the reinforcement undergoes the same displacements as the concrete elements in which it is embedded», con dos contribuciones al trabajo interno, la del sólido y la de la barra.

### 1. Cinemática de desplazamientos
Sea $\mathbf X_b$ la posición de referencia de un nodo de la barra, contenido en el elemento anfitrión $e$, con coordenadas naturales $\boldsymbol\xi_b$ (inversa del mapeo isoparamétrico, $\mathbf X_b = \sum_i N_i^e(\boldsymbol\xi_b)\,\mathbf X_i$). La adherencia perfecta es
$$\mathbf u(\mathbf X_b) = \sum_{i=1}^{n_e} N_i^e(\boldsymbol\xi_b)\,\mathbf u_i,$$
una restricción afín por componente de traslación ($u_x$, $u_y$ y, en 3D, $u_z$). Como $\boldsymbol\xi_b$ es la coordenada de un punto material, la restricción es exacta en el campo discreto del anfitrión también con desplazamientos grandes.

**Corte en los bordes del anfitrión.** La polilínea se corta en los bordes de los elementos que atraviesa (aristas en 2D, caras en 3D), de modo que cada tramo queda dentro de un solo anfitrión. Si la interpolación del anfitrión es lineal a lo largo de una recta, la deformación de la barra coincide con la del sólido proyectada sobre su eje, $\varepsilon_b = \mathbf t^\top \boldsymbol\varepsilon\,\mathbf t$. Es el caso de Tri3 y Tet4, y de Quad4 y Hex8 no distorsionados con la barra paralela a una arista. En los demás casos el tramo toma la cuerda del campo del anfitrión, es decir, su deformación media, y la diferencia con la proyección punto a punto disminuye al refinar la malla. Datos del 2026-11-10 (voladizo en esfuerzo plano, cuantía 2 %): la diferencia es ≤ 0.18 % con la malla más gruesa, frente a un error de malla de 14 %, y ≤ 0.02 % con 32×8. El parámetro `subdivisions` divide cada tramo en partes iguales.

Un nodo de la barra en el borde común de dos anfitriones conformes tiene la misma interpolación desde los dos; se ata a uno de ellos de forma determinista.

### 2. Cinemática de deformaciones
La del elemento 1D elegido: $\varepsilon = (\ell - L)/L$ en `Truss2D` y `Truss3D`, y la medida corrotacional en sus variantes.

### 3. Ecuación constitutiva
Cualquier material 1D registrado: `Elastic1D`, `Elastoplastic1D`, `ManderConcrete1D`, `CableMaterial1D`… El área de la barra $A_s$ es un parámetro del refuerzo.

**Área de concreto desplazada: no se descuenta** (decisión del usuario, 2026-11-10). El sólido ocupa todo su volumen y la barra se superpone, como en de Borst et al. (2012, ec. 6.136, que integra el concreto sobre todo $V$), Abaqus (`*EMBEDDED ELEMENT`) y DIANA. El sobrante relativo de rigidez axial elástica es
$$\frac{\rho}{(1-\rho) + n\rho},\qquad n = E_s/E_c ,$$
y el de resistencia a compresión es $\rho f'_c/[(1-\rho)f'_c + \rho f_y]$. Con $E_c$ = 30 GPa, $f'_c$ = 30 MPa y $f_y$ = 420 MPa:

| $\rho$ | rigidez elástica | resistencia a compresión |
|---|---|---|
| 1 % | +0.9 % | +0.9 % |
| 2 % | +1.8 % | +1.6 % |
| 4 % | +3.3 % | +2.6 % |

En tracción, con el concreto agrietado, el sobrante es cero. La masa del sólido tampoco se descuenta (sobrante $\rho_c A_s L_b$).

### 4. Equilibrio — forma fuerte
Sólido: $\nabla\cdot\boldsymbol\sigma + \mathbf b = \mathbf 0$ en $\Omega$. Barra: $\dfrac{d}{ds}(A_s\sigma_s) + q_s = 0$ a lo largo de $\Gamma_b$. Con adherencia perfecta, la fuerza de adherencia $q_s$ que la barra recibe del sólido la recibe el sólido con signo opuesto, y no aparece en la forma débil.

### 5. Equilibrio — forma débil
$$\int_\Omega \delta\boldsymbol\varepsilon:\boldsymbol\sigma\,d\Omega + \int_{\Gamma_b} \delta\varepsilon_b\,A_s\sigma_s\,ds = \delta W_{\text{ext}}, \qquad \delta\mathbf u_b = \delta\mathbf u(\mathbf X_b),$$
con los desplazamientos virtuales de la barra restringidos por la misma relación de §1 (de Borst et al. 2012, ec. 6.136).

---

## Formulación numérica (FEM)

### 6. Discretización
Los tramos son elementos 1D de dos nodos, o de tres con el intermedio en el punto medio del tramo, atado como los extremos. Los nodos que crea la barra son esclavos (ADR 0004): $\mathbf u = \mathbf T\,\mathbf u_{\text{libre}} + \mathbf g$, donde las filas de $\mathbf T$ son los $N_i^e(\boldsymbol\xi_b)$. El sistema reducido es $\mathbf T^\top\mathbf K\mathbf T$. En un tramo dentro de un anfitrión lineal, $\mathbf T^\top \mathbf k_t \mathbf T$ coincide con la rigidez de la barra embebida de de Borst et al. (2012, ec. 6.140).

### 7. Funciones de forma
Las del anfitrión, evaluadas en $\boldsymbol\xi_b$, y las lineales del elemento 1D a lo largo del tramo.

### 8. Localización de los puntos
Para cada nodo de la barra, y para los cortes con los bordes:
- se buscan los anfitriones candidatos por su caja envolvente;
- se resuelve la inversa del mapeo isoparamétrico por Newton;
- se comprueba la pertenencia al dominio de referencia, con tolerancia.

El punto de salida de un tramo se obtiene sobre la recta, con la misma pertenencia. Un nodo de la barra fuera de todo anfitrión es un error con mensaje claro.

### 9. Rigidez y fuerzas internas
Las del elemento 1D (su spec); la transmisión al sólido la hace $\mathbf T$. Las fuerzas de las restricciones de embebido en los nodos de la barra son las fuerzas de adherencia concentradas que el sólido transmite a la barra: salen en `constraint_forces`, separadas de las reacciones (ADR 0004 §3).

### 10. Masa
La del elemento 1D ($\rho_s A_s L$, consistente o concentrada), transmitida al sólido por $\mathbf T^\top\mathbf M\mathbf T$.

### 11. Cuadratura
La del elemento 1D. La restricción no se integra.

### 12. Adherencia-deslizamiento (opción `bond`)
Decisiones físicas del usuario del 2026-11-10: adherencia distribuida, integración nodal por omisión con Gauss como opción, barra cuadrática recomendada, leyes `BondSlipModelCode` y `BondSlipMultilinear` (normal elástica desacoplada, descarga R1).

Con `bond`, los nodos de la barra **no** se atan al anfitrión. Cada tramo, que está dentro de un solo anfitrión (la barra se corta en sus bordes, §8), se une a él por una **interfaz distribuida** a lo largo del tramo, gobernada por un material de interfaz de línea (ADR 0033).

- **Desplazamiento relativo** (barra menos concreto) en el punto de coordenada $\xi \in [0, 1]$ del tramo:
$$\Delta\mathbf u(\xi) = \sum_a N^b_a(\xi)\,\mathbf u^b_a - \sum_i N^h_i\big(\boldsymbol\xi_h(\xi)\big)\,\mathbf u^h_i,$$
  con $N^b$ las funciones de forma del elemento de barra (2 o 3 nodos) y $N^h$ las del anfitrión, evaluadas en las coordenadas naturales del punto, que da la localización de §8.
- **Sistema local del tramo**: $\mathbf t$ es el eje. En 2D, $\mathbf n = \mathbf e_z \times \mathbf t$. En 3D, $\mathbf n_1, \mathbf n_2$ es una base ortonormal cualquiera del plano normal, porque la ley es isótropa en ese plano. Componentes $\boldsymbol\Delta = \mathbf R\,\Delta\mathbf u$ con $\mathbf R = [\mathbf n\ (\mathbf n_1, \mathbf n_2);\ \mathbf t]$, normales primero.
- **Fuerza por unidad de longitud** que el concreto ejerce sobre la barra: $\mathbf q = -\mathbf R^\top\mathbf W\,\mathbf t(\boldsymbol\Delta)$, con $\mathbf W = \mathrm{diag}(d_b,\ [d_b],\ p)$. La tracción normal actúa sobre el ancho proyectado $d_b$ y la de adherencia sobre el perímetro $p$. Por omisión $p = \pi d_b$; en un paquete de barras se da su perímetro.
- **Forma débil**: la interfaz añade $\delta W_{int} = \int_{tramo} \delta\Delta\mathbf u^\top\,\mathbf R^\top\mathbf W\,\mathbf t(\boldsymbol\Delta)\,ds$. La barra recibe $-\mathbf q$ y el concreto $+\mathbf q$ (acción y reacción).
- **Discretización**, con $\mathbf B_g = [\mathbf N^b(\xi_g)\otimes\mathbf I_d,\ -\mathbf N^h(\boldsymbol\xi_{h,g})\otimes\mathbf I_d]$:
$$\mathbf F_{int} = \sum_g w_g L\ \mathbf B_g^\top\mathbf R^\top\mathbf W\,\mathbf t_g, \qquad \mathbf K = \sum_g w_g L\ \mathbf B_g^\top\mathbf R^\top\mathbf W\,\mathbf D_g\,\mathbf R\,\mathbf B_g,\qquad \mathbf D_g = \partial\mathbf t/\partial\boldsymbol\Delta.$$
- **Integración**:
  - *nodal* (por omisión, Lobatto en los nodos del tramo): pesos $\tfrac12, \tfrac12$ en barras de 2 nodos y $\tfrac16, \tfrac16, \tfrac23$ en las de 3 nodos $[i, j, m]$. La matriz de la interfaz no acopla nodos de la barra y no oscila (de Borst et al. 2012, §13.1, pp. 431–434).
  - *Gauss* (opción): 2 puntos en barras de 2 nodos y 3 en las de 3. Oscila cuando $\lambda h > \sqrt 6$, con $\lambda = \sqrt{k\,p/(E_s A_s)}$.
- **Barra recomendada**: `Truss2D3`/`Truss3D3` con integración nodal. Con los mismos nodos de barra da un deslizamiento entre 10 y 50 veces más exacto que la lineal, sin oscilar.
- **Nodo en un borde común**: un nodo de la barra en el borde de dos anfitriones recibe la contribución de los dos tramos, cada uno con su anfitrión. El campo del sólido es continuo en el borde.
- **Masa**: la de la barra; la interfaz no tiene masa.
- **Salida**: en cada punto de integración, posición, $s$, $\tau$, $\Delta u_n$, $t_n$ y la fuerza de adherencia por unidad de longitud $p\,\tau$.
- **Elementos de interfaz**: los genera el refuerzo, uno por tramo, como genera hoy las MPC. No se registran ni se citan en el YAML por `type`: su nodo es la barra con su anfitrión.


---

## Contrato de implementación

```yaml
name: EmbeddedReinforcement
kind: constraint
status: validated        # draft → implemented → validated

interface:
  dof_names: [ux, uy]    # [ux, uy, uz] en 3D; sólo traslaciones
  hosts: "todo elemento de sólido 2D/3D isoparamétrico (Tri3, Tri6, Quad4, Quad8, Quad9, Tet4, Tet10, Hex8, Hex20, Hex27) y los de módulos de usuario que declaren sus funciones de forma"
  bar_elements: "Truss2D, Truss2DCorot, Truss3D, Truss3DCorot, Cable2DCorot, Cable3DCorot, Truss2D3, Truss3D3"

admission:
  textbooks:
    - de_borst_crisfield: "§6.8.2, ecs. 6.136-6.140 (refuerzo embebido con adherencia perfecta)"
  general_codes:
    - abaqus: "*EMBEDDED ELEMENT (nodos de la barra restringidos al campo del anfitrión)"
    - opensees: "ASDEmbeddedNodeElement"
    - ansys: "REINF264 (barra integrada en el anfitrión)"
    - diana: "embedded bar reinforcement"

parameters:
  - { name: points, type: "lista de puntos", required: true, desc: "Polilínea de la barra en la configuración de referencia" }
  - { name: area, type: float, required: true, desc: "Área de la barra A_s" }
  - { name: material, type: "material 1D", required: true, desc: "Ley de la barra" }
  - { name: element, type: str, required: false, desc: "Elemento 1D de dos o tres nodos (Truss2D por omisión en 2D, Truss3D en 3D)" }
  - { name: hosts, type: "grupo o lista de elementos", required: false, desc: "Anfitriones admitidos (todos los sólidos por omisión)" }
  - { name: subdivisions, type: int, required: false, desc: "Partes iguales por tramo (1 por omisión)" }
  - { name: bond, type: "mapping", required: false, desc: "Adherencia-deslizamiento (§12): {law: id en interface_materials (GEOMETRY line), diameter: d_b, perimeter: p (por omisión π·d_b), integration: nodal|gauss (nodal por omisión)}. Sin bond, adherencia perfecta." }

conventions:
  sign: "N > 0 ⇔ tracción en la barra (la del elemento 1D)"
  voigt: "N/A"
  node_orientation: "la de la polilínea"

validity:
  - adherencia perfecta por omisión; con `bond`, adherencia-deslizamiento con una ley local τ–s (§12)
  - "el área de concreto que ocupa la barra no se descuenta: sobrante ρ/((1−ρ)+nρ) en rigidez elástica (§3)"
  - la barra sigue el campo regular del anfitrión; dentro de un anfitrión con discontinuidad embebida ve la deformación media del elemento, apertura de grieta incluida

out_of_scope:
  - degradación cíclica de la adherencia y efecto de la fluencia de la barra sobre τ (leyes de módulo de usuario)
  - barras con rigidez a flexión (elementos de marco embebidos)
  - barras curvas (se aproximan con la polilínea)
  - nodos de la barra fuera del sólido

acceptance:
  verification:
    - name: deformacion_uniforme_angulo_arbitrario
      setup: "Tri3, Quad4, Hex8 y Tet4; campo de deformación uniforme impuesto en todos los nodos del sólido; barra recta a un ángulo cualquiera que cruza varios elementos"
      expect: "ε de cada tramo = tᵀ ε t; rigidez añadida = E_s A_s t⊗t⊗t⊗t por unidad de longitud"
      tol_rel: 1.0e-10
    - name: modos_rigidos
      setup: "sólido + barra libres, 2D y 3D"
      expect: "3 (2D) o 6 (3D) modos de energía nula, todos de cuerpo rígido"
      tol_rel: 1.0e-10
    - name: malla_de_barra_independiente
      setup: "Tri3; la misma barra con nodos alineados con la malla y con nodos arbitrarios"
      expect: "mismo resultado"
      tol_rel: 1.0e-10
  specific:
    - name: placa_traccion_uniaxial
      setup: "placa en esfuerzo plano con una barra paralela a x en toda su longitud, desplazamiento impuesto en los bordes"
      expect: "rigidez axial efectiva E_c A_g + E_s A_s (área bruta: sin descuento)"
      tol_rel: 1.0e-10
    - name: fluencia_de_la_barra
      setup: "barra Elastoplastic1D en un sólido de rigidez despreciable"
      expect: "la placa cede cuando la barra alcanza A_s σ_y"
      tol_rel: 1.0e-8
    - name: viga_seccion_transformada
      setup: "viga simplemente apoyada en esfuerzo plano con refuerzo inferior, elástica, mallas sucesivas"
      expect: "flecha central → la de la sección transformada sin descuento (Euler-Bernoulli + cortante de Timoshenko), convergencia con la malla"
      tol_rel: 1.0e-2
    - name: fuerzas_de_adherencia
      setup: "barra en un sólido bajo carga cualquiera"
      expect: "Σ fuerzas de las restricciones de embebido sobre la barra + fuerzas en sus nodos = 0 (equilibrio de la barra)"
      tol_rel: 1.0e-10
    - name: masa
      setup: "modal de sólido + barra"
      expect: "masa total = la del sólido + ρ_s A_s L_b"
      tol_rel: 1.0e-12
    - name: yaml
      setup: "el mismo modelo en Python y en YAML, 2D y 3D"
      expect: "resultados idénticos"
      tol_rel: 1.0e-12
    - name: bond_extraccion_lineal
      setup: "extracción: barra en concreto de rigidez muy alta, ley lineal τ = k·s (BondSlipMultilinear de un tramo), Truss2D3 con integración nodal"
      expect: "deslizamiento igual a la solución cerrada s(x) = P/(E_s A_s λ)·cosh(λ(ℓ−x))/sinh(λℓ) con error decreciente al refinar la barra (valores medidos en el Diálogo)"
      tol_rel: 1.0e-3
    - name: bond_sin_oscilacion
      setup: "la misma extracción con λh = 3.5"
      expect: "integración nodal: τ monótona a lo largo de la barra; Gauss: τ cambia de signo (propiedad documentada)"
    - name: bond_limite_adherencia_perfecta
      setup: "k0 y k_n multiplicados por 1e4 frente a un valor de servicio, viga con refuerzo"
      expect: "desplazamientos del sólido → los de la adherencia perfecta (MPC)"
      tol_rel: 1.0e-4
    - name: bond_equilibrio
      setup: "barra con adherencia-deslizamiento bajo carga cualquiera"
      expect: "Σ fuerzas de la interfaz sobre la barra + fuerzas aplicadas en sus nodos = 0; la reacción sobre el concreto es la opuesta"
      tol_rel: 1.0e-10
    - name: bond_traslacion_rigida
      setup: "traslación rígida de barra y sólido juntos"
      expect: "sin fuerzas de interfaz"
      tol_abs: 1.0e-12
    - name: bond_anclaje_uniforme
      setup: "una interfaz con el anfitrión fijo y la barra desplazada en bloque (deslizamiento uniforme) siguiendo el ciclo de la Fig. 5.3b, confirmando el estado en cada paso"
      expect: "fuerza de extracción P = τ(s)·p·ℓ con la τ de la ley en cada instante, envolvente, descarga, fricción y recarga incluidas"
      tol_rel: 1.0e-8
    - name: bond_quiebre
      setup: "barra de dos tramos con un quiebre de la polilínea, en tracción"
      expect: "en el nodo del quiebre la fuerza de la interfaz sobre la barra equilibra el cambio de dirección del axial: N1·t1 − N2·t2 (no nula)"
      tol_rel: 1.0e-8
    - name: bond_3d_y_yaml
      setup: "extracción en un sólido Hex8 y el mismo modelo 2D y 3D en YAML"
      expect: "equilibrio de la barra y resultados del YAML idénticos a Python"
      tol_rel: 1.0e-10

references:
  - "de Borst R., Crisfield M. A., Remmers J. J. C., Verhoosel C. V. (2012). Non-linear Finite Element Analysis of Solids and Structures, 2.ª ed., Wiley, §6.8.2. doi:10.1002/9781118375938"
  - "Elwi A. E., Hrudey T. M. (1989). Finite element model for curved embedded reinforcement. J. Eng. Mech. 115(4):740-754. doi:10.1061/(asce)0733-9399(1989)115:4(740)"
  - "Ngo D., Scordelis A. C. (1967). Finite element analysis of reinforced concrete beams. ACI J. 64(3):152-163 (enlaces de adherencia)."
  - "Barzegar F., Maddipudi S. (1994). Generating reinforcement in FE modeling of concrete structures. J. Struct. Eng. 120(5):1656-1662. doi:10.1061/(asce)0733-9445(1994)120:5(1656)"
```

---

## Implementación

- **Archivo**: [solidum/bc/embedded.py](../../solidum/bc/embedded.py); localización de puntos en [solidum/mesh/point_location.py](../../solidum/mesh/point_location.py).
- **Clase**: `EmbeddedReinforcement` (registrada con `@ConstraintRegistry.register`); función `solidum.embed_reinforcement`; resultado `EmbeddedBar` (`node_ids`, `element_ids`, `host_ids`, `natural_coordinates`, `length`, `bond_forces(result)`). YAML: bloque `constraints` con `type: EmbeddedReinforcement`. Arquitectura: [ADR 0032](../adr/0032-restricciones-generadas-y-refuerzo-embebido.md).
- **Tests**: [tests/test_embedded_reinforcement.py](../../tests/test_embedded_reinforcement.py), un test por criterio de aceptación más la inversa del mapeo en los diez sólidos con mallas curvas, el localizador y los errores (barra fuera del sólido, tipo o parámetro desconocidos en el YAML).
- **Notas de traducción**:
  - El punto donde la recta sale de un anfitrión se busca por bisección (a 1e-12 del anfitrión más pequeño). Los tramos más cortos que 1e-9 de ese tamaño se funden con el vecino.
  - La viga de la aceptación se apoya con $u_y = 0$ en toda la sección extrema: el apoyo en un nodo produce una deformación local singular que crece al refinar (+2.1 % con Q8 40×8) y no es del refuerzo.
- **Adherencia-deslizamiento (§12)**: [solidum/bc/bond_interface.py](../../solidum/bc/bond_interface.py) · `BondSlipInterface` (elemento generado, sin registro; `VTK_EXPORT = False`, sus resultados se leen con `interface_state` o con `EmbeddedBar.bond_state(domain, U)`); `Bond` (dataclass de `solidum/bc/embedded.py`, también en la raíz del paquete). Tests: [tests/test_embedded_bond_slip.py](../../tests/test_embedded_bond_slip.py) (casos `bond_*`).

---

## Diálogo

- **2026-11-10** · Decisión del usuario: **formulación A**, nodos de la barra atados al campo interpolado del anfitrión y barra cortada en los bordes. La alternativa B (barra sin nodos propios, con $\varepsilon_b = \mathbf t\cdot\nabla\mathbf u\cdot\mathbf t$ integrada en el anfitrión) es idéntica en Tri3 y en Quad4 con la barra paralela a una arista. En los demás casos medidos difiere ≤ 0.18 %, frente a un error de malla de 14 %. A prepara el deslizamiento barra-sólido: basta cambiar la restricción por una ley τ–s entre los nodos de la barra y el sólido.
- **2026-11-10** · Decisión del usuario: **no se descuenta el área de concreto desplazada** (§3). Alternativas descartadas:
  - descuento exacto con el esfuerzo del anfitrión en la barra: sólo es exacto con esfuerzo constante por elemento; en los demás anfitriones el estado de un material con historia vive en sus puntos de Gauss;
  - descuento con $E_s - E_c$: tras el agrietamiento resta un esfuerzo que el concreto no lleva, y la fuerza de la barra en fluencia queda ~15 % baja.
- **2026-11-10** · Validada: los criterios de aceptación pasan, con la viga Q8 20×4 a +0.06 % de la sección transformada sin descuento (§3), lo que confirma también la decisión del área.
- **2026-11-10** · Clasificación: estándar (programa principal). Texto: de Borst et al. (2012) §6.8.2, leído (ecs. 6.136–6.140). Códigos: Abaqus y OpenSees (variante A), ANSYS y DIANA (variante B).
- **2026-11-10** · **Barras de tres nodos** (`Truss2D3`, `Truss3D3`; decisión del usuario sobre la adherencia-deslizamiento: barra cuadrática para la adherencia-deslizamiento). Cada tramo recibe un nodo intermedio en su punto medio, atado al campo del anfitrión como los extremos; `EmbeddedBar.node_ids` los intercala en el orden de la polilínea. Medido (`tests/test_embedded_reinforcement.py`): donde el campo del anfitrión es lineal a lo largo de la barra (Tri3 con cualquier inclinación, Quad4 con la barra paralela a un lado) la solución es la de la barra de dos nodos (desplazamientos a 1e-10, axial a 1e-7 por redondeo en tramos cortos); con anfitriones cuadráticos (Tri6, Tet10) y un campo cuadrático impuesto, la deformación de la barra de tres nodos es tᵀ·∇u·t a 1e-14 en todo el tramo, y la de dos nodos sólo en el punto medio de cada tramo. En Quad4 con la barra inclinada el campo bilineal es cuadrático a lo largo de la recta y la barra de tres nodos lo sigue (sin prueba de igualdad).
- **2026-11-10** · **Ampliación en BORRADOR (§12): adherencia-deslizamiento**, pendiente de la aprobación del usuario. Las decisiones físicas de la adherencia-deslizamiento son: distribuida, integración nodal con Gauss opcional, barra cuadrática recomendada, leyes `BondSlipModelCode` y `BondSlipMultilinear`. Los elementos de interfaz los genera el refuerzo (uno por tramo) y no se registran, porque su número de nodos es variable (barra + anfitrión) y nadie los cita por `type`. Los casos `bond_*` de `acceptance` son de esta ampliación.
- **2026-11-10** · §12 y casos `bond_*` **aprobados por el usuario** («si»).
- **2026-11-10** · **§12 implementada y validada** (casos `bond_*` en verde). Medido (extracción con el concreto fijo, ley lineal k = 1e10 Pa/m, d_b = 16 mm, ℓ = 1 m, error máximo del deslizamiento frente a la solución cerrada, relativo a s(0)): `Truss2D3` nodal 2.7e-5 / 1.7e-6 / 1.1e-7 con 5/10/20 divisiones del concreto (orden 4); `Truss2D` nodal 1.5e-2 / 3.9e-3 / 9.7e-4 (orden 2); con k = 1e11, 2.3e-3 / 1.6e-4 / 1.0e-5 frente a 1.3e-1 / 3.7e-2 / 9.6e-3. Con λh = 3.5 (k = 1e12): nodal τ monótona; Gauss invierte el signo de τ (−10.6 % de τ_máx). En 3D (Hex8 10×2×2, `Truss3D3`) el error es < 1e-4. Con rigideces 1e4 veces las de servicio, el sólido se mueve como con adherencia perfecta (1e-4). Hallazgo al implementar: la ley daba tangente nula en el estado virgen con s = 0 (empate de cotas), lo que la red de seguridad del ADR 0019 detectó como mecanismo; corregido en la ley (ver su spec).
