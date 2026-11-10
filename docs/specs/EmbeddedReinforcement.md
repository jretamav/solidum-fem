# EmbeddedReinforcement — refuerzo embebido en sólidos con adherencia perfecta

> Orden de trabajo. El usuario escribe la **especificación física**, la **formulación numérica** y el **contrato**. La IA rellena **implementación** y responde en **diálogo** durante el trabajo.

---

## Especificación física

### 0. Descripción general
Barras de refuerzo rectas o poligonales dentro de un sólido 2D (esfuerzo o deformación plana) o 3D, con adherencia perfecta. La barra es una cadena de elementos 1D ordinarios de Solidum (`Truss2D`, `Truss3D`, sus variantes corrotacionales y los cables), con cualquier material 1D registrado. Su malla es independiente de la del sólido: la define una polilínea, no los nodos del sólido. La adherencia perfecta se impone con restricciones afines (ADR 0004): cada nodo de la barra sigue el desplazamiento interpolado del elemento de sólido que lo contiene.

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
Los tramos son elementos 1D de dos nodos. Los nodos que crea la barra son esclavos (ADR 0004): $\mathbf u = \mathbf T\,\mathbf u_{\text{libre}} + \mathbf g$, donde las filas de $\mathbf T$ son los $N_i^e(\boldsymbol\xi_b)$. El sistema reducido es $\mathbf T^\top\mathbf K\mathbf T$. En un tramo dentro de un anfitrión lineal, $\mathbf T^\top \mathbf k_t \mathbf T$ coincide con la rigidez de la barra embebida de de Borst et al. (2012, ec. 6.140).

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

---

## Contrato de implementación

```yaml
name: EmbeddedReinforcement
kind: constraint
status: validated        # draft → implemented → validated

interface:
  dof_names: [ux, uy]    # [ux, uy, uz] en 3D; sólo traslaciones
  hosts: "todo elemento de sólido 2D/3D isoparamétrico (Tri3, Tri6, Quad4, Quad8, Quad9, Tet4, Tet10, Hex8, Hex20, Hex27) y los de módulos de usuario que declaren sus funciones de forma"
  bar_elements: "Truss2D, Truss2DCorot, Truss3D, Truss3DCorot, Cable2DCorot, Cable3DCorot"

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
  - { name: element, type: str, required: false, desc: "Elemento 1D (Truss2D por omisión en 2D, Truss3D en 3D)" }
  - { name: hosts, type: "grupo o lista de elementos", required: false, desc: "Anfitriones admitidos (todos los sólidos por omisión)" }
  - { name: subdivisions, type: int, required: false, desc: "Partes iguales por tramo (1 por omisión)" }

conventions:
  sign: "N > 0 ⇔ tracción en la barra (la del elemento 1D)"
  voigt: "N/A"
  node_orientation: "la de la polilínea"

validity:
  - adherencia perfecta (el deslizamiento barra-sólido es otro componente)
  - "el área de concreto que ocupa la barra no se descuenta: sobrante ρ/((1−ρ)+nρ) en rigidez elástica (§3)"
  - la barra sigue el campo regular del anfitrión; dentro de un anfitrión con discontinuidad embebida ve la deformación media del elemento, apertura de grieta incluida

out_of_scope:
  - deslizamiento barra-sólido
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

references:
  - "de Borst R., Crisfield M. A., Remmers J. J. C., Verhoosel C. V. (2012). Non-linear Finite Element Analysis of Solids and Structures, 2.ª ed., Wiley, §6.8.2. doi:10.1002/9781118375938"
  - "Elwi A. E., Hrudey T. M. (1989). Finite element model for curved embedded reinforcement. J. Eng. Mech. 115(4):740-754. doi:10.1061/(asce)0733-9399(1989)115:4(740)"
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

---

## Diálogo

- **2026-11-10** · Decisión del usuario: **formulación A**, nodos de la barra atados al campo interpolado del anfitrión y barra cortada en los bordes. La alternativa B (barra sin nodos propios, con $\varepsilon_b = \mathbf t\cdot\nabla\mathbf u\cdot\mathbf t$ integrada en el anfitrión) es idéntica en Tri3 y en Quad4 con la barra paralela a una arista. En los demás casos medidos difiere ≤ 0.18 %, frente a un error de malla de 14 %. A prepara el deslizamiento barra-sólido: basta cambiar la restricción por una ley τ–s entre los nodos de la barra y el sólido.
- **2026-11-10** · Decisión del usuario: **no se descuenta el área de concreto desplazada** (§3). Alternativas descartadas:
  - descuento exacto con el esfuerzo del anfitrión en la barra: sólo es exacto con esfuerzo constante por elemento; en los demás anfitriones el estado de un material con historia vive en sus puntos de Gauss;
  - descuento con $E_s - E_c$: tras el agrietamiento resta un esfuerzo que el concreto no lleva, y la fuerza de la barra en fluencia queda ~15 % baja.
- **2026-11-10** · Validada: los criterios de aceptación pasan, con la viga Q8 20×4 a +0.06 % de la sección transformada sin descuento (§3), lo que confirma también la decisión del área.
- **2026-11-10** · Clasificación: estándar (programa principal). Texto: de Borst et al. (2012) §6.8.2, leído (ecs. 6.136–6.140). Códigos: Abaqus y OpenSees (variante A), ANSYS y DIANA (variante B).
