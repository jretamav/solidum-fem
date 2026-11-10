# ELEMENTO FINITO SÓLIDO UNIDIMENSIONAL DE SEGUNDO ORDEN

> Orden de trabajo. El usuario escribe la **especificación física**, la **formulación numérica** y el **contrato**. La IA rellena **implementación** y responde en **diálogo** durante el trabajo.

---

## Especificación física

### 0. Descripción general
Elemento sólido 1D de segundo orden, recto, inmerso en el plano. Tres nodos: los extremos $i$, $j$ y un nodo intermedio $m$ en el punto medio. Transmite únicamente fuerza axial; aproximación $C^0$ cuadrática del desplazamiento axial. Es la barra de [`Truss2D`](Truss2D.md) con un grado más de interpolación.

Motivo de su incorporación (2026-11-10): la barra de refuerzo con adherencia-deslizamiento. Con los mismos nodos de barra, la barra cuadrática con integración nodal de la interfaz (Lobatto de 3 puntos) da un deslizamiento entre 10 y 50 veces más exacto que la lineal y no oscila (datos en el Diálogo). Es un elemento general: sirve en cualquier modelo de barras.

### 1. Cinemática de desplazamientos
Interpolación cuadrática del desplazamiento axial a lo largo del eje $s \in [0, L]$:
$$u(s) = a_0 + a_1 s + a_2 s^2.$$

### 2. Cinemática de deformaciones
Única componente (axial, en el sistema local de la barra), lineal en $s$:
$$\varepsilon_{ss}(s) = \frac{du}{ds}.$$

### 3. Ecuación constitutiva
La del material 1D en cada punto de integración: $\sigma_{ss} = \sigma(\varepsilon_{ss}, \text{estado})$, con módulo tangente $E_t$.

### 4. Equilibrio — forma fuerte
$$-\frac{d}{ds}\!\left(\sigma A\right) = b_s(s), \quad s \in (0, L),$$
con condiciones de contorno Dirichlet o Neumann ($\sigma A = \bar F$) en los extremos.

### 5. Equilibrio — forma débil
$$\int_0^L \frac{d\,\delta u}{ds}\, \sigma A\, ds \;=\; \int_0^L \delta u\, b_s\, A\, ds \;+\; \bigl[\delta u\, \bar F\bigr]_{\partial \Omega}, \quad \forall\, \delta u \in V_0.$$

---

## Formulación numérica (FEM)

### 6. Discretización
Nodos en el orden $[i, j, m]$ (extremos primero, intermedio al final, como los sólidos cuadráticos de Solidum, Gmsh y VTK), en $s = 0$, $s = L$ y $s = L/2$. El elemento es recto: el nodo $m$ debe estar en el punto medio de $i$–$j$ (tolerancia $10^{-8}L$; si no, `ValueError`). Cosenos directores del eje en globales:
$$\mathbf e = (c,\; s_\theta) = \frac{\mathbf x_j - \mathbf x_i}{L}.$$

### 7. Funciones de forma
Con $\xi = s/L \in [0, 1]$:
$$N_i = (1-\xi)(1-2\xi), \qquad N_j = \xi(2\xi-1), \qquad N_m = 4\xi(1-\xi).$$

### 8. Matriz deformación-desplazamiento
Local (1 DOF axial por nodo, orden $[i, j, m]$):
$$\mathbf B_{\text{loc}}(\xi) = \tfrac{1}{L}\,[\,4\xi-3,\;\; 4\xi-1,\;\; 4-8\xi\,].$$
Global ($\mathbf u_e = [u_x^{(i)}, u_y^{(i)}, u_x^{(j)}, u_y^{(j)}, u_x^{(m)}, u_y^{(m)}]^\top$): $\mathbf B = \mathbf B_{\text{loc}} \otimes \mathbf e^\top$, $\varepsilon = \mathbf B\,\mathbf u_e$.

### 9. Rigidez elemental
$$\mathbf K_e = \int_0^L \mathbf B^\top E_t A\, \mathbf B\, ds \;\approx\; \sum_{g=1}^{2} w_g L\; \mathbf B_g^\top E_{t,g} A\, \mathbf B_g .$$
Con $E_t A$ constante la regla de Gauss de 2 puntos es exacta (integrando cuadrático) y da, en el orden $[i, j, m]$:
$$\mathbf K_e = \frac{E A}{3L}\begin{bmatrix} 7 & 1 & -8\\ 1 & 7 & -8\\ -8 & -8 & 16\end{bmatrix} \otimes \mathbf e\,\mathbf e^\top .$$

### 10. Fuerzas internas nodales
$$\mathbf F_{\text{int}} = \sum_{g=1}^{2} w_g L\; \mathbf B_g^\top\, \sigma_g A .$$

### 11. Cuadratura y estado
Gauss de 2 puntos, $\xi_g = \tfrac12 \mp \tfrac{1}{2\sqrt3}$, $w_g = \tfrac12$. Cada punto lleva su estado del material. Longitud característica de cada punto (ADR 0028): $h_g = w_g L = L/2$, como en los marcos.

### 12. Fuerza interna (API pública)
El axial del elemento es lineal: $N(\xi) = \sigma(\xi)A$. Se exponen $N_i$ y $N_j$ por **extrapolación lineal** de los dos puntos de Gauss a los extremos (con material elástico es el campo exacto del elemento). No se usan las fuerzas nodales de extremo: con $N$ lineal, $F_{\text{int},i}\cdot\mathbf e = -(5N_i + N_j)/6 \neq -N_i$.

### 13. Masa
Consistente, en el orden $[i, j, m]$ y por componente traslacional:
$$\mathbf M_e = \frac{\rho A L}{30}\begin{bmatrix} 4 & -1 & 2\\ -1 & 4 & 2\\ 2 & 2 & 16\end{bmatrix} \otimes \mathbf I_2 .$$
Concentrada: $\rho A L\,(\tfrac16, \tfrac16, \tfrac23)$ por componente. La regla de Lobatto de 3 puntos (Simpson) y el reparto HRZ dan la misma diagonal.

### 14. Carga de cuerpo
Uniforme: $\mathbf f = b\,A L\,(\tfrac16, \tfrac16, \tfrac23)$ por componente global. Campo $b(\mathbf x)$: Gauss de 4 puntos (ADR 0022), como `Truss2D`.

---

## Contrato de implementación

```yaml
name: Truss2D3
kind: element
status: draft            # draft → implemented → validated

interface:
  dof_names: [ux, uy]
  n_nodes: 3
  strain_dim: 1
  n_integration_points: 2

admission:                # criterio de admisión al programa principal (Reglas.md §1)
  textbooks:
    - bathe: "§6.3.3, p. 543, fig. 6.3 (armadura de 2 a 4 nodos)"
  general_codes:
    - abaqus: "T2D3"

parameters:
  - { name: A, type: float, required: true, desc: Área de la sección transversal }

material_contract:
  signature: "compute_state(ε, state) -> (σ, E_tangent, state')"
  strain_kind: axial escalar

conventions:
  sign: "σ > 0 ⇔ ε > 0 (elongación); N > 0 tracción"
  voigt: "N/A (escalar)"
  node_orientation: "[i, j, m]: extremos y nodo intermedio; K_e invariante al intercambiar i y j"

validity:
  - "|ε| ≲ 1e-2"
  - pequeños desplazamientos y rotaciones
  - elemento recto con el nodo intermedio en el punto medio
  - cargas exclusivamente axiales

out_of_scope:
  - barras curvas (Bathe admite geometría curva; aquí no)
  - grandes desplazamientos (variante corotacional)
  - flexión, cortante, pandeo

acceptance:
  - name: soluciones manufacturadas (ADR 0022)
    setup: "línea recta inclinada; solución exacta impuesta en un extremo y su fuerza interna como carga nodal en el otro; carga de línea como fuerza de cuerpo b(x); n = 2, 4, 8, 16 elementos"
    expect: "con u* polinómica de grado ≤ 2 (carga constante) el campo es exacto en todo punto; con cargas polinómicas de grado ≤ 4, error nodal de redondeo (exactitud nodal, Tong 1969); con u* suave, orden ≥ 3 − 0.1 en L² y ≥ 2 − 0.1 en energía entre los dos últimos niveles, error decreciente"
    tol_abs: 0.1
  - name: barra axial bajo carga puntual
    setup: "empotrada en s = 0, carga F axial en s = L, 1 y 3 elementos, barra inclinada"
    expect: "u(L) = F·L/(E·A); N = F en los dos puntos de Gauss y en los extremos"
    tol_rel: 1.0e-10
  - name: carga axial uniforme
    setup: "empotrada en s = 0, libre en s = L, carga de cuerpo axial uniforme b"
    expect: "u(s) = b(Ls − s²/2)/E en todo punto y N(s) = bA(L − s) en los extremos de cada elemento, con 1 elemento"
    tol_rel: 1.0e-10
  - name: matrices
    setup: "elemento inclinado, E, A, ρ y L distintos de 1"
    expect: "K_e y M_e iguales a las formas cerradas §9 y §13; K_e simétrica de rango 2 con 4 modos de energía nula (traslaciones y movimientos transversales); masa total ρAL en cada dirección"
    tol_rel: 1.0e-12
  - name: frecuencias axiales
    setup: "barra empotrada-libre, 4 elementos, masa consistente"
    expect: "primera frecuencia axial ω₁ = (π/2)·√(E/ρ)/L con error relativo ≤ 2e-5 (medido 1.64e-5; la barra lineal con los mismos 9 nodos da 4.0e-4)"
    tol_rel: 2.0e-5
  - name: material no lineal
    setup: "barra elastoplástica perfecta bajo carga en el extremo hasta la fluencia"
    expect: "fluencia en N = A·σ_y en los dos puntos de Gauss; estado por punto confirmado al converger"
    tol_rel: 1.0e-8

references:
  - "Bathe K.-J., Finite Element Procedures, 2.ª ed., 2014, §6.3.3 (p. 543) y fig. 5.3"
  - "Abaqus Analysis User's Guide, elementos de armadura T2D3/T3D3"
```

---

## Implementación

*Rellena la IA tras programar.*

---

## Diálogo

- **2026-11-10** · Origen: decisión del usuario sobre la integración de la interfaz de adherencia (S3.2). Con los mismos nodos de barra (lineal con 2n tramos frente a cuadrática con n elementos), error máximo del deslizamiento en un ensayo de extracción (d_b = 16 mm, ley lineal τ = k_b·s): concreto rígido, k_b = 10¹⁰ Pa/m, 11 nodos: lineal nodal 1.52 %, lineal Gauss 0.52 %, cuadrática nodal (Lobatto 3) 0.04 %, cuadrática Gauss 0.02 %; concreto Tri3 10×3, k_b = 10¹²: 2.36 / 0.56 / 0.08 / 0.28 %. Oscilación de τ: cuadrática nodal ≥ −1 % de τ_máx sin inversiones; Gauss (lineal o cuadrática) hasta −20 %. **Decidido por el usuario (2026-11-10): añadir la barra de 3 nodos; para la adherencia se recomienda barra cuadrática con integración nodal; la barra lineal se conserva y Gauss queda como opción en los dos órdenes.**
- **2026-11-10** · Formas cerradas de §9, §10 y §13 deducidas con SymPy (Gauss de 2 puntos = integral exacta de K_e con E_tA constante). Frecuencias axiales de la barra empotrada-libre (consistente, error relativo del primer modo): 3.75e-3 / 2.56e-4 / 1.64e-5 / 1.03e-6 con 1/2/4/8 elementos (orden 4); lineal con los mismos nodos: 6.44e-3 / 1.61e-3 / 4.02e-4.
