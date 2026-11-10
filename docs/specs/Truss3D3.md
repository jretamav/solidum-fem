# ELEMENTO FINITO SÓLIDO UNIDIMENSIONAL DE SEGUNDO ORDEN EN EL ESPACIO

> Orden de trabajo. El usuario escribe la **especificación física**, la **formulación numérica** y el **contrato**. La IA rellena **implementación** y responde en **diálogo** durante el trabajo.

---

## Especificación física

### 0. Descripción general
Versión espacial de [`Truss2D3`](Truss2D3.md): barra recta de tres nodos $[i, j, m]$, con $m$ en el punto medio, que transmite sólo fuerza axial, con desplazamiento axial cuadrático. Toda la física (§1–§5) es la de `Truss2D3`.

---

## Formulación numérica (FEM)

Igual que `Truss2D3` (§6–§14) con tres componentes traslacionales por nodo:

- Cosenos directores $\mathbf e = (\mathbf x_j - \mathbf x_i)/L \in \mathbb R^3$.
- $\mathbf B = \mathbf B_{\text{loc}} \otimes \mathbf e^\top$, $\mathbf u_e = [\mathbf u^{(i)}, \mathbf u^{(j)}, \mathbf u^{(m)}]^\top$ (9 DOF).
- $\mathbf K_e = \dfrac{EA}{3L}\begin{bmatrix} 7 & 1 & -8\\ 1 & 7 & -8\\ -8 & -8 & 16\end{bmatrix} \otimes \mathbf e\,\mathbf e^\top$ con $E_tA$ constante; Gauss de 2 puntos en general.
- $\mathbf M_e = \dfrac{\rho A L}{30}\begin{bmatrix} 4 & -1 & 2\\ -1 & 4 & 2\\ 2 & 2 & 16\end{bmatrix} \otimes \mathbf I_3$; concentrada $\rho A L(\tfrac16, \tfrac16, \tfrac23)$ por componente.
- Fuerza interna pública: $N_i$, $N_j$ por extrapolación lineal desde los dos puntos de Gauss (convención 3D de Reglas §5: $N > 0$ tracción, en ejes locales).

---

## Contrato de implementación

```yaml
name: Truss3D3
kind: element
status: validated        # draft → implemented → validated

interface:
  dof_names: [ux, uy, uz]
  n_nodes: 3
  strain_dim: 1
  n_integration_points: 2

admission:                # criterio de admisión al programa principal (Reglas.md §1)
  textbooks:
    - bathe: "§6.3.3, p. 543, fig. 6.3 (armadura de 2 a 4 nodos en el espacio)"
  general_codes:
    - abaqus: "T3D3"

parameters:
  - { name: A, type: float, required: true, desc: Área de la sección transversal }

material_contract:
  signature: "compute_state(ε, state) -> (σ, E_tangent, state')"
  strain_kind: axial escalar

conventions:
  sign: "σ > 0 ⇔ ε > 0 (elongación); N > 0 tracción"
  voigt: "N/A (escalar)"
  node_orientation: "[i, j, m]: extremos y nodo intermedio"

validity:
  - "|ε| ≲ 1e-2"
  - pequeños desplazamientos y rotaciones
  - elemento recto con el nodo intermedio en el punto medio

out_of_scope:
  - barras curvas
  - grandes desplazamientos
  - flexión, cortante, torsión, pandeo

acceptance:
  - name: soluciones manufacturadas (ADR 0022)
    setup: "línea recta oblicua en el espacio; mismas condiciones que Truss2D3; n = 2, 4, 8, 16"
    expect: "campo exacto con u* de grado ≤ 2; error de redondeo en los nodos extremos con cualquier carga polinómica integrada exactamente y en el intermedio con u* de grado ≤ 3 (ver la corrección del 2026-11-10 en Truss2D3); orden ≥ 3 − 0.1 en L² y ≥ 2 − 0.1 en energía"
    tol_abs: 0.1
  - name: barra axial bajo carga puntual (3D)
    setup: "barra oblicua empotrada en un extremo, carga F axial en el otro, 1 y 3 elementos"
    expect: "u(L) = F·L/(E·A); N = F en los puntos de Gauss y en los extremos"
    tol_rel: 1.0e-10
  - name: matrices
    setup: "elemento oblicuo, E, A, ρ y L distintos de 1"
    expect: "K_e y M_e iguales a las formas cerradas; K_e de rango 2 con 7 modos de energía nula; masa total ρAL en cada dirección; igual a Truss2D3 cuando la barra está en el plano z = 0"
    tol_rel: 1.0e-12

references:
  - "Bathe K.-J., Finite Element Procedures, 2.ª ed., 2014, §6.3.3 (p. 543)"
  - "Abaqus Analysis User's Guide, elementos de armadura T2D3/T3D3"
```

---

## Implementación

- **Archivo**: [solidum/elements/truss_quadratic.py](../../solidum/elements/truss_quadratic.py)
- **Clase**: `Truss3D3`, sobre la misma base `_QuadraticTruss` que `Truss2D3`. Acepta nodos con 2 ó 3 coordenadas.
- **Tests**: los de [`Truss2D3`](Truss2D3.md) parametrizados en 3D en [tests/test_truss_quadratic.py](../../tests/test_truss_quadratic.py) (matrices, carga puntual con 1 y 3 elementos oblicuos, igualdad con `Truss2D3` en el plano z = 0, `internal_forces_from_nodal`), [tests/validation/test_mms_frames.py](../../tests/validation/test_mms_frames.py), [tests/test_line_elements.py](../../tests/test_line_elements.py) y el barrido de contrato.

---

## Diálogo

- **2026-11-10** · Abierta junto con [`Truss2D3`](Truss2D3.md) (decisión del usuario del mismo día; motivo y datos allí).
- **2026-11-10** · Implementado y validado: todos los casos de `acceptance` en verde. Promovido `status: validated`.
