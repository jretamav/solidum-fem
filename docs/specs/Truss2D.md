# ELEMENTO FINITO SÓLIDO UNIDIMENSIONAL

> Orden de trabajo. El usuario escribe la **especificación física**, la **formulación numérica** y el **contrato**. La IA rellena **implementación** y responde en **diálogo** durante el trabajo.

---

## Especificación física

### 0. Descripción general
Elemento sólido 1D de primer orden, inmerso en el plano. Dos nodos articulados; transmite únicamente esfuerzo axial. Aproximación $C^0$ del desplazamiento.

### 1. Cinemática de desplazamientos
Interpolación lineal del desplazamiento axial a lo largo del eje $s \in [0, L]$:
$$u(s) = a_0 + a_1 s$$

### 2. Cinemática de deformaciones
Única componente (axial, en el sistema local de la barra):
$$\varepsilon_{ss}(s) = \frac{du}{ds}$$

### 3. Ecuación constitutiva
$$\sigma_{ss} = E\, \varepsilon_{ss}$$

### 4. Equilibrio — forma fuerte
$$-\frac{d}{ds}\!\left(EA\,\frac{du}{ds}\right) = b(s), \quad s \in (0, L),$$
con condiciones de contorno Dirichlet o Neumann ($\sigma A = \bar F$) en los extremos.

### 5. Equilibrio — forma débil
$$\int_0^L \frac{d\,\delta u}{ds}\, EA\, \frac{du}{ds}\, ds \;=\; \int_0^L \delta u\, b\, ds \;+\; \bigl[\delta u\, \bar F\bigr]_{\partial \Omega}, \quad \forall\, \delta u \in V_0.$$

---

## Formulación numérica (FEM)

### 6. Discretización
Dos nodos en $s = 0$ y $s = L$. Cosenos directores del eje en globales:
$$c = \frac{x_2 - x_1}{L}, \qquad s_\theta = \frac{y_2 - y_1}{L}.$$

### 7. Funciones de forma
$$N_1(s) = 1 - \tfrac{s}{L}, \qquad N_2(s) = \tfrac{s}{L}.$$

### 8. Matriz deformación-desplazamiento
Local (1 DOF axial por nodo): $\mathbf B_{\text{loc}} = \tfrac{1}{L}[-1,\; 1]$.
Global ($\mathbf u_e = [u_x^{(1)}, u_y^{(1)}, u_x^{(2)}, u_y^{(2)}]^\top$):
$$\mathbf B = \tfrac{1}{L}\,[\,-c,\; -s_\theta,\; c,\; s_\theta\,], \qquad \varepsilon = \mathbf B\, \mathbf u_e.$$

### 9. Rigidez elemental
Integración analítica exacta ($EA$ constante):
$$\mathbf K_e = \int_0^L \mathbf B^\top EA\, \mathbf B\, ds = \frac{EA}{L}\, \mathbf d^\top \mathbf d, \qquad \mathbf d = [-c,\; -s_\theta,\; c,\; s_\theta].$$

### 10. Fuerzas internas
$$\mathbf F_{\text{int}} = \sigma A\, \mathbf d.$$

### 11. Cuadratura
No aplica (forma cerrada).

---

## Contrato de implementación

```yaml
name: Truss2D
kind: element
status: validated        # draft → implemented → validated

interface:
  dof_names: [ux, uy]
  n_nodes: 2
  strain_dim: 1
  n_integration_points: 1

parameters:
  - { name: A, type: float, required: true, desc: Área de la sección transversal }

material_contract:
  signature: "compute_state(ε, state) -> (σ, E_tangent, state')"
  strain_kind: axial escalar

conventions:
  sign: "tracción positiva (ε > 0 ⇔ σ > 0)"
  voigt: "N/A (escalar)"
  node_orientation: "libre; K_e y F_int invariantes bajo permutación"

validity:
  - "|ε| ≲ 1e-2"
  - pequeños desplazamientos y rotaciones
  - cargas exclusivamente axiales

out_of_scope:
  - flexión, cortante, momentos en extremos
  - grandes desplazamientos
  - pandeo

acceptance:
  - name: barra axial bajo carga puntual
    setup: "empotrada en s=0, carga F axial en s=L"
    expect: "u(L) = F·L / (E·A)"
    tol_rel: 1.0e-10

references:
  - "Bathe K.-J., Finite Element Procedures, §4.2.1"
```

---

## Implementación

- **Archivo**: [fenix/elements/structural.py](../../fenix/elements/structural.py)
- **Clase**: `Truss2D` (registrada vía `@ElementRegistry.register`)
- **Tests**:
  - [tests/test_structural.py](../../tests/test_structural.py) · `TestTruss2D` — verifica `L0`, cosenos directores, entradas de $\mathbf K_e$ (incluyendo simetría), y evaluación de $\mathbf F_{\text{int}}$ sobre desplazamientos conocidos. Los valores numéricos esperados coinciden con $\mathbf K_e = (EA/L)\,\mathbf d\mathbf d^\top$ y $\mathbf F_{\text{int}} = N\mathbf d$ para geometría de prueba $L=5$, $c=0.6$, $s=0.8$.
  - [tests/test_integration.py](../../tests/test_integration.py) · `TestSolversIntegration` — resuelve end-to-end el caso de aceptación (barra empotrada en un extremo, carga axial en el otro) con `NonlinearSolver` y `ArcLengthSolver`; en régimen elástico precedente a la fluencia reproduce $u(L)=FL/(EA)$.
- **Notas de traducción**:
  - La clase incluye una bandera `large_strains: bool` que activa una formulación corotacional (Updated Lagrangian) — **funcionalmente equivalente a `Truss2DCorot` pero empaquetada como opción interna**. Esta duplicación conceptual es deuda del refactor previo a spec-first; conservar como está hasta decidir si se unifica.
  - `L0`, `c`, `s` se calculan una vez en `__init__` sobre la configuración inicial; son los que usa el camino lineal. En `large_strains=True` se recalculan en cada evaluación de estado.
  - El contrato con el material es el estándar del proyecto: `material.compute_state(ε, state) → (σ, Eₜ, state')`.

---

## Diálogo

- **2026-04-21** · Cierre de ciclo retroactivo. La clase `Truss2D` preexistía al flujo spec-first; la spec se creó documentando la formulación ya implementada. Los tests unitarios y de integración satisfacen el único criterio de `acceptance` declarado (analíticamente, vía los valores de $\mathbf K_e$ y end-to-end vía el solver). Promovido `status: validated`.
