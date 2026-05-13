# ELEMENTO MARCO/VIGA 2D TIMOSHENKO

> Orden de trabajo. El usuario escribe **especificación física**, **formulación numérica** y **contrato**; la IA rellena **implementación** y responde en **diálogo**.

---

## Especificación física

### 0. Descripción general
Elemento 1D inmerso en el plano 2D que modela una **viga** según la teoría de **Timoshenko**. Dos nodos rígidamente conectados, 3 DOFs por nodo ($u_x, u_y, r_z$). Transmite **esfuerzo axial**, **cortante transverso** y **momento flector**. Incluye explícitamente la **deformación por cortante transversal**. Régimen de **linealidad geométrica** (pequeños desplazamientos y rotaciones).

### 1. Hipótesis cinemática de Timoshenko
- Las secciones transversales se mantienen **planas** tras la deformación (igual que Euler-Bernoulli).
- Las secciones **no** están obligadas a ser perpendiculares al eje neutro: se permite una **rotación adicional** respecto a la tangente.

Consecuencia: la rotación de la sección $\theta(s)$ y la pendiente de la elástica $dv/ds$ son **variables independientes**:
$$\gamma(s) = \frac{dv}{ds} - \theta(s) \neq 0,$$
donde $\gamma$ es la distorsión angular (deformación por cortante transverso). Euler-Bernoulli corresponde al caso $\gamma = 0$.

Apropiado para vigas **gruesas, cortas o peraltadas** ($L/h \lesssim 10$) donde la deformación por cortante no es despreciable.

### 2. Campos de desplazamiento
- Axial: $u(s)$ a lo largo del eje local $s \in [0, L]$.
- Transverso: $v(s)$ perpendicular al eje local.
- Rotación de la sección: $\theta(s)$, **independiente** de $v$.

### 3. Medidas de deformación
- Axial: $\varepsilon_{axial}(s) = du/ds$.
- Curvatura (flexión): $\kappa(s) = d\theta/ds$.
- Cortante: $\gamma(s) = dv/ds - \theta(s)$.

### 4. Ecuaciones constitutivas
- Axial: $N = EA\,\varepsilon_{axial}$ (el material entrega $\sigma$ y $E_t$).
- Flexión: $M = EI\,\kappa$.
- Cortante: $V = G A_s\,\gamma$, con $G = E/[2(1+\nu)]$ (módulo de corte) y $A_s$ el **área efectiva de cortante** (menor que $A$ por la distribución no uniforme de la tensión tangencial en la sección).

### 5. Factor de Phi (shear locking)
La formulación estándar con interpolaciones lineales de $v$ y $\theta$ sufre **shear locking** para vigas esbeltas (la rigidez por cortante domina numéricamente aunque $\gamma$ físicamente sea pequeño). La matriz local incorpora un factor correctivo:
$$\Phi = \frac{12\,E\,I}{G\,A_s\,L^2}.$$
Para $\Phi \to 0$ (viga esbelta) la matriz se reduce a la de Euler-Bernoulli.

### 6. Equilibrio — forma débil
Principio de trabajos virtuales (PTV). El elemento calcula solo el lado interno:
$$\delta W_{\text{int}} = \int_0^L (N\,\delta\varepsilon_{axial} + M\,\delta\kappa + V\,\delta\gamma)\,ds = \delta\mathbf u_e^\top\,\mathbf F_{\text{int}}.$$

---

## Formulación numérica (FEM)

### 7. Discretización
Dos nodos ($s = 0$ y $s = L$), 3 DOFs por nodo en coordenadas globales:
$$\mathbf u_e = [u_x^{(1)}, u_y^{(1)}, r_z^{(1)}, u_x^{(2)}, u_y^{(2)}, r_z^{(2)}]^\top.$$

### 8. Sistema local y transformación
Cosenos directores del eje (configuración inicial, fija): $c = (x_2 - x_1)/L$, $s_\theta = (y_2 - y_1)/L$. Matriz de transformación ortogonal $6\times 6$:
$$\mathbf T = \begin{bmatrix}
 c & s_\theta & 0 & 0 & 0 & 0 \\
-s_\theta & c & 0 & 0 & 0 & 0 \\
 0 & 0 & 1 & 0 & 0 & 0 \\
 0 & 0 & 0 & c & s_\theta & 0 \\
 0 & 0 & 0 & -s_\theta & c & 0 \\
 0 & 0 & 0 & 0 & 0 & 1
\end{bmatrix}, \qquad \mathbf u_{\text{local}} = \mathbf T\,\mathbf u_e.$$

### 9. Rigidez local con corrección por shear locking
Coeficientes:
$$a = \frac{12\,EI}{L^3(1+\Phi)}, \quad b = \frac{6\,EI}{L^2(1+\Phi)}, \quad c_m = \frac{(4+\Phi)\,EI}{L(1+\Phi)}, \quad d = \frac{(2-\Phi)\,EI}{L(1+\Phi)}.$$

$$\mathbf K_{\text{local}} = \begin{bmatrix}
 EA/L & 0 & 0 & -EA/L & 0 & 0 \\
 0 & a & b & 0 & -a & b \\
 0 & b & c_m & 0 & -b & d \\
-EA/L & 0 & 0 & EA/L & 0 & 0 \\
 0 & -a & -b & 0 & a & -b \\
 0 & b & d & 0 & -b & c_m
\end{bmatrix}.$$

Para $\Phi \to 0$: $a \to 12EI/L^3$, $b \to 6EI/L^2$, $c_m \to 4EI/L$, $d \to 2EI/L$ — se recupera Euler-Bernoulli.

### 10. Fuerzas internas
$$\mathbf F_{\text{int,local}} = \mathbf K_{\text{local}}\,\mathbf u_{\text{local}}.$$

La componente axial se corrige con el esfuerzo axial del material: $F_1 = -\sigma A$, $F_4 = +\sigma A$ (permite materiales no-lineales en la rama axial).

### 11. Ensamblaje global
$$\mathbf K_{\text{global}} = \mathbf T^\top\,\mathbf K_{\text{local}}\,\mathbf T, \qquad \mathbf F_{\text{int}} = \mathbf T^\top\,\mathbf F_{\text{int,local}}.$$

### 12. Cuadratura
No aplica (forma cerrada).

---

## Contrato de implementación

```yaml
name: Frame2DTimoshenko
kind: element
status: validated

interface:
  dof_names: [ux, uy, rz]
  n_nodes: 2
  strain_dim: 1
  n_integration_points: 1

parameters:
  - { name: A, type: float, required: true, desc: "Área de la sección transversal" }
  - { name: I, type: float, required: true, desc: "Momento de inercia respecto al eje perpendicular al plano" }
  - { name: As, type: float, required: true, desc: "Área efectiva de cortante" }
  - { name: nu, type: float, required: false, default: 0.3, desc: "Coeficiente de Poisson (si el material no lo expone)" }

material_contract:
  signature: "compute_state(ε, state) -> (σ, E_tangent, state')"
  strain_kind: "axial escalar (deformación media ε = (u_2 − u_1)/L en sistema local)"
  nonlinearity_model: "E_tangent escala toda la matriz local (axial, flexión y cortante por igual)"
  poisson_source: "si material.nu existe, se usa; en su defecto el parámetro `nu` del elemento con default 0.3 y advertencia por consola"

conventions:
  sign: "σ > 0 ⇔ ε > 0 (elongación) en el eje axial"
  rotation_sign: "r_z positivo antihorario"
  node_orientation: "eje local del nodo 1 al nodo 2"
  configuration: "inicial fija — L, c, s, T, Φ (en forma funcional) se calculan sin actualizar"

validity:
  - "vigas peraltadas o cortas (L/h ≲ 10) donde el cortante no es despreciable"
  - pequeños desplazamientos y pequeñas rotaciones
  - "|ε_axial| ≲ 1e-2"

out_of_scope:
  - grandes rotaciones o desplazamientos (no hay variante corotacional en el catálogo actual)
  - plasticidad distribuida en la sección
  - pandeo
  - "fuerzas de cuerpo distribuidas: el usuario reparte carga equivalente a los nodos"

acceptance:
  - name: convergencia a Euler-Bernoulli para viga esbelta
    setup: "voladizo con carga transversal P, viga esbelta (L/h ≳ 100) ⇒ Φ → 0"
    expect: "v(L) ≈ PL³/(3EI) (solución Euler-Bernoulli) con tolerancia relativa ≲ 1e-3"
    tol_rel: 1.0e-3

  - name: respuesta axial pura desacoplada
    setup: "voladizo con carga axial F"
    expect: "u_x(L) = F·L/(E·A); momentos y cortantes nulos"
    tol_rel: 1.0e-10

  - name: simetría de K
    setup: "viga oblicua"
    expect: "K_global = K_globalᵀ"
    tol_abs: 1.0e-12

references:
  - "Reddy J.N., An Introduction to the Finite Element Method, cap. 5 (vigas de Timoshenko)"
  - "Bathe K.-J., Finite Element Procedures, §5.4 (formulación con factor Φ)"
```

---

## Implementación

- **Archivo**: [fenix/elements/frame/timoshenko.py](../../fenix/elements/frame/timoshenko.py) — submódulo del paquete `fenix/elements/frame/` que aloja también `Frame2DEuler` y `Frame2DEulerCorot`.
- **Clase**: `Frame2DTimoshenko` — hereda directamente de `Element`, sin herencia con las otras vigas 2D. Comparte `build_geometry_2d` con `Frame2DEuler` en [fenix/elements/frame/_shared.py](../../fenix/elements/frame/_shared.py); además, los helpers libres de carga de cuerpo, masa consistente y traducción a `ElementForces` también viven en `_shared.py`.
- **Tests**: [tests/test_frame.py](../../tests/test_frame.py) · `TestFrame2DTimoshenkoAcceptance`:
  - `test_acceptance_convergencia_euler_en_viga_esbelta` (criterio 1) — viga con $L/h$ muy grande, verifica que la flecha tiende a $PL^3/(3EI)$ con tolerancia relativa $10^{-3}$.
  - `test_acceptance_respuesta_axial_pura` (criterio 2) — carga axial pura, verifica $u_x = FL/(EA)$.
  - `test_acceptance_simetria_K` (criterio 3) — viga oblicua, verifica $\mathbf K = \mathbf K^\top$.
  - `test_registro_en_registry` — autodiscover.
- **Notas de traducción**:
  - Para no colisionar con la variable cinemática $c$ (coseno director) dentro de `compute_element_state`, el coeficiente $c_m$ de la matriz local se llama `c_coef` en el código y `d_coef` su análogo (antes eran `c` y `d` inline, ambiguos).
  - El Poisson $\nu$ se toma de `material.nu` si el material lo expone; en su defecto, del parámetro `nu` del elemento con default 0.3 y una advertencia por consola la primera vez que se construye con ese default — atajo a errores silenciosos de usuarios que olvidan especificarlo.
  - La matriz $\mathbf T$ se calcula una sola vez en `__init__` y se guarda como atributo (régimen geométricamente lineal).
  - El factor $\Phi$ se recalcula en cada llamada porque depende de $E_t$ devuelto por el material (en régimen no-lineal $\Phi$ cambia con el estado).

---

## Diálogo

- **2026-04-21** · Elemento movido a archivo propio `fenix/elements/frame.py` y desacoplado del helper `_frame_geometry`. Con este movimiento, el helper compartido se elimina completamente del repositorio: `Frame2DEuler` y `Frame2DTimoshenko` replican la construcción de $\mathbf T$ como método estático `_build_geometry`, cada una en su clase. La duplicación (~15 líneas) es el precio aceptado de la independencia mutua.
- **2026-04-21** · Test del criterio 1: se eligió $L/h$ grande en lugar de reproducir un caso específico de libro porque la convergencia Timoshenko → Euler es el comportamiento **esperado por construcción** (factor $\Phi \to 0$). Verificarlo directamente con el solver completo valida simultáneamente la cinemática de Timoshenko, el tratamiento del shear locking y la integración con ensamblador + solver.
- **2026-05-13** · `frame.py` se parte en paquete `fenix/elements/frame/`. La duplicación de `_build_geometry` ya no se justifica: ambas vigas comparten `build_geometry_2d` en `_shared.py`, helper interno al paquete (no flotante). Las clases siguen sin herencia entre sí.
