# ELEMENTO MARCO/VIGA 2D TIMOSHENKO

> Orden de trabajo. El usuario escribe **especificación física**, **formulación numérica** y **contrato**; la IA rellena **implementación** y responde en **diálogo**.

---

## Especificación física

### 0. Descripción general
Elemento 1D inmerso en el plano 2D que modela una **viga** según la teoría de **Timoshenko**. Dos nodos rígidamente conectados, 3 DOFs por nodo ($u_x, u_y, r_z$). Transmite **fuerza axial**, **cortante transverso** y **momento flector**. Incluye explícitamente la **deformación por cortante transversal**. Régimen de **linealidad geométrica** (pequeños desplazamientos y rotaciones).

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
- Axial: $N = EA\,\varepsilon_{axial}$.
- Flexión: $M = EI\,\kappa$.
- Cortante: $V = G A_s\,\gamma$, con $G = E/[2(1+\nu)]$ (módulo de corte) y $A_s$ el **área efectiva de cortante** (menor que $A$ por la distribución no uniforme del esfuerzo cortante $\sigma_{xy}$ en la sección).

**Sólo material elástico lineal** (`IS_LINEAR_ELASTIC = True`; [ADR 0021](../adr/0021-secciones-de-marco-y-seccion-por-fibras.md)). Un material inelástico se rechaza al construir. Hasta el 2026-09-29 se admitía y su $E_t$ axial escalaba toda la matriz local (flector y cortante nulos con fluencia axial, de signo cambiado con ablandamiento). La sección por fibras de este elemento, que necesita su propia cinemática con cortante, es la **etapa 2** del ADR 0021 y está pendiente; mientras tanto, la plasticidad distribuida en la sección está disponible en los marcos Euler-Bernoulli (`Frame2DEuler`, `Frame2DEulerCorot`, `Frame3D`) con `section`.

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

La componente axial se escribe con el esfuerzo axial del material, $F_1 = -\sigma A$, $F_4 = +\sigma A$; con material elástico lineal (el único admitido) coincide con la fila axial del producto.

### 11. Ensamblaje global
$$\mathbf K_{\text{global}} = \mathbf T^\top\,\mathbf K_{\text{local}}\,\mathbf T, \qquad \mathbf F_{\text{int}} = \mathbf T^\top\,\mathbf F_{\text{int,local}}.$$

### 12. Cuadratura
No aplica (forma cerrada).

### 13. Matriz de masa
**Con inercia rotatoria** de la sección: el giro $\theta$ es independiente de la pendiente (§1) y su energía cinética $\tfrac12\rho I\dot\theta^2$ no se desprecia, a diferencia de los marcos Euler-Bernoulli ([Frame2DEuler](Frame2DEuler.md) §13).

- **Consistente**: axial $\tfrac{\rho A L}{6}\left[\begin{smallmatrix}2&1\\1&2\end{smallmatrix}\right]$, flexión hermítica $\tfrac{\rho A L}{420}\mathbf H$ en $(v_1, \theta_1, v_2, \theta_2)$ (la de Euler-Bernoulli; la interpolación dependiente de $\Phi$ añadiría términos de orden $\Phi$ que se desprecian) e inercia rotatoria $\tfrac{\rho I L}{6}\left[\begin{smallmatrix}2&1\\1&2\end{smallmatrix}\right]$ en $(\theta_1, \theta_2)$.
- **Concentrada**: $\rho A L/2$ traslacional y $\rho I L/2$ rotacional por nodo.

Se conserva hasta la etapa 2 de la sección por fibras ([ADR 0021](../adr/0021-secciones-de-marco-y-seccion-por-fibras.md)).

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

admission:                # criterio de admisión al programa principal (Reglas.md §1)
  textbooks:
    - zienkiewicz_taylor: "vol. 2 (6.ª ed.), cap. 10: barras de Timoshenko (formulación base, material elástico)"
  general_codes:
    - opensees: "ElasticTimoshenkoBeam (2D)"

parameters:
  - { name: A, type: float, required: true, desc: "Área de la sección transversal" }
  - { name: I, type: float, required: true, desc: "Momento de inercia respecto al eje perpendicular al plano" }
  - { name: As, type: float, required: true, desc: "Área efectiva de cortante" }
  - { name: nu, type: float, required: false, default: 0.3, desc: "Coeficiente de Poisson (si el material no lo expone)" }

material_contract:
  signature: "compute_state(ε, state) -> (σ, E_tangent, state')"
  strain_kind: "axial escalar (deformación media ε = (u_2 − u_1)/L en sistema local)"
  nonlinearity_model: "sólo materiales elásticos lineales (IS_LINEAR_ELASTIC = True; un material inelástico se rechaza al construir). Sección por fibras con cortante: etapa 2 del ADR 0021, pendiente"
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
  - material inelástico y plasticidad distribuida en la sección (sección por fibras de Timoshenko: etapa 2 del ADR 0021; hoy se rechaza al construir)
  - pandeo
  - "cargas puntuales en el vano (una fuerza de cuerpo variable b(x) sí, desde el 2026-10-27, ADR 0022)"

acceptance:
  - name: soluciones manufacturadas (ADR 0022)                   # 2026-10-27
    setup: "línea recta inclinada; solución exacta impuesta en un extremo y sus fuerzas internas como cargas nodales en el otro; carga de línea como fuerza de cuerpo b(x); u* axial y, en flexión, una función ψ de la que salen θ* = ψ' y v* = ψ − (EI/(G·A_s))·ψ'' (sin momento distribuido); viga gruesa, EI/(G·A_s) = 0.47 con longitud 2; n = 4, 8, 16, 32"
    expect: "con campos polinómicos (carga de grado ≤ 4), error nodal de redondeo con 1 y 3 elementos (exactitud nodal, Tong 1969), con la carga repartida con las funciones de Friedman y Kosmatka; con campos suaves, orden ≥ teórico − 0.1 por componente: u 2, θ 3, v 2 en L²; ε 1, κ 2, γ 1 en energía; con A_s enorme (Φ ≪ 1), v orden 4 como Navier-Bernoulli"
    tol_abs: 0.1
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

  - name: guarda de material elástico lineal            # añadido 2026-09-29 (ADR 0021)
    setup: "material Elastoplastic1D"
    expect: "ValueError al construir, con un mensaje que dice que el elemento todavía no admite sección por fibras y nombra los marcos que sí la admiten"

references:
  - "Reddy J.N., An Introduction to the Finite Element Method, cap. 5 (vigas de Timoshenko)"
  - "Bathe K.-J., Finite Element Procedures, §5.4 (formulación con factor Φ)"
```

---

## Implementación

- **Archivo**: [solidum/elements/frame/timoshenko.py](../../solidum/elements/frame/timoshenko.py) — submódulo del paquete `solidum/elements/frame/` que aloja también `Frame2DEuler` y `Frame2DEulerCorot`.
- **Clase**: `Frame2DTimoshenko(FrameSectionMixin, Element)`, sin herencia con las otras vigas 2D. Del mixin ([solidum/elements/frame/_section.py](../../solidum/elements/frame/_section.py), ADR 0021) sólo usa la guarda de material elástico lineal: no declara `section` en `REFERENCE_KWARGS` (sigue siendo `{"material": MaterialRegistry}`), así que su mensaje de rechazo dice que todavía no admite sección por fibras y nombra los marcos que sí. Comparte `build_geometry_2d` con `Frame2DEuler` en [solidum/elements/frame/_shared.py](../../solidum/elements/frame/_shared.py); además, los helpers libres de carga de cuerpo, masa consistente y traducción a `ElementForces` también viven en `_shared.py`.
- **Tests**: [`tests/test_line_elements.py`](../../tests/test_line_elements.py) (campo interior, energía igual a la de la rigidez, carga como campo) y [`tests/validation/test_mms_frames.py`](../../tests/validation/test_mms_frames.py) (soluciones manufacturadas); [tests/test_frame.py](../../tests/test_frame.py) · `TestFrame2DTimoshenkoAcceptance`:
  - `test_acceptance_convergencia_euler_en_viga_esbelta` (criterio 1) — viga con $L/h$ muy grande, verifica que la flecha tiende a $PL^3/(3EI)$ con tolerancia relativa $10^{-3}$.
  - `test_acceptance_respuesta_axial_pura` (criterio 2) — carga axial pura, verifica $u_x = FL/(EA)$.
  - `test_acceptance_simetria_K` (criterio 3) — viga oblicua, verifica $\mathbf K = \mathbf K^\top$.
  - `test_registro_en_registry` — autodiscover.
  - Criterio 4 (guarda): [tests/test_fiber_section.py](../../tests/test_fiber_section.py) · `TestGuardas.test_timoshenko_solo_elastico`.
- **Notas de traducción**:
  - Para no colisionar con la variable cinemática $c$ (coseno director) dentro de `compute_element_state`, el coeficiente $c_m$ de la matriz local se llama `c_coef` en el código y `d_coef` su análogo (antes eran `c` y `d` inline, ambiguos).
  - El Poisson $\nu$ se toma de `material.nu` si el material lo expone; en su defecto, del parámetro `nu` del elemento con default 0.3 y una advertencia por consola la primera vez que se construye con ese default — atajo a errores silenciosos de usuarios que olvidan especificarlo.
  - La matriz $\mathbf T$ se calcula una sola vez en `__init__` y se guarda como atributo (régimen geométricamente lineal).
  - El factor $\Phi = 24(1+\nu)\,I/(A_s L^2)$ no depende del módulo del material (auditoría 2026-09-22) y, con material elástico lineal, tampoco la matriz local cambia con el estado.
  - Masa (§13): las funciones de `_shared.py` compartidas con los marcos Euler-Bernoulli, con el argumento `rotary_inertia` $= \rho I$ que sólo pasa este elemento. Verificada en [tests/test_frame_mass.py](../../tests/test_frame_mass.py) · `TestTimoshenkoSinCambio` y en la energía cinética del giro rígido (`test_giro_rigido_2d_consistente`, que suma $\tfrac12\rho I L\omega^2$).

---

## Diálogo

- **2026-04-21** · Elemento movido a archivo propio `solidum/elements/frame.py` y desacoplado del helper `_frame_geometry`. Con este movimiento, el helper compartido se elimina completamente del repositorio: `Frame2DEuler` y `Frame2DTimoshenko` replican la construcción de $\mathbf T$ como método estático `_build_geometry`, cada una en su clase. La duplicación (~15 líneas) es el precio aceptado de la independencia mutua.
- **2026-04-21** · Test del criterio 1: se eligió $L/h$ grande en lugar de reproducir un caso específico de libro porque la convergencia Timoshenko → Euler es el comportamiento **esperado por construcción** (factor $\Phi \to 0$). Verificarlo directamente con el solver completo valida simultáneamente la cinemática de Timoshenko, el tratamiento del shear locking y la integración con ensamblador + solver.
- **2026-05-13** · `frame.py` se parte en paquete `solidum/elements/frame/`. La duplicación de `_build_geometry` ya no se justifica: ambas vigas comparten `build_geometry_2d` en `_shared.py`, helper interno al paquete (no flotante). Las clases siguen sin herencia entre sí.
- **2026-09-22** · Auditoría global: `internal_forces(U, equivalent_load=None)`: cuando el `Assembler` ensambló peso propio o fuerza de cuerpo, `solidum.run` pasa el vector nodal equivalente del elemento (escalado por λ_final) y las fuerzas internas de extremo son `F_int − f_eq` en ejes locales. Antes se devolvía `K·u`, que con carga distribuida son fuerzas nodales, no fuerzas internas: un voladizo bajo peso propio daba V_i = qL/2 y M_i = −5qL²/12 en vez de qL y −qL²/2, contradiciendo a las reacciones.
- **2026-09-22** · Auditoría global: `Φ = 24(1+ν)·I/(A_s·L²)`, independiente de E_t: con `E_t = 0` (plasticidad perfecta plastificada) la forma `12·E_t·I/(G·A_s·L²)` dividía 0/0.
- **2026-09-29** · Sólo material elástico lineal hasta la etapa 2 del [ADR 0021](../adr/0021-secciones-de-marco-y-seccion-por-fibras.md). La auditoría de pureza verificó que el modelo no lineal anterior —la $E_t$ de la deformación axial escalando toda la matriz local, flexión y cortante incluidos— daba fuerzas internas incorrectas en silencio (nulas con fluencia axial, de signo cambiado con ablandamiento). Desaparece: un material con `IS_LINEAR_ELASTIC = False` se rechaza al construir, con un mensaje que remite a los marcos Euler-Bernoulli con sección por fibras. La sección por fibras de este elemento necesita su propia cinemática con cortante (qué hacer con $\gamma$ en las fibras, integración del cortante, tratamiento del bloqueo) y es una decisión física pendiente del usuario (etapa 2). Contrato: criterio de aceptación 4 (guarda); `nonlinearity_model` y `out_of_scope` actualizados.
- **2026-09-29** · Aclaración, sin cambio de resultados: los marcos Euler-Bernoulli pierden la inercia rotatoria de su masa (cambio aprobado por el usuario, [Frame2DEuler](Frame2DEuler.md) *Diálogo*); este elemento la **conserva** (§13), porque en Timoshenko el giro es independiente, hasta la etapa 2 de las fibras. La matriz es idéntica bit a bit a la anterior; las funciones compartidas de `_shared.py` la reciben ahora como argumento explícito.
- **2026-10-27** · Cargas variables y verificación por soluciones manufacturadas (ADR 0022, fase 4), **con la decisión física del usuario**: la carga variable se reparte con las funciones de interdependencia del propio elemento (Friedman y Kosmatka 1993), que dan exactamente su rigidez con la corrección `Φ` (comprobado con SymPy) y conservan la exactitud nodal; con las de Hermite quedaría un error nodal de 12 %, 2.8 % y 0.7 % con 2, 4 y 8 elementos. Con `b` constante coinciden (`qL/2`, `±qL²/12`) y no cambia ningún resultado. Campo interior (`line_field`: `u`, `v`, `θz`, `ε`, `κz`, `γy` constante) con `section_rigidities` (`EA`, `EI`, `G·A_s`). Órdenes medidos: u 2, θ 3, κ 2, γ 1 y v 2: al refinar `Φ ∝ 1/h²` crece y la `v` del elemento aproxima `v''` con `θ'`, que difieren en `γ'`; con `Φ ≪ 1`, v orden 4 como Navier-Bernoulli. Sin momento distribuido (decisión del usuario), la solución manufacturada sale de una función `ψ`.
