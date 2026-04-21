# Catálogo de elementos

> Referencia rápida de los elementos implementados. Una entrada por elemento. Para detalles físicos/numéricos → código fuente.
>
> **Convenciones**: `STRAIN_DIM` = dimensión Voigt esperada del material asociado (1 = axial escalar, 3 = 2D `[ε_xx, ε_yy, γ_xy]`, 6 = 3D). DOFs por nodo = `DOF_NAMES`.

---

## Truss2D — barra axial 2D

Elemento sólido 1D de primer orden, inmerso en el plano. Dos nodos; transmite exclusivamente esfuerzo axial.

### Cinemática
- Interpolación lineal del desplazamiento entre los dos nodos.
- Deformación axial uniforme en el elemento.

### Formulación
Cosenos directores del eje: `c = (x₂−x₁)/L`, `s = (y₂−y₁)/L`.

```
ε     = B · u_e,          B = (1/L) · [−c, −s, c, s]
K_e   = (E·A / L) · dᵀd,  d = [−c, −s, c, s]
F_int = σ·A · d
```

### Convenciones
- Convención de signos: `σ > 0 ⇔ ε > 0 (elongación)`.
- `STRAIN_DIM = 1`: la deformación que recibe el material es un escalar (no Voigt).
- La orientación de los nodos no afecta el resultado (`K_e` y `F_int` son invariantes bajo su permutación).

### Parámetros
- `A` — área de la sección transversal.

### Régimen de validez
- Pequeñas deformaciones (`|ε| ≲ 10⁻²`) y pequeños desplazamientos.
- Cargas exclusivamente axiales. Si la carga transversal sobre la barra no es despreciable → `Frame2DEuler` / `Frame2DTimoshenko`.

### Validación
- Tests: [tests/test_truss.py](tests/test_truss.py) · `TestTruss2D`.
- Archivo fuente: [fenix/elements/truss.py](fenix/elements/truss.py) · clase `Truss2D`.
- Spec: [docs/specs/Truss2D.md](specs/Truss2D.md).
- Referencia: Bathe, *Finite Element Procedures*, §4.2.1.

---

## Truss2DCorot — armadura 2D corotacional

Barra axial 2D en régimen de **grandes desplazamientos y rotaciones con pequeña deformación** (Updated Lagrangian). Hereda de `Truss2D`; comparte DOFs, parámetros y contrato con el material.

### Cinemática
- Longitud y cosenos directores se recalculan en configuración **corriente** en cada evaluación.
- Deformación ingenieril corotacional: `ε = (l − L₀)/L₀`.

### Formulación
```
d = [−c_θ, −s_θ, c_θ, s_θ]      (dirección corriente del eje)
n = [−s_θ,  c_θ, s_θ, −c_θ]     (perpendicular al eje)
K_M   = (E·A / L₀) · d·dᵀ
K_G   = (N / l) · n·nᵀ           (rigidez geométrica)
K_T   = K_M + K_G
F_int = N · d
```

### Régimen de validez
- `|ε| ≲ 10⁻²` (pequeña deformación axial).
- Desplazamientos y rotaciones de cualquier magnitud.
- Cargas exclusivamente axiales.
- **Fuera de alcance**: pandeo por bifurcación (se capta ablandamiento precrítico pero no el punto crítico; para snap-through usar arc-length).

### Validación
- Tests: [tests/test_truss.py](tests/test_truss.py) · `TestTruss2DCorot` (3 tests, uno por criterio de aceptación).
- Archivo fuente: [fenix/elements/truss.py](fenix/elements/truss.py) · clase `Truss2DCorot`.
- Spec: [docs/specs/Truss2DCorot.md](specs/Truss2DCorot.md).
- Referencias: Crisfield §3.3, Belytschko §4.5.

---

## Truss3D — barra axial 3D

Elemento sólido 1D de primer orden, inmerso en el espacio tridimensional. Dos nodos articulados; transmite exclusivamente esfuerzo axial. Régimen estrictamente de **linealidad geométrica**.

### Cinemática
- Interpolación lineal del desplazamiento entre los dos nodos.
- Deformación axial uniforme en el elemento.
- Configuración inicial fija: `L₀`, `cₓ`, `c_y`, `c_z` no se actualizan con los desplazamientos.

### Formulación
Cosenos directores del eje: `cₓ = (x₂−x₁)/L`, `c_y = (y₂−y₁)/L`, `c_z = (z₂−z₁)/L`.

```
ε     = B · u_e,          B = (1/L) · [−cₓ, −c_y, −c_z, cₓ, c_y, c_z]
K_e   = (E·A / L) · d·dᵀ,  d = [−cₓ, −c_y, −c_z, cₓ, c_y, c_z]    (6×6, rango 1)
F_int = σ·A · d
```

### Convenciones
- Convención de signos: `σ > 0 ⇔ ε > 0 (elongación)`.
- `STRAIN_DIM = 1`: la deformación que recibe el material es un escalar (no Voigt).
- La orientación de los nodos no afecta el resultado (`K_e` y `F_int` son invariantes bajo su permutación).
- Acepta nodos con 2 ó 3 coordenadas (completa con `z = 0` si la tercera falta).

### Parámetros
- `A` — área de la sección transversal.

### Régimen de validez
- Pequeñas deformaciones (`|ε| ≲ 10⁻²`), pequeños desplazamientos y pequeñas rotaciones.
- Cargas exclusivamente axiales.
- Para grandes rotaciones en 3D → `Truss3DCorot`.

### Validación
- Tests: [tests/test_truss.py](tests/test_truss.py) · `TestTruss3D`.
- Archivo fuente: [fenix/elements/truss.py](fenix/elements/truss.py) · clase `Truss3D`.
- Spec: [docs/specs/Truss3D.md](specs/Truss3D.md).
- Referencias: Bathe §4.2.1; Cook §2.3.

---

## Truss3DCorot — barra axial 3D corotacional

Barra axial 3D en régimen de **grandes desplazamientos y rotaciones con pequeña deformación** (Updated Lagrangian). Hereda de `Truss3D`; comparte DOFs, parámetros y contrato con el material.

### Cinemática
- Longitud y cosenos directores se recalculan en configuración **corriente** en cada evaluación.
- Deformación ingenieril corotacional: `ε = (l − L₀)/L₀`.
- La dirección perpendicular al eje es un **plano** (dimensión 2), no un vector único.

### Formulación
```
d = [−cₓ, −c_y, −c_z, cₓ, c_y, c_z]       (dirección corriente del eje)
ê = (cₓ, c_y, c_z),  P = I₃ − ê·êᵀ        (proyector 3×3 al plano perpendicular)

K_M   = (E·A / L₀) · d·dᵀ                  (6×6, rango 1)
K_G   = (N / l) · [[P, −P], [−P, P]]       (6×6, rango 2)
K_T   = K_M + K_G
F_int = N · d
```

### Régimen de validez
- `|ε| ≲ 10⁻²` (pequeña deformación axial).
- Desplazamientos y rotaciones de cualquier magnitud en el espacio.
- Cargas exclusivamente axiales.

### Validación
- Tests: [tests/test_truss.py](tests/test_truss.py) · `TestTruss3DCorot` (3 tests, uno por criterio de aceptación; el de rigidez geométrica verifica dos direcciones transversas independientes del plano perpendicular).
- Archivo fuente: [fenix/elements/truss.py](fenix/elements/truss.py) · clase `Truss3DCorot`.
- Spec: [docs/specs/Truss3DCorot.md](specs/Truss3DCorot.md).
- Referencias: Crisfield §3.3; Belytschko §4.5.

---

## Cable2DCorot — cable 2D corotacional

Elemento 1D inmerso en el plano que modela un cable: dos nodos, transmite solo tensión. Cinemática corotacional; la unilateralidad la aporta el material (típicamente `CableMaterial1D`).

### Cinemática
- Idéntica a una barra corotacional: longitud y cosenos directores se recalculan en configuración corriente en cada evaluación.
- Deformación ingenieril corotacional: `ε = (l − L₀)/L₀`.
- La respuesta física de cable aparece **solo** si el material devuelve `σ = 0` en compresión.

### Formulación
```
d = [−c_θ, −s_θ, c_θ, s_θ]
n = [−s_θ,  c_θ, s_θ, −c_θ]

K_M   = (E_t · A / L₀) · d·dᵀ        (0 si material destensado)
K_G   = (N / l) · n·nᵀ                (0 si N = 0)
K_T   = K_M + K_G
F_int = N · d                         (0 si N = 0)
```

### Régimen de validez
- `|ε| ≲ 10⁻²` en régimen tensado.
- Grandes desplazamientos y rotaciones en el plano.
- Cables completamente destensados aportan `K_T = 0` al sistema global: otros elementos deben garantizar la estabilidad numérica.

### Independencia del diseño
El elemento **no hereda** de `Truss2DCorot`. La cinemática corotacional se implementa dentro de la clase para que futuras modificaciones de las armaduras no afecten a los cables, y viceversa.

### Validación
- Tests: [tests/test_cable_elements.py](tests/test_cable_elements.py) · `TestCable2DCorot` (4 tests de aceptación: tensado, destensado, rotación rígida, cruce por cero).
- Archivo fuente: [fenix/elements/cable.py](fenix/elements/cable.py) · clase `Cable2DCorot`.
- Spec: [docs/specs/Cable2DCorot.md](specs/Cable2DCorot.md).
- Referencias: Crisfield §3.3; Irvine §1.2.

---

## Cable3DCorot — cable 3D corotacional

Elemento 1D inmerso en el espacio que modela un cable: dos nodos, transmite solo tensión, puede rotar libremente en el espacio. Cinemática corotacional 3D; unilateralidad aportada por el material (típicamente `CableMaterial1D`).

### Cinemática
- Idéntica en filosofía a una barra corotacional 3D: longitud y cosenos directores se recalculan en configuración corriente.
- Deformación ingenieril corotacional: `ε = (l − L₀)/L₀`.
- La dirección perpendicular al eje es un **plano** (dimensión 2).

### Formulación
```
d = [−cₓ, −c_y, −c_z, cₓ, c_y, c_z]
ê = (cₓ, c_y, c_z),  P = I₃ − ê·êᵀ      (proyector 3×3)

K_M   = (E_t · A / L₀) · d·dᵀ            (0 si material destensado)
K_G   = (N / l) · [[P, −P], [−P, P]]     (0 si N = 0)
K_T   = K_M + K_G
F_int = N · d                            (0 si N = 0)
```

### Régimen de validez
- `|ε| ≲ 10⁻²` en régimen tensado.
- Grandes desplazamientos y rotaciones en el espacio.
- Cables completamente destensados aportan `K_T = 0` al sistema global.

### Independencia del diseño
No hereda de `Cable2DCorot` ni de `Truss3DCorot`. La maquinaria cinemática 3D se implementa íntegra dentro de la clase.

### Validación
- Tests: [tests/test_cable_elements.py](tests/test_cable_elements.py) · `TestCable3DCorot` (4 tests de aceptación).
- Archivo fuente: [fenix/elements/cable.py](fenix/elements/cable.py) · clase `Cable3DCorot`.
- Spec: [docs/specs/Cable3DCorot.md](specs/Cable3DCorot.md).
- Referencias: Crisfield §3.3; Belytschko §4.5; Irvine §1.2.

---

## Frame2DEuler — marco/viga 2D Euler-Bernoulli

Viga esbelta 2D con hipótesis de Euler-Bernoulli: secciones planas perpendiculares al eje neutro deformado. Dos nodos rígidamente conectados; transmite axial + cortante + momento flector.

### Cinemática
- Interpolación lineal del desplazamiento axial, Hermite cúbica del desplazamiento transverso.
- Pequeños desplazamientos y rotaciones (régimen geométricamente lineal).
- Configuración inicial fija: `L₀, c, s, T` se calculan una vez y no se actualizan.

### Formulación
Matriz de transformación global→local $\mathbf T$ ($6\times 6$) con cosenos directores del eje. Matriz local analítica:
```
K_local = [[ EA/L,       0,       0, -EA/L,       0,       0 ],
           [    0, 12EI/L³,  6EI/L²,     0,-12EI/L³,  6EI/L²],
           [    0,  6EI/L²,   4EI/L,     0, -6EI/L²,   2EI/L],
           [-EA/L,       0,       0,  EA/L,       0,       0 ],
           [    0,-12EI/L³, -6EI/L²,     0, 12EI/L³, -6EI/L²],
           [    0,  6EI/L²,   2EI/L,     0, -6EI/L²,   4EI/L]]
K_global = Tᵀ · K_local · T
F_int    = Tᵀ · F_int_local    (componente axial corregida con σ·A del material)
```

### Régimen de validez
- Vigas esbeltas (`L/h ≳ 10`); peraltadas → `Frame2DTimoshenko`.
- Pequeños desplazamientos, pequeñas rotaciones.
- No captura rótulas plásticas distribuidas: `E_tangent` escala toda la matriz por igual.

### Independencia del diseño
Archivo propio. No hereda de ningún otro elemento y no comparte utilidades con otros módulos — la matriz de transformación se construye dentro de la clase.

### Validación
- Tests: [tests/test_frame.py](tests/test_frame.py) · `TestFrame2DEulerAcceptance` (3 tests de aceptación + registro, incluye flecha analítica del voladizo $PL^3/(3EI)$).
- Archivo fuente: [fenix/elements/frame.py](fenix/elements/frame.py) · clase `Frame2DEuler`.
- Spec: [docs/specs/Frame2DEuler.md](specs/Frame2DEuler.md).
- Referencias: Bathe cap. 5; Cook §2.7.

---

## Frame2DTimoshenko — marco/viga 2D Timoshenko

Viga 2D con deformación por cortante transversal explícita. Dos nodos rígidamente conectados; transmite axial + cortante + momento flector. Apropiado para vigas peraltadas o cortas (`L/h ≲ 10`).

### Cinemática
- Secciones planas, **no** perpendiculares al eje neutro deformado (rotación independiente de la pendiente de la elástica).
- Distorsión angular $\gamma = dv/ds - \theta$ tratada como campo independiente.
- Configuración inicial fija (régimen geométricamente lineal).

### Formulación
Factor de corrección contra shear locking:
```
Φ = 12 · E · I / (G · A_s · L²),      G = E / [2(1 + ν)]
```
Matriz local con coeficientes ajustados:
```
a = 12·EI / [L³(1+Φ)],   b = 6·EI / [L²(1+Φ)]
c = (4+Φ)·EI / [L(1+Φ)], d = (2-Φ)·EI / [L(1+Φ)]

K_local = [[ EA/L,  0,  0, -EA/L,  0,  0],
           [    0,  a,  b,     0, -a,  b],
           [    0,  b,  c,     0, -b,  d],
           [-EA/L,  0,  0,  EA/L,  0,  0],
           [    0, -a, -b,     0,  a, -b],
           [    0,  b,  d,     0, -b,  c]]
```
Para `Φ → 0` (viga esbelta) se recupera la matriz Euler-Bernoulli.

### Parámetros
- `A`, `I`, `As` (área efectiva de cortante).
- `nu` (opcional): se toma de `material.nu` si está disponible; en su defecto, del parámetro del elemento con default 0.3 y aviso por consola.

### Régimen de validez
- Vigas gruesas/peraltadas (`L/h ≲ 10`); para esbeltas → `Frame2DEuler` (más simple y sin shear locking).
- Pequeños desplazamientos, pequeñas rotaciones.
- No captura plasticidad distribuida: `E_tangent` escala toda la matriz.

### Independencia del diseño
Archivo propio. No hereda de ningún otro elemento. La matriz de transformación se construye dentro de la clase (`_build_geometry` como método estático).

### Validación
- Tests: [tests/test_frame.py](tests/test_frame.py) · `TestFrame2DTimoshenkoAcceptance` (convergencia a Euler en viga esbelta, axial puro, simetría).
- Archivo fuente: [fenix/elements/frame.py](fenix/elements/frame.py) · clase `Frame2DTimoshenko`.
- Spec: [docs/specs/Frame2DTimoshenko.md](specs/Frame2DTimoshenko.md).
- Referencias: Reddy cap. 5; Bathe §5.4.

---

## Frame2DEulerCorot — viga 2D Euler-Bernoulli corotacional

Viga 2D esbelta Euler-Bernoulli en formulación **corotacional** (Updated Lagrangian). Grandes desplazamientos y grandes rotaciones rígidas del elemento; deformaciones axiales pequeñas y rotaciones nodales deformacionales moderadas.

### Cinemática corotacional
- Rotación rígida $\alpha_e$ del eje separada de las rotaciones deformacionales de las secciones $\bar\theta_1, \bar\theta_2$.
- DOFs deformacionales locales: $[u_l = l - L_0, \bar\theta_1, \bar\theta_2]$.
- Todo en configuración corriente: $c, s, l, \alpha$ se recalculan por evaluación.

### Formulación
```
Rigidez local (3×3) en DOFs deformacionales:
  K_ll = diag-block([EA/L₀], [[4EI/L₀, 2EI/L₀], [2EI/L₀, 4EI/L₀]])

Matriz B (3×6) en configuración corriente acopla DOFs locales y globales.

Rigidez tangente:
  K_T = B ᵀ · K_ll · B + K_σ
  K_σ = (N/l)·z·zᵀ − ((M₁+M₂)/l²)·(r·zᵀ + z·rᵀ)

con r = [−c,−s,0,c,s,0], z = [−s,c,0,s,−c,0].
```

### Régimen de validez
- `|ε_axial| ≲ 10⁻²`.
- Rotaciones rígidas del elemento ilimitadas entre commits (hasta $|\alpha_e| < \pi$ por paso).
- Rotaciones deformacionales de las secciones moderadas (`|θ̄| ≲ 30°`).
- Rotaciones continuas > π requerirían tracking de vueltas (fuera de alcance).

### Dominios que abre
Pandeo por flexión de columnas esbeltas, post-pandeo, snap-through de arcos, brazos flexibles, cables con rigidez a flexión. Problemas que la viga lineal `Frame2DEuler` no puede atacar.

### Independencia del diseño
Archivo compartido con `Frame2DEuler` y `Frame2DTimoshenko` por familia temática (vigas 2D), pero **sin herencia**: cada clase replica su propia cinemática internamente.

### Validación
- Tests: [tests/test_frame.py](tests/test_frame.py) · `TestFrame2DEulerCorotAcceptance` (4 criterios físicos + **chequeo por diferencias finitas de $\mathbf K_T$ contra $\mathbf F_{\text{int}}$** + registro).
- Archivo fuente: [fenix/elements/frame.py](fenix/elements/frame.py) · clase `Frame2DEulerCorot`.
- Spec: [docs/specs/Frame2DEulerCorot.md](specs/Frame2DEulerCorot.md).
- Referencias: Crisfield §7.3; Belytschko §4.11; Wriggers cap. 4.

---

## Frame3D — marco/viga 3D Euler-Bernoulli

Viga 3D esbelta con 6 DOFs por nodo (`ux, uy, uz, rx, ry, rz`). Transmite axial + cortante en dos planos + flexión en dos planos + torsión. Régimen geométricamente lineal.

### Cinemática
- Secciones planas perpendiculares al eje neutro deformado (Euler-Bernoulli).
- Torsión de Saint-Venant pura (sin alabeo).
- Flexiones en ejes principales desacopladas.
- Configuración inicial fija (no hay variante corotacional 3D por ahora).

### Ejes locales y orientación de la sección
El eje $\hat{\mathbf x}$ local va del nodo 1 al nodo 2. Los ejes $\hat{\mathbf y}, \hat{\mathbf z}$ de la sección se definen proyectando un **vector de referencia** (`ref_vector`) sobre el plano perpendicular al eje. Default: `[0, 0, 1]` (eje z global); fallback a `[1, 0, 0]` con advertencia si la barra es casi vertical.

### Formulación
Matriz de transformación $\mathbf T$ 12×12 en bloques de $\boldsymbol{\lambda}$ (3×3). Matriz local 12×12 desacoplada:
```
Axial      (2×2) : EA/L
Torsión    (2×2) : GJ/L       con G = E/[2(1+ν)]
Flexión xy (4×4) : 12EIz/L³, 6EIz/L², 4EIz/L, 2EIz/L
Flexión xz (4×4) : 12EIy/L³, 6EIy/L², 4EIy/L, 2EIy/L  (signos invertidos)

K_global = Tᵀ K_local T
```

### Parámetros
- `A`, `Iy`, `Iz` (área y momentos de inercia respecto a ejes locales y, z).
- `J` (constante torsional de Saint-Venant).
- `nu` (opcional, del material si lo expone, sino del elemento con default 0.3).
- `ref_vector` (opcional, default `[0, 0, 1]`).

### Régimen de validez
- Vigas esbeltas (`L/h ≳ 10`).
- Pequeños desplazamientos, pequeñas rotaciones.
- Sin alabeo (secciones macizas o cerradas).
- Para grandes rotaciones en 3D: pendiente `Frame3DCorot`.

### Independencia del diseño
Archivo propio. No hereda ni comparte helpers con `Frame2DEuler`, `Frame2DTimoshenko`, ni con los trusses 3D.

### Validación
- Tests: [tests/test_frame3d.py](tests/test_frame3d.py) · `TestFrame3DAcceptance` (5 criterios físicos + validación de `ref_vector` + registro).
- Archivo fuente: [fenix/elements/frame3d.py](fenix/elements/frame3d.py) · clase `Frame3D`.
- Spec: [docs/specs/Frame3D.md](specs/Frame3D.md).
- Referencias: Przemieniecki Tabla 11.3; Cook §2.8; Bathe cap. 5.

---

## Quad4 — cuadrilátero bilineal 2D

- **Propósito**: continuo 2D isoparamétrico; problema mecánico plano.
- **DOFs por nodo**: `['ux', 'uy']` · 4 nodos (antihorario) · `STRAIN_DIM = 3`.
- **Cinemática**: matriz `B(ξ, η)` calculada por mapeo isoparamétrico estándar; deformación Voigt `[ε_xx, ε_yy, γ_xy]`.
- **Integración**: por defecto Gauss 2×2 (4 puntos). Configurable vía `quadrature` desde `QuadratureRegistry`. **Aviso**: integración reducida 1×1 → riesgo de hourglassing.
- **Parámetros**: `thickness` (espesor para estado plano), `quadrature` (opcional).
- **Implementación**: kernels `_compute_kinematics` y `_compute_integrands` con `@njit` (Numba) — JIT en primera llamada.
- **Limitaciones**: bloqueo volumétrico con materiales casi-incompresibles (ν → 0.5); en ese régimen usar formulación mixta (no implementada).
- **Archivo**: [fenix/elements/solid_2d.py](fenix/elements/solid_2d.py)

---

## Tri3 — triángulo lineal 2D (CST)

- **Propósito**: triángulo de deformación constante; útil para transiciones de malla.
- **DOFs por nodo**: `['ux', 'uy']` · 3 nodos · `STRAIN_DIM = 3`.
- **Cinemática**: matriz `B` constante (deformación uniforme dentro del elemento).
- **Integración**: 1 punto central, peso 0.5 (área triángulo en coordenadas naturales).
- **Parámetros**: `thickness`.
- **Limitaciones**: shear locking severo; convergencia lenta. **Preferir `Quad4`** salvo en transiciones donde Quad4 no encaja geométricamente.
- **Implementación**: kernel `_compute_kinematics_tri3` con `@njit`.
- **Archivo**: [fenix/elements/solid_2d.py](fenix/elements/solid_2d.py)

---

## Cómo añadir un elemento nuevo

1. **Spec primero** — el usuario crea `docs/specs/<Nombre>.md` a partir de `docs/specs/_template_element.md` (especificación física + formulación + contrato YAML). Sin spec, la IA no escribe código.
2. **Scaffolding** — `/fenix-new element <Nombre>` genera archivo en `fenix/elements/`, decorador `@ElementRegistry.register` y esqueleto de test.
3. **Implementación + validación** — la IA codifica contra la spec; los tests cubren los casos de `acceptance` declarados.
4. **Catálogo** — cuando la spec pasa a `status: validated`, se añade aquí una entrada breve siguiendo el formato de arriba (la spec sigue siendo la referencia detallada).
