# Catálogo de elementos

> Referencia rápida de los elementos implementados. Una entrada por elemento. Para detalles físicos/numéricos → código fuente.
>
> **Convenciones**: `STRAIN_DIM` = dimensión Voigt esperada del material asociado (1 = axial escalar, 3 = 2D `[ε_xx, ε_yy, γ_xy]`, 6 = 3D). DOFs por nodo = `DOF_NAMES`.

---

## Truss2D — barra axial 2D

Elemento sólido 1D de primer orden, inmerso en el plano. Dos nodos; transmite exclusivamente fuerza axial.

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

### Cargas distribuidas
- `compute_body_load(b)` reparte `b · A · L₀ / 2` por nodo en cada componente global (forma cerrada exacta para `b` uniforme). Útil para peso propio: `b = (0, -ρ·g)`.

### Régimen de validez
- Pequeñas deformaciones (`|ε| ≲ 10⁻²`) y pequeños desplazamientos.
- Cargas exclusivamente axiales. Si la carga transversal sobre la barra no es despreciable → `Frame2DEuler` / `Frame2DTimoshenko`.

### Validación
- Tests: [tests/test_truss.py](tests/test_truss.py) · `TestTruss2D`.
- Archivo fuente: [solidum/elements/truss.py](solidum/elements/truss.py) · clase `Truss2D`.
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

### Cargas distribuidas
- Heredado de `Truss2D`: `compute_body_load(b)` reparte `b · A · L₀ / 2` por nodo, evaluado en geometría de referencia (aproximación estándar para cargas conservadoras en formulaciones corotacionales).

### Régimen de validez
- `|ε| ≲ 10⁻²` (pequeña deformación axial).
- Desplazamientos y rotaciones de cualquier magnitud.
- Cargas exclusivamente axiales.
- **Fuera de alcance**: pandeo por bifurcación (se capta ablandamiento precrítico pero no el punto crítico; para snap-through usar arc-length).

### Validación
- Tests: [tests/test_truss.py](tests/test_truss.py) · `TestTruss2DCorot` (3 tests, uno por criterio de aceptación).
- Archivo fuente: [solidum/elements/truss.py](solidum/elements/truss.py) · clase `Truss2DCorot`.
- Spec: [docs/specs/Truss2DCorot.md](specs/Truss2DCorot.md).
- Referencias: Crisfield §3.3, Belytschko §4.5.

---

## Truss3D — barra axial 3D

Elemento sólido 1D de primer orden, inmerso en el espacio tridimensional. Dos nodos articulados; transmite exclusivamente fuerza axial. Régimen estrictamente de **linealidad geométrica**.

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

### Cargas distribuidas
- `compute_body_load(b)` reparte `b · A · L₀ / 2` por nodo en cada componente global (forma cerrada exacta). Útil para peso propio: `b = (0, 0, -ρ·g)`.

### Régimen de validez
- Pequeñas deformaciones (`|ε| ≲ 10⁻²`), pequeños desplazamientos y pequeñas rotaciones.
- Cargas exclusivamente axiales.
- Para grandes rotaciones en 3D → `Truss3DCorot`.

### Validación
- Tests: [tests/test_truss.py](tests/test_truss.py) · `TestTruss3D`.
- Archivo fuente: [solidum/elements/truss.py](solidum/elements/truss.py) · clase `Truss3D`.
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

### Cargas distribuidas
- Heredado de `Truss3D`: `compute_body_load(b)` reparte `b · A · L₀ / 2` por nodo, evaluado en geometría de referencia.

### Régimen de validez
- `|ε| ≲ 10⁻²` (pequeña deformación axial).
- Desplazamientos y rotaciones de cualquier magnitud en el espacio.
- Cargas exclusivamente axiales.

### Validación
- Tests: [tests/test_truss.py](tests/test_truss.py) · `TestTruss3DCorot` (3 tests, uno por criterio de aceptación; el de rigidez geométrica verifica dos direcciones transversas independientes del plano perpendicular).
- Archivo fuente: [solidum/elements/truss.py](solidum/elements/truss.py) · clase `Truss3DCorot`.
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

### Cargas distribuidas
- `compute_body_load(b)` reparte `b · A · L₀ / 2` por nodo en cada componente global, evaluado en geometría de referencia. Útil para peso propio.

### Régimen de validez
- `|ε| ≲ 10⁻²` en régimen tensado.
- Grandes desplazamientos y rotaciones en el plano.
- Cables completamente destensados aportan `K_T = 0` al sistema global: otros elementos deben garantizar la estabilidad numérica.

### Independencia del diseño
El elemento **no hereda** de `Truss2DCorot`. La cinemática corotacional se implementa dentro de la clase para que futuras modificaciones de las armaduras no afecten a los cables, y viceversa.

### Validación
- Tests: [tests/test_cable_elements.py](tests/test_cable_elements.py) · `TestCable2DCorot` (4 tests de aceptación: tensado, destensado, rotación rígida, cruce por cero).
- Archivo fuente: [solidum/elements/cable.py](solidum/elements/cable.py) · clase `Cable2DCorot`.
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

### Cargas distribuidas
- `compute_body_load(b)` reparte `b · A · L₀ / 2` por nodo en cada componente global, evaluado en geometría de referencia.

### Régimen de validez
- `|ε| ≲ 10⁻²` en régimen tensado.
- Grandes desplazamientos y rotaciones en el espacio.
- Cables completamente destensados aportan `K_T = 0` al sistema global.

### Independencia del diseño
No hereda de `Cable2DCorot` ni de `Truss3DCorot`. La maquinaria cinemática 3D se implementa íntegra dentro de la clase.

### Validación
- Tests: [tests/test_cable_elements.py](tests/test_cable_elements.py) · `TestCable3DCorot` (4 tests de aceptación).
- Archivo fuente: [solidum/elements/cable.py](solidum/elements/cable.py) · clase `Cable3DCorot`.
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

### Cargas distribuidas
- `compute_body_load(b)` distribuye `q = A · b` (carga por unidad de longitud) con la fórmula consistente Hermite cúbica: axial mitad-mitad por nodo, transversal `q·L/2` en fuerza y `±q·L²/12` en momento (signo positivo en nodo i, negativo en j). Transformación local↔global vía la misma `T` del elemento. Útil para peso propio.

### Régimen de validez
- Vigas esbeltas (`L/h ≳ 10`); peraltadas → `Frame2DTimoshenko`.
- Pequeños desplazamientos, pequeñas rotaciones.
- No captura **plasticidad por flexión** distribuida en la sección: `E_tangent` escala toda la matriz de rigidez por igual. La fluencia que ocurre primero en la fibra inferior bajo flexión no se modela hasta que entre `FiberSection` (deuda técnica documentada).

### Independencia del diseño
Vive en el paquete `solidum/elements/frame/` junto a `Frame2DTimoshenko` y `Frame2DEulerCorot`, sin herencia entre ellas. La construcción de longitud, cosenos directores y matriz de transformación 6×6 se delega a `build_geometry_2d` en `solidum/elements/frame/_shared.py`, compartida con `Frame2DTimoshenko` (la versión corotacional reconstruye T desde `alpha0` por su lógica propia).

### Validación
- Tests: [tests/test_frame.py](tests/test_frame.py) · `TestFrame2DEulerAcceptance` (3 tests de aceptación + registro, incluye flecha analítica del voladizo $PL^3/(3EI)$).
- Archivo fuente: [solidum/elements/frame/euler.py](solidum/elements/frame/euler.py) · clase `Frame2DEuler`.
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

### Cargas distribuidas
- `compute_body_load(b)` usa la misma fórmula consistente Hermite cúbica que `Frame2DEuler` (idéntica para carga uniforme — la diferencia entre formulaciones aparece solo con cargas no uniformes o concentradas).

### Régimen de validez
- Vigas gruesas/peraltadas (`L/h ≲ 10`); para esbeltas → `Frame2DEuler` (más simple y sin shear locking).
- Pequeños desplazamientos, pequeñas rotaciones.
- No captura plasticidad distribuida: `E_tangent` escala toda la matriz.

### Independencia del diseño
Vive en el paquete `solidum/elements/frame/`, sin herencia con las otras vigas 2D. Comparte `build_geometry_2d` con `Frame2DEuler` en `solidum/elements/frame/_shared.py` — la construcción de la matriz de transformación 6×6 ya no se duplica.

### Validación
- Tests: [tests/test_frame.py](tests/test_frame.py) · `TestFrame2DTimoshenkoAcceptance` (convergencia a Euler en viga esbelta, axial puro, simetría).
- Archivo fuente: [solidum/elements/frame/timoshenko.py](solidum/elements/frame/timoshenko.py) · clase `Frame2DTimoshenko`.
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

### Cargas distribuidas
- `compute_body_load(b)` evalúa la integral consistente en la configuración de referencia (misma fórmula Hermite cúbica que `Frame2DEuler`). Para grandes rotaciones con cargas conservadoras (gravedad) la diferencia frente a la integración sobre la geometría corriente es de segundo orden y se asume aceptable.

### Régimen de validez
- `|ε_axial| ≲ 10⁻²`.
- Rotaciones rígidas del elemento ilimitadas entre commits (hasta $|\alpha_e| < \pi$ por paso).
- Rotaciones deformacionales de las secciones moderadas (`|θ̄| ≲ 30°`).
- Rotaciones continuas > π requerirían tracking de vueltas (fuera de alcance).

### Dominios que abre
Pandeo por flexión de columnas esbeltas, post-pandeo, snap-through de arcos, brazos flexibles, cables con rigidez a flexión. Problemas que la viga lineal `Frame2DEuler` no puede atacar.

### Independencia del diseño
Vive en el paquete `solidum/elements/frame/` junto a las otras vigas 2D, **sin herencia**: la cinemática corotacional se implementa íntegra dentro de la clase (reconstruye `T` desde `alpha0` en lugar de usar `build_geometry_2d` compartido). Sí comparte con `Frame2DEuler` y `Frame2DTimoshenko` los helpers libres de `_shared.py` para carga de cuerpo (`_frame2d_consistent_body_load`), masa (`_frame2d_consistent_mass_local`) y traducción a `ElementForces` (`_frame2d_forces_from_local`).

### Validación
- Tests: [tests/test_frame.py](tests/test_frame.py) · `TestFrame2DEulerCorotAcceptance` (4 criterios físicos + **chequeo por diferencias finitas de $\mathbf K_T$ contra $\mathbf F_{\text{int}}$** + registro).
- Archivo fuente: [solidum/elements/frame/euler_corot.py](solidum/elements/frame/euler_corot.py) · clase `Frame2DEulerCorot`.
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

### Cargas distribuidas
- `compute_body_load(b)` proyecta `q = A · b` sobre los ejes locales y reparte con la fórmula Hermite cúbica en ambos planos de flexión: `q · L / 2` en fuerzas, `±q · L² / 12` en momentos (con signos `+` en nodo i y `−` en nodo j para `Mz`, y signos opuestos para `My` por la regla de la mano derecha). No genera torsión (carga uniforme sin excentricidad respecto al centro de cortadura). Útil para peso propio: `b = (0, 0, -ρ·g)`.

### Régimen de validez
- Vigas esbeltas (`L/h ≳ 10`).
- Pequeños desplazamientos, pequeñas rotaciones.
- Sin alabeo (secciones macizas o cerradas).
- Para grandes rotaciones en 3D: pendiente `Frame3DCorot`.

### Masa lumped: limitación en orientación oblicua
La masa lumped es estrictamente diagonal en ejes locales (`ρAL/2` traslacional, `ρI·L/2` rotacional con `I` específica por DOF). Tras la rotación `T^T·M·T` a globales, el bloque traslacional 3×3 por nodo permanece diagonal porque `m_t·I₃` es invariante bajo SO(3), pero el bloque rotacional 3×3 queda **lleno** porque `m_rx = ρJp·L/2`, `m_ry = ρIy·L/2`, `m_rz = ρIz·L/2` son valores distintos para una sección real. `M_global` resulta entonces **bloque-diagonal por nodo** (no estrictamente diagonal). Es la limitación estándar del lumping en frames 3D (Cook-Malkus-Plesha §11.4); aceptable para Newmark/HHT pero `CentralDifferenceSolver` rechaza con `ValueError` y orienta al usuario hacia Newmark o hacia un eje del elemento alineado con un eje global.

### Independencia del diseño
Archivo propio. No hereda ni comparte helpers con `Frame2DEuler`, `Frame2DTimoshenko`, ni con los trusses 3D.

### Validación
- Tests: [tests/test_frame3d.py](tests/test_frame3d.py) · `TestFrame3DAcceptance` (5 criterios físicos + validación de `ref_vector` + registro).
- Archivo fuente: [solidum/elements/frame3d.py](solidum/elements/frame3d.py) · clase `Frame3D`.
- Spec: [docs/specs/Frame3D.md](specs/Frame3D.md).
- Referencias: Przemieniecki Tabla 11.3; Cook §2.8; Bathe cap. 5.

---

## Quad4 — cuadrilátero bilineal 2D

- **Propósito**: continuo 2D isoparamétrico; problema mecánico plano.
- **DOFs por nodo**: `['ux', 'uy']` · 4 nodos (antihorario) · `STRAIN_DIM = 3`.
- **Cinemática**: matriz `B(ξ, η)` calculada por mapeo isoparamétrico estándar; deformación Voigt `[ε_xx, ε_yy, γ_xy]`.
- **Integración**: por defecto Gauss 2×2 (4 puntos). Configurable vía `quadrature` desde `QuadratureRegistry`. **Aviso**: integración reducida 1×1 → riesgo de hourglassing.
- **Parámetros**: `thickness` (espesor para estado plano), `quadrature` (opcional).
- **Cargas distribuidas**: `compute_body_load(b)` integra $\int \mathbf N^\top \mathbf b\, dV$ sobre el elemento; `compute_edge_traction(edge, t̄)` reparte una tracción uniforme sobre un borde (índices 0..3 según conectividad nodal). Tracción siempre en globales; sin soporte aún para tracción variable o presión normal.
- **Salida por Gauss**: `compute_gauss_state(U)` devuelve σ, ε y coordenadas (naturales y globales) en cada punto de Gauss. Habilita el suavizado nodal del exporter VTK (`Sigma_*_nodal`).
- **Implementación**: kernels `_compute_kinematics` y `_compute_integrands` con `@njit` (Numba) — JIT en primera llamada. Viven en `solidum/elements/solid_2d/_shared.py`, compartidos con Tri3.
- **Limitaciones**: bloqueo volumétrico con materiales casi-incompresibles (ν → 0.5); en ese régimen usar formulación mixta (no implementada).
- **Archivo**: [solidum/elements/solid_2d/quad4.py](../solidum/elements/solid_2d/quad4.py)

---

## Tri3 — triángulo lineal 2D (CST)

- **Propósito**: triángulo de deformación constante; útil para transiciones de malla.
- **DOFs por nodo**: `['ux', 'uy']` · 3 nodos · `STRAIN_DIM = 3`.
- **Cinemática**: matriz `B` constante (deformación uniforme dentro del elemento).
- **Integración**: 1 punto central, peso 0.5 (área triángulo en coordenadas naturales).
- **Parámetros**: `thickness`.
- **Cargas distribuidas**: `compute_body_load(b)` reparte 1/3 a cada nodo (exacto para b uniforme); `compute_edge_traction(edge, t̄)` reparte L/2 a cada uno de los 2 nodos del borde (índices 0..2). Tracción en globales.
- **Salida por Gauss**: `compute_gauss_state(U)` devuelve el único σ, ε y centroide del elemento (mismo formato que Quad4).
- **Limitaciones**: shear locking severo; convergencia lenta. **Preferir `Quad4`** salvo en transiciones donde Quad4 no encaja geométricamente.
- **Implementación**: kernel `_compute_kinematics_tri3` con `@njit` en `solidum/elements/solid_2d/_shared.py`.
- **Archivo**: [solidum/elements/solid_2d/tri3.py](../solidum/elements/solid_2d/tri3.py)

---

## Quad8 — cuadrilátero serendípito 2D de orden 2

- **Propósito**: continuo 2D cuadrático sin shear locking severo en flexión; reproduce campos cuadráticos exactamente.
- **DOFs por nodo**: `['ux', 'uy']` · 8 nodos (4 vértices antihorarios + 4 medios de borde) · `STRAIN_DIM = 3`.
- **Funciones de forma**: serendípitas (sin término ξ²η²).
- **Integración**: Gauss 3×3 (9 puntos) por defecto. Reducida 2×2 disponible vía `quadrature` (con riesgo de modos espurios).
- **Cargas**: body load por cuadratura; tracción de borde uniforme reparte 1/6, 4/6, 1/6 (vértice, medio, vértice).
- **Salida por Gauss**: `compute_gauss_state(U)` con 9 puntos por defecto.
- **Spec**: [docs/specs/Quad8.md](specs/Quad8.md).
- **Archivo**: [solidum/elements/solid_2d/quad8.py](../solidum/elements/solid_2d/quad8.py) (subclase de la base interna `_HigherOrderSolid2D` definida en `_shared.py`, compartida con Quad9 y Tri6).

---

## Quad9 — cuadrilátero Lagrangiano 2D de orden 2

- **Propósito**: continuo 2D cuadrático con espacio polinómico completo $Q_2$; añade un noveno nodo central interior al Quad8.
- **DOFs por nodo**: `['ux', 'uy']` · 9 nodos · `STRAIN_DIM = 3`.
- **Funciones de forma**: producto tensorial Lagrange 1D-1D.
- **Integración**: Gauss 3×3.
- **Cargas y salida por Gauss**: idénticas a Quad8 (el nodo 8 es interior y no participa en bordes).
- **Spec**: [docs/specs/Quad9.md](specs/Quad9.md).
- **Archivo**: [solidum/elements/solid_2d/quad9.py](../solidum/elements/solid_2d/quad9.py) (subclase de la base interna `_HigherOrderSolid2D`, compartida con Quad8 y Tri6).

---

## Tri6 — triángulo 2D cuadrático completo P₂

- **Propósito**: triángulo isoparamétrico cuadrático que cura el shear locking del Tri3; reproduce campos cuadráticos exactamente.
- **DOFs por nodo**: `['ux', 'uy']` · 6 nodos (3 vértices + 3 medios de borde) · `STRAIN_DIM = 3`.
- **Integración**: 3 puntos en los puntos medios (cuadratura `tri_3`).
- **Cargas**: body load por cuadratura; tracción de borde reparte 1/6, 4/6, 1/6.
- **Spec**: [docs/specs/Tri6.md](specs/Tri6.md).
- **Archivo**: [solidum/elements/solid_2d/tri6.py](../solidum/elements/solid_2d/tri6.py) (subclase de la base interna `_HigherOrderSolid2D` en `_shared.py`, compartida con Quad8 y Quad9; declara `_MASS_QUADRATURE = "tri_6"` para integrar exactamente la masa consistente, que la cuadratura del elemento subintegraría).

---

# Elementos sólidos 3D (ADR 0012, Etapa 7)

Sólidos tridimensionales con convención Voigt 6D del proyecto (`Reglas.md §5`):
``[ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz]`` con ``γ_ij = 2·ε_ij``. Esta familia
**no expone `internal_forces`** — la primitiva natural de salida en sólidos
es ``compute_gauss_state(U)``, que devuelve ``{stress, strain, points_*}``
por punto de integración (ADR 0012: cierre del contrato `internal_forces`
por dominio explícito).

## Hex8 — hexaedro trilineal 3D

- **Propósito**: sólido 3D isoparamétrico de primer orden, espejo natural de Quad4.
- **DOFs por nodo**: `['ux', 'uy', 'uz']` · 8 nodos (orden VTK_HEXAHEDRON) · `STRAIN_DIM = 6` · `N_INTEGRATION_POINTS = 8` (default Gauss 2×2×2).
- **Cinemática**: matriz `B(ξ, η, ζ)` 6×24 por mapeo isoparamétrico estándar; deformación Voigt 3D.
- **Integración**: por defecto Gauss 2×2×2 (8 puntos). Configurable vía `quadrature` (`hex_3x3x3` para no-lineales severos, `hex_1x1x1` reducido con riesgo de hourglass).
- **Caras (ADR 0012)**: 6 caras numeradas con normal saliente — 0 (−ζ), 1 (+ζ), 2 (−η), 3 (+ξ), 4 (+η), 5 (−ξ). Ver `Hex8.FACE_NODES`.
- **Cargas distribuidas**: `compute_body_load(b)` integra `∫ Nᵀb dV`; `compute_face_traction(face, t̄)` reparte tracción uniforme sobre una cara con cuadratura 2D Gauss 2×2 (4 puntos por cara). `t̄` en globales.
- **Salida por Gauss**: `compute_gauss_state(U)` devuelve `{points_natural, points_global, strain (n_g, 6), stress (n_g, 6)}`. Sin `internal_forces` (ADR 0012).
- **Implementación**: kernels `_compute_kinematics_hex8`, `_det_jacobian_hex8`, `_shape_functions_hex8` con `@njit` en `solidum/elements/solid_3d/_shared.py`. Mass lumping vía `lump_hrz` con `n_translational_dirs=3` y `_expand_scalar_mass_3d`.
- **Limitaciones declaradas**:
  - **Locking volumétrico** con ν → 0.5 — sin mitigación (política idéntica al Quad4 2D). Blindado por `tests/test_volumetric_locking_3d.py`.
  - **Hourglass** con integración reducida `hex_1x1x1`: 12 modos espurios; sin estabilización Flanagan-Belytschko.
  - **Shear locking** con malla coarse 1 capa en sección: documentado en `tests/validation/test_macneal_beam_3d.py`.
- **Spec**: [docs/specs/Hex8.md](specs/Hex8.md).
- **Archivo**: [solidum/elements/solid_3d/hex8.py](../solidum/elements/solid_3d/hex8.py).

---

## Tet4 — tetraedro lineal 3D (CST 3D)

- **Propósito**: sólido 3D isoparamétrico de primer orden, espejo natural de Tri3. Útil para mallas no estructuradas y transiciones.
- **DOFs por nodo**: `['ux', 'uy', 'uz']` · 4 nodos (orden VTK_TETRA, volumen positivo) · `STRAIN_DIM = 6` · `N_INTEGRATION_POINTS = 1`.
- **Cinemática**: `B` constante 6×12 (deformación uniforme dentro del elemento — CST 3D). Reproduce campos lineales exactamente.
- **Integración**: 1 punto baricéntrico con peso 1/6 (exacto por linealidad de B).
- **Caras (ADR 0012)**: 4 caras triangulares numeradas con normal saliente; cara `i` opuesta al nodo `i`. Ver `Tet4.FACE_NODES`.
- **Cargas distribuidas**: `compute_body_load(b)` reparte `b·V_e/4` a cada nodo (exacto); `compute_face_traction(face, t̄)` reparte `A_cara·t̄/3` a cada uno de los 3 nodos de la cara, cero al nodo opuesto.
- **Salida por Gauss**: `compute_gauss_state(U)` con 1 punto (centroide); σ y ε constantes sobre el elemento.
- **Implementación**: kernel `_compute_kinematics_tet4` con `@njit` en `_shared.py`. Masa consistente analítica `M_ij = ρ·V_e·(1+δ_ij)/20`; lumped HRZ canónico (= row-sum por simetría).
- **Limitaciones declaradas**:
  - **Shear locking severo** en flexión — peor que Hex8 por la pobreza del espacio lineal. Recomendación: usar Hex8 en mallas hexaédricas; Tet4 para transiciones de malla no estructurada o cuando geometría compleja impide hexaedros.
  - **Locking volumétrico** con ν → 0.5 — aún peor que en Hex8.
- **Spec**: [docs/specs/Tet4.md](specs/Tet4.md).
- **Archivo**: [solidum/elements/solid_3d/tet4.py](../solidum/elements/solid_3d/tet4.py).

---

# Elementos con discontinuidad embebida (ADR 0010)

Subfamilia de elementos con **DOFs enriquecidos elementales** y **condensación estática local**: el ensamblador no ve los grados de libertad del salto. Cuando se cumple un criterio de activación (Rankine en fase 1), el elemento materializa una discontinuidad interna `Γ_d` y enriquece su cinemática con un salto `[[u]]` gobernado por un material cohesivo (`CohesiveMaterial`, ver el [catálogo de materiales](catalogo_materiales.md) §"Materiales cohesivos"). La activación se evalúa en el hook `Element.prepare_step(U_committed)` que los solvers no lineales invocan **una vez por paso**, antes del Newton (anti-chattering, ADR 0010 §5).

## CST_Embedded2D — triángulo CST con discontinuidad interior embebida (KOS)

- **Propósito**: introducir fractura computacional en aproximación discreta sobre el CST padre. Fiel a Retama (2010), Caps. 2, 5, 6 y 7: cinemática KOS, condensación estática local, longitud efectiva `l_d = (A_e/h)·cos(θ−α)`.
- **DOFs por nodo**: `['ux', 'uy']` · 3 nodos · `STRAIN_DIM = 3` · `N_INTEGRATION_POINTS = 1`.
- **DOFs enriquecidos (elementales, no globales)**: `[[u]] ∈ ℝ²` en frame local `(n, s)` de `Γ_d`. Se condensan dentro de `compute_element_state` y nunca llegan al ensamblador.
- **Estado intacto**: bit-exact con [Tri3](#tri3--triángulo-lineal-2d-cst) (mismo `B`, mismo material, misma cuadratura).
- **Estado agrietado**: Newton local sobre `[[u]]` hasta `R^{[[u]]} = 0`; después la condensación `K_cond = K_dd − K_du · K_{[[u]][[u]]}⁻¹ · K_du^T` se devuelve al ensamblador. El bulk descarga elásticamente conforme `[[u]]` crece (discrete approach); la disipación va al cohesivo en `Γ_d`.
- **Activación**: criterio Rankine (`σ_I > σ_t0` del cohesivo) en el centroide del CST, con el estado convergido del paso anterior. **Irreversible**: una vez activada, la discontinuidad persiste aunque el siguiente paso descargue.
- **Materiales aceptados**: bulk **solo `Elastic2D`** en fase 1 (la discrete approach presupone bulk elástico, ADR 0010); cohesivo cualquier `CohesiveMaterial` con `JUMP_DIM = 2`.
- **Estado de la discontinuidad**: `DiscontinuityState` (en [solidum/core/discontinuity_state.py](../solidum/core/discontinuity_state.py)) — paralela a `ElementState`, semántica trial/commit.
- **YAML**:
  ```yaml
  cohesive_materials:
    - {id: 1, type: CohesiveDamageIsotropic, sigma_t0: 2.5e6, G_f: 100.0,
       K_e: 1.0e13, softening: linear}
  elements:
    - {id: 1, type: CST_Embedded2D, nodes: [1, 2, 3],
       material: 1, cohesive_material: 1}
  ```
- **Post-procesamiento**: `compute_gauss_state(U)` añade clave `'discontinuity'` con `{normal, tangent, centroid, solitary_node, l_d, jump, traction, damage}` cuando el elemento está agrietado.
- **Limitaciones declaradas** (`out_of_scope` en la spec): modo mixto I-II (fase G del ADR 0010), reorientación de `n` tras activación (tracking no trivial, fase F), múltiples discontinuidades por elemento, contacto unilateral en compresión, bulks no elásticos, orden superior (`Tri6_Embedded`), 3D (`Tet4_Embedded`).
- **Spec**: [docs/specs/CST_Embedded2D.md](specs/CST_Embedded2D.md).
- **Archivo**: [solidum/elements/solid_2d/embedded_cst.py](../solidum/elements/solid_2d/embedded_cst.py).

---

## Cómo añadir un elemento nuevo

1. **Spec primero** — el usuario crea `docs/specs/<Nombre>.md` a partir de `docs/specs/_template_element.md` (especificación física + formulación + contrato YAML). Sin spec, la IA no escribe código.
2. **Scaffolding** — `/solidum-new element <Nombre>` genera archivo en `solidum/elements/`, decorador `@ElementRegistry.register` y esqueleto de test.
3. **Implementación + validación** — la IA codifica contra la spec; los tests cubren los casos de `acceptance` declarados.
4. **Catálogo** — cuando la spec pasa a `status: validated`, se añade aquí una entrada breve siguiendo el formato de arriba (la spec sigue siendo la referencia detallada).
