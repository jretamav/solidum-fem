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
- Tracción positiva: `ε > 0 ⇔ elongación ⇔ σ > 0`.
- `STRAIN_DIM = 1`: la deformación que recibe el material es un escalar (no Voigt).
- La orientación de los nodos no afecta el resultado (`K_e` y `F_int` son invariantes bajo su permutación).

### Parámetros
- `A` — área de la sección transversal.

### Régimen de validez
- Pequeñas deformaciones (`|ε| ≲ 10⁻²`) y pequeños desplazamientos.
- Cargas exclusivamente axiales. Si la carga transversal sobre la barra no es despreciable → `Frame2DEuler` / `Frame2DTimoshenko`.

### Validación
- Tests: [tests/test_structural.py](tests/test_structural.py) · `TestTruss2D`.
- Archivo fuente: [fenix/elements/structural.py](fenix/elements/structural.py) · clase `Truss2D`.
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
- Tests: [tests/test_structural.py](tests/test_structural.py) · `TestTruss2DCorot` (3 tests, uno por criterio de aceptación).
- Archivo fuente: [fenix/elements/structural.py](fenix/elements/structural.py) · clase `Truss2DCorot`.
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
- Tracción positiva: `ε > 0 ⇔ elongación ⇔ σ > 0`.
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
- Tests: [tests/test_structural.py](tests/test_structural.py) · `TestTruss3D`.
- Archivo fuente: [fenix/elements/structural.py](fenix/elements/structural.py) · clase `Truss3D`.
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
- Tests: [tests/test_structural.py](tests/test_structural.py) · `TestTruss3DCorot` (3 tests, uno por criterio de aceptación; el de rigidez geométrica verifica dos direcciones transversas independientes del plano perpendicular).
- Archivo fuente: [fenix/elements/structural.py](fenix/elements/structural.py) · clase `Truss3DCorot`.
- Spec: [docs/specs/Truss3DCorot.md](specs/Truss3DCorot.md).
- Referencias: Crisfield §3.3; Belytschko §4.5.

---

## Cable2DCorot — cable 2D corotacional

Elemento 1D inmerso en el plano que modela un cable: dos nodos, transmite solo tracción. Cinemática corotacional; la unilateralidad la aporta el material (típicamente `CableMaterial1D`).

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

## Frame2DEuler — pórtico/viga 2D Euler-Bernoulli

- **Propósito**: viga esbelta 2D que transmite axial + cortante + momento flector.
- **DOFs por nodo**: `['ux', 'uy', 'rz']` · 2 nodos · `STRAIN_DIM = 1`.
- **Cinemática**: hipótesis Euler-Bernoulli (secciones planas, perpendiculares al eje neutro deformado).
- **Integración**: matriz `K_local 6×6` analítica; estado constitutivo evaluado a partir de la deformación axial media `(u₃ − u₀)/L`.
- **Parámetros**: `A` (área), `I` (inercia respecto Z).
- **Limitaciones**: la no-linealidad material escala toda la rigidez por `E_t` (no captura rótulas plásticas distribuidas en sección); válido solo para vigas esbeltas (L/h ≳ 10).
- **Referencia**: Bathe, *Finite Element Procedures*, cap. 5.
- **Archivo**: [fenix/elements/structural.py](fenix/elements/structural.py)

---

## Frame2DTimoshenko — pórtico/viga 2D Timoshenko

- **Propósito**: viga 2D con deformación por cortante transversal incluida explícitamente.
- **DOFs por nodo**: `['ux', 'uy', 'rz']` · 2 nodos · `STRAIN_DIM = 1`.
- **Cinemática**: Timoshenko (secciones planas pero **no** perpendiculares al eje deformado). Factor `Φ = 12·E·I / (G·A_s·L²)` corrige la rigidez para evitar shear locking.
- **Integración**: matriz `K_local 6×6` analítica con corrección por `Φ`.
- **Parámetros**: `A`, `I`, `As` (área efectiva de cortante), `nu` (Poisson, opcional si el material lo expone).
- **Limitaciones**: misma que Euler en no-linealidad material (E_t escala todo).
- **Cuándo usarlo**: vigas gruesas (L/h ≲ 10), peraltadas, o cuando la flecha por cortante no es despreciable.
- **Referencia**: Reddy, *An Introduction to the Finite Element Method*, cap. 5.
- **Archivo**: [fenix/elements/structural.py](fenix/elements/structural.py)

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
