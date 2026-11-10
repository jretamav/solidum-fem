# BondSlipMultilinear — adherencia-deslizamiento con envolvente multilineal

> Orden de trabajo. El usuario escribe la **especificación física**, la **formulación numérica** y el **contrato**. La IA rellena **implementación** y responde en **diálogo** durante el trabajo.
>
> **Variante** de [`BondSlipModelCode`](BondSlipModelCode.md) (Reglas §4): la misma ley salvo la envolvente monótona, que da el usuario por puntos. Aquí sólo se documenta lo que cambia.

---

## Especificación física

### 0. Descripción general
Ley de adherencia-deslizamiento de línea (ADR 0033), con una envolvente que el usuario da por puntos `(s_k, τ_k)`. Sirve para una curva medida en un ensayo de arrancamiento o para la de un código que no tenga función auxiliar. Igual que en `BondSlipModelCode`:
- componente normal elástica lineal desacoplada (§3 de aquella);
- descarga y recarga R1 sin degradación (§4 de aquella).

### 1. Envolvente
$$\tau_{env}(a) = \text{interpolación lineal de } \{(0, 0), (s_1, \tau_1), \dots, (s_m, \tau_m)\}, \qquad \tau_{env}(a) = \tau_m \ \text{para}\ a > s_m$$
con `0 < s_1 < s_2 < … < s_m` y `τ_k ≥ 0`. La rigidez inicial y de descarga es la del primer tramo, `k0 = τ_1/s_1`. La fricción de descarga es `τ_f`, por omisión la residual `τ_m`, con `τ_f ≤ max τ_k`. La envolvente puede ablandarse (`τ_{k+1} < τ_k`).

---

## Formulación numérica
La de `BondSlipModelCode` (§4 a §6) con esta `τ_env`. La tangente en la envolvente es la pendiente del tramo; en un vértice se toma la del tramo hacia el que avanza el deslizamiento.

---

## Contrato de implementación

```yaml
name: BondSlipMultilinear
kind: interface_material
status: draft

interface:
  geometry: line
  jump_dim: null
  primary_state_var: slip_max
  is_symmetric: true

admission:
  textbooks:
    - de_borst_crisfield: "§13.1 elementos de interfaz (marco); el criterio no decide la adherencia por los textos y lo decidió el usuario"
  general_codes:
    - diana: "BONDSL 3, ley multilineal (DISTAU)"

parameters:
  - { name: points, type: "list[[float, float]]", required: true,  desc: "envolvente [(s_k [m], τ_k [Pa])], s creciente desde s_1 > 0" }
  - { name: tau_f,  type: float,                  required: false, desc: "fricción de descarga [Pa]; por omisión τ_m" }
  - { name: k_n,    type: float,                  required: true,  desc: "rigidez normal [Pa/m], > 0" }

acceptance:
  verification:
    - name: igual_a_model_code_con_sus_puntos
      setup: "envolvente de Eligehausen et al. 1983 con la rama potencial muestreada en 400 puntos; historia del ciclo de la Fig. 5.3b"
      expect: "misma respuesta que BondSlipModelCode con esos parámetros, salvo el error de interpolación de la rama potencial"
      tol_rel: 1.0e-3
    - name: envolvente_por_puntos
      setup: "tres puntos con ablandamiento, carga monótona en ambos sentidos"
      expect: "τ = interpolación lineal; τ_m más allá del último punto; τ(−s) = −τ(s)"
      tol_rel: 1.0e-12
    - name: tangente_y_disipacion
      setup: "como en BondSlipModelCode"
      expect: "tangente por diferencias finitas en cada tramo y en descarga, fricción y recarga; disipación ≥ 0 en historias aleatorias"
      tol_rel: 1.0e-6
  specific:
    - name: rechazos
      setup: "s no creciente, s_1 ≤ 0, τ_k < 0, τ_f > max τ_k, k_n ≤ 0, menos de un punto"
      expect: "ValueError al construir"

references:
  - "DIANA FEA, Material Library §26.3, bond-slip multilineal (BONDSL 3)."
  - "Eligehausen, Popov, Bertero (1983), Fig. 5.3 (reglas de descarga)."
```

---

## Implementación

*Rellena la IA tras programar.*

---

## Diálogo

- **2026-11-10** · Aprobada en la decisión 2 del usuario junto con `BondSlipModelCode`.
- **2026-11-10** · Spec **aprobada por el usuario** («si»).
