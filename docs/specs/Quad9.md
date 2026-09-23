# Quad9 — cuadrilátero Lagrangiano 2D de orden 2

> Variante Lagrangiana del Quad8: añade un noveno nodo central interior.
> Espacio polinómico completo $Q_2$ (incluye el término $\xi^2\eta^2$ que el serendípito omite).

---

## Especificación física

Idéntica a `Quad8` salvo por las funciones de forma.

### 7. Funciones de forma
Producto tensorial de los polinomios de Lagrange 1D de orden 2:
$$L_0(\xi) = \tfrac{1}{2}\xi(\xi - 1), \quad L_1(\xi) = 1 - \xi^2, \quad L_2(\xi) = \tfrac{1}{2}\xi(\xi + 1).$$
$N_i(\xi, \eta) = L_a(\xi)\,L_b(\eta)$ con $(a, b)$ asociado al nodo $i$ según la numeración:

- 0..3: vértices (-1,-1), (+1,-1), (+1,+1), (-1,+1).
- 4: medio del borde 0-1 (η = -1).
- 5: medio del borde 1-2 (ξ = +1).
- 6: medio del borde 2-3 (η = +1).
- 7: medio del borde 3-0 (ξ = -1).
- 8: nodo central interior (ξ, η) = (0, 0) — *bubble*.

### 11. Cuadratura
Gauss $3 \times 3$ por defecto.

### 12. Cargas distribuidas
Body load con la cuadratura del elemento. Tracción en bordes idéntica a `Quad8` (1/6, 4/6, 1/6 — el nodo 8 no participa, no toca bordes).

---

## Contrato de implementación

```yaml
name: Quad9
kind: element
status: validated

interface:
  dof_names: [ux, uy]
  n_nodes: 9
  strain_dim: 3
  n_integration_points: 9

acceptance:
  verification:
    - name: patch_cuadratico
      setup: "u = (x², 0) impuesto en los 9 nodos"
      expect: "ε_xx = 2x exacto en todos los Gauss"
      tol_rel: 1.0e-10
    - name: Cook's membrane (Cook 1974)                              # 2026-05-19
      setup: "trapezoid (0,0)-(48,44)-(48,60)-(0,44), E=1, ν=1/3, plane stress, cortante total F=1 distribuido en borde derecho; malla 4×4 Q9"
      expect: "u_y en (48,52) ≈ 23.91 ± 0.5"
      tol_abs: 0.5
```

---

## Implementación

- **Archivo**: [solidum/elements/solid_2d/quad9.py](../../solidum/elements/solid_2d/quad9.py) · clase `Quad9` (subclase de la base interna `_HigherOrderSolid2D` en [_shared.py](../../solidum/elements/solid_2d/_shared.py), compartida con Quad8 y Tri6).
- **Tests**: [tests/test_higher_order_solid_2d.py](../../tests/test_higher_order_solid_2d.py).

## Diálogo

- **2026-09-22** · Auditoría global: `quadrature` acepta clave del registro o tupla `(points, weights)` en todos los sólidos (`resolve_quadrature` valida la familia por medida de referencia y dimensión). `_MASS_QUADRATURE = "3x3"` explícita: la masa se integra con la regla completa aunque K use 2×2 (antes M de rango 8/18). tolerancia del jacobiano ahora **relativa** (`det J ≤ JACOBIAN_RTOL · Π‖fila_i(J)‖`, cota de Hadamard, adimensional): la absoluta `1e-10` sobre `det J` rechazaba mallas finas en metros sin distorsión.
- **2026-09-22** · Ensamblaje por lotes (ADR 0014): la cinemática del elemento pasa a un núcleo compilado compartido (`_kin2d_from_dN`) y se expone como `BATCH_KINEMATICS`; `compute_element_state` no cambia de formulación. El `Assembler` agrupa los elementos en familias y los evalúa en un kernel único; equivalencia con el camino por elemento verificada a 1e-14 en `K` y `F_int` sobre todos los materiales compatibles (`tests/test_batch_assembly.py`).
- **2026-09-22** · Cierre del ADR 0014 (§9): la cinemática por lotes (`BATCH_KINEMATICS`) ya **no lanza** con jacobiano degenerado —devuelve `det J ≤ 0` y el ensamblador lanza el `ValueError` con el id del elemento—, porque el kernel de familia tiene ahora variante paralela (`prange`, bit a bit igual a la serie) y una excepción dentro de la región paralela se pierde; el camino por elemento sigue lanzando como siempre. Declara las funciones de forma para el post-proceso por familia (`Family.gauss_state`, mismo resultado que `compute_gauss_state`, verificado en `tests/test_batch_postprocess.py`). Formulación sin cambios.
