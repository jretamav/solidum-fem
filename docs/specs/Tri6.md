# Tri6 — triángulo 2D cuadrático completo $P_2$

> Triángulo isoparamétrico de 6 nodos. Cura el shear locking severo del `Tri3` y reproduce campos cuadráticos exactamente. Útil cuando la malla es triangular por geometría.

---

## Especificación física

Idéntica al `Tri3` salvo por la cinemática y la cuadratura.

### 7. Funciones de forma
Coordenadas baricéntricas $L_1 = 1 - \xi - \eta$, $L_2 = \xi$, $L_3 = \eta$:
$$N_i^{vert} = L_i(2 L_i - 1), \quad i = 1, 2, 3$$
$$N_4 = 4 L_1 L_2, \quad N_5 = 4 L_2 L_3, \quad N_6 = 4 L_3 L_1.$$

Numeración: 0,1,2 vértices en (0,0), (1,0), (0,1); 3 medio del borde 0-1, 4 medio del 1-2, 5 medio del 2-0.

### 11. Cuadratura
Cuadratura triangular de 3 puntos (puntos medios de los lados; peso 1/6 cada uno) — exacta para polinomios hasta grado 2 sobre el triángulo de referencia.

### 12. Cargas distribuidas
- Body load con la cuadratura del elemento.
- Edge traction uniforme: reparto 1/6, 4/6, 1/6 a (vértice, medio, vértice).

---

## Contrato de implementación

```yaml
name: Tri6
kind: element
status: validated

interface:
  dof_names: [ux, uy]
  n_nodes: 6
  strain_dim: 3
  n_integration_points: 3

acceptance:
  verification:
    - name: patch_cuadratico
      setup: "u = (x², 0) impuesto en los 6 nodos"
      expect: "ε_xx = 2x exacto en cada Gauss"
      tol_rel: 1.0e-10
    - name: cooks_membrane_4x4
      setup: "Trapezoide de Cook triangulado 4×4 (32 Tri6, diagonal c1→c3 con center node compartido)"
      expect: "u_y(48, 52) ≈ 23.91 (Bathe; Hughes)"
      tol_abs: 1.5
    - name: cooks_membrane_refinamiento
      setup: "Mallas 2×2 → 4×4 → 6×6 con la misma triangulación"
      expect: "Error |u_y − 23.91| decreciente monótonamente"
      tol_rel: 0.0
```

---

## Implementación

- **Archivo**: [fenix/elements/solid_2d/tri6.py](../../fenix/elements/solid_2d/tri6.py) · clase `Tri6` (subclase de la base interna `_HigherOrderSolid2D` en [_shared.py](../../fenix/elements/solid_2d/_shared.py), compartida con Quad8 y Quad9). Declara `_MASS_QUADRATURE = "tri_6"` (Dunavant 6 puntos, orden 4) porque la cuadratura del elemento (`tri_3`, orden 2) subintegra el producto cuadrático×cuadrático de la masa consistente.
- **Tests**: [tests/test_higher_order_solid_2d.py](../../tests/test_higher_order_solid_2d.py) (patch cuadrático), [tests/test_cooks_membrane.py](../../tests/test_cooks_membrane.py) (benchmark Bathe/Hughes con malla triangulada 4×4 + refinamiento monótono, añadido Fase B 2026-05-19).
