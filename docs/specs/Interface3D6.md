# INTERFACE 3D6 — elemento de interfaz de espesor nulo con caras triangulares (extensión de `Interface2D4`)

> Spec corta tipo extensión (Reglas §4): la formulación de la familia está en [`Interface2D4`](Interface2D4.md) y la arquitectura en el ADR 0043. Aquí sólo lo que cambia.

## Lo que cambia respecto de `Interface2D4`

- **Caras**: triángulos lineales de 3 nodos. Cara A = (1, 2, 3), cara B = (4, 5, 6), con 4 frente a 1, 5 frente a 2 y 6 frente a 3 (orden de COH3D6 de Abaqus; celda VTK `wedge`).
- **Ejes** (§2 de la familia): `n = (x_2 − x_1) × (x_3 − x_1)/‖·‖` sobre la cara A, `s1 = (x_2 − x_1)/‖·‖`, `s2 = n × s1`; constantes en la cara.
- **Salto**: 3 componentes `(g_n, g_s1, g_s2)`; la ley debe aceptar `("surface", 3)`.
- **Integración**: nodal, los tres vértices con peso 1/6 (área del triángulo de referencia 1/2); Gauss, 3 puntos interiores `(1/6, 1/6), (2/3, 1/6), (1/6, 2/3)` con peso 1/6.
- **Medida**: `J = ‖(x_2 − x_1) × (x_3 − x_1)‖` de la cara A; sin espesor.

## Contrato de implementación

```yaml
name: Interface3D6
kind: element
status: implemented

interface:
  field: interface
  dof_names: [ux, uy, uz]
  n_nodes: 6
  jump_dim: 3
  n_integration_points: 3

admission:
  # Variante 3D de la familia admitida en Interface2D4 (informes/bibliotecario/busqueda_contacto_2026-11-26.md).
  # COH3D6 de Abaqus: de memoria del bibliotecario, sin abrir el manual; confirmar antes de status: validated.
  textbooks:
    - de_borst_crisfield: "§13.1 Interface Elements (formulación general de la familia)"
  general_codes:
    - abaqus: "COH3D6 (cohesivo de 6 nodos; por confirmar en el manual)"

parameters:
  - { name: material,    type: int, required: true,  desc: "Id en interface_materials; ley de superficie con 3 componentes" }
  - { name: integration, type: str, required: false, desc: "'nodal' | 'gauss' (Diálogo Q1 de Interface2D4)" }

acceptance:
  verification:
    - name: parcela_de_contacto_3d
      setup: "dos bloques de Tet4 con mallas coincidentes en la interfaz, presión uniforme, ContactPenalty sin fricción; nodal y Gauss"
      expect: "presión de contacto p en todos los puntos; σ_zz = −p uniforme"
      tol_rel: 1.0e-10
  specific:
    - name: cuerpo_rigido_3d
      setup: "traslación y giro rígido pequeño con nodos coincidentes"
      expect: "g = 0 y f_e = 0"
      tol_rel: 1.0e-12
    - name: tangente_del_elemento_3d
      setup: "diferencias centradas en adherencia y deslizamiento 3D"
      expect: "coinciden"
      tol_rel: 1.0e-6

references:
  - "Spec Interface2D4 (formulación de la familia); ADR 0043."
```

## Implementación

- Archivo: `solidum/elements/interface/interface_3d.py`; base común en `_shared.py` (ver `Interface2D4`).
- Clase: `Interface3D6` (`FACE_CELL = "tri"`, celda VTK `wedge`).
- Tests: `tests/test_interface_elements.py`, `test_interface_generation.py`, `test_interface_output.py`, barrido de contrato.
- Verificación: `informes/teorico/verificacion_teorica_2026-11-26_contacto.md` (parcela 3D con Tet4, nodal y Gauss, `g_0 = 0` y `> 0`: ≤ 1.6e-14; *medido* por el teórico).
- Notas de traducción: orden de nodos de COH3D6 (Abaqus). Sin comprobar: si la base de la celda `wedge` de VTK tiene la normal al revés que este orden (sólo visualización).

## Diálogo

- *(vacío; las decisiones de la familia están en el Diálogo de `Interface2D4`)*
