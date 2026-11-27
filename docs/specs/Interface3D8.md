# INTERFACE 3D8 — elemento de interfaz de espesor nulo con caras cuadriláteras (extensión de `Interface2D4`)

> Spec corta tipo extensión (Reglas §4): la formulación de la familia está en [`Interface2D4`](Interface2D4.md) y la arquitectura en el ADR 0043. Aquí sólo lo que cambia.

## Lo que cambia respecto de `Interface2D4`

- **Caras**: cuadriláteros bilineales de 4 nodos. Cara A = (1, 2, 3, 4), cara B = (5, 6, 7, 8), con 5 frente a 1, 6 frente a 2, etc. (orden de COH3D8 de Abaqus; celda VTK `hexahedron`).
- **Ejes** (§2 de la familia): `a_ξ = ∂x_A/∂ξ`, `a_η = ∂x_A/∂η`, `n = a_ξ × a_η/‖·‖`, `s1 = a_ξ/‖a_ξ‖`, `s2 = n × s1`, evaluados en cada punto de Gauss; con la regla nodal, los del nodo `i` se evalúan en `(ξ_i/3, η_i/3)` (abajo). Varían dentro de una cara alabeada.
- **Salto**: 3 componentes `(g_n, g_s1, g_s2)`; la ley debe aceptar `("surface", 3)`.
- **Integración**: nodal, en las cuatro esquinas, con el **vector de área tributaria consistente** de cada nodo, `m_i n_i = ∫N_i (a_ξ × a_η) dξdη`, que en la cara bilineal es exactamente `(a_ξ × a_η)` en `(ξ_i/3, η_i/3)` (identidad comprobada con SymPy por el teórico): la medida del nodo es su norma y la normal del nodo, su dirección (decisión del usuario del 2026-11-26, Diálogo). En una cara plana `J` es lineal y la medida es el área tributaria `∫N_i J dξdη` de la decisión Q3 de `Interface2D4` (Newton–Cotes reparte mal el área en caras que no son paralelogramos); en una alabeada, además, la normal es la que equilibra a los sólidos. El salto del nodo sigue siendo el de su par: la regla no acopla pares. Gauss, 2×2.
- **Medida**: `J = ‖a_ξ × a_η‖` de la cara A; sin espesor.

## Contrato de implementación

```yaml
name: Interface3D8
kind: element
status: implemented

interface:
  field: interface
  dof_names: [ux, uy, uz]
  n_nodes: 8
  jump_dim: 3
  n_integration_points: 4

admission:
  # Variante 3D de la familia admitida en Interface2D4 (informes/bibliotecario/busqueda_contacto_2026-11-26.md).
  # COH3D8 de Abaqus: de memoria del bibliotecario, sin abrir el manual; confirmar antes de status: validated.
  textbooks:
    - de_borst_crisfield: "§13.1 Interface Elements (formulación general de la familia)"
  general_codes:
    - abaqus: "COH3D8 (cohesivo de 8 nodos; por confirmar en el manual)"

parameters:
  - { name: material,    type: int, required: true,  desc: "Id en interface_materials; ley de superficie con 3 componentes" }
  - { name: integration, type: str, required: false, desc: "'nodal' (por omisión) | 'gauss' (Diálogo Q1 de Interface2D4)" }

acceptance:
  verification:
    - name: parcela_de_contacto_3d
      setup: "dos bloques de Hex8 con mallas coincidentes (distorsionadas, caras planas) en la interfaz, presión uniforme, ContactPenalty sin fricción; nodal y Gauss"
      expect: "presión de contacto p en todos los puntos; σ_zz = −p uniforme"
      tol_rel: 1.0e-10
  specific:
    - name: cuerpo_rigido_3d
      setup: "traslación y giro rígido pequeño con nodos coincidentes"
      expect: "g = 0 y f_e = 0"
      tol_rel: 1.0e-12
    - name: tangente_del_elemento_3d
      setup: "diferencias centradas en adherencia y deslizamiento 3D, cara alabeada"
      expect: "coinciden"
      tol_rel: 1.0e-6
    - name: area_tributaria_en_cara_trapezoidal
      setup: "tracción uniforme sobre la cara plana trapezoidal (0,0), (2,0), (1.6,1.1), (0.2,0.7) del teórico; integración nodal"
      expect: "fuerzas nodales = ∫N_i J (0.375, 0.4417, 0.400, 0.3333 por unidad de tracción)"
      tol_rel: 1.0e-12
    - name: cara_alabeada
      setup: "cara de 4 nodos no coplanares"
      expect: "n(ξ) normal unitaria de a_ξ × a_η en cada punto; ejes dextrógiros"
      tol_rel: 1.0e-12
    - name: parcela_hidrostatica_con_cara_alabeada
      setup: "dos bloques de Hex8 del mismo material separados por un paraboloide hiperbólico (7.2° de alabeo por cara), presión hidrostática p, K_n grande; nodal y Gauss"
      expect: "t_n = −p y σ = −p·I (error de la penalización); con la normal de la esquina la regla nodal daba 1.8e-2 y 8.7e-2"
      tol_rel: 1.0e-4

references:
  - "Spec Interface2D4 (formulación de la familia); ADR 0043."
```

## Implementación

- Archivo: `solidum/elements/interface/interface_3d.py`; base común en `_shared.py` (ver `Interface2D4`).
- Clase: `Interface3D8` (`FACE_CELL = "quad"`, celda VTK `hexahedron`; ejes evaluados en cada punto de Gauss o, con la regla nodal, en `(ξ_i/3, η_i/3)`: `axes_pts` en `InterfaceElement._build_geometry`).
- Tests: `tests/test_interface_patch.py` (parcela de contacto 3D con caras planas no paralelogramo, `g_0 = 0` y `> 0`, nodal y Gauss; parcela hidrostática con cara alabeada; medida nodal = `∫N_i J` en cara plana), `tests/test_interface_elements.py`, `test_interface_generation.py`, `test_interface_output.py`, barrido de contrato. Mutaciones del paso 7: 7/7 detectadas (`informes/prototipos_2026-11-26/contacto/mutar_paso7.py`).
- Verificación: `informes/teorico/verificacion_teorica_2026-11-26_contacto.md` (parcela 3D con Hex8 de caras planas no paralelogramo, nodal y Gauss: ≤ 2.6e-14; la regla de Newton–Cotes la rompe, 4.3e-2; cara alabeada con la regla nodal de la esquina: hallazgo 1, corregido).
- Notas de traducción: orden de nodos de COH3D8 (Abaqus).

## Diálogo

- Las decisiones de la familia están en el Diálogo de `Interface2D4`.
- **2026-11-26 — Normal y medida nodales en la cara alabeada. DECIDIDO por el usuario («ok la 1»): opción 1, el vector de área tributaria consistente.** Hallazgo de la verificación de `teorico-solidum` (`informes/teorico/verificacion_teorica_2026-11-26_contacto.md`, hallazgo 1): con la normal de la esquina y la medida `∫N_i J`, una cara alabeada no transmite una presión uniforme (error en los esfuerzos de los sólidos de 3 % con 2.4° de alabeo y 8.7 % con 7.2°, que no baja con `K_n`; con torsión y curvatura, 17–18 % en la tracción puntual, que no baja al refinar; *medido*), porque el vector que equilibra a los sólidos es `∫N_i n J`, no `n(ξ_i)·∫N_i J`. Opciones: (1) normal y medida del nodo `i` en `(ξ_i/3, η_i/3)`, idéntica en caras planas, error de la penalización en la alabeada (2.3e-6 en `t_n`, 1.5e-8 en σ; en la superficie con torsión y curvatura, residuo de segundo orden que converge: 1.5e-3 con 3 divisiones, 3.8e-4 con 6); (2) mantenerla y documentar el límite, con Gauss para caras alabeadas. Recomendación de la IA y del teórico: (1), análogo 3D de la decisión Q1 de `Interface2D6`. No se buscó si algún código general la usa así. Implementado el mismo día (`tests/test_interface_patch.py`: alabeada 7.2° con error < 1e-4 en `t_n` y < 1e-6 en σ, *medido*).
