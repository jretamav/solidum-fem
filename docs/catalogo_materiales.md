# Catálogo de materiales

> Referencia rápida de los modelos constitutivos implementados. Una entrada por material. Para detalles del algoritmo → código fuente.
>
> **Convenciones**: `STRAIN_DIM` = dimensión Voigt de la deformación de entrada (1 = escalar, 3 = 2D `[ε_xx, ε_yy, γ_xy]`, 6 = 3D). `PRIMARY_STATE_VAR` = variable interna que el `VtkExporter` lleva al output.

---

## Elastic1D — elástico lineal 1D

- **Ley**: `σ = E · ε` (Hooke 1D).
- **STRAIN_DIM**: 1.
- **Parámetros**: `E` (módulo de Young).
- **Variables internas**: ninguna (sin historia).
- **Tangente**: constante = `E`.
- **Compatible con**: `Truss2D`, `Truss3D`, `Frame2DEuler`, `Frame2DTimoshenko`.
- **Archivo**: [fenix/materials/elastic.py](fenix/materials/elastic.py)

---

## Elastic2D — elástico lineal 2D

- **Ley**: `σ = C · ε`, con `C` el tensor constitutivo isotrópico en notación Voigt `[xx, yy, xy]`.
- **STRAIN_DIM**: 3.
- **Parámetros**: `E`, `nu`, `hypothesis ∈ {'plane_stress', 'plane_strain'}`.
- **Variables internas**: ninguna.
- **Tangente**: constante = `C`.
- **Compatible con**: `Quad4`, `Tri3`.
- **Archivo**: [fenix/materials/elastic_2d.py](fenix/materials/elastic_2d.py)

---

## IsotropicDamage1D — daño isotrópico 1D con softening exponencial

- **Modelo**: daño escalar `d ∈ [0, 1)`, módulo secante `E_sec = (1 − d) · E`, `σ = E_sec · ε`.
- **Evolución del daño** (Kuhn-Tucker):
  - Deformación equivalente: `ε_eq = |ε|`.
  - Variable histórica: `κ = max(κ_old, ε_eq)`.
  - Si `κ ≤ κ_0` → `d = 0`. Si no: `d = 1 − (κ_0/κ) · exp(−α·(κ − κ_0))`, saturado en `DAMAGE_MAX`.
- **STRAIN_DIM**: 1 · **PRIMARY_STATE_VAR**: `'damage'`.
- **Parámetros**: `E`, `kappa_0` (umbral elástico), `alpha` (velocidad de degradación).
- **Variables internas**: `kappa`, `damage`.
- **Tangente**: secante (`E_sec`); no consistente con tangente algorítmica del daño — adecuado para Newton-Raphson amortiguado, no óptimo para arc-length agresivo.
- **Compatible con**: `Truss2D`, `Truss3D`.
- **Referencia**: Lemaitre & Chaboche, *Mechanics of Solid Materials*; Oliver, "A consistent characteristic length for smeared cracking models" (IJNME 1989).
- **Archivo**: [fenix/materials/damage_1d.py](fenix/materials/damage_1d.py)

---

## IsotropicDamage2D — daño isotrópico 2D con softening exponencial

- **Modelo**: extensión 2D de `IsotropicDamage1D`. Tensor constitutivo secante `C_sec = (1 − d) · C_e`.
- **Deformación equivalente**: `ε_eq = √(ε_xx² + ε_yy² + ½·γ_xy²)` (norma simple sin distinguir tracción/compresión).
- **Evolución**: idéntica a 1D (κ máxima histórica + ley exponencial).
- **STRAIN_DIM**: 3 · **PRIMARY_STATE_VAR**: `'damage'`.
- **Parámetros**: `E`, `nu`, `kappa_0`, `alpha`, `hypothesis` (delegado a `Elastic2D` interno).
- **Variables internas**: `kappa`, `damage`.
- **Tangente**: secante.
- **Limitación**: la `ε_eq` simétrica no distingue daño en tracción vs compresión; para hormigón usar modelo de Mazars o split tracción/compresión (no implementado).
- **Compatible con**: `Quad4`, `Tri3`.
- **Archivo**: [fenix/materials/damage_2d.py](fenix/materials/damage_2d.py)

---

## Elastoplastic1D — plasticidad J2 1D con endurecimiento isotrópico lineal

- **Modelo**: plasticidad asociativa con criterio `f = |σ_trial| − (σ_y + H·α) ≤ 0`.
- **Algoritmo**: return mapping clásico 1D.
  - Predictor elástico: `σ_trial = E · (ε − ε_p_old)`.
  - Corrector: `Δγ = f / (E + H)`, `σ = σ_trial − Δγ·E·sign(σ_trial)`.
  - Tangente algorítmica consistente: `E_t = E·H / (E + H)`.
- **STRAIN_DIM**: 1 · **PRIMARY_STATE_VAR**: `'alpha'`.
- **Parámetros**: `E`, `sigma_y` (fluencia inicial), `H` (módulo de endurecimiento; 0 = perfectamente plástico).
- **Variables internas**: `eps_p` (deformación plástica), `alpha` (acumulada equivalente).
- **Compatible con**: `Truss2D`, `Truss3D`, `Frame2DEuler`, `Frame2DTimoshenko` (en estos últimos solo se aplica al esfuerzo axial).
- **Referencia**: Simo & Hughes, *Computational Inelasticity*, cap. 1.
- **Archivo**: [fenix/materials/plastic_1d.py](fenix/materials/plastic_1d.py)

---

## VonMises2D — plasticidad J2 2D con endurecimiento isotrópico lineal

- **Modelo**: J2 (Von Mises) en deformación plana; criterio `f = ‖s_trial‖ − √(2/3)·(σ_y + H·α) ≤ 0`.
- **Algoritmo**: return mapping radial sobre la parte desviadora.
  - Descomposición volumétrica/desviadora: `p = K·tr(ε)`, `s = 2G·e_dev`.
  - Predictor elástico desviador, corrector radial `Δγ = f_trial / (2G + ⅔·H)`.
  - Actualización: `s_new = s_trial − 2G·Δγ·N`, con `N = s_trial / ‖s_trial‖`.
  - Tangente algorítmica consistente derivada por linealización del return mapping.
- **STRAIN_DIM**: 3 · **PRIMARY_STATE_VAR**: `'alpha'`.
- **Parámetros**: `E`, `nu`, `sigma_y`, `H`, `hypothesis` (**obligatorio** `'plane_strain'`).
- **Variables internas**: `eps_p` (tensorial 4 componentes `[xx, yy, zz, xy]`), `alpha` (acumulada equivalente).
- **Restricción crítica**: solo deformación plana. Usar con `hypothesis='plane_stress'` lanza `NotImplementedError` (el return mapping para plane stress requiere algoritmo diferente).
- **Implementación**: kernel `_compute_j2_plasticity` con `@njit`.
- **Compatible con**: `Quad4`, `Tri3` (configurados en plane strain).
- **Referencia**: Simo & Hughes, *Computational Inelasticity*, cap. 3 (return mapping J2).
- **Archivo**: [fenix/materials/von_mises_2d.py](fenix/materials/von_mises_2d.py)

---

## Cómo añadir un material nuevo

`/fenix-new material <Name>` — genera archivo en `fenix/materials/`, decorador `@MaterialRegistry.register`, esqueleto de test.
Declarar **`STRAIN_DIM`** y, si tiene historia, **`PRIMARY_STATE_VAR`** (la variable que aparecerá en el VTK).
Tras implementar el modelo constitutivo, **añadir una entrada a este catálogo** siguiendo el formato de arriba.
