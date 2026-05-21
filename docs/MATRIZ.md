# Matriz de combinaciones

> Cuadrícula navegable de qué combinaciones del catálogo funcionan, cuáles están testeadas, y cuáles están vetadas (con razón).
>
> **Para el usuario** — antes de modelar, mira aquí qué combinación elemento × material × solver es segura.
> **Para la IA** — al añadir un componente nuevo, esta matriz dice qué celdas hay que cubrir con test.
>
> Leyenda:
> - **✓** = combinación válida y cubierta por al menos un test **de sistema** (elemento + material + solver acoplados, no test unitario aislado del material).
> - **○** = combinación válida según contratos (`STRAIN_DIM` compatible) pero sin test de sistema. Un test unitario del material aislado puede existir y no escala a ✓.
> - **(nota)** = combinación válida con restricción semántica (ver pie de la tabla). Se combina con ✓ u ○ para indicar si hay test.
> - *vacío* = combinación inválida (incompatibilidad de `STRAIN_DIM` o de semántica del modelo).

---

## 1. Elemento × Material

Compatibilidad determinada por `STRAIN_DIM` (1 = axial escalar, 3 = 2D Voigt `[ε_xx, ε_yy, γ_xy]`, 6 = 3D Voigt `[ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz]` — ADR 0012) y por la semántica del elemento.

| Elemento \ Material           | Elastic1D | Elastoplastic1D | IsotropicDamage1D | CableMaterial1D | Elastic2D | VonMises2D | DruckerPrager2D | IsotropicDamage2D | CohesiveDamageIsotropic | Elastic3D |
|------------------------------|:---------:|:---------------:|:-----------------:|:---------------:|:---------:|:----------:|:---------------:|:-----------------:|:-----------------------:|:---------:|
| `Truss2D`                    | ✓         | ✓               | ○                 |                 |           |            |                 |                   |                         |           |
| `Truss2DCorot`               | ✓         | ○               | ○                 |                 |           |            |                 |                   |                         |           |
| `Truss3D`                    | ✓         | ○               | ○                 |                 |           |            |                 |                   |                         |           |
| `Truss3DCorot`               | ✓         | ○               | ○                 |                 |           |            |                 |                   |                         |           |
| `Cable2DCorot`               |           |                 |                   | ✓               |           |            |                 |                   |                         |           |
| `Cable3DCorot`               |           |                 |                   | ✓               |           |            |                 |                   |                         |           |
| `Frame2DEuler`               | ✓         | (a) ○           |                   |                 |           |            |                 |                   |                         |           |
| `Frame2DTimoshenko`          | ✓         | (a) ○           |                   |                 |           |            |                 |                   |                         |           |
| `Frame2DEulerCorot`          | ✓         | (a) ○           |                   |                 |           |            |                 |                   |                         |           |
| `Frame3D`                    | ✓         | (a) ○           |                   |                 |           |            |                 |                   |                         |           |
| `Quad4`                      |           |                 |                   |                 | ✓         | ✓          | (b) ✓           | ✓                 |                         |           |
| `Tri3`                       |           |                 |                   |                 | ✓         | ✓          | (b) ✓           | ✓                 |                         |           |
| `Quad8`                      |           |                 |                   |                 | ✓         | ✓          | (b) ✓           | ✓                 |                         |           |
| `Quad9`                      |           |                 |                   |                 | ✓         | ✓          | (b) ✓           | ✓                 |                         |           |
| `Tri6`                       |           |                 |                   |                 | ✓         | ✓          | (b) ✓           | ✓                 |                         |           |
| `CST_Embedded2D`             |           |                 |                   |                 | (c) ✓     | (c) ○      | (c) ○           | (c) ○             | (c) ✓                   |           |
| `Hex8`                       |           |                 |                   |                 |           |            |                 |                   |                         | ✓         |
| `Tet4`                       |           |                 |                   |                 |           |            |                 |                   |                         | ✓         |

**Notas semánticas**:

- **(a)** Frames + plasticidad 1D: válido, pero la plasticidad se aplica **sólo al esfuerzo axial $\sigma$**. La fluencia por flexión no está modelada (espera `FiberSection`, ver deuda técnica de STATUS.md).
- **(b)** `DruckerPrager2D` está validado únicamente en `plane_strain`. `plane_stress` declarado out-of-scope (proyección con σ_zz=0 acoplada a flujo dilatante es notoriamente delicada).
- **(c)** `CST_Embedded2D` requiere **dos materiales**: un *bulk* 2D estándar (Elastic2D, VonMises2D, DruckerPrager2D, IsotropicDamage2D) que gobierna el continuo, y un *cohesivo* de la familia `CohesiveMaterial` (hoy sólo `CohesiveDamageIsotropic`) que gobierna el salto en la discontinuidad embebida. La fila refleja la compatibilidad con cada uno por separado; la combinación válida se especifica en YAML con dos campos (`bulk_material`, `cohesive_material`). Validado con Elastic2D + CohesiveDamageIsotropic (Etapa 5, ADR 0010).

**Casillas ○ (válidas no testeadas)**: los corotacionales (`Truss2DCorot`, `Truss3DCorot`) aceptan `IsotropicDamage1D` por contrato (`STRAIN_DIM=1`) pero los tests existentes lo combinan sólo con materiales lineales y plásticos. Cubrir si entra un caso de uso o si se decide cerrar el catálogo formalmente.

---

## 2. Solver × Tipo de problema

La elección del solver es **ortogonal al elemento**: depende de la linealidad del problema y del régimen de análisis.

| Solver                       | Pipeline kind  | Tipo de problema                                                  | Pre-condición sobre materiales/elementos                              |
|------------------------------|----------------|-------------------------------------------------------------------|------------------------------------------------------------------------|
| `LinearSolver`               | `static`       | Estático lineal (`K·U = F` en un paso).                           | Todos los materiales con tangente constante; sin corotacional; sin cable. |
| `NonlinearSolver`            | `static`       | Estático no lineal con control de carga (Newton-Raphson).         | Cualquier no-linealidad material o geométrica suave (sin snap-back).  |
| `ArcLengthSolver`            | `static`       | Estático no lineal con snap-through / snap-back / softening.      | Igual que `NonlinearSolver`, además captura puntos límite.            |
| `ModalSolver`                | `modal`        | Autovalor generalizado `K·φ = ω²M·φ` (frecuencias y modos).       | Todos los materiales (lineales en `u = 0`); `density` declarada; `compute_mass_matrix`. |
| `NewmarkSolver`              | `transient`    | Transitorio lineal `M·ü + C·u̇ + K·u = F(t)`.                      | Mismas restricciones que `LinearSolver` + `density` declarada.        |
| `HHTSolver`                  | `transient`    | Variante de Newmark con disipación numérica controlada (HHT-α).   | Idénticas a `NewmarkSolver` (es subclase).                            |
| `NewtonNewmarkSolver`        | `transient`    | Transitorio no lineal (Newton dentro de cada paso temporal).      | Cualquier material con historia + `density` declarada.                |
| `NewtonHHTSolver`            | `transient`    | Variante no lineal con disipación numérica HHT-α + Newton interno.| Idénticas a `NewtonNewmarkSolver` (es subclase).                      |
| `CentralDifferenceSolver`    | `transient`    | Transitorio explícito por diferencias centradas (leapfrog).       | Masa lumped (`lumping="lumped"`); CFL `Δt < 2/ω_max`; Frame3D oblicuo rechazado. |
| `HarmonicSolver`             | `harmonic`     | Respuesta forzada armónica `(-ω²M + iωC + K)·û = F̂` con barrido en `ω`. | Lineal; `density` declarada.                                          |
| `ResponseSpectrumSolver`     | `spectrum`     | Análisis sísmico por combinación modal SRSS/CQC contra espectro.  | Lineal; `density` declarada; suficientes modos (verificar `cumulative_effective_mass_ratio`). |

**Compatibilidad cruzada solver → elementos**: todos los elementos del catálogo implementan `compute_mass_matrix(lumping)` con `lumping ∈ {"consistent", "lumped"}` (ADR 0009 fases 1 y 2, cerradas 2026-05-18). Cualquier elemento es compatible con los 11 solvers anteriores siempre que el material declare `density`. **Excepción**: `CentralDifferenceSolver` requiere `lumping="lumped"` y rechaza Frame3D con eje oblicuo a los ejes globales (el bloque rotacional 3×3 del lumping no es estrictamente diagonal cuando `ρJp ≠ ρIy ≠ ρIz` — limitación documentada estándar).

---

## 3. Casos test representativos por combinación

Selección de tests "canónicos" que cubren combinaciones clave. La intención no es enumerar los 558 tests sino apuntar al fichero de referencia para cada celda no trivial.

| Combinación                                                    | Test representativo                                                                 |
|----------------------------------------------------------------|-------------------------------------------------------------------------------------|
| `Truss2D` + `Elastic1D` + `LinearSolver`                       | [`test_truss.py`](../tests/test_truss.py) · `TestTruss2D`                           |
| `Truss2DCorot` + `Elastic1D` + `NonlinearSolver`               | [`test_truss.py`](../tests/test_truss.py) · `TestTruss2DCorot`                      |
| `Truss3D` / `Truss3DCorot` + `Elastic1D`                       | [`test_truss.py`](../tests/test_truss.py) · `TestTruss3D` / `TestTruss3DCorot`      |
| `Truss2D` + `Elastoplastic1D` (Newton-Raphson + ArcLength)     | [`test_integration.py`](../tests/test_integration.py)                               |
| `Cable2D/3DCorot` + `CableMaterial1D`                          | [`test_cable_elements.py`](../tests/test_cable_elements.py) + [`test_cable_material.py`](../tests/test_cable_material.py) |
| `Frame2DEuler` + `Elastic1D`                                   | [`test_frame.py`](../tests/test_frame.py) · `TestFrame2DEulerAcceptance`            |
| `Frame2DTimoshenko` + `Elastic1D`                              | [`test_frame.py`](../tests/test_frame.py) · `TestFrame2DTimoshenkoAcceptance`       |
| `Frame2DEulerCorot` + `Elastic1D`                              | [`test_frame.py`](../tests/test_frame.py) · `TestFrame2DEulerCorotAcceptance`       |
| `Frame3D` + `Elastic1D`                                        | [`test_frame3d.py`](../tests/test_frame3d.py) · `TestFrame3DAcceptance`             |
| Fuerzas internas `Frame2D` / `Frame3D`                         | [`test_frame2d_internal_forces.py`](../tests/test_frame2d_internal_forces.py) · [`test_frame3d_internal_forces.py`](../tests/test_frame3d_internal_forces.py) |
| `Quad4` / `Tri3` + `Elastic2D`                                 | [`test_solid_2d.py`](../tests/test_solid_2d.py)                                     |
| `Quad8` / `Quad9` / `Tri6` + `Elastic2D`                       | [`test_higher_order_solid_2d.py`](../tests/test_higher_order_solid_2d.py) (patch cuadrático, isoparametría, traction de borde) |
| `Quad4` + `VonMises2D` (plane_strain + plane_stress)           | [`test_solid_2d_plasticity.py`](../tests/test_solid_2d_plasticity.py)               |
| `Quad4` + `DruckerPrager2D`                                    | [`test_solid_2d_drucker_prager.py`](../tests/test_solid_2d_drucker_prager.py)       |
| `Quad4` + `IsotropicDamage2D`                                  | [`test_solid_2d_damage.py`](../tests/test_solid_2d_damage.py)                       |
| `Tri3` / `Quad8` / `Quad9` / `Tri6` + materiales no lineales 2D (cobertura de cableado) | [`test_solid_2d_nonlinear_higher_order.py`](../tests/test_solid_2d_nonlinear_higher_order.py) |
| `Hex8` / `Tet4` + `Elastic3D` (ADR 0012)                       | [`test_solid_3d.py`](../tests/test_solid_3d.py) (unitarios) · [`test_cube_lame_3d.py`](../tests/validation/test_cube_lame_3d.py) (tracción uniaxial + hidrostática exactas) · [`test_macneal_beam_3d.py`](../tests/validation/test_macneal_beam_3d.py) (locking + convergencia h) |
| Modal — barra axial / viga simplemente apoyada / álgebra       | [`test_modal.py`](../tests/test_modal.py)                                           |
| Modal — Truss3D / Frame3D / Frame2DTimoshenko / Corot / Cables / Solid2D | [`test_modal_catalog.py`](../tests/test_modal_catalog.py)                 |
| Transitorio Newmark lineal                                     | [`test_newmark.py`](../tests/test_newmark.py)                                       |
| Transitorio Newton-Newmark (no lineal)                         | [`test_newmark_nonlinear.py`](../tests/test_newmark_nonlinear.py)                   |
| HHT-α (lineal y no lineal) — disipación numérica controlada    | [`test_hht.py`](../tests/test_hht.py)                                               |
| Mass lumping (HRZ canónico + nodal directo) — todos los elementos | [`test_mass_lumping.py`](../tests/test_mass_lumping.py)                          |
| Transitorio explícito por diferencias centradas                | [`test_central_difference.py`](../tests/test_central_difference.py)                 |
| Respuesta forzada armónica en frecuencia                       | [`test_harmonic.py`](../tests/test_harmonic.py)                                     |
| Análisis sísmico por combinación modal SRSS/CQC                | [`test_response_spectrum.py`](../tests/test_response_spectrum.py)                   |
| Plasticidad sólido 2D (unitario del material)                  | [`test_materials_unit.py`](../tests/test_materials_unit.py)                         |
| `solidum.run` y `solidum.run_yaml` end-to-end (estático + dinámico) | [`test_entry.py`](../tests/test_entry.py)                                          |
| Peso propio (`assemble_self_weight`, ADR 0008)                 | [`test_density_self_weight.py`](../tests/test_density_self_weight.py) · [`test_body_force_pipeline.py`](../tests/test_body_force_pipeline.py) · [`test_body_load_truss_frame.py`](../tests/test_body_load_truss_frame.py) |

---

## 4. Huecos visibles

Casillas **○** (válidas no testeadas) que el barrido sistemático revela. Priorizadas por valor de cobertura:

1. ~~**Sólidos 2D cuadráticos + materiales no lineales**~~ — **cerrado el 2026-05-14** con [`test_solid_2d_nonlinear_higher_order.py`](../tests/test_solid_2d_nonlinear_higher_order.py): 16 tests de cobertura del cableado para las 4 combinaciones de elemento (`Tri3`, `Quad8`, `Quad9`, `Tri6`) × 4 escenarios materiales (`VonMises2D` plane_strain + plane_stress, `DruckerPrager2D`, `IsotropicDamage2D`). Todos verdes a la primera; no había bug latente.
2. **Truss/Frame + plasticidad 1D acoplada con dinámica no lineal**: hoy `test_newmark_nonlinear.py` lo cubre sobre un sistema 1DOF; falta integración con elementos `Frame*` reales.
3. **Corotacional 1D + `IsotropicDamage1D`** (`Truss2DCorot`, `Truss3DCorot`): válido por contrato, sin test. **Coste**: bajo.
4. **`Truss3D` / `Truss3DCorot` + `Elastoplastic1D`**: hoy sólo `Truss2D` lo cubre (`test_integration.py`). **Coste**: bajo.
5. **Frames + `Elastoplastic1D`** (los 4 frames): válido sólo axialmente (nota a), sin test específico que ejercite plasticidad en un frame. **Coste**: bajo. La plasticidad por flexión es hueco *físico* que espera `FiberSection` (deuda #3 de STATUS.md) — no se confunde con este.
6. ~~**Sólidos 3D**: no figuran porque no existen aún.~~ — **abierto el 2026-05-19** con la Etapa 7 (ADR 0012): `Hex8`, `Tet4`, `Elastic3D` entran y la columna 3D aparece en la tabla §1. Próximos materiales 3D no lineales (`VonMises3D`, `DruckerPrager3D`, `IsotropicDamage3D`) abrirán columnas adicionales en sub-etapas posteriores; los elementos 3D cuadráticos (`Hex20`, `Hex27`, `Tet10`) abrirán filas.

**No-huecos** (combinaciones que parecen ausentes pero son decisiones documentadas):

- **`DruckerPrager2D` en plane_stress**: declarado out-of-scope.
- **`Cable*` + plasticidad o daño**: extensión futura, no hueco del catálogo actual.
- **Frames + materiales 2D**: incompatibilidad estructural por `STRAIN_DIM`, no hueco.

---

## Cómo se mantiene este documento

- **Componente nuevo**: añadir su columna o fila; marcar celdas con ✓ donde haya test y ○ donde sea compatible sin test. Si una compatibilidad requiere restricción, añadirla como nota al pie.
- **Test nuevo que cubre celda ○**: cambiar a ✓.
- **Combinación que se descubre incompatible** (durante implementación): marcar como nota explícita, no como ✓ ni vacío.
- **Etapa cerrada**: barrer la matriz para asegurar que todas las nuevas combinaciones quedaron reflejadas.

---

*Última actualización: 2026-05-18 — sincronización tras cierre completo de la Etapa 6 (ADR 0009 completo: 11 solvers, mass lumping en todos los elementos, dispatch declarativo `PIPELINE_KIND`, reglas C y D de auditoría aplicadas).*
