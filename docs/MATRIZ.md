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

> Los **elementos y materiales térmicos** (Etapa 8) no aparecen en esta tabla: forman una familia paralela cuya compatibilidad la fija `FLUX_DIM`, no `STRAIN_DIM`. Ver §1.b.

**Materiales 1D y 2D**:

| Elemento \ Material 1D/2D     | Elastic1D | Elastoplastic1D | IsotropicDamage1D | CableMaterial1D | Elastic2D | VonMises2D | DruckerPrager2D | IsotropicDamage2D | CohesiveDamageIsotropic |
|------------------------------|:---------:|:---------------:|:-----------------:|:---------------:|:---------:|:----------:|:---------------:|:-----------------:|:-----------------------:|
| `Truss2D`                    | ✓         | ✓               | ○                 |                 |           |            |                 |                   |                         |
| `Truss2DCorot`               | ✓         | ○               | ○                 |                 |           |            |                 |                   |                         |
| `Truss3D`                    | ✓         | ○               | ○                 |                 |           |            |                 |                   |                         |
| `Truss3DCorot`               | ✓         | ○               | ○                 |                 |           |            |                 |                   |                         |
| `Cable2DCorot`               |           |                 |                   | ✓               |           |            |                 |                   |                         |
| `Cable3DCorot`               |           |                 |                   | ✓               |           |            |                 |                   |                         |
| `Frame2DEuler`               | ✓         | (a) ○           |                   |                 |           |            |                 |                   |                         |
| `Frame2DTimoshenko`          | ✓         | (a) ○           |                   |                 |           |            |                 |                   |                         |
| `Frame2DEulerCorot`          | ✓         | (a) ○           |                   |                 |           |            |                 |                   |                         |
| `Frame3D`                    | ✓         | (a) ○           |                   |                 |           |            |                 |                   |                         |
| `Quad4`                      |           |                 |                   |                 | ✓         | ✓          | (b) ✓           | ✓                 |                         |
| `Tri3`                       |           |                 |                   |                 | ✓         | ✓          | (b) ✓           | ✓                 |                         |
| `Quad8`                      |           |                 |                   |                 | ✓         | ✓          | (b) ✓           | ✓                 |                         |
| `Quad9`                      |           |                 |                   |                 | ✓         | ✓          | (b) ✓           | ✓                 |                         |
| `Tri6`                       |           |                 |                   |                 | ✓         | ✓          | (b) ✓           | ✓                 |                         |
| `CST_Embedded2D`             |           |                 |                   |                 | (c) ✓     | (c) ○      | (c) ○           | (c) ○             | (c) ✓                   |

**Materiales 3D** (sub-etapa A.bis + Etapa 7):

| Elemento \ Material 3D | Elastic3D | VonMises3D | DruckerPrager3D | IsotropicDamage3D |
|------------------------|:---------:|:----------:|:---------------:|:-----------------:|
| `Hex8`                 | ✓         | ✓          | (d) ✓           | ✓                 |
| `Hex20`                | ✓         | ✓          | (d) ✓           | ✓                 |
| `Hex27`                | ✓         | ✓          | (d) ✓           | ✓                 |
| `Tet4`                 | ✓         | ✓          | (d) ✓           | ✓                 |
| `Tet10`                | ✓         | ✓          | (d) ✓           | ✓                 |

**Notas semánticas**:

- **(a)** Frames + plasticidad 1D: válido, pero la plasticidad se aplica **sólo al esfuerzo axial $\sigma$**. La fluencia por flexión no está modelada (espera `FiberSection`, ver deuda técnica de STATUS.md).
- **(b)** `DruckerPrager2D` está validado únicamente en `plane_strain`. `plane_stress` declarado out-of-scope (proyección con σ_zz=0 acoplada a flujo dilatante es notoriamente delicada).
- **(c)** `CST_Embedded2D` requiere **dos materiales**: un *bulk* 2D estándar (Elastic2D, VonMises2D, DruckerPrager2D, IsotropicDamage2D) que gobierna el continuo, y un *cohesivo* de la familia `CohesiveMaterial` (hoy sólo `CohesiveDamageIsotropic`) que gobierna el salto en la discontinuidad embebida. La fila refleja la compatibilidad con cada uno por separado; la combinación válida se especifica en YAML con dos campos (`bulk_material`, `cohesive_material`). Validado con Elastic2D + CohesiveDamageIsotropic (Etapa 5, ADR 0010).
- **(d)** `DruckerPrager3D` se ofrece con variantes `outer`/`inner` del cono; el ajuste `plane_strain_matched` del `DruckerPrager2D` **no aplica en 3D** y se rechaza explícitamente. La cobertura de los cuadráticos (`Hex20`, `Hex27`, `Tet10`) con materiales 3D no lineales quedó cerrada en la sub-fase 5 de A.ter: las 20 celdas de la tabla están en ✓.

**Casillas ○ (válidas no testeadas)**: los corotacionales (`Truss2DCorot`, `Truss3DCorot`) aceptan `IsotropicDamage1D` por contrato (`STRAIN_DIM=1`) pero los tests existentes lo combinan sólo con materiales lineales y plásticos. Cubrir si entra un caso de uso o si se decide cerrar el catálogo formalmente.

---

## 2. Solver × Tipo de problema

La elección del solver es **ortogonal al elemento**: depende de la linealidad del problema y del régimen de análisis.

| Solver                       | Pipeline kind  | Tipo de problema                                                  | Pre-condición sobre materiales/elementos                              |
|------------------------------|----------------|-------------------------------------------------------------------|------------------------------------------------------------------------|
| `LinearSolver`               | `static`       | Estático lineal (`K·U = F` en un paso).                           | Todos los materiales con tangente constante; sin corotacional; sin cable. |
| `NonlinearSolver`            | `static`       | Estático no lineal con control de carga (Newton-Raphson).         | Cualquier no-linealidad material o geométrica suave (sin snap-back).  |
| `ArcLengthSolver`            | `static`       | Estático no lineal con snap-through / snap-back / softening.      | Igual que `NonlinearSolver`, además captura puntos límite.            |
| `DissipationArcLengthSolver` | `static`       | Estático no lineal con softening severo, controlando la disipación de energía por paso (Gutiérrez 2004). Restricción lineal en vez de cuadrática; switching automático cilíndrico↔disipación. | Igual que `ArcLengthSolver` (es subclase). Validado para daño continuo bulk 1D/2D; **no** para cohesivo+embedded con penalty `K_e` rígido (deuda técnica #4). |
| `ModalSolver`                | `modal`        | Autovalor generalizado `K·φ = ω²M·φ` (frecuencias y modos).       | Todos los materiales (lineales en `u = 0`); `density` declarada; `compute_mass_matrix`. |
| `NewmarkSolver`              | `transient`    | Transitorio lineal `M·ü + C·u̇ + K·u = F(t)`.                      | Mismas restricciones que `LinearSolver` + `density` declarada.        |
| `HHTSolver`                  | `transient`    | Variante de Newmark con disipación numérica controlada (HHT-α).   | Idénticas a `NewmarkSolver` (es subclase).                            |
| `NewtonNewmarkSolver`        | `transient`    | Transitorio no lineal (Newton dentro de cada paso temporal).      | Cualquier material con historia + `density` declarada.                |
| `NewtonHHTSolver`            | `transient`    | Variante no lineal con disipación numérica HHT-α + Newton interno.| Idénticas a `NewtonNewmarkSolver` (es subclase).                      |
| `CentralDifferenceSolver`    | `transient`    | Transitorio explícito por diferencias centradas (leapfrog).       | Masa lumped (`lumping="lumped"`); CFL `Δt < 2/ω_max`; Frame3D oblicuo rechazado. |
| `HarmonicSolver`             | `harmonic`     | Respuesta forzada armónica `(-ω²M + iωC + K)·û = F̂` con barrido en `ω`. | Lineal; `density` declarada.                                          |
| `ResponseSpectrumSolver`     | `spectrum`     | Análisis sísmico por combinación modal SRSS/CQC contra espectro.  | Lineal; `density` declarada; suficientes modos (verificar `cumulative_effective_mass_ratio`). |

**Compatibilidad cruzada solver → elementos**: todos los elementos del catálogo implementan `compute_mass_matrix(lumping)` con `lumping ∈ {"consistent", "lumped"}` (ADR 0009 fases 1 y 2, cerradas 2026-05-18). Cualquier elemento **mecánico** es compatible con los 12 solvers anteriores siempre que el material declare `density`. Los elementos térmicos van por su propia vía (§2.b). **Excepción**: `CentralDifferenceSolver` requiere `lumping="lumped"` y rechaza Frame3D con eje oblicuo a los ejes globales (el bloque rotacional 3×3 del lumping no es estrictamente diagonal cuando `ρJp ≠ ρIy ≠ ρIz` — limitación documentada estándar).

---

### 1.b Elementos térmicos × materiales térmicos

Familia **paralela**, no una fila más de la tabla anterior. La compatibilidad no la determina `STRAIN_DIM` sino **`FLUX_DIM`** (2 ó 3): el campo es un escalar `T` por nodo y el gradiente `∇T` es un vector genuino, no un tensor simétrico en notación Voigt. Un material mecánico pasado a un elemento térmico se rechaza en construcción con `TypeError` explícito, y viceversa.

| Elemento | `FLUX_DIM` | `ThermalConduction` (isótropo) | `ThermalConduction` (tensorial) |
|---|---|---|---|
| `Quad4Thermal` | 2 | ✓ | ✓ |
| `Hex8Thermal`  | 3 | ✓ | ✓ |

- **(t1)** El mismo material `ThermalConduction` cubre el caso isótropo y el anisótropo: el constructor acepta un escalar `k` y lo expande a `k·I`, o un tensor completo. **La anisotropía no requiere un material nuevo** — decisión de alcance del usuario en la Etapa 8.
- **(t2)** La dimensión del material la fija el tensor con que se construyó, no la clase: `FLUX_DIM` es **propiedad de instancia**, no `ClassVar`. Construir el material con `dim=2` y pasarlo a un `Hex8Thermal` se rechaza con `ValueError` nombrando ambas dimensiones.
- **(t3)** `ρ` y `c` son **opcionales**. Sólo el análisis transitorio los exige; el estacionario funciona sin ellos, y el mensaje de error lo dice explícitamente cuando faltan.

---

### 2.b Solver térmico

| Solver | Pipeline kind | Tipo de problema | Pre-condición |
|---|---|---|---|
| `LinearSolver` | `static` | **Conducción estacionaria** `K·T = F`. | Ninguna adicional: el solver mecánico existente resuelve el térmico **sin modificación**. No requiere `ρ` ni `c`. |
| `ThetaMethodSolver` | `thermal_transient` | **Conducción transitoria** `C·Ṫ + K·T = F(t)`, integración θ de primer orden. | `ρ` y `c` declaradas; al menos un Dirichlet de temperatura. |

**El estacionario no tiene solver propio.** Que el `LinearSolver` resuelva el problema térmico sin un cambio de línea es el hallazgo arquitectural de la Etapa 8: la infraestructura de ensamblaje, imposición de Dirichlet (ADR 0004) y despacho algebraico (ADR 0003) resultó agnóstica al campo físico.

**El transitorio sí lo tiene, y no es una variante de Newmark.** La ecuación térmica es de **primer orden** en el tiempo; la familia Newmark integra la de segundo orden mediante hipótesis sobre la aceleración, que aquí no tienen sobre qué aplicarse. `ThetaMethodSolver` no admite `rayleigh` (modelo de disipación de la ecuación de segundo orden) ni condición inicial de velocidad, y devuelve `ThermalTransientResult`, no `TransientResult`.

**Defaults invertidos respecto a la dinámica estructural**, ambos por el mismo motivo físico —el principio del máximo de la ecuación de difusión, que un esquema oscilante viola produciendo temperaturas fuera del rango de los datos:

| Parámetro | Default mecánico | Default térmico |
|---|---|---|
| `lumping` de la matriz de masa/capacidad | `"consistent"` | **`"lumped"`** |
| Esquema temporal | Newmark β=1/4, γ=1/2 (sin disipación numérica) | **`θ = 1`** (Euler implícito, L-estable) |

**Combinaciones vetadas**: los solvers `ModalSolver`, `NewmarkSolver`, `HHTSolver`, `NewtonNewmarkSolver`, `NewtonHHTSolver`, `CentralDifferenceSolver`, `HarmonicSolver` y `ResponseSpectrumSolver` **no aplican** a elementos térmicos: todos integran o diagonalizan la ecuación de segundo orden. Los solvers estáticos no lineales (`NonlinearSolver`, `ArcLengthSolver`, `DissipationArcLengthSolver`) tampoco tienen uso hoy porque el único material térmico es lineal; lo tendrán cuando entre `k(T)` o radiación.

---

### 2.c Sintaxis YAML del análisis térmico

La familia térmica es accesible desde YAML como cualquier otro componente del catálogo. Tres bloques específicos, todos verificados por [`tests/test_thermal_yaml.py`](../tests/test_thermal_yaml.py):

| Bloque | Papel | Notas |
|---|---|---|
| `thermal_materials` | Paralelo a `materials`, con espacio de nombres de `type` propio | Un elemento resuelve su `material: <id>` contra el bloque donde ese id se declaró; no hay que indicar la familia |
| `thermal_loads.body_source` | Fuente volumétrica `Q` [W/m³] | `elements: [...]` la acota a un subconjunto; sin ese campo va a todo el dominio |
| `thermal_loads.boundary_flux` | Flujo prescrito `q̄` [W/m²] | `edge:` en 2D, `face:` en 3D. **`q̄ > 0` = saliente** (enfriamiento) |

Las condiciones Dirichlet y las fuentes nodales concentradas **no necesitan sintaxis nueva**: `boundary_conditions` y `point_loads` se aplican por *nombre* de DOF, así que `T: <valor>` funciona por el mecanismo genérico ya existente. Fue lo único de la cadena térmica que no hubo que cablear.

```yaml
thermal_materials:
  - {id: 1, type: ThermalConduction, k: 45.0, c: 460.0, density: 7850.0}

elements:
  - {id: 1, type: Quad4Thermal, nodes: [1,2,3,4], material: 1, thickness: 0.3}

boundary_conditions:
  - {node_id: 1, T: 100.0}

thermal_loads:
  body_source:
    - {Q: 1.0e5}
  boundary_flux:
    - {element: 1, edge: 2, q: 500.0}

solver:
  type: ThetaMethodSolver      # o LinearSolver para el estacionario
  dt: 60.0
  n_steps: 500
  T_initial: 20.0
```

**Limitación de salida**: el `VtkExporter` escribe únicamente campos mecánicos (desplazamientos y rotaciones) y **no exporta el campo de temperatura**. Es la deuda práctica más visible del subsistema térmico; el resultado se consulta hoy desde `SolveResult.U` o `ThermalTransientResult.T_history`.

**Carga variable en el tiempo**: el YAML deriva un `F_func` constante de `thermal_loads`. Para una carga o un Dirichlet que varíen en el tiempo hay que construir el solver desde código y pasar `F_func` o `dirichlet_func` propios.

---

## 3. Casos test representativos por combinación

Selección de tests "canónicos" que cubren combinaciones clave. La intención no es enumerar la suite completa (ver el recuento vigente en [STATUS.md](STATUS.md)) sino apuntar al fichero de referencia para cada celda no trivial.

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
| `Hex20` + `Elastic3D` (A.ter sub-fase 1)                       | [`test_solid_3d_higher_order.py`](../tests/test_solid_3d_higher_order.py) (patch lineal + cuadrático, simetría K, body load, face traction serendípito, masa, hourglass reducido) · [`test_cube_lame_3d.py`](../tests/validation/test_cube_lame_3d.py) (Hex20 uniaxial e hidrostática exactas) · [`test_macneal_beam_3d.py`](../tests/validation/test_macneal_beam_3d.py) (Hex20 6×1×1 alcanza 97% u_EB vs Hex8 < 55%) · [`test_volumetric_locking_3d.py`](../tests/test_volumetric_locking_3d.py) (locking atenuado vs Hex8) |
| `Hex27` + `Elastic3D` (A.ter sub-fase 2)                       | [`test_solid_3d_higher_order.py`](../tests/test_solid_3d_higher_order.py) (incluye patch test triquadrático completo `u_x=c·x²y²z²` — capacidad distintiva del Hex27 vs Hex20; conteo de 27 hourglass modes en reducción 2×2×2) · [`test_cube_lame_3d.py`](../tests/validation/test_cube_lame_3d.py) (Hex27 uniaxial e hidrostática exactas) · [`test_macneal_beam_3d.py`](../tests/validation/test_macneal_beam_3d.py) (Hex27 vs Hex20 difieren < 3% en flexión simple) · [`test_volumetric_locking_3d.py`](../tests/test_volumetric_locking_3d.py) (locking comparable a Hex20) |
| `Tet10` + `Elastic3D` (A.ter sub-fase 3)                       | [`test_solid_3d_higher_order.py`](../tests/test_solid_3d_higher_order.py) (`TestTet10Element` — patch lineal + cuadrático, body load V_e = 1/6, face traction triangular Tri6 con `tri_3`, masa consistente con `tet_15`, lumped HRZ positivo) · [`test_cube_lame_3d.py`](../tests/validation/test_cube_lame_3d.py) (malla 5-Tet10 del cubo con mid-edges compartidos; uniaxial e hidrostática exactas) |
| `Hex20`/`Hex27`/`Tet10` + materiales 3D no lineales (A.ter sub-fase 5) | [`test_solid_3d_higher_order_nonlinear.py`](../tests/test_solid_3d_higher_order_nonlinear.py) — 9 smoke tests cross-check (3 elementos × 3 materiales): `VonMises3D` plastifica (α > 0), `DruckerPrager3D` plastifica en cortante (α > 0), `IsotropicDamage3D` daña (ω > 0, κ > κ_0) bajo BCs de confinamiento uniaxial idénticas. |
| `Hex8` / `Tet4` + materiales 3D no lineales (A.bis)            | [`test_solid_3d_plasticity.py`](../tests/test_solid_3d_plasticity.py) (VM3D uniaxial libre y confinado, DP3D regular/ápice, Damage3D elástico/post-pico; smoke tests con Tet4) |
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
| `ThermalConduction` (unitario del material térmico)            | [`test_thermal_material.py`](../tests/test_thermal_material.py) — Fourier isótropo y tensorial, validación de simetría con tolerancia escalada y definición positiva por autovalores, `ρc` opcional con mensaje accionable |
| `Quad4Thermal` + `ThermalConduction` (estacionario)            | [`test_quad4_thermal.py`](../tests/test_quad4_thermal.py) — conductividad, capacidad consistente/lumped, fuente volumétrica, flujo por borde, pared plana vs perfil lineal analítico con el `LinearSolver` sin modificar |
| `Hex8Thermal` + `ThermalConduction` (estacionario + cross-check) | [`test_hex8_thermal.py`](../tests/test_hex8_thermal.py) — rango 7/8 con el modo nulo de temperatura uniforme, 4 hourglass con cuadratura reducida, flujo en las 6 caras, **cross-check 2D↔3D** contra `Quad4Thermal` |
| Transitorio térmico θ-method (orden, estabilidad, L-estabilidad) | [`test_theta_method.py`](../tests/test_theta_method.py) — orden temporal 1 / **2** / 1 para θ = 1 / 0.5 / 2/3 contra la solución exacta del sistema semidiscreto; principio del máximo; Carslaw-Jaeger semi-infinito; Dirichlet variable en el tiempo con factorización única; cross-check 2D↔3D paso a paso |
| Vía YAML del análisis térmico (bloques, cargas, despacho, no-regresión mecánica) | [`test_thermal_yaml.py`](../tests/test_thermal_yaml.py) — `thermal_materials`, `thermal_loads` con fuente y flujo contra solución analítica, despacho a `thermal_transient`, y blindaje de que un modelo mecánico no cambia de comportamiento |
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
6. ~~**Sólidos 3D**: no figuran porque no existen aún.~~ — **abierto el 2026-05-19** con la Etapa 7 (ADR 0012). **A.bis cerrado 2026-05-21**: `VonMises3D`, `DruckerPrager3D`, `IsotropicDamage3D` (✓ para `Hex8`/`Tet4`). **A.ter cerrada por completo 2026-05-27**: filas `Hex20`, `Hex27` y `Tet10` con ✓ contra todos los materiales 3D (Elastic + tres no lineales). Matriz §1 sección 3D cerrada con ✓ en las 15 celdas. **Sub-fase 4 (validación externa NAFEMS 3D) cerrada 2026-05-27**: NAFEMS LE10 thick plate pressure (Hex20+Hex27 vs σ_yy(D) canónico) y Lamé thick cylinder 3D (Hex20+Hex27 vs solución analítica cerrada en cada Gauss; primera demostración cuantitativa de capacidad isoparamétrica curva). Tet10 sobre superficie curva diferido (descomposición hex→5tets diverge; deuda técnica #8).

7. **Validación térmica sobre geometría curva** (Etapa 8): el cilindro hueco con perfil logarítmico —2D y 3D— y el balance energético global quedaron **diferidos** por requerir mallar una corona circular, capacidad que el proyecto aún no tiene. Es la misma carencia que difirió NAFEMS LE10 en su momento y que hoy bloquea también al `Tet10` sobre superficie curva (deuda técnica #8). La validación térmica vigente se apoya en pared plana analítica, sólido semi-infinito de Carslaw-Jaeger, convergencia al estacionario del `LinearSolver` y cross-check 2D↔3D. **Coste**: medio — el mallador desbloquearía tres frentes a la vez.

**No-huecos** (combinaciones que parecen ausentes pero son decisiones documentadas):

- **`DruckerPrager2D` en plane_stress**: declarado out-of-scope.
- **`Cable*` + plasticidad o daño**: extensión futura, no hueco del catálogo actual.
- **Frames + materiales 2D**: incompatibilidad estructural por `STRAIN_DIM`, no hueco.
- **Elementos térmicos × solvers dinámicos** (`ModalSolver`, `NewmarkSolver`, `HHTSolver`, `CentralDifferenceSolver`, `HarmonicSolver`, `ResponseSpectrumSolver`): incompatibilidad **física**, no hueco — todos integran o diagonalizan la ecuación de segundo orden, y la conducción es de primer orden. Ver §2.b.
- **Elementos térmicos × solvers estáticos no lineales**: sin uso *hoy* porque el único material térmico es lineal. Dejará de ser un no-hueco cuando entre `k(T)` o radiación.
- **Materiales térmicos × elementos mecánicos** (y viceversa): rechazado en construcción con `TypeError` explícito. Son familias con contratos distintos (`compute_flux` vs `compute_stress`), no una celda vacía.

---

## Cómo se mantiene este documento

- **Componente nuevo**: añadir su columna o fila; marcar celdas con ✓ donde haya test y ○ donde sea compatible sin test. Si una compatibilidad requiere restricción, añadirla como nota al pie.
- **Test nuevo que cubre celda ○**: cambiar a ✓.
- **Combinación que se descubre incompatible** (durante implementación): marcar como nota explícita, no como ✓ ni vacío.
- **Etapa cerrada**: barrer la matriz para asegurar que todas las nuevas combinaciones quedaron reflejadas.

---

*Última actualización: 2026-08-25 — **Etapa 8 (análisis térmico)**: añadidas §1.b (elementos × materiales térmicos, compatibilidad por `FLUX_DIM` en vez de `STRAIN_DIM`) y §2.b (solver térmico, con la tabla de defaults invertidos respecto a la dinámica estructural y las combinaciones vetadas por incompatibilidad física). La familia térmica va en secciones propias y no como filas de la tabla mecánica porque el campo es escalar y el gradiente no se comprime en notación Voigt. Cuatro filas nuevas en §3, hueco #7 (validación sobre geometría curva, diferida por falta de mallador de corona circular) y tres no-huecos nuevos en §4.*

*Anterior 2026-05-27 — **A.ter sub-fases 1, 2 y 3**: añade filas `Hex20`, `Hex27` y `Tet10` a la tabla 3D; centralización en `_HigherOrderSolid3D` aplicada al entrar el `Hex27` (regla de los dos casos reales) y reutilizada por `Tet10` con cara triangular Tri6 + cuadraturas `tet_4`/`tet_15` nuevas. Las columnas 3D se triplicaron con `VonMises3D`/`DruckerPrager3D`/`IsotropicDamage3D` (✓ para `Hex8`/`Tet4`). Anterior 2026-05-18: cierre Etapa 6 (ADR 0009 completo).*
