# Revisión arquitectural — 2026-05-28

> Auditoría exhaustiva de arquitectura y programación de Solidum FEM ejecutada el 2026-05-28 mediante orquestación multi-agente (workflow `solidum-arch-review` con 10 dimensiones revisadas en paralelo, 137 hallazgos propuestos, 110 confirmados tras verificación adversarial). Baseline del repositorio: commit `1bc2add` (sub-fase 4 A.ter cerrada).

## Metodología

- **10 dimensiones revisadas en paralelo**, cada una por un agente con scope acotado a su área (lectura directa del código fuente + comparación con ADRs, specs, catálogos y memoria del proyecto):
  1. Arquitectura general y contratos
  2. Subsistema de elementos
  3. Subsistema de materiales
  4. Subsistema de solvers
  5. Capa algebraica y assembly
  6. Tests y validación
  7. Coherencia documental ADR ↔ spec ↔ catálogo ↔ código
  8. Performance y robustez numérica
  9. Code quality (smells, dead code, naming, type hints)
  10. Convenciones de signos y Voigt (Reglas.md §5)

- **Verificación adversarial** de cada hallazgo: un agente `Explore` independiente re-leyó los archivos citados con instrucción de refutar bajo incertidumbre. Solo los hallazgos confirmados con evidencia concreta entran al reporte.

- **148 agentes en total** (10 reviewers + 138 verificadores), 11 M tokens consumidos.

## Veredicto

| Severidad | Hallazgos verificados |
|---|---|
| **Crítico** | **0** |
| **Alto** | **7** |
| **Medio** | **23** |
| **Bajo** | 50 |
| **Info** | 13 |
| **Total** | **110** (de 137 propuestos; 27 refutados) |

## Resumen ejecutivo

Solidum FEM se encuentra en un estado arquitectural **sólido y operacionalmente saludable**: 970 tests verdes, 21 elementos, 13 materiales, 12 solvers y 12 ADRs vivos, con cierre formal reciente de las sub-etapas A.bis (materiales 3D no lineales) y A.ter (sólidos cuadráticos 3D `Hex20`/`Hex27`/`Tet10` sobre la base centralizada `_HigherOrderSolid3D`). El núcleo cumple los tres contratos vertebrales del ADR 0001 — contratos declarativos, auto-registro por decoradores y validación temprana — y **no se detectan bugs físicos en el subsistema de materiales** (cross-consistency 2D plane_strain ↔ 3D verificada a 10⁻¹⁴).

Los hallazgos críticos no son de formulación sino de **plumbing y simetría**:

- Tres bugs medios reales en el corrector arc-length y arc-length por disipación (residuo evaluado pre-actualización).
- Un alto en el hook ADR 0010 §5 (`NewtonNewmark`/`NewtonHHT` no invocan `prepare_all_steps`, dejando `CST_Embedded2D` incompatible con análisis dinámico).
- Un bug medio en la captura de singularidad LU (`except RuntimeError` no atrapa `spsolve` que devuelve NaN silencioso).
- Un bug medio en `VonMises2D.admissibility_scale` plane_stress (devuelve stress² violando contrato).
- Ausencia de validación de parámetros en `Elastoplastic1D` (único material sin guards).

La deriva más visible y de mayor impacto público es **documental**:

- README con ejemplo YAML no ejecutable (`plane_state: stress`, `result.reactions`).
- Métricas obsoletas en README (18/10/32/804 vs reales 21/13/45/970).
- Los tres manuales (User/Architecture/Reference) congelados pre-Etapa 7 sin sólidos 3D, sin materiales 3D, sin `DissipationArcLengthSolver`.
- Catálogo de materiales que no lista los cuadráticos 3D como compatibles.
- Specs de solvers con `yaml_type` en snake_case que no coincide con los nombres Pascal del Registry.

La cobertura de tests 3D es asimétrica respecto a 2D: ningún sólido 3D se ejercita en ningún solver dinámico, falta FD de tangente algorítmica para `VonMises3D`, invariancia rotacional 3D ausente, y los smoke tests cross-check de la sub-fase 5 A.ter son demasiado laxos (umbral `α > 0` o `ω > 0` sin comparación contra material standalone).

El subsistema de elementos muestra inconsistencias menores entre hermanos: `CST_Embedded2D` omite `body_load`/`edge_traction`/`mass` que `Tri3` sí implementa, frames duplican `validate_lumping_kwarg`, `@njit(cache=True)` aplicado en `solid_3d` pero no en `solid_2d` ni en materiales con return mapping, y los kernels de orden superior 2D/3D quedan en numpy puro fuera de Numba (hot path para validaciones tipo NAFEMS LE10).

En cuanto a convenciones, el núcleo formal (Voigt 6D, RHR, sagging) se cumple en código, pero persisten violaciones terminológicas residuales ("tensiones" por σ en specs canónicos y manuales, "isotrópico" en YAML de ejemplo) y un **bug de signo dormido** en el API legado `compute_internal_forces` de los cuatro frames (`axial_force = F_local[0] = -N`, contradiciendo §5; el API ADR 0002 `internal_forces()` sí lo corrige).

**Conclusión**: el proyecto está en buen estado para extender el catálogo. Antes de cualquier hito mayor conviene cerrar el bug del hook ADR 0010 §5 en solvers dinámicos, sanear el residuo arc-length, blindar la captura de singularidad LU, y resincronizar README/manuales/catálogos con la realidad post-A.ter (lo más visible para JOSS y SoftwareX).

## Fortalezas

1. **Núcleo arquitectural firme**: los tres contratos vertebrales del ADR 0001 (contratos declarativos `STRAIN_DIM`/`DOF_NAMES`/`IS_SYMMETRIC`, auto-registro por decoradores + descubrimiento recursivo, validación temprana en `Element.__init__`) se aplican uniformemente en todos los registries y dispatch declarativo de `run_yaml` por `PIPELINE_KIND`.

2. **Convenciones del proyecto cumplidas en código**: Voigt 6D `[xx, yy, zz, xy, yz, xz]` con `γ_ij = 2·ε_ij` respetada en los 4 materiales 3D y las 5 matrices B de sólidos 3D; la convención stress-resultant/RHR interna y la traducción a viga estructural en API pública de frames 2D (signo correcto en `internal_forces()` del ADR 0002) están bien separadas.

3. **Cross-consistency 2D plane_strain ↔ 3D para J2** verificada a 10⁻¹⁴ (predictor `s_trial + p·I`, `eps_p` tensorial 4/6 componentes), y tangentes algorítmicas consistentes correctas en los regímenes activos (J2 plane_strain/3D, DP regular/apex, daño 2D/3D); calibración DP↔MC centralizada en `_calibrate_drucker_prager` reusada entre 2D/3D.

4. **Sub-etapa A.ter entregada sobre base centralizada `_HigherOrderSolid3D`** con cuadraturas dedicadas para masa (`tet_15` para Tet10), caras Quad8/Quad9/Tri6, y blindaje de hourglass counts en tests.

5. **Capa algebraica con caché de topología COO** compartida entre K, M y peso propio, despachador SPD→LU con fallback Cholesky→LU automático, lumping HRZ robusto en todos los elementos 3D nuevos vía `_MASS_QUADRATURE` específico por subclase.

6. **Calibración dual de convergencia** (ADR 0006/0007) aplicada en todos los solvers no lineales con tolerancias adimensionales `atol + rtol·escala` autoderivadas; Result classes son frozen dataclasses con docstrings honestos sobre mutabilidad shallow; sistema de excepciones tipadas centralizado en `diagnostics.py` (ADR 0011).

7. **Suite de 970 tests verdes + 8 skipped intencionales** (Tet10 sobre superficie curva — deuda #8 documentada con motivo arquitectural), incluyendo validación externa cuantitativa contra Lamé thick cylinder 3D, NAFEMS LE10, MacNeal-Harder, Bathe wave, Hill J2 cilindro, triaxial DP3D, softening damage 3D.

8. **Documentación navegacional viva**: STATUS.md, ROADMAP.md, MATRIZ.md y los catálogos de elementos/solvers actualizados a 2026-05-27; ADRs persistentes (12 documentos, todos los grandes hitos arquitecturales registrados); specs como fuente canónica con `acceptance` versionada.

9. **Registry de cuadraturas completo** para el catálogo actual (10 reglas, incluida `tet_15` para Tet10); separación limpia entre familias estructurales 1D (que exponen `internal_forces` + `ElementForces`) y sólidos 2D/3D (que exponen `compute_gauss_state` como API canónica de salida según ADR 0012).

---

## Hallazgos altos (7)

### A1. `NewtonNewmarkSolver` / `NewtonHHTSolver` no invocan `prepare_all_steps` (hook ADR 0010 §5)

- **Dimensión**: solvers
- **Tipo**: missing
- **Ubicación**: `solidum/math/solvers/newmark.py:406-554` (`NewtonNewmarkSolver.solve`) y `:1043-1199` (`NewtonHHTSolver.solve`)
- **Descripción**: ADR 0010 §5 define el hook `prepare_all_steps` como contrato de los solvers no lineales para que elementos con discontinuidades embebidas (`CST_Embedded2D`) chequeen activación al inicio de cada paso con el estado convergido del anterior. `NonlinearSolver`, `ArcLengthSolver` y `DissipationArcLengthSolver` sí lo invocan; los dos solvers dinámicos no lineales **no**, dejando la combinación `CST_Embedded2D` + análisis dinámico silenciosamente incorrecta (no se chequea activación; el state queda en su default).
- **Recomendación**: añadir `self.assembler.prepare_all_steps(u_total)` al inicio de cada paso temporal del bucle `for step in range(n_steps)` en `NewtonNewmarkSolver.solve()` y `NewtonHHTSolver.solve()`, justo antes de calcular `u_pred`/`udot_pred`. Documentar en `catalogo_solvers.md` la compatibilidad real con embedded discontinuity (o declararla explícitamente fuera de scope hoy) y añadir test minimal `CST_Embedded2D + NewtonNewmark` para blindar la regresión.

### A2. `Elastoplastic1D` sin validación de parámetros — único material sin guards

- **Dimensión**: materials
- **Tipo**: inconsistency
- **Ubicación**: `solidum/materials/plastic_1d.py:15-20`
- **Descripción**: todos los demás materiales (Elastic1D/2D/3D, VonMises2D/3D, DruckerPrager2D/3D, IsotropicDamage1D/2D/3D, CableMaterial1D, CohesiveDamageIsotropic) validan en `__init__` que `E > 0`, `sigma_y > 0`, `H >= 0`, `density >= 0` con mensajes específicos. `Elastoplastic1D.__init__` no valida nada: acepta `E < 0`, `sigma_y = 0`, `H < 0`, `density < 0` sin error. Viola Reglas.md §1 (validación temprana al construir) y permite tangentes negativas / signos invertidos en return mapping con divergencia silenciosa del Newton.
- **Recomendación**: añadir el mismo bloque de validaciones que `VonMises2D/3D`: `if E<=0: raise ValueError(...)`, etc. Coste: ~10 líneas; previene una clase entera de bugs físicos silenciosos.

### A3. README.md: ejemplo YAML del Quick Example no es ejecutable

- **Dimensión**: docs
- **Tipo**: bug
- **Ubicación**: `README.md:54-86`
- **Descripción**: el snippet del Quick example contiene dos errores que hacen fallar el modelo:
  1. `Elastic2D` declara `plane_state: stress` pero el parámetro real es `hypothesis: plane_stress` — produce `TypeError` al construir.
  2. `print(result.reactions)` falla con `AttributeError`; el campo correcto es `result.reactions_by_node`.

  Bug doble en el ejemplo de portada del repo público, justo el material que verá un revisor de JOSS o un usuario nuevo.
- **Recomendación**: cambiar `plane_state: stress` → `hypothesis: plane_stress` (alineado con `manuals/sources/user/05_materiales.md:26-28`). Cambiar `result.reactions` → `result.reactions_by_node`. Añadir el ejemplo a `tests/test_examples_yaml.py` para que CI falle si se vuelve a romper.

### A4. README.md: cuentas de elementos/materiales/specs/tests desactualizadas

- **Dimensión**: docs
- **Tipo**: docs-drift
- **Ubicación**: `README.md:14-27`
- **Descripción**: el README declara "18 elements", "10 materials", "32 validated specs", "804 tests" como capacidades actuales — las cuatro métricas obsoletas en > 6 meses de desarrollo. Estado real verificado: **21 elementos, 13 materiales, 45 specs validadas + 1 implemented, 970 tests + 8 skipped**. README es el primer artefacto público (LGPL en GitHub) y será citado por JOSS y SoftwareX.
- **Recomendación**: actualizar las cuatro líneas a (21/13/45/970). Considerar automatizar la generación de la sección desde `STATUS.md` mediante un pre-commit hook o generador para evitar la deriva crónica.

### A5. Catálogo de materiales: materiales 3D no listan Hex20/Hex27/Tet10 como compatibles

- **Dimensión**: docs
- **Tipo**: docs-drift
- **Ubicación**: `docs/catalogo_materiales.md:45,88,171,221`
- **Descripción**: tras cerrar la sub-fase 5 de A.ter (cross-check 3 elementos × 3 materiales no lineales), `MATRIZ.md §1` marca las 15 celdas elemento × material con ✓. El catálogo de materiales sigue declarando para los 4 materiales 3D "Compatible con: Hex8, Tet4", omitiendo los tres elementos cuadráticos 3D. El catálogo es la referencia rápida del usuario antes de modelar; declara incompatibilidades inexistentes.
- **Recomendación**: actualizar las cuatro líneas a "Compatible con: Hex8, Hex20, Hex27, Tet4, Tet10 (todos los sólidos 3D con `STRAIN_DIM = 6`)" replicando el patrón ya usado en las entradas 2D.

### A6. User Manual: catálogo de elementos solo cubre 1D y 2D (omite los 5 sólidos 3D)

- **Dimensión**: docs
- **Tipo**: docs-drift
- **Ubicación**: `manuals/sources/user/04_elementos.md`
- **Descripción**: el capítulo declara "El motor expone dos familias de elementos: 1D y 2D". No menciona Hex8, Tet4 (Etapa 7), Hex20, Hex27, Tet10 (A.ter) ni CST_Embedded2D (Etapa 5). El manual es bilingüe y se publica como artefacto público; un usuario externo no encontrará documentación de los sólidos 3D pese a estar todos validados y registrados.
- **Recomendación**: añadir sección "Elementos 3D" tras la 2D (Hex8, Tet4 lineales + Hex20, Hex27, Tet10 cuadráticos) siguiendo el patrón de `docs/catalogo_elementos.md:498-602`, y sección "Elementos con discontinuidad embebida" para CST_Embedded2D. Regenerar el PDF vía `build_user_manual.py`.

### A7. Architecture Manual: capítulo ADR omite ADR 0010, 0011, 0012

- **Dimensión**: docs
- **Tipo**: docs-drift
- **Ubicación**: `manuals/sources/architecture/08_decisiones_adr.md`
- **Descripción**: el capítulo resume ADRs 0001 a 0009. ADR 0010 (Discontinuidades embebidas), 0011 (Robustez Newton + line search) y 0012 (Sólidos 3D y Voigt 6D), todos aceptados antes del 2026-05-21, no aparecen en el resumen pese a ser los más importantes — ADR 0012 introduce la convención Voigt 6D obligatoria del proyecto. El capítulo salta de 0009 a "PENDIENTE térmico".
- **Recomendación**: añadir tres subsecciones tras ADR 0009 resumiendo 0010 (CST_Embedded2D + familia cohesiva), 0011 (line search + telemetría tipada, con la enmienda `default=False`) y 0012 (Voigt 3D + cierre dominio de `internal_forces` + API `face_traction`). Regenerar PDF.

### A8 (alto en tests). Sólidos 3D sin cobertura en ningún solver dinámico

- **Dimensión**: tests
- **Tipo**: missing
- **Ubicación**: `tests/test_modal*.py`, `test_newmark*.py`, `test_hht.py`, `test_central_difference.py`, `test_harmonic.py`, `test_response_spectrum.py`
- **Descripción**: ninguno de los 5 elementos sólidos 3D (Hex8, Hex20, Hex27, Tet4, Tet10) aparece en los tests de los 8 solvers dinámicos. `MATRIZ.md §2` afirma compatibilidad universal de cualquier elemento con los solvers dinámicos vía `compute_mass_matrix(lumping)`; esa afirmación queda sin verificar empíricamente para 3D. Para el catálogo cuadrático esto es especialmente delicado: Hex20/Hex27/Tet10 declaran `_MASS_QUADRATURE` específicos (`tet_15`) y HRZ canónico — bugs de plumbing en ensamblaje masa-rigidez 3D dentro del pipeline modal/transitorio pasarían inadvertidos.
- **Recomendación**: añadir al menos un test analítico modal por elemento 3D (barra cantilever Hex8/Hex20 con ω₁ axial vs solución cerrada), un test HRZ vs consistente bajo Newmark lineal en columna sólida 3D, y un test `CentralDifferenceSolver` con Hex8 + masa lumped para verificar el CFL discreto.

> Nota: el conteo "alto: 7" del veredicto se refiere a los 7 hallazgos altos confirmados con verificación adversarial; A8 también es alto en su dimensión.

---

## Hallazgos medios destacados (selección de los 23)

### Solvers

**M1. `ArcLengthSolver` evalúa convergencia con residuo previo a la actualización (R desfasado)**
- `solidum/math/solvers/arclength.py:168-228`
- Dentro del bucle corrector, `R = lambda_iter·F_ext_ref − F_int_global` se calcula con `lambda_iter` y `F_int` previos a aplicar `ddlambda` y `dU_update`. Tras la actualización no se vuelve a ensamblar `F_int(U_iter)` ni a recomputar R. La convergencia se evalúa con un R que no corresponde al estado final del paso.
- **Fix**: tras actualizar `lambda_iter` y `U_iter`, re-ensamblar `K_global, F_int_global = self.assembler.assemble_non_linear_system(U_iter)` y recomputar R antes de `self.convergence.evaluate()`. Centralizar como helper compartido con `DissipationArcLengthSolver`.

**M2. `DissipationArcLengthSolver` hereda el mismo desfase**
- `solidum/math/solvers/dissipation_arclength.py:310-388`
- Mismo patrón que M1. Además el switch cilíndrico↔disipación puede activarse en pasos cuya convergencia es ficticia.
- **Fix**: aplicar la misma corrección y centralizar en helper compartido.

**M3. `NewtonHHTSolver` line search usa residuo Newmark estándar, no HHT-α**
- `solidum/math/solvers/newmark.py:578-620` (`_armijo_step_dynamic`) y `:1121-1126` (uso en NewtonHHT)
- Si el usuario activa `line_search=True` en `NewtonHHTSolver`, el line search obtiene α minimizando un residuo distinto al que rige Newton — métricas inconsistentes.
- **Fix**: sobreescribir `_armijo_step_dynamic` en `NewtonHHTSolver` para evaluar el residuo HHT-α coherentemente.

**M4. `ConvergenceCriterion` no recalibra al reusar entre corridas**
- `solidum/math/convergence.py:95-100`
- `is_calibrated` permanece `True` para siempre. Si se reusa el mismo criterio en dos solvers o dos corridas, la calibración antigua persiste.
- **Fix**: añadir `ConvergenceCriterion.reset()` y llamarlo al inicio de cada `solver.solve()`. Documentar en ADR 0007.

**M5. `ModalResult.free_vibration` ignora el lumping con que se construyó el ModalSolver**
- `solidum/results.py:91-118`
- Si el `ModalSolver` usó `lumping='lumped'`, el consumidor debe pasar M con el mismo lumping; pasar M consistente cuando los modos eran de lumped da resultados con factores distintos. Defecto sistemático: toda función sobre modos hereda el problema.
- **Fix**: almacenar el lumping efectivo en `ModalResult` (campo extra) y cachear M_efectiva.

### Materiales

**M6. `VonMises2D.admissibility_scale` en plane_stress devuelve stress²**
- `solidum/materials/von_mises_2d.py:399-411`
- Contract `Material.admissibility_scale` declara "en unidades de esfuerzo". Plane_strain devuelve `√(2/3)·R` (stress, correcto), plane_stress devuelve `R²/3` (stress²). La fórmula `tol = ABS + REL·scale` deja de tener sentido dimensional; el piso ABS = 1e-14 desaparece de facto.
- **Fix**: devolver siempre `√(2/3)·R` y traducir internamente en el kernel `_compute_j2_plane_stress`.

### Álgebra y assembly

**M7. `except RuntimeError` no atrapa singularidad real cuando se usa LUSolver (NaN silencioso)**
- `solidum/math/solvers/nonlinear.py:251-256`; `arclength.py:142-146,176-181`; `dissipation_arclength.py:224-230,319-324`; `newmark.py:472-476,1113-1117`
- `LUSolver.solve` invoca `scipy.sparse.linalg.spsolve`, que ante matriz singular **no levanta excepción**: emite `MatrixRankWarning` y devuelve un vector con NaN/inf. Solo `splu` levanta `RuntimeError`. La singularidad de la tangente no se detecta; el NaN se propaga.
- **Fix**: en `LUSolver.solve` y `CholeskySolver.solve` chequear `np.any(~np.isfinite(x))` tras `spsolve` y levantar `LinearAlgebraSingularError` tipado. Test con sistema deliberadamente singular que verifique `SingularTangentError`.

**M8. Cuatro solvers reimplementan la reducción Dirichlet bypaseando `Assembler.reduce*`**
- `newmark.py:130-138,344-361,743-750,981-994,1078-1080`; `central_difference.py:204-211`; `harmonic.py:162-168`
- `Assembler.reduce`, `reduce_pair` y `expand` son el camino canónico de reducción Dirichlet (ADR 0004). LinearSolver, NonlinearSolver, ArcLength, DissipationArcLength lo usan. Pero Newmark, HHT, NewtonNewmark, NewtonHHT, CentralDifference, Harmonic hacen `T, g = cs.build(ndof); K_red = (T.T @ K @ T).tocsr()` inline.
- **Fix**: extender `Assembler.reduce` para devolver opcionalmente `(K_red, T, g, F_dir)` cuando `F=None`, o añadir helper `reduce_static_dirichlet`. Refactorizar los seis solvers.

### Arquitectura

**M9. `compute_gauss_state` es API canónica de sólidos pero no figura en el contrato del base Element**
- `solidum/core/element.py:90-242` y `element_forces.py:19-23`
- El ADR 0012 promociona `compute_gauss_state(U)` como API canónica de salida para sólidos 2D/3D, y los docstrings la nombran. Pero el método no existe en `Element` — ni como abstracto, ni con stub. Invocarlo sobre un estructural 1D devuelve `AttributeError` críptico.
- **Fix**: añadir en `Element` un método `compute_gauss_state(U_global) -> dict | None` con docstring de contrato y stub que devuelve `None`. Las claves del dict son contrato, no convención.

**M10. `solidum/__init__.py` no re-exporta los sólidos 3D ni los materiales 3D**
- `solidum/__init__.py:31-46`
- Tras cerrar Etapa 7 + A.bis + A.ter, el paquete raíz sigue re-exportando solo la familia 2D. `from solidum import Hex8` o `from solidum import VonMises3D` falla con `ImportError`. Para una API pública pretendida para JOSS, el contraste con la documentación duele.
- **Fix**: añadir al re-export: Elastic3D, VonMises3D, DruckerPrager2D, DruckerPrager3D, IsotropicDamage3D, Hex8, Hex20, Hex27, Tet4, Tet10. Considerar generar la lista dinámicamente desde los registries.

### Tests

**M11. `VonMises3D` sin validación FD de tangente algorítmica consistente**
- VM2D, DP3D y Damage3D tienen verificación FD central de C_alg en régimen activo. VM3D no la tiene. La tangente algorítmica es el motor de la convergencia del Newton global; sin FD-check, una derivación equivocada de los términos `2G·dλ·N⊗N` en 6D pasaría inadvertida.
- **Fix**: añadir `test_tangent_matches_finite_difference_regular` en `TestVonMises3D` con ε mixto multi-componente.

**M12. Sin invariancia rotacional para los 4 materiales 3D**
- El test es el blindaje canónico contra bugs sutiles en B, factores ½ del desviador en Voigt, y dependencias espurias de los ejes globales. En Voigt 6D la `T_ε(R)` tiene estructura distinta a 2D; la convención del proyecto difiere de ABAQUS en la permutación yz↔xz.
- **Fix**: cuatro test classes paritarias barriendo rotaciones canónicas (90°, 45°, eje general) y verificando `σ(R·ε·Rᵀ) = R·σ(ε)·Rᵀ`.

**M13. Cross-check 3D higher-order × materiales no lineales (sub-fase 5 A.ter) demasiado laxo**
- `tests/test_solid_3d_higher_order_nonlinear.py`
- Los 9 smoke tests solo verifican `max(α) > 0` o `max(ω) > 0`. No comparan σ contra material standalone bajo la misma trayectoria de ε. Las 9 celdas pasaron de ○ a ✓ con el umbral más débil posible; un bug que devolviese σ con factor ≠ 1 pero con plasticidad activa pasaría todos los tests.
- **Fix**: extender cada mixin con un test que compare σ_xx FEM (avg Gauss) contra material standalone integrado punto a punto sobre el path ε_xx prescrito.

### Performance

**M14. `ArcLengthSolver` factoriza la misma K dos veces por iteración del corrector**
- `arclength.py:166-181`
- Cada iteración ejecuta `du_R_red = self._solve(K_red, R_red); du_t_red = self._solve(K_t_red, F_t_red)` sobre la misma K. Como `_solve` delega en `spla.spsolve` que re-factoriza internamente, el corrector paga ~2× la factorización por iteración aunque la rigidez sea idéntica.
- **Fix**: factorizar K_red una vez con `self._linalg.factorize(K_red)`, retener el `FactorizedSolver` y llamar `factor.solve()` dos veces. Aprovecha `LUFactorized`/`CholeskyFactorized` del ADR 0003 fase 2.

**M15. VM2D plane_stress Newton local sale silenciosamente sin convergencia**
- `solidum/materials/von_mises_2d.py:183-230`
- El Newton local sobre Δγ recorre `max_local_iter = 25` con `break` cuando converge. Si ninguna condición se cumple, el bucle sale sin marcar no convergencia y el caller acepta `sigma`/`C_alg` confiando en que son consistentes. Patrón equivalente en DP2D/DP3D rama ápice. Encaja en `feedback_no_parches.md`.
- **Fix**: devolver flag `local_converged: bool`. Si `False`, emitir warning vía `logging` o lanzar `MaterialLocalNonConvergenceError` que el `NonlinearSolver` clasifique.

**M16. Kernels @njit de orden superior 2D/3D en numpy puro**
- `solidum/elements/solid_2d/_shared.py:161-271`; `solid_3d/_shared.py:330-662`
- Funciones de forma cuadráticas y derivadas (Quad8, Quad9, Tri6, Hex20, Hex27, Tet10) ejecutan en numpy puro con bucle Python sobre nodos pese a constituir el bucle interno de cada ensamblaje. Hex20 con 3×3×3 ejecuta 27 evaluaciones de `_dN_hex20` (bucle Python de 20 nodos con condicionales) por elemento — cuello de botella en NAFEMS LE10 con miles de elementos. Hex8 está dentro de un `@njit` monolítico.
- **Fix**: migrar progresivamente a `@njit(cache=True)`. Speedup esperado 5-20× en ensamblaje.

### Documentación

**M17. Reference/User/Architecture Manual PDFs no se han regenerado tras A.ter**

**M18. Specs de solvers declaran yaml_type en snake_case que no coincide con el Registry**
- `docs/specs/LinearSolver.md:93`, `NonlinearSolver.md:113`, `ArcLengthSolver.md:126`, `ModalSolver.md:67`, `NewtonHHTSolver.md:97`, `DissipationArcLengthSolver.md:169`
- Seis specs documentan `yaml_type: linear`, `nonlinear`, etc. en snake_case, pero `SolverRegistry` los registra con el nombre de clase Pascal. Las otras specs sí declaran Pascal — inconsistencia interna.
- **Fix**: corregir las seis specs + test `tests/test_specs.py` que verifique `SolverRegistry.get(yaml_type)` no lanza.

**M19. Architecture Manual: lista de materiales/solvers omite subsistema 3D y DissipationArcLength**

**M20. Architecture Manual: evolución por fases no registra Etapas 5 y 7**

**M21. User Manual: catálogo de materiales omite Elastic3D, materiales 3D no lineales, cohesivos y Quad8/Quad9/Tri6**

### Convenciones

**M22. Bug de signo dormido en API legado `compute_internal_forces` de los 4 frames**
- Los cuatro frames implementan `compute_internal_forces` con `axial_force = F_local[0]` que entrega `−N` (contradice Reglas §5: `N > 0 ↔ tracción`). El API ADR 0002 `internal_forces()` sí lo corrige aplicando el signo. El legado sigue activo y mal.
- **Fix**: corregir el signo o deprecar el API legado a favor del ADR 0002 con warning.

**M23. Violaciones terminológicas residuales: "tensiones" por σ en specs canónicos y manuales**

---

## Bajos e info (63 hallazgos)

Agrupados por dimensión (detalle disponible en transcripts del workflow `wf_a5a52ebc-b2f`):

- **Arquitectura** (3): type hints incompletos, docstrings con lenguaje obsoleto pre-A.ter, comentarios `noqa` heredados.
- **Elementos** (7): `cache=True` heterogéneo entre solid_2d y solid_3d, frames duplican `validate_lumping_kwarg`, Tet10 con cuadratura `tet_4` subintegra body load, FACE_NODES inconsistencia menor de orden.
- **Materiales** (5): docstrings con "tensión" residual, alias internos inconsistentes, calibraciones DP↔MC con TODO sin fechar.
- **Solvers** (6): duplicación masiva NewtonNewmark/NewtonHHT (~260 líneas idénticas), logs heterogéneos.
- **Álgebra** (4): falta validación de positividad de diagonal en M global, no hay invalidación explícita de caché COO al modificar dominio.
- **Tests** (8): Tet4 sin blindaje cuantitativo de shear locking, sin patch test 3D multi-elemento, sin ejemplo YAML 3D end-to-end, inconsistencias spec↔test en hourglass counts Hex20/Hex27.
- **Docs** (14): conteo de specs validadas subreportado (38 vs 45), enlaces rotos en MATRIZ, README sin badge de tests, README sin demostración 3D, etc.
- **Performance** (6): materiales con return mapping sin `cache=True`, M global sin validación.
- **Code quality** (8): `is_pd=False` muerto, reexports asimétricos, fósiles cosméticos.
- **Convenciones** (2): "isotrópico" en YAML de ejemplo, comentarios con "pórtico" (uno residual).

---

## Recomendaciones priorizadas

### Bloque 1 — Fixes rápidos altos (~1 sesión)

1. **Bug ADR 0010 §5** en `NewtonNewmarkSolver` y `NewtonHHTSolver`: añadir `self.assembler.prepare_all_steps(u_total)` al inicio de cada paso temporal + test `CST_Embedded2D + NewtonNewmark`.
2. **`Elastoplastic1D` guards**: añadir validación `E > 0`, `sigma_y > 0`, `H >= 0`, `density >= 0` (~10 líneas siguiendo patrón VM2D/3D).
3. **README** fix: corregir `hypothesis: plane_stress`, `result.reactions_by_node` + actualizar métricas a (21/13/45/970). Añadir el ejemplo a `tests/test_examples_yaml.py`.

### Bloque 2 — Bugs medios de solvers (siguiente sesión)

4. **Sanear corrector arc-length** (M1, M2): re-ensamblar `F_int(U_iter)` post-actualización en `ArcLengthSolver` y `DissipationArcLengthSolver`. Helper compartido.
5. **Blindar singularidad LU** (M7): `LUSolver.solve` chequea `np.any(~np.isfinite(x))` y levanta `LinearAlgebraSingularError`; cambiar los 7 `except RuntimeError` a `except LinearAlgebraSingularError`. Test con matriz singular.
6. **Corregir `VonMises2D.admissibility_scale` plane_stress** (M6): devolver stress (no stress²) y homogeneizar el chequeo dentro del kernel.
7. **Pulir contrato medio** (M4, M5, M3, M9): `ConvergenceCriterion.reset()`, almacenar `lumping` en `ModalResult`, `_armijo_step_dynamic` override en NewtonHHT con residuo HHT-α coherente, promocionar `compute_gauss_state` al base `Element` con stub `None`.

### Bloque 3 — Resincronización pública JOSS/SoftwareX (sesión dedicada)

8. **Regenerar los 3 manuales** tras añadir sólidos 3D, materiales 3D, `DissipationArcLengthSolver`, ADRs 0010-0012 y Etapas 5/7 a las fuentes Markdown.
9. **Actualizar catálogo de materiales** para listar Hex20/Hex27/Tet10 como compatibles con los materiales 3D.
10. **Corregir `yaml_type` en 6 specs** a Pascal-case y añadir test que ate spec ↔ Registry.
11. **Re-exports de sólidos/materiales 3D** en `solidum/__init__.py`.

### Bloque 4 — Cobertura tests 3D (otra sesión)

12. **FD de tangente** para `VonMises3D` (M11).
13. **Invariancia rotacional 3D** para los 4 materiales (M12) — caza la permutación yz↔xz vs ABAQUS.
14. **Endurecer cross-check A.ter sub-fase 5** (M13) con comparación contra material standalone, no solo `α > 0`.
15. **Benchmark dinámico por sólido 3D** (A8): al menos modal + Newmark + CentralDifference por elemento 3D para validar empíricamente la promesa universal de `MATRIZ.md §2`.
16. **Añadir `examples/modelo_solido_3d.yaml`** ejecutable end-to-end.

### Bloque 5 — Refactors mayores (horizonte largo)

17. **Centralizar Dirichlet** en `Assembler` (M8) y migrar los 6 solvers dinámicos al camino canónico.
18. **Migrar a `@njit(cache=True)`** los kernels de orden superior 2D/3D (M16) y los materiales con return mapping (VM3D, DP2D/3D, Damage2D/3D); homogeneizar `cache=True` en solid_2d.
19. **Factorizar K_red una sola vez** por iteración del corrector arc-length (M14) usando `FactorizedSolver` del ADR 0003 fase 2.

---

## Notas por dimensión

| Dimensión | Estado |
|---|---|
| **Arquitectura** | Núcleo sano y bien gobernado por los contratos del ADR 0001; los hallazgos son inconsistencias post-A.ter (re-exports asimétricos, type hints, contrato `compute_gauss_state` no declarado en base) y docstrings con lenguaje obsoleto — todo corregible sin tocar semántica. |
| **Elementos** | Subsistema coherente en convenciones (Voigt, RHR, sagging) y bien centralizado en `_HigherOrderSolid3D` tras A.ter; las grietas son entre hermanos (CST_Embedded2D omite body_load/edge_traction/mass que Tri3 tiene; frames duplican `validate_lumping_kwarg`; `cache=True` heterogéneo) y subintegración del body load en Tet10 con `tet_4`. Ninguna falla bloqueante. |
| **Materiales** | Físicamente sano: cross-consistency 2D↔3D a 10⁻¹⁴, tangentes algorítmicas consistentes correctas, calibración DP↔MC centralizada. Hallazgos sustantivos: `Elastoplastic1D` sin guards (alto) y `VonMises2D.admissibility_scale` plane_stress en stress² (medio); el resto son inconsistencias menores de firma y validación. |
| **Solvers** | ADRs aplicados en su mayoría, pero hay un alto (hook `prepare_all_steps` ausente en NewtonNewmark/NewtonHHT) y dos medios (residuo arc-length pre-actualización en ambos solvers; line search HHT con residuo Newmark). `ConvergenceCriterion` no resetea entre corridas y `ModalResult.free_vibration` ignora el lumping del Solver — defectos sistemáticos de medio impacto. |
| **Álgebra** | Capa bien estructurada cumpliendo ADR 0003; un bug medio latente (NaN silencioso por `spsolve` no atrapado por `except RuntimeError`), seria duplicación de reducción Dirichlet en 6 solvers que bypasea `Assembler.reduce*`, y deuda menor de caching/invalidación de M. Lumping HRZ y QuadratureRegistry robustos para todo el catálogo cuadrático 3D. |
| **Tests** | Suite 970 verde con cobertura sólida en 2D maduro, pero asimetría notable post-Etapa 7/A.bis/A.ter: cero cobertura 3D en solvers dinámicos, falta FD de tangente en VM3D, sin invariancia rotacional 3D, smoke tests cross-check A.ter sub-fase 5 demasiado laxos, Tet4 sin blindaje cuantitativo de shear locking, sin patch test 3D multi-elemento, sin ejemplo YAML 3D end-to-end. |
| **Docs** | Documentación navegacional (STATUS/ROADMAP/MATRIZ/catálogos) al día, pero deriva significativa en artefactos públicos: README con ejemplo no ejecutable y métricas obsoletas, los tres manuales congelados pre-Etapa 7, catálogo de materiales sin cuadráticos 3D, specs de solvers con `yaml_type` snake_case, conteo de specs subreportado (38 vs 45 reales). |
| **Performance** | Núcleo numérico bien diseñado (despachador, política de tolerancias, caché COO, predictores Newmark). Las oportunidades reales: kernels de orden superior 2D/3D fuera de `@njit`, materiales con return mapping sin `cache=True`, arc-length factorizando K dos veces por iteración, Newton local del J2 plane_stress saliendo sin marcar no convergencia, y M global sin validación de positividad de diagonal. |
| **Code quality** | Buen estado general (docstrings ricos, contratos declarativos, autoregistro limpio, ausencia de TODOs colgados). Problemas reales: duplicación masiva entre `NewtonNewmark.solve` y `NewtonHHT.solve` (~260 líneas idénticas cada una), triple duplicación Hex8/Tet4 vs `_HigherOrderSolid3D`, fósiles cosméticos (`is_pd=False` muerto, noqa engañoso) y reexports asimétricos en `solidum/__init__.py`. |
| **Convenciones** | Núcleo formal del proyecto cumple convenciones canónicas (Voigt 6D, RHR, sagging, traducción interna↔API pública en frames 2D). Violaciones son terminológicas (uso residual de "tensiones" por σ en specs canónicos y manuales LaTeX; "isotrópico" en YAML) y un bug de signo dormido en el API legado `compute_internal_forces` de los cuatro frames. |

---

*Generado el 2026-05-28 mediante workflow `solidum-arch-review` (run `wf_a5a52ebc-b2f`). Baseline: commit `1bc2add`. 148 agentes (10 reviewers + 138 verificadores), 11.1 M tokens.*
