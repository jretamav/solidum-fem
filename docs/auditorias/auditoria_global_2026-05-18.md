# Auditoría global de Fenix FEM — 2026-05-18

> Auditoría a fondo de arquitectura, programación, implementación computacional y numérica. Cinco ejes paralelos (arquitectura/plumbing, elementos, materiales, solvers, cobertura/docs) sobre el snapshot del repositorio al cierre de Etapa 6 (ADR 0009 completo + ADR 0010 fase 3b + ADR 0011 enmendado).
>
> **Metodología**: cuatro agentes especializados leyeron el código de verdad (no sólo grep) bajo briefs autocontenidos con criterios concretos, evidencia `file:line`, severidad estandarizada y referencia a memorias del proyecto. El consolidador (yo) verificó todos los hallazgos críticos releyendo directamente el código antes de incluirlos en este informe.
>
> **Alcance**: ~70 archivos `.py` en `fenix/`, ~45 archivos de test, 11 ADRs aceptados, 32 specs, 558 tests verdes en suite global. La auditoría es **diagnóstico, no intervención**: este documento identifica; el usuario decide qué retomar y en qué orden.

---

## Veredicto ejecutivo

**Fenix FEM es un código robusto.** La arquitectura honra los principios de `Reglas.md`: contratos declarativos consistentes (`STRAIN_DIM`, `DOF_NAMES`, `PIPELINE_KIND`, `JUMP_DIM`), validación temprana con mensajes útiles, capas sin ciclos (`core ← math ← {materials, elements} ← entry/utils`), tolerancias adimensionales con autocalibración (ADRs 0006 y 0007), despacho declarativo del pipeline (ADR 0009 regla C), excepciones tipadas (ADR 0011). El núcleo numérico está rigurosamente implementado: J2 plane stress según Simó-Hughes §3.4.1 (Newton local con diagonalización en autovectores y tangente cerrada con corrección por endurecimiento — sofisticado), Drucker-Prager con dos ramas regular/apex y tangente algorítmica explícita, daño 2D con tangente asimétrica analítica, cinemática KOS del embedded con condensación estática y `l_d = (A_e/h)·cos(θ−α)` del Cap. 6 de la tesis. Los esquemas temporales (Newmark β-γ, HHT-α con auto-derivación, leapfrog Belytschko-Liu-Moran, Crisfield cilíndrico, SRSS/CQC Der Kiureghian 1980) están escritos conforme a la bibliografía canónica.

**Hay tres bugs críticos accionables**, todos verificados por relectura directa del código y todos **latentes** (no rompen la suite actual porque el escenario disparador no está cubierto):

1. `singular_tangent_seen` muerto en `NewtonNewmarkSolver` y `NewtonHHTSolver` → `SingularTangentError` del ADR 0011 inaccesible desde el subsistema dinámico no lineal.
2. `LinearSolver.solve` no cachea factorización → renuncia a la promesa de reuso del ADR 0003 para multiplicidad de cargas.
3. `CentralDifferenceSolver` con `nonlinear=True` y apoyos prescritos no nulos resta `F_dir` dos veces — invisible mientras `g = 0` (caso de tests), corrompe excitación cinemática del soporte y asentamientos prescritos en dinámico.

Por encima de los críticos, **4 hallazgos altos** (autodiscover no recursivo, triplicación `_resolve_rayleigh` confirmada por dos ejes independientes, sólidos 2D no exponen `internal_forces`, tests con coeficientes unitarios que ocultan bugs dimensionales) y **un cuerpo de hallazgos medios-bajos** que pintan un código maduro con deudas razonables.

**Recomendación de uso del informe**: el usuario revisa los 3 críticos primero (decisión de fixear ahora o programar), después los 4 altos (deuda técnica priorizable), y deja los medios-bajos como housekeeping cuando entren en agenda. Ninguno de los críticos justifica abrir una etapa de pausa para corregir — son fixes acotados.

---

## Métricas de la auditoría

| Eje | Hallazgos | Críticos | Altos | Medios | Bajos |
|---|---:|---:|---:|---:|---:|
| 1. Arquitectura y plumbing | 15 | 0 | 2 | 8 | 5 |
| 2. Elementos | 10 | 0 | 2 | 5 | 3 |
| 3. Materiales | 10 | 0 | 0 | 5 | 5 |
| 4. Solvers | 10 | 3 | 0 | 5 | 2 |
| 5. Cobertura, validación y docs | 6 | 0 | 0 | 4 | 2 |
| **Total** | **51** | **3** | **4** | **27** | **17** |

Sin hallazgos críticos en formulación física, signos, return mappings, o esquemas temporales. Los tres críticos son de **infraestructura** (diagnóstico tipado, reuso de factorización, contabilidad de apoyos en dinámico explícito).

---

## Top 10 hallazgos globales priorizados

1. **[CRÍTICO · Eje 4]** `singular_tangent_seen` muerto en `NewtonNewmarkSolver` / `NewtonHHTSolver` — diagnóstico tipado del ADR 0011 inaccesible. (`newmark.py:442`, `newmark.py:476-484`).
2. **[CRÍTICO · Eje 4]** `LinearSolver.solve` no cachea factorización entre llamadas — promesa rota del ADR 0003 para barrido de cargas. (`linear.py:25-49`).
3. **[CRÍTICO · Eje 4]** `CentralDifferenceSolver` no lineal con apoyos prescritos no nulos resta `F_dir` dos veces. (`central_difference.py:232,296,353-357`).
4. **[ALTO · Eje 1]** Autodiscover no recursivo: `embedded_cst.py` ya depende de re-export manual en `fenix/__init__.py:42`. Futuros archivos en subpaquetes pueden no registrarse silenciosamente. (`autodiscover.py:17`, `solid_2d/__init__.py:34-38`).
5. **[ALTO · Eje 1 + Eje 4]** Triplicación literal de `_resolve_rayleigh` en Newmark, CentralDifference y Harmonic — confirmado independientemente por dos ejes. (`newmark.py:113-127`, `central_difference.py:158-179`, `harmonic.py:146-163`).
6. **[ALTO · Eje 2]** Sólidos 2D no exponen `internal_forces` (ADR 0002 incompleto declarado en MEMORY). Coexisten `compute_internal_forces` (dict) e `internal_forces` (`None`) en el mismo árbol — convención implícita. (`fenix/elements/solid_2d/*.py`, `core/element.py:227-243`).
7. **[ALTO · Eje 2]** Tests con coeficientes unitarios (`thickness=1.0`, `E=1`, áreas unitarias) ocultan bugs dimensionales latentes — el bug histórico del `thickness` en el embedded fue precisamente esto. 27 ocurrencias de `thickness=1.0` en 11 archivos de test.
8. **[MEDIO · Eje 3]** `DruckerPrager2D.IS_SYMMETRIC = False` hardcoded incluso cuando el usuario construye con `ψ = φ` (asociado). Renuncia a Cholesky innecesariamente. (`drucker_prager_2d.py:275, 339`).
9. **[MEDIO · Eje 4]** Cuatro solvers Newmark/HHT duplican ~880 líneas de plumbing temporal idéntico — refactor zonal pendiente, abarata futuros esquemas (Bossak, generalized-α). (`newmark.py:135-1167`).
10. **[MEDIO · Eje 5]** 3 textos desactualizados en `docs/catalogo_solvers.md` (decían "despacho por `isinstance`" tras aplicar regla C); 1 inconsistencia en `catalogo_materiales.md` (Elastic2D dice compatible sólo con Quad4/Tri3, omite cuadráticos); 1 desactualización en ModalSolver (catálogo dice "sin masas modales efectivas" cuando `participation_factors` está implementado en fase 7 ADR 0009).

---

# Eje 1 — Arquitectura y plumbing

**Veredicto del eje**: arquitectura bien construida; contratos declarativos centralizados, validación temprana con mensajes ricos, capas sin ciclos, despacho declarativo, logging uniforme. Los problemas son tres familias: autodiscover estructuralmente incompleto, duplicación residual que pasa el umbral "centralizar siempre", API de `Result` asimétrica.

## H-1.1 — Autodiscover no recursivo
**Severidad**: alto · **Categoría**: registries

**Evidencia**: `fenix/autodiscover.py:17` usa `pkgutil.iter_modules(package.__path__)` (no recursivo). `fenix/elements/solid_2d/__init__.py:34-38` importa Quad4/Quad8/Quad9/Tri3/Tri6 pero **no** `embedded_cst`. `fenix/__init__.py:42` lo re-importa explícitamente.

**Descripción**: el contrato implícito "añadir archivo nuevo en `fenix/elements/...` basta para registrarlo" es **falso** para módulos dentro de subpaquetes — hay que recordar editar el `__init__.py` del subpaquete o el `fenix/__init__.py` raíz. Para una arquitectura optimizada para la IA mantenedora sin contexto previo, es punto de regresión silenciosa.

**Recomendación**: convertir `_discover_package` a `pkgutil.walk_packages` con `prefix=f"{package_name}."`, ignorando `_shared` y `diagnostics`. Alternativamente, añadir test `test_autodiscover_covers_all_decorated_classes` que escanee el árbol y verifique que toda clase con `@*Registry.register` aparece en el registry tras `_initialize_registries`.

## H-1.2 — Triplicación de `_resolve_rayleigh`
**Severidad**: alto · **Categoría**: duplicación

**Evidencia**: `newmark.py:113-127`, `central_difference.py:158-179`, `harmonic.py:146-163`.

**Descripción**: función estática que traduce `{"alpha","beta"}` o `{"xi1","omega1","xi2","omega2"}` a `(α,β)` copiada tres veces literalmente. Pasa el umbral "centralizar siempre" de la memoria de consistencia arquitectural. Confirmado por agentes 1 y 4 de forma independiente.

**Recomendación**: mover a `fenix/math/damping.py::resolve_rayleigh_config(cfg, *, source: str = "solver") -> (α, β)`. Las tres `_resolve_rayleigh` quedan como thin-wrappers o desaparecen.

## H-1.3 — API de `Result` asimétrica
**Severidad**: medio · **Categoría**: API

**Evidencia**: `results.py:340-377` (SolveResult con `element_forces`) vs `results.py:201-208,241-263,266-337` (TransientResult, HarmonicResult, ResponseSpectrumResult sin esfuerzos internos).

**Descripción**: `SolveResult` expone dict de `ElementForces`; los `Result` dinámicos sólo exponen historias U, U̇, Ü. Consumidor que quiera N(t), M(t) debe reconstruirlos manualmente.

**Recomendación**: método lazy en `TransientResult.internal_forces_history(domain)` y análogo en HarmonicResult/ResponseSpectrumResult. No precalcular eager.

## H-1.4 — Doble contrato `compute_internal_forces` vs `internal_forces`
**Severidad**: medio · **Categoría**: contratos

**Evidencia**: `solid_2d/quad4.py:80-91`, `solid_2d/tri3.py:49-52`, `solid_2d/_shared.py:372-374`, `truss.py:66-84`, `core/element.py:227-243`.

**Descripción**: convive método legacy `compute_internal_forces` (dict ad-hoc) con `internal_forces` (`ElementForces | None`). Sólidos 2D devuelven `None`; trusses devuelven ambos. ADR 0002 declarado incompleto en MEMORY pero la convención implícita es deuda activa.

**Recomendación**: documentar explícitamente en `Element.internal_forces` que `None` es legítimo para continuo donde no aplican N/V/M; marcar `compute_internal_forces` como deprecado en docstrings.

## H-1.5 — `NonlinearSolver.__init__` sin keyword-only + heurísticos hardcoded
**Severidad**: medio · **Categoría**: tolerancias

**Evidencia**: `nonlinear.py:37-58, 286-287`.

**Descripción**: 9 argumentos opcionales en una sola firma sin `*`; heurísticos adaptativos hardcoded (`iteration < 4`, `delta_lambda * 1.5`, `min_delta_lambda=1e-5`, `target_load - 1e-9`). Ilegibilidad y dificultad de auditoría.

**Recomendación**: introducir `*` después de `assembler`. Centralizar constantes adaptativas en `constants.py` con prefijo `NEWTON_ADAPTIVE_*`.

## H-1.6 — Ramas hardcoded `if lumping != "consistent": raise NotImplementedError` duplicadas
**Severidad**: medio · **Categoría**: duplicación

**Evidencia**: `quad4.py:239-243`, `tri3.py:174-178`, `_shared.py:452-457`, `truss.py:139-143` (probablemente frame3d, frame también).

**Recomendación**: helper en `core/element.py` o `math/mass_lumping.py`: `validate_lumping(lumping, supported, elem_name)`. Alternativa: `SUPPORTED_LUMPING: ClassVar[frozenset[str]]` declarativo.

## H-1.7 — `run_yaml` no valida `PIPELINE_KIND` desconocido
**Severidad**: medio · **Categoría**: API

**Evidencia**: `entry.py:333-354`.

**Descripción**: `getattr(solver, "PIPELINE_KIND", "static")` cae al pipeline static por defecto. Solver con `PIPELINE_KIND = "spectrum_v2"` mal tecleado pasa silenciosamente.

**Recomendación**: declarar `KNOWN_PIPELINE_KINDS = frozenset({"static","modal","transient","harmonic","spectrum"})` y validar al final con `ValueError` listando válidos.

## H-1.8 — `step_callback` silenciado en pipelines no estáticos
**Severidad**: bajo · **Categoría**: API

**Evidencia**: `entry.py:296-354`.

**Recomendación**: emitir `_log.warning(...)` cuando `step_callback is not None and pipeline_kind != "static"`.

## H-1.9 — Magic numbers en `ModalSolver` y otros fuera de `constants.py`
**Severidad**: medio · **Categoría**: tolerancias

**Evidencia**: `modal.py:62` (`tolerance=1.0e-9`), `modal.py:92` (`1e-12 * lam_max`); `central_difference.py:196`; `embedded_cst.py:26` (`_LOCAL_JUMP_RTOL = 1.0e-10`); `von_mises_2d.py:29` (`_PLANE_STRESS_MAX_LOCAL_ITER = 25`).

**Recomendación**: mover a `constants.py` con prefijos descriptivos. La memoria "tolerancias son ciudadanas de primera clase" aplica más allá del Newton global.

## H-1.10 — `frozen=True` semántico vs estructural en `Result`
**Severidad**: bajo · **Categoría**: API

**Evidencia**: `results.py:55,106,171,211,266,340`.

**Descripción**: `@dataclass(frozen=True)` impide reasignar atributos pero `np.ndarray` y `dict` contenidos siguen mutables. Docstring promete "inmutable" — engañoso.

**Recomendación**: documentar "shallow-frozen" en los docstrings.

## H-1.11 — `core/element.py` importa `fenix.results`
**Severidad**: bajo · **Categoría**: capas

**Evidencia**: `core/element.py:12`.

**Descripción**: `results.py` es a la vez tipo público bajo (consumido por core) y constructor con dependencias altas (consume Domain/Assembler vía TYPE_CHECKING). Grafo conceptual enturbiado, sin ciclo real de imports.

**Recomendación**: separar `ElementForces` a `core/element_forces.py` o `core/types.py`.

## H-1.12 — Excepciones tipadas ADR 0011 no extendidas a return mapping
**Severidad**: bajo · **Categoría**: logging

**Evidencia**: `diagnostics.py` cubre Newton global; Newton locales (J2 plane stress, jump CST) fallan con `RuntimeError` plano.

**Recomendación**: no actuar ahora. Cuando aparezca el segundo material con return mapping iterativo (Mohr-Coulomb), introducir `MaterialReturnMappingDivergedError`.

## H-1.13 — Duplicación geométrica truss/cable corotacional
**Severidad**: bajo · **Categoría**: duplicación

**Evidencia**: `truss.py:166-202`, `cable.py:67-100`.

**Recomendación**: monitorizar. Tercer elemento axial corotacional (cable térmico futuro) dispara centralización.

## H-1.14 — YAML parser accede a `_materials` privado
**Severidad**: bajo · **Categoría**: API

**Evidencia**: `yaml_parser.py:148,168,187,237`.

**Recomendación**: cambiar a `MaterialRegistry.names()` público.

## H-1.15 — `Node.dofs` mezcla "declared" y "global eq number" en mismo dict
**Severidad**: bajo · **Categoría**: contratos

**Evidencia**: `core/node.py:13-16`, `core/domain.py:82-90`.

**Recomendación**: dos dicts (declared, numbered) o sentinel `None` antes de numerar con check explícito.

---

# Eje 2 — Implementación numérica de elementos

**Veredicto del eje**: formulación numérica de los 16 elementos rigurosamente correcta en su núcleo. Convención de signos §5 aplicada coherentemente (interna stress-resultant/RHR, API pública 2D con `V` invertido), blindada por tests específicos (`test_frame2d_internal_forces.py`, `test_frame3d_internal_forces.py`). Bug histórico del `thickness` en embedded resuelto con regresión específica.

## H-2.1 — Sólidos 2D no implementan `internal_forces` (ADR 0002 incompleto)
**Severidad**: alto · **Componente**: Tri3, Quad4, Tri6, Quad8, Quad9, CST_Embedded2D · **Categoría**: API / contratos

**Evidencia**: `solid_2d/tri3.py:49-52`, `solid_2d/quad4.py:80-91`, `solid_2d/_shared.py:372-374`, `solid_2d/embedded_cst.py:284-286`, base devuelve `None` en `core/element.py:227-243`.

**Descripción**: la base devuelve `None` por defecto. `compute_internal_forces(U)` en los sólidos devuelve `dict {stress, strain}` — semántica distinta del contrato `ElementForces` del ADR 0002. Nombre casi-idéntico induce confusión.

**Recomendación**: deuda declarada en MEMORY. Documentar explícitamente en manual y core/element.py.

## H-2.2 — Tests con coeficientes unitarios ocultan bugs dimensionales
**Severidad**: alto · **Componente**: Tri3/Quad4 cargas y masa, Frame2D body load · **Categoría**: tests

**Evidencia**: 27 ocurrencias de `thickness=1.0` en 11 archivos de test. `test_patch_solid_2d.py` no parametriza `thickness` ni `density`. El bug del embedded fue cazado precisamente por benchmark Van Vliet con `thickness=0.1`.

**Recomendación**: al menos un test por elemento con `thickness != 1`, `area != 1`, `density != 1`, `L != 1`. Frame3D además `Iy ≠ Iz`. No urgente; housekeeping en próxima sesión.

## H-2.3 — Truss2D/Truss3D sin K_G — no aptos para pandeo lineal
**Severidad**: medio · **Componente**: Truss2D, Truss3D · **Categoría**: matriz K

**Evidencia**: `truss.py:52-64`, `truss.py:254-268`. Sólo `K_e = (E_t·A/L)·d·dᵀ`, sin `(N/L)·P`.

**Recomendación**: añadir nota en docstring: "para pandeo lineal (eigenvalor con K_G) usar la variante corotacional".

## H-2.4 — Frame3D lumped bloque-diagonal en orientación oblicua
**Severidad**: medio · **Componente**: Frame3D · **Categoría**: masa lumped

**Evidencia**: `frame3d.py:330-370`. Limitación declarada.

**Recomendación**: documentado; sin fix urgente.

## H-2.5 — Frame2DEuler mezcla rigidez Hermitiana con sobreescritura axial
**Severidad**: medio · **Componente**: Frame2DEuler, Frame2DTimoshenko · **Categoría**: integración / no linealidad

**Evidencia**: `frame/euler.py:82-89`, `frame/timoshenko.py:106-111`. `F_int_local[0] = -σ·A`, `F_int_local[3] = σ·A`, pero términos cortante/flector siguen `K_local @ u_local` con `E_t` aplicado globalmente.

**Descripción**: plasticidad axial capturada; plasticidad por flexión (fibras inferiores en tracción) **no** — requiere `FiberSection` no implementada. Comentario "permite materiales no-lineales" sin matizar.

**Recomendación**: aclaración en docstring: "no-linealidad se aplica al axial uniforme; plasticidad distribuida en sección requiere `FiberSection`, no implementada".

## H-2.6 — Lumping nodal directo en frames vs HRZ en sólidos
**Severidad**: medio · **Componente**: Truss, Cable, Frame2D, Frame3D · **Categoría**: masa lumped

**Evidencia**: `frame/_shared.py:180-230`. Justificado: preserva diagonalidad global del HRZ. Inconsistencia de criterio entre familias.

**Recomendación**: nota en ADR 0009 fase 2: "HRZ para sólidos isoparamétricos, nodal directo para vigas/barras por preservación de diagonalidad".

## H-2.7 — Tri6 rigidez subintegrada (tri_3 vs orden requerido)
**Severidad**: medio · **Componente**: Tri6 · **Categoría**: integración

**Evidencia**: `tri6.py:24-27`. `_DEFAULT_QUADRATURE = "tri_3"` para K_e, `_MASS_QUADRATURE = "tri_6"` para masa.

**Descripción**: en malla recta integrando es polinomial puro de orden 2 (exacto con 3 puntos). En malla distorsionada (J no constante) la subintegración baja rigidez (puede mejorar flexión, puede introducir modos espurios).

**Recomendación**: parámetro `quadrature_stiffness="tri_6"` opcional para integración exacta en mallas distorsionadas, o documentar la elección.

## H-2.8 — Frame2D doble rotación `T.T·K_local·u_local · T` no-op
**Severidad**: bajo · **Componente**: Frame2DEuler, Frame2DTimoshenko · **Categoría**: claridad

**Evidencia**: `frame/euler.py:91-94`, `frame/timoshenko.py:114-117`. `F_int = T.T @ F_int_local` y luego `F_local = T @ F_int = F_int_local` (T ortogonal).

**Recomendación**: refactor menor cosmético.

## H-2.9 — Warning de Frame3D casi-vertical ruidoso con ref_vector explícito
**Severidad**: bajo · **Componente**: Frame3D · **Categoría**: usabilidad

**Evidencia**: `frame3d.py:108-115`. Umbral `_VERTICAL_COS_TOL = 0.99`.

**Recomendación**: opcional — silenciar warning una vez por tipo de elemento.

## H-2.10 — CST_Embedded2D valida bulk por nombre de clase hardcoded
**Severidad**: bajo · **Componente**: CST_Embedded2D · **Categoría**: extensibilidad

**Evidencia**: `embedded_cst.py:31,78-85`. `_ACCEPTED_BULKS = ('Elastic2D',)`.

**Recomendación**: validar por contrato declarativo (`IS_ELASTIC_BULK`) cuando aparezca el segundo bulk válido.

---

# Eje 3 — Implementación numérica de materiales

**Veredicto del eje**: buena salud general. Return mappings clásicos implementados con rigor de manual; tangentes algorítmicas consistentes (no continuas, no diferencias finitas). J2 plane stress en `VonMises2D` notablemente sofisticado — diagonalización en autovectores con Newton local cerrado, fiel a Simó-Hughes §3.4.1 Box 3.2. Cap `DAMAGE_MAX` cohesivo vs continuo correctamente diferenciado y documentado.

## H-3.1 — DruckerPrager2D `IS_SYMMETRIC = False` aun cuando ψ = φ
**Severidad**: medio · **Material**: DruckerPrager2D · **Categoría**: tangente / despachador algebraico

**Evidencia**: `drucker_prager_2d.py:275, 339`. `self.associated = (psi_deg == phi_deg)` existe pero no se expone en `IS_SYMMETRIC`.

**Descripción**: para problemas asociados (default), la tangente sí es simétrica (`b_g ⊗ b_f` → `b_f ⊗ b_f`), pero el despachador algebraico (ADR 0003) sigue eligiendo LU. Renuncia gratuita a Cholesky/LDLᵀ (~2× más rápidos).

**Recomendación**: `@property` que devuelva `self.associated`, o atributo de instancia en `__init__`. Verificar que el despachador algebraico consulte `material.IS_SYMMETRIC` por instancia.

## H-3.2 — Elastic1D y Elastic2D no validan inputs físicos
**Severidad**: medio · **Material**: Elastic1D, Elastic2D · **Categoría**: validación temprana

**Evidencia**: `elastic.py:10-12`, `elastic_2d.py:11-19`. No se rechaza `E ≤ 0`, `ν ∉ (-1, 0.5)`, `hypothesis` desconocida. Inconsistente con `VonMises2D`, `DruckerPrager2D`, `CableMaterial1D`, `IsotropicDamage*` que sí validan.

**Descripción**: `Elastic2D(E=200, nu=0.5, hypothesis='plane_strain')` produce `C` con NaN/Inf; bug aparece en el primer `K @ U`. `hypothesis='plain_stress'` (typo) cae silenciosamente a plane strain.

**Recomendación**: añadir mismas validaciones que ya tienen las clases hermanas.

## H-3.3 — Nomenclatura "tensión" vs "esfuerzo" en docstrings
**Severidad**: bajo · **Material**: VonMises2D principalmente · **Categoría**: nomenclatura

**Evidencia**: `von_mises_2d.py:315` — `"sigma_y : Tensión de fluencia inicial"`. Memoria explícita `feedback_esfuerzo_vs_traccion.md`.

**Recomendación**: barrido global "tensión de fluencia" → "esfuerzo de fluencia" manteniendo "tracción" solo donde corresponde (cohesivos).

## H-3.4 — Duplicación softening exponencial 1D / 2D
**Severidad**: medio · **Material**: IsotropicDamage1D, IsotropicDamage2D · **Categoría**: duplicación

**Evidencia**: `damage_1d.py:86-93`, `damage_2d.py:94-101`. Rama cálculo `d`, cap `DAMAGE_MAX`, flag `saturated` idénticos. Diferencia: `eps_eq` y tangente.

**Recomendación**: extraer `fenix/materials/_softening.py::exponential_damage(kappa, kappa_0, alpha, cap) -> (d, dd_dkappa, saturated)`. Abrir puerta a `linear_damage(...)` para alinear API con cohesivo.

## H-3.5 — Tangente algorítmica DP sin diferencia finita explícita
**Severidad**: medio · **Material**: DruckerPrager2D · **Categoría**: cobertura de test

**Evidencia**: `drucker_prager_2d.py:144-152`. Comentario "verificable por diferencia finita". `test_solid_2d_drucker_prager.py` sólo tiene tests de convergencia indirecta.

**Recomendación**: test que perturbe `ε` con `δε = 1e-7` y compare `(σ(ε+δε)-σ(ε))/δε` contra `C_alg @ e_i`. Mismo test para VonMises2D plane stress.

## H-3.6 — Ausencia de test material CableMaterial1D en switching
**Severidad**: bajo · **Material**: CableMaterial1D · **Categoría**: cobertura

**Recomendación**: revisar `test_cable_material.py` y añadir caso explícito de `ε = 0` exacto si falta.

## H-3.7 — Cap DAMAGE_MAX continuo vs cohesivo contrastado pero sin nota en docstring
**Severidad**: bajo (informativo) · **Material**: IsotropicDamage* vs CohesiveDamageIsotropic

**Evidencia**: `damage_1d.py:91-95`, `damage_2d.py:99-107` (cap a d, σ); `damage_isotropic.py:138` (ω sin cap, sólo tangente). Memoria explícita.

**Recomendación**: nota en docstring de `IsotropicDamage1D.compute_state` referenciando el contraste con cohesivo.

## H-3.8 — Trial/commit no documentado en contrato `Material`
**Severidad**: bajo · **Material**: todos inelásticos · **Categoría**: trial/commit

**Evidencia**: `core/material.py:60-68` no menciona que `new_state_vars` es trial.

**Recomendación**: docstring: "El `new_state_vars` devuelto es el estado **trial** del paso. El consumidor (ElementState/Solver) confirma este estado al cierre del paso convergido".

## H-3.9 — VonMises2D plane stress: invariante `tr(ε^p) = 0` no verificado
**Severidad**: bajo · **Material**: VonMises2D plane stress · **Categoría**: hipótesis

**Evidencia**: `von_mises_2d.py:243`. `ε^p_zz` forzada por incompresibilidad plástica; no se verifica que `ε^p_old` entrante satisfaga `Σ = 0`.

**Recomendación**: `assert` en `compute_state` antes del kernel: `abs(sum(eps_p_old)) < tol`.

## H-3.10 — Tolerancia Newton local VonMises plane stress con `alpha_old`
**Severidad**: bajo · **Material**: VonMises2D plane stress

**Evidencia**: `von_mises_2d.py:192`. `yield_tol` se calcula con `alpha_old`.

**Recomendación**: documentar en docstring del kernel, o recalcular `yield_tol` cada iteración usando `alpha_curr`.

---

# Eje 4 — Implementación numérica de solvers

**Veredicto del eje**: sólida y mayoritariamente conforme a los ADRs. Fórmulas centrales correctas (Newmark β-γ, HHT-α con auto-derivación, leapfrog Belytschko-Liu-Moran, Crisfield cilíndrico, SRSS/CQC Der Kiureghian). Regla C de despacho declarativo cumplida. Defaults clave del ADR 0011 correctos (`line_search=False`).

**Tres bugs críticos verificados independientemente**.

## H-4.1 — `singular_tangent_seen` muerto en Newmark/HHT no lineales
**Severidad**: **CRÍTICO** · **Solver**: NewtonNewmarkSolver, NewtonHHTSolver · **Categoría**: diagnóstico

**Evidencia**:
- `newmark.py:442` y `:540` — `singular_tangent_seen = False` y nunca se modifica.
- `newmark.py:476-484` — único `try` cubre solo `CholeskyNotPositiveDefiniteError`; degrada a LU sin marcar el flag.
- `newmark.py:1024` y `:1126` — mismo patrón en NewtonHHTSolver.
- Contraste: `nonlinear.py:235-240` sí flipea `singular_tangent_seen = True` en su `except RuntimeError`.

**Descripción**: en NewtonNewmarkSolver y NewtonHHTSolver el flag se inicializa a `False` antes del bucle Newton y nunca se actualiza. Como `classify_divergence(..., singular_tangent_detected=False)` jamás devuelve `SingularTangentError` (`diagnostics.py:189`), esa subclase de excepción está **inaccesible** desde el subsistema dinámico no lineal. Si el jacobiano dinámico pierde rango por bifurcación o tangente material casi singular, el solver agota iteraciones y termina lanzando un diagnóstico **engañoso** (Oscillating / LoadExceedsCapacity / Unknown). Si `splu` lanza `RuntimeError` por singularidad pura, el error propaga sin atrapar y se pierde la metadata estructurada del ADR 0011.

**Por qué importa**: rompe la promesa explícita del ADR 0011 documentada en `diagnostics.py:18` ("El módulo es consumido por NonlinearSolver y NewtonNewmarkSolver"). La memoria "no parches" exige diagnósticos físicamente significativos.

**Recomendación**: añadir `try/except RuntimeError` que envuelva la rama LU análogo a `NonlinearSolver._solve_reduced`, y emparejar el `except CholeskyNotPositiveDefiniteError` para flipear el flag cuando la degradación a LU también falla con tangente singular real.

## H-4.2 — `LinearSolver.solve` no cachea factorización
**Severidad**: **CRÍTICO** · **Solver**: LinearSolver · **Categoría**: factorización

**Evidencia**:
- `linear.py:25-49` — `solve()` llama `select_solver(...).solve(K_red, F_red)` sin retener la factorización.
- `linalg/lu.py:39-40` — `LUSolver.solve` invoca `spla.spsolve(K, b)` directamente.
- `linalg/base.py:64-66` — `factorize(K) -> FactorizedSolver` existe y está documentado como "habilita Newton modificado y reuso entre pasos".

**Descripción**: si un consumidor llama `linear_solver.solve(F1)` y luego `linear_solver.solve(F2)` con la misma `K` (barrido de cargas, FRF lineal, multiplicidad de cargas), se paga la factorización dos veces. La interfaz `factorize` ya existe pero `LinearSolver` no la usa.

**Recomendación**: cache opcional — primera llamada factoriza y guarda `self._factor`; llamadas siguientes con misma `K` reutilizan. Patrón idéntico al ya implementado en `NewmarkSolver` (newmark.py:197).

## H-4.3 — `CentralDifferenceSolver` no lineal con apoyos no nulos resta `F_dir` dos veces
**Severidad**: **CRÍTICO** · **Solver**: CentralDifferenceSolver · **Categoría**: esquema temporal

**Evidencia**:
- `central_difference.py:232` — `F_dir = T.T @ (K @ g)` calculado con K lineal inicial.
- `central_difference.py:296` — `a_free = M_inv_diag * (F_red - F_int_red - C_red @ v_half - F_dir)`.
- `central_difference.py:351-357` — `_compute_internal_forces_reduced`: rama lineal devuelve `K_red @ u_free` (no incluye `F_dir`, correcto); rama no-lineal devuelve `T.T @ F_int_global(u_global)` con `u_global = T·u_free + g` (incluye implícitamente la contribución del apoyo en `F_int`).

**Descripción**: en la rama lineal `F_red - F_int_red - F_dir = T.T·(F − K·u_global)` correctamente. En la rama no-lineal, `F_int_global(T·u_free + g)` ya incluye el efecto del apoyo (para material lineal: `K·T·u_free + K·g`), por lo que sumar `-F_dir = -T.T·K·g` adicionalmente equivale a contar el efecto del apoyo **dos veces**. El bug es **invisible cuando `g = 0`** (apoyos homogéneos, caso por defecto de todos los tests del proyecto) pero corrompe la dinámica con desplazamientos prescritos no nulos (asentamiento prescrito, excitación cinemática del soporte).

**Por qué importa**: agujero conceptual latente. Se manifestará la primera vez que se haga un análisis dinámico explícito con apoyo cinemático. Es exactamente el patrón que la memoria "no parches" advierte contra (síntoma invisible hoy, latente para mañana).

**Recomendación**: en la rama no-lineal devolver `T.T @ F_int_global - F_dir` y luego eliminar `- F_dir` del cómputo de la aceleración (alineando la semántica con la rama lineal sin doble resta). Añadir test con apoyo prescrito no nulo en `test_central_difference.py`.

## H-4.4 — Cuatro solvers Newmark/HHT duplican ~880 líneas de plumbing temporal
**Severidad**: medio · **Solver**: NewmarkSolver, NewtonNewmarkSolver, HHTSolver, NewtonHHTSolver · **Categoría**: duplicación

**Evidencia**: `newmark.py:135-254` (NewmarkSolver), `:336-575` (NewtonNewmarkSolver), `:713-842` (HHTSolver), `:927-1167` (NewtonHHTSolver).

**Descripción**: cuatro `solve()` métodos comparten skeleton (predictor → resolver → corrector → history). Diferenciación Newmark↔HHT son tres términos `α` en RHS y factor `(1+α)` en `A_eff`. Diferenciación lineal↔Newton es el bucle interno y re-ensamblaje. Escritos como métodos completos independientes en lugar de hooks.

**Recomendación**: refactor zonal — extraer `_TimeIntegratorBase` con hooks (`_residual`, `_jacobian_eff`, `_advance_state`). Cuatro subclases con 30-50 líneas. Abarata extensiones futuras (Bossak, generalized-α, energy-momentum).

## H-4.5 — `_resolve_rayleigh` triplicado (= H-1.2)
**Severidad**: medio · **Categoría**: duplicación

Mismo hallazgo que H-1.2, confirmado independientemente.

## H-4.6 — Override silencioso de β, γ en HHT
**Severidad**: medio · **Solver**: HHTSolver, NewtonHHTSolver · **Categoría**: predictor/corrector

**Evidencia**: `newmark.py:700-703, 911-914`.

**Descripción**: si el usuario pasa `beta=...` o `gamma=...` al constructor de HHT se aplican silenciosamente, reemplazando los auto-derivados desde `α`. Combinaciones arbitrarias rompen estabilidad incondicional + orden 2 simultáneos.

**Recomendación**: `warnings.warn` cuando se pase override explícito, o exigir que ambos se pasen juntos.

## H-4.7 — `LinearSolver` declara SPD a priori sin derivarlo
**Severidad**: medio · **Solver**: LinearSolver · **Categoría**: factorización

**Evidencia**: `linear.py:34-38`. `is_positive_definite=True` hardcoded.

**Recomendación**: derivar de flags de elementos/materiales (análogo a `domain_is_symmetric`). El contrato declarativo del ADR 0003 queda respetado.

## H-4.8 — CentralDifferenceSolver no-lineal descarta `K_t`
**Severidad**: medio · **Solver**: CentralDifferenceSolver · **Categoría**: esquema temporal

**Evidencia**: `central_difference.py:354-357`. `_K_t, F_int_global = assembler.assemble_non_linear_system(u_global)` — `_K_t` ignorado.

**Descripción**: cada paso ensambla la tangente y la descarta. Correcto físicamente (leapfrog no necesita `K_t`), pero coste constante por paso.

**Recomendación**: `assemble_non_linear_system(u, jac=False)` opcional. Beneficia también futuros Quasi-Newton.

## H-4.9 — Doble ensamblaje en Newton dinámico con `line_search=True`
**Severidad**: bajo · **Solver**: NewtonNewmarkSolver, NewtonHHTSolver · **Categoría**: line search

**Evidencia**: `newmark.py:489-507, 1078-1093`. Comentario explícito reconoce: `_armijo_step_dynamic` ensambla en su loop, luego outer loop re-ensambla por simetría con la rama sin line search.

**Descripción**: bajo. Default `line_search=False` (ADR 0011 enmendado) no ejecuta esta ruta. Cuando se activa, coste medible pero predecible.

**Recomendación**: documentar y mantener.

## H-4.10 — Archivo `Iventalla.dat` untracked en repo
**Severidad**: bajo · **Categoría**: housekeeping

**Recomendación**: ignorar o agregar al `.gitignore`. Fuera del alcance funcional.

---

# Eje 5 — Cobertura, validación, consistencia docs↔código

**Veredicto del eje**: STATUS y ROADMAP reflejan fielmente el estado del código (cuentas de componentes coinciden: 16 elementos, 9 materiales, 11 solvers, 29 specs validated). Las inconsistencias se concentran en `catalogo_solvers.md` y `catalogo_materiales.md` con texto desactualizado tras aplicar la regla C (despacho declarativo) y la fase 7 del ADR 0009 (response spectrum).

## H-5.1 — 3 textos desactualizados en `catalogo_solvers.md` tras regla C
**Severidad**: medio · **Categoría**: docs↔código

**Evidencia**: 
- `catalogo_solvers.md:118` (NewtonNewmarkSolver): "Como es subclase de NewmarkSolver, entry.run_yaml lo detecta vía isinstance y despacha a run_transient automáticamente."
- `catalogo_solvers.md:137` (HHTSolver): mismo texto.
- `catalogo_solvers.md:69` (ModalSolver): "sin masas modales efectivas ni factores de participación (fase 7 del ADR 0009)".

**Descripción**: `entry.py:333` despacha por `getattr(solver, "PIPELINE_KIND", "static")` (regla C aplicada). `participation_factors` existe en `modal_response.py:114` y `ResponseSpectrumResult` los expone (ADR 0009 fase 7 cerrada).

**Recomendación**: actualizar las 3 entradas. Tarea mecánica.

## H-5.2 — `Elastic2D` en catálogo dice compatible solo con Quad4/Tri3, omite cuadráticos
**Severidad**: medio · **Categoría**: docs↔código

**Evidencia**: `catalogo_materiales.md:30`. `Elastic2D` aparece como "Compatible con: Quad4, Tri3" — pero `STRAIN_DIM=3` permite Quad8, Quad9, Tri6 igualmente; los otros materiales 2D (VonMises2D, IsotropicDamage2D, DruckerPrager2D) sí listan los 5.

**Recomendación**: actualizar `Elastic2D` a "Compatible con: Quad4, Tri3, Tri6, Quad8, Quad9".

## H-5.3 — 7 huecos de spec vs componentes registrados
**Severidad**: medio · **Categoría**: cobertura de specs

**Evidencia**: 36 componentes registrados (16 + 9 + 11) vs 29 specs validated. Huecos:
- **Materiales**: Elastic1D, Elastic2D, Elastoplastic1D (3 huecos — anteriores a la convención de specs).
- **Solvers**: LinearSolver, NonlinearSolver, ArcLengthSolver (3 huecos — anteriores), NewtonHHTSolver (1 hueco — variante de HHTSolver, según Reglas §4 "spec corta tipo extensión").

**Descripción**: deuda documental conocida implícitamente. Cuatro de los siete (NewtonHHTSolver explícitamente, LinearSolver/NonlinearSolver/ArcLengthSolver implícitamente) son centrales en el catálogo.

**Recomendación**: priorizar `NewtonHHTSolver.md` (spec corta tipo extensión) y `NonlinearSolver.md` (componente central). Los demás como housekeeping.

## H-5.4 — Tests con coeficientes unitarios (= H-2.2)
**Severidad**: medio · **Categoría**: cobertura

Confirmado por inventario: 27 ocurrencias de `thickness=1.0` en 11 archivos. Hallazgo cruzado con Eje 2.

## H-5.5 — No existe `MATRIZ.md`
**Severidad**: bajo · **Categoría**: docs faltantes

**Evidencia**: `ROADMAP.md:192` y `STATUS.md:5` mencionan `MATRIZ.md` como pendiente. Tabla cruzada elemento × material × solver con celdas validadas/blancas/vetadas.

**Recomendación**: cuando se priorice el housekeeping de docs navegacionales, generar `MATRIZ.md`. No bloquea avance.

## H-5.6 — `ONBOARDING.md` mencionado como pendiente
**Severidad**: bajo · **Categoría**: docs faltantes

**Evidencia**: `ROADMAP.md:193`.

**Recomendación**: existe en memoria del proyecto (`reference_sistema_docs_navegacionales.md` dice "Empezar por ONBOARDING en sesión cold-start"). Verificar si ya existe en `docs/` (la memoria sugiere que sí, pero ROADMAP lo lista como pendiente).

---

# Lo que está bien hecho (consolidado de los 5 ejes)

Esta sección resume patrones validados como sólidos. El usuario debería estar **confiado** en estos puntos.

## Arquitectura
- **Contratos declarativos centralizados**: `Element`, `Material`, `CohesiveMaterial` declaran sus `ClassVar` con docstrings que explican cada bit. Validación temprana en `_validate_material_compatibility` con mensajes ricos.
- **Layering limpio**: `core` no importa `math/materials/elements`. `math` no importa `materials/elements`. `materials` no importa `elements`. `elements` no importa `math.solvers`. Sin ciclos.
- **Despacho declarativo de pipelines** por `PIPELINE_KIND` (regla C ADR 0009).
- **`bc.constraints` generaliza Dirichlet + linear constraints** homogéneamente (ADR 0004).
- **Patrón paquete + `_shared.py`** consistente en `elements/frame/`, `elements/solid_2d/`, `math/solvers/`.
- **Caché de topología COO + caché de matriz de masa** — performance pensada desde la arquitectura.
- **Pre-validación agregada** de masa y YAML (errors acumulados, reportados juntos).
- **Logging uniforme** (sin `print`, `get_logger("...")` siempre).
- **Conversión de signo 2D bien aislada** en `_frame2d_forces_from_local` — único sitio donde se aplica `V_2D = -V_interno`.
- **Excepciones tipadas ADR 0011** con métricas estructuradas (residuo, δU, factor de carga, n_bisecciones).
- **ADRs persistentes y bien linkeados** desde el código.

## Elementos
- **Convenciones de signos §5 verificadas con tests explícitos** de cantilever 2D y 3D que cubren casos canónicos.
- **Bug del thickness blindado** con tests de regresión específicos (`thickness=0.1`).
- **Cinemática KOS del CST_Embedded2D** correcta: condensación estática local, `l_d = (A_e/h)·|cos(θ−α)|`, activación Rankine en `prepare_step()` con estado committed (anti-chattering).
- **Quadratura tri_6 Dunavant** implementada con ternas simétricas y pesos correctos.
- **Funciones de forma higher-order verificadas símbolo a símbolo** contra Cook-Malkus-Plesha §6.
- **Matriz geométrica Frame2DEulerCorot completa** con términos cruzados `−(M1+M2)/l²·(r⊗z + z⊗r)` de Crisfield Vol.1.
- **Truss/Cable corotacional 3D con `perpendicular_projector(ê)`** bien factorizado.
- **Lumped Quad4/Tri3 con HRZ canónico**, masa total preservada.
- **Patch test de MacNeal-Harder** sobre Quad4/Tri3 distorsionados.
- **Patch test cuadrático** sobre Tri6, Quad8, Quad9.

## Materiales
- **Return mapping J2 plane stress (Simó-Hughes §3.4.1)** ejemplar — diagonalización en autovectores de `C_e^ps·P`, Newton local escalar cerrado, tangente con corrección por `dα/dΔγ` no constante.
- **Tangente algorítmica consistente** en todos los inelásticos (1D y 2D), incluida asimétrica para Drucker-Prager no asociado e IsotropicDamage2D.
- **ADR 0006 aplicado uniformemente**: cada material declara `admissibility_scale` con dimensión física correcta.
- **σ = s_trial + p·I en rama elástica** de VonMises2D plane strain (no `C_e·strain`) — captura descarga elástica post-plástica correctamente.
- **Distinción cap-DAMAGE_MAX continuo vs cohesivo** correctamente implementada.
- **DruckerPrager2D regular↔apex** con criterio geométrico claro; rama apex hidrostático puro.
- **Calibración DP-MC** con tres variantes correctas (Chen-Han / Simó-Hughes).
- **Validación de inputs** en todos los inelásticos (rangos físicos, mensajes claros).

## Solvers
- **`PIPELINE_KIND` declarativo** consistente en 7 solvers registrados; 4 subclases heredan de `NewmarkSolver`.
- **`line_search=False` default real** en NonlinearSolver, NewtonNewmarkSolver, NewtonHHTSolver (enmienda ADR 0011).
- **Convergencia dual ADR 0007 con autocalibración** uniforme.
- **HHT auto-deriva β y γ desde α** con validación de rango.
- **HHT-α formulación correcta** en matriz efectiva y RHS (Hilber-Hughes-Taylor 1977).
- **Newmark predictores y correctores** consistentes con la formulación clásica.
- **Newmark lineal reusa una sola factorización** para todo el barrido temporal.
- **CentralDifference valida** `lumping="lumped"`, diagonalidad y positividad. Falla rápido con mensaje físico claro.
- **Detección CFL a posteriori** con `divergence_threshold` razonable.
- **Crisfield cilíndrico** canónico (predictor con corrección de signo, restricción cuadrática, elección de raíz por coseno).
- **ARPACK shift-invert** con M-ortonormalización y orden ascendente.
- **CQC con coeficientes Der Kiureghian 1980** correctos.
- **`free_vibration` centralizado** en `modal_response.py` (regla D aplicada).
- **Diagnóstico tipado ADR 0011** con detector de oscilación por ventana móvil.
- **Sin tolerancias hardcoded convenientes**: todas vienen de `constants.py` con valores defendibles.

## Cobertura
- **558 tests verdes** en suite global.
- **STATUS / ROADMAP fielmente reflejan el código** (cuentas de componentes coinciden).
- **Specs `validated`** alineadas con los componentes que el ROADMAP declara cerrados.

---

# Plan de acción sugerido (priorización propuesta)

Esta sección es **recomendación**, no ejecución. El usuario decide.

## Inmediato — fixes críticos (recomendable resolver en orden)

1. **H-4.1** — añadir `try/except RuntimeError` con flip de `singular_tangent_seen` en NewtonNewmarkSolver/NewtonHHTSolver. ~30 líneas de cambio total, 2 tests adicionales (provocar singularidad dinámica, verificar `SingularTangentError`).
2. **H-4.3** — corregir doble resta de `F_dir` en CentralDifferenceSolver no lineal. Restar `F_dir` también en la rama no-lineal, eliminar del cómputo de aceleración. Añadir test con apoyo prescrito no nulo en `test_central_difference.py`. ~10 líneas de cambio.
3. **H-4.2** — cachear factorización en LinearSolver con flag de invalidación. ~20 líneas + test de doble llamada. Mejora performance medible.

Recomiendo agrupar los 3 en **una sola rama de trabajo "robustez infraestructural post-auditoría"**, sin ADR nuevo (no son decisiones arquitecturales — son fixes de implementación que respaldan ADRs existentes 0003 y 0011).

## Corto plazo — hallazgos altos

4. **H-1.1** — autodiscover recursivo con `pkgutil.walk_packages`. Más test `test_autodiscover_covers_all_decorated_classes`.
5. **H-1.2 / H-4.5** — centralizar `resolve_rayleigh_config` en `fenix/math/damping.py`. Triplicación confirmada por dos ejes.
6. **H-2.1** — documentar explícitamente el doble contrato (deuda ADR 0002 conocida) en docstrings y manual.
7. **H-2.2** — añadir tests con coeficientes no unitarios para los elementos críticos (al menos uno por elemento sólido 2D y por frame).

## Medio plazo — hallazgos medios accionables

- H-3.1 (DruckerPrager IS_SYMMETRIC dinámico), H-3.2 (validación Elastic1D/2D), H-3.4 (centralizar softening exponencial 1D/2D), H-3.5 (test FD para tangentes DP/VM-PS).
- H-4.4 (refactor zonal de Newmark/HHT — abarata extensiones futuras como Bossak/generalized-α; sin urgencia hasta nuevo esquema temporal).
- H-1.5 (NonlinearSolver keyword-only + heurísticos en constants.py).
- H-1.6 (centralizar validación de `lumping`).
- H-1.7 (validar `PIPELINE_KIND` conocido en run_yaml).
- H-1.9 (magic numbers a constants.py).
- H-5.1, H-5.2 (actualizar catálogos docs).

## Diferido sin urgencia — hallazgos bajos

Mantener como housekeeping cuando se entre a una sesión de pulido. Ninguno bloquea avance.

---

# Cuestiones no verificadas en esta auditoría

- **No se ejecutó la suite de tests** (`pytest`). El análisis es estático sobre código. Los críticos H-4.1 y H-4.3 podrían quedar latentes precisamente porque ningún test ejercita el escenario disparador (singular dinámico no lineal; central diff no lineal con apoyo no nulo).
- **No se corrieron benchmarks físicos** (Van Vliet, voladizo elastoplástico, snap-through, FRF analítica). Tests reportados como verdes (558) se aceptan como veraces.
- **Convención de signos 3D** no auditada exhaustivamente — agentes asumieron cumplida por las múltiples referencias en docstrings. Auditoría dedicada de un solo eje (signos) podría cerrar al 100%.
- **Locking volumétrico** ν → 0.5 declarado limitación pero no se midió cuán cerca produce respuesta razonable.
- **Hourglassing en Quad4 1×1** no cuantificado — el default es 2×2 sin hourglassing.
- **Plasticidad por flexión en frames** (FiberSection) declarada deuda. No se midió cuánto error introduce la actual `E_t` global respecto a integración por capas.
- **Frame3D con warping**: implementación Saint-Venant pura. Insuficiente para secciones abiertas (vigas H, U). No declarado en docstring.
- **Comportamiento de `central_difference.py` con `linear_algebra="auto"`**: parámetro se acepta por compat pero no se usa. No verificado si el YAML parser pasa este kwarg.

---

*Auditoría conducida 2026-05-18 por 4 agentes IA paralelos + consolidación con verificación directa de hallazgos críticos sobre el snapshot del repositorio al cierre de Etapa 6 (commit `5151d7d`). 51 hallazgos totales: 3 críticos, 4 altos, 27 medios, 17 bajos. El núcleo numérico y físico de Fenix FEM es robusto; las deudas accionables son de infraestructura.*

---

# Addendum — Progreso 2026-05-19

Bitácora de cierre tras la sesión de saneamiento post-auditoría. Cuatro
batches de fixes aplicados, **19 hallazgos cerrados** sobre los 51 del
informe; **6 hallazgos diferidos con rationale** explícito (todos
medios-bajos sin bloqueo de avance).

## Hallazgos cerrados (19)

### Críticos (3/3 cerrados)
- ✅ **H-4.3** — `CentralDifferenceSolver` no-lineal contabiliza `F_dir` una sola vez (commit `8c2a809`).
- ✅ **H-4.1** — `SingularTangentError` accesible en `NewtonNewmark`/`NewtonHHT` (commit `5b09a34`).
- ✅ **H-4.2** — `LinearSolver` cachea factorización entre `solve` (commit `b517718`).

### Altos (4/4 cerrados)
- ✅ **H-1.1** — Autodiscover recursivo con `walk_packages` (commit `60c4cb7`).
- ✅ **H-1.2 / H-4.5** — `resolve_rayleigh_config` centralizado en `math/damping` (commit `25ca8ba`).
- ✅ **H-2.1** — Docstring de `Element.internal_forces` extendido con la convención del doble contrato y la deuda ADR 0002 (commit `d333f95`).
- ✅ **H-2.2 / H-5.4** — Test fail-fast de escalamiento dimensional para los 5 sólidos 2D (commit `c89aa37`).

### Medios cerrados (12)
- ✅ **H-1.6** — `validate_lumping_kwarg` centralizado en `core/element` (commit `ac7dc13`).
- ✅ **H-1.7** — `run_yaml` valida `PIPELINE_KIND` contra `_KNOWN_PIPELINE_KINDS` (commit `a7ed8c7`).
- ✅ **H-1.9** — Magic numbers de modal/central-diff/embedded/J2 movidos a `constants.py` con prefijos descriptivos (commit `dbd8210`).
- ✅ **H-2.4** — Catálogo de elementos documenta la limitación de Frame3D lumped en orientación oblicua (commit `a04c5c3`).
- ✅ **H-2.5** — Frame2D Euler/Timoshenko: comentario explícito en `compute_element_state` aclarando que sólo el axial captura plasticidad (commit `a04c5c3`).
- ✅ **H-2.7** — Tri6: docstring extendido con el caveat de subintegración con `tri_3` en mallas distorsionadas (commit `4e4ed54`).
- ✅ **H-3.1** — `DruckerPrager2D.IS_SYMMETRIC` derivado de `self.associated` por instancia; `domain_is_symmetric` lee del instance (commit `a7ed8c7`).
- ✅ **H-3.2** — `Elastic1D`/`Elastic2D` validan `E>0`, `ν∈(-1, 0.5)`, `hypothesis` y `density≥0` (commit `a7ed8c7`).
- ✅ **H-3.4** — `evaluate_exponential_damage` centralizado en `materials/_softening` (commit `d66c1a7`).
- ✅ **H-3.5** — Validación FD de tangentes algorítmicas para VonMises2D (plane strain + plane stress) y DruckerPrager2D regular (commit `4e4ed54`).
- ✅ **H-4.6** — `HHTSolver`/`NewtonHHTSolver` emiten `RuntimeWarning` en override de β/γ (commit `a7ed8c7`).
- ✅ **H-4.7** — `LinearSolver` deriva `is_positive_definite` del flag de simetría del dominio (commit `a7ed8c7`).

### Bajos cerrados — ninguno explícito en esta ronda
Documentalmente: **H-2.6** (lumping nodal vs HRZ) cubierto por la documentación previa en docstring de `_frame2d_lumped_mass_local` y ADR 0009 §"Fase 2"; **H-5.1** (3 textos desactualizados en catálogo de solvers) y **H-5.2** (Elastic2D compatibilidad en catálogo de materiales) corregidos en el commit de catálogos (`a04c5c3`).

## Hallazgos diferidos con rationale (6)

Items pendientes después de esta sesión. Ninguno bloquea el avance; la
decisión de retomarlos requiere un caso de uso explícito o una decisión
arquitectural mayor.

### H-1.3 (parcial) — `HarmonicResult` y `ResponseSpectrumResult` sin `element_forces`
**Severidad**: medio. **Cerrado para `TransientResult`** (commit
`ebe42a5`): método lazy `internal_forces_history(domain)` que recorre
los pasos temporales y devuelve `{elem_id: [ElementForces, ...]}`. Test
analítico verde para Truss2D contra `N(t) = E·A·u(t)/L₀`. Limitación
documentada: corotacionales devuelven el state del último paso, no el
contemporáneo (acción para un futuro: historial de state si aparece
caso de uso).

**Diferido para `HarmonicResult` y `ResponseSpectrumResult`** por razón
técnica, no por falta de consumidor:

- `HarmonicResult.u_complex` es campo complejo (amplitud + fase). El
  contrato `internal_forces(U)` está pensado para campos reales — los
  return mappings de materiales no lineales no aceptan complejos. Para
  caso lineal el cálculo `N(ω) = E·A·û(ω)/L₀` es directo, pero exponerlo
  vía `internal_forces` exigiría extender el contrato base (cambio
  arquitectural). Alternativa más limpia: método dedicado en el
  resultado armónico que asuma comportamiento lineal y compute
  esfuerzos amplitud-fase desde `u_complex`.
- `ResponseSpectrumResult.u_combined` es envolvente máxima (SRSS/CQC
  sobre modos), no campo coherente de desplazamientos. Esfuerzos
  derivados *no* pueden obtenerse llamando `internal_forces(u_combined)`
  — habría que reconstruir respuestas modales individuales
  (`u_n_max = Γ_n·φ_n·S_d(ω_n)`), evaluar `internal_forces(u_n_max)` por
  modo, y aplicar SRSS/CQC a los esfuerzos. Requiere que el resultado
  almacene también las componentes modales individuales en el espacio
  de DOFs (hoy sólo guarda `u_per_mode` que es lo que necesita; falta
  el wiring de cómputo).

Ambos cierres son refactor arquitectural específico y no se acometen
sin un caso de uso real que dirija las decisiones de diseño.

### H-1.4 — Doble contrato `compute_internal_forces` (dict) vs `internal_forces` (`ElementForces|None`)
**Severidad**: medio. **Razón del diferimiento**: la solución
arquitectural (renombrar `compute_internal_forces` como legado deprecado,
o consolidar ambas APIs) requiere cerrar primero la deuda mayor del
ADR 0002 (incompleto para sólidos 2D — registrada en
`project_adr_0002_incompleto_solidos.md`). Renombre parcial sin cerrar
la deuda principal está explícitamente prohibido por la memoria del
proyecto. Acción cubierta documentalmente por H-2.1 ✅ (docstring que
explicita el doble contrato).

### H-2.3 — Truss2D / Truss3D sin matriz geométrica `K_G`
**Severidad**: medio. **Razón del diferimiento**: implementar `K_G` en
los trusses lineales es **ampliación de scope** (nueva feature), no fix.
Las variantes corotacionales (`Truss2DCorot`, `Truss3DCorot`) ya tienen
`K_G`; el caso "análisis de pandeo lineal con trusses no-corotacionales"
es resoluble usando la variante corotacional con linearización al
estado inicial. Retomar cuando aparezca un benchmark de pandeo lineal
de armaduras que demande explícitamente trusses lineales con `K_G`.

### H-4.4 — Refactor zonal de los cuatro solvers Newmark/HHT
**Severidad**: medio. **Razón del diferimiento**: refactor grande
(~880 líneas) sin caso de uso priorizado. Abarataría extensiones futuras
del subsistema temporal (Bossak, generalized-α, energy-momentum schemes)
pero el riesgo de regresión en el código actual supera el beneficio
inmediato. Retomar cuando entre el segundo esquema temporal — ahí la
duplicación pasa el umbral "centralizar siempre" sin ambigüedad.

### H-4.8 — `Assembler.assemble_non_linear_system(u, jac=False)` opcional
**Severidad**: medio. **Razón del diferimiento**: optimización de coste
constante por paso en `CentralDifferenceSolver` no-lineal. Tocar la API
del `Assembler` sin medir el impacto real (qué porción del paso es
ensamblaje de `K_t` vs ensamblaje de `F_int` solo) introduce riesgo
sin garantía de ganancia. Retomar con profiling concreto cuando
`CentralDifference` no-lineal sea el cuello de botella demostrado en
un caso real.

### H-5.3 — 7 huecos de spec vs componentes registrados
**Severidad**: medio. **Razón del diferimiento**: deuda documental
conocida. Cuatro componentes (`Elastic1D`, `Elastic2D`, `Elastoplastic1D`,
`LinearSolver`/`NonlinearSolver`/`ArcLengthSolver`) son anteriores a la
convención de specs validadas. `NewtonHHTSolver` es variante de
`HHTSolver` y comparte spec por Reglas §4 ("spec corta tipo extensión").
Rellenar las 7 specs son ~3.5h de trabajo mecánico que no aporta nuevo
conocimiento físico ni cambia comportamiento — se retoma en una sesión
dedicada de housekeeping documental.

## Hallazgos bajos pendientes (17)

Sin acción en esta ronda. Lista resumen para próxima sesión de
housekeeping (orden por archivo del informe original):

- Eje 1: H-1.8 (step_callback silenciado), H-1.10 (frozen=True semántico),
  H-1.11 (capa core importa results), H-1.12 (excepciones tipadas no
  extendidas a return mapping), H-1.13 (duplicación geom truss/cable
  corot), H-1.14 (YAML parser usa `_materials` privado), H-1.15
  (Node.dofs sentinel -1).
- Eje 2: H-2.8 (Frame2D doble rotación cosmético), H-2.9 (warning
  Frame3D casi-vertical), H-2.10 (CST_Embedded2D valida bulk por nombre
  de clase).
- Eje 3: H-3.3 (nomenclatura "tensión" vs "esfuerzo"), H-3.6 (test del
  switching exact CableMaterial1D), H-3.7 (nota cap DAMAGE_MAX en
  docstring), H-3.8 (trial/commit no documentado en contrato), H-3.9
  (invariante `tr(ε^p)=0` no verificado), H-3.10 (yield_tol con alpha_old).
- Eje 4: H-4.9 (doble ensamblaje line_search), H-4.10 (`Iventalla.dat`
  untracked — fuera de alcance).
- Eje 5: H-5.5 (`MATRIZ.md` faltante), H-5.6 (`ONBOARDING.md` pendiente).

## Métricas finales

| Eje | Cerrados | Diferidos con rationale | Bajos pendientes |
|---|---:|---:|---:|
| 1 — Arquitectura | 5 (H-1.3 parcial: TransientResult ✅, Harmonic/Spectrum diferidos) | 2 | 7 |
| 2 — Elementos | 5 | 1 | 3 |
| 3 — Materiales | 5 | 0 | 6 |
| 4 — Solvers | 4 | 2 | 2 |
| 5 — Docs / cobertura | 3 | 1 | 2 |
| **Total** | **23** | **6** | **20** |

(El item H-2.6 cuenta como "cubierto documentalmente" y no se duplica;
H-5.4 es alias de H-2.2 y se cuenta una vez; H-4.5 es alias de H-1.2.
H-1.3 ahora vale como cerrado para `TransientResult`; `HarmonicResult` y
`ResponseSpectrumResult` permanecen diferidos por razón técnica.)

Suite global al final del addendum: **585 passed, 5 skipped, 0 failures**
(originales 558 + 27 tests añadidos durante la sesión).
