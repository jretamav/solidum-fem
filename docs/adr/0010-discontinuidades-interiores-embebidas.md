# ADR 0010 — Discontinuidades interiores embebidas

- **Estado**: Aceptado
- **Fecha**: 2026-05-13
- **Alcance**: nueva familia paralela `CohesiveMaterial` (clase base independiente de `Material`); nueva subjerarquía de elementos con DOFs enriquecidos a nivel elemental y condensación estática local; nueva semántica de "elemento que cambia de modo" durante la simulación; catálogos de elementos y materiales; manuales.

## Terminología

- **Embedded discontinuity (ED-FEM)**: familia de métodos donde un elemento finito alberga internamente una discontinuidad de desplazamientos `[[u]]` sobre una superficie `Γ_d` que cruza el elemento. La continuidad del salto entre elementos vecinos **no** se impone (a diferencia de XFEM). Los DOFs del salto son **elementales**, no globales — se condensan estáticamente antes de devolver `K_e` y `f_e` al ensamblador.
- **Discrete approach** (Retama 2010): el bulk se modela con una ley continua `σ–ε` elástica; la discontinuidad con una ley cohesiva `t–[[u]]` con softening. Movimiento de cuerpo rígido relativo entre las dos partes del elemento. Esta es la aproximación adoptada.
- **Continuum approach** (alternativa no adoptada): una sola ley continua en todo el dominio, con la grieta apareciendo como un caso límite (`weak/strong discontinuity`). Queda fuera de este ADR.
- **Clasificación de formulaciones (Jirásek 2000)**: SOS (Statically Optimal Symmetric), **KOS (Kinematically Optimal Symmetric)**, SKON (Statically/Kinematically Optimal Nonsymmetric). La **fase 1 de Fenix implementa KOS**, fiel a la formulación principal de Retama (2010). KSON queda diferida (fase H).
- **Traction-jump material**: material constitutivo que opera sobre el salto `[[u]] ∈ ℝ²` y devuelve tracciones `t ∈ ℝ²` definidas en `Γ_d`. Es la familia análoga a `Material` (que opera sobre `ε` y devuelve `σ`) pero con dimensión, escala y unidades distintas. Sus parámetros típicos son `σ_t0` (resistencia a tracción) y `G_F` (energía de fractura), no `E, ν`.
- **`DiscontinuityState`**: estructura de datos elemental que recoge `(n, l_d, [[u]], κ, ω)` para un elemento agrietado. Análoga a `ElementState` pero específica de la discontinuidad.
- **Activación**: transición de un elemento desde "intacto" (sin discontinuidad, comportamiento estándar) a "agrietado" (con discontinuidad activa, modos enriquecidos contribuyendo). Una vez activado, no revierte.

## Resumen ejecutivo

Fenix abre la línea de fractura computacional adoptando **discontinuidades interiores embebidas en aproximación discreta** (Retama 2010), por las razones discutidas con el usuario: maximiza el aprovechamiento de su experticia doctoral y construye abstracciones que el proyecto necesitará tarde o temprano (materiales cohesivos, DOFs enriquecidos elementales, condensación estática local). La fase 1 implementa **KOS sobre CST 2D, Modo-I, fiel a la tesis 2010**, incluyendo explícitamente la aportación del Cap. 6 (`l_d = (A/h)·cos(θ−α)`) y validándola por primera vez con el experimento numérico que la tesis no contiene. Este ADR fija cinco decisiones arquitecturales cross-cutting (familia `CohesiveMaterial` paralela, elementos con DOFs enriquecidos elementales, condensación estática local, estado de discontinuidad, semántica de activación) y deja registrada la hoja de ruta hacia extensiones posteriores (sin tracking estilo Brank/Zhang, modo mixto, KSON, 3D, orden superior). Las decisiones se toman ahora para que las fases posteriores no reabran debates.

## Contexto

Hoy Fenix tiene una línea de daño continuo madura (`IsotropicDamage1D/2D`, softening exponencial, tangente algorítmica consistente). Lo que no tiene — y lo que la línea de embedded discontinuity introduce — son **discontinuidades físicas reales** del campo de desplazamientos: grietas que abren, descargan el bulk vecino, dejan de transmitir tracción más allá de un umbral. Tu tesis doctoral (UNAM 2010, Cap. 8) muestra que esto es esencial en materiales cuasi-frágiles donde el bulk descarga elásticamente mientras la grieta dispara la disipación.

**Relación con ADRs previos**:
- **ADR 0002** (API de resultados): `SolveResult` admite campos por elemento sin cambios. El `DiscontinuityState` se expone como un campo elemental en post-proceso (`crack_opening`, `damage_at_crack`, `crack_normal`).
- **ADR 0003** (despachador algebraico): no afectado. Las matrices que llegan al ensamblador son ya las condensadas `K_cond`. Para KOS son simétricas, así que el path Cholesky/LU sigue eligiéndose igual. KSON (fase H) introduciría matrices asimétricas; entonces se evaluará si forzar LU global o mantener mixto.
- **ADR 0006/0007** (tolerancias): el Newton-Raphson global no cambia. La tolerancia adimensional ya parametrizada (atol + rtol·escala física) sirve sin retocar.
- **ADR 0008** (density): no aplica — los materiales cohesivos no tienen densidad propia (la inercia es del bulk).
- **ADR 0009** (modal y dinámico): tampoco afectado en fase 1 (sólo cuasi-estática). Una extensión futura a dinámica con embedded discontinuity se acomoda sin cambios estructurales por la decisión 5 de este ADR (estado de la discontinuidad fuera de `Node`).

**Reubicación dentro del plan**: la decisión de implementar esto ya hoy se toma con autorización explícita del usuario, por encima de la Etapa 5 pendiente del ROADMAP (sólidos 3D / placas / térmico / completar ADR 0009 / Mohr-Coulomb). El ROADMAP se actualizará tras aceptación de este ADR.

## Decisión

### 1. Familia `CohesiveMaterial` paralela a `Material`, no subclase

Se introduce una jerarquía nueva, paralela e independiente, para materiales que operan sobre saltos en una superficie y devuelven tracciones:

```python
class CohesiveMaterial(ABC):
    """Material cohesivo traction-jump. Define t = T([[u]], state) en Γ_d.

    Distinto en dominio, dimensión y unidades de Material:
      - Material: σ = C·ε,           ε ∈ ℝ^{STRAIN_DIM},  σ ∈ ℝ^{STRAIN_DIM},  [σ] = Pa.
      - CohesiveMaterial: t = T·[[u]],  [[u]] ∈ ℝ^{ndim},  t ∈ ℝ^{ndim},        [t] = Pa.
                          (parámetros típicos: σ_t0 [Pa], G_F [N/m].)
    """
```

**Razón de no unificar con `Material`**:

- ε es un tensor de deformación (adimensional), [[u]] es un vector de desplazamiento (longitud). Sus unidades, sus operaciones algebraicas (gradiente vs evaluación directa) y sus escalas físicas son distintas.
- σ es esfuerzo (fuerza/área); t también lo es **pero sobre una superficie**, con escala asociada a `G_F` (fuerza/longitud, energía por unidad de área de grieta) — no a `E` (rigidez volumétrica).
- Los tests de aceptación son fundamentalmente distintos: para `Material` se cicla `ε`; para `CohesiveMaterial` se aplica un historial de `[[u]]` y se verifica que la integral `∫t·d[[u]] = G_F` cierre la grieta correctamente.
- Reusabilidad: `CohesiveMaterial` también lo consumirán **elementos de interfaz cohesiva** (CZM clásico entre dos cuerpos) si entran en el futuro — no sólo embedded discontinuity. La abstracción tiene vida propia.

**Registro**: se reutiliza el patrón de auto-registro vía decorador, pero con su propio registry — `CohesiveMaterialRegistry`, paralelo a `MaterialRegistry`. Razón: mezclar continuos y cohesivos en un único registry obligaría al parser YAML a discriminar por tipo en cada uso, ensuciando el contrato.

### 2. Elementos con DOFs enriquecidos elementales

Los DOFs del salto `[[u]] ∈ ℝ²` viven como **atributo del elemento**, no como nodos del modelo global. El contrato del `Element` se preserva: el ensamblador llama a `compute_stiffness_and_internal_forces()` y recibe `(K_e, f_e)` de dimensiones estándar (`n_nodos · n_DOFs_por_nodo`); el elemento se encarga internamente de la condensación.

```python
class CST_Embedded2D(Element):
    discontinuity_state: DiscontinuityState | None  # None = elemento intacto

    def compute_stiffness_and_internal_forces(self):
        if self.discontinuity_state is None:
            return self._compute_standard()        # CST regular
        return self._compute_enriched_condensed()  # con condensación local
```

**Justificación**: el ensamblador no necesita saber que el elemento tiene modos enriquecidos. Esta es exactamente la promesa de la condensación estática local (Retama 2010 Cap. 7): los modos elementales no llegan al sistema global. Preserva todas las invariantes del `Assembler` actual — topología COO cacheada, despachador algebraico, expansión Dirichlet — sin tocar nada.

### 3. Condensación estática local en el elemento

Implementación canónica (Retama 2010 §7.1):

```
Sistema local antes de condensar:
  [K_ûû  K_ûũ] {dû}   {R_û}
  [K_ũũ^T K_ũũ] {dũ} = {R_ũ}

K_cond = K_ûû - K_ûũ · K_ũũ^{-1} · K_ûũ^T          ← devuelta al ensamblador
R_cond = R_û  - K_ûũ · K_ũũ^{-1} · R_ũ              ← devuelto al ensamblador

Tras resolución global y obtener dû:
  dũ = K_ũũ^{-1} · (R_ũ - K_ũũ^T · dû)              ← recuperación local
```

La variante alternativa de Retama 2010 §7.2 (almacenar sólo `ũ` actual, re-evaluar `K_ũũ` por iteración) queda **diferida**. La elegirá una optimización posterior si la memoria del estado se vuelve problema en mallas grandes — no es decisión arquitectural sino implementacional.

### 4. `DiscontinuityState` — estado del elemento agrietado

Dataclass paralela a `ElementState`, específica de la discontinuidad:

```python
@dataclass
class DiscontinuityState:
    normal: np.ndarray        # n ∈ ℝ², perpendicular a σ_I al activarse
    centroid: np.ndarray      # punto por el que pasa Γ_d (centroide en CST)
    l_d: float                # (A/h)·cos(θ−α), Cap. 6 Retama 2010
    jump_committed: np.ndarray   # [[u]] al cierre del último paso convergido
    jump_trial: np.ndarray       # [[u]] dentro de la iteración Newton
    kappa: float              # max ⟨[[u]]_eq⟩ histórico
    damage: float             # ω(κ) actual
```

Semántica **trial / commit** consistente con `ElementState` actual: `jump_trial` evoluciona en cada iteración Newton; `jump_committed` se actualiza sólo cuando el paso converge. Mismo patrón que en plasticidad J2 actual.

### 5. Activación — semántica de elemento que cambia de modo

Antes de cada paso del solver no lineal:

```
para cada elemento intacto:
    σ_I = max principal stress en el elemento (centroide para CST)
    si σ_I > σ_t0:
        n = autovector de σ_I (Rankine: perpendicular a la tracción máxima)
        l_d = (A/h)·cos(θ−α)  con θ, α derivados de n y la geometría
        instanciar DiscontinuityState
        marcar elemento como agrietado (irreversiblemente)
```

A partir de esa activación, el elemento contribuye con su `K_cond` enriquecida.

**El elemento no cambia de clase** — sigue siendo `CST_Embedded2D` antes y después. Lo que cambia es su estado interno (`discontinuity_state: None → instancia`). Esto evita "elementos que mutan de tipo" como concepto arquitectural (que rompería invariantes de assembly, registry, post-proceso).

La verificación de activación se ejecuta **al principio del paso, no dentro del Newton local**. Razón: dentro del Newton, σ_I oscila por el predictor lineal y puede cruzar y descruzar el umbral varias veces antes de converger; activar dentro del Newton lleva a chattering. Al principio del paso, σ_I se evalúa con el estado convergido del paso anterior — decisión estable.

### 6. `l_d = (A/h)·cos(θ−α)` — aportación del Cap. 6, decisión cerrada

La longitud de la discontinuidad se calcula con la fórmula de Retama (2010) §6.1:

```
l_d = (A / h) · cos(θ − α)
```

donde `A` es el área del CST, `h` la altura desde el nodo solitario, `θ` el ángulo de la grieta con respecto a `x` global, `α` el ángulo del lado opuesto al nodo solitario. La derivación geométrica está en la tesis y demuestra que ésta es **la única longitud que satisface continuidad fuerte de tracciones** sobre `Γ_d` para CST con KOS.

Esta fórmula es **decisión cerrada en el ADR** — no se relitiga en specs ni en fases posteriores. La alternativa ingenua `l_d = A/h` (paralela al lado opuesto) sólo coincide cuando la grieta cae exactamente paralela a un lado; en cualquier otra orientación introduce stress locking. Esto se verifica numéricamente en la fase 3 (ver hoja de ruta) como test de aceptación del elemento — el experimento que la tesis no incluyó.

### 7. Tratamiento en YAML

```yaml
# Material continuo del bulk (existente)
materials:
  - name: hormigon_elastico
    type: Elastic2D
    parameters: { E: 39.8e9, nu: 0.2, mode: plane_stress }

# Material cohesivo de la discontinuidad (nuevo)
cohesive_materials:
  - name: concreto_cohesivo
    type: CohesiveDamageIsotropic
    parameters:
      sigma_t0: 2.57e6      # Pa, resistencia a tracción
      G_f: 121.9            # N/m, energía de fractura
      softening: exponential  # linear | exponential
      penalty_stiffness: 1.0e15  # rigidez elástica del salto, Retama 2010 §3.1

# Elemento con discontinuidad embebida (nuevo)
elements:
  - type: CST_Embedded2D
    nodes: [...]
    bulk_material: hormigon_elastico
    cohesive_material: concreto_cohesivo
    activation_criterion: rankine  # único soportado fase 1
```

Decisión: `bulk_material` y `cohesive_material` se declaran por separado en el elemento. El bulk se evalúa siempre; el cohesivo entra al sistema sólo tras activación. El parser YAML acepta una sección nueva `cohesive_materials` para no mezclar con la lista de continuos (consecuencia directa de la decisión 1).

## Hoja de ruta — fases

| Fase | Componente | Entrega | Estado |
|---|---|---|---|
| 1 | Spec + implementación `CohesiveDamageIsotropic` (Mode-I, daño isótropo, softening lineal/exponencial) | Material testeado en aislamiento con historial prescrito de `[[u]]`; verifica `∫t·d[[u]] = G_F` | Completado (commit `c7e68de`, 2026-05-18) |
| 2 | Spec + implementación `CST_Embedded2D` (KOS, fiel a Retama 2010 Caps. 2, 5, 7) | Elemento aislado en tracción uniaxial reproduce curva σ–ε esperada | Completado (commit `2e85a70`, 2026-05-18) |
| 3 | Validación del Cap. 6 — test `l_d ingenuo` vs `l_d = (A/h)·cos(θ−α)` | Test integrado en spec del elemento documenta stress locking en la versión ingenua | Completado (`tests/test_ld_chapter_6_validation.py`, 2026-05-18) |
| 3b | Bring-up de integración end-to-end + caza del bug dimensional `thickness` en términos cohesivos | 4 tests de cableado + 4 tests de regresión dimensional (`TestThicknessDimensionalConsistency`); ADR §"Caveats y lecciones aprendidas" | Completado (commit `e54b5ce`, 2026-05-18) |
| 4 | Benchmark Van Vliet faithful (Retama 2010 §8.1) — comparativa con la curva experimental | Test end-to-end, ambas leyes de softening, rama de descarga completa contra experimental | **Diferida**. Bloqueante: los solvers actuales (`NonlinearSolver`, `ArcLengthSolver` cilíndrico) no atraviesan la rama post-pico del embedded discontinuity con penalty `K_e` stiff (transición elástico→softening abrupta hace `K_t` indefinida). El modelo físico está verificado en aislamiento por los tests del subsistema (incluido el blindaje dimensional). El límite es solver, no formulación. Retoma: cuando se priorice el mini-ADR 0011 de solvers para softening severo (dissipation arc-length, CMOD/CTOD control vía MPC con reacción en el cabezal, sign-of-pivot tracking). |
| 5 | Catálogos + manuales | Entradas en `catalogo_elementos.md` y `catalogo_materiales.md` añadidas en fases 1 y 2; manuales auto-regeneran del MD existente | Cubierto parcialmente; regeneración explícita de manuales con la sección "Materiales cohesivos / Elementos embedded" diferida con la fase 4. |
| F | Cracking-elements-style sin crack tracking (Brank 2021, Zhang 2018+) | Variante del elemento que elimina la necesidad de tracking trivial | Diferida |
| G | Modo mixto I-II | Extensión del material cohesivo con componente tangencial | Diferida |
| H | KSON nonsymmetric | Implementación de la variante no simétrica para comparativa académica con KOS | Diferida |
| I | 3D embedded discontinuity (Tet4_Embedded, etc.) | Cuando entren sólidos 3D al proyecto | Diferida |
| J | Orden superior (Tri6_Embedded, Quad8_Embedded) | Posterior a fase I | Diferida |

Cada fase es un commit cerrado con tests verdes, no bloquea las siguientes y deja Fenix en estado funcional.

## Consecuencias

**Positivas**

- Fenix gana fractura computacional con una formulación que el usuario domina a nivel doctoral. Coste de aprendizaje cero.
- Las cuatro abstracciones nuevas (`CohesiveMaterial`, elementos con DOFs enriquecidos elementales, condensación estática local, `DiscontinuityState`) son **inversiones reutilizables**: sirven a futuros CZM clásicos, futuros incompatible-modes, futuras formulaciones EAS, eventual XFEM-style si llega.
- La aportación original de la tesis (`l_d` correcto) queda implementada, validada numéricamente por primera vez, y atribuida explícitamente en el ADR y los tests.
- Se preservan todos los contratos existentes: `Assembler`, despachador algebraico, Newton-Raphson, parser YAML, sistema de manuales. La cirugía está localizada en clases nuevas y un atributo nuevo en `Element` base.

**Negativas / costes**

- Subsistema grande. Esfuerzo de fase 1 estimado en 4–6 sesiones (material, elemento, validación Cap. 6, Van Vliet, catálogo). No es trabajo de fin de semana.
- Introduce el concepto de "elemento que cambia de modo" durante la simulación. Aunque se contiene como bandera interna (decisión 5), conceptualmente es semántica nueva en Fenix.
- La condensación estática local añade un Newton interno por elemento agrietado (resolución de `K_ũũ · dũ = ...`). Para mallas con muchos elementos agrietados simultáneos el coste es no trivial; mitigación habitual: pocos elementos activos en cada paso (típico en problemas de fractura localizada).
- Fase 1 hereda dos limitaciones reconocidas en la tesis: sólo Modo-I (modo mixto en fase G) y tracking trivial heredado de la dirección principal (sin tracking en fase F). Se documentan como `out_of_scope` en las specs.

**Alternativas consideradas**

- *`CohesiveMaterial` como subclase de `Material` con `STRAIN_DIM` reinterpretado*: rechazada. La semántica (operar sobre `[[u]]` no sobre `ε`; unidades de `G_F` no de `E`; tests sobre integral de tracción no sobre ciclo de deformación) es lo suficientemente distinta para que la unificación introdujese más confusión que ahorro de código.
- *Cohesive zone model clásico (elemento de interfaz entre dos cuerpos) en lugar de embedded discontinuity*: rechazada para fase 1, no para el futuro. CZM requiere remallado para introducir interfaces — la promesa del embedded discontinuity es precisamente evitarlo. El `CohesiveMaterial` definido aquí es reutilizable por un CZM futuro si entra al proyecto.
- *Phase-field fracture como línea alternativa*: rechazada en la discusión previa con el usuario. Es la corriente dominante hoy en fractura computacional pero requiere expertise que no aprovecha la tesis y un campo adicional con resolución fina (coste computacional alto). Phase-field como línea propia es decisión futura, independiente de ésta.
- *Empezar por la versión sin crack tracking (Brank 2021, Zhang cracking elements)*: rechazada por decisión explícita del usuario — fiel a la tesis 2010 primero, validar el trabajo original, después extender.
- *Implementar KOS y KSON en paralelo en fase 1 para replicar la comparativa de la tesis*: rechazada — infla el alcance ~30% sin valor inmediato. KSON entra como fase H si se quiere reproducir la comparativa académica.

## Caveats y lecciones aprendidas

### Consistencia dimensional de los términos cohesivos (descubierto 2026-05-18)

**Síntoma**: durante un ensayo *diagnóstico* del subsistema (probeta Van Vliet en tracción centrada — no la versión faithful con excentricidad de la fase 4), la curva carga-deformación nunca alcanzaba pico. El cohesivo no degradaba aunque `σ_bulk` superase `σ_t0` en los elementos del crack. En diagnóstico por elemento: `jump_n ≈ 3×10⁻¹⁰ m` (10× menor que el equilibrio físico requeriría, ~3×10⁻⁹ m), `damage = 0` siempre, y la tracción cohesiva era `t_n ≈ σ_bulk/10`.

**Causa raíz**: en `CST_Embedded2D._element_state_cracked`, el residual y la rigidez locales del Newton del salto sumaban dos términos en **unidades diferentes**:

- Bulk: `-G^T · B^φ^T · σ · vol` con `vol = A_e · thickness` → unidades de N (fuerza 3D).
- Cohesivo: `l_d · t` con `l_d` en m y `t` en Pa → unidades de N/m (fuerza por unidad de espesor).

Sumar `N + N/m` produce un residuo numéricamente válido pero físicamente inconsistente: el cohesivo entraba sobre-pesado por `1/thickness` en el balance local. Newton convergía a un salto `thickness×` **menor** del correcto, y el cohesivo nunca alcanzaba `κ_0` para activar damage.

**Fix**: multiplicar los términos cohesivos por `self.thickness` para cerrar la consistencia dimensional:

```python
R_jump = -G.T @ (B_phi.T @ sigma) * vol + ds.l_d * t_local * self.thickness
K_jj   = G.T @ (B_phi.T @ C_tan_bulk @ B_phi) @ G * vol + ds.l_d * T_coh * self.thickness
```

**Por qué fue invisible**: todos los tests existentes del subsistema (fases 1–3 + bring-up de integración + 24 tests del elemento) usaban `thickness = 1.0`. Con espesor unitario, el factor faltante es la identidad y el bug no se manifiesta. El primer caso real con `thickness ≠ 1.0` fue Van Vliet (`thickness = 0.1 m`).

**Blindaje**: cuatro tests nuevos en `tests/test_cst_embedded.py::TestThicknessDimensionalConsistency`, todos ejercitan `thickness = 0.1`:
1. `test_cohesive_engages_with_non_unit_thickness` — verifica que `damage > 0` cuando la carga rebasa el umbral.
2. `test_traction_balances_bulk_stress` — verifica continuidad `t_n ≈ σ·n` con tolerancia `1e-3 · σ_t0`.
3. `test_residual_with_correct_formula_vanishes` — reconstruye `R^{[[u]]}` con el factor `thickness` correcto y confirma anulación.
4. `test_jump_scales_independently_of_thickness` — el salto físico es invariante bajo cambio de `thickness` con cinemática idéntica (antes del fix, escalaba con `thickness`).

El test antiguo `TestLocalRecovery.test_residual_vanishes_after_loading` también se actualizó para usar la fórmula corregida (era *accidentalmente* correcto con `thickness=1.0` y habría dado falso negativo con cualquier otro espesor).

**Lección general**: tests de unidad con coeficientes "1.0" (espesores unitarios, módulos unitarios, áreas unitarias) son **trampa de cobertura**. Validan la lógica simbólica del operador pero hacen invisible cualquier inconsistencia dimensional. Política para el subsistema: tests con al menos un parámetro físico de escala ≠ 1.0 cuando ese parámetro aparezca en la formulación.

## Referencias

- **Retama Velasco, J. (2010).** *Formulation and Approximation to Problems in Solids by Embedded Discontinuity Models*. Tesis Doctoral, Instituto de Ingeniería, UNAM. Director: Dr. A. Gustavo Ayala Milián. **Referencia primaria; el autor de la tesis es el usuario de Fenix FEM.** PDF en [`docs/referencias/PhD_Thesis_Retama_2010.pdf`](../referencias/PhD_Thesis_Retama_2010.pdf).
- Jirásek, M. (2000). Comparative study on finite elements with embedded discontinuities. *Computer Methods in Applied Mechanics and Engineering* 188, 307-330. — Clasificación SOS/KOS/SKON adoptada.
- Oliver, J., Huespe, A.E., Sánchez, P.J. (2006). A comparative study on finite elements for capturing strong discontinuities. *CMAME* 195, 4732-4752.
- Linder, C., Armero, F. (2007). Finite elements with embedded strong discontinuities for the modeling of failure in solids. *International Journal for Numerical Methods in Engineering* 72, 1391-1433.
- Brank, B., Stanić, A., Ibrahimbegović, A. (2021). Crack propagation simulation without crack tracking algorithm: embedded discontinuity formulation with incompatible modes. *CMAME* 386. — Referencia para la fase F futura.
- Zhang, Y. et al. (2018-2024). Cracking elements method y derivados. — Referencia para la fase F futura.
- ADR 0002 — API de resultados (`SolveResult` admite `DiscontinuityState` como campo elemental).
- ADR 0003 — Despachador algebraico (no afectado en fase 1; KSON fase H podría requerir LU global).
- ADR 0007 — Convergencia y tolerancias (Newton global hereda sin cambios).
