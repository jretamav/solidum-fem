# Estado de Fenix FEM

> Foto de una página del estado actual. **Se regenera o actualiza tras cada PR significativa** (cierre de etapa, nuevo componente validado, nuevo ADR aceptado, deuda saldada).
>
> Para visión por etapas y bifurcaciones: ver [ROADMAP.md](ROADMAP.md). Para combinaciones validadas: ver `MATRIZ.md` (pendiente).

---

## Métricas

| Indicador | Valor |
|---|---|
| **Tests** | 420 pasan / 5 skipped / 0 fallos (425 colectados) |
| **Elementos** | 16 (10 estructurales 1D + 5 sólidos 2D + 1 sólido 2D con discontinuidad embebida) |
| **Materiales** | 9 (8 continuos + 1 cohesivo traction-jump) |
| **Solvers** | 6 (3 estáticos + 1 modal + 2 transitorios) |
| **ADRs aceptados** | 10 (0001–0010) |
| **Specs `validated`** | 25 |
| **Etapas cerradas** | 4 completas + 1 parcial (etapa 4 vía ADR 0009 fases 1, 3, 4) |

---

## Capacidad de análisis hoy

**Análisis estáticos**
- Lineal (`LinearSolver`).
- No lineal con control de carga (`NonlinearSolver`, Newton-Raphson incremental con paso adaptativo).
- Trazado de curvas de equilibrio con snap-through/snap-back (`ArcLengthSolver`, Crisfield cilíndrico).

**Análisis dinámicos**
- Modal — autovalores generalizados con shift-invert ARPACK (`ModalSolver`).
- Transitorio lineal Newmark con amortiguamiento Rayleigh, β-γ parametrizables (`NewmarkSolver`).
- Transitorio no lineal Newton-Newmark, mismo amortiguamiento (`NewtonNewmarkSolver`).

**Geometrías cubiertas**
- Estructuras 1D: barras (`Truss2D/3D` lineales y corotacionales), cables (`Cable2D/3D` corotacionales), vigas (`Frame2DEuler`, `Frame2DTimoshenko`, `Frame2DEulerCorot`, `Frame3D`).
- Sólidos 2D: lineales (`Quad4`, `Tri3`) y cuadráticos (`Quad8`, `Quad9`, `Tri6`).
- **Sólidos 2D con discontinuidad interior embebida** (Etapa 5, ADR 0010): `CST_Embedded2D` — CST con cinemática KOS enriquecida, DOFs del salto condensados a nivel elemental, activación Rankine en `prepare_step` con `l_d = (A_e/h)·cos(θ−α)` (Cap. 6 Retama 2010).

**Materiales cubiertos**
- Elásticos: `Elastic1D`, `Elastic2D` (plane stress + plane strain).
- Cable unilateral: `CableMaterial1D`.
- Plasticidad J2: `Elastoplastic1D`, `VonMises2D` (kernels Numba especializados para plane strain y plane stress).
- Plasticidad friccional: `DruckerPrager2D` (plane strain, dos ramas regular/apex, plasticidad no asociada).
- Daño isótropo con softening exponencial: `IsotropicDamage1D`, `IsotropicDamage2D` (tangente algorítmica consistente, asimétrica en 2D).
- **Material cohesivo traction-jump** (familia paralela `CohesiveMaterial`, ADR 0010): `CohesiveDamageIsotropic` — Mode-I, daño isótropo con softening lineal/exponencial parametrizado por `σ_t0` y `G_F`. Penalty `K_e` separado de la rigidez del bulk.

**Infraestructura común**
- Auto-registro vía decoradores en `ElementRegistry`, `MaterialRegistry`, `SolverRegistry`.
- Despachador algebraico automático Cholesky/LU según simetría/SPD (ADR 0003).
- Caché de topología COO compartida entre `K`, `M`, masa, peso propio.
- Tolerancias adimensionales con escala física por material (ADR 0006) y convergencia dual desplazamiento + residuo (ADR 0007).
- Parser YAML con despacho automático estático / modal / transitorio.
- Tres manuales `.tex/.pdf` generados desde markdown: Reference (FF-MR), User (FF-MU), Architecture (FF-MA).

---

## Limitaciones declaradas (out-of-scope hoy)

- **Sólidos 3D**: pendiente (`Hex8`, `Tet4`, etc.).
- **Placas y láminas**: pendiente.
- **Análisis térmico**: pendiente (decoupled → coupled).
- **Mohr-Coulomb 2D**, **FiberSection para frames no-lineales**: pendientes.
- **HHT-α, central differences, harmonic, response spectrum**: fases 5–7 del ADR 0009 diferidas.
- **Contacto mecánico**: horizonte largo.
- **Grandes deformaciones en sólidos** (lagrangiano total/actualizado): horizonte largo.
- **Hiperelasticidad**, **plasticidad anisótropa**, **daño con regularización** (gradient damage, phase-field): horizonte largo.
- **Locking volumétrico** en `Quad4`/`Quad8` con ν → 0.5: declarado limitación, mitigaciones B-bar/F-bar diferidas.
- **Plasticidad por flexión en frames**: hoy sólo plastifican axialmente (espera `FiberSection`).
- **Trazado completo de la rama post-pico con embedded discontinuity**: `NonlinearSolver` y `ArcLengthSolver` cilíndrico no atraviesan la transición elástico→softening con penalty cohesivo stiff (`K_t` casi singular en `κ_0`). El modelo físico está verificado en aislamiento (material cohesivo + elemento + condensación); el límite es solver. Retoma vía mini-ADR de "solvers para softening" (dissipation arc-length, control CMOD/CTOD, sign-of-pivot tracking) cuando se priorice.
- **Modo mixto I-II en cohesivo**, **KSON no-simétrico**, **embedded sin tracking**, **3D embedded**, **orden superior**: fases F-J del ADR 0010 diferidas.

---

## Deuda técnica priorizada

| # | Item | Criterio de retoma |
|---|---|---|
| 1 | `internal_forces` devuelve `None` en sólidos 2D (ADR 0002 incompleto). | Cuando entren sólidos 3D, post-proceso avanzado o consumidor externo que pida `ElementForces` para sólidos. |
| 2 | `FiberSection` para plasticidad por flexión en frames. | Si la próxima etapa se decide por la opción E (Mohr-Coulomb + FiberSection). |
| 3 | Reglas de disparo C y D arquitecturales pendientes ([memoria](../../../.claude/projects/g--Mi-unidad-Proyectos-IA-fenix-fem/memory/project_reglas_disparo_pendientes.md)). | Cuando ocurra el evento que cada regla espera. |
| 4 | Solver para softening severo con embedded discontinuity (fase 4 ADR 0010 + curva descarga Van Vliet). | Cuando se priorice un mini-ADR de "solvers para softening": dissipation arc-length, control CMOD/CTOD del crack, sign-of-pivot tracking. |

Ninguno de los tres bloquea el avance. Todos están documentados con su contexto en la memoria del proyecto o en [MATRIZ.md](MATRIZ.md), y no requieren acción proactiva.

**Deuda saldada 2026-05-14**: cobertura asimétrica de sólidos 2D no-Quad4 con materiales no lineales — cerrada con [`test_solid_2d_nonlinear_higher_order.py`](../tests/test_solid_2d_nonlinear_higher_order.py) (16 tests verdes a la primera; no había bug latente).

**Etapa cerrada 2026-05-18**: discontinuidades interiores embebidas (Etapa 5, [ADR 0010](adr/0010-discontinuidades-interiores-embebidas.md)). Subsistema completo de fractura computacional vía embedded discontinuity discreta, fiel a Retama 2010: material cohesivo `CohesiveDamageIsotropic`, elemento `CST_Embedded2D` con condensación local, validación numérica del `l_d` del Cap. 6 y bring-up end-to-end. Bug arquitectural cazado durante el bring-up (factor `thickness` faltante en términos cohesivos, invisible con tests de espesor unitario) documentado en el ADR §"Caveats y lecciones aprendidas" y blindado con cuatro tests de regresión. Fases 4 (benchmark Van Vliet faithful) y 5 (manuales explícitos) diferidas — modelo físico verificado por curva analítica 1D, límite en solver no en formulación.

---

## Próximo hito

**Elección de Etapa 6** entre las opciones diferidas A-E del ROADMAP:
- A. Sólidos 3D (`Hex8`, `Tet4`, …).
- B. Placas y láminas.
- C. Análisis térmico desacoplado.
- D. Completar ADR 0009 (HHT-α, central differences, harmonic).
- E. Mohr-Coulomb 2D + `FiberSection`.

El argumentario completo de cada opción está en [ROADMAP.md](ROADMAP.md). La decisión la toma el usuario con base en la dirección que quiera dar al proyecto tras esta etapa.

---

## Cómo se regenera este documento

1. Actualizar la tabla de métricas contando ADRs, specs, componentes registrados.
2. Si entró un componente nuevo: añadir a "Capacidad de análisis hoy" en la sección correspondiente.
3. Si saldó deuda: tachar el item de la tabla con `~~strikethrough~~` y mover la línea bajo el cuerpo o eliminarla en la siguiente actualización.
4. Si cambió el próximo hito: reescribir la sección "Próximo hito".
5. Actualizar la fecha al pie.

---

*Última actualización: 2026-05-18 — cierre de Etapa 5 (discontinuidades interiores embebidas, ADR 0010); fases 4 y 5 diferidas con justificación; próxima decisión = Etapa 6 entre opciones A-E.*
