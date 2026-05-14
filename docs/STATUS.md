# Estado de Fenix FEM

> Foto de una página del estado actual. **Se regenera o actualiza tras cada PR significativa** (cierre de etapa, nuevo componente validado, nuevo ADR aceptado, deuda saldada).
>
> Para visión por etapas y bifurcaciones: ver [ROADMAP.md](ROADMAP.md). Para combinaciones validadas: ver `MATRIZ.md` (pendiente).

---

## Métricas

| Indicador | Valor |
|---|---|
| **Tests** | 401 pasan / 5 skipped / 0 fallos (406 colectados) |
| **Elementos** | 15 (10 estructurales 1D + 5 sólidos 2D) |
| **Materiales** | 8 (2 elásticos + 2 plasticidad + 2 daño + 1 cable + 1 friccional) |
| **Solvers** | 6 (3 estáticos + 1 modal + 2 transitorios) |
| **ADRs aceptados** | 10 (0001–0010) |
| **Specs `validated`** | 23 |
| **Etapas cerradas** | 3 completas + 1 parcial (etapa 4 vía ADR 0009 fases 1, 3, 4) |

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

**Materiales cubiertos**
- Elásticos: `Elastic1D`, `Elastic2D` (plane stress + plane strain).
- Cable unilateral: `CableMaterial1D`.
- Plasticidad J2: `Elastoplastic1D`, `VonMises2D` (kernels Numba especializados para plane strain y plane stress).
- Plasticidad friccional: `DruckerPrager2D` (plane strain, dos ramas regular/apex, plasticidad no asociada).
- Daño isótropo con softening exponencial: `IsotropicDamage1D`, `IsotropicDamage2D` (tangente algorítmica consistente, asimétrica en 2D).

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

---

## Deuda técnica priorizada

| # | Item | Criterio de retoma |
|---|---|---|
| 1 | `internal_forces` devuelve `None` en sólidos 2D (ADR 0002 incompleto). | Cuando entren sólidos 3D, post-proceso avanzado o consumidor externo que pida `ElementForces` para sólidos. |
| 2 | `FiberSection` para plasticidad por flexión en frames. | Si la Etapa 5 se decide por la opción E (Mohr-Coulomb + FiberSection). |
| 3 | Reglas de disparo C y D arquitecturales pendientes ([memoria](../../../.claude/projects/g--Mi-unidad-Proyectos-IA-fenix-fem/memory/project_reglas_disparo_pendientes.md)). | Cuando ocurra el evento que cada regla espera. |

Ninguno de los tres bloquea el avance. Todos están documentados con su contexto en la memoria del proyecto o en [MATRIZ.md](MATRIZ.md), y no requieren acción proactiva.

**Deuda saldada 2026-05-14**: cobertura asimétrica de sólidos 2D no-Quad4 con materiales no lineales — cerrada con [`test_solid_2d_nonlinear_higher_order.py`](../tests/test_solid_2d_nonlinear_higher_order.py) (16 tests verdes a la primera; no había bug latente).

---

## Próximo hito

**Fase 1 del ADR 0010** — discontinuidades interiores embebidas. Subsistema nuevo (fractura computacional vía embedded discontinuity en aproximación discreta, fiel a Retama 2010). Arranca con la spec + implementación de `CohesiveDamageIsotropic` (material traction-jump, Modo-I, daño isótropo, softening lineal/exponencial), seguida del elemento `CST_Embedded2D` (KOS), la validación numérica del `l_d` correcto (Cap. 6 de la tesis) y el benchmark de Van Vliet.

Esta decisión sustituye al "decidir Etapa 5" anterior; las 5 opciones del ROADMAP original quedan diferidas en favor de esta línea por autorización explícita del usuario (autor de la formulación).

---

## Cómo se regenera este documento

1. Actualizar la tabla de métricas contando ADRs, specs, componentes registrados.
2. Si entró un componente nuevo: añadir a "Capacidad de análisis hoy" en la sección correspondiente.
3. Si saldó deuda: tachar el item de la tabla con `~~strikethrough~~` y mover la línea bajo el cuerpo o eliminarla en la siguiente actualización.
4. Si cambió el próximo hito: reescribir la sección "Próximo hito".
5. Actualizar la fecha al pie.

---

*Última actualización: 2026-05-13 — tras aceptación de ADR 0010 (discontinuidades interiores embebidas).*
