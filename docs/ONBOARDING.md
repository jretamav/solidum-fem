# Onboarding — primer día en Solidum FEM

> Documento de entrada para agente IA o desarrollador humano que arranca sesión sin contexto previo. Lo lees una vez, te orientas, y vas a la fuente concreta para profundizar.
>
> **No duplica contenido**: apunta a `Reglas.md`, `CLAUDE.md`, ADRs, specs, catálogos y memoria. Si entran en conflicto, **`Reglas.md` manda** (es el contrato firmado entre usuario y agente).

---

## 0. Cinco minutos para situarse

Solidum FEM es un programa de elementos finitos en aproximación de desplazamientos para **investigación en mecánica de sólidos** (mecánico ± térmico, acoplado o desacoplado). Hoy resuelve:

- Estática lineal y no lineal (material y geométrica) sobre **1D estructural**, **sólidos 2D** y **sólidos 3D lineales** (Hex8, Tet4 — Etapa 7 cerrada 2026-05-19, ADR 0012).
- **Subsistema modal/dinámico/espectral completo** (ADR 0009 cerrado 2026-05-18): modal por autovalores generalizados, transitorio implícito Newmark/HHT-α (lineal y no lineal con Newton), transitorio explícito por diferencias centradas, respuesta forzada armónica en frecuencia y análisis sísmico por combinación modal espectral (SRSS/CQC).
- Catálogo de materiales con plasticidad J2 (plane strain + plane stress), Drucker-Prager, daño isótropo y cohesivo traction-jump. Elásticos en 1D, 2D y 3D.
- Fractura computacional vía discontinuidades embebidas (CST_Embedded2D con condensación local).

Lo que **no** resuelve aún: sólidos 3D cuadráticos (Hex20/Hex27/Tet10), materiales 3D no lineales (VonMises3D/DruckerPrager3D/IsotropicDamage3D), placas/láminas, térmico, contacto, Mohr-Coulomb, FiberSection. Ver [STATUS.md](STATUS.md) §"Limitaciones declaradas".

### Estado de publicación (2026-05-20)

El repo es **público en GitHub** bajo `https://github.com/jretamav/solidum-fem` con licencia **LGPL-3.0** (ver `LICENSE` y `LICENSE.GPL`). Decisiones macro de publicación están en memoria:

- **Paper 1a → JOSS**. Draft completo en `paper/joss/paper.md` (~897 palabras, dentro del techo JOSS 1000) + `paper/joss/paper.bib` con 16 entradas auditadas contra CrossRef/JOSS/sitios oficiales. Sección no-estándar "Statement on the use of generative artificial intelligence" añadida en línea con la práctica del paper CAEE-Winkler del autor. Commits `155558e` (redacción) y `bba6d99` (enriquecimiento + auditoría bib) acumulados, **2 ahead de `origin/main`, push diferido para relectura en frío**. Checklist operativo de submission (push → tag → Zenodo DOI → formulario JOSS) en `[[project_paper_joss_estado]]`.
- **Paper 1b → SoftwareX** (software paper integral). Manuscrito vivirá fuera de Solidum (Overleaf o repo separado). Pendiente sesión dedicada tras cerrar el JOSS.
- **Papers 2+ → metodológicos** en IJNME / Computational Mechanics / Engineering Fracture Mechanics / AES / Computer Physics Communications (embedded KOS, dissipation arc-length, etc.).
- AES fue evaluada y **descartada como destino del paper 1** tras consulta del scope real (la revista viró hacia PINN/ML; ya no publica framework papers).

Estrategia detallada en `[[project_estrategia_publicaciones]]`; estado vivo del paper JOSS en `[[project_paper_joss_estado]]`; estado del repo público en `[[project_repo_publico_decidido]]`; columna vertebral conceptual del proyecto para futuras comunicaciones en `[[project_solidum_motivacion_fundacional]]`.

---

## 1. Orden de lectura recomendado

Para arrancar sin contexto previo, lee en este orden — **no necesitas más para empezar**:

1. **[Reglas.md](../Reglas.md)** (~5 min): contrato Usuario ↔ IA. Identidad del proyecto, división de dominios, convenciones de signos (Voigt 2D y 3D), salvaguardas. **Volver siempre que dudes sobre qué decisión te toca**.
2. **[STATUS.md](STATUS.md)** (~2 min): foto del estado actual. Métricas, capacidades, deuda técnica, próximo hito.
3. **[ROADMAP.md](ROADMAP.md)** (~5 min): etapas cerradas y bifurcación pendiente. Hoy: Etapas 1-7 cerradas; Etapa 8 abierta entre opciones B (placas/láminas), C (térmico), E (Mohr-Coulomb + FiberSection) y las ramificaciones naturales de la Etapa 7 (A.bis materiales 3D no lineales, A.ter cuadráticos 3D).
4. **[MATRIZ.md](MATRIZ.md)** (~3 min): qué combinaciones elemento × material son válidas y testeadas.
5. **Último ADR aceptado** ([ADR 0012 — Sólidos 3D y Voigt 6D](adr/0012-solidos-3d-y-voigt-6d.md)): para entender la última decisión arquitectural grande. Convención Voigt 3D del proyecto, cierre del contrato `internal_forces` por dominio explícito (sólidos exponen `compute_gauss_state`), API de caras 3D con normal saliente. Anteriores: [ADR 0011](adr/0011-robustez-newton-line-search.md) (robustez Newton), [ADR 0010](adr/0010-discontinuidades-interiores-embebidas.md) (embedded discontinuities), [ADR 0009](adr/0009-analisis-modal-y-dinamico.md) (subsistema dinámico).

Si la sesión va sobre un componente concreto: ir directo a su spec en `docs/specs/<Nombre>.md` y su entrada en `docs/catalogo_<elementos|materiales|solvers>.md`.

---

## 2. Estructura del repositorio

```
solidum_fem/
├── Reglas.md, CLAUDE.md          ← Contrato y guía operativa
├── solidum/                        ← Código fuente
│   ├── core/                     ← Domain, Node, Element base, Material base, Assembler
│   ├── elements/                 ← truss, cable, frame/, frame3d, solid_2d/, solid_3d/
│   ├── materials/                ← elastic, elastic_2d, elastic_3d, plastic_1d, von_mises_2d, drucker_prager_2d, damage_*
│   ├── math/
│   │   ├── solvers/              ← linear, nonlinear, arclength, modal, newmark (+HHT),
│   │   │                            central_difference, harmonic, response_spectrum
│   │   ├── linalg/               ← dispatcher, eigen
│   │   ├── mass_lumping.py       ← HRZ canónico (ADR 0009 fase 2)
│   │   ├── modal_response.py     ← free_vibration + SRSS/CQC + helpers de espectros
│   │   ├── integration.py        ← QuadratureRegistry
│   │   └── geometry.py, damping.py, convergence.py
│   ├── results.py                ← SolveResult, ModalResult, TransientResult,
│   │                                HarmonicResult, ResponseSpectrumResult
│   ├── registry.py, constants.py, logging.py
│   └── utils/                    ← YAML parser, gmsh parser, VTK exporter
├── tests/                        ← 804 verdes + 5 skipped (pytest); 38 en tests/validation/ contra benchmarks publicados (Lamé 2D, NAFEMS LE1, MacNeal-Harder, Bathe wave, Hill J2) + 5 nuevos 3D en tests/validation/ (cubo Lamé 3D, MacNeal 3D)
├── docs/
│   ├── adr/                      ← 0001-0012: decisiones arquitecturales
│   ├── specs/                    ← una por componente: contrato + acceptance
│   ├── referencias/              ← PDFs/papers citados desde ADRs y specs (tesis Retama 2010, etc.)
│   ├── catalogo_*.md             ← índices navegables por dominio
│   ├── ROADMAP.md, STATUS.md,
│   ├── MATRIZ.md, ONBOARDING.md  ← documentos navegacionales (este set)
├── manuals/
│   ├── sources/{reference,user,architecture}/  ← markdown fuente
│   ├── build_*_manual.py         ← un builder por manual; salida en LaTeX/PDF
│   └── glossary.md               ← (pendiente) glosario ES↔EN para regenerar manuales en inglés
├── paper/joss/                   ← paper.md + paper.bib para JOSS (skeleton)
├── examples/                     ← YAMLs de ejemplo (estático, modal, transitorio)
├── README.md, CONTRIBUTING.md, CITATION.cff, LICENSE, LICENSE.GPL
└── pyproject.toml
```

---

## 3. Comandos esenciales

```bash
# Tests
python -m pytest tests/ -q                # suite completa (~30 s en Windows)
python -m pytest tests/test_modal.py -q   # un solo módulo

# Manuales (regenerar tras cambios en sources/)
python manuals/build_reference_manual.py
python manuals/build_user_manual.py
python manuals/build_architecture_manual.py

# Pipeline YAML end-to-end
python -c "import solidum; print(solidum.run_yaml('examples/<archivo>.yaml'))"
```

---

## 4. Memoria del proyecto

Hay dos memorias que persisten entre sesiones:

- **`C:\Users\jreta\.claude\projects\g--Mi-unidad-Proyectos-IA-solidum-fem\memory\MEMORY.md`** y los `*.md` que indexa: feedback del usuario, estado de proyecto, referencias externas. **Cárgalas al arranque cuando sea relevante**. Tipos: `user`, `feedback`, `project`, `reference`. Detalle en el system prompt; protocolo en CLAUDE.md global.
- **`docs/specs/<Nombre>.md`**: memoria *física* y *numérica* de cada componente. Donde se guarda la formulación rigurosa.

**Antes de actuar**: si la memoria menciona una decisión, **verificar que sigue vigente** consultando el código (los memos envejecen).

---

## 5. División de dominios — qué decide el agente vs el usuario

Resumen ágil de Reglas.md §3 (no la sustituye; léela entera al menos una vez).

- **Decide el agente, reporta una frase**: plumbing puro — registry, decoradores, scaffolding, type hints, refactor zonal, vectorización, caché, manejo sparse.
- **Decide el usuario, el agente propone con ecuaciones + test analítico**: formulación física-matemática — matriz `B`, return mapping, criterios de convergencia, constantes con significado físico, traducción de convenciones bibliográficas.
- **Zona gris (el agente propone, el usuario valida concepto)**: clases base, semántica trial/commit, API pública (`SolveResult`, `solidum.run*`), decisiones con ADR en curso. **Ante la duda, eleva a zona gris**.

---

## 6. Protocolo de cambios — qué tipo de cambio merece qué artefacto

Resumen ágil de Reglas.md §4 (la fuente es ese párrafo).

| Tipo de cambio                                          | Artefacto                                                       |
|---------------------------------------------------------|-----------------------------------------------------------------|
| Plumbing pequeño                                        | Una frase en el commit.                                         |
| Refactor zonal (centralización, renombre estructural)   | Diff + rationale en commit + actualización de specs/manuales.   |
| Cambio físico sobre componente existente                | Ecuaciones + diff + test analítico. Visto bueno antes de commit si toca formulación. |
| Componente nuevo (elemento / material / solver)         | Spec primero (`docs/specs/<Nombre>.md`) → implementación → catálogo. |
| Variante de componente existente                        | Spec corta (solo lo que cambia) + nota en catálogo si es seleccionable por YAML. Sin ADR. |
| Cambio arquitectural grande                             | ADR en `docs/adr/000N-titulo.md` antes o junto al commit.       |

---

## 7. Convenciones críticas a no olvidar

- **Signos**: convención única del proyecto (Reglas.md §5). 2D interno = stress-resultant; 2D público (`internal_forces`) = sagging-positivo de vigas. La traducción B↔A vive dentro de `internal_forces()` de cada elemento 2D.
- **Idioma**: política híbrida por capa (Reglas.md §6 — política revisada 2026-05-20 al abrir el repo al público). Externa pública (README, CITATION, CONTRIBUTING, paper.md, manuscritos) en **inglés**; manuales **bilingües ES+EN** (español canónico, inglés generado por la IA con glosario terminológico — glosario pendiente en `manuals/glossary.md`); documental navegacional (STATUS, ROADMAP, ONBOARDING, MATRIZ, catálogos) en **español**; interna técnica (specs, ADRs, docstrings, commits) en **español**; identificadores de código en **inglés**. Ver `[[project_politica_idioma]]` para el detalle.
- **Unidades**: el sistema no convierte. El usuario es responsable de la consistencia (Reglas.md §5 + ADR 0008).
- **Tolerancias**: patrón `f ≤ atol + rtol · escala` con escala adimensional autoderivada (ADR 0006 / 0007). No introducir constantes hardcoded.

---

## 8. Antes de commitear

1. `python -m pytest tests/ -q` verde (741 pasan, 5 skipped es normal).
2. Si el cambio afecta a un componente con spec: actualizar la spec en el mismo commit si la formulación cambió, o subir `status: validated` si el componente se acaba de validar.
3. Si el cambio renombra/elimina símbolos públicos: barrer specs, catálogos y manuales que los mencionen, en el mismo commit.
4. **Pre-commit hooks**: no usar `--no-verify`. Si un hook falla, investigar la causa raíz.
5. **Commit message**: en español, descriptivo, una línea + cuerpo si hace falta. Co-autoría con Claude al pie es opcional.

---

## 9. Si te bloqueas

- **Falta contexto sobre un componente**: ir a su spec, no a chat con el usuario.
- **Falta contexto histórico sobre una decisión**: buscar en `docs/adr/`. Si no está, es señal de que la decisión no se ha tomado formalmente (zona gris).
- **El usuario delega "decide tú"**: **responsabilidad amplificada, no descarga**. Aplicar reglas + memoria + juicio arquitectural. Saltarse cualquiera produce decisiones plausibles pero mal alineadas. Ver memoria `feedback-rol-experto-arquitectura`.
- **Una decisión sale del catálogo de los 6 tipos del §6**: probablemente es zona gris — proponer en chat antes de actuar.

---

## 10. Lo que este documento NO es

- **No es un tutorial de FEM**: presupone formación del lector (el usuario tiene doctorado en mecánica del medio continuo).
- **No es la documentación del API**: el manual de referencia (`manuals/Reference_manual.pdf`) cumple ese rol.
- **No es un changelog**: para eso está `git log` y los ADRs.
- **No reemplaza Reglas.md**: si discrepan, manda Reglas.md.

---

*Última actualización: 2026-05-20 (segunda sesión del día) — **paper JOSS draft completo y verificado** en `paper/joss/paper.md` (~897 palabras, 16 entradas bibliográficas auditadas, Statement of need enriquecido con la inversión de dirección + paisaje de 3 niveles + prueba de existencia embedded discontinuity). Commits `155558e` + `bba6d99` acumulados sin pushear (push diferido para relectura en frío). **Próximo hito estratégico**: flujo de submission JOSS — relectura → push → git tag de release → archivo Zenodo (verificar/configurar linkage) → formulario en https://joss.theoj.org/papers/new. Ver `[[project_paper_joss_estado]]` para el checklist operativo.*
