# Onboarding — primer día en Fenix FEM

> Documento de entrada para agente IA o desarrollador humano que arranca sesión sin contexto previo. Lo lees una vez, te orientas, y vas a la fuente concreta para profundizar.
>
> **No duplica contenido**: apunta a `Reglas.md`, `CLAUDE.md`, ADRs, specs, catálogos y memoria. Si entran en conflicto, **`Reglas.md` manda** (es el contrato firmado entre usuario y agente).

---

## 0. Cinco minutos para situarse

Fenix FEM es un programa de elementos finitos en aproximación de desplazamientos para **investigación en mecánica de sólidos** (mecánico ± térmico, acoplado o desacoplado). Hoy resuelve:

- Estática lineal y no lineal (material y geométrica) sobre **1D estructural** y **sólidos 2D**.
- **Subsistema modal/dinámico/espectral completo** (ADR 0009 cerrado 2026-05-18): modal por autovalores generalizados, transitorio implícito Newmark/HHT-α (lineal y no lineal con Newton), transitorio explícito por diferencias centradas, respuesta forzada armónica en frecuencia y análisis sísmico por combinación modal espectral (SRSS/CQC).
- Catálogo de materiales con plasticidad J2 (plane strain + plane stress), Drucker-Prager, daño isótropo y cohesivo traction-jump.
- Fractura computacional vía discontinuidades embebidas (CST_Embedded2D con condensación local).

Lo que **no** resuelve aún: sólidos 3D, placas/láminas, térmico, contacto, Mohr-Coulomb, FiberSection. Ver [STATUS.md](STATUS.md) §"Limitaciones declaradas".

---

## 1. Orden de lectura recomendado

Para arrancar sin contexto previo, lee en este orden — **no necesitas más para empezar**:

1. **[Reglas.md](../Reglas.md)** (~5 min): contrato Usuario ↔ IA. Identidad del proyecto, división de dominios, convenciones de signos, salvaguardas. **Volver siempre que dudes sobre qué decisión te toca**.
2. **[STATUS.md](STATUS.md)** (~2 min): foto del estado actual. Métricas, capacidades, deuda técnica, próximo hito.
3. **[ROADMAP.md](ROADMAP.md)** (~5 min): etapas cerradas y bifurcación pendiente. Hoy: Etapas 1-6 cerradas; Etapa 7 abierta entre las opciones A, B, C, E del catálogo de bifurcaciones.
4. **[MATRIZ.md](MATRIZ.md)** (~3 min): qué combinaciones elemento × material son válidas y testeadas.
5. **Último ADR aceptado** ([ADR 0011 — Robustez de Newton-Raphson con line search](adr/0011-robustez-newton-line-search.md)): para entender la última decisión arquitectural grande. Línea de fractura computacional articulada por [ADR 0010](adr/0010-discontinuidades-interiores-embebidas.md); subsistema dinámico/modal/espectral articulado por [ADR 0009](adr/0009-analisis-modal-y-dinamico.md) (cerrado en su totalidad 2026-05-18, fases 1-7 + variante HHT-α + reglas C y D aplicadas).

Si la sesión va sobre un componente concreto: ir directo a su spec en `docs/specs/<Nombre>.md` y su entrada en `docs/catalogo_<elementos|materiales|solvers>.md`.

---

## 2. Estructura del repositorio

```
fenix_fem/
├── Reglas.md, CLAUDE.md          ← Contrato y guía operativa
├── fenix/                        ← Código fuente
│   ├── core/                     ← Domain, Node, Element base, Material base, Assembler
│   ├── elements/                 ← truss, cable, frame/, frame3d, solid_2d/
│   ├── materials/                ← elastic, plastic_1d, von_mises_2d, drucker_prager_2d, damage_*
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
├── tests/                        ← 703 verdes + 5 skipped (pytest)
├── docs/
│   ├── adr/                      ← 0001-0010: decisiones arquitecturales
│   ├── specs/                    ← una por componente: contrato + acceptance
│   ├── referencias/              ← PDFs/papers citados desde ADRs y specs (tesis Retama 2010, etc.)
│   ├── catalogo_*.md             ← índices navegables por dominio
│   ├── ROADMAP.md, STATUS.md,
│   ├── MATRIZ.md, ONBOARDING.md  ← documentos navegacionales (este set)
├── manuals/
│   ├── sources/{reference,user,architecture}/  ← markdown fuente
│   └── build_*_manual.py         ← un builder por manual; salida en LaTeX/PDF
└── examples/                     ← YAMLs de ejemplo (estático, modal, transitorio)
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
python -c "import fenix; print(fenix.run_yaml('examples/<archivo>.yaml'))"
```

---

## 4. Memoria del proyecto

Hay dos memorias que persisten entre sesiones:

- **`C:\Users\jreta\.claude\projects\g--Mi-unidad-Proyectos-IA-fenix-fem\memory\MEMORY.md`** y los `*.md` que indexa: feedback del usuario, estado de proyecto, referencias externas. **Cárgalas al arranque cuando sea relevante**. Tipos: `user`, `feedback`, `project`, `reference`. Detalle en el system prompt; protocolo en CLAUDE.md global.
- **`docs/specs/<Nombre>.md`**: memoria *física* y *numérica* de cada componente. Donde se guarda la formulación rigurosa.

**Antes de actuar**: si la memoria menciona una decisión, **verificar que sigue vigente** consultando el código (los memos envejecen).

---

## 5. División de dominios — qué decide el agente vs el usuario

Resumen ágil de Reglas.md §3 (no la sustituye; léela entera al menos una vez).

- **Decide el agente, reporta una frase**: plumbing puro — registry, decoradores, scaffolding, type hints, refactor zonal, vectorización, caché, manejo sparse.
- **Decide el usuario, el agente propone con ecuaciones + test analítico**: formulación física-matemática — matriz `B`, return mapping, criterios de convergencia, constantes con significado físico, traducción de convenciones bibliográficas.
- **Zona gris (el agente propone, el usuario valida concepto)**: clases base, semántica trial/commit, API pública (`SolveResult`, `fenix.run*`), decisiones con ADR en curso. **Ante la duda, eleva a zona gris**.

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
- **Idioma**: docstrings y mensajes al usuario en **español**; identificadores en **inglés** (Reglas.md §6).
- **Unidades**: el sistema no convierte. El usuario es responsable de la consistencia (Reglas.md §5 + ADR 0008).
- **Tolerancias**: patrón `f ≤ atol + rtol · escala` con escala adimensional autoderivada (ADR 0006 / 0007). No introducir constantes hardcoded.

---

## 8. Antes de commitear

1. `python -m pytest tests/ -q` verde (703 pasan, 5 skipped es normal).
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

*Última actualización: 2026-05-14 — primer redactado tras formalizar el set de documentos navegacionales (ROADMAP + STATUS + MATRIZ + ONBOARDING).*
