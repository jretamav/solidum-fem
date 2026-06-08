# Solidum FEM — instrucciones para la IA

Lee siempre el contrato del proyecto antes de actuar:

@Reglas.md

**Memoria del agente del proyecto** vive en `.claude/memory/` dentro del repo (sincronizada por Drive, excluida de git por `.gitignore`). Al iniciar sesión léela siguiendo el índice `@.claude/memory/MEMORY.md`. No crees memoria del agente en `~/.claude/projects/.../memory/` (el usuario trabaja en varias PCs vía Drive; la memoria fuera del repo no se sincroniza).

Resumen operativo:
- **Modelo de colaboración**: la IA escribe el código; el usuario revisa diffs y resultados físicos. Ver Reglas.md §2.
- **División de dominios**: ver Reglas.md §3 para qué se discute con detalle (modelos físicos/matemáticos) y qué se reporta brevemente (plumbing).
- **Cambios arquitecturales grandes**: registrar un ADR en `docs/adr/000N-titulo.md`. Ver Reglas.md §4.
- **Componentes nuevos**: usar la skill `/solidum-new <kind> <Name>` para scaffolding (kind ∈ material/element/solver). Ver `.claude/skills/solidum-new/SKILL.md`.
