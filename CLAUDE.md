# Solidum FEM — instrucciones para la IA

Lee siempre el contrato del proyecto antes de actuar:

@Reglas.md

**Memoria del agente del proyecto** vive en `.claude/memory/` dentro del repo (sincronizada por Drive, excluida de git por `.gitignore`). Al iniciar sesión léela siguiendo el índice `@.claude/memory/MEMORY.md`. No crees memoria del agente en `~/.claude/projects/.../memory/` (el usuario trabaja en varias PCs vía Drive; la memoria fuera del repo no se sincroniza).

Resumen operativo:
- **Preguntas al usuario** (pedido del usuario, 2026-11-12): cada pregunta debe entenderse sola, aunque repita lo dicho antes. En 3-4 líneas de lenguaje llano: qué se decide, por qué importa y la recomendación con su razón. Nunca depender de un mensaje largo anterior («¿apruebas el paso 1?» sin decir qué es el paso 1). Una sola pregunta por mensaje.
- **Verificación y validación** (Reglas §6, 2026-11-12): antes de programar, revisión de la spec, y antes de hacer commit de un cambio de física, verificación, ambas con `teorico-solidum` (también para consultas teóricas); después de implementar, `validador-solidum` (escalera de validación); literatura y criterio de admisión, `bibliotecario-solidum`. Son agentes exclusivos del proyecto (`.claude/agents/`, fuera del repositorio); no usar los generales (`teorico_fem`, `bibliotecario`) en Solidum. Sus informes van a `informes/` y se resumen al usuario.
- **Modelo de colaboración**: la IA decide y escribe la arquitectura y el código (delegación total desde 2026-09-29); el usuario decide la física y valida los resultados. No presentarle decisiones de software. Ver Reglas.md §2.
- **División de dominios**: ver Reglas.md §3 para qué se discute con detalle (modelos físicos/matemáticos) y qué se reporta brevemente (plumbing).
- **Cambios arquitecturales grandes**: registrar un ADR en `docs/adr/000N-titulo.md`. Ver Reglas.md §4.
- **Componentes nuevos**: usar la skill `/solidum-new <kind> <Name>` para scaffolding (kind ∈ material/element/solver). Ver `.claude/skills/solidum-new/SKILL.md`.
- **Componentes no estándar** (formulación de investigación que no es FEM clásico): van a un módulo de usuario `solidum/user/<nombre>/`, fuera del programa principal, con carga explícita (`user_modules` en YAML); scaffolding con `/solidum-new <kind> <Name> --user <nombre>`. Ver Reglas.md §4 y ADR 0020.
