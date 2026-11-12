# Solidum FEM — instrucciones para la IA

Lee siempre el contrato del proyecto antes de actuar:

@Reglas.md

**Memoria del agente del proyecto** vive en `.claude/memory/` dentro del repo (sincronizada por OneDrive, excluida de git por `.gitignore`). No crees memoria del agente en `~/.claude/projects/.../memory/` (el usuario trabaja en varias PCs; la memoria fuera del repo no se sincroniza). Su índice se carga al iniciar sesión (la importación sólo funciona fuera de comillas de código; mantenerlo corto, una línea por nota):

@.claude/memory/MEMORY.md

Resumen operativo:
- **Preguntas al usuario** (pedido del usuario, 2026-11-12): cada pregunta debe entenderse sola, aunque repita lo dicho antes. En 3-4 líneas de lenguaje llano: qué se decide, por qué importa y la recomendación con su razón. Nunca depender de un mensaje largo anterior («¿apruebas el paso 1?» sin decir qué es el paso 1). Una sola pregunta por mensaje, en texto libre, nunca con la ventana de opciones (`AskUserQuestion`). Si surgen varias subdecisiones, proponer el orden y preguntar sólo la primera; un «sí» no cierra las siguientes.
- **Mensajes al usuario** (2026-11-12): su atención es el recurso más escaso. Estructura fija: qué hice, qué encontré, qué necesito de ti. Sólo las cifras que sirven para decidir; sin jerga de software (commit, hook, plumbing) si no la necesita. Cada afirmación marcada según su origen: *medido*, *calculado*, *leído en la fuente* o *hipótesis*. No preguntarle lo que Reglas §2-§3 delegan a la IA: hacerlo y reportarlo en una línea.
- **Un tema por sesión** (2026-11-12): las sesiones largas se resumen automáticamente y se pierden detalles. Al cambiar de tema, proponer `/cerrar-sesion` y una sesión nueva. Delegar a los agentes del proyecto lo que llena el contexto (leer artículos, revisar registros o diffs largos) y quedarse con sus conclusiones.
- **Orden de archivos** (2026-11-11): nada suelto en la raíz. Reportes, PDF generados e informes de agentes van a `informes/` (fuera de Git); trabajos intermedios, al scratchpad; artículos y tesis, a `docs/referencias/`. Si no hay carpeta adecuada, proponerla.
- **Al correr ejemplos o validaciones** (2026-11-11): la IA usa Solidum como un usuario más (`.geo`/YAML, `python -m solidum`, salida y figuras estándar) y lo menos posible como calculista; lo que Solidum no haga se propone como mejora en vez de rodearlo con un guion. Al implementar, el papel de Reglas §2-§4.
- **Dónde vive cada decisión** (2026-11-12): física → la spec (su Diálogo); software → el ADR; forma de trabajar → `Reglas.md` y este archivo; estado del trabajo → `ESTADO.md`; la memoria sólo para preferencias y referencias que no estén en otro lado. Al cambiar una decisión, actualizar su lugar canónico y buscar copias desactualizadas.
- **Verificación y validación** (Reglas §6, 2026-11-12): antes de programar, revisión de la spec, y antes de hacer commit de un cambio de física, verificación, ambas con `teorico-solidum` (también para consultas teóricas); después de implementar, `validador-solidum` (escalera de validación); literatura y criterio de admisión, `bibliotecario-solidum`; en cambios de arquitectura grandes (ADR, subsistema, refactorización transversal, rendimiento), `arquitecto-solidum` revisa el plan antes de programar y el diff antes del commit (en plumbing menor basta `/code-review`). Son agentes exclusivos del proyecto (`.claude/agents/`, fuera del repositorio); no usar los generales (`teorico_fem`, `bibliotecario`) en Solidum. Sus informes van a `informes/` y se resumen al usuario.
- **Modelo de colaboración**: la IA decide y escribe la arquitectura y el código (delegación total desde 2026-09-29); el usuario decide la física y valida los resultados. No presentarle decisiones de software. Ver Reglas.md §2.
- **División de dominios**: ver Reglas.md §3 para qué se discute con detalle (modelos físicos/matemáticos) y qué se reporta brevemente (plumbing).
- **Cambios arquitecturales grandes**: registrar un ADR en `docs/adr/000N-titulo.md`. Ver Reglas.md §4.
- **Componentes nuevos**: usar la skill `/solidum-new <kind> <Name>` para scaffolding (kind ∈ material/element/solver). Ver `.claude/skills/solidum-new/SKILL.md`.
- **Componentes no estándar** (formulación de investigación que no es FEM clásico): van a un módulo de usuario `solidum/user/<nombre>/`, fuera del programa principal, con carga explícita (`user_modules` en YAML); scaffolding con `/solidum-new <kind> <Name> --user <nombre>`. Ver Reglas.md §4 y ADR 0020.
