# ADR 0001 — Mecanismos anti-caja-negra

**Fecha**: 2026-04-19
**Estado**: Aceptado

## Contexto

Reglas.md §2 establece el modelo "calculadora con manual disponible": el usuario delega la escritura del código pero no su comprensión. Para que ese principio se opere — y no se quede en declaración de intenciones — el proyecto necesita mecanismos prácticos que prevengan la deriva hacia opacidad a medida que el código crece y los conceptos se acumulan.

Sin estos mecanismos, el riesgo realista es: tras 6 meses y 30 commits, el usuario no puede explicar la mitad de los conceptos sostenedores del sistema. La calculadora se ha convertido en caja negra por inercia.

## Decisión

Se introducen tres mecanismos vivos, baratos de mantener y mutuamente reforzantes:

1. **Mapa arquitectural** — `docs/arquitectura.md`. Una página con diagrama de bloques + flujo de datos. Función: índice navegable del sistema. Mantenimiento: la IA lo actualiza cuando cambia la topología (no la implementación interna).

2. **Conceptos clave** — `docs/conceptos_clave.md`. Inventario explícito (10–20 entradas) de los conceptos sostenedores del sistema, cada uno con 2–3 frases conceptuales. Función: define el "test de no-caja-negra" — el usuario debe poder explicar cualquier entrada en términos generales en cualquier momento. Mantenimiento: la IA añade entradas al introducir piezas conceptuales nuevas; el usuario solicita refresco si alguna se le vuelve opaca.

3. **Modo "explícame" bajo demanda** — protocolo conversacional. El usuario invoca "explícame X" sobre archivo, función, concepto o decisión, y la IA responde al nivel solicitado (conceptual, intermedio, detalle). Función: válvula de escape contra opacidad puntual sin esperar a sesiones programadas. No requiere artefacto en el repo.

## Consecuencias

**Positivas**:
- El sistema sigue siendo entendible para el usuario aunque no escriba el código.
- El usuario futuro (6 meses, 2 años) y eventuales revisores externos tienen mapa y glosario para navegar el proyecto.
- La carga de explicación es asíncrona: el usuario no se ve obligado a leer todo en el momento, pero tiene material consultable cuando lo necesite.

**Negativas / costes**:
- Mantener `arquitectura.md` y `conceptos_clave.md` añade overhead a cada cambio estructural. Mitigado: la IA los actualiza junto al código, no como tarea separada.
- Riesgo de drift entre documentos y código si la IA olvida actualizar. Mitigado: revisar como parte del checklist mental al cerrar cualquier refactor mediano.

**Implicación para Reglas.md §4**: los cambios arquitecturales grandes ya generan ADR; ahora también deben revisar `arquitectura.md` y `conceptos_clave.md` por si necesitan actualización.
