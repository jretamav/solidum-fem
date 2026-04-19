# Reglas de Fenix FEM — contrato de colaboración Usuario ↔ IA

Este documento establece los principios del proyecto y el protocolo de trabajo entre el usuario humano y la IA. La IA debe leerlo al inicio de cada sesión.

---

## 0. Identidad del proyecto

**Fenix FEM** es un programa de elementos finitos en aproximación de desplazamientos para investigación en mecánica de sólidos: simulación de problemas mecánico y térmico, acoplados o desacoplados. Su desarrollo debe ser riguroso tanto en arquitectura de software como en métodos numéricos. La arquitectura debe permitir el crecimiento futuro **sin grandes cambios en el código fuente original**.

## 1. Arquitectura optimizada para extensión por IA

El programa debe estar optimizado para que la IA implemente con bajo coste futuras incorporaciones de **elementos**, **materiales** y **métodos de solución** derivados de la investigación del usuario en mecánica del medio continuo. Cada decisión arquitectural se evalúa también contra el coste de añadir el componente N+1: si un cambio no abarata las extensiones futuras, no se justifica solo por estética. Equilibrio constante entre velocidad de ejecución y consumo de RAM.

Concretamente, esto se traduce en preferir:
- Contratos declarativos (`STRAIN_DIM`, `DOF_NAMES`, `N_INTEGRATION_POINTS`) sobre métodos abstractos cuando el comportamiento puede inferirse del dato.
- Auto-registro vía decoradores + descubrimiento automático sobre listas manuales de imports.
- Validación temprana al construir con mensajes claros sobre fallos crípticos en runtime.
- Resistir abstracciones especulativas no justificadas por al menos dos casos reales del dominio.

## 2. Modelo de colaboración: la IA como calculadora con manual disponible

Toda modificación o ampliación del código de Fenix FEM se realiza a través de la IA. El flujo es **Usuario → IA → Código → resultados → Usuario**: la IA produce y modifica el código; el usuario cierra el lazo de validación sobre los resultados físicos y las decisiones arquitecturales.

El usuario delega la **escritura** del código (sintaxis, plumbing, patrones de software) pero **no la comprensión conceptual**. La analogía operativa es la de una calculadora científica: el usuario sabe qué es el seno, qué propiedades tiene y qué orden de magnitud esperar, pero no memoriza el algoritmo CORDIC interno; sí detecta resultados absurdos. Trasladado al proyecto, el usuario debe tener un modelo mental completo del **qué** y del **por qué** de cada pieza arquitectural — sin necesidad de escribir su implementación.

La IA, por tanto, **nunca trata el sistema como caja negra para el usuario**. Cuando introduce una pieza arquitectural nueva (decorador, autodiscover, contrato declarativo, cache, etc.), explica en lenguaje conceptual y en pocas frases: qué es, por qué resuelve el problema, dónde vive en el repositorio.

## 3. División de dominios

| Dominio | Quién decide | Qué se comunica |
|---|---|---|
| **Modelos físicos y matemáticos** — materiales (constitutivos, return mapping, criterios de fluencia/daño, módulos tangentes), elementos (cinemática, matriz B, formulación variacional, integración, transformaciones), solvers (criterios de convergencia, predictores/correctores, paso adaptativo, longitud de arco), constantes con significado físico (tolerancias de fluencia, factor de penalidad). | **Usuario revisa con detalle.** | Ecuaciones (LaTeX/ASCII), referencia bibliográfica si aplica, diff acotado, propuesta de test contra solución analítica o benchmark conocido. |
| **Plumbing de software** — registry, autodiscover, decoradores, contratos abstractos, ensamblaje sparse, cache de topología COO, vectorización, Numba, parsers (YAML, gmsh, VTK), scaffolding, estructura de directorios, type hints. | **IA decide.** | Una frase de qué cambió y por qué; tests verdes como evidencia. No se pide al usuario que lea la implementación. |
| **Zona gris** — clases base (`Element`, `Material`), semántica trial/commit de variables internas, esquema de `ElementState`. | IA propone, usuario valida concepto. | Aviso conceptual breve; lectura detallada solo si el usuario la pide. |

## 4. Protocolo de presentación de cambios

- **Cambio físico/matemático** → ecuaciones, diff acotado, justificación física, propuesta de test. La IA pide visto bueno antes de commitear si el cambio afecta a la formulación, no solo a la implementación.
- **Cambio pequeño de plumbing** → una frase descriptiva, sin requerir lectura del código.
- **Cambio arquitectural grande** (refactor transversal, nuevo subsistema, ruptura de contratos) → la IA crea un **ADR breve** (Architecture Decision Record, ~1 página) en `docs/adr/000N-titulo.md` antes o junto al commit. Los ADRs son persistentes y consultables; sustituyen a las explicaciones efímeras de chat para que el sistema siga siendo entendible meses o años después.
- **Bajo demanda**: si el usuario pide profundizar en cualquier pieza, la IA explica al nivel solicitado, desde lo conceptual hasta el detalle de implementación.

## 5. Salvaguardas

- **Tests blindan física no obvia.** Toda formulación nueva entra acompañada de un test de validación contra solución analítica o benchmark conocido. La IA puede equivocarse en sutilezas físicas (signos en Voigt, factores de ½, hipótesis plane stress vs plane strain) que no rompen la compilación; los tests son la red de seguridad.
- **El usuario lee diffs y commits.** No escribe, sí supervisa. Esta es la salvaguarda primaria contra errores numéricos sutiles.
- **ADRs preservan memoria arquitectural** para el yo futuro del usuario y para revisores externos.
- **Idioma**: docstrings y mensajes al usuario en español; identificadores de código en inglés.
