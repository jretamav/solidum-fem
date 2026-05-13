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

- **Componente nuevo (elemento / material / solver)** → el usuario abre una **spec** en `docs/specs/<Nombre>.md` a partir de la plantilla `docs/specs/_template_<kind>.md`. La spec contiene: especificación física, formulación numérica, contrato (YAML con `interface`, `parameters`, `acceptance`, `references`). La IA **no escribe código hasta que la spec está completa**; si falta algo, pregunta en la sección *Diálogo* de la spec antes de implementar. Tras validar (tests verdes contra los casos de `acceptance`), la IA actualiza `status: validated` en la spec y añade una entrada al catálogo correspondiente (`docs/catalogo_<elementos|materiales|solvers>.md`). Relación: **spec = orden de trabajo y referencia detallada por componente**; **catálogo = índice navegable** del conjunto.
- **Variante de componente existente** (subclase de un solver / elemento / material que reusa decisiones de un ADR/spec ya aceptados — p. ej. Newton-Newmark sobre `NewmarkSolver`, HHT-α sobre `NewmarkSolver`, una calibración alternativa de Drucker-Prager) → **spec corta tipo extensión** en `docs/specs/<Variante>.md` que documenta **solo lo que cambia** respecto al componente padre (cambios de formulación, parámetros nuevos, comportamientos distintos) y enlaza el padre. **Sin ADR nuevo** mientras reuse decisiones del ADR padre. **Sin entrada estructural nueva** en el manual si la variante es invisible para el usuario del API (típicamente lo es: la subclase comparte el mismo entrypoint). Sí actualizar el catálogo (entrada propia, o nota en la del padre) cuando la variante sea seleccionable en YAML. Sí test analítico/numérico cuando la variante cambia física, no solo plumbing. El objetivo de esta regla es evitar la sobreingeniería del ritual completo cuando la decisión arquitectural ya está tomada.
- **Cambio físico/matemático** sobre componente existente → ecuaciones, diff acotado, justificación física, propuesta de test. La IA pide visto bueno antes de commitear si el cambio afecta a la formulación, no solo a la implementación. Si la formulación de la spec cambia, se actualiza la spec en el mismo commit.
- **Cambio pequeño de plumbing** → una frase descriptiva, sin requerir lectura del código.
- **Cambio arquitectural grande** (refactor transversal, nuevo subsistema, ruptura de contratos) → la IA crea un **ADR breve** (Architecture Decision Record, ~1 página) en `docs/adr/000N-titulo.md` antes o junto al commit. Los ADRs son persistentes y consultables; sustituyen a las explicaciones efímeras de chat para que el sistema siga siendo entendible meses o años después.
- **Bajo demanda**: si el usuario pide profundizar en cualquier pieza, la IA explica al nivel solicitado, desde lo conceptual hasta el detalle de implementación.

## 5. Convenciones de signos

Convención única del proyecto. Obligatoria en toda formulación, test, API pública (`internal_forces`, `SolveResult`), documentación y catálogo. Si una referencia bibliográfica usa otra convención, se traduce al implementar y se anota la traducción en el test.

**Ejes y giros — 2D estática**
- `x` positivo a la derecha, `y` positivo hacia arriba.
- Giro/momento positivo en sentido antihorario (regla de la mano derecha con `z` saliendo del plano).

**Ejes y giros — 3D estática**
- Regla de la mano derecha para todos los ejes (locales y globales) y todos los momentos (`Mx`, `My`, `Mz`, `T`).

**Esfuerzos internos sobre elemento diferencial 2D — signos por deformación**
Convención clásica de vigas: `N` positivo ↔ tracción; `V` positivo tiende a rotar el diferencial en sentido horario; `M` positivo ↔ tracción en la fibra inferior (sagging).

- Cara izquierda (normal saliente en `−x`): `N` apunta en `−x`, `V` apunta en `+y`, `M` actúa en sentido horario.
- Cara derecha (normal saliente en `+x`): `N` apunta en `+x`, `V` apunta en `−y`, `M` actúa en sentido antihorario.

**Esfuerzos internos en 3D — convención stress-resultant / RHR pura**
En la cara con normal saliente `+x_local`, los esfuerzos positivos son:
- `N` en `+x_local` (tracción positiva).
- `Vy` en `+y_local`, `Vz` en `+z_local`.
- `T ≡ Mx`, `My`, `Mz`: vectores momento en `+x_local`, `+y_local`, `+z_local` respectivamente (regla de la mano derecha).

En la cara con normal saliente `−x_local`, todos los sentidos se invierten (Newton 3ª ley).

Justificación: es la convención que sale de integrar directamente el tensor de tensiones sobre la sección (`N=∫σxx dA`, `Vy=∫σxy dA`, `Vz=∫σxz dA`, `T=∫(y·σxz − z·σxy) dA`, `My=−∫z·σxx dA`, `Mz=∫y·σxx dA`). Es la que adoptan Bathe, Crisfield, Cook-Malkus-Plesha, SAP2000, OpenSees, ANSYS (Beam188/189), Abaqus (B31/B32). La convención estructural de "sagging positivo" no tiene extensión canónica a 3D — se usa solo en 2D.

Convención de fibra para flectores (consecuencia del signo de arriba):
- `Mz > 0` ⇒ tracción en fibras con `y < 0` (equivalente a sagging en el plano xy con `+y=arriba`).
- `My > 0` ⇒ tracción en fibras con `z > 0`.

**Capa interna vs API pública**
- *Interna (formulación)*: todos los elementos — 2D y 3D — trabajan internamente en la convención stress-resultant/RHR. Matrices `B`, fuerzas internas `f_int = ∫Bᵀσ dΩ`, jacobianos, residuales, etc., en esta convención. Evita signos especiales por dimensión.
- *API pública (`internal_forces`, diagramas, catálogo)*: los elementos 3D la exponen tal cual. Los elementos 2D exponen la convención de viga estructural de la sección anterior — `V` con signo opuesto al interno — por ser la intuición clásica de diagramas. La traducción B↔A es un simple cambio de signo en `V`, aplicado dentro de `internal_forces()` del elemento 2D.

Relación 2D↔3D en el API público:
- `M_2D` ≡ `Mz_3D` (mismo signo, ambos sagging positivo con `+y=arriba`).
- `V_2D` y `Vy_3D` difieren en signo en el API público por la razón anterior.

Los valores se exponen siempre en **ejes locales del elemento**; la transformación a globales es una capa separada.

## 6. Salvaguardas

- **Tests blindan física no obvia.** Toda formulación nueva entra acompañada de un test de validación contra solución analítica o benchmark conocido. La IA puede equivocarse en sutilezas físicas (signos en Voigt, factores de ½, hipótesis plane stress vs plane strain) que no rompen la compilación; los tests son la red de seguridad.
- **El usuario lee diffs y commits.** No escribe, sí supervisa. Esta es la salvaguarda primaria contra errores numéricos sutiles.
- **ADRs preservan memoria arquitectural** para el yo futuro del usuario y para revisores externos.
- **Idioma**: docstrings y mensajes al usuario en español; identificadores de código en inglés.
