# Filosofía del proyecto

## Identidad del programa

Fenix FEM es un programa de Método de Elementos Finitos (MEF) en aproximación de desplazamientos, orientado a la investigación en mecánica de sólidos: simulación de problemas mecánico y térmico, acoplados o desacoplados. Su público objetivo es la comunidad de investigación en mecánica computacional que requiere incorporar materiales, elementos o esquemas de solución sin verse limitada por la rigidez interna del programa.

## Hipótesis fundacional

La mayor parte del coste de un programa de MEF de investigación no reside en la implementación del primer elemento, sino en el mantenimiento de un coste bajo para la incorporación del elemento N+1. Todo proyecto que crece sin disciplina de arquitectura acaba penalizando cada nueva incorporación: el primer material se implementa en una sesión; el décimo requiere modificar varios archivos; el vigésimo se desincentiva. Fenix FEM se diseña explícitamente para que la pendiente de coste de extensión se mantenga plana.

## Modelo de colaboración

El desarrollo de Fenix FEM se rige por un contrato explícito entre el usuario (investigador con dominio profundo en MEF) y un asistente de inteligencia artificial responsable de la escritura del código. El flujo es: usuario → asistente → código → resultados → usuario. El asistente produce y modifica el código; el usuario cierra el lazo de validación sobre los resultados físicos y las decisiones de arquitectura.

La analogía operativa es la de una calculadora científica: se conoce qué representa la función seno, sus propiedades y el orden de magnitud esperado de su resultado, pero no se memoriza el algoritmo CORDIC interno; sí se detectan resultados absurdos. Trasladado al proyecto: la escritura del código (sintaxis, infraestructura, patrones de software) se delega; la comprensión conceptual no. El usuario mantiene un modelo mental completo del *qué* y del *por qué* de cada pieza de la arquitectura sin necesidad de escribir su implementación.

## Principios rectores

**Optimización para extensión.** Cada decisión de arquitectura se evalúa también contra el coste de añadir el componente N+1. Un cambio que no abarata las extensiones futuras no se justifica únicamente por motivos estéticos. La consecuencia práctica es la preferencia por contratos declarativos sobre métodos abstractos cuando el comportamiento puede inferirse del dato, y por mecanismos de auto-registro sobre listas manuales de importaciones.

**Equilibrio entre velocidad de cómputo y consumo de memoria.** Las decisiones de implementación interna (estructuras dispersas, almacenamiento en formato Coordinate (COO) en caché, vectorización, compilación Just-In-Time (JIT)) buscan el compromiso entre tiempo de ejecución y consumo de memoria, y no la optimización aislada de uno solo de los dos ejes.

**Mecanismos de transparencia conceptual.** El sistema mantiene tres artefactos vivos para preservar la comprensión del usuario a medida que el código crece: un mapa de arquitectura (este manual y `docs/arquitectura.md`), un inventario explícito de los conceptos sostenedores (`docs/conceptos_clave.md`) y un protocolo conversacional para profundizar bajo demanda. Estos mecanismos se formalizan en el Architecture Decision Record (ADR) 0001.

**Validación física mediante pruebas.** Las sutilezas físicas (signos en notación de Voigt, factores de un medio, hipótesis de tensión plana frente a deformación plana) no se manifiestan en errores de compilación. Toda formulación nueva se acompaña de una prueba de validación contra solución analítica o frente a un caso de referencia conocido.

**Convención de signos única.** Toda formulación, prueba, interfaz pública y documentación utiliza la misma convención de signos. Si una referencia bibliográfica adopta otra, la traducción se realiza al implementar y se anota en la prueba correspondiente. Esta convención se detalla en el capítulo dedicado.

## Cuándo usar Fenix FEM y cuándo no

La identidad de un programa se define tanto por su alcance como por sus límites. Esta sección delimita explícitamente lo que Fenix FEM pretende ser y lo que no.

**Casos para los que Fenix FEM es la herramienta indicada.** Investigación en mecánica de sólidos donde la prioridad es la **transparencia y la modificabilidad** del código: prototipado de modelos constitutivos nuevos antes de portarlos a un programa industrial; ensayo de variantes algorítmicas en solvers no lineales (predictores, correctores, criterios de longitud de arco); evaluación de formulaciones de elementos finitos contra soluciones analíticas y casos de referencia conocidos; trabajo docente y de tesis donde el estudiante debe comprender, modificar y validar cada pieza del cálculo.

**Casos para los que Fenix FEM no es la herramienta indicada.** Análisis de producción de piezas industriales con mallas de millones de elementos: el rendimiento absoluto y la paralelización masiva no son objetivos del proyecto. Análisis dinámicos transitorios o problemas multifísicos hoy: están previstos como evolución futura pero no implementados. Cálculos sometidos a normativa de certificación que requieren un programa con trazabilidad de validación industrial: se recomienda Code_Aster, Abaqus o ANSYS según el contexto. Modelado geomecánico avanzado o modelado de hormigón con plasticidad acoplada al daño: existen programas especializados (OpenSees, ATENA, DIANA) más maduros para estos dominios.

**Posición frente a programas de referencia.** Fenix FEM no compite en escala con los programas industriales: complementa el flujo de investigación. La hipótesis es que el código que el investigador puede leer, modificar y validar línea a línea es más valioso, en su contexto, que un programa cerrado del que se obtiene un resultado sin posibilidad de auditoría interna. Cuando el modelo investigado madura y debe aplicarse a producción, la traducción a un programa industrial es la trayectoria natural.
