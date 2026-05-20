# Flujo de un análisis

Esta sección describe, paso a paso, la secuencia de operaciones que se ejecuta sobre un caso desde su lanzamiento hasta la generación de resultados. La descripción es de alto nivel y omite detalles de implementación.

## 1. Arranque

El caso se lanza desde un *script* de usuario o desde la interfaz pública de entrada (`solidum.run` / `solidum.run_yaml`). Antes de que se abra el archivo YAML, se ejecuta la importación del paquete `solidum`. Dicha importación dispara el procedimiento de descubrimiento automático (*autodiscover*), que recorre las carpetas canónicas del proyecto e importa todos los módulos. Como cada material, elemento y solver lleva el decorador `@register` correspondiente, los registros quedan poblados sin enumeración manual.

## 2. Lectura del caso

El intérprete del caso abre el archivo YAML y construye los objetos del problema. Para cada material, elemento o solver mencionado en el archivo, se localiza su nombre en el registro correspondiente, se recupera la clase, se introspeccionan los argumentos requeridos por su constructor y se le transmiten los campos del YAML. El resultado es un grafo de objetos: nodos, materiales, elementos asociados a sus materiales, condiciones de contorno, cargas y un solver seleccionado.

En esta etapa se ejecuta la primera comprobación de consistencia: cuando un elemento se construye, su clase base verifica que el material recibido sea dimensionalmente compatible. Un elemento `Truss2D` con material `Elastic2D` produce un error en la fase de construcción del caso, no doscientas iteraciones después en forma de fallo opaco durante la solución.

## 3. Numeración de grados de libertad

El objeto `Domain` recorre nodos y elementos. Cada elemento declara los DOF requeridos en cada uno de sus nodos mediante el atributo `DOF_NAMES`. El dominio asigna índices globales de forma consistente con esa información y obtiene así el tamaño y la topología de la matriz de rigidez global. Esta numeración rige el ensamblaje durante toda la simulación.

## 4. Bucle del solver

El solver toma el control. Su comportamiento depende del tipo seleccionado (lineal, no lineal con Newton-Raphson, longitud de arco), pero todos los tipos comparten una estructura común: una secuencia de pasos de carga, una o varias iteraciones dentro de cada paso y, dentro de cada iteración, la siguiente secuencia interna:

1. Se solicita a cada elemento su matriz de rigidez tangente local `K_local` y su vector de fuerzas internas `f_internal`. El elemento las obtiene a partir del estado *trial* de sus puntos de Gauss y mediante una llamada al material con la deformación correspondiente.
2. El ensamblador construye la matriz global `K_global` y el vector global utilizando topología COO en caché: en la primera invocación se calculan los pares de índices; en las invocaciones siguientes se reescribe únicamente el vector de datos sobre la misma topología.
3. Se imponen las condiciones de Dirichlet por eliminación directa (ADR 0004). Toda restricción se expresa en forma afín `u_s = g_s + Σ α_si · u_mi` y se acumula en un `ConstraintSet`; el ensamblador produce el par `(T, g)` tal que `u = T · u_libre + g` y entrega al solver el sistema reducido `K_red = TᵀKT`, `F_red = Tᵀ(F − K · g)`. La imposición es exacta a redondeo y preserva la simetría y la positividad de la matriz original.
4. El subsistema algebraico resuelve el sistema reducido `K_red · δU_libre = F_red`. Aquí interviene el despachador (ADR 0003): se inspeccionan las propiedades de `K_red` (simetría, posible definición positiva) y se selecciona el algoritmo numérico adecuado (factorización de Cholesky cuando aplica, factorización LU en el caso general). El ensamblador reconstruye después el incremento completo mediante `expand`.
5. Se evalúa la convergencia mediante un criterio dual: norma del residuo de fuerza y norma del incremento de desplazamientos, ambas comparadas contra una tolerancia de la forma `atol + rtol · escala(estado)` con la escala adaptativa al estado corriente del problema (ADR 0007). Los términos absolutos `atol` se autoderivan de las escalas del primer ensamblaje, lo que garantiza invariancia bajo cambio de unidades: el mismo problema planteado en N/m, kN/mm o MN/m converge en el mismo número de iteraciones hasta paridad de bits. El residuo de fuerza se evalúa únicamente sobre los DOF libres, dado que en los DOF prescritos el residuo contiene la reacción física, no un fallo de equilibrio. La convergencia del paso requiere que ambos criterios se cumplan simultáneamente (semántica AND, no OR).

## 5. Cierre del paso

Tras la convergencia de un paso, el solver invoca `commit_state()` sobre cada elemento. Esta operación promueve el estado *trial* a estado comprometido: las variables internas (deformaciones plásticas, daño, historia) que los puntos de Gauss han recorrido durante las iteraciones se consolidan. Si el paso no convergiera y se redujese el incremento, el estado *trial* se descartaría sin contaminar el estado comprometido.

## 6. Salida

A la finalización del solver, el objeto `Domain` construye un `SolveResult` inmutable que contiene los desplazamientos finales, las cargas aplicadas, las reacciones y los esfuerzos internos por elemento en la convención de signos del proyecto. Este `SolveResult` constituye la interfaz pública (Application Programming Interface, API) que consume cualquier herramienta externa (GUI, *scripts* de post-proceso). De forma paralela, el módulo `VtkExporter` escribe un archivo VTK con los desplazamientos, las tensiones y la variable principal del material declarada por cada material mediante el atributo `PRIMARY_STATE_VAR`.

## Visión sintética

A grandes rasgos, la secuencia es: importación → descubrimiento automático → interpretación del caso → numeración de DOF → bucle del solver { ensamblaje → resolución algebraica → evaluación de convergencia } → consolidación → salida. Cada flecha corresponde a un contrato declarado y verificado; cada paso intermedio es sustituible sin afectar al resto.
