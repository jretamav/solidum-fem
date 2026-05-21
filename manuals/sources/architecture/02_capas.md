# Estructura en capas

El programa se organiza en seis capas con responsabilidades disjuntas. Cada capa se comunica con las contiguas mediante contratos explícitos; ninguna capa accede al interior de otra que no sea la inmediatamente inferior.

## Mapa de capas

```latex
\begin{figure}[h]
\centering
\begin{tikzpicture}[
  node distance = 5mm and 0mm,
  layer/.style = {
    rectangle, draw=darkgray, thick, rounded corners=2pt,
    fill=bglight, text width=12.5cm, align=left,
    inner sep=6pt, font=\small
  },
  title/.style = {font=\small\bfseries\color{yamlkey}},
  arr/.style   = {-{Stealth[length=3mm]}, thick, draw=darkgray}
]
\node[layer] (entrada) {
  {\title Entrada del usuario}\\
  \texttt{case.yaml} (descripción declarativa del problema)\\
  Mallas \texttt{.msh} (geometría desde Gmsh)
};
\node[layer, below=of entrada] (init) {
  {\title Capa de inicialización} \hfill \texttt{solidum/\_\_init\_\_.py}\\
  Al importar el paquete se invoca \texttt{autodiscover.initialize()};
  los decoradores \texttt{@*Registry} pueblan
  \texttt{MaterialRegistry}, \texttt{ElementRegistry} y
  \texttt{SolverRegistry}.
};
\node[layer, below=of init] (parse) {
  {\title Capa de interpretación del caso} \hfill \texttt{solidum/utils/yaml\_parser.py}\\
  Lee el archivo YAML y construye objetos consultando los registros
  por introspección del constructor.
};
\node[layer, below=of parse] (dominio) {
  {\title Capa de dominio} \hfill \texttt{solidum/core/}\\
  \texttt{Domain} -- \texttt{Node} -- DOF (numeración global). Bases
  \texttt{Element} (con \texttt{DOF\_NAMES}, \texttt{STRAIN\_DIM},
  \texttt{N\_INTEGRATION\_POINTS}, \texttt{ElementState} en estados
  \emph{trial} y comprometido) y \texttt{Material} (con
  \texttt{STRAIN\_DIM}, \texttt{PRIMARY\_STATE\_VAR}).
};
\node[layer, below=of dominio] (numerica) {
  {\title Capa numérica} \hfill \texttt{solidum/math/}\\
  \texttt{assembly.py} (ensamblaje disperso con caché COO),
  \texttt{integration.py} (cuadraturas de Gauss),
  paquete \texttt{solvers/} con un módulo por solver
  (\texttt{LinearSolver}, \texttt{NonlinearSolver},
  \texttt{ArcLengthSolver}, \texttt{ModalSolver},
  \texttt{NewmarkSolver}, \texttt{NewtonNewmarkSolver}) y subsistema algebraico
  \texttt{linalg/} (\texttt{LUSolver}, \texttt{CholeskySolver},
  \texttt{LDLTSolver} reservado, despachador
  \texttt{select\_solver}, \texttt{EigenSolver}).
};
\node[layer, below=of numerica] (salida) {
  {\title Capa de salida} \hfill \texttt{solidum/utils/vtk\_exporter.py} + \texttt{solidum.results}\\
  Exporta $\mathbf U$, $\sigma$ y \texttt{PRIMARY\_STATE\_VAR} en
  formato VTK. \texttt{SolveResult} expone $\mathbf U$,
  $\mathbf F_\text{aplicada}$ y las fuerzas internas por elemento
  como interfaz pública para consumidores externos.
};

\draw[arr] (entrada) -- (init);
\draw[arr] (init) -- (parse);
\draw[arr] (parse) -- (dominio);
\draw[arr] (dominio) -- (numerica);
\draw[arr] (numerica) -- (salida);
\end{tikzpicture}
\caption{Estructura en capas de Solidum FEM. Cada capa se comunica con la inmediatamente inferior mediante contratos explícitos.}
\end{figure}
```

## Responsabilidad de cada capa

**Entrada del usuario.** Un caso de Solidum FEM se describe de forma declarativa en un archivo en formato YAML (acrónimo recursivo de *YAML Ain't Markup Language*): nodos, materiales, elementos, condiciones de contorno, cargas y solver. La malla puede embeberse en el YAML o leerse desde un archivo `.msh` generado por Gmsh. El YAML constituye contrato público; cualquier consumidor externo (Interfaz Gráfica de Usuario o GUI, *scripts* de automatización) opera contra él.

**Capa de inicialización.** Su única responsabilidad consiste en poblar los registros (`MaterialRegistry`, `ElementRegistry`, `SolverRegistry`) sin que el resto del programa enumere manualmente las clases existentes. Esta tarea se ejecuta importando todos los módulos bajo carpetas canónicas; los decoradores `@register` se ejecutan en el momento de la importación e insertan cada clase en el registro correspondiente.

**Capa de interpretación del caso.** Traduce el archivo YAML en objetos. No contiene ramificaciones del tipo `if material_type == "Elastic1D"` por cada material: introspecciona los argumentos formales del constructor de la clase recuperada del registro y le transmite los campos del YAML. La incorporación de un material nuevo no requiere modificación del intérprete.

**Capa de dominio.** El objeto `Domain` mantiene la numeración global de grados de libertad (DOF) y la conectividad. Las clases base `Element` y `Material` definen el contrato que cumple cada componente concreto. El objeto `ElementState` encapsula las dos copias de las variables internas (tensiones, variables de estado, historia) sobre las que opera el solver no lineal: una copia *trial* y otra comprometida.

**Capa numérica.** Concentra los componentes propios de la mecánica computacional: ensamblaje de matrices dispersas, cuadraturas de integración, estrategia de paso (lineal, no lineal incremental, longitud de arco) y el subsistema algebraico que resuelve los sistemas lineales que aparecen dentro de cada iteración.

**Capa de salida.** Provee dos productos. Por un lado, `VtkExporter` para visualización en programas compatibles con VTK (ParaView). Por otro, `SolveResult` con los desplazamientos, las cargas aplicadas, las reacciones y las fuerzas internas por elemento, concebido para consumidores externos como FenixBAR (ADR 0002).

## Justificación de la separación

La frontera entre capas no es decorativa: define qué cambia con coste bajo y qué cambia con coste alto. La incorporación de un material nuevo afecta exclusivamente a la capa de dominio (`solidum/materials/<nombre>.py`) y no requiere modificar el intérprete del caso, el ensamblador ni el solver. La modificación del algoritmo del ensamblaje disperso no repercute sobre ningún material ni elemento concretos. La estanqueidad entre capas es el mecanismo que mantiene plana la pendiente de coste de extensión.
