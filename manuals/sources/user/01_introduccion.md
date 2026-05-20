# Introducción y Configuración del Entorno

## Filosofía del Proyecto

Solidum FEM nace como herramienta de investigación en mecánica del medio continuo y análisis no lineal de estructuras. Su diseño persigue dos objetivos simultáneos:

1. **Rigor numérico y físico.** Cada formulación incorporada está validada contra solución analítica o benchmark establecido. La separación estricta `Elemento` ↔ `Material` permite combinar cualquier cinemática con cualquier ley constitutiva siempre que sus protocolos sean compatibles.
2. **Crecimiento sin fricción.** La arquitectura está optimizada para que añadir un elemento, material o solver nuevo sea barato. Auto-registro vía decoradores, contratos declarativos (`STRAIN_DIM`, `DOF_NAMES`) y validación temprana del input son los mecanismos centrales.

### Modelo Data-Driven

Toda la orquestación de un análisis se concentra en un archivo `.yaml` legible por humanos. El usuario no necesita escribir código Python para definir un modelo; la API programática queda reservada para casos avanzados (optimización paramétrica, mallas generadas proceduralmente, acoplamiento con otros sistemas).

## Requisitos del Sistema

Solidum FEM se ejecuta sobre Python 3.8 o superior en Windows, Linux y macOS sin compilación nativa.

### Dependencias

```bash
pip install numpy scipy pyyaml meshio numba
```

- **`numpy` / `scipy`**: álgebra lineal densa y dispersa; `spsolve` sobre matrices CSR.
- **`pyyaml`**: análisis sintáctico del archivo de entrada.
- **`meshio`**: importación de mallas Gmsh y exportación a VTK.
- **`numba`**: compilación JIT de los kernels críticos en elementos 2D y plasticidad.

## Herramientas Externas Recomendadas

- **Gmsh**: generador de mallas; produce archivos `.msh` consumibles directamente desde el YAML.
- **ParaView**: visualización interactiva de los archivos `.vtu` generados.

## Ejecución de un Análisis

El runner principal vive en `examples/ejecutar_yaml.py`. Para ejecutar un modelo:

```bash
python examples/ejecutar_yaml.py ruta/al/modelo.yaml
```

El script carga el YAML, valida su contenido (acumulando *todos* los errores en una sola pasada antes de abortar), construye el dominio, ensambla el solver pedido, ejecuta el análisis y exporta los resultados según el bloque `output`.
