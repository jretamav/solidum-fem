# Ejemplos de Solidum FEM

Modelos ejecutables que demuestran las capacidades del programa. A diferencia
de `tests/`, que verifica corrección, estos existen para **ser leídos y
modificados** por quien aprende a usar Solidum.

---

## Convención de organización

Decidida el 2026-09-10 al crecer el directorio a ~28 entradas planas, donde los
seis archivos de un mismo ejemplo quedaban intercalados alfabéticamente con
YAMLs sueltos de temas distintos.

**La regla es una sola:**

> Un ejemplo con **más de un archivo** vive en su propia carpeta, con un
> `README.md` corto. Un ejemplo de **un solo archivo** vive suelto en
> `basicos/`.

La carpeta gana cuando hay piezas que pertenecen juntas —malla `.geo`, modelo
`.yaml`, scripts de post-proceso, figuras generadas—. Con un único YAML sólo
añade un nivel de navegación sin agrupar nada.

```
examples/
├── README.md                  ← este índice
├── ejecutar_yaml.py           ← runner genérico, transversal a todos
├── api_resultados_minimal.py  ← demo mínima de la API de resultados
│
├── basicos/                   ← un archivo = un ejemplo
│   └── modelo_*.yaml
│
└── <nombre_del_ejemplo>/      ← ejemplo compuesto
    ├── README.md              ← qué demuestra, cómo correrlo, qué esperar
    ├── modelo.yaml
    ├── malla.geo
    ├── figura_*.py
    └── figura_*.png
```

### Al crear un ejemplo nuevo

1. Si genera figuras, necesita malla o tiene scripts de post-proceso →
   **carpeta propia**.
2. `README.md` dentro, con: qué capacidad demuestra, comando para ejecutarlo,
   y qué resultado esperar (una cifra citable si la hay).
3. Añadir una fila a la tabla de abajo.
4. Los `.png` generados **se versionan** — el manual y la charla los consumen,
   y regenerarlos exige matplotlib + malla, que no todo lector tiene.

### Migración de lo existente

Los ejemplos actuales siguen planos **deliberadamente**. La convención se
aplica desde los nuevos y se migrará lo anterior cuando demuestre valer la
pena, siguiendo la regla de los dos casos reales de `Reglas.md §1`. La
migración no es gratuita: [`tests/test_examples_yaml.py`](../tests/test_examples_yaml.py)
construye rutas planas (`_path('modelo_marco.yaml')`) y habría que ajustarlo,
igual que las referencias en [`docs/ONBOARDING.md`](../docs/ONBOARDING.md) y en
`manuals/sources/user/01_introduccion.md`.

(`tests/test_entry.py` recorre el directorio con `os.walk`, así que ése ya
soporta subdirectorios sin cambios.)

---

## Índice de ejemplos

### Compuestos

| Ejemplo | Qué demuestra |
|---|---|
| [`bambu_ortotropo/`](bambu_ortotropo/) | Acoplamiento tracción–cortante en un material ortótropo: una tira de bambú con la fibra fuera de eje se distorsiona al traccionarla, y su rigidez cae por debajo del 10 % de E₁ a 45°. Usa `Orthotropic2D`. |

### Modelos sueltos (pendientes de migrar)

| Archivo | Análisis | Componentes |
|---|---|---|
| `modelo_elastico_2d.yaml` | Estático lineal | `Quad4` + `Elastic2D` + `LinearSolver` |
| `modelo_marco.yaml` | Estático no lineal | `Frame2DEulerCorot` + `NonlinearSolver` |
| `modelo_plasticidad.yaml` | Estático no lineal | plasticidad J2 |
| `modelo_placa.yaml` | Estático lineal | sólido 2D |
| `placa_gmsh.yaml` | Estático lineal | malla importada de gmsh (`placa.geo`, `placa.msh`) |
| `modelo_modal.yaml` | Modal | `ModalSolver` |
| `modelo_dinamico_plastico.yaml` | Transitorio no lineal | `NewtonNewmarkSolver` |
| `modelo_central_difference.yaml` | Transitorio explícito | `CentralDifferenceSolver` |
| `modelo_harmonic.yaml` | Armónico | `HarmonicSolver` |
| `modelo_response_spectrum.yaml` | Espectral | `ResponseSpectrumSolver` |
| `charla_cook.*` | Cook's membrane J2 | conjunto de 6 archivos — candidato natural a carpeta |
| `charla_placa_perfecta.yaml` | plasticidad perfecta | + `charla_placa.geo` |

---

## Cómo ejecutar

```bash
python examples/ejecutar_yaml.py examples/<ruta>/<modelo>.yaml
```

o desde Python:

```python
import solidum
resultado = solidum.run_yaml('examples/<ruta>/<modelo>.yaml')
```

## Garantías de la suite

[`tests/test_examples_yaml.py`](../tests/test_examples_yaml.py) protege este
directorio en tres niveles:

1. **Cifras fijadas** para los ejemplos históricos, enumerados a mano.
2. **`TestTodosLosEjemplosEjecutan`** descubre los YAML de forma **recursiva**
   —subdirectorios incluidos— y verifica que todos siguen ejecutándose de
   extremo a extremo. Es lo que hace segura la organización en carpetas: un
   ejemplo nuevo entra en la regresión sin que nadie lo registre.
3. **`test_los_ejemplos_compuestos_tienen_readme`** comprueba que toda carpeta
   de ejemplo lleva su `README.md`, para que la convención no se erosione.

Excepción declarada: los YAML que consumen una malla `.msh` **no versionada**
(los `charla_*`, que requieren correr gmsh sobre su `.geo`) se saltan con
motivo explícito en vez de fallar. Distinguir "falta una dependencia externa"
de "el ejemplo está roto" evita que el test se vuelva ruido y acabe
silenciado.
