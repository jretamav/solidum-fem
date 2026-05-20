# Solidum FEM

[![License: LGPL v3](https://img.shields.io/badge/license-LGPL%20v3-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-%E2%89%A5%203.10-blue.svg)](pyproject.toml)
[![Tests](https://img.shields.io/badge/tests-804%20passed-brightgreen.svg)](#validación)
<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.PENDIENTE.svg)](https://doi.org/10.5281/zenodo.PENDIENTE) -->
<!-- [![JOSS](https://joss.theoj.org/papers/PENDIENTE/status.svg)](https://joss.theoj.org/papers/PENDIENTE) -->

Programa de elementos finitos en aproximación de desplazamientos para investigación
en mecánica de sólidos: simulación de problemas mecánico y térmico, acoplados o
desacoplados, con arquitectura optimizada para extender mediante asistencia de IA.

## Estado actual

- **18 elementos**: estructurales 1D (truss, cable, frame 2D Euler/Timoshenko/corotacional,
  frame 3D), sólidos 2D (Quad4, Tri3, Quad8, Quad9, Tri6), sólidos 2D con discontinuidad
  embebida (CST\_Embedded2D), sólidos 3D (Hex8, Tet4).
- **10 materiales**: elásticos 1D/2D/3D, cable unilateral, plasticidad J2 (1D, plane strain,
  plane stress), Drucker-Prager 2D, daño isótropo continuo 1D/2D, cohesivo de daño isótropo.
- **12 solvers**: estático lineal y no lineal, arc-length cilíndrico (Crisfield) y por
  disipación (Gutiérrez 2004 + switching automático), modal con shift-invert ARPACK,
  transitorio Newmark/HHT lineal y no lineal, diferencia central, armónico, espectro
  de respuesta (SRSS, CQC).
- **12 ADRs aceptados** documentando las decisiones arquitecturales.
- **32 specs `validated`** con criterios de aceptación numéricos.
- **804 tests** verdes incluyendo 5 benchmarks canónicos publicados (Lamé 2D, NAFEMS LE1,
  MacNeal-Harder, Bathe wave, Hill 1950).

Estado completo en [docs/STATUS.md](docs/STATUS.md). Onboarding para sesiones nuevas en
[docs/ONBOARDING.md](docs/ONBOARDING.md).

## Instalación

```bash
git clone https://github.com/jretamav/solidum-fem.git
cd solidum-fem
pip install -e .
```

Requiere Python ≥ 3.10. Dependencias: `numpy`, `scipy`, `numba`, `pyyaml`, `meshio`.

Para correr la suite de tests:

```bash
pip install -e .[dev]
pytest -q
```

## Ejemplo rápido

Modelo en YAML — placa cuadrada elástica con desplazamiento impuesto:

```yaml
# minimal.yaml
nodes:
  - {id: 1, coords: [0.0, 0.0]}
  - {id: 2, coords: [1.0, 0.0]}
  - {id: 3, coords: [1.0, 1.0]}
  - {id: 4, coords: [0.0, 1.0]}

materials:
  - {id: 1, type: Elastic2D, E: 210e9, nu: 0.3, plane_state: stress}

elements:
  - {id: 1, type: Quad4, material: 1, thickness: 0.01, nodes: [1, 2, 3, 4]}

boundary_conditions:
  - {node_id: 1, ux: 0.0, uy: 0.0}
  - {node_id: 4, ux: 0.0, uy: 0.0}
  - {node_id: 2, ux: 0.001}
  - {node_id: 3, ux: 0.001}

solver:
  type: LinearSolver
```

Ejecución desde Python:

```python
import solidum

result = solidum.run_yaml("minimal.yaml")
print(result.U)            # desplazamientos nodales
print(result.reactions)    # reacciones en nodos restringidos
```

Más ejemplos en [examples/](examples/): plasticidad 2D, frames 2D/3D, modal,
transitorio (Newmark, diferencia central), armónico, espectro de respuesta.

## Documentación

- **Reference manual** (`manuals/Reference_manual.pdf`) — especificación formal de cada
  componente: ecuaciones, matrices B, contratos YAML, criterios de aceptación.
- **User manual** (`manuals/User_manual.pdf`) — guía de uso: sintaxis YAML, capítulos
  por tipo de análisis, ejemplos completos.
- **Architecture manual** (`manuals/Architecture_manual.pdf`) — visión arquitectural:
  capas, bloques funcionales, ADRs, evolución prevista.

Catálogos navegables cortos:

- [docs/catalogo\_elementos.md](docs/catalogo_elementos.md)
- [docs/catalogo\_materiales.md](docs/catalogo_materiales.md)
- [docs/catalogo\_solvers.md](docs/catalogo_solvers.md)

ADRs en [docs/adr/](docs/adr/). Specs por componente en [docs/specs/](docs/specs/).

## Validación

La suite cubre desde tests unitarios hasta benchmarks externos con valor canónico:

- Tracción uniaxial Lamé 2D plane strain (12 tests, convergencia O(h²) Q8).
- NAFEMS LE1 elliptic membrane (10 tests, Q4/Q8/Tri3/Tri6 con malla 32×32).
- MacNeal-Harder slender beam (8 tests, Frame Euler/Timoshenko exacto con 1 elemento).
- Bathe wave propagation con diferencia central (4 tests, *c\_num* ≈ *c* analítico).
- Hill 1950 cilindro elastoplástico J2 (4 tests, contra solución cerrada).
- Cubo Lamé 3D y MacNeal-Harder 3D para Hex8/Tet4 (Etapa 7).

Detalle en [tests/validation/README.md](tests/validation/README.md).

## Convenciones del proyecto

Documentadas en [Reglas.md](Reglas.md). Resumen:

- Convención única de signos para esfuerzos internos y giros (RHR 3D / sagging-positivo 2D).
- Voigt 2D `[ε_xx, ε_yy, γ_xy]`, Voigt 3D `[ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz]` (ADR 0012).
- Docstrings y mensajes en español; identificadores de código en inglés.

## Citación

Pendiente: cuando se publique el paper en
*Journal of Open Source Software* y el release `v1.0.0` se archive en Zenodo, esta
sección referenciará el archivo [CITATION.cff](CITATION.cff) y el DOI permanente.

## Licencia

GNU Lesser General Public License v3.0 (LGPL-3.0). Ver [LICENSE](LICENSE) (LGPL-3.0)
y [LICENSE.GPL](LICENSE.GPL) (GPL-3.0, requerida porque la LGPL-3.0 se define como
permisos adicionales sobre la GPL-3.0).

## Autor

**J. Retama-Velasco** — [ORCID 0000-0001-6451-5597](https://orcid.org/0000-0001-6451-5597)
Civil Engineering Department, Facultad de Estudios Superiores Aragón
Universidad Nacional Autónoma de México (UNAM)
