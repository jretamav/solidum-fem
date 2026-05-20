# Solidum FEM

[![License: LGPL v3](https://img.shields.io/badge/license-LGPL%20v3-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-%E2%89%A5%203.10-blue.svg)](pyproject.toml)
[![Tests](https://img.shields.io/badge/tests-804%20passed-brightgreen.svg)](#validation)
<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.PENDING.svg)](https://doi.org/10.5281/zenodo.PENDING) -->
<!-- [![JOSS](https://joss.theoj.org/papers/PENDING/status.svg)](https://joss.theoj.org/papers/PENDING) -->

A displacement-based finite element framework for research in solid mechanics,
with an architecture optimized for extension via AI-assisted development.

## Current capabilities

- **18 elements**: 1D structural (truss, cable, frame 2D Euler/Timoshenko/corotational,
  frame 3D), 2D solids (Quad4, Tri3, Quad8, Quad9, Tri6), 2D solids with embedded
  discontinuity (CST\_Embedded2D), 3D solids (Hex8, Tet4).
- **10 materials**: elastic 1D/2D/3D, unilateral cable, J2 plasticity (1D, plane
  strain, plane stress), 2D Drucker-Prager, isotropic continuum damage 1D/2D,
  isotropic cohesive damage.
- **12 solvers**: linear and nonlinear static, cylindrical arc-length (Crisfield)
  and dissipation arc-length (Gutiérrez 2004 with automatic switching), modal
  via shift-invert ARPACK, linear and nonlinear Newmark/HHT time integration,
  central difference, harmonic response, response spectrum (SRSS, CQC).
- **12 accepted ADRs** documenting architectural decisions.
- **32 validated specs** with quantitative acceptance criteria.
- **804 tests** green, including 5 published canonical benchmarks (Lamé 2D,
  NAFEMS LE1, MacNeal-Harder, Bathe wave propagation, Hill 1950).

Project status in [docs/STATUS.md](docs/STATUS.md) (in Spanish; see the
[Documentation](#documentation) section for the language policy).

## Installation

```bash
git clone https://github.com/jretamav/solidum-fem.git
cd solidum-fem
pip install -e .
```

Requires Python ≥ 3.10. Dependencies: `numpy`, `scipy`, `numba`, `pyyaml`,
`meshio`.

To run the test suite:

```bash
pip install -e .[dev]
pytest -q
```

## Quick example

YAML model — square elastic plate with imposed displacement:

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

Execution from Python:

```python
import solidum

result = solidum.run_yaml("minimal.yaml")
print(result.U)            # nodal displacements
print(result.reactions)    # reactions at constrained nodes
```

Further examples in [examples/](examples/): 2D plasticity, 2D/3D frames, modal
analysis, time integration (Newmark, central difference), harmonic response,
response spectrum analysis.

## Documentation

The project documentation follows a layered language policy:

| Layer | Language | Artifacts |
|---|---|---|
| Public-facing | **English** | This README, `CITATION.cff`, JOSS paper, contributor docs |
| Manuals | **Bilingual (ES + EN)** | Reference, User, Architecture (English versions auto-generated from Spanish source via a controlled glossary) |
| Navigational | **Spanish** | `docs/STATUS.md`, `docs/ROADMAP.md`, `docs/ONBOARDING.md`, `docs/MATRIZ.md`, catalogs |
| Internal / technical | **Spanish** | Component specs, ADRs, docstrings, code comments, governance docs |
| Code identifiers | **English** | Classes, functions, variables, file names |

The three manuals (in `manuals/`):

- **Reference manual** — formal specification of each component: equations,
  $\mathbf{B}$ matrices, YAML contracts, acceptance criteria.
- **User manual** — usage guide: YAML syntax, chapters by analysis type,
  complete examples.
- **Architecture manual** — architectural overview: layers, functional blocks,
  ADRs, planned evolution.

Component catalogs (short navigational entries):

- [docs/catalogo\_elementos.md](docs/catalogo_elementos.md) (elements)
- [docs/catalogo\_materiales.md](docs/catalogo_materiales.md) (materials)
- [docs/catalogo\_solvers.md](docs/catalogo_solvers.md) (solvers)

ADRs live in [docs/adr/](docs/adr/). Per-component specs in
[docs/specs/](docs/specs/).

## Validation

The test suite spans from unit tests to external benchmarks with canonical
reference values:

- Lamé 2D plane strain uniaxial tension (12 tests, O(h²) convergence for Q8).
- NAFEMS LE1 elliptic membrane (10 tests, Q4/Q8/Tri3/Tri6 on a 32×32 mesh).
- MacNeal-Harder slender beam (8 tests, Euler/Timoshenko frame exact with one
  element).
- Bathe wave propagation with central difference (4 tests, *c\_num* ≈ *c*
  analytical).
- Hill 1950 elastic-plastic J2 cylinder (4 tests, against the closed-form
  solution).
- 3D Lamé cube and 3D MacNeal-Harder for Hex8/Tet4 (Stage 7).

Details in [tests/validation/README.md](tests/validation/README.md).

## Project conventions

Documented in [Reglas.md](Reglas.md) (in Spanish). Highlights:

- Unified sign convention for internal forces and rotations (3D right-hand
  rule; 2D sagging-positive bending).
- Voigt 2D `[ε_xx, ε_yy, γ_xy]`; Voigt 3D `[ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz]`
  (per ADR 0012).
- Docstrings and end-user messages in Spanish; code identifiers in English.

## How to cite

Pending: once the paper is published in the *Journal of Open Source Software*
and the `v1.0.0` release is archived on Zenodo, this section will point to
[CITATION.cff](CITATION.cff) and the permanent DOI.

## License

GNU Lesser General Public License v3.0 (LGPL-3.0). See [LICENSE](LICENSE)
(LGPL-3.0) and [LICENSE.GPL](LICENSE.GPL) (GPL-3.0, required because the
LGPL-3.0 is defined as a set of additional permissions on top of the GPL-3.0).

## Author

**J. Retama-Velasco** — [ORCID 0000-0001-6451-5597](https://orcid.org/0000-0001-6451-5597)
Civil Engineering Department, Facultad de Estudios Superiores Aragón
Universidad Nacional Autónoma de México (UNAM)
