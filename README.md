# Solidum FEM

[![License: LGPL v3](https://img.shields.io/badge/license-LGPL%20v3-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-%E2%89%A5%203.10-blue.svg)](pyproject.toml)
[![Tests](https://github.com/jretamav/solidum-fem/actions/workflows/tests.yml/badge.svg)](https://github.com/jretamav/solidum-fem/actions/workflows/tests.yml)
<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.PENDING.svg)](https://doi.org/10.5281/zenodo.PENDING) -->
<!-- [![JOSS](https://joss.theoj.org/papers/PENDING/status.svg)](https://joss.theoj.org/papers/PENDING) -->

A displacement-based finite element framework for research in solid mechanics,
with an architecture optimized for extension via AI-assisted development.

## Current capabilities

- **23 elements**: 1D structural (truss, cable, frame 2D Euler/Timoshenko/corotational,
  frame 3D), 2D solids (Quad4, Tri3, Quad8, Quad9, Tri6), 2D solids with embedded
  discontinuity (CST\_Embedded2D), 3D solids linear (Hex8, Tet4) and quadratic
  (Hex20, Hex27, Tet10), and thermal conduction elements (Quad4Thermal, Hex8Thermal).
- **14 materials** across three parallel families: elastic 1D/2D/3D, unilateral
  cable, J2 plasticity (1D, plane strain, plane stress, 3D), Drucker-Prager (2D
  and 3D), isotropic continuum damage 1D/2D/3D, isotropic cohesive damage
  (traction-jump), and Fourier thermal conduction with tensor conductivity.
- **13 solvers**: linear and nonlinear static, cylindrical arc-length (Crisfield)
  and dissipation arc-length (Gutiérrez 2004 with automatic switching), modal
  via shift-invert ARPACK, linear and nonlinear Newmark/HHT time integration,
  central difference, harmonic response, response spectrum (SRSS, CQC), and
  θ-method integration for transient heat conduction.
- **Thermal analysis** (uncoupled): steady-state and transient heat conduction
  in 2D and 3D, with Dirichlet and Neumann boundary conditions. The steady-state
  regime needs no dedicated solver — the existing linear solver handles it
  unchanged, since assembly, constraint elimination and algebraic dispatch make
  no assumption about the physical meaning of a degree of freedom.
- **12 accepted ADRs** documenting architectural decisions.
- **49 validated specs** with quantitative acceptance criteria.
- **1176 tests** green, including 8 published canonical benchmarks (Lamé 2D and
  3D, NAFEMS LE1 and LE10, MacNeal-Harder 2D and 3D, Bathe wave propagation,
  Hill 1950, Carslaw-Jaeger semi-infinite solid).

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
  - {id: 1, type: Elastic2D, E: 210e9, nu: 0.3, hypothesis: plane_stress}

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
print(result.U)                    # nodal displacements
print(result.reactions_by_node)    # reactions at constrained nodes
```

The same entry point runs every analysis type; the solver declared in the YAML
determines the pipeline and the result type. A thermal model, for instance,
declares its material under `thermal_materials` and its temperatures as `T`
boundary conditions. In the steady state it returns a `SolveResult` whose `U`
holds nodal temperatures; the transient solver below returns a
`ThermalTransientResult` with the full time history:

```yaml
thermal_materials:
  - {id: 1, type: ThermalConduction, k: 45.0, c: 460.0, density: 7850.0}

elements:
  - {id: 1, type: Quad4Thermal, material: 1, thickness: 0.3, nodes: [1, 2, 3, 4]}

boundary_conditions:
  - {node_id: 1, T: 100.0}
  - {node_id: 2, T: 20.0}

solver:
  type: ThetaMethodSolver    # or LinearSolver for the steady state
  dt: 60.0
  n_steps: 500
  T_initial: 20.0
```

Further examples in [examples/](examples/): 2D plasticity, 2D/3D frames, modal
analysis, time integration (Newmark, central difference), harmonic response,
response spectrum analysis.

## Documentation

The project documentation follows a layered language policy:

| Layer | Language | Artifacts |
|---|---|---|
| Public-facing | **English** | This README, `CITATION.cff`, JOSS paper, contributor docs |
| Manuals | **Spanish today; bilingual planned** | Reference, User, Architecture. English versions are to be generated from the Spanish source via a controlled glossary; **neither the glossary nor the English PDFs exist yet** |
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
- 3D Lamé cube and 3D MacNeal-Harder for Hex8/Tet4 (Stage 7), extended to the
  quadratic elements Hex20/Hex27/Tet10 (Stage A.ter).
- NAFEMS LE10 thick plate under pressure, against the canonical σ\_yy(D) =
  −5.38 MPa for Hex20 and Hex27.
- Lamé thick cylinder in 3D, against the closed-form Timoshenko-Goodier
  solution evaluated at every Gauss point — the quantitative demonstration of
  curved isoparametric capability for the quadratic elements.
- Carslaw-Jaeger semi-infinite solid for transient heat conduction, against the
  error-function solution.

Thermal time integration is verified by measured convergence rates rather than
a single tolerance: 0.99 for backward Euler and 2.00 for Crank-Nicolson, taken
against the exact solution of the semidiscrete system so that spatial error does
not mask the temporal rate.

Details in [tests/validation/README.md](tests/validation/README.md).

## Project conventions

Documented in [Reglas.md](Reglas.md) (in Spanish). Highlights:

- Unified sign convention for internal forces and rotations (3D right-hand
  rule; 2D sagging-positive bending).
- Voigt 2D `[ε_xx, ε_yy, γ_xy]`; Voigt 3D `[ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz]`
  (per ADR 0012). The thermal problem uses **no Voigt notation**: ∇T is a
  genuine vector, not a compressed symmetric tensor, so the thermal
  $\mathbf{B}$ matrix is the raw shape-function gradient.
- Prescribed heat flux is **positive outwards** (cooling); an undeclared edge
  or face is adiabatic, by analogy with the traction-free boundary.
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
