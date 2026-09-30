# Contributing to Solidum FEM

Thank you for your interest in Solidum FEM. This document describes how to
report issues, ask for help, suggest changes, and submit code.

> Solidum FEM is research software developed primarily by a single author with
> AI-assisted development. External contributions are welcome, but the project
> moves at the pace of a research workflow rather than a commercial product.
> Please calibrate response-time expectations accordingly.

## Code of conduct

We expect all contributors to follow standard open-source community norms:
respectful communication, focus on technical content, no personal attacks,
and a welcoming stance toward newcomers regardless of background.

## How to report a bug

1. Search the [issue tracker](https://github.com/jretamav/solidum-fem/issues)
   to see whether the bug has already been reported.
2. If not, open a new issue and include:
   - A descriptive title.
   - Steps to reproduce, ideally a minimal YAML model file or a short Python
     snippet.
   - Expected behavior vs. actual behavior, with numeric values, error
     messages, or stack traces.
   - Your environment: Python version, OS, and Solidum version
     (`pip show solidum-fem`).

For numerical discrepancies, please include the analytical or reference value
against which you compared, with its source.

## How to suggest a feature

Open an issue describing:

- The physical or numerical problem you want to solve.
- Why the current capabilities don't cover it.
- A bibliographic reference if the feature corresponds to a published
  formulation.

Feature requests aligned with the project's roadmap (see
[docs/ROADMAP.md](docs/ROADMAP.md), in Spanish) have the highest chance of
being prioritized.

## How to ask for support

For usage questions:

- Consult the [User manual](manuals/User_manual.pdf) (currently Spanish; an
  English version is planned — see the language policy in
  [Reglas.md §6](Reglas.md)) and the runnable [examples/](examples/).
- Search closed issues.
- Open a new issue.

## How to submit a pull request

1. Fork the repository.
2. Create a feature branch from `main`.
3. Make your changes following the project conventions (next section).
4. Run the test suite (`pytest -q`) and verify it passes.
5. Add tests for any new physical formulation (see *Tests* below).
6. Commit with a clear message in the format `type(scope): description` —
   see recent commits for examples.
7. Open a pull request against `main` describing what changes and why.

By submitting a pull request you agree to license your contribution under
the LGPL-3.0-or-later license (the project license).

## Development setup

```bash
git clone https://github.com/jretamav/solidum-fem.git
cd solidum-fem
pip install -e .[dev]
pytest -q
```

Python ≥ 3.10 required. Numba's JIT compilation warms up on first run;
subsequent runs are faster.

## Project conventions

The project follows strict conventions documented in [Reglas.md](Reglas.md)
(in Spanish). Key points for contributors:

- **Sign conventions**: unified across the project — right-hand rule in 3D,
  sagging-positive bending in 2D (Reglas.md §5).
- **Voigt notation**: 2D `[ε_xx, ε_yy, γ_xy]`; 3D `[ε_xx, ε_yy, ε_zz, γ_xy,
  γ_yz, γ_xz]`.
- **Language policy** (Reglas.md §6): code identifiers in English; docstrings
  and code comments in Spanish; public-facing docs (this file, README,
  CITATION) in English. This hybrid is deliberate — see Reglas.md §6 for the
  full rationale.
- **Spec-first workflow**: new elements, materials, or solvers require a spec
  in `docs/specs/` before implementation. The spec sets the physical contract
  (equations, parameters, acceptance criteria) and is validated against the
  resulting tests. Its first step is the admission criterion below, declared
  in the spec's `admission` block (`textbooks` and `general_codes`, written
  with the canonical identifiers of `ADMISSION_TEXTBOOKS` and
  `ADMISSION_GENERAL_CODES` in `solidum/tools/spec.py`).
- **Admission criterion** (Reglas.md §1): a formulation belongs to the main
  program (standard Solidum) only if (i) it appears in the reference texts of
  the field (Zienkiewicz-Taylor, Hughes, Bathe, Belytschko-Liu-Moran,
  de Souza Neto-Perić-Owen, de Borst-Crisfield-Remmers-Verhoosel, Chopra) and (ii) it
  is implemented in at least one established general-purpose code (Abaqus,
  ANSYS, DIANA, Code_Aster, OpenSees). When the criterion is inconclusive,
  the maintainer decides.
- **User modules** (ADR 0020, `docs/adr/0020-modulos-de-usuario.md`): a
  research formulation that does not meet the criterion (for example, an
  embedded strong discontinuity with KOS kinematics, crack-path tracking, or
  Gutiérrez's energy-release path-following control) goes to a user module,
  in the spirit of FEAP's user elements, materials and macros:
  code in `solidum/user/<name>/`, specs in `docs/user/<name>/specs/`, tests in
  `tests/user/<name>/`. A user module is never loaded automatically: a model
  requests it with `user_modules: [<name>]` in its YAML file, or with
  `solidum.load_user_module("<name>")` in Python. The main program never
  imports or names a user module; if a module needs something the main
  program lacks, it is added as a generic mechanism, never as a special case.
  `tests/test_main_program_purity.py` enforces this boundary with three
  checks: a structural one (no main-program file imports a user module in any
  form, whether absolute, relative, lazy, by string, or through
  `load_user_module` with a literal name; only the generic loader
  `load_user_module` / `available_user_modules` is allowed), a lexical one (no
  main-program file mentions what a module provides; the watched names are
  derived automatically from what each module registers when loaded, so a new
  module does not have to be listed in the test) and a runtime one
  (`import solidum` loads no user module and every registered class belongs to
  the main program). A user module may only add registry entries: if it
  reassigns an existing one, `load_user_module` raises `ValueError` and undoes
  the whole load. In the other direction,
  `tests/test_user_module_public_imports.py` checks that a user module uses
  only the public API of the main program (no module path or name with a
  leading underscore).
- **ADRs for architectural decisions**: cross-cutting changes require an
  Architecture Decision Record in `docs/adr/`. See [Reglas.md §4](Reglas.md).

Contributors uncomfortable with Spanish-language internal docs are welcome
to raise this as an issue; we are actively working toward bilingual manuals
(Spanish canonical, English generated; see Reglas.md §6).

## Tests

Solidum FEM has a comprehensive test suite covering unit tests and canonical
benchmark validation:

```bash
pytest -q                            # full suite
pytest tests/validation -q           # canonical benchmark validation only
pytest tests/test_solid_3d.py -q     # a single module
```

Every new physical formulation must come with a validation test against an
analytical solution or a published benchmark. See [tests/validation/](tests/validation/)
for examples (Lamé 2D/3D, NAFEMS LE1, MacNeal-Harder, Bathe wave, Hill 1950).

## Maintainer contact

For matters not appropriate for public issues (security disclosures,
sensitive collaboration inquiries), contact the maintainer at
[jretamav@comunidad.unam.mx](mailto:jretamav@comunidad.unam.mx).

For everything else, the [issue tracker](https://github.com/jretamav/solidum-fem/issues)
is the right channel.

## License

Solidum FEM is licensed under the GNU Lesser General Public License v3.0
(LGPL-3.0-or-later). Contributions are licensed under the same terms.
See [LICENSE](LICENSE) and [LICENSE.GPL](LICENSE.GPL).
