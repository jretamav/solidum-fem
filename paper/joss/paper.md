---
title: 'Solidum FEM: An extension-oriented finite element framework for research in solid mechanics'
tags:
  - Python
  - finite element method
  - computational solid mechanics
  - nonlinear mechanics
  - embedded discontinuity
  - plasticity
  - damage
authors:
  - name: J. Retama-Velasco
    orcid: 0000-0001-6451-5597
    affiliation: 1
affiliations:
  - name: Civil Engineering Department, Facultad de Estudios Superiores Aragón, Universidad Nacional Autónoma de México (UNAM), Ciudad Nezahualcóyotl, Estado de México, Mexico
    index: 1
date: 20 May 2026
bibliography: paper.bib
---

# Summary

Solidum FEM is a Python finite element library for research in computational solid mechanics. It supports linear and nonlinear static analysis, modal and transient dynamics (Newmark, HHT-α, central differences), harmonic response, and seismic response-spectrum analysis (SRSS, CQC). Its current catalogue comprises 18 elements (1D structural, 2D and 3D solids, including a 2D solid with embedded discontinuity for fracture), 10 material models (elastic, J2 plasticity, Drucker-Prager, isotropic damage with softening, traction-jump cohesive), and 12 solvers — including a dissipation-based arc-length continuation method.

The library is deliberately structured so that AI coding assistants, working alongside a domain researcher, can add new elements, materials, and solvers without modifying the framework's core. Each component is a self-registering Python module guided by a structured specification document that captures the physical formulation, the software contract, and the acceptance tests in one place. Architectural decisions are persisted as ADRs (Architecture Decision Records), so a fresh agent session can re-acquire the design rationale without dialogue. The catalogue is validated quantitatively against canonical FEM benchmarks — Lamé, NAFEMS LE1, MacNeal-Harder, Bathe wave propagation, and Hill — and includes formulations from the author's own research that are not commonly available off-the-shelf in existing Python frameworks, in particular embedded strong discontinuities [@Retama:2010; @Retama:2020] and dissipation-based arc-length continuation [@Gutierrez:2004; @Verhoosel:2009].

# Statement of need

Computational solid mechanics research routinely produces formulations — return-mapping algorithms for new yield criteria, enriched element kinematics for fracture, custom path-following solvers — that need to be implemented and validated as part of the publication of the underlying theory. Researchers face two unappealing paths. General-purpose finite element frameworks in Python — FEniCS [@Logg:2012; @Baratta:2023], SfePy [@Cimrman:2019], scikit-fem [@Gustafsson:2020] — define problems through their weak form and scale to large simulations, but their abstractions are far from the level at which constitutive models and elements are formulated; the researcher ends up rewriting the formulation to fit the framework's interfaces, which in practice means solving a problem subtly different from the one originally posed. Research-oriented frameworks such as OpenSees [@McKenna:2010; @Zhu:2018] and FEAP [@Taylor:2020] provide well-defined extension points — abstract base classes for new elements and materials in OpenSees, user-element and user-material subroutines in FEAP — so that researchers can plug new formulations into a stable kernel; the kernel itself, however, is not intended to be reshaped by the end-user, and modifying it requires navigating a large, design-pattern-heavy C++ or Fortran codebase whose architecture is fixed by its maintainers. This raises the entry cost for formulations that do not fit the existing extension contracts. In both cases, the software constrains the science it is supposed to serve.

Solidum FEM inverts that relationship: the framework is structured so that new elements, materials and solvers can be added without touching the core, and the formulation drives the implementation rather than being reshaped by it. The primary unit of extension is an element-level Python class (with explicit shape functions, strain–displacement matrix, and integration scheme) or a material-level class (with explicit stress and tangent calculation). Each component is a self-registering module, paired with a structured specification document (`docs/specs/<Name>.md`) that fixes the physical formulation, the parameters, and the acceptance criteria in one place, and with an automated test against an analytical solution or a published benchmark.

The framework is also deliberately designed for collaboration with AI coding assistants. Declarative contracts, self-registering modules, per-component specifications, and persisted Architecture Decision Records together allow an AI agent — operating in independent sessions with no memory of past dialogue — to re-acquire enough context to implement and validate a new component from its specification alone. Solidum FEM has been developed primarily through this workflow: the author defines formulations and acceptance criteria, the AI implements, the author reviews diffs and validates results against benchmarks. The architecture is the artefact that makes this division of labour viable.

The catalogue reflects this division of labour in practice. It includes formulations from the author's research that are not commonly available as off-the-shelf components in existing Python frameworks — embedded strong discontinuities for fracture, with a cohesive traction–jump model condensed locally onto a CST element [@Retama:2010; @Retama:2020]; dissipation-based arc-length continuation with automatic switching between cylindrical and dissipation constraints [@Gutierrez:2004; @Verhoosel:2009]; time-domain integration via Hilber–Hughes–Taylor α [@Hilber:1977] — non-trivial components whose implementation under the AI-assisted workflow serves as a working test of the methodology proposed here. The catalogue is validated against canonical benchmarks of the FEM literature — Lamé closed-form solutions, NAFEMS LE1 elliptic membrane [@NAFEMS:1990], MacNeal–Harder slender beam [@MacNeal:1985], Bathe wave propagation [@Bathe:1996], and Hill's 1950 elastoplastic cylinder [@Hill:1950] — with 38 dedicated validation tests in `tests/validation/` and 804 tests in total covering both software contracts and physical correctness.

# Acknowledgements

This work was made possible by the institutional support of the Universidad Nacional Autónoma de México (UNAM) through the PAPIIT project IT101926 and the PAPIME project PE101626.

# Statement on the use of generative artificial intelligence

During the development of Solidum FEM, the author used Anthropic's Claude language models (Claude Sonnet and Claude Opus, accessed through the Claude Code command-line interface) as coding assistants for the implementation of elements, materials, solvers, tests, and documentation, following formulations, architectural decisions, and acceptance criteria defined by the author in the project's rules document, Architecture Decision Records, and per-component specifications. After each assisted change, the author reviewed the resulting diffs and validated the behaviour against analytical solutions and published benchmarks, taking full responsibility for the scientific rigor, mathematical accuracy, and physical integrity of the results.

# References
