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

The library is deliberately structured so that AI coding assistants, working alongside a domain researcher, can add new elements, materials, and solvers without modifying the framework's core. Each component is a self-registering Python module guided by a structured specification document that captures the physical formulation, the software contract, and the acceptance tests in one place. Architectural decisions are persisted as ADRs (Architecture Decision Records), so a fresh agent session can re-acquire the design rationale without dialogue. The catalogue is validated quantitatively against canonical FEM benchmarks — Lamé, NAFEMS LE1, MacNeal-Harder, Bathe wave propagation, and Hill — and includes formulations from the author's own research that are not commonly available off-the-shelf in existing Python frameworks, in particular embedded strong discontinuities [@Retama:2010] and dissipation-based arc-length continuation [@Gutierrez:2004; @Verhoosel:2009].

# Statement of need

Computational solid mechanics research routinely produces formulations — return-mapping algorithms for new yield criteria, enriched element kinematics for fracture, custom path-following solvers — that need to be implemented and validated as part of the publication of the underlying theory. The implementation effort is large enough to deter many researchers from writing reusable code, and often produces software that is not preserved beyond the lifetime of the project. General-purpose finite element frameworks in Python — FEniCS [@Logg:2012], SfePy [@Cimrman:2019], scikit-fem [@Gustafsson:2020] — are well suited to defining problems through their weak form and solving them at scale, but their abstraction layer is not where the researcher's contribution typically lives. The primary artefacts in solid mechanics research are constitutive models and elements; weak-form abstractions, however powerful, are infrastructure layered above them.

Solidum FEM addresses this gap. It is a Python framework in which the primary unit of extension is an element-level Python class (with explicit shape functions, strain–displacement matrix, and integration scheme) or a material-level class (with explicit stress and tangent calculation). New components are added by writing a single self-registering module, with no changes to the core. Each component is accompanied by a structured specification document (`docs/specs/<Name>.md`) that fixes the physical formulation, the parameters, and the acceptance criteria in one place, and by an automated test against an analytical solution or a published benchmark.

The framework is also deliberately designed for collaboration with AI coding assistants. The combination of declarative contracts, self-registering modules, per-component specifications, and persisted Architecture Decision Records allows an AI agent — operating in independent sessions with no memory of past dialogue — to re-acquire enough context to implement and validate a new component. Solidum FEM has been developed primarily through AI-assisted programming, with the author defining formulations, reviewing diffs, and validating results against benchmarks; the architecture is the artefact that makes this division of labour viable.

The catalogue includes formulations relevant to the author's research that are not commonly available as off-the-shelf components in existing Python frameworks: embedded strong discontinuities for fracture, with a cohesive traction–jump model condensed locally onto a CST element [@Retama:2010]; dissipation-based arc-length continuation with automatic switching between cylindrical and dissipation constraints [@Gutierrez:2004; @Verhoosel:2009]; and time-domain integration via Hilber–Hughes–Taylor α [@Hilber:1977]. The catalogue is validated against canonical benchmarks of the FEM literature — Lamé closed-form solutions, NAFEMS LE1 elliptic membrane [@NAFEMS:1990], MacNeal–Harder slender beam [@MacNeal:1985], Bathe wave propagation [@Bathe:1996], and Hill's 1950 elastoplastic cylinder [@Hill:1950] — with 38 dedicated validation tests in `tests/validation/` and 804 tests in total covering both software contracts and physical correctness.

# Acknowledgements

This work was made possible by the institutional support of the Universidad Nacional Autónoma de México (UNAM) through the PAPIIT project IT101926 and the PAPIME project PE101626.

# Statement on the use of generative artificial intelligence

During the development of Solidum FEM, the author used Anthropic's Claude language models (Claude Sonnet and Claude Opus, accessed through the Claude Code command-line interface) as coding assistants for the implementation of elements, materials, solvers, tests, and documentation, following formulations, architectural decisions, and acceptance criteria defined by the author in the project's rules document, Architecture Decision Records, and per-component specifications. After each assisted change, the author reviewed the resulting diffs and validated the behaviour against analytical solutions and published benchmarks, taking full responsibility for the scientific rigor, mathematical accuracy, and physical integrity of the results.

# References
