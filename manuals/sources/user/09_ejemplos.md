# Ejemplos Completos

## Marco 2D — Viga en Voladizo

Combinación de un tramo Euler-Bernoulli y un tramo Timoshenko, empotrada en el extremo y cargada en la punta.

```yaml
nodes:
  - {id: 1, coords: [0.0, 0.0]}
  - {id: 2, coords: [2.0, 0.0]}
  - {id: 3, coords: [4.0, 0.0]}

materials:
  - id: 1
    type: Elastic1D
    E: 200.0e9
    density: 7850.0

elements:
  - id: 1
    type: Frame2DEuler
    material: 1
    nodes: [1, 2]
    A: 0.01
    I: 8.33e-6

  - id: 2
    type: Frame2DTimoshenko
    material: 1
    nodes: [2, 3]
    A: 0.01
    I: 8.33e-6
    As: 0.00833

boundary_conditions_by_node:
  - {node_id: 1, ux: 0.0, uy: 0.0, rz: 0.0}

point_loads_by_node:
  - {node_id: 3, uy: -15000.0}

solver:
  type: LinearSolver

output:
  file_name: "viga_voladizo_marco"
  pre_process: {export: true}
  results:
    frequency: "last_step"
    nodal_results: ["Displacements"]
```

## Placa Continua con Plasticidad J2 (Gmsh)

Placa con perforación central importada desde `.msh`, sometida a carga axial creciente hasta superar el límite elástico del acero. Solver no lineal con paso adaptativo.

```yaml
mesh: "placa_agujero.msh"

mesh_physical_groups:
  "Superficie_Placa":
    material: 1
    thickness: 0.05
    quadrature: "2x2"

materials:
  - id: 1
    type: VonMises2D
    E: 200.0e9
    nu: 0.3
    sigma_y: 250.0e6
    H: 1.5e9
    hypothesis: plane_strain
    density: 7850.0

boundary_conditions_by_group:
  "Borde_Fijo": {ux: 0.0, uy: 0.0}

point_loads_by_group:
  "Borde_Carga": {ux: 1500000.0}

solver:
  type: NonlinearSolver
  num_steps: 30
  max_iter: 10
  adaptive: true
  convergence:
    rtol_force: 1.0e-5
    rtol_disp: 1.0e-5

output:
  file_name: "estudio_placa"
  pre_process: {export: true}
  results:
    frequency: "all_steps"
    nodal_results: ["Displacements"]
    element_results: ["Von_Mises", "Internal_State"]
```

**Interpretación**: durante los primeros pasos (régimen elástico) Newton-Raphson converge en 1–2 iteraciones. Cuando la concentración de esfuerzos en la periferia del orificio supera $\sigma_y = 250$ MPa, comienzan iteraciones del *return mapping*; `Internal_State` (la deformación plástica acumulada $\alpha$) crece visiblemente alrededor del agujero.

## Cable Pretensado en Régimen Corotacional

Cable 2D con material unilateral, cargado transversalmente en su nodo central. Demuestra el acoplamiento elemento corotacional + material unilateral.

```yaml
nodes:
  - {id: 1, coords: [0.0, 0.0]}
  - {id: 2, coords: [5.0, 0.0]}
  - {id: 3, coords: [10.0, 0.0]}

materials:
  - id: 1
    type: CableMaterial1D
    E: 150.0e9
    density: 7850.0

elements:
  - id: 1
    type: Cable2DCorot
    material: 1
    nodes: [1, 2]
    A: 5.0e-5

  - id: 2
    type: Cable2DCorot
    material: 1
    nodes: [2, 3]
    A: 5.0e-5

boundary_conditions_by_node:
  - {node_id: 1, ux: 0.0, uy: 0.0}
  - {node_id: 3, ux: 0.0, uy: 0.0}

point_loads_by_node:
  - {node_id: 2, uy: -500.0}

solver:
  type: NonlinearSolver
  num_steps: 50
  max_iter: 25
  adaptive: true
  convergence:
    rtol_force: 1.0e-5
    rtol_disp: 1.0e-5

output:
  file_name: "cable_carga_transversal"
  pre_process: {export: true}
  results:
    frequency: "all_steps"
    nodal_results: ["Displacements"]
```
