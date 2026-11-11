# Catálogo de restricciones generadas

> Referencia rápida de los generadores de restricciones ([ADR 0032](adr/0032-restricciones-generadas-y-refuerzo-embebido.md)). Un generador, aplicado a un dominio ya mallado, crea a partir de su geometría los nodos, los elementos y las restricciones afines (MPC, [ADR 0004](adr/0004-imposicion-de-condiciones-de-dirichlet.md)) que necesita. En Python se aplica con su función; en el YAML se declara en el bloque `constraints` (lista de `{type, …}`), que se aplica después de crear los elementos. Las MPC escritas una a una siguen en el bloque `linear_constraints`.

---

## EmbeddedReinforcement — refuerzo embebido con adherencia perfecta o con adherencia-deslizamiento

Barra de refuerzo recta o poligonal dentro de sólidos 2D o 3D. La barra es una cadena de elementos 1D ordinarios (`Truss2D`, `Truss3D`, sus variantes corrotacionales, cables, y las de tres nodos `Truss2D3` y `Truss3D3`), con cualquier material 1D. Su malla no depende de la del sólido. Spec: [EmbeddedReinforcement](specs/EmbeddedReinforcement.md).

### Formulación
- **Adherencia perfecta** (por omisión): cada nodo de la barra sigue el desplazamiento interpolado del elemento de sólido que lo contiene: `u(X_b) = Σ N_i(ξ_b)·u_i`, una MPC por traslación. Es exacta también con desplazamientos grandes, porque `ξ_b` es un punto material.
- **Adherencia-deslizamiento** (`bond`, spec §12): sin MPC; cada tramo de la barra se une a su anfitrión con una interfaz distribuida (`BondSlipInterface`, generada) cuyo desplazamiento relativo `Δu = u_barra − u_sólido` se descompone en deslizamiento `s` a lo largo de la barra y componentes normales; la ley de interfaz de línea da `τ(s)` sobre el perímetro y `t_n = k_n·Δu_n` sobre el diámetro.
- La polilínea se corta en los bordes de los anfitriones. Con anfitriones de interpolación lineal a lo largo de la recta (Tri3, Tet4; Quad4 y Hex8 sin distorsión con la barra paralela a una arista), la deformación de la barra es `tᵀ ε t`. En los demás casos es la deformación media del tramo, y converge al refinar la malla.
- El área de concreto que ocupa la barra **no se descuenta** (de Borst et al. 2012, ec. 6.136; como Abaqus y DIANA). El sobrante de rigidez elástica es `ρ/((1−ρ)+nρ)`, con `n = E_s/E_c`: +0.9/+1.8/+3.3 % con ρ = 1/2/4 %. Con el concreto agrietado en tracción es cero.

### Uso
```python
bar = solidum.embed_reinforcement(domain, [[0.0, 0.05], [4.0, 0.05]], area=1.6e-3,
                                  material=solidum.Elastic1D(E=200e9))
r = solidum.run(domain, F_applied=F)
bar.bond_forces(r)      # fuerza de adherencia en cada nodo de la barra
r.element_forces[bar.element_ids[0]].components["N"]
```
```yaml
constraints:
  - type: EmbeddedReinforcement
    points: [[0.0, 0.05], [4.0, 0.05]]
    area: 1.6e-3
    material: 2          # id de un material 1D de 'materials'
    element: Truss2D     # opcional (Truss3D por omisión en 3D)
    hosts: [1, 2, 3]     # opcional: ids de los anfitriones
    subdivisions: 1      # opcional: partes por tramo
```

### Parámetros
- `points`: polilínea de la barra, con todos sus puntos dentro de los anfitriones.
- `area`: área `A_s`.
- `material`: material 1D.
- `element`: elemento 1D de dos nodos o de tres (con el intermedio de cada tramo en su punto medio, atado como los extremos).
- `hosts`: anfitriones admitidos.
- `subdivisions`: partes por tramo.
- `element_kwargs`: argumentos extra del elemento 1D (p. ej. `L0` de un cable).
- `bond`: adherencia-deslizamiento (spec §12): `{law: id en interface_materials, diameter: d_b, perimeter: (π·d_b), integration: nodal|gauss}`. Cada tramo se une a su anfitrión con una interfaz distribuida (`BondSlipInterface`, generada, una por tramo) en lugar de atarse. Barra recomendada: `Truss2D3`/`Truss3D3` con integración nodal (deslizamiento de orden 4 frente al orden 2 de la lineal, sin oscilación). Leyes: `BondSlipModelCode`, `BondSlipMultilinear` ([catálogo de materiales](catalogo_materiales.md), «Materiales de interfaz»).

### Resultados
- Los nodos y elementos de la barra son normales: desplazamientos, fuerza axial `N` en `element_forces` y exportación VTK.
- Con `bond`: `EmbeddedBar.bond_state(domain, result.U)` da, en cada punto de integración de cada tramo, posición, deslizamiento `s`, `τ`, `Δu_n`, `t_n` y la fuerza de adherencia por unidad de longitud `p·τ`. Las interfaces no se exportan a VTK (no tienen celda propia).
- Con adherencia perfecta, las fuerzas de las MPC de embebido salen en `constraint_forces`, no en las reacciones: son las fuerzas de adherencia concentradas, `EmbeddedBar.bond_forces(result)` (con `bond` no hay MPC y se usa `bond_state`).
- Los nodos de la barra quedan en `Domain.generated_nodes`. Los selectores `*_by_coord` del YAML no los toman.

### Régimen de validez
- Adherencia perfecta, o adherencia-deslizamiento distribuida con una ley de interfaz de línea (`bond`).
- La barra sigue el campo regular del anfitrión. En un anfitrión con discontinuidad embebida ve la deformación media del elemento, apertura de grieta incluida.

### Fuera de alcance
- Degradación cíclica de la adherencia más allá de la descarga y recarga de las leyes del catálogo (regla R1 de `BondSlipModelCode`, sin degradación).
- Barras con rigidez a flexión.
- Barras curvas (se aproximan con la polilínea).
- Nodos fuera del sólido.
