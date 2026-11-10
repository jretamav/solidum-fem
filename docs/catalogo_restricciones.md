# Catálogo de restricciones generadas

> Referencia rápida de los generadores de restricciones ([ADR 0032](adr/0032-restricciones-generadas-y-refuerzo-embebido.md)). Un generador, aplicado a un dominio ya mallado, crea a partir de su geometría los nodos, los elementos y las restricciones afines (MPC, [ADR 0004](adr/0004-imposicion-de-condiciones-de-dirichlet.md)) que necesita. En Python se aplica con su función; en el YAML se declara en el bloque `constraints` (lista de `{type, …}`), que se aplica después de crear los elementos. Las MPC escritas una a una siguen en el bloque `linear_constraints`.

---

## EmbeddedReinforcement — refuerzo embebido con adherencia perfecta

Barra de refuerzo recta o poligonal dentro de sólidos 2D o 3D. La barra es una cadena de elementos 1D ordinarios (`Truss2D`, `Truss3D`, sus variantes corrotacionales, cables), con cualquier material 1D. Su malla no depende de la del sólido. Spec: [EmbeddedReinforcement](specs/EmbeddedReinforcement.md).

### Formulación
- Cada nodo de la barra sigue el desplazamiento interpolado del elemento de sólido que lo contiene: `u(X_b) = Σ N_i(ξ_b)·u_i`, una MPC por traslación. Es exacta también con desplazamientos grandes, porque `ξ_b` es un punto material.
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
- `element`: elemento 1D de dos nodos.
- `hosts`: anfitriones admitidos.
- `subdivisions`: partes por tramo.
- `element_kwargs`: argumentos extra del elemento 1D (p. ej. `L0` de un cable).

### Resultados
- Los nodos y elementos de la barra son normales: desplazamientos, fuerza axial `N` en `element_forces` y exportación VTK.
- Las fuerzas de las MPC de embebido salen en `constraint_forces`, no en las reacciones: son las fuerzas de adherencia concentradas, `EmbeddedBar.bond_forces(result)`.
- Los nodos de la barra quedan en `Domain.generated_nodes`. Los selectores `*_by_coord` del YAML no los toman.

### Régimen de validez
- Adherencia perfecta.
- La barra sigue el campo regular del anfitrión. En un anfitrión con discontinuidad embebida ve la deformación media del elemento, apertura de grieta incluida.

### Fuera de alcance
- Deslizamiento barra-sólido.
- Barras con rigidez a flexión.
- Barras curvas (se aproximan con la polilínea).
- Nodos fuera del sólido.
