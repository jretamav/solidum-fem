# Post-Procesamiento y Visualización

## Formato VTU (VTK XML)

Mediante `meshio`, Solidum FEM exporta los resultados en archivos `.vtu` compatibles con ParaView. La información se clasifica en dos diccionarios estándar:

### Datos Nodales (Point Data)

- `Displacements`: campo vectorial $(u_x, u_y, u_z)$ por nodo.
- `Rotations`: giros $(r_x, r_y, r_z)$, sólo si algún nodo tiene grados de libertad de giro (marcos).
- `External_Forces`: distribución nodal del vector de cargas aplicado en el paso actual.
- `Supports` (y `Supports_Rot` con giros): indicador de DOFs restringidos, 1.0 en cada componente con Dirichlet.
- `Temperature`: temperatura nodal, sólo si algún nodo tiene el DOF `T`.
- `Nodal_Mass`, `Nodal_Inertia`: masas de nudo $(m_x, m_y, m_z)$ e inercias rotacionales $(J_x, J_y, J_z)$ declaradas (`nodal_masses`, ADR 0024), sólo si el modelo tiene alguna; `Nodal_Inertia`, si tiene giros. Sirven para revisar la entrada de masas en ParaView.
- `Sigma_XX_nodal` … `Tau_XZ_nodal`, `Von_Mises_nodal`: suavizado nodal de orden 0 de los esfuerzos por celda (media simple sobre las celdas adyacentes a cada nodo). Sólo promedian las celdas de **sólidos mecánicos**: las barras, cables y marcos (celdas `line`) y los elementos térmicos no participan, tampoco cuando un elemento térmico se superpone a un sólido sobre los mismos nodos en un dominio multicampo $(u_x, u_y, T)$. Un modelo sin sólidos mecánicos (sólo térmico o sólo estructural 1D) no exporta estos campos.

#### Restricciones multipunto y suavizado nodal

Con restricciones multipunto (bloque `linear_constraints` en YAML, `add_linear_constraint` en Python) el archivo recibe la solución completa, con los DOF esclavos ya recuperados: `Displacements`, `Temperature` y todos los campos por celda son los del modelo resuelto también en los nodos esclavos. El suavizado nodal tiene una limitación: dos nodos coincidentes atados por una restricción ($u_s = u_m$, p. ej. para unir mallas no conformes) siguen siendo dos nodos distintos, y cada uno promedia sólo las celdas de su lado. En la interfaz, `Sigma_*_nodal` del nodo maestro y del esclavo difieren entre sí y del valor que daría la malla con el nodo fundido (que, si cada lado aporta el mismo número de celdas, es la media de las dos copias). Para leer el esfuerzo en la interfaz, use los campos por celda o promedie las dos copias.

### Datos de Elemento (Cell Data)

- `Von_Mises`: esfuerzo equivalente J2, calculado como invariante 3D con las seis componentes; en plane strain incluye $\sigma_{zz}$ (que aporta el material vía `out_of_plane_stress`), en plane stress $\sigma_{zz} = 0$.
- `Internal_State`: variable principal del material según su `PRIMARY_STATE_VAR` (acumulada equivalente $\alpha$ en plasticidad, $d$ en daño, etc.).
- `N_i`, `N_j`, `V_i`, `M_i`… (en 3D, `Vy`, `Vz`, `T`, `My`, `Mz`): **fuerzas internas de extremo** de barras, cables y marcos (`internal_forces`, ADR 0002), en ejes locales del elemento, en el nodo i y el nodo j, con la convención de signos del proyecto; `NaN` en las celdas que no tienen esa componente. Con carga distribuida, `VtkExporter.export(..., element_forces=resultado.element_forces)` descuenta la carga equivalente.
- `Sigma_XX`, `Sigma_YY`, `Sigma_ZZ`, `Tau_XY`, `Tau_YZ`, `Tau_XZ`: las seis componentes del tensor de esfuerzos (en 2D las tres últimas son $\sigma_{zz}$ del material y ceros). Se evalúan en el campo `U` exportado vía `compute_gauss_state` y se promedian por elemento: el exportador no depende de qué solver produjo `U` ni de que el estado esté comiteado. Cubre sólidos 2D y 3D (`quad`, `triangle`, `quad8`, `quad9`, `triangle6`, `hexahedron`, `tetra`, `hexahedron20`, `hexahedron27`, `tetra10`), elementos térmicos (campo nodal `Temperature` y `Flux` por celda) y los elementos de los módulos de usuario, como el `CST_Embedded2D` del módulo `discontinuities`: el exportador deriva la celda VTK del número de nodos (o de `VTK_CELL_TYPE`, si el elemento lo declara), sin conocer la clase. Las barras y cables de dos nodos y los marcos son celdas `line`; las barras de tres nodos `Truss2D3` y `Truss3D3`, celdas `line3` (segmento cuadrático, con el nodo intermedio).

Las interfaces de adherencia que genera el refuerzo embebido con `bond` (`BondSlipInterface`) no tienen celda propia y no se exportan: unen los nodos de la barra con los de su anfitrión. Su estado (deslizamiento, esfuerzo de adherencia) se consulta con `EmbeddedBar.bond_state` (capítulo «Uso Avanzado: API Programática»).

## Animación con Archivo PVD

Cuando se exporta con `frequency: "all_steps"`, además de un `.vtu` por paso se genera un archivo maestro `.pvd` (XML colección) que ParaView abre como animación temporal. El "tiempo" en el archivo PVD coincide con el factor de carga $\lambda$.

## Salida en Texto Plano (Estilo FEAP)

Opcional para verificación manual o post-procesamiento con scripts. Configurable por nodos y elementos seleccionados:

```yaml
output:
  text_export:
    export: true
    nodes: [1, 5, 10]                # o "all"
    elements: all
    nodal_results: ["Displacements", "External_Forces"]
    element_results: ["Von_Mises", "Internal_State"]
```

## Workflow ParaView

1. `File > Open` y seleccionar el archivo `.pvd` (o el grupo de `.vtu`).
2. *Apply* en el panel *Properties*.
3. Reproducir la animación con los controles temporales de la barra superior.
4. Aplicar filtro **Warp By Vector** con `Displacements` para visualizar la deformación física (ajustar `Scale Factor` para amplificar).
5. Cambiar el coloreado de `Solid Color` a la variable de interés (`Von_Mises`, `Internal_State`, `Sigma_XX`, …).
