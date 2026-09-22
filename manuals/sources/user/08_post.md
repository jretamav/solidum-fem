# Post-Procesamiento y Visualización

## Formato VTU (VTK XML)

Mediante `meshio`, Solidum FEM exporta los resultados en archivos `.vtu` compatibles con ParaView. La información se clasifica en dos diccionarios estándar:

### Datos Nodales (Point Data)

- `Displacements`: campo vectorial $(u_x, u_y, u_z)$ por nodo.
- `External_Forces`: distribución nodal del vector de cargas aplicado en el paso actual.
- `Supports`: indicador booleano de DOFs restringidos.

### Datos de Elemento (Cell Data)

- `Von_Mises`: esfuerzo equivalente J2, calculado como invariante 3D con las seis componentes; en plane strain incluye $\sigma_{zz}$ (que aporta el material vía `out_of_plane_stress`), en plane stress $\sigma_{zz} = 0$.
- `Internal_State`: variable principal del material según su `PRIMARY_STATE_VAR` (acumulada equivalente $\alpha$ en plasticidad, $d$ en daño, etc.).
- `Sigma_XX`, `Sigma_YY`, `Sigma_ZZ`, `Tau_XY`, `Tau_YZ`, `Tau_XZ`: las seis componentes del tensor de esfuerzos (en 2D las tres últimas son $\sigma_{zz}$ del material y ceros). Se evalúan en el campo `U` exportado vía `compute_gauss_state` y se promedian por elemento: el exportador no depende de qué solver produjo `U` ni de que el estado esté comiteado. Cubre sólidos 2D y 3D (`quad`, `triangle`, `quad8`, `quad9`, `triangle6`, `hexahedron`, `tetra`, `hexahedron20`, `hexahedron27`, `tetra10`), elementos térmicos (campo nodal `Temperature` y `Flux` por celda) y `CST_Embedded2D`.

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
