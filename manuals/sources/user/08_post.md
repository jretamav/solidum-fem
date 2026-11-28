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
- **Todas las variables internas escalares** de los materiales, con su nombre (`damage_t`, `damage_c`, `kappa_t` del modelo de tracción y compresión; `cracked`, `w_max_n`, `damage_n`, `slip` de la grieta fija...), promediadas sobre los puntos de integración de la celda. `NaN` en las celdas cuyo material no tiene esa variable (en las figuras, gris). Si dos materiales usan el mismo nombre, cada celda lleva el valor del suyo, con el significado de su material. Los ángulos (la dirección $\theta$ de una grieta fija) no se promedian: la media de direcciones de puntos distintos, o de un punto sin grieta, no es una dirección; se leen en el archivo por puntos.
- `N_i`, `N_j`, `V_i`, `M_i`… (en 3D, `Vy`, `Vz`, `T`, `My`, `Mz`): **fuerzas internas de extremo** de barras, cables y marcos (`internal_forces`, ADR 0002), en ejes locales del elemento, en el nodo i y el nodo j, con la convención de signos del proyecto; `NaN` en las celdas que no tienen esa componente. Las que escribe `python -m solidum` descuentan la carga distribuida (peso propio, fuerza de cuerpo, cargas de línea) y son las mismas de la historia y del archivo de diagramas (hasta el 2026-11-27 el VTK daba las fuerzas nodales $K\,u$ con carga distribuida). Desde Python, `VtkExporter.export(..., element_forces=resultado.element_forces)` hace lo mismo. También las etapas estáticas de un análisis con alguna etapa en el tiempo. En los instantes en el tiempo de un análisis con carga distribuida valen `NaN`: no se conoce la carga del instante ni entra la inercia repartida, y la historia tampoco las registra.

- `Sigma_XX`, `Sigma_YY`, `Sigma_ZZ`, `Tau_XY`, `Tau_YZ`, `Tau_XZ`: las seis componentes del tensor de esfuerzos (en 2D las tres últimas son $\sigma_{zz}$ del material y ceros). Se evalúan en el campo `U` exportado vía `compute_gauss_state` y se promedian por elemento: el exportador no depende de qué solver produjo `U` ni de que el estado esté comiteado. Cubre sólidos 2D y 3D (`quad`, `triangle`, `quad8`, `quad9`, `triangle6`, `hexahedron`, `tetra`, `hexahedron20`, `hexahedron27`, `tetra10`), elementos térmicos (campo nodal `Temperature` y `Flux` por celda) y los elementos de los módulos de usuario, como el `CST_Embedded2D` del módulo `discontinuities`: el exportador deriva la celda VTK del número de nodos (o de `VTK_CELL_TYPE`, si el elemento lo declara), sin conocer la clase. Las barras y cables de dos nodos y los marcos son celdas `line`; las barras de tres nodos `Truss2D3` y `Truss3D3`, celdas `line3` (segmento cuadrático, con el nodo intermedio).

**Elementos de interfaz** (bloque `interfaces`): cada uno es una celda con los nodos de sus dos caras, `quad` (`Interface2D4`), celda 30 de VTK, cuadrilátero de 6 nodos (`Interface2D6`), `wedge` (`Interface3D6`) o `hexahedron` (`Interface3D8`). Con los nodos coincidentes la celda tiene volumen nulo; en ParaView se ve al deformar con *Warp By Vector* o al separar las caras. Sus campos por celda son la media de sus puntos de integración:

- `Jump_n`, `Jump_s` (en 3D `Jump_s1`, `Jump_s2`): desplazamiento relativo en ejes locales; `Jump_n` es la holgura actual (en el pre-proceso, la holgura inicial).
- `Traction_n`, `Traction_s` (`Traction_s1`, `Traction_s2`): tracción sobre la cara A en ejes locales; en contacto, `Traction_n` = $-p_N$.
- Las variables internas escalares de la ley (`contact_pressure`, `slip_accumulated`, `damage`…), salvo los códigos: el régimen `status` de `ContactPenalty` (0 abierta, 1 adherida, 2 deslizando) no se promedia y se lee en el archivo por puntos.

Son componentes escalares, no un vector, porque un vector en ParaView se lee en ejes globales. En las celdas de interfaz los esfuerzos de sólido (`Sigma_*`, `Von_Mises`) valen `NaN` y en las de sólido los `Jump_*` y `Traction_*` también; las interfaces no entran en el suavizado nodal.

Las interfaces de adherencia que genera el refuerzo embebido con `bond` (`BondSlipInterface`) no tienen celda propia y no se exportan al VTK: unen los nodos de la barra con los de su anfitrión. Su estado (deslizamiento, esfuerzo de adherencia) se registra en la historia (`var: jump` o `traction` con `component: s`, y sus variables) o se consulta con `EmbeddedBar.bond_state` (capítulo «Uso Avanzado: API Programática»).

## Diagramas de fuerzas internas

Con `output.results.diagrams`, Solidum escribe las fuerzas internas **a lo largo** de cada barra y cada marco (ADR 0044):

```yaml
output:
  results:
    diagrams: {stations: 11, at: [final]}    # o diagrams: true
```

- `<nombre>_diagramas_final.csv` (y `<nombre>_diagramas_NNNN.csv` por cada factor de `at`): una fila por estación de cada elemento con `element`, `s` (distancia desde el nodo i), `xi` (fracción de la longitud), `x`, `y`, `z` (posición sobre la geometría sin deformar) y las fuerzas internas (`N`, `V`, `M` en 2D; `N`, `Vy`, `Vz`, `T`, `My`, `Mz` en 3D; vacías donde no aplican).
- Estaciones: `stations` uniformes por elemento (11 por omisión) más los puntos donde se anula cada cortante dentro del elemento, donde el flector es extremo: el máximo del claro sale exacto aunque no caiga en una estación.
- Se calculan por equilibrio del tramo, a partir de las fuerzas de extremo y de la carga a lo largo del elemento: son exactas con cualquier material una vez convergido el paso, también dentro de un elemento cargado (la interpolación del elemento, lineal en el flector, no ve la carga).
- Sólo en análisis estáticos (con etapas, si ninguna es en el tiempo): en un dinámico faltarían la inercia y el amortiguamiento repartidos sobre el elemento. No en las armaduras cuadráticas.
- En los elementos corotacionales, el equilibrio se plantea sobre la cuerda deformada: no incluye el momento de la fuerza axial por la separación entre el eje deformado y la cuerda, que se recupera con varios elementos por barra.

Desde Python, `resultado.internal_forces_along(elemento, xi)` da lo mismo en los puntos `xi` (fracción de la longitud desde el nodo i): un `ElementForcesAlong` con `xi`, `s` y `components`.

## Valores en los puntos de integración

Con `output.results.integration_points: true`, Solidum escribe en los mismos instantes que el VTK por celda `<nombre>_puntos_NNNN.vtu`, `<nombre>_puntos_final.vtu` y su índice `<nombre>_puntos.pvd`: una nube de puntos (celdas `vertex`) en la posición de cada punto de integración de los sólidos y de las interfaces, sin promediar. Por punto:

- `Displacements`: el desplazamiento interpolado $\mathbf u(\boldsymbol\xi_g) = \sum_a N_a(\boldsymbol\xi_g)\,\mathbf u_a$ (para deformar la nube con *Warp By Vector*);
- `Element_Id` y `Point` (número del punto en su elemento, desde 1);
- `Sigma_XX` … `Tau_XZ` y `Von_Mises`;
- todas las variables internas escalares del material, también los ángulos; `NaN` donde el material no la tiene;
- en los puntos de los elementos de interfaz (sobre la cara A): `Jump_*`, `Traction_*` y todas las variables escalares de la ley, también `status`; los esfuerzos de sólido valen `NaN`.

Los valores son los del estado con que el análisis equilibró el paso, los mismos que registra la historia, sin volver a evaluar el material. Sirve para leer lo que el promedio por celda esconde: el ángulo de cada grieta fija, qué puntos de un elemento se agrietaron, el daño punto a punto.

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
