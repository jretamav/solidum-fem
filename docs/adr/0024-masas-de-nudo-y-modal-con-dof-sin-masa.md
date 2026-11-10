# ADR 0024 — Masas concentradas en nudos y análisis modal con grados de libertad sin masa

- **Estado**: aceptado (2026-11-09); masas sólo inerciales, decisión del usuario del mismo día
- **Fecha**: 2026-11-09
- **Alcance**: `Domain.add_nodal_mass`, YAML (`nodal_masses`), `Assembler.assemble_mass_matrix`, `EigenSolver` (condensación estática de los grados de libertad sin masa), `ModalSolver` (error tipado `MasslessDofsError`) y, por herencia, `ResponseSpectrumSolver`. Cumple lo que el ADR 0008 dejó para «un ADR propio»: masas concentradas que no salen de la densidad del material.

## Contexto

Hasta ahora la única fuente de masa era la densidad del material (ADR 0008). En un edificio, la masa sísmica de cada nivel viene sobre todo de la losa, los acabados y una fracción de la carga viva, no del peso propio de vigas y columnas, y la práctica estándar la concentra en los nudos de cada nivel: `mass` en OpenSees, masa de nudo en SAP2000, el elemento `MASS21` de ANSYS, `*MASS` en Abaqus. Sin ella no había análisis modal ni espectral de un edificio con masas de piso, ni patrones de carga lateral proporcionales a la masa.

Con masas sólo en algunos nudos, casi todos los grados de libertad quedan **sin masa** y la matriz de masa es semidefinida. Pasa lo mismo con elementos de `density = 0.0`, que el ADR 0008 admite. El `ModalSolver` no lo resolvía:

- con más de 32 grados de libertad libres, ARPACK lanzaba `ArpackError -9999` («Could not build an Arnoldi factorization»), que no dice la causa;
- con 32 o menos, el camino denso (`scipy.linalg.eigh`) lanzaba `LinAlgError`, porque exige `M` definida.

Medido en un marco de un nivel con columnas de densidad cero y masa sólo en la viga: 18 de sus 33 grados de libertad libres no tienen masa, y el problema sí tiene solución (15 frecuencias finitas, calculadas con `numpy` como valores propios de `K⁻¹M`).

## Decisión

### 1. La masa de nudo es un dato del dominio

```python
domain.add_nodal_mass(5, ux=12.0e3, uy=12.0e3)      # masa [kg]
domain.add_nodal_mass(5, rz=4.0e3)                  # inercia rotacional [kg·m²]
```

- Se declara por nodo y por **nombre de grado de libertad mecánico**: traslacional (`ux`, `uy`, `uz`) con masa, o rotacional (`rx`, `ry`, `rz`) con inercia rotacional (masa·longitud²). Las rotacionales son una opción, por decisión del usuario (2026-11-09): para un empuje lateral bastan las traslacionales, pero la dinámica de marcos con masas de piso las usa a veces.
- Varias declaraciones sobre el mismo grado de libertad **se suman**, como las cargas nodales.
- Cada valor debe ser finito y no negativo, y el nodo debe existir. Que el nodo tenga ese grado de libertad se comprueba al ensamblar, porque los grados de libertad los crean los elementos.
- Un grado de libertad no mecánico (`T`) se rechaza: su «masa» sería una capacidad calorífica concentrada, otra magnitud con otras unidades.

**Por qué en el dominio y no como elemento** (`MASS21`, `*MASS`): un elemento necesitaría material, una rigidez nula, cinemática, sitio en el ensamblaje por lotes y entrada en los catálogos, y sólo aporta una entrada en la diagonal de `M`. Las cargas nodales ya viven como datos por (nodo, grado de libertad), y OpenSees y SAP2000 tratan la masa de nudo igual. **Por qué no en el material**: no es una propiedad del medio continuo (ADR 0008), sino una masa no estructural apoyada en el nudo.

### 2. Ensamblaje

`Assembler.assemble_mass_matrix` suma cada masa de nudo a la entrada diagonal de su grado de libertad, **igual con `lumping="consistent"` que con `"lumped"`**: una masa puntual es la misma en los dos esquemas. La suma se hace sobre el patrón de dispersión ya cacheado, así que `M` conserva el patrón de `K`. La caché de `M` incluye la huella de las masas de nudo, de modo que una masa añadida después del primer ensamblaje no se ignora.

Consecuencias:
- Todo consumidor de `M` (modal, espectral, armónico, Newmark, HHT, diferencias centrales) ve las masas de nudo sin cambios propios.
- La regla del ADR 0008 no cambia: un material sin `density` sigue siendo error. Un elemento que no aporta masa se declara con `density = 0.0`.

### 3. Sólo inercia

La masa de nudo entra en `M`. **No genera peso con `gravity`**, igual que el `mass` de OpenSees. Si su peso actúa, se declara como carga.

**Decisión del usuario (2026-11-09): `gravity` no genera el peso de las masas de nudo.** Razones:
- en un análisis de empuje, la masa sísmica y la carga de gravedad salen de combinaciones distintas (ASCE 41: gravedad 1.1·(D + 0.25L); masa de la carga muerta y una fracción de la viva), y acoplarlas impondría una igualdad que las normas no piden;
- coincide con el `mass` de OpenSees;
- evita contar dos veces el peso de quien declara también las cargas de piso.

La carga de gravedad de los pisos se declara como carga (`point_loads`, o carga constante en un análisis por etapas).

### 4. YAML

```yaml
nodal_masses:
  - {node_id: 5, ux: 12.0e3, uy: 12.0e3}
  - {node_id: 9, ux: 12.0e3, uy: 12.0e3, rz: 4.0e3}
```

El bloque se valida al leer, con todos los errores juntos: forma, nombres de grados de libertad, valores y existencia de nodos y grados de libertad. Las variantes por grupo físico de gmsh y por coordenada (`nodal_masses_by_group` y `nodal_masses_by_coord`, como `point_loads_by_group` y `point_loads_by_coord`) se añadieron el 2026-11-09: cada nodo seleccionado recibe la masa indicada. No cambian nada de lo decidido aquí.

### 5. Análisis modal con grados de libertad sin masa: condensación estática exacta

Un grado de libertad libre no tiene masa si su entrada en la diagonal de `M_red` es cero exacto. En una matriz semidefinida positiva, una diagonal nula implica fila y columna nulas. Con `m` los grados de libertad con masa y `0` los que no la tienen:

$$\begin{bmatrix}\mathbf K_{mm} & \mathbf K_{m0}\\ \mathbf K_{0m} & \mathbf K_{00}\end{bmatrix}\begin{bmatrix}\boldsymbol\phi_m\\ \boldsymbol\phi_0\end{bmatrix} = \omega^2\begin{bmatrix}\mathbf M_{mm} & \mathbf 0\\ \mathbf 0 & \mathbf 0\end{bmatrix}\begin{bmatrix}\boldsymbol\phi_m\\ \boldsymbol\phi_0\end{bmatrix}$$

La segunda fila no tiene inercia: $\boldsymbol\phi_0 = -\mathbf K_{00}^{-1}\mathbf K_{0m}\,\boldsymbol\phi_m$. No es una aproximación, a diferencia de la reducción de Guyan con grados de libertad que sí tienen masa. Sustituyendo:

$$\hat{\mathbf K}\,\boldsymbol\phi_m = \omega^2\,\mathbf M_{mm}\,\boldsymbol\phi_m,\qquad \hat{\mathbf K} = \mathbf K_{mm} - \mathbf K_{m0}\,\mathbf K_{00}^{-1}\,\mathbf K_{0m}.$$

Sus valores propios son exactamente las frecuencias finitas del problema original; las otras $n_0$ son infinitas.

**Implementación**, sin formar $\hat{\mathbf K}$ en los modelos grandes:
- Con $\mathbf A = \mathbf K - \sigma\mathbf M$, el bloque $mm$ de $\mathbf A^{-1}$ es $(\hat{\mathbf K} - \sigma\mathbf M_{mm})^{-1}$: es el complemento de Schur, porque $\mathbf M$ no tiene filas $0$.
- Así, el operador de desplazamiento e inversión del problema condensado se aplica con la factorización dispersa de $\mathbf A$ del modelo completo, y ARPACK trabaja con el producto interno de $\mathbf M_{mm}$, que es definida.
- $\boldsymbol\phi_0$ se recupera con la factorización de $\mathbf K_{00}$.
- Con 32 o menos grados de libertad con masa, o pidiendo casi todos los modos, se forma $\hat{\mathbf K}$ densa y se usa `eigh`.

Los modos expandidos siguen siendo $\mathbf M$-ortonormales. Un modelo con masa en todos sus grados de libertad sigue exactamente el camino de antes.

El reintento con desplazamiento negativo de las estructuras libres (con `K` singular) se conserva en el camino condensado.

### 6. Error tipado cuando el problema no tiene solución

`MasslessDofsError` (`ValueError`, en `solidum.math.solvers.diagnostics`) se lanza en estos casos, y su mensaje nombra los grados de libertad sin masa como (nodo, dof):
- se piden más modos que grados de libertad con masa (sólo esos tienen frecuencia finita);
- ningún grado de libertad libre tiene masa;
- los grados de libertad sin masa forman un mecanismo ($\mathbf K_{00}$ singular): se mueven sin rigidez ni inercia, y su movimiento (en el modal, su frecuencia) está indeterminado; también en los transitorios implícitos (§6 bis);
- la matriz de masa de los grados de libertad con masa sigue siendo singular.

### 6 bis. Transitorios implícitos: aceleración inicial consistente (decisión del usuario, 2026-11-09)

La integración implícita no necesita $\mathbf M$ definida: $\mathbf M + \gamma\Delta t\,\mathbf C + \beta\Delta t^2\mathbf K$ es regular aunque haya grados de libertad sin masa. Sólo la aceleración inicial, $\mathbf M\ddot{\mathbf u}_0 = \mathbf r_0 = \mathbf F(t_0) - \mathbf F_\text{int}(\mathbf u_0) - \mathbf C\dot{\mathbf u}_0$, queda indeterminada en ellos. Newmark, HHT y sus variantes de Newton la toman **consistente**:

$$\ddot{\mathbf u}_m = \mathbf M_{mm}^{-1}\,\mathbf r_{0,m}, \qquad \ddot{\mathbf u}_0 = -\mathbf K_{00}^{-1}\mathbf K_{0m}\,\ddot{\mathbf u}_m ,$$

con $\mathbf K$ la rigidez (tangente) en $\mathbf u_0$. Las filas sin masa no tienen inercia: $\mathbf K_0\mathbf u = \mathbf F_0$ en todo instante (con Rayleigh, $\mathbf C_0 = a_1\mathbf K_0$ en esas filas, y la relación es la misma), y derivarla da esa aceleración. Los grados de libertad sin masa siguen de forma cuasiestática a los que la tienen, como en el análisis modal (§5).

Alternativa descartada: $\ddot{\mathbf u}_0 = \mathbf 0$ en los grados de libertad sin masa (lo que hace OpenSees al arrancar con aceleraciones nulas). Medido en el edificio de cortante de tres pisos (3 grados de libertad con masa y 18 sin masa, Newmark con $\Delta t = T_1/40$): los desplazamientos y las velocidades de todos los nudos son los mismos con las dos opciones (diferencia ≤ 5e-14), porque la ecuación algebraica de las filas sin masa fija $\mathbf u_0$ en cada paso; pero con una carga aplicada desde $t_0$ la aceleración reportada en los grados de libertad sin masa arrastra con $\ddot{\mathbf u}_0 = \mathbf 0$ un error del 50 % de la aceleración de piso que Newmark no amortigua, y con la consistente, de 3e-11. Partiendo del reposo con $\mathbf F(t_0) = \mathbf 0$ las dos coinciden.

Si el estado inicial no equilibra las filas sin masa ($\mathbf r_{0,0} \neq \mathbf 0$: una carga aplicada de golpe en un nudo sin masa), se avisa: esos grados de libertad se ajustan en el primer paso. Un mecanismo sin masa ($\mathbf K_{00}$ singular) da `MasslessDofsError` con causa `massless_mechanism`. Pruebas: `tests/test_transient_massless_dofs.py` (el edificio contra el problema 3×3 condensado integrado con el mismo esquema, Newmark y HHT-α, sin y con Rayleigh, a 1e-9) y `tests/test_transient_termination.py`.

### 7. Fuera de este ADR

- ~~Transitorio con grados de libertad sin masa~~: resuelto en el §6 bis (2026-11-09).
- Diferencias centrales exige masa en todos los grados de libertad, porque es un método explícito; eso no cambia.

## Alternativas consideradas

- **Elemento de masa puntual** (`MASS21`, `*MASS`): rechazado (§1).
- **Repartir la masa del nivel como densidad equivalente de las vigas**: distorsiona la distribución de la masa y su inercia rotacional. No es la masa que se quiere modelar.
- **Masa ficticia pequeña en los grados de libertad sin masa**, como hacen algunos programas: introduce frecuencias espurias y cambia los resultados en función de un número sin significado físico. Descartada por no ser exacta.
- **Problema inverso** $\mathbf M\boldsymbol\phi = \mu\,\mathbf K\boldsymbol\phi$ con $\mu = 1/\omega^2$: admite $\mathbf M$ semidefinida, pero exige $\mathbf K$ definida, así que falla en estructuras libres, y pierde el desplazamiento espectral.
- **ARPACK con $\mathbf M$ semidefinida directamente**: el producto interno de $\mathbf M$ es degenerado y es justo lo que fallaba (`-9999`).

## Consecuencias

- Análisis modal y espectral de edificios con masas de piso, con columnas y vigas sin masa propia o con ella.
- Las frecuencias son exactas para el modelo discreto: la condensación no aproxima nada.
- El error de ARPACK sin causa se sustituye por un error del programa que nombra los grados de libertad.
- Los modelos con masa en todos los grados de libertad dan los mismos resultados que antes.

## Referencias

- Bathe, K.-J. (2014), *Finite Element Procedures*, §10.3.1 (condensación estática en problemas de valores propios).
- Chopra, A. K., *Dynamics of Structures*, §9.3 (condensación estática de los grados de libertad sin masa).
- OpenSees, comando `mass`; ANSYS, elemento `MASS21`; Abaqus, `*MASS`; SAP2000, masas de nudo.
- ADR 0008 (densidad como propiedad del material), ADR 0009 (análisis modal y dinámico).
