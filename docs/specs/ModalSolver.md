# ModalSolver — Análisis modal por autovalores generalizados

> Orden de trabajo. La especificación física y la formulación numérica son patrimonio del usuario (revisa con detalle). La IA implementa y rellena las secciones marcadas.

---

## Especificación física

### 0. Descripción general
Análisis modal lineal de una estructura discretizada por FEM: cálculo de frecuencias naturales y modos de vibración no amortiguados. Hipótesis: linealidad geométrica y material, amortiguamiento nulo (`C = 0`), masa constante en el tiempo (lagrangeano total).

### 1. Problema físico
Vibración libre no amortiguada del sistema discreto:
$$\mathbf M\,\ddot{\mathbf u}(t) + \mathbf K\,\mathbf u(t) = \mathbf 0.$$
Solución armónica $\mathbf u(t) = \boldsymbol\phi\, e^{i\omega t}$. Sustituyendo:
$$\bigl(\mathbf K - \omega^2 \mathbf M\bigr)\,\boldsymbol\phi = \mathbf 0.$$

### 2. Problema de autovalor generalizado
Solución no trivial requiere $\det(\mathbf K - \omega^2 \mathbf M) = 0$. Equivalente:
$$\mathbf K\,\boldsymbol\phi_n = \omega_n^2\, \mathbf M\,\boldsymbol\phi_n, \qquad n = 1, \ldots, N_{\text{dof}}.$$
Autovalores $\lambda_n = \omega_n^2 \ge 0$; autovectores $\boldsymbol\phi_n$ son los modos de vibración.

### 3. Propiedades del par $(\mathbf K, \mathbf M)$
- $\mathbf K$ simétrica, positiva (semi)definida tras Dirichlet.
- $\mathbf M$ simétrica, **positiva semidefinida**. Es definida cuando todo grado de libertad libre recibe masa (de algún elemento con $\rho > 0$ o de una masa de nudo). Si no, los grados de libertad sin masa se condensan (§6 bis).
- Fuentes de masa: la densidad del material (ADR 0008) y las **masas concentradas en nudos** (ADR 0024), que se suman a la diagonal de $\mathbf M$ igual con masa consistente que concentrada.
- Espectro real no negativo: $\omega_n \in \mathbb R_{\ge 0}$.

### 4. Salidas
- Frecuencias naturales $\omega_n$ (rad/s), $f_n = \omega_n/(2\pi)$ (Hz), períodos $T_n = 1/f_n$ (s).
- Modos $\boldsymbol\phi_n$ M-ortonormales: $\boldsymbol\phi_i^\top \mathbf M\, \boldsymbol\phi_j = \delta_{ij}$.
- Ordenación ascendente: $\omega_1 \le \omega_2 \le \ldots \le \omega_{n_\text{modes}}$.
- Participación modal en cada traslación del modelo ($u_x$, $u_y$, $u_z$): factores $\Gamma_n = \boldsymbol\phi_n^\top \mathbf M\,\mathbf r$ y masas efectivas $M_n^* = \Gamma_n^2$, con $\mathbf r$ la traslación rígida unitaria (vale 1 también en los apoyos, que se mueven con el suelo; Chopra 2017, §13.1). Con modos M-ortonormales no hace falta dividir entre $\boldsymbol\phi_n^\top\mathbf M\boldsymbol\phi_n$.
- Los grados de libertad sin masa condensados (§6 bis), como (nodo, grado de libertad).

---

## Formulación numérica

### 5. Reducción por Dirichlet
Por eliminación directa (ADR 0004) con operador $\mathbf T \in \mathbb R^{N_\text{dof} \times N_\text{libre}}$ que selecciona DOFs libres:
$$\mathbf K_\text{red} = \mathbf T^\top \mathbf K\, \mathbf T, \qquad \mathbf M_\text{red} = \mathbf T^\top \mathbf M\, \mathbf T.$$
El problema reducido $\mathbf K_\text{red}\,\boldsymbol\phi_\text{red} = \omega^2\, \mathbf M_\text{red}\,\boldsymbol\phi_\text{red}$ se resuelve sobre los DOFs libres. El término $\mathbf g_\text{indep}$ de restricciones lineales no homogéneas no aplica: los modos son solución del problema homogéneo (los modos de la estructura con apoyos a desplazamiento prescrito no nulo coinciden con los de apoyos nulos; el campo prescrito solo entra en la respuesta estática). Tras resolver, expansión $\boldsymbol\phi = \mathbf T\,\boldsymbol\phi_\text{red}$; en DOFs prescritos $\phi_i = 0$.

### 6. Algoritmo numérico: Lanczos con shift-invert
ARPACK vía `scipy.sparse.linalg.eigsh`:
$$\text{eigsh}(\mathbf K_\text{red},\; k = n_\text{modes},\; \mathbf M = \mathbf M_\text{red},\; \sigma,\; \text{which}\!=\!\text{`LM'},\; \text{tol}).$$
- **Shift-invert** transforma el problema en $\bigl(\mathbf K_\text{red} - \sigma\,\mathbf M_\text{red}\bigr)^{-1}\mathbf M_\text{red}\,\boldsymbol\phi = \mu\,\boldsymbol\phi$ con $\mu = 1/(\omega^2 - \sigma)$. Buscar `which="LM"` (mayor magnitud de $\mu$) equivale a buscar autovalores cercanos al shift $\sigma$.
- $\sigma = 0$ por defecto: las frecuencias más bajas, que son las relevantes en ingeniería estructural.
- Internamente ARPACK factoriza $(\mathbf K_\text{red} - \sigma\,\mathbf M_\text{red})$ una sola vez y la reusa en cada iteración Lanczos.

### 6 bis. Grados de libertad sin masa: condensación estática exacta (ADR 0024)
Con masas sólo en algunos nudos o elementos de densidad cero, hay grados de libertad libres con diagonal de $\mathbf M_\text{red}$ nula (cero exacto; en una matriz semidefinida positiva eso implica fila y columna nulas). Con $m$ los grados de libertad con masa y $0$ los demás, la fila $0$ no tiene inercia:
$$\boldsymbol\phi_0 = -\mathbf K_{00}^{-1}\,\mathbf K_{0m}\,\boldsymbol\phi_m, \qquad \hat{\mathbf K}\,\boldsymbol\phi_m = \omega^2\,\mathbf M_{mm}\,\boldsymbol\phi_m, \qquad \hat{\mathbf K} = \mathbf K_{mm} - \mathbf K_{m0}\,\mathbf K_{00}^{-1}\,\mathbf K_{0m}.$$
Es exacta, no una reducción de Guyan aproximada: los valores propios de $\hat{\mathbf K}$ son las frecuencias finitas del problema original, tantas como grados de libertad con masa; las demás son infinitas.
- Con 32 o menos grados de libertad con masa, o pidiendo casi todos los modos: $\hat{\mathbf K}$ densa y `eigh`.
- Con más: ARPACK sobre el problema condensado sin formar $\hat{\mathbf K}$. Con $\mathbf A = \mathbf K - \sigma\mathbf M$, el bloque $mm$ de $\mathbf A^{-1}$ es $(\hat{\mathbf K} - \sigma\,\mathbf M_{mm})^{-1}$ (complemento de Schur), así que el operador de desplazamiento e inversión usa la factorización dispersa de $\mathbf A$ del modelo completo, con el producto interno de $\mathbf M_{mm}$, que es definida. Con $\sigma = 0$ y $\mathbf K$ singular se reintenta con desplazamiento negativo, como en el camino sin condensar.
- $\boldsymbol\phi_0$ se recupera con la factorización de $\mathbf K_{00}$; los modos expandidos son $\mathbf M$-ortonormales.
- Si no hay solución (más modos que grados de libertad con masa, ninguno con masa, grados de libertad sin masa sin rigidez que los fije, o masa singular más allá de ellos), se lanza `MasslessDofsError`, que los nombra como (nodo, dof).

### 7. Normalización
`eigsh` con $\mathbf M$ explícita devuelve autovectores M-ortonormales sin reescalado adicional. Signo arbitrario (autovector definido módulo signo); convención de signo se documenta solo si emerge como preferencia (no parte del contrato).

### 8. Post-procesado
- $\omega_n = \sqrt{\lambda_n}$ (clip a cero los autovalores ligeramente negativos por ruido numérico, $|\lambda| < 10^{-14}\,\lambda_\text{max}$).
- $f_n = \omega_n / (2\pi)$, $T_n = 1/f_n$ (con $T_n = \infty$ si $\omega_n = 0$, modo de cuerpo rígido).
- Reordenación ascendente por $\omega_n$.

---

## Contrato de implementación

```yaml
name: ModalSolver
kind: solver
status: validated       # draft → implemented → validated

interface:
  yaml_type: modal      # solver: { type: modal }
  output: ModalResult   # campos: frequencies_hz, frequencies_rad, periods, modes, n_modes, participation_factors, effective_masses, condensed_dofs

admission:                # criterio de admisión al programa principal (Reglas.md §1)
  textbooks:
    - bathe: "cap. 10-11 (autoproblemas generalizados)"
    - hughes: "cap. 10 (métodos para autoproblemas)"
  general_codes:
    - abaqus: "*FREQUENCY (Lanczos)"
    - ansys: "MODOPT,LANB"
    - opensees: "eigen"

parameters:
  - { name: n_modes,   type: int,   required: true,  desc: "Número de modos a calcular" }
  - { name: sigma,     type: float, required: false, default: 0.0,     desc: "Shift (rad²/s²) para shift-invert; 0 → frecuencias más bajas" }
  - { name: which,     type: str,   required: false, default: "LM",    desc: "Estrategia ARPACK: LM | SM | LA | SA | BE" }
  - { name: tolerance, type: float, required: false, default: 1.0e-9,  desc: "Tolerancia ARPACK sobre el residuo del autovalor" }
  - { name: lumping,   type: str,   required: false, default: "consistent", desc: "Discretización de masa: 'consistent' o 'lumped' (ADR 0009 fases 1 y 2, cerradas)" }
  - { name: linear_algebra, type: str, required: false, default: "auto", desc: "Backend para la factorización de (K - σM) interna a shift-invert (ADR 0003)" }

requirements:
  - "Todos los materiales del modelo declaran `density` (ADR 0008); `density = 0.0` es válida para elementos sin masa propia."
  - "Al menos un DOF libre con masa (de la densidad o de una masa de nudo, ADR 0024)."
  - "n_modes ≤ número de DOF libres con masa: sólo ellos tienen frecuencia finita. Los DOF sin masa se condensan estáticamente (§6 bis); si no se puede, `MasslessDofsError`."

conventions:
  units:        "heredadas del modelo (ADR 0008); ω en rad/s, f en Hz, T en s consistentes con kg/N/m o equivalentes."
  frequencies:  "ω y f expuestos simultáneamente; el usuario consulta cualquiera."
  modes:        "M-ortonormales (Φᵀ M Φ = I); signo arbitrario por modo."
  ordering:     "ascendente en ω: ω₁ ≤ ω₂ ≤ … ≤ ω_n_modes."
  rigid_body:   "modos de cuerpo rígido aparecen con ω≈0; T=∞ se reporta sin error."

out_of_scope:
  - "Modos complejos con amortiguamiento no proporcional."
  - "Análisis modal sobre estado deformado (modos de pequeña amplitud alrededor de configuración prestressed) — diferido."
  - "Masa lumped y factores de participación / masas efectivas — implementados (fase 2 y fase 7 del ADR 0009, cerradas 2026-05-18); accesibles vía `lumping='lumped'` y `ResponseSpectrumSolver` respectivamente."

acceptance:
  verification:
    - name: barra_axial_empotrada_libre
      setup: |
        Barra 1D elástica empotrada en x=0, libre en x=L. Discretización con
        N=20 elementos Truss2D alineados con el eje x. E=210e9 Pa, ρ=7850 kg/m³,
        A=1e-4 m², L=1.0 m.
      expect: |
        Frecuencias propias analíticas (ondas axiales en barra empotrada-libre):
          ω_n = ((2n-1)π / (2L)) · sqrt(E/ρ),    n = 1, 2, 3, …
        Los primeros 3 modos calculados coinciden con la solución analítica.
      tol_rel: 0.01
    - name: viga_bernoulli_simplemente_apoyada
      setup: |
        Viga Bernoulli-Euler simplemente apoyada en sus dos extremos.
        Discretización con N=20 elementos Frame2DEuler alineados con el eje x.
        E=210e9 Pa, ρ=7850 kg/m³, A=1e-4 m², I=8.33e-10 m⁴, L=1.0 m.
      expect: |
        Frecuencias propias analíticas (flexión transversal):
          ω_n = (nπ/L)² · sqrt(E·I / (ρ·A)),    n = 1, 2, 3, …
        Los primeros 3 modos coinciden con la solución analítica.
      tol_rel: 0.02
    - name: oscilador_un_gdl_masa_de_nudo
      setup: |
        Columna en voladizo de Frame2DEuler sin masa propia (density = 0.0),
        con 1 y 4 elementos; masa m en el ux de la cabeza (ADR 0024).
        E = 2.0e11 Pa, I = 8.5e-5 m⁴, L = 3.5 m, m = 2.5e4 kg.
      expect: |
        T = 2π·√(m/k) con k = 3EI/L³, con masa consistente y concentrada.
        Los DOF sin masa del modo siguen la elástica de la carga en la punta:
        v(y)/v(L) = y²(3L − y)/(2L³), θ(L) = −3/(2L)·v(L). Con una inercia J
        en el rz de la cabeza, T = 2π·√(J/(EI/L)).
      tol_rel: 1.0e-6
      ref: "tests/test_nodal_mass.py::TestOsciladorDeUnGradoDeLibertad"
    - name: edificio_cortante_tres_pisos
      setup: |
        Tres niveles, dos columnas biempotradas de Frame2DEuler sin masa por
        nivel (2 elementos cada una), vigas rígidas impuestas por cinemática
        (uy = rz = 0 en los nudos de piso y ux del nudo derecho igual al del
        izquierdo por restricción lineal) y la masa de cada piso repartida
        por mitades entre sus dos nudos.
      expect: |
        Periodos y modos iguales a los del problema 3×3 a mano, rigidez de
        entrepiso 2·12EI/h³: con pisos iguales, la forma cerrada
        ω_j = 2·√(k/m)·sin((2j−1)π/14), φ_j(i) ∝ sin(i(2j−1)π/7); con pisos
        distintos, el 3×3 con sus masas y rigideces.
      tol_rel: 1.0e-6
      ref: "tests/test_nodal_mass.py::TestEdificioDeCortante"
  specific:
    - name: ortonormalidad_de_masa
      setup: "Cualquier modelo bien formado tras resolver el problema modal."
      expect: "Φᵀ M Φ = I (n_modes × n_modes)."
      tol_abs: 1.0e-10
    - name: ortogonalidad_de_rigidez
      setup: "Mismo modelo."
      expect: "Φᵀ K Φ = diag(ω²₁, …, ω²_n)."
      tol_rel: 1.0e-8
    - name: dirichlet_anula_modos_en_prescritos
      setup: "Barra empotrada-libre; nodo x=0 con ux=0."
      expect: "Componente del modo en el DOF prescrito = 0 exacto."
      tol_abs: 1.0e-14
    - name: mpc_tie_matches_merged_node
      setup: "ménsula de Frame2DEuler inclinada 30° partida en su nodo central y atada con u_s = u_m en ux, uy y rz, con el esclavo a cada lado del corte; masa consistente y concentrada"
      expect: "Frecuencias y modos expandidos Φ = T·Φ_red (también en el nodo duplicado) = los del modelo de nodo fundido; ΦᵀMΦ = I con la masa global del modelo atado"
      tol_rel: 1.0e-8
      ref: "tests/test_mpc_dynamic.py::test_modal_reproduce_el_nodo_fundido"
    - name: density_faltante_lanza_value_error
      setup: "Material sin density declarada en el modelo."
      expect: "ModalSolver.solve() lanza ValueError listando los materiales afectados (mismo patrón que assemble_self_weight, ADR 0008)."
    - name: masa_de_nudo_entrada_exacta
      setup: "Columna en voladizo de Frame2DEuler con densidad; masa m en el ux del nudo de la cabeza (ADR 0024); masa consistente y concentrada"
      expect: "M_uxux aumenta exactamente en m y ninguna otra entrada cambia; lo mismo con una inercia rotacional en rz"
      tol_abs: 0.0
      ref: "tests/test_nodal_mass.py::TestEnsamblajeDeLaMasaDeNudo"
    - name: dof_sin_masa_marco_de_un_nivel
      setup: "Marco de acero de un nivel (h = 3 m, L = 5 m, columnas 0.30 × 0.50 m de 20 fibras y viga 0.40 × 0.70 m de 28 fibras, 4 elementos por miembro, 3 puntos de Gauss) con columnas de densidad cero: 18 de 33 DOF libres sin masa"
      expect: "T₁ = 0.03097 s y T₂ = 0.01295 s, iguales a la solución densa (valores propios de K⁻¹M); residuo relativo de cada modo ‖Kφ − ω²Mφ‖/‖Kφ‖ ≤ 1e-8; con 12 elementos por miembro (39 DOF con masa) el mismo criterio por el camino ARPACK, con masa consistente y concentrada"
      tol_rel: 1.0e-6
      ref: "tests/test_modal_massless_dofs.py::TestMarcoDeUnNivelSinMasaEnColumnas"
    - name: dof_sin_masa_error_tipado
      setup: "Voladizo sin masa propia con masa sólo en el ux de la cabeza y n_modes = 2; modelo sin ninguna masa; dos barras alineadas cuyo nudo intermedio no tiene masa ni rigidez transversal"
      expect: "MasslessDofsError (ValueError) con causa too_many_modes, no_mass y massless_mechanism; el mensaje nombra los DOF sin masa como nodo <id> (<dof>)"
      ref: "tests/test_modal_massless_dofs.py::TestErrorTipado"
    - name: masa_en_todos_sin_cambios
      setup: "Modelos con masa en todos los DOF libres (toda la suite modal y dinámica previa)"
      expect: "Siguen el camino de siempre, sin condensación, y dan los mismos resultados que antes del ADR 0024"
      ref: "tests/test_modal_massless_dofs.py::TestMarcoDeUnNivelSinMasaEnColumnas::test_con_masa_en_todos_sigue_el_camino_de_siempre"

references:
  - "Bathe K.-J., Finite Element Procedures (2014), cap. 9 (matriz de masa consistente) y cap. 11 (algoritmos de autovalor)."
  - "Cook, Malkus & Plesha, Concepts and Applications of Finite Element Analysis, §11.3–§11.5."
  - "Hughes T. J. R., The Finite Element Method, cap. 7."
  - "Wilson E. L., Three-Dimensional Static and Dynamic Analysis of Structures, cap. 10."
  - "Lehoucq, Sorensen, Yang — ARPACK Users' Guide (Lanczos con shift-invert)."
  - "ADR 0003 — Capa algebraica (factorización reutilizable de (K − σM))."
  - "ADR 0008 — Material.density."
  - "ADR 0024 — Masas concentradas en nudos y análisis modal con grados de libertad sin masa."
  - "Bathe K.-J., Finite Element Procedures (2014), §10.3.1 (condensación estática en problemas de valores propios)."
  - "Chopra A. K., Dynamics of Structures, §9.3 (condensación estática de los grados de libertad sin masa)."
  - "ADR 0009 — Análisis modal y dinámico (este solver es la fase 1: el análisis modal estricto, problema generalizado K·φ = ω²·M·φ; el análisis dinámico transitorio entra en fases 3-7)."
```

---

## Implementación

- **Archivo**: [solidum/math/solvers/modal.py](../../solidum/math/solvers/modal.py) (clase `ModalSolver`, registrada vía `@SolverRegistry.register`).
- **Backend algebraico**: [solidum/math/linalg/eigen.py](../../solidum/math/linalg/eigen.py) (clase `EigenSolver` envolviendo `scipy.sparse.linalg.eigsh`).
- **Tipo de resultado**: [`ModalResult`](../../solidum/results.py) — dataclass frozen con `frequencies_rad`, `frequencies_hz`, `periods`, `modes`, `n_modes`, `converged`. Expone `free_vibration(M, u0, u0_dot, t)` para reconstruir analíticamente la respuesta temporal de vibración libre sin amortiguamiento por superposición modal — solución cerrada de `M·ü + K·u = 0` proyectada sobre los modos calculados (modos rígidos tratados aparte con evolución lineal `q_n(t) = a_n + b_n·t`).
- **Entrypoint público**: [`solidum.run_modal`](../../solidum/entry.py) para uso programático; [`solidum.run_yaml`](../../solidum/entry.py) despacha automáticamente a `run_modal` cuando el YAML declara `solver: type: ModalSolver`.
- **Pipeline**: `Assembler.assemble_system()` → `Assembler.assemble_mass_matrix()` → `Assembler.reduce_pair(K, M)` → `EigenSolver.solve(K_red, M_red, n_modes)` → expansión `Φ = T · Φ_red` → `ModalResult`.
- **Caché**: `Assembler` cachea `M_global` con el lumping usado; reusos posteriores del mismo análisis no recomputan.
- **Validación de density**: pre-chequeo agregado en `Assembler.assemble_mass_matrix` con el mismo patrón de [ADR 0008 — `assemble_self_weight`](../adr/0008-densidad-propiedad-del-material.md).
- **Tests**:
  - [tests/test_modal.py](../../tests/test_modal.py) — 13 tests cubriendo:
    - `TestModalAxialBar`: barra empotrada-libre, 20 elementos Truss2D, frecuencias contra `ω_n = ((2n-1)π/(2L))·√(E/ρ)` (tol 1%).
    - `TestModalSimplySupportedBeam`: viga biapoyada, 20 elementos Frame2DEuler, frecuencias contra `ω_n = (nπ/L)²·√(EI/(ρA))` (tol 2%).
    - `TestModalAlgebraicProperties`: ortonormalidad `ΦᵀMΦ = I` (atol 1e-10) y ortogonalidad `ΦᵀKΦ = diag(ω²)` (rtol 1e-8) en ambos modelos.
    - `TestModalSolverContract`: errores agregados (`density` faltante, `run_modal` sin solver ni n_modes) y verificación de que `lumping="lumped"` corre tras ADR 0009 fase 2 (cerrada 2026-05-18).
    - `TestModalYamlPipeline`: pipeline end-to-end desde YAML.
- **Notas de traducción**:
  - El parámetro `linear_algebra` está en el constructor por coherencia con los demás solvers, pero en fase 1 no se usa: ARPACK gestiona internamente la factorización LU de `(K − σM)` para shift-invert. Cuando se enchufe a la capa algebraica (ADR 0003), se podrá inyectar Cholesky o LDLᵀ.
  - El bloque `convergence:` del YAML (ADR 0007) se omite silenciosamente para `ModalSolver` igual que para `LinearSolver`: no hay iteración Newton ni norma de residuo configurable; la tolerancia es la propia de ARPACK (`tolerance:`).
  - Matriz de masa consistente "translacional invariante rotacional" para `Truss*` y `Cable*` (`M = kron([[2,1],[1,2]] · ρAL/6, I_dim)`). Para frames se ensambla la masa local axial + Hermitiana cúbica y se rota con `T^T · M_local · T`. Para Frame3D se incluye masa torsional con momento polar geométrico `Jp = Iy + Iz` (no la constante de Saint-Venant `J`, que rige la rigidez). Para sólidos 2D, integración por la cuadratura del propio elemento (Tri3 usa fórmula analítica exacta porque su cuadratura 1-punto sería insuficiente).
  - Inercia rotacional propia de sección (`ρI·L`) **incluida** en los DOFs rotacionales de Frame2D (Euler, Timoshenko, EulerCorot) y Frame3D — coherente con Timoshenko y con la rigidez `K_local` que tiene términos rotacionales puros (ADR 0009 §1).

---

## Diálogo

- **2026-05-13** · Validado. Los 13 tests de `tests/test_modal.py` cubren los criterios de `acceptance` (verificación analítica + propiedades algebraicas + contrato de errores + pipeline YAML) y pasan con las tolerancias declaradas. Promovido `status: validated`.
- **2026-09-22** · Auditoría global: `K` singular con σ = 0 (estructura free-free) fallaba según la malla ("Factor is exactly singular" con 1, 2 o 4 elementos, éxito por redondeo con 3 o 10): reintento automático con un shift negativo relativo a tr(K)/tr(M). `ArpackNoConvergence` → `converged=False` con los modos obtenidos (`n_modes` reales); antes el flag estaba cableado a `True`. Modelos con n_free ≤ 32 o n_modes ≥ n_free − 1 usan `scipy.linalg.eigh` denso (M-ortonormal), así que un modelo de 1 GDL admite análisis modal.
- **2026-09-30** · Auditoría de las MPC (ADR 0004): las restricciones multipunto entran por `Assembler.reduce_pair` con el mismo `T` que los apoyos (`K_red = TᵀKT`, `M_red = TᵀMT`, `Φ = T·Φ_red`). Primer test con MPC: el oráculo del nodo fundido se cumple hasta el redondeo en frecuencias y modos, y los modos expandidos son M-ortonormales con la masa global (`tests/test_mpc_dynamic.py`). Nota sobre `ModalResult.free_vibration(M, u0, u0_dot, t)`: proyecta las condiciones iniciales en norma M (`a = ΦᵀM·u0`), así que su resultado no depende de qué nodo se declaró esclavo; los solvers transitorios, en cambio, conservan las componentes libres y reconstruyen los esclavos como `T·u0[libres] + g`. Con condiciones iniciales coherentes con las MPC ambos caminos coinciden; con incoherentes difieren (los transitorios lo avisan en el log). Unificar las dos proyecciones cambiaría resultados con condiciones iniciales incoherentes: queda como decisión física pendiente del usuario. Sin cambio de código.
- **2026-09-30** · Revisión del paquete dinámico de la auditoría de las MPC: el parámetro `lumping` del contrato aún decía «solo 'consistent' en fase 1», pero `ModalSolver` admite `lumping='lumped'` desde el cierre de la fase 2 del ADR 0009 (lo ejercitan `tests/test_modal.py::test_lumped_runs_after_phase2` y el caso `mpc_tie_matches_merged_node` con masa concentrada). Corregida la descripción aquí y en el catálogo. Sin cambio de código.
- **2026-11-09** · ADR 0024. Masas concentradas en nudos (`Domain.add_nodal_mass`, bloque YAML `nodal_masses`), sumadas a la diagonal de $\mathbf M$ con los dos esquemas de masa. Con grados de libertad libres sin masa, ARPACK lanzaba `ArpackError -9999` sin causa y el camino denso exigía $\mathbf M$ definida; ahora se condensan estáticamente de forma exacta (§6 bis) y, si el problema no tiene solución, `MasslessDofsError` los nombra. Criterios nuevos: `oscilador_un_gdl_masa_de_nudo`, `edificio_cortante_tres_pisos`, `masa_de_nudo_entrada_exacta`, `dof_sin_masa_marco_de_un_nivel`, `dof_sin_masa_error_tipado`, `masa_en_todos_sin_cambios`. Medido: error relativo de los periodos de 1e-12 a 1e-14 contra la solución densa, residuos de 1e-11. `ResponseSpectrumSolver` lo hereda (masa total y masas efectivas con masas de piso en `tests/test_nodal_mass.py`).
- **2026-11-09** · Participación modal en el resultado. `ModalResult` trae, para cada traslación del modelo, los factores de participación `participation_factors["ux"]` ($\Gamma_n = \boldsymbol\phi_n^\top\mathbf M\mathbf r$) y las masas efectivas `effective_masses["ux"]` ($\Gamma_n^2$), con la misma definición que ya usaba `ResponseSpectrumSolver`; y `condensed_dofs`, los grados de libertad sin masa condensados. Para un giro alrededor de un punto o una dirección oblicua: `solidum.math.rigid_body.rigid_body_vector` y `solidum.math.modal_response.participation_factors`. Sin cambio de física. La razón de masa acumulada no se añade todavía: con qué masa total se compara (la de todo el modelo o sólo la que puede moverse) es una decisión física pendiente, la misma que la del espectral. Criterios: edificio de cortante de tres pisos contra el problema 3×3 (1e-6), suma de las masas efectivas de todos los modos igual a $\mathbf b^\top\mathbf M_{LL}^{-1}\mathbf b$ con masa consistente y concentrada, igual que el espectral, $\mathbf K\mathbf r = \mathbf 0$ para los movimientos rígidos de cinco familias de elementos. Test: `tests/test_modal_participation.py`.
- **2026-11-09** · Masas de nudo por grupo físico y por coordenada en el YAML (`nodal_masses_by_group`, `nodal_masses_by_coord`): cada nodo seleccionado recibe la masa indicada, como las cargas puntuales por grupo. Sin cambio de física. Test: `tests/test_nodal_mass_selectors.py`.
- **2026-11-09** · Revisión de la masa: `solidum.math.mass_properties.mass_properties(asm, lumping)` da la masa total por traslación ($\mathbf r_d^\top\mathbf M\mathbf r_d$) y el centro de masa, por dirección como línea de acción de la resultante de inercia ($\bar x_j = c_j + \varepsilon_{dkj}\,\mathbf r_d^\top\mathbf M\mathbf r_k(c)/m_d$ con el giro rígido $\mathbf r_k$), exacto con masa consistente y concentrada también en barras inclinadas; un solo punto si la masa no depende de la dirección. Masas de nudo en el VTK. Sin cambio de física. Test: `tests/test_mass_properties.py` (masas de nudo, marco plano con barra inclinada, marco espacial, placas de densidad distinta, VTK).
