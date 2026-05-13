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
- $\mathbf M$ simétrica, **positiva definida** estrictamente para masa consistente (toda fila/columna tiene aporte de algún elemento con $\rho > 0$).
- Espectro real no negativo: $\omega_n \in \mathbb R_{\ge 0}$.

### 4. Salidas
- Frecuencias naturales $\omega_n$ (rad/s), $f_n = \omega_n/(2\pi)$ (Hz), períodos $T_n = 1/f_n$ (s).
- Modos $\boldsymbol\phi_n$ M-ortonormales: $\boldsymbol\phi_i^\top \mathbf M\, \boldsymbol\phi_j = \delta_{ij}$.
- Ordenación ascendente: $\omega_1 \le \omega_2 \le \ldots \le \omega_{n_\text{modes}}$.

---

## Formulación numérica

### 5. Reducción por Dirichlet
Por eliminación directa (ADR 0004) con operador $\mathbf T \in \mathbb R^{N_\text{dof} \times N_\text{libre}}$ que selecciona DOFs libres:
$$\mathbf K_\text{red} = \mathbf T^\top \mathbf K\, \mathbf T, \qquad \mathbf M_\text{red} = \mathbf T^\top \mathbf M\, \mathbf T.$$
El problema reducido $\mathbf K_\text{red}\,\boldsymbol\phi_\text{red} = \omega^2\, \mathbf M_\text{red}\,\boldsymbol\phi_\text{red}$ se resuelve sobre los DOFs libres. El término $\mathbf g_\text{indep}$ de restricciones lineales no homogéneas no aplica: los modos son solución del problema homogéneo (los modos de la estructura con apoyos a desplazamiento prescrito no nulo coinciden con los de apoyos nulos; el campo prescrito solo entra en la respuesta estática). Tras resolver, expansión $\boldsymbol\phi = \mathbf T\,\boldsymbol\phi_\text{red}$; en DOFs prescritos $\phi_i = 0$.

### 6. Algoritmo numérico: Lanczos con shift-invert
ARPACK vía `scipy.sparse.linalg.eigsh`:
$$\text{eigsh}(\mathbf K_\text{red},\; k = n_\text{modes},\; \mathbf M = \mathbf M_\text{red},\; \sigma,\; \text{which="LM"},\; \text{tol}).$$
- **Shift-invert** transforma el problema en $\bigl(\mathbf K_\text{red} - \sigma\,\mathbf M_\text{red}\bigr)^{-1}\mathbf M_\text{red}\,\boldsymbol\phi = \mu\,\boldsymbol\phi$ con $\mu = 1/(\omega^2 - \sigma)$. Buscar `which="LM"` (mayor magnitud de $\mu$) equivale a buscar autovalores cercanos al shift $\sigma$.
- $\sigma = 0$ por defecto: las frecuencias más bajas, que son las relevantes en ingeniería estructural.
- Internamente ARPACK factoriza $(\mathbf K_\text{red} - \sigma\,\mathbf M_\text{red})$ una sola vez y la reusa en cada iteración Lanczos.

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
  output: ModalResult   # campos: frequencies_hz, frequencies_rad, periods, modes, n_modes

parameters:
  - { name: n_modes,   type: int,   required: true,  desc: "Número de modos a calcular" }
  - { name: sigma,     type: float, required: false, default: 0.0,     desc: "Shift (rad²/s²) para shift-invert; 0 → frecuencias más bajas" }
  - { name: which,     type: str,   required: false, default: "LM",    desc: "Estrategia ARPACK: LM | SM | LA | SA | BE" }
  - { name: tolerance, type: float, required: false, default: 1.0e-9,  desc: "Tolerancia ARPACK sobre el residuo del autovalor" }
  - { name: lumping,   type: str,   required: false, default: "consistent", desc: "Discretización de masa; solo 'consistent' en fase 1" }
  - { name: linear_algebra, type: str, required: false, default: "auto", desc: "Backend para la factorización de (K - σM) interna a shift-invert (ADR 0003)" }

requirements:
  - "Todos los materiales del modelo declaran `density > 0` (ADR 0008). Densidad cero genera M_red singular → eigsh con sigma=0 falla; el solver propaga el error con mensaje físico."
  - "Al menos un DOF libre tras Dirichlet."
  - "n_modes < N_libre (ARPACK exige k < n)."

conventions:
  units:        "heredadas del modelo (ADR 0008); ω en rad/s, f en Hz, T en s consistentes con kg/N/m o equivalentes."
  frequencies:  "ω y f expuestos simultáneamente; el usuario consulta cualquiera."
  modes:        "M-ortonormales (Φᵀ M Φ = I); signo arbitrario por modo."
  ordering:     "ascendente en ω: ω₁ ≤ ω₂ ≤ … ≤ ω_n_modes."
  rigid_body:   "modos de cuerpo rígido aparecen con ω≈0; T=∞ se reporta sin error."

out_of_scope:
  - "Masa lumped (fase 2 del ADR 0009)."
  - "Masas modales efectivas y factores de participación (fase 7, análisis sísmico)."
  - "Modos complejos con amortiguamiento no proporcional."
  - "Análisis modal sobre estado deformado (modos de pequeña amplitud alrededor de configuración prestressed) — diferido."

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
    - name: density_faltante_lanza_value_error
      setup: "Material sin density declarada en el modelo."
      expect: "ModalSolver.solve() lanza ValueError listando los materiales afectados (mismo patrón que assemble_self_weight, ADR 0008)."

references:
  - "Bathe K.-J., Finite Element Procedures (2014), cap. 9 (matriz de masa consistente) y cap. 11 (algoritmos de autovalor)."
  - "Cook, Malkus & Plesha, Concepts and Applications of Finite Element Analysis, §11.3–§11.5."
  - "Hughes T. J. R., The Finite Element Method, cap. 7."
  - "Wilson E. L., Three-Dimensional Static and Dynamic Analysis of Structures, cap. 10."
  - "Lehoucq, Sorensen, Yang — ARPACK Users' Guide (Lanczos con shift-invert)."
  - "ADR 0003 — Capa algebraica (factorización reutilizable de (K − σM))."
  - "ADR 0008 — Material.density."
  - "ADR 0009 — Subsistema de análisis dinámico (este solver es la fase 1)."
```

---

## Implementación

- **Archivo**: [fenix/math/solvers.py](../../fenix/math/solvers.py) (clase `ModalSolver`, registrada vía `@SolverRegistry.register`).
- **Backend algebraico**: [fenix/math/linalg/eigen.py](../../fenix/math/linalg/eigen.py) (clase `EigenSolver` envolviendo `scipy.sparse.linalg.eigsh`).
- **Tipo de resultado**: [`ModalResult`](../../fenix/results.py) — dataclass frozen con `frequencies_rad`, `frequencies_hz`, `periods`, `modes`, `n_modes`, `converged`.
- **Entrypoint público**: [`fenix.run_modal`](../../fenix/entry.py) para uso programático; [`fenix.run_yaml`](../../fenix/entry.py) despacha automáticamente a `run_modal` cuando el YAML declara `solver: type: ModalSolver`.
- **Pipeline**: `Assembler.assemble_system()` → `Assembler.assemble_mass_matrix()` → `Assembler.reduce_pair(K, M)` → `EigenSolver.solve(K_red, M_red, n_modes)` → expansión `Φ = T · Φ_red` → `ModalResult`.
- **Caché**: `Assembler` cachea `M_global` con el lumping usado; reusos posteriores del mismo análisis no recomputan.
- **Validación de density**: pre-chequeo agregado en `Assembler.assemble_mass_matrix` con el mismo patrón de [ADR 0008 — `assemble_self_weight`](../adr/0008-densidad-propiedad-del-material.md).
- **Tests**:
  - [tests/test_modal.py](../../tests/test_modal.py) — 13 tests cubriendo:
    - `TestModalAxialBar`: barra empotrada-libre, 20 elementos Truss2D, frecuencias contra `ω_n = ((2n-1)π/(2L))·√(E/ρ)` (tol 1%).
    - `TestModalSimplySupportedBeam`: viga biapoyada, 20 elementos Frame2DEuler, frecuencias contra `ω_n = (nπ/L)²·√(EI/(ρA))` (tol 2%).
    - `TestModalAlgebraicProperties`: ortonormalidad `ΦᵀMΦ = I` (atol 1e-10) y ortogonalidad `ΦᵀKΦ = diag(ω²)` (rtol 1e-8) en ambos modelos.
    - `TestModalSolverContract`: errores agregados (`density` faltante, `lumping="lumped"` no implementado, `run_modal` sin solver ni n_modes).
    - `TestModalYamlPipeline`: pipeline end-to-end desde YAML.
- **Notas de traducción**:
  - El parámetro `linear_algebra` está en el constructor por coherencia con los demás solvers, pero en fase 1 no se usa: ARPACK gestiona internamente la factorización LU de `(K − σM)` para shift-invert. Cuando se enchufe a la capa algebraica (ADR 0003), se podrá inyectar Cholesky o LDLᵀ.
  - El bloque `convergence:` del YAML (ADR 0007) se omite silenciosamente para `ModalSolver` igual que para `LinearSolver`: no hay iteración Newton ni norma de residuo configurable; la tolerancia es la propia de ARPACK (`tolerance:`).
  - Matriz de masa consistente "translacional invariante rotacional" para `Truss*` y `Cable*` (`M = kron([[2,1],[1,2]] · ρAL/6, I_dim)`). Para frames se ensambla la masa local axial + Hermitiana cúbica y se rota con `T^T · M_local · T`. Para Frame3D se incluye masa torsional con momento polar geométrico `Jp = Iy + Iz` (no la constante de Saint-Venant `J`, que rige la rigidez). Para sólidos 2D, integración por la cuadratura del propio elemento (Tri3 usa fórmula analítica exacta porque su cuadratura 1-punto sería insuficiente).
  - Inercia rotacional propia de sección (`ρI·L`) **incluida** en los DOFs rotacionales de Frame2D (Euler, Timoshenko, EulerCorot) y Frame3D — coherente con Timoshenko y con la rigidez `K_local` que tiene términos rotacionales puros (ADR 0009 §1).

---

## Diálogo

- **2026-05-13** · Validado. Los 13 tests de `tests/test_modal.py` cubren los criterios de `acceptance` (verificación analítica + propiedades algebraicas + contrato de errores + pipeline YAML) y pasan con las tolerancias declaradas. Promovido `status: validated`.
