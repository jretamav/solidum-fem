# LinearSolver — Análisis estático lineal en un solo paso

> Spec **retroactiva**: el solver existe desde la primera fase del proyecto y está validado por la inmensa mayoría de tests del pipeline estático. Esta spec lo documenta sin cambiar comportamiento — H-5.3 de la auditoría 2026-05-18.

---

## Especificación física

### 0. Descripción general

Análisis estático **lineal** de una estructura discretizada por FEM: resolución directa del sistema algebraico `K·U = F`. Hipótesis: linealidad geométrica y material (tangente constante = rigidez secante), apoyos Dirichlet (homogéneos o no homogéneos), restricciones lineales multipunto opcionales (ADR 0004). Sin iteración. Sin historia de carga.

### 1. Ecuación de equilibrio resuelta

Forma fuerte semidiscreta:

$$\mathbf K\,\mathbf U \;=\; \mathbf F_{\text{ext}}$$

donde $\mathbf K$ es la matriz de rigidez global ensamblada, $\mathbf U$ los desplazamientos nodales (incluidos los Dirichlet) y $\mathbf F_{\text{ext}}$ las cargas externas (incluidas peso propio si está activado).

### 2. Condiciones de contorno

Apoyos Dirichlet (homogéneos $u_s = 0$ o no homogéneos $u_s = g$) eliminados directamente por ADR 0004 fase 1. Restricciones lineales multipunto admitidas vía transformación $\mathbf T$ + offset $\mathbf g$ (ADR 0004 fase 2). El sistema reducido se resuelve en DOFs libres:

$$\mathbf K_{\text{red}} = \mathbf T^\top \mathbf K\,\mathbf T, \qquad \mathbf F_{\text{red}} = \mathbf T^\top\,(\mathbf F - \mathbf K\,\mathbf g)$$

### 3. Salidas físicas

- Vector de desplazamientos global $\mathbf U \in \mathbb R^{n_{\text{dof}}}$ (libres + prescritos consistentes).
- Reacciones se recuperan a posteriori vía `domain.reactions(U)`.
- Esfuerzos en puntos de Gauss vía `domain.recompute_gauss_state(U)`.
- El resultado se envuelve en `SolveResult` por la capa pública (`fenix.run_static`).

---

## Formulación numérica

### 4. Esquema operativo

Un solo paso:

1. Ensamblar $\mathbf K$ (sparse COO → CSR) con caché topológica del Assembler.
2. Aplicar la reducción ADR 0004: $\mathbf K_{\text{red}}$, $\mathbf F_{\text{red}}$, operadores $\mathbf T$ y $\mathbf g$.
3. Calcular $\mathbf F_{\text{dir}} = \mathbf T^\top\,(\mathbf K\,\mathbf g)$ (contribución de Dirichlet no homogéneo).
4. Factorizar $\mathbf K_{\text{red}}$ con el backend algebraico seleccionado.
5. Resolver $\mathbf U_{\text{red}} = \mathbf K_{\text{red}}^{-1}\,(\mathbf F_{\text{red}} - \mathbf F_{\text{dir}})$.
6. Expandir $\mathbf U = \mathbf T\,\mathbf U_{\text{red}} + \mathbf g$.

### 5. Predictor / corrector

No aplica — solver no iterativo.

### 6. Criterio de convergencia

No aplica — la solución directa es exacta (a precisión de máquina y condicionamiento de $\mathbf K$).

### 7. Imposición de Dirichlet y MPC

ADR 0004 fase 1 (eliminación directa de apoyos) y fase 2 (transformación $\mathbf T$ + offset $\mathbf g$ para MPC). Manejados uniformemente por `Assembler.reduce(...)`.

### 8. Backend algebraico

Despacho automático ADR 0003:

- Si el dominio es **simétrico** (todos los materiales con tangente simétrica) ⇒ se asume PD y se factoriza con **Cholesky** (`scipy.sparse.linalg.splu` con `options={'SymmetricMode': True}` o equivalente Cholesky sparse).
- Si Cholesky reporta no-positividad (puede ocurrir con casi-singular o mal condicionado) ⇒ **fallback automático** a LU.
- Si el dominio declara material **no simétrico** ⇒ LU directo.

Override manual: parámetro `linear_algebra ∈ {"auto", "cholesky", "lu"}`.

**Cache de factorización (ADR 0003 fase 2)**: la primera llamada a `solve` ensambla $\mathbf K$, reduce y factoriza; la factorización se guarda en el solver. Llamadas posteriores con un $\mathbf F$ distinto reutilizan el factor sin reensamblaje — el coste se reduce a una resolución triangular barata. Si el usuario modifica el modelo entre llamadas (nuevos elementos, BCs, materiales) debe invocar `invalidate_cache()` explícitamente.

### 9. Adaptatividad / control de paso

No aplica.

### 10. Caveats numéricos

- **Suposición de PD por simetría**: el despachador asume que un dominio con materiales simétricos tiene $\mathbf K$ positiva definida. Es cierto para elasticidad estándar con suficientes apoyos; falso para casos degenerados (apoyos insuficientes, modos rígidos no restringidos, casi-singularidad geométrica). El fallback Cholesky→LU cubre estos casos.
- **Cache desactualizado**: si el modelo se modifica entre `solve` calls sin llamar a `invalidate_cache`, los resultados serán inconsistentes. La política deliberada es **cache silenciosa** (no se detectan modificaciones del modelo automáticamente) por velocidad — la responsabilidad de invalidar es del usuario.
- **No es no lineal**: si el problema es no lineal (material plástico, geometría corotacional con cargas significativas) usar `NonlinearSolver`. Este solver evaluará la tangente en estado virgen ($\mathbf U = 0$) y producirá una respuesta lineal incorrecta sin avisar.

---

## Contrato de implementación

```yaml
name: LinearSolver
kind: solver
status: validated

interface:
  yaml_type: linear
  output: SolveResult
  pipeline_kind: static

parameters:
  - { name: linear_algebra, type: str, required: false, default: "auto",
      desc: "Selección del backend algebraico: 'auto' (Cholesky→LU según simetría), 'cholesky', 'lu'" }

requirements:
  - "Modelo lineal: materiales con tangente constante, geometría no corotacional o cargas suficientemente pequeñas"
  - "Apoyos suficientes para descartar modos rígidos (estructura cinemáticamente estable)"
  - "Compatibilidad ensamblaje-elemento-material verificada por Domain"

conventions:
  units: "heredadas del modelo (ADR 0008)"
  stability: "trivial — un solo paso"

out_of_scope:
  - "Análisis no lineal ⇒ usar NonlinearSolver"
  - "Snap-through / snap-back ⇒ usar ArcLengthSolver"
  - "Análisis dinámico ⇒ usar NewmarkSolver, HHTSolver, etc."
  - "Detección automática de modos rígidos (apoyos insuficientes)"

acceptance:
  verification:
    - name: viga_en_voladizo_carga_puntual
      setup: "viga de Frame2DEuler en voladizo con carga puntual P en el extremo"
      expect: "δ = P·L³/(3·E·I) ± O(elementos)"
      tol_rel: 1.0e-8
    - name: barra_traccion_uniaxial
      setup: "Truss2D, longitud L, área A, carga axial F"
      expect: "δ = F·L/(E·A) exacto (sin discretización)"
      tol_rel: 1.0e-12
    - name: convergencia_h
      setup: "viga discretizada con N elementos, N ∈ {2, 4, 8, 16}"
      expect: "error decrece con orden teórico del elemento (4 para Hermite cúbicos)"
      tol_rel: 0.1
    - name: reactions_equilibrio
      setup: "estructura cargada en una dirección"
      expect: "Σ(reacciones) + Σ(cargas externas) ≈ 0 (Newton 3ª ley)"
      tol_abs: 1.0e-9
    - name: dirichlet_no_homogeneo
      setup: "barra con apoyo prescrito u_s = g ≠ 0 en un extremo"
      expect: "solución exacta esperada (translación rígida + alargamiento si carga)"
      tol_rel: 1.0e-12
    - name: cache_factorizacion
      setup: "dos llamadas consecutivas a solve con F distintos sin invalidar cache"
      expect: "segunda llamada reutiliza el factor (mensaje de log diferenciado); resultado igual al cómputo desde cero"
      tol_rel: 1.0e-12
    - name: fallback_cholesky_lu
      setup: "matriz simétrica pero casi-singular o no PD efectiva (mecanismo cinemático añadido y luego soltado)"
      expect: "Cholesky aborta, LU resuelve sin fallar"
    - name: mpc_simple
      setup: "Quad4 con MPC u_2 = u_1 (master/slave)"
      expect: "u_2 = u_1 en la solución, equilibrio respetado"
      tol_abs: 1.0e-12

references:
  - "Bathe K.J. (2014). Finite Element Procedures. Prentice Hall. §2."
  - "Cook R.D., Malkus D.S., Plesha M.E., Witt R.J. (2002). Concepts and Applications of FEA. Wiley. §2.4 (eliminación de DOFs Dirichlet)."
  - "ADR 0003 — selección automática de backend algebraico."
  - "ADR 0004 — manejo de Dirichlet y constraints lineales."
```

---

## Implementación

- **Archivo**: [fenix/math/solvers/linear.py](../../fenix/math/solvers/linear.py).
- **Clase**: `LinearSolver`, registrada vía `@SolverRegistry.register` con `PIPELINE_KIND = "static"`.
- **Cache**: atributos `_factor`, `_T`, `_g_full`, `_F_dir`, `_n_free` se rellenan lazy en la primera llamada a `solve`. Método público `invalidate_cache()` los pone a `None` para forzar reensamblaje.
- **Fallback Cholesky→LU**: capturado en `_build_cache` vía `CholeskyNotPositiveDefiniteError` (ADR 0003 §5); registra un warning y reintenta con LU.
- **Entrypoint público**: `fenix.run_static(model, solver="linear", ...)` (despacho declarativo por `PIPELINE_KIND`, regla C de la auditoría).
- **Tests**: cobertura masiva implícita — prácticamente todos los tests del pipeline estático lineal lo usan. Tests específicos del cache: `tests/test_solver_robustness.py::test_linear_solver_cache_reuse`.

---

## Diálogo

- **2026-05-19** · Spec creada retroactivamente para cerrar el hueco H-5.3. Solver anterior a la convención de specs validadas. La cache de factorización se añadió en la sesión de saneamiento post-auditoría (commit `b517718`, H-4.2); esta spec recoge ese comportamiento como parte del contrato actual.
