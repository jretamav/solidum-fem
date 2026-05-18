# Auditoría de robustez de solvers no lineales — fase A

**Fecha**: 2026-05-18.
**Tests**: [`tests/test_solver_robustness.py`](../../tests/test_solver_robustness.py) — 11 tests, 11 verdes.

---

## Objetivo y alcance

Mapear empíricamente el comportamiento de `NonlinearSolver`, `ArcLengthSolver` y
`NewtonNewmarkSolver` en regímenes que ejercitan sus límites de robustez:
puntos límite, snap-back, bifurcación, plasticidad perfecta, dt sobredimensionado.
Esta fase **no modifica los solvers**. Documenta dónde sostienen y dónde se rompen,
para informar el mini-ADR de la fase B.

Criterios de marcado de tests:

- **`assertEqual/assertLess/...`** verifica comportamiento esperable; tests verdes
  documentan robustez ya presente.
- **`assertRaises(RuntimeError)`** documenta un modo de fallo limpio del solver
  actual (informativo, no aspiracional). Cuando la fase C introduzca la mejora
  correspondiente, el test se reformula como caso de éxito.

---

## Resumen ejecutivo

**Resultado neto**: los tres solvers son más robustos de lo que sugerían los huecos
de validación. Los modos de fallo encontrados están concentrados en un solo patrón
arquitectural (**ausencia de globalización del Newton**); el resto de regímenes
problemáticos predichos no se reproducen en el código actual.

| Régimen | Solver | Resultado |
|---|---|---|
| Snap-through (von Mises 2-bar) | `ArcLengthSolver` | ✅ Traza la curva (24 pasos hasta λ=1, postcrítico u_y ≈ −2h) |
| Snap-through (mismo problema) | `NonlinearSolver` | ✅ **Hallazgo positivo**: el adaptativo sostiene el snap-through (20 pasos con bisección de Δλ) |
| Snap-back con daño 1D continuo | `ArcLengthSolver` | ✅ **Hallazgo positivo**: traza la rama post-pico con softening exponencial fuerte (α=500, 1 paso, 4 iter) |
| Bifurcación con imperfección (1% L) | `ArcLengthSolver` | ✅ Llega a λ ≥ 0.9·P_cr; flecha lateral amplifica como esperado |
| Control de `dl` agresivo (10× razonable) | `ArcLengthSolver` | ✅ Auto-shrink reacciona y traza el snap-through completo |
| Drucker-Prager con endurecimiento suave | `NonlinearSolver` | ✅ Plastificación moderada, 10 pasos a λ=1 |
| Drucker-Prager perfectamente plástico + sobrecarga | `NonlinearSolver` | ⚠ **Hallazgo crítico**: Newton oscila entre dos estados; bisección de Δλ no rescata; `RuntimeError` con mensaje genérico |
| Daño 2D plane stress con softening abrupto | `NonlinearSolver` | ✅ Diverge limpiamente con `RuntimeError` (sin oscilación) |
| Reversal de carga con plasticidad J2 1D | `NonlinearSolver` | ✅ Commit/rollback de history bajo bisección sin contaminación |
| dt sobredimensionado (T/3) + plasticidad perfecta | `NewtonNewmarkSolver` | ✅ **Hallazgo positivo**: la masa estabiliza el jacobiano; converge en 1-2 iter/paso |
| Daño que se activa mid-transitorio | `NewtonNewmarkSolver` | ✅ Rayleigh-K₀ no introduce divergencias en degradación de K_t |

---

## Hallazgos por solver

### `NonlinearSolver` (Newton-Raphson con paso adaptativo)

**Hallazgo principal — robustez infravalorada**. El solver sostiene snap-through
con control de carga puro cuando se le permite bisecar el incremento. En el von
Mises 2-bar con `h/L = 0.1` y carga 10× el pico estimado, completó 20 pasos
hasta λ=1 con `u_y` postcrítico ≈ −2h. Esto **contradice** la afirmación común
(reflejada en [`docs/catalogo_solvers.md`](../catalogo_solvers.md) §`NonlinearSolver` →
*Limitación*) de que "no puede atravesar puntos límite". La bisección agresiva,
junto con el cambio cinemático de la formulación corotacional, basta en este
caso. El arc-length sigue siendo necesario en otros regímenes (snap-back puro
con `du/dλ < 0`), pero la frontera de aplicabilidad de Newton+bisección es más
ancha que la documentada.

**Hallazgo crítico — oscilación de Newton sin line search.** En Drucker-Prager
plane strain con `H=0` (perfectamente plástico) y carga ~50× la cohesión
característica, el Newton **oscila entre dos pozos** con ratios alternados
2.1 ↔ 440 durante las 30 iteraciones permitidas. El paso se biseca, vuelve a
oscilar, repite, y muere al alcanzar `min_delta_lambda` con mensaje genérico
*"el solver ha divergido"*.

Patrón observado en el log de iteraciones:

```
Iteración 28 | R/tol_F: 1.03e+02 | dU/tol_d: 1.74e+02
Iteración 29 | R/tol_F: 4.97e-01 | dU/tol_d: 1.74e+02   ← casi convergido
Iteración 30 | R/tol_F: 1.03e+02 | dU/tol_d: 1.74e+02   ← rebote
```

Este es exactamente el patrón que motiva la globalización del Newton-Raphson
mediante **line search** (búsqueda lineal con Armijo o suficiente decrecimiento):
el incremento completo de Newton aterriza fuera del pozo de convergencia y
rebota. Con `α < 1` óptimo (encontrado por backtracking), el residuo descendería
monótonamente.

**Diagnóstico secundario — falta de telemetría útil**. El mensaje de fallo no
distingue entre:

1. Singularidad real de `K_t` (problema físico mal planteado).
2. Newton fuera del pozo de convergencia (rescatable con line search).
3. Carga > capacidad plástica del modelo (problema sin solución estática).
4. Tolerancia mal calibrada para la escala del problema.

Un futuro `RuntimeError` debería incluir al menos: último residuo y delta, patrón
de oscilación detectado (sí/no), número de bisecciones consumidas, último δλ.

### `ArcLengthSolver` (Crisfield cilíndrico)

**Hallazgo principal — sostén de softening de daño bulk**. STATUS.md menciona
*"NonlinearSolver y ArcLengthSolver cilíndrico no atraviesan la transición
elástico→softening con penalty cohesivo stiff"* como limitación. Esta
auditoría confirma que el límite es específico al **penalty cohesivo stiff** de
`CST_Embedded2D`: el arc-length sí traza la rama post-pico de
`IsotropicDamage1D` con softening exponencial fuerte (α=500, κ_0=1e-3) con
carga 2× el pico, convergiendo en 4 iteraciones en un único paso. El problema
del embedded no es el solver; es la mala condición de la rigidez del cono
elástico-cohesivo cerca de la activación, que será objeto del mini-ADR de
"solvers para softening severo" (deuda #4 de STATUS.md).

**Bifurcación con imperfección**: el solver converge en 2-3 iter por paso,
pero los factores de carga por paso en los primeros pasos son inestables
(λ salta entre valores cercanos a 0 e incluso negativos durante el post-pandeo
inicial). Es comportamiento esperado en columna esbelta cerca de carga crítica
con auto-shrink del `dl`; no es un bug, pero conviene documentarlo para no
sorprenderse al ver `λ` retroceder transitoriamente. La flecha lateral final
sí amplifica como predice el modelo de imperfección.

**Sin tracking de signo de pivote**: el método `_negative_pivots` retorna `None`
siempre. En la geometría con imperfección (rama única bien definida) no es
necesario; en bifurcación pura sería. No se ejercitó en esta auditoría
(decisión deliberada de incluir solo imperfección, ver discusión previa).

### `NewtonNewmarkSolver` (Newmark + Newton dinámico no lineal)

**Hallazgo principal — la masa estabiliza el jacobiano**. El régimen predicho
como problemático (dt sobredimensionado con plasticidad perfecta) **no se
manifiesta** porque el jacobiano dinámico `J = M + γΔt·C + βΔt²·K_t` mantiene
la diagonal dominante de `M` aunque `K_t → 0` en plasticidad perfecta. El test
con `dt = T/3` (un periodo natural en 3 puntos de muestreo) y plasticidad J2
perfecta convergió en 1-2 iteraciones por paso.

Consecuencia para la fase B: la mejora "paso temporal adaptativo en
`NewtonNewmarkSolver`" tiene **prioridad menor** de la que parecía. El régimen
que la motivaría (oscilación del Newton con dt grande) está cubierto por la
masa. La regla asimétrica entre `NonlinearSolver` (biseca) y `NewtonNewmark`
(falla con `RuntimeError`) sigue presente, pero el escenario donde realmente
hace daño es estrecho.

**Daño mid-transitorio**: el Rayleigh calibrado con K₀ no introduce
divergencias cuando K_t se degrada durante el transitorio. Validación
positiva de la decisión documentada en
[`docs/specs/NewtonNewmarkSolver.md`](../specs/NewtonNewmarkSolver.md).

---

## Implicaciones para la fase B (mini-ADR de estrategias de robustez)

Reordenamiento de prioridades respecto a la propuesta original:

1. **Line search Armijo en `NonlinearSolver` y `NewtonNewmarkSolver`** — **alta**.
   Es la mejora que más expansiona el régimen de aplicabilidad: rescata el modo
   de oscilación de Newton documentado en Drucker-Prager y previsiblemente
   cualquier otro régimen plástico con tangente marginal. Trabajo acotado:
   backtracking sobre `‖R(U + α·δU)‖` con condición de suficiente decrecimiento.

2. **Telemetría de divergencia útil** — **alta**.
   Antes de tocar el algoritmo, distinguir en el `RuntimeError` los modos de
   fallo (oscilación, singularidad, mal escalado, sobrecarga física). Es trabajo
   barato y multiplica el valor diagnóstico de cualquier divergencia futura.

3. **Tracking de signo de pivote en `ArcLengthSolver`** — **media-baja**.
   Útil solo en bifurcación pura sin imperfección; en problemas físicos reales
   con imperfecciones (las habituales), el predictor por proyección basta.
   Postergar hasta que aparezca un caso real con bifurcación pura.

4. **Paso temporal adaptativo en `NewtonNewmarkSolver`** — **baja**.
   El régimen que lo motivaba (oscilación con dt grande) no se reproduce.
   La asimetría con `NonlinearSolver` sigue siendo conceptualmente
   incoherente, pero el riesgo práctico es bajo.

5. **Actualizar [`docs/catalogo_solvers.md`](../catalogo_solvers.md)** —
   inmediato. La línea *"NonlinearSolver no puede atravesar puntos límite"*
   debería matizarse: "puede atravesar puntos límite suaves con bisección
   adaptativa si la cinemática del elemento es no lineal (corotacional); no
   atraviesa snap-back con `du/dλ < 0`".

El mini-ADR de la fase B debería formalizar las dos primeras decisiones, dejando
las otras tres como notas para revisión cuando emerjan los disparadores.

---

## Limitaciones de esta auditoría

- **Cobertura no exhaustiva**. Once tests no agotan el espacio de regímenes; en
  particular, no se ejercitaron: combinación arc-length + materiales no
  simétricos (daño 2D con tangente consistente asimétrica), Drucker-Prager
  con tracking de paso por el ápice (`f_trial > 0` con `s_trial → 0`),
  problemas grandes (>1000 DOFs) donde la condición de `K_t` empieza a
  importar.

- **Sin medición de coste**. No se cuentan factorizaciones LU ni tiempos
  acumulados. Para optimización de rendimiento (no objetivo de esta fase),
  haría falta otra pasada con `time` y profiler.

- **Inferencias sobre line search**. La aserción de que line search rescataría
  el modo de oscilación de Drucker-Prager es plausible pero no probada;
  validarla es objeto de la fase C, no de esta.

---

*Última actualización: 2026-05-18 — fase A cerrada con 11 tests verdes y modos
de fallo documentados. Próximo paso: fase B (mini-ADR de estrategias de
robustez de Newton).*
