# ADR 0006 — Tolerancias del criterio de admisibilidad constitutiva

- **Estado**: aceptado
- **Fecha**: 2026-05-12
- **Alcance**: clase base `Material`, todos los modelos constitutivos con check de admisibilidad (`Elastoplastic1D`, `IsotropicDamage1D`, `IsotropicDamage2D`, `VonMises2D`).

## Contexto

Cada modelo constitutivo no lineal tiene un criterio de admisibilidad — una función `f(estado)` que decide si el material está en régimen reversible (`f ≤ 0`) o disipativo (`f > 0`). En aritmética de doble precisión la comparación no puede ser exacta: tras un return mapping, el ruido de redondeo deja `f` cerca de cero pero no exactamente cero, por lo que la comparación se escribe `f ≤ tol`.

La implementación anterior usaba `PLASTIC_YIELD_TOL = 1e-9` (constante absoluta, sin unidades). Eso es frágil por dos razones:

1. **No invariante bajo cambio de unidades.** `f` tiene unidades de esfuerzo. El ruido de redondeo en `f` escala como `σ_y · ε_machine`. Si el usuario trabaja en Pa con `σ_y ~ 1e8`, el ruido es ~`1e-8` — por encima de la tolerancia, y el código entra espuriamente en la rama plástica. Si trabaja con `σ_y ~ 1e-6` (adimensional fuerte), una violación real del 0.1% (`f ~ 1e-9`) cae dentro de la tolerancia y se enmascara.
2. **No adaptativo al estado.** Modelos con endurecimiento severo (Voce, Chaboche) o que carecen de `σ_y` (daño puro, cinemático puro) no encajan en una tolerancia hardcodeada a un parámetro inicial fijo.

## Decisión

Se adopta el esquema combinado absoluto + relativo, estilo ODE solver (`scipy.integrate`, SUNDIALS, Code_Aster):

```
f ≤ ADMISSIBILITY_TOL_ABS + ADMISSIBILITY_TOL_REL · escala(estado)
```

con:
- `ADMISSIBILITY_TOL_REL = 1e-10` — número puro, banda de seis órdenes sobre `ε_machine ≈ 1e-16` en doble.
- `ADMISSIBILITY_TOL_ABS = 1e-14` — esfuerzo, piso para estados donde `escala → 0` (estado virgen, daño antes del umbral si la escala se evaluara contra `κ` corriente, modelos cinemáticos puros).
- `escala(estado)` la declara cada material vía `Material.admissibility_scale(state_vars)`, devolviendo un esfuerzo positivo característico del criterio en el estado corriente.

## Contrato `admissibility_scale`

Método nuevo de la clase base `Material`:

```python
def admissibility_scale(self, state_vars=None) -> float:
    """Escala característica del criterio de admisibilidad, en esfuerzo."""
    return 1.0   # default; solo materiales con check lo sobreescriben
```

Implementaciones actuales:

| Material | `admissibility_scale` |
|---|---|
| `Elastic1D`, `Elastic2D` | default 1.0 (no se invoca, no tienen check) |
| `CableMaterial1D` | default 1.0 (ramificación directa por signo de ε) |
| `Elastoplastic1D` | `σ_y + H·α` (fluencia corriente, adaptativa) |
| `IsotropicDamage1D` | `E·κ_0` (esfuerzo umbral inicial) |
| `IsotropicDamage2D` | `E·κ_0` (idem) |
| `VonMises2D` | `√(2/3)·(σ_y + H·α)` (norma desviadora de la superficie corriente) |

La comparación se encapsula en `Material.is_admissible(f, state_vars)` para que los modelos no la repliquen.

## Consecuencias

**Inmediatas**:
- Misma física en MPa, Pa, GPa o adimensional produce **idéntico** estado interno hasta precisión de doble. Cubierto por test `TestAdmissibilityToleranceUnitInvariance`.
- En endurecimiento severo la tolerancia crece junto con la fluencia, evitando que el ruido relativo se vuelva una fracción dominante de la escala.
- En el caso degenerado `escala → 0` (no actualmente realizable con los modelos del repo, sí en extensiones futuras como cinemático puro o phase-field con `g_c` muy pequeña), el piso absoluto `1e-14` evita comparar contra cero exacto.

**Hacia adelante**:
- Modelos energéticos (phase-field, gradient damage, cohesive) encajan sin modificar el contrato: devuelven `admissibility_scale = √(2·E·g_c)` o equivalente, derivado de su densidad de energía crítica.
- Modelos con return mapping no lineal (Voce, Chaboche) pueden migrar a semi-smooth Newton interno reutilizando `admissibility_scale` como referencia de normalización del residuo del Newton.
- `Material.is_admissible` es el único punto del código que aplica la fórmula; cambiar la política (p. ej. a piso absoluto adaptativo `max(ABS, ε_machine · escala)`) es una edición en un sitio.

**Migración**:
- Constante `PLASTIC_YIELD_TOL` retirada de `solidum/constants.py` y sus importadores. No tenía consumidores externos al proyecto.
- Tests existentes pasan sin modificación: el rango de unidades habitual (MPa) cae dentro del rango donde la nueva política coincide en valor numérico con la antigua hasta el límite de doble precisión.

## Alternativas consideradas

- **Mantener la tolerancia absoluta** y aceptar la fragilidad. Rechazada: no cumple el estándar de robustez que comparten Abaqus, ANSYS, Code_Aster, LS-DYNA.
- **Relativa a parámetro fijo `σ_y` inicial**. Rechazada: no adaptativa al endurecimiento, no extensible a modelos sin `σ_y`.
- **Relativa a la rigidez** (`E · ε_ref` estilo OpenSees). Descartada: la rigidez no siempre es la escala correcta del criterio (en plasticidad cinemática pura la rigidez no se mueve mientras la back-stress sí lo hace).
- **Formulación variacional/energética pura** (estilo Carstensen-Mielke). Más rigurosa, pero implica reescribir todos los modelos como minimización incremental. Coste desproporcionado para el beneficio inmediato; queda como camino futuro para modelos específicos.
- **Semi-smooth Newton interno** con la tolerancia viviendo en el Newton. Compatible con el contrato propuesto y se reservaría para cuando los modelos lo demanden (Voce, Chaboche).

## Deuda conocida

Identificada al cerrar el ADR; no bloqueante con los cuatro materiales actuales, a revisitar cuando el catálogo crezca.

- **Contrato por convención, sin enforcement.** Un material futuro puede tener un criterio de admisibilidad no trivial y, por descuido, no invocar `is_admissible`/`admissibility_tol` — replicando su propia comparación con tolerancia hardcoded y saltándose la política central. La clase base no tiene un punto natural de validación que detecte esto en construcción (a diferencia de `STRAIN_DIM`, que sí lo valida el elemento). Mitigación pendiente: cuando haya 6-8 materiales se justifica un mecanismo más estricto (p. ej. flag declarativo `HAS_ADMISSIBILITY_CHECK` que obligue a sobreescribir `admissibility_scale`).

- **Footgun del default `admissibility_scale = 1.0`.** Si un material con check olvida sobreescribir la escala, hereda `1.0` y la tolerancia resultante (`~1e-10` en unidades absolutas) es absurdamente apretada en cualquier sistema físico real. El fallo es silencioso, no ruidoso: produce un régimen "todo plástico/dañado siempre" sin excepción. La mitigación natural es la misma del punto anterior (declarar explícitamente la presencia del check para forzar el override).

- **Sin test del piso `ABS`.** El caso degenerado `escala → 0` está cubierto por el término absoluto `1e-14`, pero ningún material actual lo ejercita. Cuando entre el primero (cinemático puro, phase-field con `g_c` muy pequeña, etc.), escribir el test antes de la implementación.

## Referencias

- de Souza Neto, Perić, Owen (2008), *Computational Methods for Plasticity*, §7.3 (tolerancias de retorno).
- Simo & Hughes (1998), *Computational Inelasticity*, §2.7.
- Code_Aster, documentación de `RESI_INTE_RELA` / `RESI_INTE_ABSO`.
- SUNDIALS, *CVODES User Guide*, sección sobre tolerancias mixtas atol+rtol.
