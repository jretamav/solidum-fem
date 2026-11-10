# ADR 0030 — Reacciones de los apoyos en los análisis dinámicos

- **Estado**: aceptado (2026-11-09)
- **Fecha**: 2026-11-09
- **Alcance**:
  - `TransientResult`, `HarmonicResult` y `ResponseSpectrumResult`: campos `reaction_labels`, `R_history` / `R_damping_history`, `R_complex` / `R_damping_complex`, `R_per_mode` / `R_combined`, y el método `total_reaction`;
  - `TransientTrace` (registro por paso, común a los cinco integradores mecánicos);
  - `solidum.results.reaction_dofs`, `reactions_from_residual`, `reaction_labels`; `solidum.math.modal_response.combine_modal`.
- **Extiende**: el ADR 0004 §3 (reacciones y fuerzas de restricción del estático) a los análisis dinámicos.

## Contexto

Los análisis dinámicos no daban reacciones: sin ellas no hay historia del cortante basal, la magnitud con que se revisa un análisis sísmico. En un apoyo, la ecuación de movimiento tiene tres partes: la fuerza de los elementos, la inercia de la masa del apoyo o acoplada a él (masa consistente) y el amortiguamiento. Con amortiguamiento de Rayleigh, `C = a₀·M + a₁·K`, la parte `a₁·K·u̇` es la fuerza de amortiguamiento que pasa por los elementos de la base. Qué partes cuenta la reacción es una decisión física. OpenSees deja elegir: `reactions -dynamic` las cuenta todas, y `-rayleigh` quita la de Rayleigh.

## Decisión (del usuario, 2026-11-09)

1. **La reacción es la ecuación de movimiento en las filas de apoyo**, con todas sus partes:
   $$\mathbf R = \mathbf F_\text{int}(\mathbf u) + \mathbf M\ddot{\mathbf u} + \mathbf C\dot{\mathbf u} - \mathbf F .$$
   Cumple el equilibrio de la estructura entera: $\sum_x \mathbf R = \sum_x(\mathbf M\ddot{\mathbf u} + \mathbf C\dot{\mathbf u}) - \sum_x\mathbf F$, porque $\sum_x\mathbf F_\text{int} = 0$ por traslación rígida.
2. **La parte de amortiguamiento, $\mathbf C\dot{\mathbf u}$ en los apoyos, se reporta por separado**, para obtener la reacción sin Rayleigh al posprocesar (`damping=False`), sin volver a correr el análisis.
3. **Armónico**: la misma regla en amplitudes complejas, $\hat{\mathbf R} = (\mathbf K - \omega^2\mathbf M + i\omega\mathbf C)\hat{\mathbf u} - \hat{\mathbf F}$, con su parte $i\omega\mathbf C\hat{\mathbf u}$.
4. **Espectral**: la reacción de cada modo es la ecuación del modo en los apoyos, $(\mathbf K - \omega_n^2\mathbf M)\mathbf u_n$ (el amortiguamiento está en el espectro), y se combina como cualquier otra respuesta (SRSS o CQC). Una magnitud que suma reacciones, como el cortante basal, se suma en cada modo antes de combinar (`total_reaction`).

Se guardan sólo los grados de libertad de apoyo (`reaction_labels`, `(nodo, dof)`): los Dirichlet y los de las restricciones multipunto de apoyo, con el mismo reparto entre apoyos y uniones que el estático (ADR 0004 §3). En el transitorio, la reacción se evalúa en cada instante guardado con el estado de ese instante. En HHT-α, cuya ecuación se cumple en $t_{n+1+\alpha}$, el equilibrio global en $t_{n+1}$ se cumple a $O(\alpha\,\Delta t)$, no exactamente. En la diferencia central con amortiguamiento ocurre lo mismo, porque usa $\dot{\mathbf u}_{n+1/2}$.

## Medido

En un voladizo de acero con masa en la cabeza, Rayleigh al 5 % y un pulso lateral, quitar el amortiguamiento cambia la historia del cortante basal en un 0.8 % del pico. La inercia en el apoyo aporta un 0.07 % con masa consistente y nada con masa concentrada. En análisis inelásticos con Rayleigh proporcional a la rigidez inicial la parte de amortiguamiento puede ser mucho mayor (Charney, 2008): por eso se reporta aparte.

Pruebas (`tests/test_dynamic_reactions.py`):
- equilibrio global a 1e-10 con Newmark amortiguado y con la diferencia central sin amortiguar;
- cortante basal del edificio de cortante igual a $-k_1(u_1 + a_1\dot u_1)$ a 1e-9;
- armónico contra el oscilador con amortiguamiento de rigidez a 1e-10;
- cortante basal espectral de cada modo igual a $-M_n^*S_a$ a 1e-9.

## Consecuencias

- El cortante basal es `result.total_reaction("ux")`, y sin Rayleigh, `total_reaction("ux", damping=False)`.
- El costo es un producto matriz-vector por matriz y por paso guardado, y un vector por grado de libertad de apoyo.

## Referencias

- Chopra, A. K. (2017). *Dynamics of Structures*, §9.4 y cap. 13 (fuerzas de apoyo y respuesta espectral).
- Charney, F. A. (2008). Unintended consequences of modeling damping in structures. *Journal of Structural Engineering*, 134(4), 581–592.
