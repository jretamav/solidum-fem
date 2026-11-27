# ADR 0038 — Energías y magnitudes de elemento en la historia de la salida

- **Estado**: aceptado (2026-11-11)
- **Fecha**: 2026-11-11
- **Alcance**:
  - `solidum/utils/output_writer.py`: columnas `energy` y `element` de `output.history`;
  - `solidum/core/material.py`, `solidum/core/interface_material.py`: contrato `SECANT_UNLOADING`, declarado por los materiales elásticos, de daño y las cohesivas de daño;
  - `solidum/entry.py`: el escritor recibe la carga aplicada sin convertir (con sus cargas equivalentes por elemento).

## Contexto

La historia de la respuesta (ADR 0036) registraba desplazamientos y reacciones. La validación de los modelos de concreto necesita además la **energía disipada**, cuya comprobación central en la banda de grieta es `D = G_f·A`, y el **estado de los elementos** a lo largo del análisis: el daño en la zona de fisuración, el esfuerzo en una sección, las fuerzas internas de un marco. Hasta ahora esas cantidades se calculaban aparte, con scripts sobre la curva carga–desplazamiento. Los programas generales las dan como salida de historia: las energías `ALLWK`, `ALLSE` y `ALLDMD` y las variables de elemento de Abaqus.

## Decisión

1. **Energías** (`energy: external_work | elastic | dissipated`), sólo en análisis estáticos sin etapas (los mismos en que se registra la reacción). *Ampliado el 2026-11-27 (ADR 0036, «Ampliación del 2026-11-27»): también por etapas si ninguna es en el tiempo; el trabajo externo se acumula de una etapa a otra.* En equilibrio, la fuerza que el exterior (cargas y apoyos) ejerce sobre la estructura es `F_int(U)`, así que:
   - trabajo externo `W_n = W_{n−1} + ½ (F_int,n + F_int,n−1)·(U_n − U_{n−1})` (trapecio entre pasos convergidos; incluye el trabajo de los desplazamientos impuestos);
   - energía elástica recuperable `E_n = ½ U_n·F_int,n`; *ampliado por el ADR 0043 §5 (2026-11-26)*: `E_n = Σ_e E_e`, con `½ u_e·f_e` salvo en los elementos que entregan su energía almacenada (`Element.stored_energy`: las interfaces, `Σ_g m_g·ψ_g` de su ley, que cubre la fricción y la holgura inicial); sin interfaces es la misma cifra;
   - energía disipada `D_n = W_n − E_n`.

   `E = ½ U·F_int` es la energía que se recupera al descargar hasta el origen sólo si **todos** los materiales descargan por la secante (`σ = C_s·ε`): elásticos, daño, cohesivas de daño. Lo declaran con `SECANT_UNLOADING = True`; si algún material no lo declara (plasticidad, concreto con deformación residual, secciones por fibras) o hay fuerzas internas en `U = 0` (cables pretensados), pedir `elastic` o `dissipated` es un error al crear el escritor, y `external_work` sigue disponible.
2. **Magnitudes de elemento** (`element: {var, ...}`), en cualquier análisis:
   - `var: stress` o `strain` con `component`: una componente Voigt del proyecto (`xx`, `yy`, `zz`, `xy`, `yz`, `xz`; `zz` del esfuerzo en 2D es el fuera del plano del material), `mises` (esfuerzo) o `max_principal`; por punto de integración;
   - una variable interna escalar del `STATE_SCHEMA` del material (`damage`, `kappa`, `alpha`...), por punto de integración;
   - una fuerza interna de un elemento estructural 1D (`N`, `V`, `M`, `Vy`, `Vz`, `T`, `My`, `Mz`) en el extremo `end: i | j`, con la convención de Reglas §5 y la carga equivalente de la carga distribuida descontada (escalada por el factor de carga; sólo en análisis estáticos si hay carga distribuida).

   Selector: `id`, `ids`, `group` (elementos con todos sus nodos en el grupo físico) o `where` (condición sobre el centro del elemento, con el compilador de condiciones del ADR 0037); sin selector, todos los elementos que tienen esa magnitud. `reduce`: `mean` (por omisión; media de los puntos de integración o de los elementos), `min` o `max`.
3. **Estado leído, no recalculado.** El esfuerzo y las variables internas se leen del estado confirmado del paso (el que usó el análisis), por arreglos en los elementos del ensamblaje por lotes; no se vuelve a evaluar el material, lo que con IMPL-EX daría otro esfuerzo. En el análisis lineal, que no confirma el estado, se leen del de prueba que deja el ensamblaje en `U`. La deformación se calcula como `B·u`.

## Alternativas descartadas

- **Energía elástica por punto de integración** (`½ σ:ε^e` con la deformación elástica de cada material). Exige que cada material exponga su deformación recuperable; para los materiales de la secante coincide con `½ U·F_int`. Se adoptará cuando una validación con plasticidad lo pida, extendiendo el contrato.
- **Recalcular el estado en el callback.** Con IMPL-EX y con el estado ya confirmado, el material extrapolaría de nuevo; el estado guardado es el que equilibró el paso.
- **Variables de elemento por número de punto de integración.** El orden de los puntos depende de la regla; la reducción (media, mínimo, máximo) sobre un selector cubre los casos de validación sin depender de él.

## Consecuencias

- La comprobación `D/(G_f·A)` de la banda de grieta sale de la salida estándar: en la barra de objetividad por el YAML (Quad4, Δ = 10 mm), `D/(G_f·A) = 0.998`.
- Un material nuevo que descargue por la secante debe declarar `SECANT_UNLOADING = True` para admitir las energías elástica y disipada; por omisión no las admite.
- Las energías en etapas y en dinámica (cinética, amortiguamiento) quedan para cuando se necesiten. *Las de las etapas estáticas, hechas el 2026-11-27 (ADR 0036); las de dinámica siguen pendientes.*

## Ampliación (2026-11-25): energías por zona

4. **`energy` con `elements: {id | ids | group | where}`**: las mismas tres energías restringidas a un conjunto de elementos, con la fuerza interna de cada elemento `f_e = Σ_g w_g·detJ_g·s·B_gᵀσ_g` (`s`, el espesor en 2D) evaluada con los esfuerzos del estado que equilibró el paso (§3): `W_n = W_{n−1} + ½ Σ_e (f_e,n + f_e,n−1)·(u_e,n − u_e,n−1)` (el trabajo que recibe el material de la zona, `∫∫σ:dε dV`; para todo el modelo, el trabajo externo), `E_n = ½ Σ_e u_e,n·f_e,n` y `D_n = W_n − E_n`. Las energías de las zonas de una partición suman las del modelo (prueba con espesor ≠ 1). Sólo sólidos isoparamétricos (cinemática compilada, fuerza interna `∫Bᵀσ`; los marcos y los elementos con estado propio son un error) y, para `elastic` y `dissipated`, materiales de la zona con `SECANT_UNLOADING`. *Ampliado por el ADR 0043 §5 (2026-11-26)*: también elementos de interfaz (`f_e = G_e·t_e` con `G_e = [m_g·(W·R·B)_gᵀ]`), cuya energía almacenada es la de su ley.

Separa la energía por **zona** (la banda de grieta frente al concreto bajo la carga), no por **mecanismo** en un mismo punto (agrietamiento frente a aplastamiento con d⁺ y d⁻ a la vez): esa separación exige definir la disipación de cada mecanismo en cada material, que es una decisión física (con ν ≠ 0 el modelo 2 no deriva de un potencial; Reglas §3) y queda para el usuario.
