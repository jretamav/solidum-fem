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

1. **Energías** (`energy: external_work | elastic | dissipated`), sólo en análisis estáticos sin etapas (los mismos en que se registra la reacción). En equilibrio, la fuerza que el exterior (cargas y apoyos) ejerce sobre la estructura es `F_int(U)`, así que:
   - trabajo externo `W_n = W_{n−1} + ½ (F_int,n + F_int,n−1)·(U_n − U_{n−1})` (trapecio entre pasos convergidos; incluye el trabajo de los desplazamientos impuestos);
   - energía elástica recuperable `E_n = ½ U_n·F_int,n`;
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
- Las energías en etapas y en dinámica (cinética, amortiguamiento) quedan para cuando se necesiten.
