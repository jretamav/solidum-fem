# ADR 0031 — Análisis lineales sobre el estado de una etapa

- **Estado**: aceptado (2026-11-09)
- **Fecha**: 2026-11-09
- **Alcance**:
  - `LinearizedStateMixin` (`continue_from`, `_state_stiffness`), común a `ModalSolver`, `HarmonicSolver` y `ResponseSpectrumSolver`;
  - `run_modal`, `run_harmonic` y `run_response_spectrum` con `previous=`;
  - etapas `stages` del YAML;
  - `ModalResult.eigenvalues` y `unstable_modes`.
- **Extiende**: los ADR 0026 (etapas estáticas) y 0029 (continuidad en el tiempo) a los análisis lineales.

## Contexto

El modal, el armónico y el espectral usaban siempre la rigidez elástica inicial $\mathbf K(\mathbf 0)$. Un modal después de la gravedad no veía el efecto P-Δ, y uno sobre una estructura empujada no veía el agrietamiento. Por eso las etapas del YAML rechazaban un análisis modal tras otra etapa. Además, un autovalor negativo se recortaba a cero sin aviso.

## Decisión (del usuario, 2026-11-09)

1. **Rigidez tangente del estado de la etapa.** Con `continue_from(previous)` (o `previous=`), el análisis linealiza alrededor del estado final de `previous` (un `SolveResult` o el último instante de un `TransientResult`) con la rigidez que dan los elementos en esos desplazamientos, con su estado interno y sin incremento adicional:
   - la rigidez geométrica de los elementos corrotacionales (P-Δ);
   - en plasticidad, la de descarga, que es elástica, como en el paso de perturbación lineal de Abaqus;
   - en daño y en el concreto agrietado, la dañada.

   Sin etapa, sigue siendo $\mathbf K(\mathbf 0)$. Es lo que hacen el modal pretensado de ANSYS y `eigen` de OpenSees después de la gravedad. En el armónico, el amortiguamiento de Rayleigh se forma con esa misma rigidez.
2. **No cambia el estado**: una etapa posterior continúa la última que sí lo cambió (una estática o una en el tiempo).
3. **Autovalores negativos**: no se recortan. Quedan con su signo en `eigenvalues`, con frecuencia y periodo `NaN`, `unstable_modes` y un aviso: el estado es inestable en esos modos. El espectral se rechaza con un error si hay modos inestables.
4. **El ensamblador de la etapa** guarda el estado interno de los elementos. `previous=` sin el `assembler` de la etapa, ni un `solver` que lo traiga, da un error en lugar de partir del estado virgen.

## Medido

`tests/test_linearized_stage.py`:
- **Columna esbelta sin masa con masa en la cabeza y P = 0.5·P_cr** (corrotacional): la frecuencia coincide con la del oscilador de rigidez lateral exacta $k(P) = P\kappa/(\tan\kappa L - \kappa L)$ a 4e-4 con 20 elementos. El error baja de forma cuadrática con la malla, hacia un resto del orden de $P/EA$: la fórmula clásica supone la columna inextensible. Con una columna robusta ($P/EA \approx 10^{-3}$), el acortamiento axial cambia $k$ en ese mismo orden.
- **Armónico y espectral**: lo mismo con esa rigidez.
- **P = 1.2·P_cr**: autovalor negativo, que antes daba frecuencia 0.
- **Columna de concreto reforzado** (N = 600 kN) empujada a 20 mm: $T_1$ pasa de 0.269 s a 0.441 s.
- **Acero que fluyó**: el periodo no cambia (descarga elástica), y empujar después del modal da lo mismo que sin él.

## Fuera de este ADR

- El `LinearSolver` como etapa (una perturbación estática con la tangente del estado). Falta decidir si la etapa lineal cambia el estado interno de los elementos para la etapa siguiente.
- Un estático después de un transitorio.

## Referencias

- Abaqus Analysis User's Guide, §6.1.3 (*General and linear perturbation procedures*).
- Bathe, K.-J. (2014). *Finite Element Procedures*, §6.8 (estabilidad lineal) y §10.3.
