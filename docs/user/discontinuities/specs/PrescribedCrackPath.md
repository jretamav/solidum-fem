# PrescribedCrackPath — trayectoria de grieta prescrita

> Spec corta de variante (Reglas.md §4): documenta sólo lo que cambia respecto de [`GlobalCrackTracking`](GlobalCrackTracking.md). Redactada e implementada por la IA el 2026-11-10 para los modelos con grieta prescrita de la tesis.

---

## Especificación física

La trayectoria de la grieta es un dato: una poligonal dada, prolongada por sus extremos. Es el modelo de **grieta prescrita** de Retama (2010) §8.1 (ensayo de Van Vliet): una banda de elementos embebidos a lo largo de la trayectoria conocida.

Lo que cambia respecto del seguimiento global:
- $\theta$ es la **distancia con signo** a la poligonal, positiva a su izquierda, en lugar de la solución del problema auxiliar.
- Hay una sola grieta. Se activa todo elemento cortado por la trayectoria con $\sigma_I > \sigma_{t0}$: no hay propagación desde una punta ni nucleación.

No cambian:
- el criterio de Rankine;
- la normal de Rankine orientada hacia el nodo aislado;
- $l_d = (A/h)\cos\psi$;
- la irreversibilidad;
- la regla del vértice: si el nodo aislado da $\cos\psi \le 0$, cambia de lado y la trayectoria pasa por ese vértice; si ya está fijado, cambia el nodo libre.

Con una banda de un solo elemento de espesor, el seguimiento global no es aplicable: la isolínea sale de la banda en cuanto las direcciones principales giran. La variante reproduce exactamente la partición elemental anterior (arriba/abajo en la banda) y el ensayo de Van Vliet.

---

## Contrato de implementación

```yaml
name: PrescribedCrackPath
kind: crack_tracker
status: validated

interface:
  parent: GlobalCrackTracking
  region: "los CST_Embedded2D que lo referencian"

parameters:
  - { name: points, type: list, required: true, desc: "Vértices [x, y] de la poligonal (al menos dos), prolongada por sus extremos" }

validity:
  - trayectoria conocida (grieta prescrita)
  - poligonal sin quiebres bruscos (el signo es el del tramo más cercano)

out_of_scope:
  - varias grietas

acceptance:
  specific:
    - name: van_vliet_grieta_prescrita
      setup: "ensayo de Van Vliet de la tesis (§8.1) con la trayectoria y = 0"
      expect: "pico, deformación en el pico y carga a 20, 50 y 100 µm dentro de la lectura de la fig. 8.3 (las mismas tolerancias que antes del seguimiento)"
    - name: astillas
      setup: "acceptance.astillas de GlobalCrackTracking (5 mallas Delaunay con cos ψ < 0)"
      expect: "todas se separan, 0 conflictos, 0 avisos, 0 fallos del Newton local"

references:
  - "Retama Velasco J. (2010). Tesis doctoral, UNAM, §8.1 (Van Vliet, grieta prescrita)."
```

---

## Implementación

- Archivo: `solidum/user/discontinuities/crack_tracking.py`
- Clase: `PrescribedCrackPath(GlobalCrackTracking)` (registrada en `CrackTrackerRegistry`); distancia con signo en `polyline_signed_distance`.
- Redefine sólo los ganchos del padre: `_solve` (θ = distancia, con los valores fijados y los cambios de lado del paso), `_update_chain` (todo intacto cortado es punta), `_can_nucleate` (nunca), `_prepare_cracks` (una grieta), `_assemble_field` (nada que ensamblar), `_fixed_values`.
- Uso: `tests/user/discontinuities/_van_vliet.py` (banda con `crack_tracker=PrescribedCrackPath([[-D, 0], [D, 0]])`), `test_disc_crack_tracking.py::TestAstillas`. YAML: `crack_trackers: [{id: 1, type: PrescribedCrackPath, points: [[x0, y0], [x1, y1]]}]`.

---

## Diálogo

- **2026-11-10** · Creada al implementar el seguimiento global: el ensayo de Van Vliet de la tesis (banda de un elemento de espesor, grieta prescrita) se bloqueaba con el seguimiento global. Antes lo resolvía la regla elemental, que con la banda daba la partición arriba/abajo; ahora la da la trayectoria prescrita, con el mismo resultado.
