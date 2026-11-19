# PrescribedCrackPath — trayectoria de grieta prescrita

> Spec corta de variante (Reglas.md §4): documenta sólo lo que cambia respecto de [`CrackTracking`](CrackTracking.md). Redactada e implementada por la IA el 2026-11-10 para los modelos con grieta prescrita de la tesis, como variante del seguimiento global de entonces (`GlobalCrackTracking`); el 2026-11-19 su padre pasó a ser `CrackTracking` (prolongación recta) sin cambio de comportamiento: la variante reproduce bit a bit su secuencia de activaciones anterior (astillas; banda de Van Vliet en las pruebas).

---

## Especificación física

La trayectoria de la grieta es un dato: una poligonal dada, prolongada por sus extremos. Es el modelo de **grieta prescrita** de Retama (2010) §8.1 (ensayo de Van Vliet): una banda de elementos embebidos a lo largo de la trayectoria conocida.

Lo que cambia respecto de `CrackTracking`:
- Los valores nodales son la **distancia con signo** a la poligonal, positiva a su izquierda, en lugar de los de la prolongación recta; los nodos de los elementos ya agrietados guardan esa misma distancia (con el mínimo $\delta$), no la distancia a la recta de su elemento.
- Hay una sola grieta. Todo elemento intacto cortado por la trayectoria es punta y se activa si $\sigma_I > \sigma_{t0}$: no hay propagación desde una punta ni nucleación, ni salida por el borde libre.

No cambian:
- el criterio de Rankine;
- la normal de Rankine orientada hacia el nodo aislado;
- $l_d = (A/h)\cos\psi$;
- la irreversibilidad;
- la regla del vértice: si el nodo aislado da $\cos\psi \le 0$, cambia de lado y la trayectoria pasa por ese vértice; si ya está fijado, cambia el nodo libre que deje $\cos\psi > 0$ (en esta variante el aislado puede estar libre, porque todo elemento cortado es punta).

Con una banda de un solo elemento de espesor, el seguimiento no es aplicable: la grieta sale de la banda en cuanto las direcciones principales giran (medido el 2026-11-10 con el seguimiento global de entonces). La variante reproduce exactamente la partición elemental anterior (arriba/abajo en la banda) y el ensayo de Van Vliet.

---

## Contrato de implementación

```yaml
name: PrescribedCrackPath
kind: crack_tracker
status: validated

interface:
  parent: CrackTracking
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
      setup: "acceptance.astillas de CrackTracking (5 mallas Delaunay con cos ψ < 0)"
      expect: "todas se separan, 0 conflictos, 0 avisos, 0 fallos del Newton local"

references:
  - "Retama Velasco J. (2010). Tesis doctoral, UNAM, §8.1 (Van Vliet, grieta prescrita)."
```

---

## Implementación

- Archivo: `solidum/user/discontinuities/crack_tracking.py`
- Clase: `PrescribedCrackPath(CrackTracking)` (registrada en `CrackTrackerRegistry`); distancia con signo en `polyline_signed_distance`.
- Redefine sólo los puntos de extensión del padre: `_refresh` (valores = distancia a la poligonal, con los guardados y los cambios de lado de la evaluación; puntas = todo intacto cortado), `_can_nucleate` (nunca), `_prepare_cracks` (una grieta), `_adoption_crack` (la grieta 0, orientada como la poligonal) y `_fixed_values` (la distancia a la poligonal).
- Uso: `tests/user/discontinuities/_van_vliet.py` (banda con `crack_tracker=PrescribedCrackPath([[-D, 0], [D, 0]])`), `test_disc_crack_tracking.py::TestAstillas`. YAML: `crack_trackers: [{id: 1, type: PrescribedCrackPath, points: [[x0, y0], [x1, y1]]}]`.

---

## Diálogo

- **2026-11-10** · Creada al implementar el seguimiento global: el ensayo de Van Vliet de la tesis (banda de un elemento de espesor, grieta prescrita) se bloqueaba con el seguimiento global. Antes lo resolvía la regla elemental, que con la banda daba la partición arriba/abajo; ahora la da la trayectoria prescrita, con el mismo resultado.
