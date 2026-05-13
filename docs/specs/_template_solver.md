# <Nombre> — <descripción corta>

> Orden de trabajo. El usuario escribe la **especificación física**, la **formulación numérica** y el **contrato**. La IA rellena **implementación** y responde en **diálogo** durante el trabajo.

---

## Especificación física

### 0. Descripción general
*Tipo de análisis (estático lineal, estático no lineal, modal, transitorio lineal, transitorio no lineal, longitud de arco, etc.). Hipótesis cinemáticas globales (pequeñas deformaciones, corotacional…) y material (lineal, no lineal, con historia).*

### 1. Ecuación de equilibrio resuelta
*Forma fuerte semidiscreta del problema. Estática: $\mathbf K\,\mathbf u = \mathbf F$ o $\mathbf F_{\text{int}}(\mathbf u) = \mathbf F_{\text{ext}}$. Dinámica: $\mathbf M\ddot{\mathbf u} + \mathbf C\dot{\mathbf u} + \mathbf K\mathbf u = \mathbf F(t)$. Autovalores: $\mathbf K\boldsymbol\phi = \omega^2\mathbf M\boldsymbol\phi$. Longitud de arco: $\mathbf F_{\text{int}}(\mathbf u) = \lambda\,\mathbf F_{\text{ext}}^{\text{ref}}$ + restricción cuadrática.*

### 2. Condiciones iniciales / de contorno
*Apoyos Dirichlet (constantes en el tiempo si aplica), cargas externas, condiciones iniciales $\mathbf u(0)$, $\dot{\mathbf u}(0)$ para problemas transitorios. Comportamiento esperado en presencia de restricciones lineales multipunto (MPC).*

### 3. Salidas físicas
*Magnitudes que el usuario espera recuperar: campo $\mathbf u$, $\dot{\mathbf u}$, $\ddot{\mathbf u}$, esfuerzos en puntos de Gauss, factor de carga $\lambda$, frecuencias propias $\omega_i$, modos $\boldsymbol\phi_i$, historia temporal completa o solo último paso, diagnósticos de convergencia.*

---

## Formulación numérica

### 4. Esquema operativo
*Algoritmo paso a paso del solver. Estático lineal: una resolución $\mathbf K^{-1}\mathbf F$. Estático no lineal: bucle externo de incrementos $\lambda$, bucle interno Newton-Raphson. Transitorio: integrador temporal (Newmark-β, HHT-α, central differences) con predictor + corrector. Modal: shift-invert + ARPACK Lanczos.*

### 5. Predictor / corrector (si aplica)
*Fórmulas explícitas del predictor y del corrector. Para integradores temporales: $\tilde{\mathbf u}_{n+1}$, $\tilde{\dot{\mathbf u}}_{n+1}$ y sistema efectivo a resolver. Para arc-length: predictor tangente y restricción cuadrática.*

### 6. Criterio de convergencia
*Si el solver es iterativo: norma del residuo, norma del incremento, criterio dual fuerza+desplazamiento. Patrón canónico del proyecto (ADR 0007):*

$$\lVert\mathbf R[\text{libres}]\rVert \le \texttt{atol\_force} + \texttt{rtol\_force} \cdot \max(\lVert\mathbf F_{\text{ext}}\rVert, \lVert\mathbf F_{\text{int}}\rVert)$$
$$\lVert\Delta\mathbf U\rVert \le \texttt{atol\_disp} + \texttt{rtol\_disp} \cdot \lVert\mathbf U\rVert$$

*Semántica AND (ambos simultáneamente). `atol` autoderivada de escala característica para invariancia bajo cambio de unidades.*

### 7. Imposición de Dirichlet y MPC
*Mecanismo de imposición: eliminación directa (ADR 0004 fase 1) para apoyos rígidos; transformación $\mathbf T$ + offset $\mathbf g$ para restricciones lineales multipunto (ADR 0004 fase 2). Proyección a DOFs libres del sistema efectivo del solver.*

### 8. Backend algebraico
*Tipo de operación lineal requerida (factorización + resolución, eigsolver, solver iterativo). Reutilización de factorización entre pasos cuando proceda (`FactorizedSolver`, ADR 0003 fase 2). Selección automática del backend según propiedades del operador (ADR 0003).*

### 9. Adaptatividad / control de paso
*Si aplica: cómo se ajusta $\Delta\lambda$, $\Delta t$ o $dl$ según convergencia. Política de bisección al fallar. Cota inferior `min_delta_*` que provoca aborto.*

### 10. Caveats numéricos
*Limitaciones del esquema: puntos límite (Newton-Raphson no atraviesa snap-through), estabilidad condicional (central differences), oscilaciones espurias en transiciones (cable destensado), divergencia cerca de bifurcaciones, etc.*

---

## Contrato de implementación

```yaml
name: <Nombre>
kind: solver
status: draft            # draft → implemented → validated

interface:
  yaml_type: <Nombre>             # cómo se referencia desde el YAML del usuario
  output: ""                       # tipo de resultado: SolveResult, ModalResult, TransientResult, …

parameters:
  - { name: , type: , required: true, desc:  }

requirements:
  - ""                              # precondiciones sobre el modelo (density > 0 si dinámico, K constante si lineal, …)

conventions:
  units:    "heredadas del modelo (ADR 0008)."
  stability: ""                    # estabilidad incondicional/condicional; cota sobre Δt si procede

out_of_scope:
  - ""                              # variantes o capacidades diferidas (HHT-α, masa lumped, multi-support, …)

acceptance:
  # Bloque de verificación obligatorio para todo solver nuevo.
  # Sigue la disciplina V&V: comparar contra solución analítica o
  # benchmark de referencia en al menos un caso de cada régimen.
  verification:
    - name: caso_lineal_con_solucion_analitica
      setup: "<configuración mínima donde la solución se conoce en forma
             cerrada: viga en voladizo, oscilador 1 GDL, autovalores
             tabulados de una malla regular, etc.>"
      expect: "<magnitud y valor analítico esperado>"
      tol_rel: 1.0e-8
    - name: convergencia_h_o_dt
      setup: "<refinamiento de malla o paso temporal, midiendo la pendiente
             del error en log–log>"
      expect: "<orden de convergencia teórico del esquema: 2 para Newmark
              average acceleration, 2 para Quad4, etc.>"
      tol_rel: 0.1                # tolerancia sobre la pendiente
    - name: recuperacion_del_caso_lineal     # solvers no lineales
      setup: "modelo enteramente lineal (todos los materiales con tangente
             constante, sin grandes desplazamientos); ejecutar con este solver
             y con el solver lineal de referencia"
      expect: "ambos resultados coinciden a paridad de bits (o muy cerca);
              el solver no lineal converge en 1 iteración por paso"
      tol_rel: 1.0e-12
  # Tests adicionales específicos del solver (degeneración a un solver más
  # simple en límites paramétricos, robustez ante mal condicionamiento,
  # invariancia bajo cambio de unidades, etc.).
  specific:
    - name: ""
      setup: ""
      expect: ""
      tol_rel: 1.0e-10

references:
  - ""                              # Crisfield (1981), Newmark (1959), Bathe (2014), etc.
```

---

## Implementación

*Rellena la IA tras programar.*

- Archivo: —
- Clase: —
- Tests:
  - —
- Notas de traducción: —

---

## Diálogo

*Preguntas, aclaraciones y hallazgos durante la implementación. Entradas fechadas.*

- *(vacío)*
