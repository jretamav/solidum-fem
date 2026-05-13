# <Nombre> — <descripción corta>

> Orden de trabajo. El usuario escribe la **especificación física**, la **formulación numérica** y el **contrato**. La IA rellena **implementación** y responde en **diálogo** durante el trabajo.

---

## Especificación física

### 0. Descripción general
*Tipo de modelo constitutivo (elástico, plástico, daño, viscoelástico…); dimensionalidad (`STRAIN_DIM` ∈ {1, 3, 6}); hipótesis cinemáticas (small strain, finite strain); hipótesis 2D (`plane_strain` / `plane_stress`) si aplica.*

### 1. Descomposición de la deformación
*Si el modelo lo requiere: $\boldsymbol\varepsilon = \boldsymbol\varepsilon^e + \boldsymbol\varepsilon^p$ (plasticidad), $\boldsymbol\varepsilon = \boldsymbol\varepsilon^e$ (elasticidad pura), etc.*

### 2. Ley elástica (relación esfuerzo–deformación elástica)
*$\boldsymbol\sigma = \mathbf C_e : \boldsymbol\varepsilon^e$ con $\mathbf C_e$ parametrizado por las constantes elásticas del modelo (isótropo: $E,\nu$; ortótropo: $E_1, E_2, G_{12}, \nu_{12}$; …).*

### 3. Superficie/criterio de fluencia, daño o evolución
*Si aplica: función $f(\boldsymbol\sigma, q)$ que delimita régimen elástico vs inelástico; criterio de iniciación de daño; energía libre $\Psi$ para viscoelasticidad. Variables internas $q$ del modelo.*

### 4. Regla de flujo / ley de evolución
*Cómo evolucionan las variables internas. Plasticidad asociada: $\dot{\boldsymbol\varepsilon}^p = \dot\gamma\,\partial f/\partial\boldsymbol\sigma$. No asociada: potencial plástico $g \neq f$. Daño: $\dot d = \dot d(\kappa)$. Viscoplasticidad: $\dot\gamma = \langle\phi(f)\rangle/\eta$.*

### 5. Endurecimiento / ablandamiento / saturación
*Ley que controla la expansión, traslación o ablandamiento de la superficie de fluencia o la evolución del daño. Endurecimiento isótropo lineal, exponencial, Voce, etc.*

### 6. Condiciones de Kuhn-Tucker (si aplica)
*Plasticidad rate-independent y daño rate-independent imponen $\dot\gamma \ge 0$, $f \le 0$, $\dot\gamma\,f = 0$ (complementariedad).*

### 7. Variables internas
*Lista con tipo, dimensión y significado físico. Cuál es `PRIMARY_STATE_VAR` (la variable exportada al post-procesado).*

### 8. Notación Voigt y manejo de componentes auxiliares
*Convención del proyecto: $\boldsymbol\varepsilon = [\varepsilon_{xx}, \varepsilon_{yy}, \gamma_{xy}]$ con $\gamma_{xy} = 2\varepsilon_{xy}$. Si el modelo necesita componentes 3D internas (e.g. $\varepsilon^p_{zz}$ en plasticidad 2D), documentarlo aquí.*

---

## Formulación numérica

### 9. Esquema temporal — discretización implícita
*Backward Euler implícito sobre las variables internas (plasticidad, daño, viscoelasticidad). Dado el estado $\{q_n\}$ al paso $n$ y la deformación total $\boldsymbol\varepsilon_{n+1}$, resolver para $\{\boldsymbol\sigma_{n+1}, q_{n+1}\}$.*

### 10. Return mapping / actualización del estado
*Algoritmo central: estado trial elástico, comprobación de admisibilidad, corrector si activo. Cuando aplica, expresión cerrada para $\Delta\gamma$, $\Delta d$, etc. Si requiere Newton local (no asociada, plane stress proyectado…), describir el sistema de ecuaciones a resolver.*

### 11. Tangente algorítmica consistente
*$\mathbf C^{\text{alg}} = \partial\boldsymbol\sigma_{n+1}/\partial\boldsymbol\varepsilon_{n+1}$ derivada del algoritmo discretizado, no de la ley continua. Imprescindible para convergencia cuadrática del Newton global. Indicar simetría: si es no simétrica, declarar `IS_SYMMETRIC = False` para que el despachador algebraico ADR 0003 elija un backend adecuado (LU).*

### 12. Caveats numéricos
*Singularidades del modelo (e.g. ápice de Drucker-Prager, $\varepsilon = 0$ en cable, $d \to 1$ en daño saturado), heurísticas de regularización, comportamiento esperado de Newton-Raphson cerca de transiciones.*

---

## Contrato de implementación

```yaml
name: <Nombre>
kind: material
status: draft            # draft → implemented → validated

interface:
  strain_dim: 0          # 1 axial · 3 plano (ε_xx, ε_yy, γ_xy) · 6 3D
  primary_state_var: ""  # nombre de la variable interna exportada al post
  is_symmetric: true     # tangente algorítmica simétrica? (afecta al despachador algebraico ADR 0003)

parameters:
  - { name: , type: , required: true, desc:  }

state_variables:
  # Variables internas que el material guarda y actualiza entre commits.
  - { name: , type: , shape: , desc:  }

conventions:
  voigt: ""              # orden de componentes, factor en cortante
  units: ""              # unidades de cada parámetro
  sign: ""               # convención de signos (tracción positiva, etc.)

validity:
  - ""                   # regímenes de aplicación válidos (e.g. |ε| ≲ 1e-2)

out_of_scope:
  - ""                   # hipótesis no soportadas, características diferidas

acceptance:
  # Bloque de verificación obligatorio para todo material nuevo.
  # Sigue la disciplina V&V: el modelo numérico reproduce la respuesta
  # física esperada en casos canónicos.
  verification:
    - name: respuesta_elastica
      setup: "carga monotónica en régimen elástico (por debajo del umbral
             de fluencia / daño); comparar σ vs ε con la ley elástica pura"
      expect: "σ = C_e:ε exacto; tangente = C_e elástica"
      tol_rel: 1.0e-12
    - name: respuesta_inelastica_uniaxial
      setup: "<test uniaxial monotónico con solución analítica conocida:
             curva σ–ε de plasticidad J2 con endurecimiento lineal,
             evolución de daño con softening exponencial, etc.>"
      expect: "<curva analítica tabulada o forma cerrada>"
      tol_rel: 1.0e-8
    - name: ciclo_carga_descarga
      setup: "ciclo carga → descarga → recarga; verificar que las variables
             internas evolucionan solo en carga activa y permanecen
             constantes en descarga"
      expect: "histéresis correcta; estado se conserva en descarga"
      tol_rel: 1.0e-10
  # Tests adicionales específicos del modelo (degeneración a otro material
  # conocido en límites paramétricos, simetría/asimetría de la tangente,
  # invariancia bajo cambio de unidades, etc.).
  specific:
    - name: ""
      setup: ""
      expect: ""
      tol_rel: 1.0e-10

references:
  - ""                   # Simó & Hughes (1998), de Souza Neto et al. (2008), etc.
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
