# ThetaMethodSolver — integración temporal de primer orden (θ-method)

> Spec pre-redactada por la IA (2026-08-25). **El usuario revisa la física; el plumbing está resuelto.**
>
> Cuarto y último componente de la Etapa 8 — C1 núcleo mínimo. Resuelve el transitorio de conducción sobre [`Quad4Thermal`](Quad4Thermal.md) y [`Hex8Thermal`](Hex8Thermal.md).

---

## Especificación física

### 0. Descripción general y por qué es un solver nuevo

Integrador temporal para sistemas semidiscretos de **primer orden**:

$$\mathbf C\,\dot{\mathbf T} + \mathbf K\,\mathbf T = \mathbf F(t)$$

**No es una variante de `NewmarkSolver`.** La familia Newmark integra la ecuación de **segundo orden** $\mathbf M\ddot{\mathbf u} + \mathbf C\dot{\mathbf u} + \mathbf K\mathbf u = \mathbf F$ mediante hipótesis sobre la aceleración dentro del paso. En conducción no existe segunda derivada temporal: no hay aceleración, no hay inercia, y las hipótesis de Newmark no tienen sobre qué aplicarse. Forzar Newmark con $\mathbf M = \mathbf 0$ degenera el esquema.

Es, por tanto, un solver de familia propia, no una subclase — y no le corresponde spec corta de extensión sino spec completa (`Reglas.md §4`).

El **régimen estacionario** no lo usa: $\mathbf K\mathbf T = \mathbf F$ lo resuelve el `LinearSolver` existente sin modificación.

### 1. Esquema θ

Se aproxima la ecuación en un punto intermedio del paso, ponderando el estado con el parámetro $\theta \in [0,1]$:

$$\mathbf C\,\frac{\mathbf T_{n+1} - \mathbf T_n}{\Delta t} + \mathbf K\left[\theta\,\mathbf T_{n+1} + (1-\theta)\,\mathbf T_n\right] = \theta\,\mathbf F_{n+1} + (1-\theta)\,\mathbf F_n$$

Reordenando, el sistema a resolver por paso:

$$\underbrace{\left(\mathbf C + \theta\,\Delta t\,\mathbf K\right)}_{\mathbf A}\,\mathbf T_{n+1} = \left[\mathbf C - (1-\theta)\,\Delta t\,\mathbf K\right]\mathbf T_n + \Delta t\left[\theta\,\mathbf F_{n+1} + (1-\theta)\,\mathbf F_n\right]$$

$\mathbf A$ es **simétrica y definida positiva** para $\theta > 0$ (suma de $\mathbf C$ definida positiva y $\mathbf K$ semidefinida, con Dirichlet impuesto), luego el despachador algebraico (ADR 0003) elige Cholesky. Con $\Delta t$ constante, $\mathbf A$ **se factoriza una sola vez** y se reutiliza en todos los pasos.

### 2. Casos particulares de θ

| $\theta$ | Nombre | Orden | Estabilidad |
|---|---|---|---|
| $0$ | Euler explícito (forward) | 1 | **Condicional**: $\Delta t \le 2/\lambda_{\max}$ |
| $1/2$ | Crank-Nicolson | **2** | Incondicional (A-estable) |
| $2/3$ | Galerkin | 1 | Incondicional |
| $1$ | Euler implícito (backward) | 1 | Incondicional, **L-estable** |

Con $\theta \ge 1/2$ el esquema es **incondicionalmente estable**: ningún $\Delta t$ produce divergencia. Es la propiedad que distingue este solver del `CentralDifferenceSolver` mecánico, cuya estabilidad es siempre condicional (CFL).

### 3. Estabilidad no es lo mismo que ausencia de oscilaciones

Punto físico central de esta spec, y la razón de la elección de default.

Crank-Nicolson ($\theta = 1/2$) es de **segundo orden** y por tanto el más preciso de la tabla, pero **A-estable y no L-estable**: su factor de amplificación tiende a $-1$ para los modos de frecuencia alta, en lugar de a $0$. Consecuencia práctica: ante un **cambio brusco** —un escalón de temperatura impuesto en la frontera, típico en el arranque de un análisis— produce **oscilaciones amortiguadas lentamente** en los primeros pasos, con temperaturas que pueden salirse del rango de los datos.

Euler implícito ($\theta = 1$) es sólo de primer orden, pero **L-estable**: amortigua los modos altos por completo en un paso. Nunca oscila.

Esto es el mismo fenómeno físico que motivó el default `lumped` en la matriz de capacidad, atacado desde el otro lado: **el principio del máximo** de la ecuación de difusión dice que la solución exacta nunca excede los extremos de los datos iniciales y de frontera. Un esquema que lo viola no es simplemente impreciso — produce un resultado cualitativamente imposible.

**Por eso el default propuesto es $\theta = 1$** (Euler implícito), robusto ante el arranque típico. Crank-Nicolson queda disponible para quien busque precisión de segundo orden y sepa que su problema es suave. Se documenta el criterio para que la elección sea informada y no una preferencia oculta.

### 4. Condición inicial

El campo $\mathbf T_0$ debe declararse. Dos formas admitidas:

- **Uniforme**: un escalar, el caso más frecuente (cuerpo inicialmente a temperatura ambiente).
- **Por nodo**: vector completo, para arrancar de un estado no uniforme (p. ej. el resultado de un análisis estacionario previo).

A diferencia del mecánico, **no hay velocidad inicial**: la ecuación es de primer orden y $\mathbf T_0$ la determina por completo.

**Consistencia con Dirichlet**: si $\mathbf T_0$ contradice una temperatura impuesta en un nodo con Dirichlet, existe una discontinuidad en $t=0$. Es físicamente legítimo (un choque térmico es exactamente eso) pero es la situación que más excita las oscilaciones descritas en §3. El solver **avisa** cuando lo detecta, sin abortar: es modelización válida, no error.

### 5. Selección del paso de tiempo

El paso no se elige por estabilidad —con $\theta \ge 1/2$ no hay límite— sino por **precisión**. La escala física la da la difusividad $\alpha = k/(\rho c)$ [m²/s] del material y el tamaño característico de elemento $h$:

$$\Delta t_{\text{car}} \sim \frac{h^2}{\alpha}$$

Es el tiempo que tarda el frente térmico en atravesar un elemento. Un $\Delta t \gg \Delta t_{\text{car}}$ es estable pero no resuelve el transitorio: salta por encima de la física que se quiere ver.

El solver **no calcula ni impone** este valor —igual que el `CentralDifferenceSolver` no calcula el $\Delta t$ crítico—, pero sí lo **reporta como diagnóstico** al arrancar, para que el usuario contraste su elección. Es información, no restricción.

### 6. Sin amortiguamiento de Rayleigh

El amortiguamiento de Rayleigh $\mathbf C = \alpha\mathbf M + \beta\mathbf K$ del subsistema dinámico **no aplica**: es un modelo de disipación de la ecuación de segundo orden. En conducción la "disipación" es la propia conducción, ya contenida en $\mathbf K$. El parámetro no existe en este solver.

---

## Formulación numérica

### 7. Algoritmo

```
A = C + θ·Δt·K                     # una vez, si Δt constante
factor = factorize(A)              # Cholesky (A es SPD); reutilizada
B = C − (1−θ)·Δt·K                 # una vez

T = T_0
para cada paso n:
    rhs = B·T_n + Δt·[θ·F_{n+1} + (1−θ)·F_n]
    aplicar Dirichlet (eliminación, ADR 0004)
    T_{n+1} = factor.solve(rhs_reducido)
    almacenar según output_every
```

**Coste por paso**: una sustitución hacia adelante/atrás. La factorización es única mientras $\Delta t$ y las matrices no cambien — el problema es lineal y sin historia, así que $\mathbf K$ y $\mathbf C$ se ensamblan una sola vez.

### 8. Dirichlet dependiente del tiempo

Una temperatura impuesta puede variar en el tiempo, $\bar T(t)$ — es un caso de uso central (una superficie sometida a un ciclo térmico). La eliminación del ADR 0004 se aplica en cada paso con el valor correspondiente; el término de acoplamiento $\mathbf K_{fp}\bar{\mathbf T}_p$ se recalcula, la factorización de $\mathbf A_{ff}$ **no**.

### 9. Resultado

Devuelve un `ThermalTransientResult` con `t_history`, `T_history`, `n_steps` y `converged`.

**Decisión de zona gris que se señala**: no se reutiliza `TransientResult`, cuyos campos `udot_history`, `uddot_history`, `alpha_rayleigh` y `beta_rayleigh` no tienen significado térmico. Rellenarlos con `None` o ceros produciría un objeto que miente sobre su contenido. Se sigue el criterio ya aplicado en el proyecto con `ModalResult`, `HarmonicResult` y `ResponseSpectrumResult`: un tipo de análisis con semántica propia tiene su dataclass propia.

### 10. Caveats numéricos

- **Sin Dirichlet, $\mathbf A$ es singular** en el estacionario y mal condicionada en el transitorio con $\Delta t$ grande. El diagnóstico debe nombrar el modo de temperatura uniforme, no reportar un fallo algebraico genérico.
- **Oscilaciones con Crank-Nicolson ante escalón**: descrito en §3. No es un error del solver; se documenta y se cuantifica con test.
- **$\theta = 0$ (explícito)**: se admite por completitud pero requiere $\Delta t \le 2/\lambda_{\max}$, y con capacidad consistente $\lambda_{\max}$ es grande. Sin masa lumped es prácticamente inutilizable; se avisa.
- **Contraste alto de conductividad**: mal condicionamiento heredado del problema, no del esquema.

---

## Contrato de implementación

```yaml
name: ThetaMethodSolver
kind: solver
status: draft            # draft → implemented → validated

interface:
  pipeline_kind: thermal_transient
  field: temperature
  equation_order: 1      # primer orden en el tiempo (Newmark integra segundo)

parameters:
  - { name: theta, type: float, required: false, default: 1.0,
      desc: "Peso del esquema ∈ [0,1]. 1.0 Euler implícito (L-estable, default);
             0.5 Crank-Nicolson (2º orden, A-estable, puede oscilar ante escalón);
             2/3 Galerkin; 0.0 explícito (condicionalmente estable)" }
  - { name: dt, type: float, required: true,
      desc: "Paso de tiempo [s]. No limitado por estabilidad si θ ≥ 0.5, sí por
             precisión: contrastar con el Δt característico h²/α reportado" }
  - { name: n_steps, type: int, required: true, desc: "Número de pasos" }
  - { name: T_initial, type: float | array-like, required: true,
      desc: "Condición inicial. Escalar ⇒ campo uniforme; vector ⇒ valor por nodo" }
  - { name: lumping, type: str, required: false, default: "lumped",
      desc: "Forma de la matriz de capacidad. Default lumped por el mismo criterio
             físico que en los elementos térmicos (Quad4Thermal §5)" }
  - { name: output_every, type: int, required: false, default: 1,
      desc: "Almacenar el campo cada N pasos; limita la memoria en corridas largas" }
  - { name: linear_algebra, type: str, required: false, default: "auto",
      desc: "Backend algebraico (ADR 0003). A es SPD ⇒ Cholesky por defecto" }
  # Sin `rayleigh`: el amortiguamiento de Rayleigh es de la ecuación de 2º orden.

conventions:
  stability: "θ ≥ 0.5 incondicionalmente estable; θ = 1 además L-estable
              (amortigua modos altos por completo, nunca oscila)"
  dt_guidance: "Δt característico ~ h²/α con α = k/(ρc). Reportado como
                diagnóstico al arrancar; no impuesto"
  order_reporting: "el solver reporta el orden efectivo del esquema según θ
                    (2 si θ=0.5, 1 en el resto), para que el coste en precisión
                    del default robusto sea visible y no silencioso"
  factorization: "A = C + θ·Δt·K se factoriza UNA vez si Δt es constante"

validity:
  - "sistema lineal de primer orden: C·Ṫ + K·T = F"
  - "θ ∈ [0, 1]; dt > 0; n_steps ≥ 1"
  - "el modelo debe imponer al menos un Dirichlet de temperatura"
  - "material con c y density declaradas (ThermalConduction §4)"

out_of_scope:
  - "no linealidad: k(T), cambio de fase, radiación ⇒ requerirían Newton por paso"
  - "paso de tiempo adaptativo ⇒ diferido hasta caso de uso real"
  - "amortiguamiento de Rayleigh (no aplica a 1er orden)"
  - "acoplamiento termomecánico"

acceptance:
  verification:
    - name: decaimiento_exponencial_1dof
      setup: "sistema de 1 DOF c·Ṫ + k·T = 0 con T(0) = T_0; comparar contra
             la solución analítica T(t) = T_0·exp(−k·t/c)"
      expect: "error < tol para θ=1 y θ=0.5; el error de Crank-Nicolson decrece
              como O(Δt²) y el de Euler implícito como O(Δt) — verificar AMBAS
              tasas por refinamiento sucesivo de Δt"
      tol_rel: 1.0e-3
    - name: orden_de_convergencia_temporal
      setup: "mismo problema, Δt sucesivamente halvado"
      expect: "pendiente log-log ≈ 1 para θ=1 y ≈ 2 para θ=0.5.
               Es el test que distingue realmente los dos esquemas"
      tol_abs: 0.15
    - name: estabilidad_incondicional
      setup: "θ = 1 con Δt tres órdenes de magnitud mayor que h²/α"
      expect: "la solución NO diverge; converge monótonamente al estacionario.
               Con θ = 0 y el mismo Δt sí diverge (contraste explícito)"
    - name: convergencia_al_estacionario
      setup: "problema de pared plana con Dirichlet fijo en ambos extremos;
             integrar hasta t → ∞"
      expect: "el campo converge al perfil lineal que da el LinearSolver sobre
              el MISMO modelo; comparación nodo a nodo"
      tol_rel: 1.0e-8
  specific:
    - name: L_estabilidad_ante_escalon
      setup: "cuerpo a T_0 uniforme; escalón de Dirichlet en una cara;
             Δt grande; comparar θ=1 contra θ=0.5"
      expect: "θ=1 monótono, sin temperaturas fuera del rango [T_0, T_pared];
               θ=0.5 exhibe oscilación amortiguada en los primeros pasos.
               El test DOCUMENTA la diferencia, no la corrige"
    - name: principio_del_maximo
      setup: "θ=1, capacidad lumped, sin fuente; condiciones de frontera
             acotadas entre T_min y T_max"
      expect: "en TODO paso y TODO nodo, T ∈ [T_min, T_max]. Ningún valor
               fuera del rango de los datos"
      tol_abs: 1.0e-12
    - name: conduccion_transitoria_semi_infinita
      setup: "barra larga inicialmente a T_0, con la cara x=0 llevada
             súbitamente a T_s; comparar con la solución analítica de
             Carslaw-Jaeger T(x,t) = T_s + (T_0−T_s)·erf(x/(2√(αt)))"
      expect: "perfil coincidente en varios instantes mientras el frente no
              alcanza el extremo opuesto (hipótesis de medio semi-infinito)"
      tol_rel: 2.0e-2
    - name: balance_energetico_transitorio
      setup: "dominio aislado (todo Neumann homogéneo) con fuente Q constante"
      expect: "la energía acumulada ΣρcVΔT iguala a Q·V·t en cada paso —
              el esquema conserva la energía introducida"
      tol_rel: 1.0e-10
    - name: dirichlet_variable_en_el_tiempo
      setup: "temperatura impuesta con variación sinusoidal en el tiempo"
      expect: "el campo sigue la excitación; la factorización de A_ff se reutiliza
               (verificado por contador de llamadas) mientras el acoplamiento
               K_fp·T_p se recalcula por paso"
      tol_rel: 1.0e-10
    - name: condicion_inicial_escalar_y_vectorial
      setup: "arrancar con T_initial escalar y con vector equivalente"
      expect: "resultados idénticos; el escalar se expande correctamente"
      tol_abs: 1.0e-14
    - name: aviso_incompatibilidad_inicial_dirichlet
      setup: "T_initial contradice un Dirichlet impuesto en un nodo"
      expect: "el solver AVISA identificando los nodos, y continúa: un choque
               térmico es modelización legítima, no un error"
    - name: reporte_del_orden_efectivo
      setup: "construir el solver con θ = 1, 0.5 y 2/3"
      expect: "reporta orden 1, 2 y 1 respectivamente. Hace visible el coste en
               precisión del default robusto en vez de dejarlo implícito"
    - name: rechazo_inputs_invalidos
      setup: "θ fuera de [0,1]; dt ≤ 0; n_steps < 1; material sin c ni density"
      expect: "ValueError con mensaje accionable en cada caso"

references:
  - "Zienkiewicz O.C., Taylor R.L. (2000). The Finite Element Method, Vol. 1. Butterworth-Heinemann. §18 (esquemas de un paso para 1er orden; análisis de estabilidad)."
  - "Lewis R.W., Nithiarasu P., Seetharamu K.N. (2004). Fundamentals of the FEM for Heat and Fluid Flow. Wiley. §6 (θ-method, oscilaciones, criterios de Δt)."
  - "Hughes T.J.R. (2000). The Finite Element Method: Linear Static and Dynamic FEA. Dover. §8 (algoritmos para problemas parabólicos; A-estabilidad vs L-estabilidad)."
  - "Carslaw H.S., Jaeger J.C. (1959). Conduction of Heat in Solids. Oxford. §2.4 (sólido semi-infinito, solución con función error)."
  - "Wood W.L. (1990). Practical Time-Stepping Schemes. Oxford."
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

*Puntos abiertos para el usuario. Sólo física y alcance.*

**2026-08-25 — IA → usuario. Un punto:**

1. ~~**Valor por defecto de $\theta$.**~~ ✅ **Resuelto 2026-08-25 — $\theta = 1$ (Euler implícito)**. El usuario delegó la decisión en la IA por ser materia de criterio técnico; queda razonada aquí para que sea auditable.

   **La disyuntiva.** Crank-Nicolson ($\theta = 1/2$) es de orden 2 y más preciso para un mismo $\Delta t$; Euler implícito ($\theta = 1$) es de orden 1 pero **L-estable**. Ante el arranque típico —un escalón de temperatura en la frontera— Crank-Nicolson produce oscilaciones que pueden dar **temperaturas fuera del rango de los datos**, violando el principio del máximo.

   **El criterio que decide** no es el orden de convergencia sino qué debe hacer el programa cuando el usuario no elige: un **default protege a quien no eligió**. Los dos modos de fallo no son simétricos — un resultado impreciso se detecta refinando el paso y se corrige; un resultado *imposible* desconcierta, y quien no conozca la teoría de A-estabilidad frente a L-estabilidad no tiene forma de diagnosticarlo. La precisión perdida es recuperable por el usuario; la confianza en un resultado cualitativamente absurdo, no.

   **Coherencia interna**: es el mismo criterio que fijó `lumped` como default de la matriz de capacidad, aplicado ahora al eje temporal en vez del espacial. Ambos atacan el mismo fenómeno —oscilación espuria ante frente abrupto— y que apunten en la misma dirección hace el sistema predecible.

   **El default no encierra a nadie**: Crank-Nicolson queda a un parámetro de distancia para quien sepa que su problema es suave y busque orden 2.

   **Requisito añadido con la decisión**: el solver **reporta el orden efectivo del esquema** según el $\theta$ en uso, para que el coste en precisión del default sea visible y no una penalización silenciosa. Es la contrapartida honesta de elegir robustez.
