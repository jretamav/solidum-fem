# ThermalConduction — conducción de calor de Fourier con conductividad tensorial

> Spec pre-redactada por la IA (2026-08-25) sobre las seis decisiones de alcance cerradas con el usuario para la Etapa 8. **El usuario revisa la física; el plumbing está resuelto.** Los puntos abiertos de §Diálogo son físicos o de alcance, nunca arquitecturales.
>
> Primer material de la familia térmica. Etapa 8 — C1 núcleo mínimo: conducción pura, sin convección, radiación ni acoplamiento termomecánico.

---

## Especificación física

### 0. Descripción general

Modelo constitutivo de **conducción de calor** según la ley de Fourier. Relaciona el vector de flujo de calor $\mathbf q$ con el gradiente de temperatura $\nabla T$. Es el análogo térmico de `Elastic2D`/`Elastic3D`: lineal, sin variables internas, con "tangente" constante.

No pertenece a la familia `Material` mecánica. Los materiales mecánicos relacionan $\boldsymbol\sigma$ con $\boldsymbol\varepsilon$ en notación Voigt; éste relaciona dos vectores en el espacio físico, sin notación Voigt y sin tensor simétrico de por medio. Vive en una **familia paralela** `ThermalMaterial` con su propio registro, exactamente como `CohesiveMaterial` (ADR 0010) resolvió el mismo tipo de incompatibilidad semántica.

### 1. Ley constitutiva — Fourier

$$\mathbf q = -\,\mathbf k\, \nabla T$$

con:

- $\mathbf q$ — vector de flujo de calor [W/m²]. Apunta en la dirección del transporte de energía.
- $\nabla T$ — gradiente de temperatura [K/m].
- $\mathbf k$ — tensor de conductividad térmica [W/(m·K)]. Matriz $2\times2$ en 2D, $3\times3$ en 3D.

**El signo negativo es la física, no una convención elegible**: el calor fluye de mayor a menor temperatura (segunda ley de la termodinámica). De ahí se sigue que $\mathbf k$ debe ser **definida positiva** — si no lo fuera, existiría una dirección en la que el calor fluiría espontáneamente hacia la región más caliente.

### 2. Conductividad tensorial e isótropa

El contrato guarda siempre $\mathbf k$ como tensor (decisión 4 del alcance). El caso isótropo es el particular

$$\mathbf k = k\,\mathbf I$$

y el constructor acepta un escalar $k$, expandiéndolo internamente. Esto evita romper el contrato cuando entren materiales ortótropos (madera, composites, roca estratificada, hormigón fibroreforzado), habituales en el dominio del usuario.

**Requisitos sobre $\mathbf k$**:

- **Simétrica**: $\mathbf k = \mathbf k^\top$. Lo exigen las relaciones recíprocas de Onsager para el transporte en medios sin campo magnético.
- **Definida positiva**: todos los autovalores $> 0$. Equivale a exigir $k > 0$ en el caso isótropo.

### 3. Ecuación de balance que este material alimenta

El material aporta $\mathbf k$ (y, para el transitorio, $\rho c$) a la ecuación de conducción:

$$\rho c\,\frac{\partial T}{\partial t} = -\nabla\cdot\mathbf q + Q = \nabla\cdot(\mathbf k\nabla T) + Q$$

En régimen **estacionario** el término temporal se anula:

$$\nabla\cdot(\mathbf k\nabla T) + Q = 0$$

con $Q$ la fuente volumétrica de calor [W/m³].

### 4. Capacidad calorífica volumétrica

El producto $\rho c$ [J/(m³·K)] es lo que gobierna el término transitorio:

- $\rho$ — densidad [kg/m³].
- $c$ — calor específico [J/(kg·K)].

Se declaran por separado —no como producto— para que $\rho$ tenga el mismo significado y las mismas unidades que en los materiales mecánicos, y para que el usuario pueda tomar $c$ directamente de tabla.

**Ambas son opcionales al construir**, siguiendo el precedente de `density` en el ADR 0008: un análisis **estacionario** no necesita ninguna de las dos. Si un solver transitorio encuentra `None`, falla con `ValueError` identificando el material — nunca con capacidad cero silenciosa.

La **difusividad térmica** $\alpha = k/(\rho c)$ [m²/s] es magnitud derivada, no parámetro: gobierna la velocidad de propagación del frente térmico y sirve para estimar el paso de tiempo característico. Se expone como propiedad calculada, no se almacena.

### 5. Variables internas

**Ninguna.** El modelo es lineal y sin memoria: $\mathbf q$ depende sólo del $\nabla T$ instantáneo. `PRIMARY_STATE_VAR = None`, igual que en `Elastic1D`/`Elastic2D`/`Elastic3D`.

La dependencia $\mathbf k(T)$ —conductividad función de la temperatura, que sí introduciría no linealidad— queda fuera de alcance (§out_of_scope).

### 6. Unidades y consistencia

El sistema **no convierte unidades** (`Reglas.md §5`); la consistencia es responsabilidad del usuario. En SI: $[k]$ = W/(m·K), $[c]$ = J/(kg·K), $[\rho]$ = kg/m³, $[Q]$ = W/m³, $[T]$ = K.

**Sobre la escala de temperatura**: la conducción pura sólo involucra **diferencias** de temperatura, así que Kelvin y Celsius dan resultados idénticos. La distinción sí importaría con radiación (donde $T^4$ exige escala absoluta), que está fuera de esta etapa. Se documenta para que la extensión futura no herede una ambigüedad.

---

## Formulación numérica

### 7. Contribución a la matriz de conductividad

El elemento integra la forma débil de la ecuación de balance. El material aporta $\mathbf k$ al integrando:

$$\mathbf K_e = \int_\Omega \mathbf B^\top\,\mathbf k\,\mathbf B\;d\Omega, \qquad \mathbf B = \nabla\mathbf N$$

donde $\mathbf B$ es $2\times n_{nodos}$ en 2D y $3\times n_{nodos}$ en 3D — **el gradiente de las funciones de forma, sin reordenamiento de Voigt**. Es el mismo $\partial N/\partial x$ que ya calculan los kernels mecánicos; cambia sólo cómo se ensambla.

$\mathbf K_e$ es **simétrica** por serlo $\mathbf k$, luego el despachador algebraico (ADR 0003) puede elegir Cholesky. `IS_SYMMETRIC = True`.

### 8. Contribución a la matriz de capacidad

$$\mathbf C_e = \int_\Omega \rho c\;\mathbf N^\top\mathbf N\;d\Omega$$

Estructuralmente idéntica a la matriz de masa consistente, con $\rho c$ en lugar de $\rho$. Por la decisión 5 del alcance se ofrecen ambas formas, **con `lumped` por defecto** — al revés que en dinámica estructural. El motivo es físico y se detalla en la spec del elemento: la capacidad consistente produce oscilaciones espurias ante un frente térmico abrupto y puede violar el principio del máximo de la ecuación de difusión.

### 9. Interfaz de cómputo

$$\texttt{compute\_flux}(\nabla T) \;\longrightarrow\; (\mathbf q,\; \mathbf k)$$

Devuelve el flujo y el tensor de conductividad —análogo a `(σ, C_tangent)` del contrato mecánico—. Sin `state_vars`: el modelo no tiene memoria, y añadir el parámetro por simetría con `Material` sería ruido. Cuando entre un material térmico no lineal, la firma se ampliará entonces, con el caso real delante.

### 10. Caveats numéricos

- **Ninguna singularidad**: el modelo es lineal y $\mathbf k$ definida positiva por validación al construir. No hay régimen degenerado análogo al ápice de Drucker-Prager o a $d \to 1$ en daño.
- **Condicionamiento con contraste alto de conductividad**: en un dominio multimaterial con varios órdenes de magnitud entre las $k$ (p. ej. metal contra aislante), $\mathbf K$ queda mal condicionada. Es propiedad del problema físico, no del modelo; se documenta como caveat operativo.
- **Anisotropía severa**: con $k_{\max}/k_{\min}$ muy grande el problema se vuelve fuertemente direccional y la malla debe alinearse razonablemente con los ejes de anisotropía, o la solución degrada. Igual que arriba: propiedad del problema.

---

## Contrato de implementación

```yaml
name: ThermalConduction
kind: thermal_material
status: draft            # draft → implemented → validated

interface:
  field: temperature     # campo primario que gobierna (1 DOF escalar por nodo)
  flux_dim: 0            # 2 en 2D · 3 en 3D — se fija al construir según la dimensión
  primary_state_var: null
  is_symmetric: true

parameters:
  - { name: k,       type: float | array-like, required: true,
      desc: "Conductividad térmica [W/(m·K)]. Escalar ⇒ isótropo, se expande a k·I.
             Matriz 2×2 o 3×3 ⇒ ortótropo/anisótropo; debe ser simétrica y definida positiva" }
  - { name: c,       type: float, required: false, default: null,
      desc: "Calor específico [J/(kg·K)]. Opcional al construir; obligatorio solo
             en análisis transitorio (mismo criterio que density en ADR 0008)" }
  - { name: density, type: float, required: false, default: null,
      desc: "Densidad [kg/m³]. Opcional al construir; obligatoria solo en transitorio" }

signature:
  compute_flux: "(grad_T: ndarray(d,)) -> (q: ndarray(d,), k: ndarray(d,d))  con d ∈ {2, 3}"
  gradient_kind: "gradiente físico ∇T = [∂T/∂x, ∂T/∂y(, ∂T/∂z)] — sin notación Voigt"

state_variables: []

conventions:
  sign: "q = -k·∇T; el flujo va de mayor a menor temperatura (segunda ley)"
  units: "k [W/(m·K)] · c [J/(kg·K)] · density [kg/m³] · Q [W/m³] · T [K o °C]"
  temperature_scale: "sólo intervienen diferencias de T ⇒ K y °C equivalentes en
                      conducción pura. La distinción importará al entrar radiación (T⁴)"
  derived: "difusividad térmica α = k/(ρ·c) [m²/s] expuesta como propiedad calculada,
            no almacenada"

validity:
  - "k simétrica y definida positiva (autovalores > 0); escalar k > 0"
  - "c > 0 y density > 0 cuando se declaran"
  - "conductividad independiente de la temperatura (modelo lineal)"
  - "medio continuo homogéneo dentro de cada elemento"

out_of_scope:
  - "conductividad dependiente de la temperatura k(T) ⇒ introduce no linealidad; etapa futura"
  - "convección (Robin) y radiación ⇒ condiciones de frontera, no del material; fuera de C1"
  - "cambio de fase / calor latente (Stefan) ⇒ etapa futura"
  - "acoplamiento termomecánico y deformación térmica α·ΔT ⇒ decisión 2 del alcance"
  - "medios porosos, convección forzada interna, transporte de masa"

acceptance:
  verification:
    - name: ley_de_fourier_exacta
      setup: "compute_flux con ∇T arbitrario; k isótropo y k ortótropo"
      expect: "q = -k·∇T exacto componente a componente; segundo retorno = k"
      tol_rel: 1.0e-14
    - name: expansion_escalar_a_tensor
      setup: "construir con k escalar en 2D y en 3D"
      expect: "k interno == k·I (2×2 y 3×3 respectivamente)"
      tol_abs: 1.0e-14
    - name: simetria_y_positividad
      setup: "construir con k escalar y con k tensorial válido"
      expect: "k == k.T exacto; todos los autovalores estrictamente positivos"
      tol_abs: 1.0e-14
    - name: flujo_opuesto_al_gradiente
      setup: "k isótropo; ∇T unitario en varias direcciones"
      expect: "q antiparalelo a ∇T; q·∇T < 0 estricto (disipación termodinámica)"
      tol_rel: 1.0e-14
    - name: anisotropia_desvia_el_flujo
      setup: "k = diag(k1, k2) con k1 ≠ k2; ∇T a 45° de los ejes"
      expect: "q NO es antiparalelo a ∇T; se desvía hacia el eje de mayor conductividad.
               Verificación cuantitativa contra q = -k·∇T calculado a mano"
      tol_rel: 1.0e-12
    - name: difusividad_derivada
      setup: "k, c, density conocidos"
      expect: "alpha == k/(density·c) para el caso isótropo"
      tol_rel: 1.0e-14
    - name: rechazo_inputs_invalidos
      setup: "k ≤ 0; k tensorial no simétrico; k tensorial no definido positivo;
              c ≤ 0; density ≤ 0; k de dimensión inconsistente con el elemento"
      expect: "ValueError con mensaje claro identificando el parámetro en cada caso"
    - name: capacidad_ausente_falla_en_transitorio
      setup: "construir sin c ni density; solicitar matriz de capacidad"
      expect: "ValueError identificando el material y el parámetro faltante,
               nunca capacidad cero silenciosa (precedente ADR 0008)"
    - name: mensaje_de_error_accionable
      setup: "construir sin c y/o sin density; solicitar matriz de capacidad"
      expect: "el mensaje nombra el material, nombra el/los parámetro(s) que faltan
               e indica la salida ('declara ambas' | 'usa un solver estacionario').
               No basta con reportar el atributo ausente"
    - name: estacionario_no_exige_capacidad
      setup: "construir sólo con k; ensamblar la matriz de conductividad"
      expect: "K_e se ensambla sin error; c y density no se consultan en ningún punto
               del camino estacionario"
      tol_rel: 1.0e-14

references:
  - "Incropera F.P., DeWitt D.P. (2007). Fundamentals of Heat and Mass Transfer. Wiley. §2 (ley de Fourier, conductividad tensorial)."
  - "Carslaw H.S., Jaeger J.C. (1959). Conduction of Heat in Solids. Oxford. §1 (formulación general, requisitos sobre k)."
  - "Lewis R.W., Nithiarasu P., Seetharamu K.N. (2004). Fundamentals of the Finite Element Method for Heat and Fluid Flow. Wiley. §3."
  - "Zienkiewicz O.C., Taylor R.L. (2000). The Finite Element Method, Vol. 1. Butterworth-Heinemann. §7 (problemas de campo escalar)."
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

**2026-08-25 — IA → usuario. Tres puntos:**

1. ~~**Nombre del material.**~~ ✅ **Resuelto 2026-08-25 — `ThermalConduction`** (decisión del usuario). El nombre inicial `ThermalIsotropic` describía mal el contrato, que admite tensor completo; `ThermalConduction` nombra la física y no compromete la simetría del tensor. Un futuro material térmico no lineal —conductividad $\mathbf k(T)$, cambio de fase— entrará como nombre propio distinto, no como variante de éste.

2. ~~**`c` y `density` opcionales.**~~ ✅ **Resuelto 2026-08-25 — se traslada el criterio del ADR 0008** (decisión del usuario). Ambas son opcionales al construir: un análisis **estacionario** sólo usa $\mathbf k$, y exigir propiedades que no entran en ninguna ecuación obligaría a inventar valores ficticios. El fallo silencioso —el riesgo real de un parámetro opcional— está cerrado por diseño: cuando un consumidor transitorio las encuentra `None`, lanza `ValueError` identificando el material.

   **Requisito añadido en la decisión**: el mensaje debe ser **accionable**, no sólo descriptivo — decir qué hacer, no sólo qué falta. Forma canónica:

   > `ThermalConduction (id=1): el análisis transitorio requiere 'c' y 'density'. Declara ambas en el material, o usa un solver estacionario (LinearSolver).`

   Es el espíritu de "validación temprana con mensajes claros" de `Reglas.md §1` aplicado al caso.

3. ~~**Fuente volumétrica $Q$.**~~ ✅ **Resuelto 2026-08-25 — $Q$ es carga aplicada en el elemento, no propiedad del material** (decisión del usuario).

   El material define la **ley constitutiva**, no las acciones sobre el medio; es el mismo reparto que en la familia mecánica, donde `Elastic2D` no sabe nada del peso propio ni de las fuerzas de cuerpo.

   La razón que decide no es la simetría arquitectural sino la **variabilidad espacial y temporal**: la generación por hidratación del hormigón depende de la madurez, y el efecto Joule de la corriente. Ninguna es constante del medio. Modelada como carga, $Q$ puede variar por zona y por paso sin duplicar materiales; como propiedad del material, cada valor distinto exigiría un material propio.

   Consecuencia para la spec del elemento: `compute_body_source(Q)` integra $\int_\Omega Q\,\mathbf N^\top d\Omega$, paralelo a `compute_body_load(b)` del mecánico. Se especifica allí, no aquí.

   No cierra la puerta a un atajo futuro —un material que declare generación propia— si aparece el caso de uso que lo justifique.
