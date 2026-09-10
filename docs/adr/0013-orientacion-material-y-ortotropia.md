# ADR 0013 — Orientación material: dónde vive el ángulo de fibra en un material anisótropo

- **Estado**: aceptado
- **Fecha**: 2026-09-09
- **Alcance**: apertura de la familia de materiales **anisótropos**, empezando por elasticidad ortótropa 2D (`Orthotropic2D`). Toca: `solidum/materials/` (componente nuevo), catálogo de materiales, manuales. **No toca** el contrato `Material` ni el contrato `Element` — esa es precisamente la decisión que este ADR registra y acota.

## Contexto

El proyecto de investigación PAPIIT del usuario sobre **bambú** exige modelos de material ortótropos. El bambú es ortótropo de forma severa: la relación $E_L/E_T$ ronda 20:1, muy por encima de la madera. Modelarlo como isótropo no es una aproximación tosca — es otro material.

El catálogo actual tiene 14 materiales y **todos son isótropos**. La palabra "ortótropo" no aparece en el código: verificado 2026-09-09 con `grep -i "orientacion|orientation|anisotr|ortotrop|local_axes|fiber_dir"` sobre `solidum/`, cero resultados. No existe ningún concepto de orientación material.

### La matriz constitutiva no es el problema

Escribir $\mathbf C$ ortótropa es álgebra de libro: cuatro constantes independientes en 2D ($E_1$, $E_2$, $G_{12}$, $\nu_{12}$) ligadas por la reciprocidad $\nu_{12}/E_1 = \nu_{21}/E_2$. No hay variables internas, ni historia, ni return mapping. Es **más simple** que `VonMises3D` o `IsotropicDamage2D`, ya implementados y validados.

### El problema es la orientación

Un material ortótropo es inútil salvo que sus ejes principales coincidan con los globales. En ejes del material no hay acoplamiento (bloque cortante nulo); al rotar un ángulo $\theta$, la constitutiva $\mathbf C_{\text{glob}} = \mathbf T^\top \mathbf C_{\text{mat}} \mathbf T$ **se llena por completo** y aparece acoplamiento tracción–cortante ($C_{16}, C_{26} \neq 0$): una tracción pura genera distorsión angular. Ese acoplamiento no existe en ningún material isótropo.

La orientación es, por tanto, **parte de la física del problema**, no un detalle de implementación. La pregunta arquitectural es dónde vive.

### Las alternativas reales

Son **dos**, no tres — y conviene decirlo así porque en la discusión inicial se presentaron tres opciones al mismo nivel, lo cual fue impreciso:

- **A — la orientación vive en el material**: `Orthotropic2D(..., theta=30.0)`. El material precalcula $\mathbf C_{\text{glob}}$ en el constructor y `compute_state` la devuelve ya rotada, exactamente como `Elastic2D` precalcula `self.C` según la hipótesis. **Cero cambios al código existente.**
- **B — la orientación vive en el elemento**: cada elemento declara sus ejes locales de material y se los pasa al material. Es lo que hacen Abaqus (`*ORIENTATION`) y ANSYS. Requiere ampliar la firma de `compute_state`.
  - **B′ — por punto de Gauss**: variante de B donde la orientación se evalúa dentro del bucle de integración en vez de una vez por elemento. Necesaria para **ortotropía cilíndrica** (el culmo real de bambú, cuyos ejes materiales siguen la geometría). **No es una tercera opción rival**: es la misma decisión arquitectural de B con distinta granularidad, y se implementa como extensión natural de B, no como alternativa.

### Coste de intrusión, medido

Verificado 2026-09-09: hay **37 llamadas a `compute_state` en 13 archivos**. Ese número engaña. Los sitios reales que B tocaría son **seis**:

- Los seis archivos de elementos 1D (`truss.py`, `cable.py`, `frame3d.py`, `frame/euler.py`, `frame/timoshenko.py`, `frame/euler_corot.py`) tienen `STRAIN_DIM = 1`. **La ortotropía no les aplica**: un material 1D tiene un solo módulo. No se tocan.
- Los sólidos ya están centralizados: `solid_2d/_shared.py` y `solid_3d/_shared.py` (bases `_HigherOrderSolid2D` / `_HigherOrderSolid3D` de la sub-etapa A.ter). Los cuadráticos heredan de ahí.
- Quedan esas dos bases más `Quad4`, `Tri3`, `Hex8`, `Tet4`, que tienen bucle propio.

## Decisión

### 1. Se adopta **A** ahora, con **B** declarado como destino

`Orthotropic2D` recibe `theta` en el constructor y precalcula la constitutiva rotada. No se modifica ni el contrato `Material` ni el contrato `Element`. La entrega es **puramente aditiva**: un archivo nuevo en `solidum/materials/`, su decorador de registro, sus tests.

**Por qué A y no B como primer paso.** No es que B sea peor — es casi con certeza donde acabará esto (ver §2). Es que B tiene el peor perfil coste/beneficio *como punto de partida*, y esto es una particularidad de este caso, no una regla general:

- **B cuesta lo mismo que B′** (mismos seis archivos, misma ampliación de firma, mismo ADR) y **resuelve menos**: sirve para fibra constante por elemento, no para ortotropía cilíndrica.
- **B da lo mismo que A** en el caso que la investigación ataca primero — probetas y tiras a ángulo constante, donde la fibra es uniforme. Ahí A y B producen resultados idénticos.

Es decir: si el caso es probetas, A entrega lo mismo sin tocar nada. Si el caso fuera culmos completos, B se queda corto y habría que ir a B′ igualmente. B queda apretado en medio.

**El argumento de fondo** es `Reglas.md §1`: resistir abstracciones especulativas no justificadas por al menos dos casos reales. Hoy hay **un** caso real (probetas de bambú a distintos ángulos) y una incertidumbre genuina sobre si el siguiente será B o B′ — incertidumbre que depende de resultados experimentales que aún no existen. Fijar ahora el contrato `Element` sería decidir sobre una suposición acerca de la investigación, no sobre un caso real.

Precedente inverso deliberado: en la Etapa 8 la base `_ThermalSolid` se creó con el **primer** elemento y no con el segundo, contra la regla habitual, porque la estructura de $\mathbf C\dot{\mathbf T} + \mathbf K\mathbf T = \mathbf F$ era manifiestamente idéntica en 2D y 3D. Aquí ocurre lo contrario: la forma final de la abstracción **no** es evidente todavía, y ese es justamente el criterio que distingue un caso del otro.

### 2. `theta` en el material es un **atajo consciente y desechable**

Se registra explícitamente que **A es conceptualmente incorrecta a largo plazo**: la orientación es una propiedad del *elemento en la malla*, no del material. Dos elementos del mismo culmo comparten $E_1, E_2, G_{12}, \nu_{12}$ —son el mismo material— pero tienen fibra en direcciones distintas. Meter $\theta$ en el material fuerza a duplicar el material para cambiar un dato geométrico.

**Condición de migración a B**: cuando aparezca un caso real con orientación variable en la malla. La vía es aditiva y no tira trabajo:

```python
compute_state(strain, state_vars=None, orientation=None)
```

con `orientation=None` significando "sin rotación". Todos los materiales isótropos ignoran el parámetro; `Orthotropic2D` usa el `theta` del constructor como *default* cuando el elemento no impone orientación. **Ningún test existente debe cambiar** — ése es el criterio de que la migración es aditiva y no ruptura de contrato.

### 3. Alcance acotado de la primera entrega

Decidido en sesión con el usuario:

- **Sólo `plane_stress`.** La hipótesis `plane_strain` ortótropa **no cierra sin datos 3D**: requiere $E_3$, $\nu_{13}$, $\nu_{23}$. En isótropo esto no se nota porque $E_3 = E$ por definición; en ortótropo $E_3 \neq E_1$. Se declara fuera de alcance en vez de pedir constantes que la campaña experimental quizá no haya medido. Para probetas y ensayos de bambú, plane stress es la hipótesis natural.
- **Sólo 2D.** `Orthotropic3D` (nueve constantes independientes, rotación 6×6) queda para cuando lo pida un caso real — regla de los dos casos.
- **Sólo elasticidad.** Los criterios de falla ortótropos (Tsai-Wu, Hashin, Hoffman) son un componente distinto y posterior, con el agravante de que el bambú falla de modo muy distinto a lo largo y a través de la fibra (el *split* longitudinal es su modo característico).
- **Sin material gradado.** La variación radial de densidad de haces vasculares en la pared del culmo es una decisión distinta y más ambiciosa (propiedad función de la posición, no constante por elemento).

### 4. Admisibilidad: **no** se puede copiar la validación de `Elastic2D`

`Elastic2D` valida $-1 < \nu < 0.5$ ([`elastic_2d.py:20`](../../solidum/materials/elastic_2d.py#L20)). Ese criterio es **incorrecto para ortótropo**. La condición correcta es que $\mathbf C$ sea definida positiva:

$$E_1, E_2, G_{12} > 0, \qquad 1 - \nu_{12}\nu_{21} > 0 \iff |\nu_{12}| < \sqrt{E_1/E_2}$$

Con $E_1/E_2 = 20$ el límite es $|\nu_{12}| < 4.47$: **$\nu_{12} > 0.5$ es legítimo** en un material muy ortótropo. Un validador heredado de `Elastic2D` rechazaría datos físicamente válidos de bambú. No se modifica `Elastic2D` — su criterio es correcto para isótropo; simplemente no se reutiliza.

### 5. Convención de subíndices de Poisson, fijada explícitamente

Existen dos convenciones opuestas en la literatura para $\nu_{ij}$, y confundirlas **invierte los papeles de $E_1$ y $E_2$** sin que nada falle ruidosamente. Se adopta la de Jones, Tsai y la mayoría de textos de composites:

$$\nu_{12} = -\frac{\varepsilon_{22}}{\varepsilon_{11}} \quad \text{bajo carga uniaxial en la dirección 1}$$

**Primer índice = dirección de carga; segundo = dirección de la contracción medida.** Debe quedar en el docstring de la clase y en la spec, no sólo aquí.

## Consecuencias

**Positivas**

- Entrega puramente aditiva: cero riesgo de regresión, los 1176 tests siguen válidos por construcción y no por verificación. Diff revisable en minutos, coherente con el modelo de colaboración de `Reglas.md §2` donde el usuario es la salvaguarda primaria leyendo diffs.
- Desbloquea de inmediato la calibración experimental del PAPIIT (probetas a distintos ángulos) — que es precisamente lo que generará la información para decidir entre B y B′.
- El isótropo queda como **caso particular verificable**: con $E_1 = E_2 = E$, $\nu_{12} = \nu$, $G_{12} = E/[2(1+\nu)]$ debe recuperarse exactamente la $\mathbf C$ de `Elastic2D` en plane stress. Es el test de regresión más valioso de la entrega.
- Abarata una eventual Etapa de **placas y láminas** (opción B del ROADMAP): las láminas ortótropas —composites laminados— son el caso de uso dominante de esa familia.

**Negativas, aceptadas**

- Un material por orientación. En una malla con fibra variable haría falta una instancia por ángulo. Es el límite conocido de A, y la condición de migración de §2.
- La decisión correcta a largo plazo (B) queda diferida, con el riesgo de que el atajo se fosilice. Mitigación: §2 lo declara desechable y fija la condición de retoma; el ROADMAP y `docs/STATUS.md` lo registran.

**Neutras**

- `IS_SYMMETRIC = True` se mantiene: la constitutiva ortótropa es simétrica incluso rotada (lo garantiza la existencia de energía de deformación), así que el despachador algebraico del ADR 0003 la sigue tratando como candidata a Cholesky.
- No hay variables internas: `PRIMARY_STATE_VAR = None`, como en el resto de materiales puramente elásticos.

## Referencias

- Jones R.M. (1999). *Mechanics of Composite Materials*, 2ª ed. Taylor & Francis. Cap. 2 (constitutiva ortótropa 2D, reciprocidad, transformación de ejes).
- Lekhnitskii S.G. (1963). *Theory of Elasticity of an Anisotropic Elastic Body*. Holden-Day. (Restricciones de admisibilidad por definición positiva.)
- Ting T.C.T. (1996). *Anisotropic Elasticity: Theory and Applications*. Oxford University Press.
