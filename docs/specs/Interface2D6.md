# INTERFACE 2D6 — elemento de interfaz de espesor nulo 2D con caras cuadráticas (extensión de `Interface2D4`)

> Spec corta tipo extensión (Reglas §4): la formulación de la familia está en [`Interface2D4`](Interface2D4.md) y la arquitectura en el ADR 0043. Aquí sólo lo que cambia. Incorporada al bloque de contacto por decisión del usuario del 2026-11-26 («si la 1»): sin ella, `pair` y `split` no se pueden usar con sólidos cuadráticos (Quad8, Quad9, Tri6), que el usuario usa en sus validaciones de fractura.

## Lo que cambia respecto de `Interface2D4`

### 1. Caras y orden de los nodos

Caras de **3 nodos** (dos extremos y uno intermedio), isoparamétricas: pueden ser curvas, como los lados de los sólidos cuadráticos que unen.

Orden: **primero las esquinas, después los intermedios**, como los sólidos cuadráticos de Solidum (Quad8: esquinas y después nodos de lado):

- nodos 1, 2: extremos de la cara A; nodos 3, 4: extremos de la cara B, con **4 frente a 1 y 3 frente a 2** (los cuatro primeros, igual que `Interface2D4`);
- nodo 5: intermedio de la cara A (entre 1 y 2); nodo 6: intermedio de la cara B, frente a 5.

Es el orden de la celda VTK `VTK_QUADRATIC_LINEAR_QUAD` (**tipo 30**, `quad6` en meshio: cuatro esquinas antihorarias, el intermedio del lado 1–2 y el del lado 3–4), así que la celda se escribe sin permutar (paso de salida, ADR 0043 §7). Corregido por el teórico: la 26 es la cuña cuadrática de 15 nodos.

### 2. Funciones de forma, ejes y medida

En la cara A, con `ξ ∈ [−1, 1]` (nodo 1 en `ξ = −1`, nodo 2 en `ξ = +1`, nodo 5 en `ξ = 0`):

$$N_1 = \tfrac12\,\xi(\xi-1), \qquad N_2 = \tfrac12\,\xi(\xi+1), \qquad N_5 = 1-\xi^2$$

y las mismas en la cara B, asignadas a los nodos frente a cada nodo de A (4, 3, 6).

- `a(ξ) = ∂x_A/∂ξ = Σ_i (∂N_i/∂ξ)·x_i`; `J(ξ) = ‖a(ξ)‖`; `s = a/J`, `n = e_z × s`. En una cara curva los ejes **varían de punto a punto** y se evalúan en cada punto de integración (como en el cuadrilátero alabeado 3D).
- Holgura inicial `g_0(ξ) = n(ξ)·(x_B(ξ) − x_A(ξ))`, con `x_A`, `x_B` interpoladas con las `N_i` (familia §3). Con nodos coincidentes, 0.
- Medida `dΓ = J·e·dξ`, con `e` el espesor de los sólidos adyacentes (familia §4 y §6).

### 3. Salto, fuerzas y rigidez

Sin cambios de fórmula (familia §5–§6): `g = g_0·e_1 + R(ξ)·Σ_i N_i(ξ)·(u_B,par(i) − u_A,i)`, `f_e = Σ_g m_g·B_gᵀ·t_g`, `K_e = Σ_g m_g·B_gᵀ·T_g·B_g`, con `m_g` la medida del punto.

### 4. Integración (3 puntos en las dos reglas)

- **Nodal** (Newton–Cotes / Lobatto de 3 puntos, en los nodos `ξ = −1, +1, 0`): cada punto ve un solo par de nodos (contacto nodo a nodo). En un lado recto, con el intermedio centrado o no, el peso de cada nodo es el área tributaria consistente `m_i = e·∫N_i J dξ` (decisión del usuario del 2026-11-26 para la familia, Diálogo Q3 de `Interface2D4`); con el intermedio en el centro, `m = e·L·(1/6, 1/6, 2/3)` para los nodos (1, 2, 5) (*calculado*): la regla de Simpson. **En cualquier lado, recto o curvo: Simpson con el jacobiano del nodo**, `m_i = e·w_i·J(ξ_i)` con `w = (1/3, 1/3, 4/3)` (decisión del usuario del 2026-11-26, Diálogo Q1). En un lado recto es lo mismo que `∫N_i J`; en uno curvo es la norma del vector de área tributaria `∫N_i n J dξ` y reproduce una presión uniforme.
- **Gauss**: 3 puntos (`ξ = 0, ±√(3/5)`, pesos 8/9, 5/9, 5/9). En un lado recto el integrando de `K_e` con ley lineal es de grado 4 en `ξ` y la regla es exacta hasta grado 5.

Por omisión, `nodal` (decisión del usuario del 2026-11-26 para la familia).

Literatura (*leído en la fuente* por el bibliotecario): con Gauss, un gradiente grande de tracción dentro del elemento y rigidez alta hacen oscilar la tracción entre puntos; Newton–Cotes/Lobatto de 3 puntos no acopla los pares y lo evita (Schellekens y de Borst 1993, pp. 57–58; de Borst et al. 2012, p. 432). En la línea cuadrática Newton–Cotes, Lobatto y la integración concentrada (*lumping*) dan los mismos pesos 1/6, 2/3, 1/6 de la longitud (Schellekens y de Borst 1993, Fig. 2(b)). DIANA (CL12I) integra con Newton–Cotes de 3 puntos por omisión; ANSYS (INTER203), con 2 puntos de Gauss.

### 5. Generación desde la malla (ADR 0043 §6)

`pair` y `split` crean un `Interface2D6` en cada lado de 3 nodos de un sólido cuadrático (Quad8, Quad9, Tri6; `EDGE_NODES` = extremo, intermedio, extremo). La cara A se recorre en sentido contrario al del lado en el contorno antihorario del sólido A (familia §2), con el intermedio en la posición 5. En `pair`, cada nodo intermedio de A se empareja como los extremos (sobre la normal de A en ese nodo). Con la integración nodal, la cara A debe ser recta o casi: con pares a distancia sobre una cara curva el par no queda sobre la normal (familia §2, «qué cara es A»).

## Advertencias y límites (además de los de la familia)

- **Cara A curva con holgura**: igual que en la familia, la cara A con pares a distancia debe ser la recta.
- **Intermedio descentrado**: si el nodo 5 no está en el centro del lado, `J` varía y los pesos nodales dejan de ser 1/6, 1/6, 2/3 (en lado recto, `m = (L/6 + 2δ/3, L/6 − 2δ/3, 2L/3)` con el descentramiento `δ`, *calculado* por el teórico). **Cara degenerada** (corregido por el teórico: `J = ‖a‖` nunca es negativo, así que «`J ≤ 0`» no sirve): es error al construir si algún peso nodal no es positivo o si la tangente de la cara se invierte (`a(±1)·a(0) ≤ 0`). El nodo a un cuarto del elemento singular de Barsoum da peso nulo en el extremo (`m = (0, 1/3, 2/3)·L`) y también es error.
- **3D cuadrático** (caras de 6 y 8 nodos): fuera de este bloque. En la cara serendípita de 8 nodos la fracción consistente de una presión uniforme es −1/12 en las esquinas, y 0 en las esquinas del triángulo de 6 nodos (*calculado*); con integración nodal esos nodos no transmitirían el contacto. Necesita una decisión del usuario sobre la integración.

## Contrato de implementación

```yaml
name: Interface2D6
kind: element
status: implemented

interface:
  field: interface
  dof_names: [ux, uy]
  n_nodes: 6
  jump_dim: 2
  n_integration_points: 3       # nodal o Gauss, los dos con 3 puntos

admission:                      # evidencia en informes/bibliotecario/busqueda_interfaz_cuadratica_2026-11-26.md
  textbooks:
    - de_borst_crisfield: "§13.1 Interface Elements (pp. 428–436): formulación para n pares de nodos (ecs. 13.1–13.10); interpolación cuadrática y Newton–Cotes (pp. 432–433); ejemplo con interfaces de seis nodos y Newton–Cotes de 3 puntos (p. 434)"
    - zienkiewicz_taylor: "vol. 2, 7.ª ed., §4.11.2 elementos de junta isoparamétricos (p. 128) [genérico]"
  general_codes:
    - diana: "CL12I (línea 2D, 3+3 nodos, 2.º orden; Newton–Cotes de 3 puntos por omisión), Theory §33.7.5"
    - ansys: "INTER203 (interfaz 2D de 6 nodos, cuadrática; 2 puntos, Gauss), Element Reference y Theory Reference"
    - code_aster: "PLAN_JOINT con malla QUAD8 (V6.03.162); R3.06.09 §1.2; PLAN_INTERFACE QUAD8 (U3.13.14, R3.06.13)"

parameters:
  - { name: material,    type: int,   required: true,  desc: "Id en interface_materials; ley de superficie con 2 componentes" }
  - { name: integration, type: str,   required: false, desc: "'nodal' (por omisión) | 'gauss'" }
  - { name: thickness,   type: float, required: false, desc: "Espesor fuera del plano; lo pone el generador desde los sólidos adyacentes" }

conventions:
  node_orientation: "A = (1, 2, 5), B = (3, 4, 6) con 4 frente a 1, 3 frente a 2 y 6 frente a 5; esquinas primero (VTK 30, quad6)"

out_of_scope:
  - "Caras cuadráticas 3D (6 y 8 nodos): decisión de integración pendiente"

acceptance:
  verification:
    - name: parcela_de_contacto_cuadratica
      setup: "dos bloques de Quad8 (y de Quad9 y Tri6) con mallas coincidentes en la interfaz, materiales distintos, presión uniforme p, ContactPenalty real (unilateral: hace visible el signo de la normal) sin fricción, K_n de 10 a 100·E/h, espesor ≠ 1; nodal y Gauss; creados con pair y con split"
      expect: "t_n = −p en todos los puntos; σ_yy = −p uniforme; interpenetración p/K_n"
      tol_rel: 1.0e-9
    - name: salto_uniforme_cuadratico
      setup: "una interfaz aislada recta e inclinada de longitud L; cara A fija, cara B trasladada; ley lineal con K_n ≠ K_s"
      expect: "fuerzas nodales e·L·(1/6, 1/6, 2/3)·t en (1, 2, 5) con nodal y con Gauss (detecta pesos de Gauss desordenados)"
      tol_rel: 1.0e-12
  specific:
    - name: cuerpo_rigido_cara_curva
      setup: "cara curva (intermedio fuera de la cuerda), giro rígido pequeño (obligatorio: la traslación sola no detecta pares mal asignados) con nodos coincidentes"
      expect: "g = 0 y f_e = 0"
      tol_rel: 1.0e-12
    - name: ejes_en_cara_curva
      setup: "cara parabólica conocida"
      expect: "n y s en los tres puntos iguales a los de a(ξ) = Σ N_i'·x_i, con signo (detecta la normal de la cuerda y el signo cambiado)"
      tol_rel: 1.0e-12
    - name: pesos_en_cara_curva
      setup: "cara en arco y con intermedio descentrado; salto normal local uniforme"
      expect: "pesos y fuerzas nodales según la decisión Q1 (con Simpson: e·t·e_z × (a(−1)/3, a(+1)/3, 4a(0)/3)); con ∫N_i J, cuadratura fina de 24 puntos"
      tol_rel: 1.0e-12
    - name: cara_degenerada
      setup: "intermedio que da un peso nodal no positivo, o la tangente invertida; nodo a un cuarto"
      expect: "error al construir"
    - name: rango_de_la_rigidez
      setup: "ley lineal, cara curva, nodal y Gauss"
      expect: "6 valores propios no nulos (detecta una regla de 2 puntos)"
      tol_rel: 1.0e-10
    - name: pruebas_de_la_familia
      setup: "las de Interface2D4 repetidas para este elemento: nodal desacopla pares, escalado con el espesor, simetría de la rigidez, orientación invertida (también en cara curva), cohesiva con holgura rechazada"
      expect: "como en la familia"
    - name: holgura_inicial_cuadratica
      setup: "cara A recta y cara B paralela a distancia g_0 (los tres pares sobre la normal); ContactPenalty"
      expect: "g_0 en los tres puntos; cierre con desplazamiento impuesto como en la familia"
      tol_rel: 1.0e-12
    - name: tangente_del_elemento_cuadratico
      setup: "diferencias centradas con ContactPenalty con fricción (adherencia y deslizamiento) en una cara curva"
      expect: "coinciden"
      tol_rel: 1.0e-6

references:
  - "Spec Interface2D4 (formulación de la familia); ADR 0043."
  - "Schellekens, J. C. J., de Borst, R. (1993). On the numerical integration of interface elements. IJNME 36(1), 43–66. DOI 10.1002/nme.1620360104 (docs/referencias/Schellekens_deBorst_1993_integracion_interfaces.pdf). Línea cuadrática: Newton–Cotes = Lobatto = lumping = 1/6, 2/3, 1/6 de la longitud (Fig. 2(b), p. 49); sin acople entre pares, a diferencia de Gauss (Fig. 13, p. 58)."
  - "de Borst, R., Crisfield, M. A., Remmers, J. J. C., Verhoosel, C. V. (2012). Nonlinear Finite Element Analysis of Solids and Structures, 2.ª ed., §13.1. Wiley."
```

## Implementación

- Archivo: `solidum/elements/interface/interface_2d.py` (clase); lo común en `solidum/elements/interface/_shared.py` (funciones de forma `line3`, reglas de 3 puntos, rechazo de la cara degenerada). La base de la familia pasó de `B_PARTNERS` a `A_NODES`/`B_NODES` (posiciones de la cara A en el orden de sus funciones de forma y de su par en B) para admitir el orden «esquinas primero»; los tres elementos lineales dan lo mismo (pruebas de la familia y barrido verdes).
- Clase: `Interface2D6` (`FACE_CELL = "line3"`, `A_NODES = (0, 1, 4)`, `B_NODES = (3, 2, 5)`, `VTK_CELL_TYPE = "quad6"`).
- Salida (paso 5, ADR 0043 §7): meshio 5.3.5 conoce el tipo 30 pero no su dimensión ni su número de nodos y no puede escribir ni leer la celda; `ensure_meshio_cells()` de `solidum/utils/vtk_exporter.py` le añade las dos entradas. Ida y vuelta con el orden del elemento en `tests/test_interface_output.py`.
- Generador (`solidum/mesh/interfaces.py`): los lados de 3 nodos de Quad8, Quad9 y Tri6 dan `Interface2D6`, con el intermedio al final de la cara.
- Tests:
  - `tests/test_interface_2d6.py` (30): todos los casos de `acceptance`, con el convenio de nodos escrito a mano en la prueba (no leído de la clase); parcela exacta desde el YAML con `split` para Quad8, Quad9 y Tri6, nodal y Gauss (tolerancia 1e-9), y con `pair` en Python.
  - `tests/test_interface_generation.py`: `split` con sólidos cuadráticos y el error de las caras cuadráticas 3D.
  - Barrido de contrato con la fábrica de `Interface2D6`.
  - Mutaciones: 11/11 detectadas (`informes/prototipos_2026-11-26/contacto/mutar_interface2d6.py`), con los errores típicos de la tabla del teórico.

## Diálogo

- **2026-11-26 — Origen**: el mensaje de error de `split` con Quad8 recomendaba «use elementos sólidos lineales»; el usuario lo objetó («a mayor grado del elemento es mejor su solución»). La recomendación se quitó (Reglas §3: elegir el elemento es del usuario) y el usuario decidió añadir la interfaz cuadrática 2D a este bloque (opción 1, «si la 1»); la 3D cuadrática queda para cuando haya un caso, por la pregunta de integración de §«Advertencias».
- Integración nodal por omisión: heredada de la familia (Diálogo Q1 de `Interface2D4`).
- **2026-11-26 — Admisión** (`bibliotecario-solidum`, `informes/bibliotecario/busqueda_interfaz_cuadratica_2026-11-26.md`): **estándar** (de Borst et al. 2012 §13.1; DIANA CL12I, ANSYS INTER203, Code_Aster PLAN_JOINT QUAD8; Abaqus no tiene cohesivo cuadrático). El orden «esquinas primero» coincide con ANSYS y Code_Aster; DIANA numera cara por cara. Dato sobre Q1 (*leído*): Schellekens y de Borst (1993, p. 65) observan, con elementos planos lineales distorsionados, que Newton–Cotes calcula mal `det J` en los puntos y el *lumping* `∫N_i dA` no; trata superficies lineales, no líneas cuadráticas curvas, y no contradice el cálculo del teórico para este elemento (en lado recto las dos reglas coinciden).
- **2026-11-26 — Revisión de la spec por `teorico-solidum`** (`informes/teorico/revision_spec_interface2d6_2026-11-26.md`, «cambios antes de programar»): formulación, orden de nodos y las dos reglas en lado recto correctos; parcela exacta con Quad8, Quad9 y Tri6 (≤ 1.4e-14, *calculado*). Aplicado: tipo VTK 30, criterio de cara degenerada por peso no positivo o tangente invertida (y nodo a un cuarto), pruebas nuevas y ajustadas. Hallazgo → Q1.
- **Q1 — DECIDIDO por el usuario el 2026-11-26 («tu recomendación»): opción (b), Simpson con el jacobiano del nodo.** Pregunta: pesos nodales en un lado curvo. En lado recto las dos opciones coinciden con la decisión Q3 de la familia. En lado curvo: (a) `e·∫N_i J dξ` (la letra de Q3): bajo presión uniforme la tracción nodal se desvía +1.3 % en los extremos y −0.7 % en el intermedio con arcos de 30° por elemento (+0.15 % / −0.08 % con 10°), no baja al subir `K_n` y desaparece al refinar; (b) Simpson con el jacobiano del nodo `e·w_i·J(ξ_i)`, `w = (1/3, 1/3, 4/3)`: reproduce la presión uniforme (2e-8 con `K_n` grande), porque es la norma del vector de área tributaria `∫N_i n J dξ` que equilibra a los sólidos (*calculado* por el teórico). Recomendación del teórico y de la IA: (b).
