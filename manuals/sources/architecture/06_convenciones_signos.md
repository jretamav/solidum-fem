# Convenciones de signos

La convención de signos es única para el conjunto del proyecto y de aplicación obligatoria en toda formulación, prueba, interfaz pública (`internal_forces`, `SolveResult`), documentación y catálogo. Si una referencia bibliográfica adopta otra convención, la traducción se realiza al implementar y se anota en la prueba correspondiente.

## Ejes y giros

**Estática 2D.** Eje `x` positivo a la derecha; eje `y` positivo hacia arriba. Giro y momento positivos en sentido antihorario (regla de la mano derecha (RHR) con `z` saliente del plano).

**Estática 3D.** Aplicación de la RHR a todos los ejes (locales y globales) y a todos los momentos (`Mx`, `My`, `Mz`, `T`).

## Esfuerzos internos en 2D — convención de viga estructural

Sobre un elemento diferencial, con `N` positivo en tracción, `V` positivo cuando tiende a rotar el diferencial en sentido horario y `M` positivo en flexión sagitada (tracción en la fibra inferior con `+y` hacia arriba):

- *Cara izquierda* (normal saliente en `−x`): `N` apunta en `−x`, `V` apunta en `+y`, `M` actúa en sentido horario.
- *Cara derecha* (normal saliente en `+x`): `N` apunta en `+x`, `V` apunta en `−y`, `M` actúa en sentido antihorario.

Esta es la convención clásica de los diagramas de esfuerzos en vigas.

## Esfuerzos internos en 3D — convención de resultantes de tensión bajo RHR

En la cara con normal saliente `+x_local`, los esfuerzos positivos se definen como:

- `N` en `+x_local` (tracción positiva).
- `Vy` en `+y_local`, `Vz` en `+z_local`.
- `T ≡ Mx`, `My`, `Mz`: vectores de momento en `+x_local`, `+y_local`, `+z_local` respectivamente, conforme a la RHR.

En la cara con normal saliente `−x_local`, todos los sentidos se invierten (tercera ley de Newton).

Justificación: esta es la convención que resulta de integrar directamente el tensor de tensiones sobre la sección. Es la adoptada por Bathe, Crisfield, Cook-Malkus-Plesha, SAP2000, OpenSees, ANSYS y Abaqus para vigas 3D. La convención estructural de "flexión sagitada positiva" no admite extensión canónica a 3D y se utiliza únicamente en 2D.

Convención de fibra para los flectores en 3D, consecuencia directa del signo:

- `Mz > 0` ⇒ tracción en fibras con `y < 0` (equivalente a flexión sagitada en el plano `xy` con `+y` hacia arriba).
- `My > 0` ⇒ tracción en fibras con `z > 0`.

## Capa interna y capa pública

Internamente, todos los elementos (2D y 3D) operan en la convención de resultantes de tensión bajo RHR. Las matrices `B`, las fuerzas internas `f_int = ∫Bᵀσ dΩ`, los jacobianos y los residuos se calculan en dicha convención. Esta uniformidad evita signos especiales por dimensión dentro del código de formulación.

En la interfaz pública (`internal_forces`, diagramas, catálogo), los elementos 3D la exponen sin transformación. Los elementos 2D, en cambio, exponen la convención de viga estructural (`V` con signo opuesto al interno) por ser la representación clásica en los diagramas de esfuerzos. La traducción consiste en un simple cambio de signo en `V`, aplicado dentro del método `internal_forces()` del elemento 2D.

Relación 2D ↔ 3D en la interfaz pública:

- `M_2D` ≡ `Mz_3D` (mismo signo; ambos corresponden a flexión sagitada positiva con `+y` hacia arriba).
- `V_2D` y `Vy_3D` difieren en signo en la interfaz pública por la razón anterior.

Los valores se exponen siempre en ejes locales del elemento; la transformación a ejes globales constituye una capa separada.
