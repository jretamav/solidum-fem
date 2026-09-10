# Tira de bambú con fibra fuera de eje

Demuestra el efecto que **distingue a un material ortótropo de uno isótropo**:
el acoplamiento tracción–cortante. Al traccionar una tira cuya fibra no está
alineada con la carga, la tira no sólo se alarga — se distorsiona.

Primer ejemplo del catálogo que usa [`Orthotropic2D`](../../docs/specs/Orthotropic2D.md),
el material añadido para el proyecto PAPIIT sobre bambú.

---

## Archivos

| Archivo | Qué es |
|---|---|
| `tira_45_grados.yaml` | Cuadrado unitario, malla 4×4 de `Quad4`, fibra a 45°, tracción de 10 MPa |
| `barrido_orientacion.py` | Ejecuta el modelo variando `theta` y tabula la respuesta |

## Cómo ejecutarlo

```bash
python examples/ejecutar_yaml.py examples/bambu_ortotropo/tira_45_grados.yaml
python examples/bambu_ortotropo/barrido_orientacion.py
```

---

## Qué esperar

Desplazamientos de la esquina superior derecha bajo 10 MPa de tracción en `x`:

| θ | u_x [mm] | u_y [mm] | u_y/u_x | E_x [GPa] | E_x/E₁ |
|---|---|---|---|---|---|
| 0°  |  0.667 | −0.233 | **−0.350** | 15.000 | 1.000 |
| 15° |  1.500 | −3.374 | −2.249 |  6.666 | 0.444 |
| 30° |  3.747 | −5.621 | −1.500 |  2.669 | 0.178 |
| 45° |  6.746 | −6.313 | −0.936 |  1.482 | 0.099 |
| 60° |  9.664 | −5.338 | −0.552 |  1.035 | 0.069 |
| 75° | 11.748 | −3.091 | −0.263 |  0.851 | 0.057 |
| 90° | 12.500 | −0.233 | **−0.019** |  0.800 | 0.053 |

Material: `E₁ = 15 GPa`, `E₂ = 0.8 GPa`, `G₁₂ = 0.7 GPa`, `ν₁₂ = 0.35`.

### Tres cosas que leer en esa tabla

**1. En los ejes principales sólo hay Poisson.** En θ=0° el cociente `u_y/u_x`
vale exactamente **−0.350 = −ν₁₂**; en θ=90°, **−0.0187 = −ν₂₁**. La
reciprocidad `ν₁₂/E₁ = ν₂₁/E₂` obliga a que el segundo sea ~19 veces menor —
un Poisson diminuto es lo normal en un material muy ortótropo, no un error.

**2. Fuera de eje aparece algo que no es Poisson.** A 15° el cociente llega a
**−2.25**, muy por encima de cualquier coeficiente de Poisson admisible. No es
contracción transversal: es el **acoplamiento tracción–cortante**, imposible en
un isótropo.

Nótese que el máximo de distorsión *relativa* no cae en 45° sino cerca de 15°,
donde la tira todavía es rígida (`E_x = 6.7 GPa`) y el acoplamiento ya es
fuerte. A 45° el material se ha ablandado tanto que `u_x` crece y el cociente
baja, aunque `u_y` alcance allí su máximo absoluto.

**3. La rigidez se desploma más rápido de lo que sugiere E₁/E₂.** A 45° queda
por debajo del **10 %** de `E₁`, pese a que `E₁/E₂ = 18.75`. Ya a 15° se ha
perdido más de la mitad. El mínimo no se obtiene interpolando entre los dos
módulos principales: el término cruzado `(1/G₁₂ − 2ν₁₂/E₁)` domina la zona
intermedia y, con `G₁₂` pequeño como en bambú, hunde la curva por debajo de la
media armónica de `E₁` y `E₂`.

Para un culmo real esto significa que una desviación de fibra de sólo 15° ya
cuesta más de la mitad de la rigidez axial.

---

## Relación con la validación

Los cocientes `E_x/E₁` de la tabla coinciden con la solución cerrada de
**Jones (1999) §2.8** para las constantes aparentes fuera de eje, verificada a
1e-12 en [`tests/validation/test_off_axis_orthotropic.py`](../../tests/validation/test_off_axis_orthotropic.py).

Este ejemplo llega a las mismas cifras por un camino distinto —modelo FEM
completo con malla, ensamblaje y solver, en vez de la matriz constitutiva
sola—, así que sirve de comprobación cruzada además de demostración.

---

## Advertencias

**Los valores del material son de literatura general**, no de una especie
concreta de bambú. Sustitúyelos por los medidos en la campaña experimental; la
spec anota en su §Diálogo qué queda por determinar (las cuatro constantes de la
especie y si el plano de análisis es *L–R* o *L–T*).

**La orientación es uniforme en toda la malla.** `theta` vive en el material
por decisión del [ADR 0013](../../docs/adr/0013-orientacion-material-y-ortotropia.md),
lo que basta para probetas y tiras pero **no** para un culmo completo, donde los
ejes materiales siguen la geometría cilíndrica. Ese caso requiere orientación
por punto de Gauss y está declarado fuera de alcance con su condición de
retoma.

**Las condiciones de borde están elegidas para no coartar el efecto.** Sólo se
restringe `u_x` en el borde izquierdo, más `u_y` en una esquina para eliminar
el sólido rígido. Empotrar el borde entero impediría la distorsión que el
ejemplo quiere mostrar.
