# Referencias bibliográficas externas

Las referencias canónicas usadas por Solidum FEM no se redistribuyen en este repositorio.
Si una decisión arquitectural o una formulación apunta a una de ellas, se cita aquí
con la fuente oficial donde obtenerla.

## Retama (2010) — Discontinuidades interiores y daño cohesivo

- **Autor**: J. Retama Velasco
- **Tipo**: Tesis Doctoral en Ingeniería Civil
- **Institución**: Universidad Nacional Autónoma de México (UNAM), Programa de Maestría y Doctorado en Ingeniería
- **Año**: 2010
- **Contenido relevante para Solidum**: formulación KOS / KSON, modelos cohesivos $t$–$[[u]]$, condensación estática del salto, expresión `l_d = (A_e/h)·cos(θ−α)` (Cap. 6), criterio de activación Rankine.
- **Vinculado a**: ADR 0010 (discontinuidades interiores embebidas), `solidum/elements/cst_embedded2d.py`, `solidum/materials/cohesive_damage_isotropic.py`.
- **Fuente oficial**: Repositorio de Tesis Digitales UNAM, [tesiunam.dgb.unam.mx](https://tesiunam.dgb.unam.mx/) (búsqueda por autor).

## de Borst & Sluys (1999) — Métodos computacionales en mecánica de sólidos no lineal

- **Autores**: René de Borst, Lambertus J. Sluys
- **Tipo**: Apuntes de curso / monografía académica
- **Institución**: Delft University of Technology, Faculty of Civil Engineering and Geosciences
- **Año**: 1999 (con ediciones posteriores)
- **Contenido relevante para Solidum**: solvers no lineales (Newton-Raphson, arc-length, line search), plasticidad clásica y de superficies múltiples, daño isótropo y anisótropo, criterios de bifurcación. La Tabla 3.1 (§3.4) sustenta la enmienda al ADR 0011 sobre line search con tangente consistente.
- **Vinculado a**: ADR 0011 (robustez de solvers no lineales), `solidum/math/solvers/`, varios materiales.
- **Fuente oficial**: distribuido por TU Delft como material académico; existen versiones publicadas posteriores por los mismos autores como referencia académica equivalente (p. ej. de Borst, Crisfield, Remmers, Verhoosel, *Nonlinear Finite Element Analysis of Solids and Structures*, 2nd ed., Wiley, 2012, ISBN 978-0-470-66644-9, que retoma y extiende el contenido).

## Otras referencias citadas en ADRs y specs

Los ADRs (`docs/adr/000N-titulo.md`) y las specs (`docs/specs/`) citan referencias bibliográficas adicionales en línea. Cada cita lleva información suficiente (autor, año, revista o editorial) para localizar la fuente original a través de buscadores académicos (Google Scholar, Scopus, Web of Science) o catálogos institucionales. No se mantiene un listado consolidado de todas ellas para evitar duplicación con las citas en contexto.
