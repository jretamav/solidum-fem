# Referencias bibliográficas externas

Las referencias canónicas usadas por Solidum FEM no se redistribuyen en este repositorio.
Si una decisión arquitectural o una formulación apunta a una de ellas, se cita aquí
con la fuente oficial donde obtenerla.

## Retama (2010) — Discontinuidades interiores y daño cohesivo

- **Autor**: J. Retama-Velasco
- **Tipo**: Tesis Doctoral en Ingeniería Civil
- **Institución**: Universidad Nacional Autónoma de México (UNAM), Programa de Maestría y Doctorado en Ingeniería
- **Año**: 2010
- **Contenido relevante para Solidum**: formulación KOS / KSON, modelos cohesivos $t$–$[[u]]$, condensación estática del salto, expresión `l_d = (A_e/h)·cos(θ−α)` (Cap. 6), criterio de activación Rankine.
- **Vinculado a**: ADR 0010 (discontinuidades interiores embebidas) y ADR 0020 (módulos de usuario): la formulación es el módulo de usuario `discontinuities`, fuera del programa principal — [`solidum/user/discontinuities/embedded_cst.py`](../../solidum/user/discontinuities/embedded_cst.py) (`CST_Embedded2D`), sus specs en [`docs/user/discontinuities/specs/`](../user/discontinuities/specs/); la ley cohesiva `CohesiveDamageIsotropic` pasó el 2026-11-10 al programa principal ([`solidum/materials/interface/cohesive_damage_isotropic.py`](../../solidum/materials/interface/cohesive_damage_isotropic.py), ADR 0033).
- **Fuente oficial**: Repositorio de Tesis Digitales UNAM, [tesiunam.dgb.unam.mx](https://tesiunam.dgb.unam.mx/) (búsqueda por autor).

## de Borst & Sluys (1999) — Métodos computacionales en mecánica de sólidos no lineal

- **Autores**: René de Borst, Lambertus J. Sluys
- **Tipo**: Apuntes de curso / monografía académica
- **Institución**: Delft University of Technology, Faculty of Civil Engineering and Geosciences
- **Año**: 1999 (con ediciones posteriores)
- **Contenido relevante para Solidum**: solvers no lineales (Newton-Raphson, arc-length, line search), plasticidad clásica y de superficies múltiples, daño isótropo y anisótropo, criterios de bifurcación. La Tabla 3.1 (§3.4) sustenta la enmienda al ADR 0011 sobre line search con tangente consistente.
- **Vinculado a**: ADR 0011 (robustez de solvers no lineales), `solidum/math/solvers/`, varios materiales.
- **Fuente oficial**: distribuido por TU Delft como material académico; existen versiones publicadas posteriores por los mismos autores como referencia académica equivalente (p. ej. de Borst, Crisfield, Remmers, Verhoosel, *Nonlinear Finite Element Analysis of Solids and Structures*, 2nd ed., Wiley, 2012, ISBN 978-0-470-66644-9, que retoma y extiende el contenido).

## Oliver, Huespe & Cante (2008) — Integración IMPL-EX

- **Autores**: J. Oliver, A. E. Huespe, J. C. Cante
- **Tipo**: artículo, *Computer Methods in Applied Mechanics and Engineering* 197(21–24), 1865–1889, DOI [10.1016/j.cma.2007.11.027](https://doi.org/10.1016/j.cma.2007.11.027)
- **Contenido relevante para Solidum**: esquema implícito-explícito para modelos de daño y plasticidad (Tabla 3: extrapolación de la variable interna, ec. 35; esfuerzo y tangente secante, ecs. 38–39), análisis del error (§3.5) y control adaptativo del paso (§3.6, ec. 78).
- **Vinculado a**: specs `IsotropicDamageCrackBand2D/3D` (§8 y Diálogo, con las dos diferencias de la implementación de Solidum aprobadas por el usuario) y ADR 0034 (pseudo-tiempo de los solvers incrementales).
- **Copia local**: `Oliver_Huespe_Cante_2008_CMAME_IMPLEX.pdf` (no versionada).

## Cervera, Chiumenti, Codina et al. (2006–2015) — Localización de deformaciones con fisuración distribuida

Serie sobre el sesgo de malla de los modelos de fisuración distribuida en triángulos lineales y sus remedios; consultada al diagnosticar la banda de grieta (`IsotropicDamageCrackBand2D`) en la viga SEN (2026-11-11). DOIs verificados en Crossref; copias locales no versionadas.

- Cervera, M., Chiumenti, M. (2006). Smeared crack approach: back to the original track. *Int. J. Numer. Anal. Meth. Geomech.* 30, 1173–1199. DOI [10.1002/nag.518](https://doi.org/10.1002/nag.518) — fisuración distribuida guiada por seguimiento de grieta. `Cervera_Chiumenti_2006_IJNAMG_smeared_crack_tracking.pdf`
- Cervera, M., Chiumenti, M., Codina, R. (2010). Mixed stabilized finite element methods in nonlinear solid mechanics. Part I: Formulation. *CMAME* 199, 2559–2570. DOI [10.1016/j.cma.2010.04.006](https://doi.org/10.1016/j.cma.2010.04.006). `Cervera_Chiumenti_Codina_2010_CMAME_mixed_I_formulation.pdf`
- Cervera, M., Chiumenti, M., Codina, R. (2010). … Part II: Strain localization. *CMAME* 199, 2571–2589. DOI [10.1016/j.cma.2010.04.005](https://doi.org/10.1016/j.cma.2010.04.005) — daño isótropo de Rankine con formulación mixta desplazamiento–deformación estabilizada en triángulos lineales. `Cervera_Chiumenti_Codina_2010_CMAME_mixed_II_localization.pdf`
- Cervera, M., Chiumenti, M., Codina, R. (2011). Mesh objective modeling of cracks using continuous linear strain and displacement interpolations. *Int. J. Numer. Meth. Engng* 87, 962–987. DOI [10.1002/nme.3148](https://doi.org/10.1002/nme.3148). `Cervera_Chiumenti_Codina_2011_IJNME_mesh_objective_cracks.pdf`
- Cervera, M., Chiumenti, M., Benedetti, L., Codina, R. (2015). … Part III: Compressible and incompressible plasticity. *CMAME* 285, 752–775. DOI [10.1016/j.cma.2014.11.040](https://doi.org/10.1016/j.cma.2014.11.040). `Cervera_Chiumenti_Benedetti_Codina_2015_CMAME_mixed_III_plasticity.pdf`

## Otras referencias citadas en ADRs y specs

Los ADRs (`docs/adr/000N-titulo.md`) y las specs (`docs/specs/` en el programa principal, `docs/user/<módulo>/specs/` en los módulos de usuario) citan referencias bibliográficas adicionales en línea. Cada cita lleva información suficiente (autor, año, revista o editorial) para localizar la fuente original a través de buscadores académicos (Google Scholar, Scopus, Web of Science) o catálogos institucionales. No se mantiene un listado consolidado de todas ellas para evitar duplicación con las citas en contexto.
