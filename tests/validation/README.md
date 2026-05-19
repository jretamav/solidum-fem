# Benchmarks de validación externa — Fenix FEM

Tests que validan Fenix contra **soluciones analíticas clásicas** y
**benchmarks publicados** (NAFEMS, MacNeal-Harder, Timoshenko-Goodier).
Complementan los tests unitarios y de integración del directorio padre
con cifras citables contra referencias externas.

## Cobertura actual

| # | Benchmark | Archivo | Elementos / solver | Resultado |
|---|-----------|---------|--------------------|-----------|
| 1 | Cilindro grueso de Lamé (plane strain) | `test_lame_cylinder.py` | Quad4, Tri3, Quad8, Quad9, Tri6 | 12/12 |
| 2 | NAFEMS LE1 — elliptic membrane (plane stress) | `test_nafems_le1.py` | Quad4, Tri3, Quad8, Quad9, Tri6 | 10/10 |
| 3 | MacNeal-Harder slender cantilever | `test_slender_cantilever.py` | Frame2DEuler, Frame2DTimoshenko, Quad4, Quad8 | 8/8 |
| 4 | Bathe wave propagation (barra 1D) | `test_bathe_wave_propagation.py` | Truss2D + CentralDifferenceSolver | 4/4 |
| 5 | Cilindro grueso de Hill — J2 perfecta (plane strain) | `test_hill_cylinder_j2.py` | Quad4 + VonMises2D + NonlinearSolver | 4/4 |

**Total: 38/38 tests verde**. Última actualización: 2026-05-19.

## Resultados cuantitativos por benchmark

### 1. Cilindro de Lamé (plane strain)

Cuadrante de corona circular Rᵢ=1, Rₑ=2, presión interna p=1, E=1000, ν=0.3.

Error L²-relativo de los esfuerzos σ_rr y σ_θθ en Gauss interiores:

| Elemento | Malla | e_L²(σ_rr) | e_L²(σ_θθ) |
|----------|-------|-----------|-----------|
| Quad4 | 8×8   | 17.8%     | 2.7%      |
| Quad4 | 16×16 |  8.8%     | 1.3%      |
| Tri3  | 8×8   | 17.1%     | 8.2%      |
| Quad8 | 4×4   |  3.8%     | 1.2%      |
| Quad8 | 8×8   |  0.7%     | 0.3%      |
| Quad9 | 4×4   |  3.5%     | 1.1%      |
| Tri6  | 4×4   |  3.4%     | 1.3%      |
| Tri6  | 8×8   |  0.9%     | 0.4%      |

- **Convergencia O(h¹)** verificada para Quad4/Tri3 (σ ∈ P0 por elemento).
- **Convergencia O(h²)** verificada para Quad8 (ratio 4×4→8×8 = 5.1).
- **Anisotropía radial vs circunferencial** capturada correctamente —
  habría destapado errores de signo en B/C.

### 2. NAFEMS LE1 — elliptic membrane

Cuadrante elíptico entre elipses concéntricas (2,1) y (3.25,2.75), presión
exterior 10 MPa. Valor de referencia: **σ_yy(D=(2,0)) = 92.7 MPa**.

σ_yy en el Gauss más cercano a D (sin stress recovery sofisticado):

| Elemento | Malla   | σ_yy(D, nearest Gauss) | Error |
|----------|---------|------------------------|-------|
| Quad4    | 32×32   | 90.1                   |  2.8% |
| Quad8    | 8×8     | 85.4                   |  7.9% |
| Quad9    | 8×8     | 85.8                   |  7.4% |
| Tri3     | 32×32   | 77.6                   | 16.3% |
| Tri6     | 16×16   | 78.4                   | 15.4% |

Convergencia h monótona verificada para los 5 elementos.

**Limitación documentada**: para alcanzar el valor canónico 92.7 con malla
coarse se necesita **recovery superconvergente (Zienkiewicz SPR)** no
implementado en el post-proceso. Los Tri tienen Gauss en el centroide, lo
que aleja inherentemente el "best sampling point" del nodo de esquina D.

### 3. MacNeal-Harder slender cantilever

Cantilever L=6, h=0.2, t=0.1 (esbeltez L/h=30), E=10⁷, ν=0.3, P_tip=1.
Referencias analíticas:

- **Euler-Bernoulli**: u_tip = PL³/(3EI) = 0.10800
- **Timoshenko**:      u_tip = 0.10809 (corrección de cortante ~0.09%)

| Caso                            | Resultado | Error |
|---------------------------------|-----------|-------|
| Frame2DEuler 1 elemento         | exacto    | < 1e-12 |
| Frame2DTimoshenko 1 elemento    | exacto    | < 1%    |
| Frame2DTimoshenko 8 elementos   | exacto    | < 1e-6  |
| Quad4 30×1 (shear locking severo) | 0.0729 | -32.6%  |
| Quad4 60×4                      | 0.0979    | -9.4%   |
| Quad4 120×8                     | 0.1053    | -2.6%   |
| Quad8 12×1                      | 0.1073    | -0.7%   |

- **Frames analíticos** dan la solución exacta con 1 elemento.
- **Shear locking en Q4** documentado cuantitativamente — Q4 1-capa
  predice ~67% del valor real. Necesita malla 120×8 para alcanzar error
  <5%. Esta es una limitación conocida del elemento bilinear sin
  reduced-integration / B-bar.
- **Q8 12×1** sin locking severo: -0.7% con malla coarse.

### 4. Bathe wave propagation (barra 1D + central difference)

Barra empotrada-libre L=10, E=1, ρ=1, A=1. Velocidad analítica de onda
c = √(E/ρ) = 1.0. Escalón de fuerza F₀ = 10⁻³ aplicado en x=L en t=0⁺.
Discretización con Truss2D, integración explícita con
`CentralDifferenceSolver` (masa lumped, requisito del integrador).

Velocidad de onda medida vs analítica (Δt = 0.5·CFL crítico):

| N elementos | Δt        | c_numérico | Error |
|-------------|-----------|------------|-------|
| 50          | 0.14142   | 1.01563    | 1.56% |
| 100         | 0.07071   | 1.00645    | 0.65% |
| 200         | 0.03536   | 1.00169    | 0.17% |

- **Tiempo de llegada del frente** al midpoint (t=5) con N=100: error <5%.
- **Velocidad inferida de 3 estaciones (x=0.25L, 0.5L, 0.75L)**: error <5%.
- **Convergencia O(h)** del error de dispersión verificada.

Este benchmark cierra el hueco identificado en `[[project_validacion_fase_b]]`
sobre validación de wave propagation en central difference — ningún test
anterior medía la velocidad de propagación.

### 5. Cilindro grueso de Hill (J2 perfecta, plane strain)

Cuadrante de corona Rᵢ=1, Rₑ=2, presión interna p. Material von Mises
perfecto (σ_y=1, H=0), E=1000, ν=0.3, plane strain.

Presiones críticas analíticas (Hill 1950):

    p_e = (σ_y/√3)·(1 − Rᵢ²/Rₑ²) = 0.4330      (yield onset)
    p_L = (2σ_y/√3)·ln(Rₑ/Rᵢ)    = 0.8005      (límite plástico)

Cuatro tests:

| Test | Régimen | Verifica |
|------|---------|----------|
| `elastic_regime_no_plasticity` | p = 0.30 < p_e | α = 0 en todos los Gauss; σ_rr y σ_θθ coinciden con Lamé puro (e_L² < 20% / 6%) |
| `elastoplastic_regime_yield_zone` | p = 0.70 (r_c ≈ 1.49) | Gauss en r < r_c−margen plastificados (α>0); Gauss en r > r_c+margen elásticos (α=0); σ_θθ − σ_rr ≈ 2σ_y/√3 en zona plástica (err < 5%) |
| `stress_in_plastic_zone` | p = 0.70 | σ_rr en zona plástica coincide con fórmula de Hill (e_L² < 18%) |
| `displacement_grows_approaching_limit` | p ∈ {0.50, 0.65, 0.78} | u_r(Rᵢ) crece monótonamente; ratio (p≈0.975·p_L)/(p<p_e) > 2 |

Este benchmark valida cuantitativamente J2 perfecta en problema con
solución analítica cerrada — el más citable de la familia de plasticidad
de Hill. Cierra hueco identificado en la auditoría de Tanda 14: ningún
test anterior comparaba J2 contra Hill 1950.

## Decisiones de diseño

### Métricas elegidas

- **Lamé**: error L²-relativo sobre Gauss interiores (descartando bordes
  donde la BC discreta introduce error de frontera). Métrica estándar de
  Bathe §4.3, Hughes §3.10.
- **LE1**: σ_yy en el Gauss más cercano a D — limitada por la ausencia
  de stress recovery; ver "Deuda técnica" más abajo.
- **Cantilever**: u_tip — escalar simple, físicamente intuitivo.

### Tolerancias

Calibradas empíricamente con margen. **No son fruto de bajar el listón
hasta que un test pase** — reflejan el orden de convergencia teórico del
elemento (O(h¹) para Q4/Tri3, O(h²) para los cuadráticos). El test de
h-refinement complementa la verificación cuantitativa.

### Convención de mallas estructuradas

Cada benchmark construye su malla en coordenadas paramétricas (s, t) y
mapea al dominio físico. La construcción de mallas con midnodes para
Q8/Q9/Tri6 reusa el patrón establecido en `tests/test_cooks_membrane.py`.

## Deuda técnica identificada

1. **Stress recovery superconvergente (SPR de Zienkiewicz)** no
   implementado. Tendría dos beneficios:
   - σ_yy(D) en NAFEMS LE1 alcanzaría ≈92.7 con malla coarse.
   - Cierres de algunos huecos abiertos de la Fase B (cohesivo G_F).
2. **Q4 con shear locking**: alternativas como Quad4 con integración
   reducida selectiva (SRI) o B-bar permitirían malla coarse en slender
   beam sin locking. La spec de un Quad4-SRI sería un componente nuevo
   (no variante) por su impacto en la rigidez.
3. **Sólidos 3D**: el cimiento de validación es solo 2D. NAFEMS LE10
   (placa gruesa 3D), LE2 (cilindro tapado), LE3 (hemisferio) requieren
   H8/H20/T4/T10 antes de poder validar.

## Cómo añadir un benchmark nuevo

1. Crear `tests/validation/test_<benchmark>.py`.
2. Module-docstring con: referencia bibliográfica completa, definición
   física del problema, BCs, cantidad medida, valor de referencia.
3. Helpers de malla y carga al estilo de los tres existentes.
4. Marcador `@pytest.mark.parametrize` por elemento cuando aplique.
5. **Al menos un test de h-convergencia** además de los tests con malla
   fija (verifica la propiedad cualitativa más fuerte de la formulación).
6. Tolerancias calibradas empíricamente — **no las hagas pasar bajando
   el listón**, ajusta la malla o el método de medida.
7. Actualizar la tabla de cobertura de este README.
