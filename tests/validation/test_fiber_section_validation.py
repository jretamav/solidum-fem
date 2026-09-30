"""Validación ciega de la sección por fibras de los marcos Euler-Bernoulli.

Materializa los ocho casos del bloque ``acceptance`` de
``docs/specs/FiberSection.md`` contra la API pública fijada por el ADR 0021
(``solidum.sections`` y el argumento ``section`` de ``Frame2DEuler``,
``Frame2DEulerCorot`` y ``Frame3D``).

Test ciego
----------
Escrito desde la spec y la física, **sin leer** la implementación de
``solidum/sections/`` ni el camino con sección de los marcos: es una
comprobación independiente. Todas las cifras esperadas salen de formas
cerradas (Neal 1977) o de identidades exactas del modelo discreto (estática
de los elementos, simetrías de la sección discreta, análisis de límite), nunca
de una ejecución de la implementación que se valida.

Referencias
-----------
- Neal B.G. (1977). *The Plastic Methods of Structural Analysis*, 3.ª ed.:
  M–κ de la sección rectangular elastoplástica perfecta, interacción N–M,
  cargas de colapso del voladizo y de la viga biempotrada.
- Spacone E., Filippou F.C., Taucer F.F. (1996). Fibre beam–column model for
  non-linear analysis of R/C frames: Part I. Formulation. EESD 25:711–725.
- Bathe K.-J. (2014). *Finite Element Procedures*, 2.ª ed., cap. 6.

Datos comunes
-------------
E = 200 GPa, σ_y = 250 MPa, b = 0.3 m (ancho, eje z local), h = 0.5 m
(peralte, eje y local), L = 3 m. Coeficientes distintos de 1 a propósito: los
unitarios esconden errores dimensionales. Para la flexión alrededor de z::

    ε_y = σ_y/E = 1.25e-3               κ_y = 2σ_y/(E h) = 5.0e-3 1/m
    M_y = σ_y b h²/6 = 3.125e6 N·m      M_p = σ_y b h²/4 = 4.6875e6 N·m
    N_p = σ_y b h   = 3.75e7 N

Alrededor de y (peralte b) se intercambian b y h: κ_y = 8.333e-3 1/m,
M_y = 1.875e6 N·m, M_p = 2.8125e6 N·m.

Convenciones (Reglas §5, spec §2): ε_i = ε0 − y_i κ_z + z_i κ_y;
M_z = −Σ y σ A > 0 tracciona las fibras con y < 0; M_y = Σ z σ A > 0
tracciona las fibras con z > 0.

Tolerancias — derivación
------------------------
(a) **Sección rectangular de n capas** (fibra en el centro de cada capa de
espesor Δ = h/n; regla del punto medio). Con N = 0, M = σ_y b ∫ g(y) dy con
g = y²/c en el núcleo elástico |y| < c = ε_y/κ y g = |y| fuera:

  - núcleo: error del punto medio por capa Δ³g''/24 con g'' = 2/c, sumado
    sobre las 2c/Δ capas del núcleo: σ_y b Δ²/6;
  - zona plástica: g lineal en cada capa (con n par ninguna capa cruza
    y = 0): error nulo;
  - las dos capas que contienen la frontera |y| = c: g tiene un quiebre de
    pendiente de salto J = 1 (de 2 a 1); el error del punto medio con un
    quiebre dentro de la celda es a lo sumo J Δ²/8: en total σ_y b Δ²/4.

  Cota: σ_y b Δ² (1/6 + 1/4) = **(5/3) M_p/n²**. En régimen elástico el error
  es exactamente M_y/n² = (2/3) M_p/n² (EI de n capas = b h³/12 (1 − 1/n²)).
  Con **N impuesta** el error de N desplaza el eje neutro; el efecto en M es
  el error del momento respecto del eje neutro (misma cota) más un término de
  segundo orden σ_y b δy_n² con δy_n ≤ Δ/2: cota **(8/3) M_p/n²**.
  Con n = 200: 4.2e-5 y 6.7e-5 de M_p.

(b) **Sección circular** de n_rad × n_circ sectores (fibra en el centroide del
sector): el núcleo aporta Σ A_s Var_s(y)/c ≤ M_p [(Δr/R)²/4 + Δφ²/12]; los
sectores cortados por la frontera elastoplástica, con el quiebre en la peor
posición, M_p (3/8)[(Δr/R)² + Δφ²]. Cota: **M_p [(5/8)(Δr/R)² + (11/24)Δφ²]**
(6.0e-3 de M_p con n_rad = 20, n_circ = 64). Es holgada; la precisión la dan
las comprobaciones exactas del punto (c).

(c) **Identidades exactas** del modelo discreto (reproducción elástica con
n_ip ≥ 2, M_p discreto, simetrías de la sección discreta, estática del
elemento): sólo redondeo y la tolerancia de admisibilidad del retorno plástico
de cada fibra (``ADMISSIBILITY_TOL_REL = 1e-10``): de 1e-12 a 1e-9.

(d) **Colapsos**: razonamiento completo en el bloque de los casos 4 y 5.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

import solidum  # noqa: F401  — autodescubrimiento de los registros
from solidum.core.domain import Domain
from solidum.elements.frame import Frame2DEuler, Frame2DEulerCorot
from solidum.elements.frame3d import Frame3D
from solidum.materials.elastic import Elastic1D
from solidum.materials.plastic_1d import Elastoplastic1D
from solidum.math.assembly import Assembler
from solidum.math.convergence import ConvergenceCriterion
from solidum.math.solvers import NonlinearSolver
from solidum.sections import BarLayer, CircularPatch, Fiber, FiberSection, RectangularPatch


# =============================================================================
# Datos físicos
# =============================================================================

E = 200.0e9            # [Pa]
SIGMA_Y = 250.0e6      # [Pa]
B = 0.3                # ancho de la sección, eje z local [m]
H = 0.5                # peralte de la sección, eje y local [m]
L = 3.0                # longitud [m]
NU = 0.3
EPS_Y = SIGMA_Y / E
G = E / (2.0 * (1.0 + NU))

# Torsión: sólo la exige Frame3D y ningún caso la ejercita (spec §1: elástica
# y desacoplada). J de Saint-Venant aproximado de un rectángulo 0.5 × 0.3.
J_RECT = 0.196 * H * B**3
GJ_RECT = G * J_RECT

N_LAYERS = 200         # capas en la dirección de flexión (M–κ, interacción)

# Identidades exactas: redondeo + admisibilidad del retorno plástico (1e-10).
TOL_EXACT = 1.0e-9

# Elemento de 3D alineado con +x global: con ref_vector = +y global los ejes
# locales coinciden con los globales (ŷ = proyección de ref_vector, ẑ = x̂ × ŷ).
REF_ALIGNED = [0.0, 1.0, 0.0]


# =============================================================================
# Soluciones analíticas
# =============================================================================

def rect_limits(depth: float, width: float):
    """(κ_y, M_y, M_p, N_p) de la sección rectangular con peralte ``depth``."""
    kappa_y = 2.0 * SIGMA_Y / (E * depth)
    M_y = SIGMA_Y * width * depth**2 / 6.0
    M_p = SIGMA_Y * width * depth**2 / 4.0
    N_p = SIGMA_Y * width * depth
    return kappa_y, M_y, M_p, N_p


def m_rect(kappa: float, depth: float, width: float) -> float:
    """M–κ de la sección rectangular elastoplástica perfecta con N = 0 (Neal 1977).

    M = EI κ para κ ≤ κ_y; M = M_p [1 − (κ_y/κ)²/3] para κ ≥ κ_y.
    """
    kappa_y, _, M_p, _ = rect_limits(depth, width)
    if kappa <= kappa_y:
        return E * width * depth**3 / 12.0 * kappa
    return M_p * (1.0 - (kappa_y / kappa) ** 2 / 3.0)


def m_rect_axial(kappa: float, n_ratio: float, depth: float, width: float) -> float:
    """M(κ) con N = n·N_p impuesta y las dos caras plastificadas.

    Núcleo elástico de semiancho c = ε_y/κ centrado en el eje neutro
    y_n = N/(2σ_y b); fuera del núcleo |σ| = σ_y. Integrando::

        N = 2 σ_y b y_n,      M = σ_y b (h²/4 − y_n² − c²/3)
                                = M_p [1 − (N/N_p)² − (κ_y/κ)²/3],

    válida mientras el núcleo no toque las caras, κ ≥ κ_y/(1 − |N|/N_p).
    Para κ → ∞ tiende a la interacción plástica M_p [1 − (N/N_p)²].
    La misma expresión vale para la flexión alrededor de y (peralte b).
    """
    kappa_y, _, M_p, _ = rect_limits(depth, width)
    return M_p * (1.0 - n_ratio**2 - (kappa_y / kappa) ** 2 / 3.0)


def m_circle(kappa: float, radius: float) -> float:
    """M–κ de la sección circular maciza elastoplástica perfecta, N = 0.

    Núcleo elástico |y| < c = ε_y/κ, ancho de la sección 2√(R² − y²)::

        M = (4/3) σ_y (R² − c²)^{3/2}
          + (4σ_y/c) [ (c/8)(2c² − R²)√(R² − c²) + (R⁴/8) arcsen(c/R) ]

    En c = R da M_y = π σ_y d³/32; para c → 0, M_p = σ_y d³/6.
    """
    kappa_y = EPS_Y / radius
    if kappa <= kappa_y:
        return E * math.pi * radius**4 / 4.0 * kappa
    c = EPS_Y / kappa
    root = math.sqrt(radius**2 - c**2)
    return (4.0 / 3.0 * SIGMA_Y * (radius**2 - c**2) ** 1.5
            + 4.0 * SIGMA_Y / c * (c / 8.0 * (2.0 * c**2 - radius**2) * root
                                   + radius**4 / 8.0 * math.asin(c / radius)))


def tol_rect(n: int) -> float:
    """Cota del error de discretización de M con N = 0, relativa a M_p (docstring, a)."""
    return (5.0 / 3.0) / n**2


def tol_rect_axial(n: int) -> float:
    """Cota del error de discretización de M con N impuesta, relativa a M_p (a)."""
    return (8.0 / 3.0) / n**2


def tol_circle(n_rad: int, n_circ: int) -> float:
    """Cota del error de discretización de M de la sección circular, relativa a M_p (b)."""
    dphi = 2.0 * math.pi / n_circ
    return 5.0 / 8.0 / n_rad**2 + 11.0 / 24.0 * dphi**2


# =============================================================================
# Construcción de secciones y elementos
# =============================================================================

def rect_section(material, ny: int, nz: int, depth: float = H, width: float = B,
                 GJ: float | None = None) -> FiberSection:
    """Rectángulo centrado en el eje del elemento, rejilla ny × nz."""
    patch = RectangularPatch(material, y=(-depth / 2.0, depth / 2.0),
                             z=(-width / 2.0, width / 2.0), ny=ny, nz=nz)
    return FiberSection(patches=[patch], GJ=GJ)


def axis_section(axis: str, material, n: int, GJ: float | None = None):
    """Sección rectangular con ``n`` capas en la dirección de flexión del eje.

    ``2d``: flexión en xy (peralte h). ``3d_z``: Frame3D, flexión alrededor de
    z (peralte h, 2 columnas en z). ``3d_y``: flexión alrededor de y (peralte b
    a lo largo de z, 2 filas en y). Devuelve (sección, peralte, ancho).
    """
    if axis == "2d":
        return rect_section(material, n, 1, GJ=GJ), H, B
    if axis == "3d_z":
        return rect_section(material, n, 2, GJ=GJ), H, B
    if axis == "3d_y":
        return rect_section(material, 2, n, GJ=GJ), B, H
    raise ValueError(axis)


def section_response(sec, axis: str, eps0: float, kappa: float, states):
    """(N, M, M_otro, k_NN, estados') con la curvatura ``kappa`` en el eje dado."""
    if axis == "2d":
        s, k, new = sec.response_2d(np.array([eps0, kappa]), states)
        return s[0], s[1], 0.0, k[0, 0], new
    if axis == "3d_z":
        s, k, new = sec.response_3d(np.array([eps0, 0.0, kappa]), states)
        return s[0], s[2], s[1], k[0, 0], new
    s, k, new = sec.response_3d(np.array([eps0, kappa, 0.0]), states)
    return s[0], s[1], s[2], k[0, 0], new


def element_2d(cls, p1, p2, **kwargs):
    """Un marco 2D aislado (sin condiciones de contorno) entre ``p1`` y ``p2``."""
    dom = Domain()
    n1 = dom.add_node(1, list(p1))
    n2 = dom.add_node(2, list(p2))
    elem = cls(1, [n1, n2], **kwargs)
    dom.add_element(elem)
    dom.generate_equation_numbers(verbose=False)
    return elem


def element_3d(p1, p2, **kwargs):
    """Un Frame3D aislado (sin condiciones de contorno) entre ``p1`` y ``p2``."""
    dom = Domain()
    n1 = dom.add_node(1, list(p1))
    n2 = dom.add_node(2, list(p2))
    elem = Frame3D(1, [n1, n2], **kwargs)
    dom.add_element(elem)
    dom.generate_equation_numbers(verbose=False)
    return elem


def uniform_curvature_dofs(axis: str, kappa: float):
    """u_e de curvatura uniforme κ: θ1 = −κL/2, θ2 = +κL/2, desplazamientos nulos.

    Con N2'' = (6x/L − 4)/L y N4'' = (6x/L − 2)/L (spec §5):
    κ(x) = N2''θ1 + N4''θ2 = κ en todo el elemento. En Frame3D, plano xz,
    κ_y = N2''θ_y1 + N4''θ_y2 con w = 0. Con κ uniforme las fuerzas nodales
    son F_θ2 = ∫N4'' M dx = +M y F_θ1 = ∫N2'' M dx = −M, cortante nulo.
    Devuelve (u_e, índice de F_θ2, índice de F_θ1).
    """
    half = kappa * L / 2.0
    if axis == "2d":
        u = np.zeros(6)
        u[2], u[5] = -half, half
        return u, 5, 2
    u = np.zeros(12)
    i1, i2 = (5, 11) if axis == "3d_z" else (4, 10)
    u[i1], u[i2] = -half, half
    return u, i2, i1


def assert_same_response(K_test, F_test, K_ref, F_ref, rtol: float, label: str) -> None:
    """K y F_int iguales a ``rtol`` con escalas dimensionalmente coherentes.

    K se compara normalizada por √(K_ii K_jj): cada entrada queda adimensional
    y O(1), y los términos de flexión (EI/L³) no se esconden detrás de los
    axiales (EA/L). F_int con atol = rtol·‖F_ref‖∞ (fuerzas y momentos del
    mismo orden: los u_e se eligen para que axial y flexión sean comparables).
    El redondeo de las sumas de fibras es ~n_f·eps ≈ 1e-14: cualquier
    diferencia por encima de 1e-10 es de formulación.
    """
    d = np.sqrt(np.abs(np.diag(K_ref)))
    assert np.all(d > 0.0), f"{label}: diagonal de K de referencia con ceros"
    dK = np.abs(K_test - K_ref) / np.outer(d, d)
    assert dK.max() <= rtol, (
        f"{label}: max |ΔK_ij|/√(K_ii K_jj) = {dK.max():.3e} > {rtol:.0e}")
    np.testing.assert_allclose(F_test, F_ref, rtol=rtol, atol=rtol * np.abs(F_ref).max(),
                               err_msg=f"{label}: F_int distinto")


# Secciones elásticas del caso 1 ------------------------------------------------

E_CONCRETE = 30.0e9
A_BAR = 3.1e-4
Y_BAR, Z_BAR = 0.2, 0.1
NY_EL, NZ_EL = 10, 6


def elastic_section(kind: str):
    """Sección elástica y sus rigideces integradas en forma cerrada.

    Rectángulo de ny × nz celdas con la fibra en el centro de cada celda:
    ∫ y² sobre una celda = Δ y_c² + Δ³/12 ⇒ Σ y² A = b h³/12 (1 − 1/ny²),
    Σ z² A = h b³/12 (1 − 1/nz²). ``mixed``: rectángulo de concreto y dos
    capas de 3 barras de acero simétricas (y = ±0.2, z ∈ {−0.1, 0, 0.1});
    rigideces de la sección transformada (sin descontar el concreto que ocupan
    las barras: la sección suma fibras). Simétrica en y y en z ⇒ sin
    acoplamiento N–M. Devuelve (sección, E_ref, EA, EI_z, EI_y, materiales, GJ).
    """
    Iz_n = B * H**3 / 12.0 * (1.0 - 1.0 / NY_EL**2)
    Iy_n = H * B**3 / 12.0 * (1.0 - 1.0 / NZ_EL**2)
    if kind == "rect":
        mat = Elastic1D(E=E)
        sec = rect_section(mat, NY_EL, NZ_EL, GJ=GJ_RECT)
        return sec, E, E * B * H, E * Iz_n, E * Iy_n, [mat], GJ_RECT
    concrete = Elastic1D(E=E_CONCRETE)
    steel = Elastic1D(E=E)
    GJ = E_CONCRETE / (2.0 * (1.0 + NU)) * J_RECT
    sec = FiberSection(
        patches=[
            RectangularPatch(concrete, y=(-H / 2, H / 2), z=(-B / 2, B / 2), ny=NY_EL, nz=NZ_EL),
            BarLayer(steel, start=(-Y_BAR, -Z_BAR), end=(-Y_BAR, Z_BAR), n=3, area=A_BAR),
            BarLayer(steel, start=(Y_BAR, -Z_BAR), end=(Y_BAR, Z_BAR), n=3, area=A_BAR),
        ],
        GJ=GJ,
    )
    EA = E_CONCRETE * B * H + E * 6 * A_BAR
    EIz = E_CONCRETE * Iz_n + E * 6 * A_BAR * Y_BAR**2
    EIy = E_CONCRETE * Iy_n + E * 2 * A_BAR * (Z_BAR**2 + 0.0 + Z_BAR**2)
    return sec, E_CONCRETE, EA, EIz, EIy, [concrete, steel], GJ


# =============================================================================
# Caso 1 — reproducción elástica
# =============================================================================
# Con fibras elásticas y sección simétrica, k_s = diag(EA, EI) (spec §4). El
# integrando de la rigidez del elemento, B_gᵀ k_s B_g, es cuadrático en x
# (B axial constante, curvatura lineal): Gauss-Legendre con n_ip puntos
# integra exactamente grado 2 n_ip − 1, así que n_ip = 2 (el mínimo, el caso
# de la spec) y n_ip = 3 (el valor por defecto) reproducen la forma cerrada.

@pytest.mark.parametrize("ny, nz", [(1, 1), (4, 3), (10, 6), (200, 1)])
def test_caso1_capas_en_el_centro_y_EI_de_n_capas(ny, nz):
    """Fibra en el centro de cada celda; EI de n capas = b h³/12 (1 − 1/n²)."""
    sec = rect_section(Elastic1D(E=E), ny, nz)
    dy, dz = H / ny, B / nz
    assert sec.n_fibers == ny * nz
    iy = (sec.y + H / 2.0) / dy - 0.5
    iz = (sec.z + B / 2.0) / dz - 0.5
    np.testing.assert_allclose(iy, np.round(iy), atol=1e-9, err_msg="y fuera del centro de la celda")
    np.testing.assert_allclose(iz, np.round(iz), atol=1e-9, err_msg="z fuera del centro de la celda")
    cells = {(int(round(a)), int(round(b))) for a, b in zip(iy, iz)}
    assert cells == {(i, j) for i in range(ny) for j in range(nz)}, "celdas repetidas o ausentes"
    np.testing.assert_allclose(sec.area, dy * dz, rtol=1e-14)

    Iz_n = B * H**3 / 12.0 * (1.0 - 1.0 / ny**2)
    Iy_n = H * B**3 / 12.0 * (1.0 - 1.0 / nz**2)
    scale = B * H**3 / 12.0
    np.testing.assert_allclose(np.sum(sec.area * sec.y**2), Iz_n, rtol=1e-12, atol=1e-14 * scale)
    np.testing.assert_allclose(np.sum(sec.area * sec.z**2), Iy_n, rtol=1e-12, atol=1e-14 * scale)
    _, k, _ = sec.response_2d(np.array([1.0e-5, 2.0e-4]), sec.initial_states())
    np.testing.assert_allclose(k[1, 1], E * Iz_n, rtol=1e-12, atol=1e-14 * E * scale,
                               err_msg="EI de la sección ≠ E b h³/12 (1 − 1/n²)")
    np.testing.assert_allclose(k[0, 0], E * B * H, rtol=1e-12)


@pytest.mark.parametrize("kind", ["rect", "mixed"])
def test_caso1_rigidez_de_seccion_elastica(kind):
    """k_s = diag(EA, EI_z) en 2D y diag(EA, EI_y, EI_z) en 3D; s = k_s e."""
    sec, _E_ref, EA, EIz, EIy, mats, GJ = elastic_section(kind)
    assert list(sec.materials()) == mats, "materials() no respeta el orden de aparición"
    assert sec.GJ == GJ
    if kind == "mixed":
        assert sec.n_fibers == NY_EL * NZ_EL + 6
        bars = np.isclose(sec.area, A_BAR, rtol=1e-12)
        assert bars.sum() == 6
        got = sorted(zip(np.round(sec.y[bars], 12), np.round(sec.z[bars], 12)))
        want = sorted((y, z) for y in (-Y_BAR, Y_BAR) for z in (-Z_BAR, 0.0, Z_BAR))
        np.testing.assert_allclose(np.array(got), np.array(want), atol=1e-12,
                                   err_msg="barras de BarLayer mal colocadas (n = 3, extremos incluidos)")

    e2 = np.array([1.0e-5, 2.0e-4])
    s2, k2, _ = sec.response_2d(e2, sec.initial_states())
    np.testing.assert_allclose(np.diag(k2), [EA, EIz], rtol=1e-12)
    assert abs(k2[0, 1]) <= 1e-12 * math.sqrt(EA * EIz)
    np.testing.assert_allclose(k2, k2.T, rtol=0.0, atol=1e-12 * math.sqrt(EA * EIz))
    np.testing.assert_allclose(s2, [EA * e2[0], EIz * e2[1]], rtol=1e-12)

    e3 = np.array([1.0e-5, 1.0e-4, 2.0e-4])
    s3, k3, _ = sec.response_3d(e3, sec.initial_states())
    diag = np.array([EA, EIy, EIz])
    np.testing.assert_allclose(np.diag(k3), diag, rtol=1e-12)
    off = np.abs(k3 - np.diag(np.diag(k3))) / np.sqrt(np.outer(diag, diag))
    assert off.max() <= 1e-12, f"acoplamiento espurio en k_s 3D: {off.max():.2e}"
    np.testing.assert_allclose(s3, diag * e3, rtol=1e-12)


# Elemento oblicuo (L = 3, cosenos 0.6 y 0.8) y u_e con axial y flexión comparables.
P1_2D = (0.4, -0.2)
P2_2D = (0.4 + 0.6 * L, -0.2 + 0.8 * L)
U_2D = np.array([2.0e-5, -3.0e-4, 1.0e-4, -1.0e-5, 2.0e-4, -2.0e-4])


@pytest.mark.parametrize("kind", ["rect", "mixed"])
@pytest.mark.parametrize("n_ip", [2, 3])
@pytest.mark.parametrize("cls", [Frame2DEuler, Frame2DEulerCorot], ids=lambda c: c.__name__)
def test_caso1_reproduccion_elastica_marcos_2d(cls, n_ip, kind):
    """K y F_int del marco con fibras = forma cerrada con EA, EI integrados (1e-10)."""
    sec, E_ref, EA, EIz, _EIy, _mats, _GJ = elastic_section(kind)
    ref = element_2d(cls, P1_2D, P2_2D, material=Elastic1D(E=E_ref), A=EA / E_ref, I=EIz / E_ref)
    fib = element_2d(cls, P1_2D, P2_2D, section=sec, n_ip=n_ip)
    K_r, F_r = ref.compute_element_state(U_2D)
    K_f, F_f = fib.compute_element_state(U_2D)
    assert_same_response(K_f, F_f, K_r, F_r, 1.0e-10, f"{cls.__name__} n_ip={n_ip} {kind}")


# Frame3D oblicuo: dirección (2, 3, 6)/7, ref_vector z global (no paralelo).
_DIR_3D = np.array([2.0, 3.0, 6.0]) / 7.0
P1_3D = (0.5, -0.3, 0.2)
P2_3D = tuple(np.array(P1_3D) + L * _DIR_3D)
U_3D = np.array([2.0e-5, -3.0e-4, 1.5e-4, 2.0e-4, 1.0e-4, -2.0e-4,
                 -1.0e-5, 2.0e-4, -1.0e-4, -1.0e-4, -2.0e-4, 1.5e-4])


@pytest.mark.parametrize("kind", ["rect", "mixed"])
@pytest.mark.parametrize("n_ip", [2, 3])
def test_caso1_reproduccion_elastica_frame3d(n_ip, kind):
    """Frame3D con fibras = forma cerrada con EA, EI_y, EI_z integrados y GJ (1e-10)."""
    sec, E_ref, EA, EIz, EIy, _mats, GJ = elastic_section(kind)
    G_ref = E_ref / (2.0 * (1.0 + NU))
    ref = element_3d(P1_3D, P2_3D, material=Elastic1D(E=E_ref), A=EA / E_ref,
                     Iy=EIy / E_ref, Iz=EIz / E_ref, J=GJ / G_ref, nu=NU,
                     ref_vector=[0.0, 0.0, 1.0])
    fib = element_3d(P1_3D, P2_3D, section=sec, n_ip=n_ip, ref_vector=[0.0, 0.0, 1.0])
    K_r, F_r = ref.compute_element_state(U_3D)
    K_f, F_f = fib.compute_element_state(U_3D)
    assert_same_response(K_f, F_f, K_r, F_r, 1.0e-10, f"Frame3D n_ip={n_ip} {kind}")


# =============================================================================
# Caso 2 — momento–curvatura de la sección rectangular elastoplástica perfecta
# =============================================================================
# Elastoplastic1D con H = 0, N = 0 (ε0 = 0; N resulta nula por simetría),
# curvatura monótona creciente. En 1D el retorno plástico es exacto para
# cualquier tamaño de incremento si la deformación de cada fibra es monótona,
# así que el barrido incremental reproduce la M–κ analítica.
# Comprobaciones exactas del modelo discreto (n = 200 capas):
#   - en κ = κ_y ninguna fibra ha fluido (la extrema está en h/2 − Δ/2), así
#     que M(κ_y) = E I_n κ_y = M_y (1 − 1/n²) exactamente;
#   - para κ > n κ_y fluyen todas (la más interior está en Δ/2) y
#     M = σ_y Σ |y_i| A_i = M_p exactamente (n par: |y| lineal en cada capa).
#   - descarga: con carga monótona la M–κ no distingue una sección que guarda
#     la historia de una que evalúa cada fibra desde el estado virgen (el
#     retorno de un solo paso es exacto). Una descarga Δκ = −1.5 κ_y sí: cada
#     fibra cambia su esfuerzo en E·y·Δκ, |Δσ| ≤ 1.5 σ_y < 2 σ_y, así que
#     ninguna fluye en sentido contrario y M' = M − E I_n · 1.5 κ_y exactamente.

KAPPA_SWEEP = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0, 100.0]
UNLOAD = 1.5           # descarga elástica, en unidades de κ_y


def _plastic_material():
    return Elastoplastic1D(E=E, sigma_y=SIGMA_Y, H=0.0)


def _EI_layers(depth: float, width: float, n: int) -> float:
    """E I de n capas en la dirección de flexión: E b h³/12 (1 − 1/n²)."""
    return E * width * depth**3 / 12.0 * (1.0 - 1.0 / n**2)


@pytest.mark.parametrize("axis", ["2d", "3d_z", "3d_y"])
def test_caso2_momento_curvatura_seccion(axis):
    """M(κ) de response_2d / response_3d frente a Neal: M_y, M_p, κ_y y M_p[1 − (κ_y/κ)²/3]."""
    sec, depth, width = axis_section(axis, _plastic_material(), N_LAYERS)
    kappa_y, M_y, M_p, N_p = rect_limits(depth, width)
    tol = tol_rect(N_LAYERS) * M_p
    states = sec.initial_states()
    for r in KAPPA_SWEEP:
        kappa = r * kappa_y
        N, M, M_other, _, states = section_response(sec, axis, 0.0, kappa, states)
        assert abs(N) <= TOL_EXACT * N_p, f"κ={r}κ_y: N = {N:.3e} ≠ 0 con ε0 = 0"
        assert abs(M_other) <= TOL_EXACT * M_p, f"κ={r}κ_y: momento en el otro eje {M_other:.3e}"
        M_an = m_rect(kappa, depth, width)
        assert abs(M - M_an) <= tol, (
            f"{axis} κ={r}κ_y: M = {M:.8e}, Neal = {M_an:.8e}, "
            f"|ΔM|/M_p = {abs(M - M_an) / M_p:.2e} > (5/3)/n² = {tol_rect(N_LAYERS):.2e}")
        if r == 1.0:
            np.testing.assert_allclose(M, M_y * (1.0 - 1.0 / N_LAYERS**2), rtol=1e-12,
                                       err_msg="M(κ_y) ≠ M_y (1 − 1/n²): núcleo aún elástico")
    kappa_full = 4.0 * N_LAYERS * kappa_y
    _, M, _, _, states = section_response(sec, axis, 0.0, kappa_full, states)
    assert abs(M - M_p) <= TOL_EXACT * M_p, (
        f"{axis}: con todas las fibras plastificadas M = {M:.10e} ≠ M_p = {M_p:.10e}")
    _, M_unl, _, _, _ = section_response(sec, axis, 0.0, kappa_full - UNLOAD * kappa_y, states)
    M_expected = M - _EI_layers(depth, width, N_LAYERS) * UNLOAD * kappa_y
    assert abs(M_unl - M_expected) <= TOL_EXACT * M_p, (
        f"{axis}: descarga Δκ = −{UNLOAD}κ_y: M = {M_unl:.10e}, descarga elástica = {M_expected:.10e} "
        "(¿la sección no conserva la historia de las fibras?)")


@pytest.mark.parametrize("axis", ["2d", "3d_z", "3d_y"])
def test_caso2_momento_curvatura_elemento(axis):
    """Curvatura uniforme impuesta en el elemento: F_θ2 = M de la sección = Neal.

    En los puntos de Gauss la curvatura es la misma, así que
    F_θ2 = M_s Σ w_g (L/2) N4''(x_g) = M_s exactamente (N4'' lineal).
    """
    sec, depth, width = axis_section(axis, _plastic_material(), N_LAYERS, GJ=GJ_RECT)
    if axis == "2d":
        elem = element_2d(Frame2DEuler, (0.0, 0.0), (L, 0.0), section=sec)
    else:
        elem = element_3d((0.0, 0.0, 0.0), (L, 0.0, 0.0), section=sec, ref_vector=REF_ALIGNED)
    kappa_y, _M_y, M_p, N_p = rect_limits(depth, width)
    tol = tol_rect(N_LAYERS) * M_p
    states = sec.initial_states()
    for r in KAPPA_SWEEP:
        kappa = r * kappa_y
        u, i2, i1 = uniform_curvature_dofs(axis, kappa)
        _, F = elem.compute_element_state(u)
        elem.commit_state()
        _, M_sec, _, _, states = section_response(sec, axis, 0.0, kappa, states)
        M_el = F[i2]
        assert abs(M_el - M_sec) <= 1e-12 * M_p, (
            f"{axis} κ={r}κ_y: F_θ2 = {M_el:.10e} ≠ M de la sección {M_sec:.10e}")
        assert abs(F[i1] + M_el) <= 1e-12 * M_p, f"{axis}: F_θ1 ≠ −F_θ2"
        others = np.delete(F, [i1, i2])
        assert np.abs(others).max() <= TOL_EXACT * N_p, (
            f"{axis} κ={r}κ_y: fuerzas nodales espurias {others}")
        M_an = m_rect(kappa, depth, width)
        assert abs(M_el - M_an) <= tol, (
            f"{axis} κ={r}κ_y: M elemento = {M_el:.8e}, Neal = {M_an:.8e}")
    # Descarga elástica (ver el encabezado del caso): prueba que el elemento
    # conserva y confirma el estado de las fibras de cada punto de Gauss.
    u, i2, _ = uniform_curvature_dofs(axis, (KAPPA_SWEEP[-1] - UNLOAD) * kappa_y)
    _, F = elem.compute_element_state(u)
    M_expected = M_el - _EI_layers(depth, width, N_LAYERS) * UNLOAD * kappa_y
    assert abs(F[i2] - M_expected) <= TOL_EXACT * M_p, (
        f"{axis}: descarga del elemento: F_θ2 = {F[i2]:.10e}, descarga elástica = {M_expected:.10e} "
        "(¿el elemento no confirma el estado de las fibras en commit_state?)")


# =============================================================================
# Caso 3 — interacción axial–flector
# =============================================================================
# N impuesta primero (κ = 0, respuesta elástica) y curvatura creciente con N
# constante: en cada incremento se resuelve N(ε0; κ) = N con Newton
# salvaguardado por bisección (N es monótona en ε0 porque E_t ≥ 0) y se
# aceptan los estados como hace el elemento. dM/dκ con N constante es el
# complemento de Schur de k_s, que es semidefinida positiva: M no decrece.

KAPPA_AXIAL = np.concatenate([[0.0], np.geomspace(0.05, 100.0, 60)])


def _solve_axial(sec, axis, kappa, N_target, states, eps0_guess, N_p):
    """ε0 tal que N(ε0; κ) = N_target partiendo de ``states`` (sin modificarlos)."""
    def evaluate(eps0):
        N, M, _, k_nn, new = section_response(sec, axis, eps0, kappa, states)
        return N - N_target, k_nn, N, M, new

    res, k_nn, N, M, new = evaluate(eps0_guess)
    lo = hi = eps0_guess
    res_lo = res_hi = res
    step = EPS_Y
    while res_lo > 0.0:
        lo -= step
        step *= 2.0
        res_lo = evaluate(lo)[0]
    step = EPS_Y
    while res_hi < 0.0:
        hi += step
        step *= 2.0
        res_hi = evaluate(hi)[0]
    eps0 = eps0_guess
    for _ in range(200):
        if abs(res) <= 1e-11 * N_p:
            return eps0, N, M, new
        if res > 0.0:
            hi = eps0
        else:
            lo = eps0
        cand = eps0 - res / k_nn if k_nn > 0.0 else 0.5 * (lo + hi)
        if not lo < cand < hi:
            cand = 0.5 * (lo + hi)
        eps0 = cand
        res, k_nn, N, M, new = evaluate(eps0)
    raise AssertionError(f"sin convergencia en ε0 para κ = {kappa:.3e} (residuo {res:.3e})")


@pytest.mark.parametrize("axis, n_ratio", [("2d", 0.3), ("2d", 0.6), ("2d", -0.5),
                                           ("3d_y", 0.4), ("3d_z", -0.3)])
def test_caso3_interaccion_axial_flector(axis, n_ratio):
    """Con N = n·N_p fija, M(κ) = M_p[1 − n² − (κ_y/κ)²/3] → M_p[1 − (N/N_p)²]."""
    sec, depth, width = axis_section(axis, _plastic_material(), N_LAYERS)
    kappa_y, _M_y, M_p, N_p = rect_limits(depth, width)
    N_target = n_ratio * N_p
    kappa_both = kappa_y / (1.0 - abs(n_ratio))
    tol = tol_rect_axial(N_LAYERS) * M_p
    states = sec.initial_states()
    eps0, M_prev, n_checked = 0.0, -np.inf, 0
    for r in KAPPA_AXIAL:
        kappa = r * kappa_y
        eps0, N, M, states = _solve_axial(sec, axis, kappa, N_target, states, eps0, N_p)
        assert abs(N - N_target) <= 1e-10 * N_p
        assert M >= M_prev - TOL_EXACT * M_p, f"κ={r:.3g}κ_y: M decrece con N constante"
        M_prev = M
        if kappa >= 1.001 * kappa_both:
            M_an = m_rect_axial(kappa, n_ratio, depth, width)
            assert abs(M - M_an) <= tol, (
                f"{axis} N/N_p={n_ratio} κ={r:.3g}κ_y: M = {M:.8e}, cerrada = {M_an:.8e}, "
                f"|ΔM|/M_p = {abs(M - M_an) / M_p:.2e} > (8/3)/n² = {tol_rect_axial(N_LAYERS):.2e}")
            n_checked += 1
    assert n_checked >= 20
    M_lim = M_p * (1.0 - n_ratio**2)
    gap = M_p * ((kappa_y / kappa) ** 2 / 3.0 + tol_rect_axial(N_LAYERS))
    assert abs(M - M_lim) <= gap, (
        f"{axis}: en κ = {KAPPA_AXIAL[-1]:.0f}κ_y, M = {M:.8e} lejos de M_p[1 − (N/N_p)²] = {M_lim:.8e}")


# =============================================================================
# Casos 4 y 5 — cargas de colapso: por qué el elemento en desplazamientos
# sobreestima la carga
# =============================================================================
# Voladizo rígido-plástico: rótula en el empotramiento (x = 0), P_c = M_p/L.
# El elemento sólo evalúa la sección en los puntos de Gauss, que nunca están
# en el extremo (el primero, en x_1 = 0.113 h_e con n_ip = 3), y su curvatura
# es lineal: para que κ_1 crezca sin límite con κ_2 (x = h_e/2) acotada haría
# falta curvatura de signo contrario en x_3, así que en el límite los puntos 1
# y 2 alcanzan juntos M_p. La rótula discreta queda repartida hacia el interior
# del elemento raíz y el brazo efectivo es menor que L: el modelo sobreestima
# la carga de colapso.
#
# Análisis de límite del modelo discreto. Las fuerzas nodales de giro del
# elemento raíz son F_θ1 = Σ c_g M_g y F_θ2 = Σ d_g M_g, con
# c_g = (w_g/2)(6ξ_g − 4) y d_g = (w_g/2)(6ξ_g − 2) (N2'' y N4'' de la spec en
# ξ = x/h_e). Las fuerzas de cada elemento están autoequilibradas, de modo que
# la estática del modelo discreto es exacta: F_θ1 = −P L (reacción) y
# F_θ2 = P (L − h_e) (lo que transmite el resto del voladizo). Eliminando M_3,
#
#     P [L + (C/A)(L − h_e)] = (A − C²/A) M_1 + B (1 + C/A) M_2,
#
# con A = −c_1 = 0.9233, B = −c_2 = 4/9, C = c_3 = 0.3677: exacta para
# cualesquiera M_1, M_2. Con |M_g| ≤ M_p (fibras con |σ| ≤ σ_y y número par de
# capas, así que Σ σ_y |y_i| A_i = M_p) da la cota superior del modelo discreto
#
#     P ≤ P_∞(n) = P_c · 1.398276/(1.398276 − 0.398276/n),
#
# 7.7 %, 3.7 %, 1.8 % y 0.90 % por encima de P_c con 4, 8, 16 y 32 elementos:
# unas 2.5 veces el estimado ingenuo ξ_1/n (sólo el punto de Gauss desplazado)
# y de primer orden en h_e. Con desplazamiento finito, M_g = M(κ_g) < M_p y la
# carga queda algo por debajo de P_∞(n). κ_1 y κ_2 del elemento raíz se leen de
# los DOF nodales (curvatura de Hermite, spec §5), así que la misma identidad
# predice P con la M–κ analítica: la diferencia es el error de discretización
# de las fibras, (5/3)/n_f² de M_p.
#
# Viga biempotrada con carga central: M = 0 en L/4 y 3L/4 por la simetría
# respecto del centro y la antisimetría de cada semivano (giro nulo en sus dos
# extremos) respecto de su punto medio; cada cuarto es un voladizo de longitud
# L/4 con n/4 elementos y carga P/2. Mismas fórmulas con L → L/4, P → P/2;
# P_c = 8 M_p/L (rótulas en los empotramientos y en el centro).
#
# Desplazamiento impuesto final: 5 δ_y (δ_y = κ_y L²/3 en el voladizo,
# κ_y L²/24 en la biempotrada, desplazamiento elástico en la primera
# fluencia). Suficiente para que la curvatura en los puntos 1–2 del elemento
# raíz sea varias veces κ_y incluso con la malla gruesa (el déficit
# (κ_y/κ)²/3 queda muy por debajo de la separación entre mallas) y, con la
# malla fina, por debajo de n_f κ_y: las fibras centrales siguen elásticas y la
# tangente no es singular. La historia es monótona en cada fibra, así que el
# resultado no depende del número de pasos.

N_FIB_COLLAPSE = 200
MESHES = (4, 8, 16, 32)
DELTA_MULT = 5.0
_XG, _WG = np.polynomial.legendre.leggauss(3)
_XI = (_XG + 1.0) / 2.0
C_G = _WG / 2.0 * (6.0 * _XI - 4.0)     # peso de M_g en F_θ1 (w_g · L/2 · N2'')
D_G = _WG / 2.0 * (6.0 * _XI - 2.0)     # peso de M_g en F_θ2 (w_g · L/2 · N4'')


def root_element_load(M1: float, M2: float, arm: float, h: float):
    """(P, M_3) de la estática exacta del elemento raíz dados M_1 y M_2.

    Σ c_g M_g = −P·arm y Σ d_g M_g = P (arm − h), resueltas en (P, M_3).
    """
    A = np.array([[arm, C_G[2]], [-(arm - h), D_G[2]]])
    rhs = np.array([-(C_G[0] * M1 + C_G[1] * M2), -(D_G[0] * M1 + D_G[1] * M2)])
    P, M3 = np.linalg.solve(A, rhs)
    return float(P), float(M3)


def root_curvatures(v2: float, theta2: float, h: float) -> np.ndarray:
    """κ en los tres puntos de Gauss del elemento raíz (v1 = θ1 = 0), spec §5."""
    return (6.0 - 12.0 * _XI) / h**2 * v2 + (6.0 * _XI - 2.0) / h * theta2


def run_collapse(support: str, n_el: int):
    """Marco de n_el Frame2DEuler con fibras, control por desplazamiento.

    Devuelve (reacción en el DOF controlado por paso, U final, nodos).
    """
    sec = rect_section(_plastic_material(), N_FIB_COLLAPSE, 1)
    dom = Domain()
    nodes = [dom.add_node(i + 1, [L * i / n_el, 0.0]) for i in range(n_el + 1)]
    for i in range(n_el):
        dom.add_element(Frame2DEuler(i + 1, [nodes[i], nodes[i + 1]], section=sec, n_ip=3))
    ends = [nodes[0]] if support == "cantilever" else [nodes[0], nodes[-1]]
    for node in ends:
        for dof in ("ux", "uy", "rz"):
            node.fix_dof(dof, 0.0)
    kappa_y = rect_limits(H, B)[0]
    if support == "cantilever":
        ctrl, delta_y = nodes[-1], kappa_y * L**2 / 3.0
    else:
        ctrl, delta_y = nodes[n_el // 2], kappa_y * L**2 / 24.0
    ctrl.fix_dof("uy", DELTA_MULT * delta_y)
    dom.generate_equation_numbers(verbose=False)
    asm = Assembler(dom)
    loads: list[float] = []

    def record(_step, U, _lam):
        _, F_int = asm.assemble_non_linear_system(U)
        loads.append(float(F_int[ctrl.dofs["uy"]]))

    solver = NonlinearSolver(asm, convergence=ConvergenceCriterion(rtol_force=1e-9, rtol_disp=1e-9),
                             num_steps=10, max_iter=40)
    U = solver.solve(np.zeros(dom.total_dofs), step_callback=record)
    return np.array(loads), U, nodes


def _check_collapse(support: str) -> None:
    kappa_y, _M_y, M_p, _N_p = rect_limits(H, B)
    if support == "cantilever":
        P_c, arm_factor, scale = M_p / L, 1.0, 1.0
    else:
        P_c, arm_factor, scale = 8.0 * M_p / L, 0.25, 2.0
    tol_pred = tol_rect(N_FIB_COLLAPSE) + 1e-6          # fibras + residuo de Newton
    P_max, P_inf = {}, {}
    for n in MESHES:
        loads, U, nodes = run_collapse(support, n)
        h, arm = L / n, arm_factor * L
        assert np.all(np.diff(loads) >= -TOL_EXACT * P_c), f"n={n}: la carga decrece en la historia"
        P_max[n] = loads.max()

        P_q_inf, M3_inf = root_element_load(M_p, M_p, arm, h)
        assert abs(M3_inf) <= M_p, "la cota de límite exige |M_3| ≤ M_p"
        P_inf[n] = scale * P_q_inf
        assert P_max[n] <= P_inf[n] * (1.0 + 1e-6), (
            f"{support} n={n}: P_max/P_c = {P_max[n] / P_c:.6f} supera la cota de límite "
            f"del modelo discreto P_∞/P_c = {P_inf[n] / P_c:.6f}")

        kap = root_curvatures(U[nodes[1].dofs["uy"]], U[nodes[1].dofs["rz"]], h)
        assert kap[0] > 0.0 and kap[1] > 0.0, f"n={n}: curvatura de signo inesperado {kap}"
        assert kap[0] < N_FIB_COLLAPSE * kappa_y, (
            f"n={n}: κ_1 = {kap[0] / kappa_y:.1f} κ_y ≥ n_f κ_y: sin núcleo elástico el test no es válido")
        P_pred = scale * root_element_load(m_rect(kap[0], H, B), m_rect(kap[1], H, B), arm, h)[0]
        assert abs(P_max[n] - P_pred) <= tol_pred * P_inf[n], (
            f"{support} n={n}: P = {P_max[n]:.8e} y la estática del elemento raíz con la M–κ "
            f"analítica en κ_1 = {kap[0] / kappa_y:.2f} κ_y, κ_2 = {kap[1] / kappa_y:.2f} κ_y "
            f"predice {P_pred:.8e} (|Δ|/P_∞ = {abs(P_max[n] - P_pred) / P_inf[n]:.2e})")

    ratios = [P_max[n] / P_c for n in MESHES]
    assert all(a > b for a, b in zip(ratios, ratios[1:])), (
        f"{support}: la carga máxima no decrece con la malla: P_max/P_c = {ratios}")
    assert ratios[-1] > 1.0, f"{support}: la malla fina queda por debajo de P_c ({ratios[-1]:.6f})"
    e_final = ratios[-1] - 1.0
    e_bound = P_inf[MESHES[-1]] / P_c - 1.0
    assert e_final <= e_bound * (1.0 + 1e-6), (
        f"{support}: error final {e_final:.3%} > cota razonada {e_bound:.3%}")


def test_caso4_colapso_voladizo_carga_en_el_extremo():
    """Voladizo, 4/8/16/32 elementos: P_max decreciente → P_c = M_p/L (error final ≤ 0.90 %)."""
    _check_collapse("cantilever")


def test_caso5_colapso_viga_biempotrada_carga_central():
    """Biempotrada, 4/8/16/32 elementos: P_max decreciente → P_c = 8 M_p/L (error final ≤ 3.7 %)."""
    _check_collapse("fixed")


# =============================================================================
# Caso 6 — sección circular maciza en flexión biaxial a 45° (Frame3D)
# =============================================================================
# Con n_circ múltiplo de 8 y el primer sector empezando en el ángulo 0, un
# giro de 45° (n_circ/8 sectores) lleva la sección discreta sobre sí misma:
# la curva M–κ a 45° es **idéntica** a la uniaxial, no sólo parecida. La
# reflexión respecto de la bisectriz (y ↔ z) cambia el signo de ε, así que
# M_y = M_z exactamente. Además la bisectriz y el eje z (y = 0) son fronteras
# de sector, así que |y| y |z − y| son lineales en cada sector y la regla del
# centroide da M_p exacto con todas las fibras plastificadas.

D_CIRC = 0.4
R_CIRC = D_CIRC / 2.0
N_RAD, N_CIRC = 20, 64
GJ_CIRC = G * math.pi * R_CIRC**4 / 2.0
M_Y_CIRC = math.pi * SIGMA_Y * D_CIRC**3 / 32.0      # 1.5708e6 N·m
M_P_CIRC = SIGMA_Y * D_CIRC**3 / 6.0                 # 2.6667e6 N·m
KAPPA_Y_CIRC = 2.0 * SIGMA_Y / (E * D_CIRC)          # 6.25e-3 1/m
N_P_CIRC = SIGMA_Y * math.pi * R_CIRC**2
KAPPA_SWEEP_CIRC = [0.25, 0.5, 1.0, 1.25, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0]


def circular_section(material) -> FiberSection:
    patch = CircularPatch(material, r_outer=R_CIRC, r_inner=0.0, center=(0.0, 0.0),
                          n_rad=N_RAD, n_circ=N_CIRC)
    return FiberSection(patches=[patch], GJ=GJ_CIRC)


def circle_layout():
    """Radios del centroide y áreas por anillo (división radial uniforme).

    Centroide de un sector anular [r1, r2] × [φ1, φ1 + Δφ]:
    r_c = (2/3)(r2³ − r1³)/(r2² − r1²) · sen(Δφ/2)/(Δφ/2), en φ1 + Δφ/2.
    """
    r = np.linspace(0.0, R_CIRC, N_RAD + 1)
    r1, r2 = r[:-1], r[1:]
    dphi = 2.0 * math.pi / N_CIRC
    area = 0.5 * dphi * (r2**2 - r1**2)
    rc = (2.0 / 3.0) * (r2**3 - r1**3) / (r2**2 - r1**2) * math.sin(dphi / 2.0) / (dphi / 2.0)
    return rc, area, dphi


def _discrete_inertia() -> float:
    """Σ A y² de la disposición de la API: Σ_j cos²((j + ½)Δφ) = n_circ/2."""
    rc, area, _ = circle_layout()
    return float(np.sum(area * rc**2) * N_CIRC / 2.0)


def test_caso6_disposicion_de_la_seccion_circular():
    """Fibras en el centroide de sectores anulares, área exacta, primer sector en φ = 0."""
    sec = circular_section(Elastic1D(E=E))
    rc, area_ring, dphi = circle_layout()
    assert sec.n_fibers == N_RAD * N_CIRC
    np.testing.assert_allclose(sec.area.sum(), math.pi * R_CIRC**2, rtol=1e-12)
    r = np.hypot(sec.y, sec.z)
    ring = np.argmin(np.abs(r[:, None] - rc[None, :]), axis=1)
    np.testing.assert_allclose(r, rc[ring], rtol=1e-12, err_msg="radio ≠ radio del centroide del sector")
    assert np.all(np.bincount(ring, minlength=N_RAD) == N_CIRC)
    np.testing.assert_allclose(sec.area, area_ring[ring], rtol=1e-12, err_msg="área ≠ área del sector")
    phi = np.mod(np.arctan2(sec.z, sec.y), 2.0 * math.pi)
    j = np.floor(phi / dphi).astype(int)
    np.testing.assert_allclose(phi - (j + 0.5) * dphi, 0.0, atol=1e-9,
                               err_msg="los sectores no empiezan en φ = 0 (eje +y hacia +z)")
    assert np.all(np.bincount(j, minlength=N_CIRC) == N_RAD)

    I_d = _discrete_inertia()
    np.testing.assert_allclose(np.sum(sec.area * sec.y**2), I_d, rtol=1e-12)
    np.testing.assert_allclose(np.sum(sec.area * sec.z**2), I_d, rtol=1e-12)
    # Déficit de la regla del centroide: Σ A Var_s(y) ≈ I (1/6)[(Δr/R)² + Δφ²/2];
    # el resto es O(h²) relativo (~1 %): margen del 10 %.
    deficit = 1.0 - I_d / (math.pi * R_CIRC**4 / 4.0)
    estimate = ((1.0 / N_RAD) ** 2 + 0.5 * dphi**2) / 6.0
    assert deficit > 0.0 and abs(deficit / estimate - 1.0) <= 0.10, (deficit, estimate)


def test_caso6_momento_curvatura_circular_uniaxial_y_45():
    """Misma M–κ a 45° que uniaxial; M_y = πσ_y d³/32, M_p = σ_y d³/6 (response_3d)."""
    sec = circular_section(_plastic_material())
    tol = tol_circle(N_RAD, N_CIRC) * M_P_CIRC
    I_d = _discrete_inertia()
    st_u, st_45 = sec.initial_states(), sec.initial_states()
    for r in KAPPA_SWEEP_CIRC:
        kappa = r * KAPPA_Y_CIRC
        k45 = kappa / math.sqrt(2.0)
        (N_u, My_u, Mz_u), _, st_u = sec.response_3d(np.array([0.0, 0.0, kappa]), st_u)
        (N_45, My_45, Mz_45), _, st_45 = sec.response_3d(np.array([0.0, k45, k45]), st_45)
        assert abs(N_u) <= TOL_EXACT * N_P_CIRC and abs(N_45) <= TOL_EXACT * N_P_CIRC
        assert abs(My_u) <= TOL_EXACT * M_P_CIRC, f"κ={r}κ_y: M_y ≠ 0 en flexión uniaxial"
        assert abs(My_45 - Mz_45) <= TOL_EXACT * M_P_CIRC, (
            f"κ={r}κ_y: a 45° M_y = {My_45:.10e} ≠ M_z = {Mz_45:.10e}")
        M_45 = math.hypot(My_45, Mz_45)
        assert abs(M_45 - Mz_u) <= TOL_EXACT * M_P_CIRC, (
            f"κ={r}κ_y: |M| a 45° = {M_45:.10e} ≠ uniaxial {Mz_u:.10e}")
        M_an = m_circle(kappa, R_CIRC)
        assert abs(Mz_u - M_an) <= tol, (
            f"κ={r}κ_y: M = {Mz_u:.8e}, analítica = {M_an:.8e}, "
            f"|ΔM|/M_p = {abs(Mz_u - M_an) / M_P_CIRC:.2e} > {tol_circle(N_RAD, N_CIRC):.2e}")
        if r == 1.0:
            # Ninguna fibra ha fluido (todas dentro de r < R): M = E I_d κ_y.
            np.testing.assert_allclose(Mz_u, E * I_d * KAPPA_Y_CIRC, rtol=1e-12)
            assert abs(Mz_u / M_Y_CIRC - 1.0) <= 1.1 * ((1.0 / N_RAD) ** 2
                                                        + 0.5 * (2.0 * math.pi / N_CIRC) ** 2) / 6.0

    d_min = min(np.abs(sec.y).min(), (np.abs(sec.z - sec.y) / math.sqrt(2.0)).min())
    assert d_min > 0.0, "hay fibras sobre y = 0 o sobre y = z: los sectores no empiezan en φ = 0"
    kappa_full = 2.0 * EPS_Y / d_min
    k45 = kappa_full / math.sqrt(2.0)
    s_u, _, _ = sec.response_3d(np.array([0.0, 0.0, kappa_full]), st_u)
    s_45, _, _ = sec.response_3d(np.array([0.0, k45, k45]), st_45)
    assert abs(s_u[2] - M_P_CIRC) <= TOL_EXACT * M_P_CIRC, (
        f"plastificación completa uniaxial: M = {s_u[2]:.10e} ≠ M_p = {M_P_CIRC:.10e}")
    assert abs(math.hypot(s_45[1], s_45[2]) - M_P_CIRC) <= TOL_EXACT * M_P_CIRC, (
        "plastificación completa a 45°: |M| ≠ M_p")


def test_caso6_frame3d_flexion_biaxial_45():
    """Frame3D con sección circular y curvatura uniforme a 45°: M_y = M_z y |M| = M uniaxial."""
    sec = circular_section(_plastic_material())
    elem = element_3d((0.0, 0.0, 0.0), (L, 0.0, 0.0), section=sec, ref_vector=REF_ALIGNED)
    tol = tol_circle(N_RAD, N_CIRC) * M_P_CIRC
    st_u = sec.initial_states()
    for r in [0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0]:
        kappa = r * KAPPA_Y_CIRC
        half = kappa / math.sqrt(2.0) * L / 2.0
        u = np.zeros(12)
        u[4], u[10] = -half, half          # θ_y: κ_y = κ/√2
        u[5], u[11] = -half, half          # θ_z: κ_z = κ/√2
        _, F = elem.compute_element_state(u)
        elem.commit_state()
        s_u, _, st_u = sec.response_3d(np.array([0.0, 0.0, kappa]), st_u)
        My, Mz = F[10], F[11]
        assert abs(My - Mz) <= TOL_EXACT * M_P_CIRC, f"κ={r}κ_y: M_y = {My:.10e} ≠ M_z = {Mz:.10e}"
        assert abs(F[4] + My) <= TOL_EXACT * M_P_CIRC and abs(F[5] + Mz) <= TOL_EXACT * M_P_CIRC
        others = F[[0, 1, 2, 3, 6, 7, 8, 9]]
        assert np.abs(others).max() <= TOL_EXACT * N_P_CIRC, f"κ={r}κ_y: fuerzas espurias {others}"
        M = math.hypot(My, Mz)
        assert abs(M - s_u[2]) <= TOL_EXACT * M_P_CIRC, (
            f"κ={r}κ_y: |M| elemento a 45° = {M:.10e} ≠ M uniaxial de la sección {s_u[2]:.10e}")
        assert abs(M - m_circle(kappa, R_CIRC)) <= tol


# =============================================================================
# Caso 7 — convenciones de signo (Reglas §5)
# =============================================================================

def test_caso7_signos_deformaciones_de_fibra():
    """κ_z > 0 tracciona y < 0 y da M_z > 0; κ_y > 0 tracciona z > 0 y da M_y > 0."""
    sec = rect_section(Elastic1D(E=E), 10, 6)     # sin fibras en y = 0 ni en z = 0
    kappa = 1.0e-3
    atol = 1e-15 * kappa * H
    states = sec.initial_states()

    eps = sec.strains_2d(np.array([0.0, kappa]))
    np.testing.assert_allclose(eps, -sec.y * kappa, rtol=0.0, atol=atol)
    assert np.all(eps[sec.y < 0.0] > 0.0) and np.all(eps[sec.y > 0.0] < 0.0)
    s, _, _ = sec.response_2d(np.array([0.0, kappa]), states)
    assert s[1] > 0.0, "κ_z > 0 debe dar M_z > 0 (tracción en y < 0)"

    eps = sec.strains_3d(np.array([0.0, kappa, 0.0]))
    np.testing.assert_allclose(eps, sec.z * kappa, rtol=0.0, atol=atol)
    assert np.all(eps[sec.z > 0.0] > 0.0) and np.all(eps[sec.z < 0.0] < 0.0)
    s, _, _ = sec.response_3d(np.array([0.0, kappa, 0.0]), states)
    assert s[1] > 0.0 and abs(s[2]) <= 1e-12 * abs(s[1]), "κ_y > 0 debe dar M_y > 0 y M_z = 0"

    eps = sec.strains_3d(np.array([0.0, 0.0, kappa]))
    np.testing.assert_allclose(eps, -sec.y * kappa, rtol=0.0, atol=atol)
    assert np.all(eps[sec.y < 0.0] > 0.0) and np.all(eps[sec.y > 0.0] < 0.0)
    s, _, _ = sec.response_3d(np.array([0.0, 0.0, kappa]), states)
    assert s[2] > 0.0 and abs(s[1]) <= 1e-12 * abs(s[2]), "κ_z > 0 debe dar M_z > 0 y M_y = 0"

    e = np.array([3.0e-4, -7.0e-4, 1.1e-3])
    np.testing.assert_allclose(sec.strains_3d(e), e[0] - sec.y * e[2] + sec.z * e[1],
                               rtol=0.0, atol=1e-15)
    np.testing.assert_allclose(sec.strains_2d(e[[0, 2]]), e[0] - sec.y * e[2], rtol=0.0, atol=1e-15)


def test_caso7_fibra_excentrica_fuerzas_internas_y_tangente():
    """Una fibra en (y, z) = (0.1, 0.05) traccionada: N > 0, M_y > 0, M_z < 0 (spec §4)."""
    y0, z0, a0 = 0.1, 0.05, 1.0e-4
    sec = FiberSection(fibers=[Fiber(y0, z0, a0, Elastic1D(E=E))])
    eps0 = 5.0e-4
    force = E * eps0 * a0

    s3, k3, _ = sec.response_3d(np.array([eps0, 0.0, 0.0]), sec.initial_states())
    np.testing.assert_allclose(s3, force * np.array([1.0, z0, -y0]), rtol=1e-14)
    k3_ref = E * a0 * np.array([[1.0, z0, -y0], [z0, z0**2, -y0 * z0], [-y0, -y0 * z0, y0**2]])
    np.testing.assert_allclose(k3, k3_ref, rtol=1e-14)
    assert s3[1] > 0.0, "tracción en z > 0 ⇒ M_y > 0"
    assert s3[2] < 0.0, "tracción en y > 0 ⇒ M_z < 0"

    s2, k2, _ = sec.response_2d(np.array([eps0, 0.0]), sec.initial_states())
    np.testing.assert_allclose(s2, force * np.array([1.0, -y0]), rtol=1e-14)
    np.testing.assert_allclose(k2, E * a0 * np.array([[1.0, -y0], [-y0, y0**2]]), rtol=1e-14)


# =============================================================================
# Caso 8 — corrotacional elástico con giros grandes
# =============================================================================
# Frame2DEulerCorot: la sección trabaja en el sistema básico
# v = [Δl, θ̄1, θ̄2]; con fibras elásticas las fuerzas básicas y su tangente
# coinciden con las de A e I integrados para n_ip ≥ 2 (integrando cuadrático),
# y la cinemática corrotacional es la misma: K y F_int iguales a 1e-10.

@pytest.mark.parametrize("alpha", [math.pi / 3.0, 5.0 * math.pi / 6.0, -2.0 * math.pi / 3.0])
@pytest.mark.parametrize("n_ip", [2, 3])
def test_caso8_corotacional_elemento_giro_grande(n_ip, alpha):
    """Giro de sólido rígido α grande + alargamiento y giros deformacionales."""
    sec, E_ref, EA, EIz, _EIy, _mats, _GJ = elastic_section("rect")
    ref = element_2d(Frame2DEulerCorot, (0.0, 0.0), (L, 0.0),
                     material=Elastic1D(E=E_ref), A=EA / E_ref, I=EIz / E_ref)
    fib = element_2d(Frame2DEulerCorot, (0.0, 0.0), (L, 0.0), section=sec, n_ip=n_ip)
    stretch = 1.0e-4
    u = np.array([0.0, 0.0, alpha + 0.05,
                  (1.0 + stretch) * L * math.cos(alpha) - L,
                  (1.0 + stretch) * L * math.sin(alpha),
                  alpha - 0.08])
    K_r, F_r = ref.compute_element_state(u)
    K_f, F_f = fib.compute_element_state(u)
    assert_same_response(K_f, F_f, K_r, F_r, 1.0e-10, f"corot n_ip={n_ip} α={math.degrees(alpha):.0f}°")


def test_caso8_corotacional_arco_de_Bathe_Bolourchi():
    """Voladizo de 10 elementos con momento en el extremo hasta θ = π/2: misma solución.

    Con K y F_int idénticos a redondeo, las dos corridas de Newton recorren los
    mismos iterados; convergidas con rtol_disp = 1e-11, la diferencia de las
    soluciones está acotada por la última corrección, ≤ 1e-11 ‖U‖. Se exige
    1e-10 ‖U‖ (la tolerancia del contrato).
    """
    sec, E_ref, EA, EIz, _EIy, _mats, _GJ = elastic_section("rect")
    n_el = 10
    M_tip = (math.pi / 2.0) * EIz / L

    def solve(**kwargs):
        dom = Domain()
        nodes = [dom.add_node(k + 1, [k * L / n_el, 0.0]) for k in range(n_el + 1)]
        for k in range(n_el):
            dom.add_element(Frame2DEulerCorot(k + 1, [nodes[k], nodes[k + 1]], **kwargs))
        for dof in ("ux", "uy", "rz"):
            nodes[0].fix_dof(dof, 0.0)
        dom.generate_equation_numbers(verbose=False)
        F_ext = np.zeros(dom.total_dofs)
        F_ext[nodes[-1].dofs["rz"]] = M_tip
        solver = NonlinearSolver(Assembler(dom),
                                 convergence=ConvergenceCriterion(rtol_force=1e-11, rtol_disp=1e-11),
                                 num_steps=8, max_iter=30)
        return solver.solve(F_ext), nodes[-1]

    U_ref, tip_ref = solve(material=Elastic1D(E=E_ref), A=EA / E_ref, I=EIz / E_ref)
    U_fib, _ = solve(section=sec)
    assert abs(U_ref[tip_ref.dofs["rz"]] - math.pi / 2.0) <= 1e-4, "la referencia no llegó al cuarto de círculo"
    diff = np.linalg.norm(U_fib - U_ref)
    assert diff <= 1e-10 * np.linalg.norm(U_ref), (
        f"‖U_fibras − U_A,I‖/‖U‖ = {diff / np.linalg.norm(U_ref):.2e} > 1e-10")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
