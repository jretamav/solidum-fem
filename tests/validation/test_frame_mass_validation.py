"""Validación ciega de la masa de los marcos Euler-Bernoulli (hallazgo clas-4).

Física aprobada por el usuario (2026-09-29) para ``Frame2DEuler``,
``Frame2DEulerCorot`` y la flexión de ``Frame3D``, por los dos caminos (con
material y con sección, en el que ρA = Σρ_i A_i de ``mass_properties``):

- **Consistente** clásica **sin inercia rotatoria** (la hipótesis de
  Euler-Bernoulli la desprecia): axial ρAL/6·[[2, 1], [1, 2]] y flexión
  hermítica en (v1, θ1, v2, θ2)::

      ρAL/420 · [[ 156,   22L,   54,  −13L ],
                 [ 22L,   4L²,  13L,  −3L² ],
                 [  54,   13L,  156,  −22L ],
                 [−13L,  −3L², −22L,   4L² ]]

  En el plano xz de ``Frame3D``, (w1, θy1, w2, θy2), la cinemática w' = −θy
  (regla de la mano derecha) equivale a sustituir θz → −θy: la matriz es
  S·M_xy·S con S = diag(1, −1, 1, −1), es decir, los acoplamientos
  traslación-giro cambian de signo y los bloques w-w y θ-θ no.
- **Concentrada HRZ** (Hinton, Rock y Zienkiewicz 1976; Cook et al.): diagonal
  de la consistente escalada para conservar la masa traslacional. Por nodo:
  axial ρAL/2, transversal ρAL/2 (156/312) y giro de flexión
  m_θ = ρAL·4L²/312 = ρAL³/78.
- **Torsión de Frame3D** (sin cambio): consistente ρJ_pL/6·[[2, 1], [1, 2]] y
  concentrada ρJ_pL/2 por nodo, con J_p = I_y + I_z (con sección,
  ρJ_p = ρI_y + ρI_z = Σρ_i A_i (y_i² + z_i²)). El J de Saint-Venant sólo
  entra en la rigidez.
- Consecuencia: la matriz de masa de flexión **no depende de I**; en
  ``Frame3D`` depende de (I_y, I_z) sólo a través de I_y + I_z, en la
  torsión.

Test ciego
----------
Escrito desde la física anterior y las specs públicas (``Frame2DEuler`` §7 para
la transformación 2D, ``Frame3D`` §6-§7 para los ejes locales), **sin leer**
``solidum/elements/frame/`` ni ``solidum/elements/frame3d.py``. Todas las
cifras esperadas salen de fórmulas cerradas o de identidades exactas del
modelo discreto construidas aquí mismo, nunca de una ejecución de la
implementación que se valida.

Datos (distintos de 1 a propósito: los unitarios esconden errores
dimensionales): E = 2.1e11 Pa, ρ = 7850 kg/m³, A = 0.02 m², I = 3e-4 m⁴,
L = 2.5 m (elemento aislado y longitud total de las vigas), ν = 0.3,
G = E/(2(1 + ν)) = 8.0769e10 Pa. Frame3D: I_y = 1e-4, I_z = 3e-4, J = 2.2e-4
(J ≠ I_y + I_z para distinguir J de J_p). Sección por fibras: rectángulo de
peralte h = √(12 I/A) = 0.42426 m y ancho b = A/h = 0.047140 m (b h = A),
20 × 6 celdas: I_z,n = b h³/12 (1 − 1/20²) = 2.9925e-4 m⁴,
I_y,n = h b³/12 (1 − 1/6²) = 3.6008e-6 m⁴, GJ = G·1.4e-5.

Valores analíticos de referencia (medio continuo)
-------------------------------------------------
- Viga biapoyada Euler-Bernoulli: ω_k = (kπ/L)²·√(EI/(ρA)):
  ω_1 = 1000.3228 rad/s, ω_2 = 4001.2913 rad/s (camino con sección, con
  I_z,n: 999.0716 y 3996.2865 rad/s).
- Voladizo: ω_1 = (β_1 L)²·√(EI/(ρA L⁴)), β_1 L = 1.8751040687 (raíz de
  cos x cosh x = −1): 356.3618 rad/s (con sección 355.9161 rad/s).
- Frame3D biapoyado: plano xz (I_y) 577.5366 y 2310.1466 rad/s; plano xy (I_z)
  1000.3228 y 4001.2913 rad/s.
- Torsión empotrada-libre: ω_1 = (π/2L)·√(GJ/(ρJ_p)) = 1494.6839 rad/s (con
  sección 433.3289 rad/s).

Identidades exactas del modelo discreto
---------------------------------------
(1) **Símbolo de la viga biapoyada.** En una malla uniforme de n elementos de
    longitud h, el modo v_j = V sin(jφ), θ_j = Θ cos(jφ), φ = kπ/n, cumple
    exactamente las ecuaciones de todos los nodos (en los apoyos, por la
    extensión especular v_{−j} = −v_j, θ_{−j} = θ_j, la ecuación interior es el
    doble de la de contorno). Con X = Θh y t = sen²(φ/2), el problema queda en
    el generalizado 2×2 K̂ − λ̃ M̂ con λ = λ̃ (EI/h³)/m::

        K̂ = [[48t, −12 sen φ], [−12 sen φ, 12 − 8t]]
        M̂_cons = [[420 − 216t, 26 sen φ], [26 sen φ, 2 + 12t]],  m = ρAh/420
        M̂_HRZ  = diag(1, 1/39),                                  m = ρAh

    y det K̂ = 192 t² (sin cancelación). Desarrollando la raíz menor:
    consistente λ^h/λ = 1 + φ⁴/720 + O(φ⁶), es decir
    **ω^h/ω = 1 + φ⁴/1440 + O(φ⁶)** (orden 4, por encima); HRZ
    **ω^h/ω = 1 − φ²/78 + O(φ⁴)** (orden 2, por debajo): la masa de giro
    m_θ actúa como una inercia rotatoria espuria de radio de giro² = h²/39, que
    tiende a cero con la malla.
(2) **Torsión empotrada-libre** (cadena de elementos lineales, θ_j = sen(jθ),
    θ = (2j − 1)π/(2n), con la misma extensión especular en el extremo libre):
    consistente ω^h = (c/h)·√(12t/(3 − 2t)), HRZ ω^h = (2c/h)·√t, con
    t = sen²(θ/2) y c = √(GJ/(ρJ_p)). Desarrollos: consistente
    ω^h/ω = 1 + θ²/24 + 0.00052 θ⁴, HRZ ω^h/ω = sen(θ/2)/(θ/2) = 1 − θ²/24 + θ⁴/1920.
(3) **Voladizo**: sin forma cerrada; referencia discreta ensamblada aquí con las
    fórmulas aprobadas (rigidez hermítica clásica y masa de arriba).
(4) **Rayleigh-Ritz**: la masa consistente con el espacio hermítico conforme
    (C¹) da ω_k^h ≥ ω_k (min-max de Courant-Fischer) y, con mallas anidadas
    (n → 2n), ω_k^h decrece con n. La HRZ no es una cota: que quede por debajo
    lo dice el desarrollo (1) y lo comprueba el voladizo.

Por qué estos casos detectan la masa anterior: el bloque ρIL/6·[[2, 1], [1, 2]]
en (θ1, θ2) es la inercia rotatoria de una viga de Rayleigh; la consistente
convergería a ω_R = ω/√(1 + (kπr/L)²) con r² = I/A, −1.16 % (k = 1) y
−4.4 % (k = 2) por DEBAJO de Euler-Bernoulli, y violaría la cota de
Rayleigh-Ritz en las mallas finas. La concentrada anterior (ρIL/2 por nodo)
tiende al mismo límite.

Tolerancias — derivación
------------------------
- Matrices elementales: fórmulas de pocas operaciones más TᵀMT con T ortogonal:
  redondeo ~1e-15; se exige |ΔM_ij|/√(M_ii M_jj) ≤ 1e-12 (normalización
  dimensionalmente coherente: kg, kg·m y kg·m² en la misma matriz).
- Frecuencias discretas frente a las identidades (1)-(3): problema
  generalizado denso resuelto como M φ = μ K φ (μ = 1/ω², Cholesky de K), cuya
  μ mayor sale con error relativo ~ε; el condicionamiento de K (∝ n⁴) lo sube a
  ≤ 1e-11 en n = 32. Tolerancia 1e-10.
- Órdenes medidos p = log₂(e_n/e_2n), e = ω^h/ω − 1. La desviación de p respecto
  del orden nominal es O(φ²) y se calcula con el símbolo exacto (1): consistente
  0.074 como máximo (k = 2, par 4 → 8) y 0.005 en el par 16 → 32; HRZ 0.121 y
  0.006. Bandas: |p − 4| ≤ 0.1 y |p − 2| ≤ 0.15 en todos los pares y ≤ 0.02 en
  el par fino. El redondeo (≤ 1e-11 en ω frente a e ≥ 6e-8) mueve p menos de
  1e-3. Voladizo: 3.973/3.993/3.998 y 1.963/1.991/1.998 con la referencia (3):
  mismas bandas.
- Coeficientes principales en n = 32: e·1440/φ⁴ − 1 = −3e-4 (k = 1) y −1.1e-3
  (k = 2); e·78/φ² + 1 = −3e-4 y −1.4e-3 (O(φ²) de (1)): banda 0.01. Voladizo
  consistente, con φ = β_1 h: 0.975, 0.994, 0.998, 1.000 para n = 4…32 (el
  término de dispersión domina; los contornos no añaden términos de menor
  orden porque el elemento hermítico es nodalmente exacto): banda 0.05.
  Torsión: desviaciones ≤ 0.002 en n = 4 por (2): banda 0.01.
"""
from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.linalg import eigh

import solidum  # noqa: F401  — autodescubrimiento de los registros
from solidum.core.domain import Domain
from solidum.elements.frame import Frame2DEuler, Frame2DEulerCorot
from solidum.elements.frame3d import Frame3D
from solidum.materials.elastic import Elastic1D
from solidum.math.assembly import Assembler
from solidum.math.solvers import ModalSolver
from solidum.sections import BarLayer, FiberSection, RectangularPatch


# =============================================================================
# Datos físicos
# =============================================================================

E = 2.1e11             # [Pa]
RHO = 7850.0           # [kg/m³]
A = 0.02               # [m²]
I = 3.0e-4             # [m⁴]
L = 2.5                # [m] longitud del elemento aislado y de las vigas
NU = 0.3
G = E / (2.0 * (1.0 + NU))
RHO_A = RHO * A

# Frame3D: I_y ≠ I_z y J ≠ I_y + I_z.
IY = 1.0e-4
IZ = 3.0e-4
J = 2.2e-4
RHO_JP = RHO * (IY + IZ)

# Sección por fibras de área A (b h = A) y, en el continuo, I_z = I.
H_SEC = math.sqrt(12.0 * I / A)
B_SEC = A / H_SEC
NY_SEC, NZ_SEC = 20, 6
IZ_SEC = B_SEC * H_SEC**3 / 12.0 * (1.0 - 1.0 / NY_SEC**2)   # Σ y² A de las celdas
IY_SEC = H_SEC * B_SEC**3 / 12.0 * (1.0 - 1.0 / NZ_SEC**2)   # Σ z² A de las celdas
J_SEC = 1.4e-5
GJ_SEC = G * J_SEC

# Sección compuesta (densidades distintas): concreto 0.3 × 0.5 y dos capas de
# 3 barras de acero en y = ±0.2, z ∈ {−0.1, 0, 0.1}.
RHO_C, E_C = 2400.0, 30.0e9
B_C, H_C = 0.3, 0.5
NY_C, NZ_C = 10, 6
A_BAR = 3.1e-4
Y_BAR, Z_BAR = 0.2, 0.1

LUMPINGS = ("consistent", "lumped")
MESHES = (4, 8, 16, 32)

TOL_MATRIX = 1.0e-12       # matrices elementales (redondeo de TᵀMT)
TOL_DISCRETE = 1.0e-10     # frecuencias frente a las identidades exactas del modelo discreto

# Geometrías del elemento aislado (longitud L en todos los casos).
GEOM_2D = {
    "alineado": ((0.7, -0.4), (0.7 + L, -0.4)),
    "oblicuo": ((0.4, -0.2), (0.4 + 0.6 * L, -0.2 + 0.8 * L)),
}
_DIR_3D = np.array([2.0, 3.0, 6.0]) / 7.0
_P1_3D = (0.5, -0.3, 0.2)
# Alineado con +x y ref_vector = +y global: ejes locales = globales (spec §6).
REF_ALIGNED = [0.0, 1.0, 0.0]
GEOM_3D = {
    "alineado": (_P1_3D, (0.5 + L, -0.3, 0.2), REF_ALIGNED),
    "oblicuo": (_P1_3D, tuple(np.array(_P1_3D) + L * _DIR_3D), [0.0, 0.0, 1.0]),
}


# =============================================================================
# Matrices de referencia (física aprobada)
# =============================================================================

def _check_lumping(lumping: str) -> None:
    if lumping not in LUMPINGS:
        raise ValueError(lumping)


def mass_axial(rho_a: float, length: float, lumping: str) -> np.ndarray:
    """Barra de 2 nodos (también la torsión con ρJ_p en lugar de ρA)."""
    _check_lumping(lumping)
    if lumping == "consistent":
        return rho_a * length / 6.0 * np.array([[2.0, 1.0], [1.0, 2.0]])
    return rho_a * length / 2.0 * np.eye(2)


def mass_bending(rho_a: float, length: float, lumping: str) -> np.ndarray:
    """Flexión en (v1, θ1, v2, θ2), plano xy: hermítica sin inercia rotatoria o HRZ."""
    _check_lumping(lumping)
    h = length
    if lumping == "consistent":
        return rho_a * h / 420.0 * np.array([
            [156.0, 22.0 * h, 54.0, -13.0 * h],
            [22.0 * h, 4.0 * h**2, 13.0 * h, -3.0 * h**2],
            [54.0, 13.0 * h, 156.0, -22.0 * h],
            [-13.0 * h, -3.0 * h**2, -22.0 * h, 4.0 * h**2],
        ])
    m_t = rho_a * h / 2.0          # 156/312 de ρAh
    m_r = rho_a * h**3 / 78.0      # 4h²/312 de ρAh
    return np.diag([m_t, m_r, m_t, m_r])


# (w1, θy1, w2, θy2): w' = −θy ⇒ misma matriz con θ → −θ.
S_XZ = np.diag([1.0, -1.0, 1.0, -1.0])


def local_mass_2d(rho_a: float, length: float, lumping: str) -> np.ndarray:
    """6×6 local en [u1, v1, θ1, u2, v2, θ2]."""
    M = np.zeros((6, 6))
    M[np.ix_([0, 3], [0, 3])] = mass_axial(rho_a, length, lumping)
    bend = [1, 2, 4, 5]
    M[np.ix_(bend, bend)] = mass_bending(rho_a, length, lumping)
    return M


def local_mass_3d(rho_a: float, rho_jp: float, length: float, lumping: str) -> np.ndarray:
    """12×12 local en [u, v, w, θx, θy, θz] por nodo."""
    M = np.zeros((12, 12))
    axial, torsion = [0, 6], [3, 9]
    plane_xy, plane_xz = [1, 5, 7, 11], [2, 4, 8, 10]
    M[np.ix_(axial, axial)] = mass_axial(rho_a, length, lumping)
    M[np.ix_(torsion, torsion)] = mass_axial(rho_jp, length, lumping)
    M_b = mass_bending(rho_a, length, lumping)
    M[np.ix_(plane_xy, plane_xy)] = M_b
    M[np.ix_(plane_xz, plane_xz)] = S_XZ @ M_b @ S_XZ
    return M


def transformation_2d(p1, p2) -> np.ndarray:
    """T de la spec Frame2DEuler §7: u_local = T u_global."""
    d = np.asarray(p2, float) - np.asarray(p1, float)
    c, s = d / np.linalg.norm(d)
    R = np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])
    return np.kron(np.eye(2), R)


def local_axes_3d(p1, p2, ref) -> np.ndarray:
    """λ de la spec Frame3D §6-§7: filas x̂, ŷ, ẑ en coordenadas globales."""
    p1, p2, ref = (np.asarray(v, float) for v in (p1, p2, ref))
    x = (p2 - p1) / np.linalg.norm(p2 - p1)
    y = ref - ref.dot(x) * x
    y /= np.linalg.norm(y)
    return np.vstack([x, y, np.cross(x, y)])


def transformation_3d(p1, p2, ref) -> np.ndarray:
    """T 12×12 = blkdiag(λ, λ, λ, λ)."""
    return np.kron(np.eye(4), local_axes_3d(p1, p2, ref))


def assert_same_mass(M, M_ref, rtol: float, label: str) -> None:
    """M = M_ref con |ΔM_ij|/√(M_ii M_jj) ≤ rtol y M simétrica."""
    M = np.asarray(M, dtype=float)
    assert M.shape == M_ref.shape, f"{label}: forma {M.shape} ≠ {M_ref.shape}"
    d = np.sqrt(np.abs(np.diag(M_ref)))
    assert np.all(d > 0.0), f"{label}: diagonal de referencia con ceros"
    scale = np.outer(d, d)
    asym = np.abs(M - M.T) / scale
    assert asym.max() <= rtol, f"{label}: M no simétrica ({asym.max():.3e})"
    dM = np.abs(M - M_ref) / scale
    i, j = np.unravel_index(np.argmax(dM), dM.shape)
    assert dM.max() <= rtol, (
        f"{label}: max |ΔM_ij|/√(M_ii M_jj) = {dM.max():.3e} > {rtol:.0e} en ({i}, {j}): "
        f"M = {M[i, j]:.10e}, fórmula = {M_ref[i, j]:.10e}")


# =============================================================================
# Construcción de elementos y secciones
# =============================================================================

def steel() -> Elastic1D:
    return Elastic1D(E=E, density=RHO)


def element_2d(cls, p1, p2, **kwargs):
    """Un marco 2D aislado entre ``p1`` y ``p2``."""
    dom = Domain()
    n1 = dom.add_node(1, list(p1))
    n2 = dom.add_node(2, list(p2))
    elem = cls(1, [n1, n2], **kwargs)
    dom.add_element(elem)
    dom.generate_equation_numbers(verbose=False)
    return elem


def element_3d(p1, p2, **kwargs):
    """Un Frame3D aislado entre ``p1`` y ``p2``."""
    dom = Domain()
    n1 = dom.add_node(1, list(p1))
    n2 = dom.add_node(2, list(p2))
    elem = Frame3D(1, [n1, n2], **kwargs)
    dom.add_element(elem)
    dom.generate_equation_numbers(verbose=False)
    return elem


def rect_section(nz: int, GJ: float | None = None) -> FiberSection:
    """Rectángulo b × h de acero (b h = A), NY_SEC × nz celdas."""
    patch = RectangularPatch(steel(), y=(-H_SEC / 2.0, H_SEC / 2.0),
                             z=(-B_SEC / 2.0, B_SEC / 2.0), ny=NY_SEC, nz=nz)
    return FiberSection(patches=[patch], GJ=GJ)


def composite_section(GJ: float | None = None):
    """Concreto + 6 barras de acero; devuelve (sección, ρA, ρJ_p) en forma cerrada.

    ρA = ρ_c b h + ρ_s 6 A_b. ρJ_p = Σρ_i A_i (y_i² + z_i²): celdas del
    rectángulo con la fibra en el centro, Σ y² A = b h³/12 (1 − 1/n_y²) y
    Σ z² A = h b³/12 (1 − 1/n_z²); barras Σ y² = 6·0.2², Σ z² = 4·0.1².
    """
    concrete = Elastic1D(E=E_C, density=RHO_C)
    bars = steel()
    sec = FiberSection(
        patches=[
            RectangularPatch(concrete, y=(-H_C / 2, H_C / 2), z=(-B_C / 2, B_C / 2), ny=NY_C, nz=NZ_C),
            BarLayer(bars, start=(-Y_BAR, -Z_BAR), end=(-Y_BAR, Z_BAR), n=3, area=A_BAR),
            BarLayer(bars, start=(Y_BAR, -Z_BAR), end=(Y_BAR, Z_BAR), n=3, area=A_BAR),
        ],
        GJ=GJ,
    )
    rho_a = RHO_C * B_C * H_C + RHO * 6 * A_BAR
    rho_jp = (RHO_C * (B_C * H_C**3 / 12.0 * (1.0 - 1.0 / NY_C**2)
                       + H_C * B_C**3 / 12.0 * (1.0 - 1.0 / NZ_C**2))
              + RHO * A_BAR * (6 * Y_BAR**2 + 4 * Z_BAR**2))
    return sec, rho_a, rho_jp


# Caminos 2D: (ρA, EI) de cada uno.
PATHS_2D = ("euler", "corot", "seccion")
PATH_PROPS_2D = {
    "euler": (RHO_A, E * I),
    "corot": (RHO_A, E * I),
    "seccion": (RHO * B_SEC * H_SEC, E * IZ_SEC),
}


def make_frame_2d(path: str, eid: int, nodes, section=None):
    if path == "euler":
        return Frame2DEuler(eid, nodes, material=steel(), A=A, I=I)
    if path == "corot":
        return Frame2DEulerCorot(eid, nodes, material=steel(), A=A, I=I)
    return Frame2DEuler(eid, nodes, section=section if section is not None else rect_section(1), n_ip=3)


# =============================================================================
# (a) Matrices elementales frente a las fórmulas
# =============================================================================

@pytest.mark.parametrize("lumping", LUMPINGS)
@pytest.mark.parametrize("geometry", list(GEOM_2D))
@pytest.mark.parametrize("cls", [Frame2DEuler, Frame2DEulerCorot], ids=lambda c: c.__name__)
def test_a_matriz_elemental_2d(cls, geometry, lumping):
    """M global = Tᵀ M_local T con la consistente sin inercia rotatoria o la HRZ."""
    p1, p2 = GEOM_2D[geometry]
    elem = element_2d(cls, p1, p2, material=steel(), A=A, I=I)
    M = elem.compute_mass_matrix(lumping=lumping)
    T = transformation_2d(p1, p2)
    M_ref = T.T @ local_mass_2d(RHO_A, L, lumping) @ T
    assert_same_mass(M, M_ref, TOL_MATRIX, f"{cls.__name__} {geometry} {lumping}")


@pytest.mark.parametrize("lumping", LUMPINGS)
@pytest.mark.parametrize("cls", [Frame2DEuler, Frame2DEulerCorot], ids=lambda c: c.__name__)
def test_a_bloque_de_giro_2d_sin_inercia_rotatoria(cls, lumping):
    """Bloque (θ1, θ2) = ρAL³/420·[[4, −3], [−3, 4]] (o ρAL³/78 por nodo) e independiente de I.

    La masa anterior sumaba ρIL/6·[[2, 1], [1, 2]] (o ρIL/2 por nodo): con
    I = 3e-4 es un 8 % de M_θθ consistente, y dependería de I.
    """
    p1, p2 = GEOM_2D["alineado"]
    M = element_2d(cls, p1, p2, material=steel(), A=A, I=I).compute_mass_matrix(lumping=lumping)
    k = RHO_A * L**3
    expected = (np.array([[4.0, -3.0], [-3.0, 4.0]]) * k / 420.0 if lumping == "consistent"
                else np.eye(2) * k / 78.0)
    block = M[np.ix_([2, 5], [2, 5])]
    np.testing.assert_allclose(block, expected, rtol=TOL_MATRIX, atol=TOL_MATRIX * k / 78.0,
                               err_msg=f"{cls.__name__} {lumping}: bloque de giro con inercia rotatoria")
    M_other_I = element_2d(cls, p1, p2, material=steel(), A=A, I=I / 100.0).compute_mass_matrix(
        lumping=lumping)
    assert_same_mass(M_other_I, M, TOL_MATRIX, f"{cls.__name__} {lumping}: M depende de I")


@pytest.mark.parametrize("lumping", LUMPINGS)
@pytest.mark.parametrize("geometry", list(GEOM_2D))
@pytest.mark.parametrize("cls", [Frame2DEuler, Frame2DEulerCorot], ids=lambda c: c.__name__)
def test_a_movimientos_rigidos_2d(cls, geometry, lumping):
    """uᵀMu de una traslación = ρAL|d|²; de un giro rígido unitario alrededor del origen:

    - consistente: ρA ∫|r|² ds = ρA (L|r1|² + L² r1·x̂ + L³/3), exacta porque el
      campo rígido (u lineal, v lineal con v' = θ) está en el espacio de
      interpolación y la masa no tiene inercia rotatoria;
    - HRZ: ρAL/2 (|r1|² + |r2|²) + 2·ρAL³/78.
    """
    p1, p2 = (np.asarray(p, float) for p in GEOM_2D[geometry])
    M = element_2d(cls, tuple(p1), tuple(p2), material=steel(), A=A, I=I).compute_mass_matrix(
        lumping=lumping)
    d = np.array([math.cos(0.3), math.sin(0.3)])
    u_t = np.array([d[0], d[1], 0.0, d[0], d[1], 0.0])
    np.testing.assert_allclose(u_t @ M @ u_t, RHO_A * L, rtol=1e-13,
                               err_msg="energía cinética de traslación ≠ ρAL")
    u_r = np.array([-p1[1], p1[0], 1.0, -p2[1], p2[0], 1.0])
    x_hat = (p2 - p1) / L
    if lumping == "consistent":
        expected = RHO_A * (L * p1.dot(p1) + L**2 * p1.dot(x_hat) + L**3 / 3.0)
    else:
        expected = RHO_A * L / 2.0 * (p1.dot(p1) + p2.dot(p2)) + 2.0 * RHO_A * L**3 / 78.0
    np.testing.assert_allclose(u_r @ M @ u_r, expected, rtol=1e-12,
                               err_msg=f"{cls.__name__} {geometry} {lumping}: energía de giro rígido")


def frame3d_material(p1, p2, ref, **overrides):
    props = dict(material=steel(), A=A, Iy=IY, Iz=IZ, J=J, nu=NU, ref_vector=ref)
    props.update(overrides)
    return element_3d(p1, p2, **props)


@pytest.mark.parametrize("lumping", LUMPINGS)
@pytest.mark.parametrize("geometry", list(GEOM_3D))
def test_a_matriz_elemental_3d(geometry, lumping):
    """Frame3D: axial, flexión en xy y xz (acoplamientos con signo opuesto) y torsión con ρJ_p."""
    p1, p2, ref = GEOM_3D[geometry]
    M = frame3d_material(p1, p2, ref).compute_mass_matrix(lumping=lumping)
    T = transformation_3d(p1, p2, ref)
    M_ref = T.T @ local_mass_3d(RHO_A, RHO_JP, L, lumping) @ T
    assert_same_mass(M, M_ref, TOL_MATRIX, f"Frame3D {geometry} {lumping}")


def test_a_signo_del_acoplamiento_en_el_plano_xz():
    """Elemento alineado (ejes locales = globales): entradas del acoplamiento traslación-giro.

    Plano xy: M[v1, θz1] = +22L, M[v1, θz2] = −13L, M[θz1, v2] = +13L (× ρAL/420).
    Plano xz (w' = −θy): M[w1, θy1] = −22L, M[w1, θy2] = +13L, M[θy1, w2] = −13L.
    Torsión: M[θx1, θx1] = ρJ_pL/3, M[θx1, θx2] = ρJ_pL/6.
    """
    p1, p2, ref = GEOM_3D["alineado"]
    M = frame3d_material(p1, p2, ref).compute_mass_matrix(lumping="consistent")
    k = RHO_A * L / 420.0
    checks = {
        (1, 5): 22.0 * L * k, (1, 11): -13.0 * L * k, (5, 7): 13.0 * L * k,
        (2, 4): -22.0 * L * k, (2, 10): 13.0 * L * k, (4, 8): -13.0 * L * k,
        (4, 4): 4.0 * L**2 * k, (5, 5): 4.0 * L**2 * k, (4, 10): -3.0 * L**2 * k,
        (3, 3): RHO_JP * L / 3.0, (3, 9): RHO_JP * L / 6.0,
        (1, 4): 0.0, (2, 5): 0.0, (3, 4): 0.0, (4, 5): 0.0,
    }
    for (i, j), value in checks.items():
        atol = TOL_MATRIX * math.sqrt(M[i, i] * M[j, j])
        assert abs(M[i, j] - value) <= atol + TOL_MATRIX * abs(value), (
            f"M[{i}, {j}] = {M[i, j]:.10e}, esperado {value:.10e}")


@pytest.mark.parametrize("lumping", LUMPINGS)
def test_a_masa_3d_solo_depende_de_rhoA_y_de_Iy_mas_Iz(lumping):
    """Sin inercia rotatoria de flexión, M no depende de la orientación de la sección
    alrededor del eje (ref_vector), ni de J, ni de (I_y, I_z) salvo por I_y + I_z.

    Invariancia: el bloque traslacional es m_ax x̂x̂ᵀ + m_tr (1 − x̂x̂ᵀ), el
    acoplamiento v θz − w θy = x̂·(u × θ) es un producto mixto y los giros de
    flexión tienen la misma masa en los dos planos: nada depende de (ŷ, ẑ).
    """
    p1, p2, ref = GEOM_3D["oblicuo"]
    M_ref = frame3d_material(p1, p2, ref).compute_mass_matrix(lumping=lumping)
    variants = {
        "ref_vector = (1, −1, 0.3)": frame3d_material(p1, p2, [1.0, -1.0, 0.3]),
        "I_y = I_z = 2e-4": frame3d_material(p1, p2, ref, Iy=2.0e-4, Iz=2.0e-4),
        "J = 5e-5": frame3d_material(p1, p2, ref, J=5.0e-5),
    }
    for label, elem in variants.items():
        assert_same_mass(elem.compute_mass_matrix(lumping=lumping), M_ref, TOL_MATRIX,
                         f"Frame3D {lumping}: M cambia con {label}")


@pytest.mark.parametrize("lumping", LUMPINGS)
@pytest.mark.parametrize("geometry", list(GEOM_3D))
def test_a_movimientos_rigidos_3d(geometry, lumping):
    """Giro rígido ω alrededor del origen: u_i = ω × r_i, θ_i = ω.

    - consistente: ρA ∫|ω × r|² ds + ρJ_p L (ω·x̂)², con
      ∫|ω × (r1 + s x̂)|² ds = L|ω×r1|² + L² (ω×r1)·(ω×x̂) + L³/3 |ω×x̂|²;
    - HRZ: ρAL/2 (|ω×r1|² + |ω×r2|²) + ρJ_p L (ω·x̂)² + ρAL³/39 (|ω|² − (ω·x̂)²).
    """
    p1, p2, ref = GEOM_3D[geometry]
    r1, r2 = np.asarray(p1, float), np.asarray(p2, float)
    M = frame3d_material(p1, p2, ref).compute_mass_matrix(lumping=lumping)
    d = np.array([0.36, -0.48, 0.8])
    u_t = np.concatenate([d, np.zeros(3), d, np.zeros(3)])
    np.testing.assert_allclose(u_t @ M @ u_t, RHO_A * L, rtol=1e-13,
                               err_msg="energía cinética de traslación ≠ ρAL")
    w = np.array([0.3, -0.5, 0.8])
    u_r = np.concatenate([np.cross(w, r1), w, np.cross(w, r2), w])
    x_hat = (r2 - r1) / L
    wx = w.dot(x_hat)
    a, b = np.cross(w, r1), np.cross(w, x_hat)
    if lumping == "consistent":
        expected = RHO_A * (L * a.dot(a) + L**2 * a.dot(b) + L**3 / 3.0 * b.dot(b)) + RHO_JP * L * wx**2
    else:
        c2 = np.cross(w, r2)
        expected = (RHO_A * L / 2.0 * (a.dot(a) + c2.dot(c2)) + RHO_JP * L * wx**2
                    + RHO_A * L**3 / 39.0 * (w.dot(w) - wx**2))
    np.testing.assert_allclose(u_r @ M @ u_r, expected, rtol=1e-12,
                               err_msg=f"Frame3D {geometry} {lumping}: energía de giro rígido")


# =============================================================================
# (b) Igualdad entre caminos (material, corotacional, sección de igual ρA)
# =============================================================================

@pytest.mark.parametrize("lumping", LUMPINGS)
def test_b_igualdad_euler_corot_seccion_2d(lumping):
    """Frame2DEuler ↔ Frame2DEulerCorot ↔ Frame2DEuler con sección de ρA = ρ b h = ρA."""
    p1, p2 = GEOM_2D["oblicuo"]
    M_euler = element_2d(Frame2DEuler, p1, p2, material=steel(), A=A, I=I).compute_mass_matrix(
        lumping=lumping)
    M_corot = element_2d(Frame2DEulerCorot, p1, p2, material=steel(), A=A, I=I).compute_mass_matrix(
        lumping=lumping)
    M_sec = element_2d(Frame2DEuler, p1, p2, section=rect_section(1), n_ip=3).compute_mass_matrix(
        lumping=lumping)
    assert_same_mass(M_corot, M_euler, TOL_MATRIX, f"Corot ≠ Euler ({lumping})")
    assert_same_mass(M_sec, M_euler, TOL_MATRIX, f"sección ≠ material ({lumping})")


@pytest.mark.parametrize("lumping", LUMPINGS)
@pytest.mark.parametrize("cls", [Frame2DEuler, Frame2DEulerCorot], ids=lambda c: c.__name__)
def test_b_seccion_compuesta_2d_rhoA_suma_de_fibras(cls, lumping):
    """Densidades distintas por fibra: ρA = Σρ_i A_i = 374.601 kg/m en las fórmulas."""
    sec, rho_a, _ = composite_section()
    p1, p2 = GEOM_2D["oblicuo"]
    M = element_2d(cls, p1, p2, section=sec, n_ip=3).compute_mass_matrix(lumping=lumping)
    T = transformation_2d(p1, p2)
    assert_same_mass(M, T.T @ local_mass_2d(rho_a, L, lumping) @ T, TOL_MATRIX,
                     f"{cls.__name__} sección compuesta {lumping}")


@pytest.mark.parametrize("lumping", LUMPINGS)
def test_b_igualdad_material_seccion_3d(lumping):
    """Frame3D con sección = con material de A = b h, I_y = I_y,n, I_z = I_z,n, y las fórmulas
    con ρJ_p = ρ(I_y,n + I_z,n); sección compuesta con ρA y ρJ_p de las fibras."""
    p1, p2, ref = GEOM_3D["oblicuo"]
    T = transformation_3d(p1, p2, ref)
    M_sec = element_3d(p1, p2, section=rect_section(NZ_SEC, GJ=GJ_SEC), n_ip=3,
                       ref_vector=ref).compute_mass_matrix(lumping=lumping)
    M_mat = frame3d_material(p1, p2, ref, A=B_SEC * H_SEC, Iy=IY_SEC, Iz=IZ_SEC,
                             J=J_SEC).compute_mass_matrix(lumping=lumping)
    M_ref = T.T @ local_mass_3d(RHO * B_SEC * H_SEC, RHO * (IY_SEC + IZ_SEC), L, lumping) @ T
    assert_same_mass(M_sec, M_ref, TOL_MATRIX, f"Frame3D sección {lumping}")
    assert_same_mass(M_mat, M_sec, TOL_MATRIX, f"Frame3D material ≠ sección {lumping}")

    sec, rho_a, rho_jp = composite_section(GJ=GJ_SEC)
    M_comp = element_3d(p1, p2, section=sec, n_ip=3, ref_vector=ref).compute_mass_matrix(lumping=lumping)
    assert_same_mass(M_comp, T.T @ local_mass_3d(rho_a, rho_jp, L, lumping) @ T, TOL_MATRIX,
                     f"Frame3D sección compuesta {lumping}")


# =============================================================================
# Herramientas modales
# =============================================================================

def _dense(X) -> np.ndarray:
    return X.toarray() if hasattr(X, "toarray") else np.asarray(X, dtype=float)


def model_frequencies(dom: Domain, lumping: str, n_modes: int, return_modes: bool = False):
    """ω más bajas de K φ = ω² M φ tras eliminar los DOF prescritos.

    Ensamblaje y reducción como ``ModalSolver`` (Assembler.assemble_system,
    assemble_mass_matrix, reduce_pair); el problema denso se resuelve como
    M φ = μ K φ (μ = 1/ω²): con Cholesky de K la μ mayor sale con error
    relativo ~ε, lo que el orden medido en n = 32 necesita (e ~ 6e-8).
    """
    asm = Assembler(dom)
    asm.assemble_system()
    K = asm.K_global
    M = asm.assemble_mass_matrix(lumping=lumping)
    K_red, M_red, T = asm.reduce_pair(K, M)
    mu, vec = eigh(_dense(M_red), _dense(K_red))
    order = np.argsort(mu)[::-1][:n_modes]
    omegas = 1.0 / np.sqrt(mu[order])
    if return_modes:
        return omegas, np.asarray(T @ vec[:, order])
    return omegas


def ss_symbol_frequency(k: int, n: int, EI: float, rho_a: float, length: float, lumping: str) -> float:
    """ω exacta del modo k de la viga biapoyada de n elementos (identidad (1) del encabezado).

    Raíz menor de det(K̂ − λ̃ M̂) = a2 λ̃² + a1 λ̃ + a0 = 0 en la forma estable
    λ̃ = 2 a0/(−a1 + √(a1² − 4 a2 a0)), con a0 = det K̂ = 192 t².
    """
    _check_lumping(lumping)
    h = length / n
    phi = k * math.pi / n
    t = math.sin(phi / 2.0) ** 2
    sin2 = 4.0 * t * (1.0 - t)                      # sen² φ
    Kvv, Kxx = 48.0 * t, 12.0 - 8.0 * t
    a0 = 192.0 * t * t
    if lumping == "consistent":
        Mvv, Mxx = 420.0 - 216.0 * t, 2.0 + 12.0 * t
        KvxMvx = -12.0 * 26.0 * sin2                 # (−12 sen φ)(26 sen φ)
        Mvx2 = 676.0 * sin2
        m = rho_a * h / 420.0
    else:
        Mvv, Mxx, KvxMvx, Mvx2 = 1.0, 1.0 / 39.0, 0.0, 0.0
        m = rho_a * h
    a2 = Mvv * Mxx - Mvx2
    a1 = -(Kvv * Mxx + Kxx * Mvv - 2.0 * KvxMvx)
    lam = 2.0 * a0 / (-a1 + math.sqrt(a1 * a1 - 4.0 * a2 * a0))
    return math.sqrt(lam * EI / h**3 / m)


def reference_beam_frequencies(n: int, lumping: str, support: str, EI: float, rho_a: float,
                               length: float, n_modes: int) -> np.ndarray:
    """Viga de n elementos ensamblada con las fórmulas aprobadas (DOF v, θ por nodo)."""
    h = length / n
    ke = EI / h**3 * np.array([
        [12.0, 6.0 * h, -12.0, 6.0 * h],
        [6.0 * h, 4.0 * h**2, -6.0 * h, 2.0 * h**2],
        [-12.0, -6.0 * h, 12.0, -6.0 * h],
        [6.0 * h, 2.0 * h**2, -6.0 * h, 4.0 * h**2],
    ])
    me = mass_bending(rho_a, h, lumping)
    nd = 2 * (n + 1)
    K, M = np.zeros((nd, nd)), np.zeros((nd, nd))
    for e in range(n):
        idx = np.arange(2 * e, 2 * e + 4)
        K[np.ix_(idx, idx)] += ke
        M[np.ix_(idx, idx)] += me
    fixed = {0, 2 * n} if support == "biapoyada" else {0, 1}
    free = [i for i in range(nd) if i not in fixed]
    mu = eigh(M[np.ix_(free, free)], K[np.ix_(free, free)], eigvals_only=True)
    return 1.0 / np.sqrt(mu[::-1][:n_modes])


def torsion_chain_frequency(j: int, n: int, c: float, length: float, lumping: str) -> float:
    """ω exacta del modo j de la barra empotrada-libre de n elementos (identidad (2))."""
    _check_lumping(lumping)
    h = length / n
    t = math.sin((2 * j - 1) * math.pi / (4.0 * n)) ** 2    # sen²(θ/2)
    if lumping == "consistent":
        return c / h * math.sqrt(12.0 * t / (3.0 - 2.0 * t))
    return 2.0 * c / h * math.sqrt(t)


def cantilever_root() -> float:
    """β_1 L: primera raíz de cos x cosh x + 1 = 0 (Newton)."""
    x = 1.875
    for _ in range(100):
        f = math.cos(x) * math.cosh(x) + 1.0
        df = -math.sin(x) * math.cosh(x) + math.cos(x) * math.sinh(x)
        dx = f / df
        x -= dx
        if abs(dx) < 1e-15:
            break
    return x


def observed_orders(errors) -> list[float]:
    """p = log₂(e_n/e_2n) sobre mallas que se duplican."""
    e = np.abs(np.asarray(errors, dtype=float))
    return [float(math.log2(e[i] / e[i + 1])) for i in range(len(e) - 1)]


def check_convergence(errors, lumping: str, p_nominal: float, band_all: float, band_fine: float,
                      label: str) -> None:
    """Signo del error (consistente > 0, HRZ < 0), monotonía y orden medido."""
    errors = np.asarray(errors, dtype=float)
    if lumping == "consistent":
        assert np.all(errors > 0.0), f"{label}: la consistente no es cota superior: e = {errors}"
    else:
        assert np.all(errors < 0.0), f"{label}: la HRZ no queda por debajo: e = {errors}"
    assert np.all(np.diff(np.abs(errors)) < 0.0), f"{label}: convergencia no monótona: e = {errors}"
    p = observed_orders(errors)
    assert all(abs(q - p_nominal) <= band_all for q in p), (
        f"{label}: órdenes medidos {np.round(p, 4)} fuera de {p_nominal} ± {band_all}")
    assert abs(p[-1] - p_nominal) <= band_fine, (
        f"{label}: orden del par fino {p[-1]:.4f} fuera de {p_nominal} ± {band_fine}")


# =============================================================================
# (c) Viga biapoyada: convergencia de ω_1 y ω_2
# =============================================================================

def build_beam_2d(path: str, n: int, support: str) -> Domain:
    """Viga de longitud L a lo largo de x con n marcos 2D.

    ``ux = 0`` en todos los nodos: en una viga recta la parte axial está
    desacoplada de la flexión en K y en M, así que eliminarla no cambia las
    frecuencias de flexión y evita modos axiales intercalados.
    """
    dom = Domain()
    nodes = [dom.add_node(i + 1, [L * i / n, 0.0]) for i in range(n + 1)]
    section = rect_section(1) if path == "seccion" else None
    for i in range(n):
        dom.add_element(make_frame_2d(path, i + 1, [nodes[i], nodes[i + 1]], section))
    for node in nodes:
        node.fix_dof("ux", 0.0)
    if support == "biapoyada":
        for node in (nodes[0], nodes[-1]):
            node.fix_dof("uy", 0.0)
    else:
        nodes[0].fix_dof("uy", 0.0)
        nodes[0].fix_dof("rz", 0.0)
    dom.generate_equation_numbers(verbose=False)
    return dom


@pytest.mark.parametrize("lumping", LUMPINGS)
@pytest.mark.parametrize("path", PATHS_2D)
def test_c_viga_biapoyada_convergencia(path, lumping):
    """n = 4, 8, 16, 32: identidad con el símbolo, cota de Rayleigh-Ritz, monotonía y orden.

    Consistente: por encima, orden 4, e → φ⁴/1440. HRZ: por debajo, orden 2,
    e → −φ²/78 (φ = kπ/n).
    """
    rho_a, EI = PATH_PROPS_2D[path]
    omega_exact = [(k * math.pi / L) ** 2 * math.sqrt(EI / rho_a) for k in (1, 2)]
    errors = {1: [], 2: []}
    for n in MESHES:
        omegas = model_frequencies(build_beam_2d(path, n, "biapoyada"), lumping, 2)
        for k in (1, 2):
            w_sym = ss_symbol_frequency(k, n, EI, rho_a, L, lumping)
            assert abs(omegas[k - 1] / w_sym - 1.0) <= TOL_DISCRETE, (
                f"{path} {lumping} n={n} k={k}: ω = {omegas[k - 1]:.12e}, "
                f"símbolo exacto del modelo discreto = {w_sym:.12e}")
            errors[k].append(omegas[k - 1] / omega_exact[k - 1] - 1.0)

    band_all = 0.1 if lumping == "consistent" else 0.15
    p_nominal = 4.0 if lumping == "consistent" else 2.0
    for k in (1, 2):
        check_convergence(errors[k], lumping, p_nominal, band_all, 0.02, f"{path} {lumping} k={k}")
        phi = k * math.pi / MESHES[-1]
        e = errors[k][-1]
        if lumping == "consistent":
            ratio, target = e * 1440.0 / phi**4, 1.0
        else:
            ratio, target = e * 78.0 / phi**2, -1.0
        assert abs(ratio - target) <= 0.01, (
            f"{path} {lumping} k={k}: coeficiente principal del error en n = {MESHES[-1]}: "
            f"{ratio:.5f} (esperado {target})")


@pytest.mark.parametrize("lumping", LUMPINGS)
def test_c_modalsolver_reproduce_las_frecuencias(lumping):
    """El pipeline de ``ModalSolver`` da las mismas ω_1..ω_3 que el símbolo exacto (n = 16)."""
    n = 16
    rho_a, EI = PATH_PROPS_2D["euler"]
    result = ModalSolver(Assembler(build_beam_2d("euler", n, "biapoyada")), n_modes=3,
                         lumping=lumping, tolerance=1.0e-12).solve()
    assert result.converged and result.n_modes == 3
    expected = [ss_symbol_frequency(k, n, EI, rho_a, L, lumping) for k in (1, 2, 3)]
    np.testing.assert_allclose(result.frequencies_rad, expected, rtol=1.0e-9,
                               err_msg=f"ModalSolver {lumping}: ω ≠ símbolo exacto")


# =============================================================================
# (d) Voladizo: primera frecuencia
# =============================================================================

@pytest.mark.parametrize("lumping", LUMPINGS)
@pytest.mark.parametrize("path", PATHS_2D)
def test_d_voladizo_primera_frecuencia(path, lumping):
    """ω_1 → 1.875104²·√(EI/(ρA L⁴)): consistente por encima con orden 4, HRZ por debajo con orden 2."""
    beta1 = cantilever_root()
    assert abs(beta1 - 1.875104) < 1e-6, beta1
    rho_a, EI = PATH_PROPS_2D[path]
    omega_exact = beta1**2 * math.sqrt(EI / (rho_a * L**4))
    np.testing.assert_allclose(omega_exact, 1.875104**2 * math.sqrt(EI / (rho_a * L**4)), rtol=1e-6)
    errors = []
    for n in MESHES:
        w1 = model_frequencies(build_beam_2d(path, n, "voladizo"), lumping, 1)[0]
        w_ref = reference_beam_frequencies(n, lumping, "voladizo", EI, rho_a, L, 1)[0]
        assert abs(w1 / w_ref - 1.0) <= TOL_DISCRETE, (
            f"{path} {lumping} n={n}: ω_1 = {w1:.12e}, modelo discreto de referencia = {w_ref:.12e}")
        errors.append(w1 / omega_exact - 1.0)
    p_nominal = 4.0 if lumping == "consistent" else 2.0
    check_convergence(errors, lumping, p_nominal, 0.1, 0.02, f"voladizo {path} {lumping}")
    if lumping == "consistent":
        for n, e in zip(MESHES, errors):
            ratio = e * 1440.0 / (beta1 / n) ** 4
            assert abs(ratio - 1.0) <= 0.05, (
                f"voladizo {path} n={n}: e·1440/(β₁h)⁴ = {ratio:.4f} (término de dispersión)")


# =============================================================================
# (e) Frame3D: flexión en los dos planos y torsión
# =============================================================================

def build_beam_3d(n: int):
    """Frame3D biapoyado a lo largo de x, ejes locales = globales (ref_vector = +y).

    ``ux`` y ``rx`` nulos en todos los nodos (axial y torsión desacoplados de la
    flexión en una barra recta); ``uy = uz = 0`` en los apoyos. Devuelve
    (dominio, nodos).
    """
    dom = Domain()
    nodes = [dom.add_node(i + 1, [L * i / n, 0.0, 0.0]) for i in range(n + 1)]
    mat = steel()
    for i in range(n):
        dom.add_element(Frame3D(i + 1, [nodes[i], nodes[i + 1]], material=mat, A=A, Iy=IY, Iz=IZ,
                                J=J, nu=NU, ref_vector=REF_ALIGNED))
    for node in nodes:
        node.fix_dof("ux", 0.0)
        node.fix_dof("rx", 0.0)
    for node in (nodes[0], nodes[-1]):
        node.fix_dof("uy", 0.0)
        node.fix_dof("uz", 0.0)
    dom.generate_equation_numbers(verbose=False)
    return dom, nodes


@pytest.mark.parametrize("lumping", LUMPINGS)
@pytest.mark.parametrize("n", [8, 16])
def test_e_frame3d_flexion_en_los_dos_planos(n, lumping):
    """I_y = 1e-4 ≠ I_z = 3e-4: plano xz (w, θy) con EI_y y plano xy (v, θz) con EI_z.

    Primeros cuatro modos: xz k=1 (577.5), xy k=1 (1000.3), xz k=2 (2310.1),
    xy k=2 (4001.3) rad/s en el continuo. Cada uno igual al símbolo exacto con
    su EI, del lado correcto del valor analítico, y con el movimiento sólo en
    su plano (los planos están desacoplados en K y en M).
    """
    dom, nodes = build_beam_3d(n)
    omegas, modes = model_frequencies(dom, lumping, 4, return_modes=True)
    expected = sorted(
        [(ss_symbol_frequency(k, n, E * EI_, RHO_A, L, lumping), plane, k, EI_)
         for plane, EI_ in (("xz", IY), ("xy", IZ)) for k in (1, 2)])
    dofs = {name: [node.dofs[name] for node in nodes] for name in ("uy", "uz", "ry", "rz")}
    for i, (w_sym, plane, k, EI_) in enumerate(expected):
        assert abs(omegas[i] / w_sym - 1.0) <= TOL_DISCRETE, (
            f"n={n} {lumping} modo {i + 1} (plano {plane}, k={k}): ω = {omegas[i]:.12e}, "
            f"símbolo = {w_sym:.12e}")
        w_exact = (k * math.pi / L) ** 2 * math.sqrt(E * EI_ / RHO_A)
        e = omegas[i] / w_exact - 1.0
        assert (e > 0.0) if lumping == "consistent" else (e < 0.0), (
            f"n={n} {lumping} plano {plane} k={k}: e = {e:.3e} del lado equivocado")
        phi = modes[:, i]
        in_plane = ("uy", "rz") if plane == "xy" else ("uz", "ry")
        out_plane = ("uz", "ry") if plane == "xy" else ("uy", "rz")
        amp_in = np.abs(phi[dofs[in_plane[0]]]).max()
        amp_out = max(np.abs(phi[dofs[name]]).max() for name in out_plane)
        assert amp_in > 0.0 and amp_out <= 1e-8 * amp_in, (
            f"n={n} {lumping}: el modo {i + 1} (plano {plane}) mezcla planos: "
            f"{amp_out:.3e} fuera del plano frente a {amp_in:.3e}")
    # Mismo φ, misma ρA: ω_xy/ω_xz = √(I_z/I_y) exactamente en el modelo discreto.
    np.testing.assert_allclose(omegas[1] / omegas[0], math.sqrt(IZ / IY), rtol=TOL_DISCRETE)


def build_torsion_bar(path: str, n: int):
    """Frame3D empotrado-libre a lo largo de x con sólo ``rx`` libre (salvo en el empotramiento).

    Devuelve (dominio, c = √(GJ/(ρJ_p))).
    """
    dom = Domain()
    nodes = [dom.add_node(i + 1, [L * i / n, 0.0, 0.0]) for i in range(n + 1)]
    if path == "material":
        props = dict(material=steel(), A=A, Iy=IY, Iz=IZ, J=J, nu=NU)
        c = math.sqrt(G * J / RHO_JP)
    else:
        props = dict(section=rect_section(NZ_SEC, GJ=GJ_SEC), n_ip=3)
        c = math.sqrt(GJ_SEC / (RHO * (IY_SEC + IZ_SEC)))
    for i in range(n):
        dom.add_element(Frame3D(i + 1, [nodes[i], nodes[i + 1]], ref_vector=REF_ALIGNED, **props))
    for node in nodes:
        for dof in ("ux", "uy", "uz", "ry", "rz"):
            node.fix_dof(dof, 0.0)
    nodes[0].fix_dof("rx", 0.0)
    dom.generate_equation_numbers(verbose=False)
    return dom, c


@pytest.mark.parametrize("lumping", LUMPINGS)
@pytest.mark.parametrize("path", ["material", "seccion"])
def test_e_frame3d_torsion_empotrada_libre(path, lumping):
    """ω_1 → (π/2L)·√(GJ/(ρJ_p)), J_p = I_y + I_z: identidad exacta (2), signo, monotonía y orden 2."""
    errors = []
    for n in MESHES:
        dom, c = build_torsion_bar(path, n)
        w1 = model_frequencies(dom, lumping, 1)[0]
        w_chain = torsion_chain_frequency(1, n, c, L, lumping)
        assert abs(w1 / w_chain - 1.0) <= TOL_DISCRETE, (
            f"torsión {path} {lumping} n={n}: ω_1 = {w1:.12e}, cadena exacta = {w_chain:.12e} "
            f"(¿ρJ_p ≠ ρ(I_y + I_z)?)")
        errors.append(w1 / (math.pi / (2.0 * L) * c) - 1.0)
    check_convergence(errors, lumping, 2.0, 0.01, 0.01, f"torsión {path} {lumping}")
    target = 1.0 if lumping == "consistent" else -1.0
    for n, e in zip(MESHES, errors):
        theta = math.pi / (2.0 * n)
        ratio = e * 24.0 / theta**2
        assert abs(ratio - target) <= 0.01, (
            f"torsión {path} {lumping} n={n}: e·24/θ² = {ratio:.5f} (esperado {target})")
