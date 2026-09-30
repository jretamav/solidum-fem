"""Matrices de masa de los marcos (ADR 0009; cambio aprobado el 2026-09-29).

Marcos Euler-Bernoulli (``Frame2DEuler``, ``Frame2DEulerCorot`` y la flexión
de ``Frame3D``), por los caminos con material y con sección, **sin inercia
rotatoria** de la sección (la hipótesis de Euler-Bernoulli la desprecia):

- Consistente: axial ``ρAL/6·[[2,1],[1,2]]`` y flexión hermítica
  ``ρAL/420·H`` en ``(v₁, θ₁, v₂, θ₂)``; en el plano ``xz`` de ``Frame3D``,
  con ``w' = −θ_y``, los términos de acoplamiento traslación-giro cambian de
  signo.
- Concentrada HRZ (Hinton, Rock y Zienkiewicz 1976): ``ρAL/2`` traslacional
  y ``ρAL³/78`` rotacional de flexión por nodo.
- Torsión de ``Frame3D`` sin cambio: ``ρJ_pL/6·[[2,1],[1,2]]`` y ``ρJ_pL/2``,
  ``J_p = I_y + I_z``.

``Frame2DTimoshenko`` conserva la masa con inercia rotatoria (``ρIL/6`` y
``ρIL/2``).

Los valores esperados se escriben desde las fórmulas (y la hermítica se
comprueba además integrando ``∫ρA NᵀN dx`` con las funciones de Hermite), no
desde la salida del código. Coeficientes físicos distintos de 1 y elementos
oblicuos a propósito: un error dimensional o de transformación no se esconde
tras ``L = 1`` o ejes canónicos.
"""
import math

import numpy as np
import pytest

import solidum  # noqa: F401 — autodescubrimiento
from solidum.core.domain import Domain
from solidum.core.node import Node
from solidum.elements.frame import Frame2DEuler, Frame2DEulerCorot, Frame2DTimoshenko
from solidum.elements.frame3d import Frame3D
from solidum.materials.elastic import Elastic1D
from solidum.math.assembly import Assembler
from solidum.math.solvers import ModalSolver
from solidum.sections import Fiber, FiberSection, RectangularPatch

RHO_S = 7850.0          # acero [kg/m³]
RHO_C = 2400.0          # concreto [kg/m³]
E_S = 2.0e11
A_S, I_S = 1.5e-3, 3.2e-6          # sección de acero (camino con material)
IY_S, IZ_S = 1.1e-6, 3.2e-6

P1_2D, P2_2D = (0.3, -0.2), (2.1, 1.4)          # elemento 2D oblicuo
P1_3D, P2_3D = (0.2, -0.1, 0.4), (1.3, 1.9, 2.2)  # elemento 3D oblicuo
REF_3D = [0.3, -0.5, 0.8]                        # ref_vector genérico


def _acero():
    return Elastic1D(E=E_S, density=RHO_S)


def _nodos(p1, p2):
    return [Node(1, list(p1)), Node(2, list(p2))]


def _longitud(p1, p2):
    return math.dist(p1, p2)


# ---------------------------------------------------------------------------
# Fórmulas esperadas
# ---------------------------------------------------------------------------

def _hermitica(rhoA, L):
    """``ρAL/420·H`` en ``(v₁, θ₁, v₂, θ₂)`` con ``v' = θ`` (Cook et al.)."""
    return rhoA * L / 420.0 * np.array([
        [156.0,     22 * L,     54.0,     -13 * L],
        [22 * L,    4 * L**2,   13 * L,   -3 * L**2],
        [54.0,      13 * L,     156.0,    -22 * L],
        [-13 * L,   -3 * L**2,  -22 * L,  4 * L**2],
    ])


def _hermitica_xz(rhoA, L):
    """Plano ``xz`` de ``Frame3D`` en ``(w₁, θ_y1, w₂, θ_y2)``: con
    ``w' = −θ_y``, los términos que llevan exactamente un giro cambian de
    signo respecto de :func:`_hermitica`."""
    return rhoA * L / 420.0 * np.array([
        [156.0,     -22 * L,    54.0,     13 * L],
        [-22 * L,   4 * L**2,   -13 * L,  -3 * L**2],
        [54.0,      -13 * L,    156.0,    22 * L],
        [13 * L,    -3 * L**2,  22 * L,   4 * L**2],
    ])


def _barra(m_L):
    """Masa lineal de barra ``m_L/6·[[2,1],[1,2]]`` (``m_L`` = masa del tramo)."""
    return m_L / 6.0 * np.array([[2.0, 1.0], [1.0, 2.0]])


def _consistente_2d(rhoA, L, rhoI=None):
    M = np.zeros((6, 6))
    M[np.ix_([0, 3], [0, 3])] = _barra(rhoA * L)
    M[np.ix_([1, 2, 4, 5], [1, 2, 4, 5])] = _hermitica(rhoA, L)
    if rhoI is not None:                       # inercia rotatoria (Timoshenko)
        M[np.ix_([2, 5], [2, 5])] += _barra(rhoI * L)
    return M


def _concentrada_2d(rhoA, L, rhoI=None):
    m_t = rhoA * L / 2
    m_r = rhoA * L**3 / 78 if rhoI is None else rhoI * L / 2
    return np.diag([m_t, m_t, m_r, m_t, m_t, m_r])


def _consistente_3d(rhoA, rhoJp, L):
    M = np.zeros((12, 12))
    M[np.ix_([0, 6], [0, 6])] = _barra(rhoA * L)
    M[np.ix_([3, 9], [3, 9])] = _barra(rhoJp * L)
    M[np.ix_([1, 5, 7, 11], [1, 5, 7, 11])] = _hermitica(rhoA, L)
    M[np.ix_([2, 4, 8, 10], [2, 4, 8, 10])] = _hermitica_xz(rhoA, L)
    return M


def _concentrada_3d(rhoA, rhoJp, L):
    m_t, m_x, m_th = rhoA * L / 2, rhoJp * L / 2, rhoA * L**3 / 78
    return np.diag([m_t, m_t, m_t, m_x, m_th, m_th] * 2)


def _local(el, M_global):
    """``T·M_global·Tᵀ``: la matriz en ejes locales del elemento."""
    T = el._reference_transform() if isinstance(el, Frame2DEulerCorot) else el.T
    return T @ M_global @ T.T


def _cerca(M, esperada, rtol=1e-13):
    np.testing.assert_allclose(M, esperada, rtol=0, atol=rtol * np.abs(esperada).max())


# ---------------------------------------------------------------------------
# La fórmula hermítica frente a la integral ∫ρA NᵀN dx
# ---------------------------------------------------------------------------

def test_hermitica_es_la_integral_de_las_funciones_de_hermite():
    """``ρA·L·∫₀¹ NᵀN dξ`` con Gauss-Legendre de 4 puntos (exacta para el
    integrando de grado 6) reproduce ``ρAL/420·H``."""
    rhoA, L = RHO_S * A_S, 2.7
    xi, w = np.polynomial.legendre.leggauss(4)
    xi, w = 0.5 * (xi + 1), 0.5 * w
    M = np.zeros((4, 4))
    for x, wg in zip(xi, w):
        N = np.array([1 - 3 * x**2 + 2 * x**3, L * (x - 2 * x**2 + x**3),
                      3 * x**2 - 2 * x**3, L * (-x**2 + x**3)])
        M += wg * rhoA * L * np.outer(N, N)
    _cerca(M, _hermitica(rhoA, L))


def test_hrz_escala_la_diagonal_consistente():
    """HRZ: la diagonal consistente ``ρAL/420·(156, 4L², 156, 4L²)`` escalada
    para que la traslación sume ``ρAL`` da ``ρAL/2`` y ``ρAL³/78``."""
    rhoA, L = RHO_S * A_S, 2.7
    d = np.diag(_hermitica(rhoA, L))
    alpha = rhoA * L / (d[0] + d[2])
    np.testing.assert_allclose(alpha * d, [rhoA * L / 2, rhoA * L**3 / 78,
                                           rhoA * L / 2, rhoA * L**3 / 78], rtol=1e-14)


# ---------------------------------------------------------------------------
# Marcos 2D Euler-Bernoulli, elemento a elemento
# ---------------------------------------------------------------------------

MARCOS_EB_2D = [Frame2DEuler, Frame2DEulerCorot]


@pytest.mark.parametrize("cls", MARCOS_EB_2D, ids=lambda c: c.__name__)
class TestMarco2DEulerBernoulli:

    def _elemento(self, cls, I=I_S):
        return cls(1, _nodos(P1_2D, P2_2D), _acero(), A=A_S, I=I)

    def test_consistente(self, cls):
        el = self._elemento(cls)
        L = _longitud(P1_2D, P2_2D)
        _cerca(_local(el, el.compute_mass_matrix("consistent")),
               _consistente_2d(RHO_S * A_S, L))

    def test_concentrada_hrz(self, cls):
        el = self._elemento(cls)
        L = _longitud(P1_2D, P2_2D)
        M = el.compute_mass_matrix("lumped")
        # Diagonal también en globales (bloque traslacional m_t·I₂).
        _cerca(M, np.diag(np.diag(M)))
        _cerca(_local(el, M), _concentrada_2d(RHO_S * A_S, L))

    @pytest.mark.parametrize("lumping", ["consistent", "lumped"])
    def test_no_depende_de_I(self, cls, lumping):
        """Sin inercia rotatoria, ``I`` sólo entra en la rigidez."""
        M1 = self._elemento(cls, I=I_S).compute_mass_matrix(lumping)
        M2 = self._elemento(cls, I=50 * I_S).compute_mass_matrix(lumping)
        np.testing.assert_array_equal(M1, M2)

    def test_corrotacional_igual_al_lineal(self, cls):
        """En la configuración de referencia las dos masas coinciden."""
        ref = Frame2DEuler(1, _nodos(P1_2D, P2_2D), _acero(), A=A_S, I=I_S)
        el = self._elemento(cls)
        for lumping in ("consistent", "lumped"):
            _cerca(el.compute_mass_matrix(lumping), ref.compute_mass_matrix(lumping))


# ---------------------------------------------------------------------------
# Frame3D, elemento a elemento
# ---------------------------------------------------------------------------

class TestFrame3D:

    def _elemento(self, Iy=IY_S, Iz=IZ_S):
        return Frame3D(1, _nodos(P1_3D, P2_3D), _acero(), A=A_S, Iy=Iy, Iz=Iz,
                       J=0.7 * (Iy + Iz), nu=0.3, ref_vector=REF_3D)

    def test_consistente_en_los_dos_planos_y_torsion(self):
        el = self._elemento()
        L = _longitud(P1_3D, P2_3D)
        _cerca(_local(el, el.compute_mass_matrix("consistent")),
               _consistente_3d(RHO_S * A_S, RHO_S * (IY_S + IZ_S), L))

    def test_signo_del_plano_xz(self):
        """El bloque ``(uz_i, ry_i, uz_j, ry_j)`` es la hermítica con los
        acoplamientos traslación-giro de signo opuesto al plano ``xy``."""
        el = self._elemento()
        L = _longitud(P1_3D, P2_3D)
        M = _local(el, el.compute_mass_matrix("consistent"))
        xy = M[np.ix_([1, 5, 7, 11], [1, 5, 7, 11])]
        xz = M[np.ix_([2, 4, 8, 10], [2, 4, 8, 10])]
        S = np.diag([1.0, -1.0, 1.0, -1.0])
        _cerca(xz, S @ xy @ S)
        # Acoplamiento v₁–θ_z1 = +22L·ρAL/420; w₁–θ_y1 = −22L·ρAL/420.
        c = RHO_S * A_S * L / 420
        assert xy[0, 1] == pytest.approx(22 * L * c, rel=1e-12)
        assert xz[0, 1] == pytest.approx(-22 * L * c, rel=1e-12)

    def test_concentrada_hrz(self):
        el = self._elemento()
        L = _longitud(P1_3D, P2_3D)
        _cerca(_local(el, el.compute_mass_matrix("lumped")),
               _concentrada_3d(RHO_S * A_S, RHO_S * (IY_S + IZ_S), L))

    def test_concentrada_diagonal_con_eje_global(self):
        el = Frame3D(1, _nodos((0.0, 0.0, 0.0), (2.7, 0.0, 0.0)), _acero(), A=A_S,
                     Iy=IY_S, Iz=IZ_S, J=IY_S + IZ_S, nu=0.3)
        M = el.compute_mass_matrix("lumped")
        _cerca(M, _concentrada_3d(RHO_S * A_S, RHO_S * (IY_S + IZ_S), 2.7))

    @pytest.mark.parametrize("lumping", ["consistent", "lumped"])
    def test_Iy_Iz_solo_entran_en_la_torsion(self, lumping):
        """Cambiar ``I_y``, ``I_z`` sólo mueve el bloque torsional ``(rx_i, rx_j)``,
        a través de ``J_p = I_y + I_z``."""
        a = self._elemento(IY_S, IZ_S)
        b = self._elemento(3 * IY_S, 5 * IZ_S)
        Ma = _local(a, a.compute_mass_matrix(lumping))
        Mb = _local(b, b.compute_mass_matrix(lumping))
        resto = [i for i in range(12) if i not in (3, 9)]
        _cerca(Mb[np.ix_(resto, resto)], Ma[np.ix_(resto, resto)])
        L = _longitud(P1_3D, P2_3D)
        rhoJp_b = RHO_S * (3 * IY_S + 5 * IZ_S)
        esperada = _barra(rhoJp_b * L) if lumping == "consistent" else np.eye(2) * rhoJp_b * L / 2
        _cerca(Mb[np.ix_([3, 9], [3, 9])], esperada)


# ---------------------------------------------------------------------------
# Camino con sección
# ---------------------------------------------------------------------------

B, H = 0.3, 0.5
A_BARRA, Y_BARRA, Z_BARRA = 1.0e-3, -0.2, 0.08


def _seccion_mixta():
    """Concreto ``B × H`` (6 × 4 fibras) más una barra de acero excéntrica."""
    concreto = Elastic1D(E=3.0e10, density=RHO_C)
    return FiberSection(
        patches=[RectangularPatch(concreto, y=(-H / 2, H / 2), z=(-B / 2, B / 2), ny=6, nz=4)],
        fibers=[Fiber(Y_BARRA, Z_BARRA, A_BARRA, _acero())], GJ=1.0e7)


def _rho_a_mixta():
    return RHO_C * B * H + RHO_S * A_BARRA


def _rho_jp_mixta():
    """``Σρ_i A_i (y_i² + z_i²)`` de la discretización: las fibras de un
    parche ``n × m`` están en los centros de las celdas, cuyo momento de
    segundo orden es ``b h³/12·(1 − 1/n²)`` (regla del punto medio)."""
    Iz_c = B * H**3 / 12 * (1 - 1 / 6**2)
    Iy_c = H * B**3 / 12 * (1 - 1 / 4**2)
    return RHO_C * (Iy_c + Iz_c) + RHO_S * A_BARRA * (Y_BARRA**2 + Z_BARRA**2)


class TestCaminoConSeccion:

    @pytest.mark.parametrize("cls", MARCOS_EB_2D, ids=lambda c: c.__name__)
    def test_marcos_2d_con_rho_A_de_las_fibras(self, cls):
        el = cls(1, _nodos(P1_2D, P2_2D), section=_seccion_mixta())
        L = _longitud(P1_2D, P2_2D)
        rhoA = _rho_a_mixta()
        _cerca(_local(el, el.compute_mass_matrix("consistent")), _consistente_2d(rhoA, L))
        _cerca(_local(el, el.compute_mass_matrix("lumped")), _concentrada_2d(rhoA, L))

    def test_frame3d_con_rho_A_y_rho_Jp_de_las_fibras(self):
        el = Frame3D(1, _nodos(P1_3D, P2_3D), section=_seccion_mixta(), ref_vector=REF_3D)
        L = _longitud(P1_3D, P2_3D)
        rhoA, rhoJp = _rho_a_mixta(), _rho_jp_mixta()
        _cerca(_local(el, el.compute_mass_matrix("consistent")), _consistente_3d(rhoA, rhoJp, L))
        _cerca(_local(el, el.compute_mass_matrix("lumped")), _concentrada_3d(rhoA, rhoJp, L))

    @pytest.mark.parametrize("lumping", ["consistent", "lumped"])
    def test_igual_al_camino_con_material(self, lumping):
        """Densidad uniforme: la masa con sección es la del marco con material
        con ``A = ΣA_i``, ``I_z = ΣA_i y_i²``, ``I_y = ΣA_i z_i²``."""
        mat = _acero()
        sec = FiberSection(
            patches=[RectangularPatch(mat, y=(-H / 2, H / 2), z=(-B / 2, B / 2), ny=8, nz=5)],
            GJ=2.0e7)
        A = B * H
        Iz = B * H**3 / 12 * (1 - 1 / 8**2)
        Iy = H * B**3 / 12 * (1 - 1 / 5**2)
        pares = [
            (Frame2DEuler(1, _nodos(P1_2D, P2_2D), mat, A=A, I=Iz),
             Frame2DEuler(2, _nodos(P1_2D, P2_2D), section=sec)),
            (Frame2DEulerCorot(1, _nodos(P1_2D, P2_2D), mat, A=A, I=Iz),
             Frame2DEulerCorot(2, _nodos(P1_2D, P2_2D), section=sec)),
            (Frame3D(1, _nodos(P1_3D, P2_3D), mat, A=A, Iy=Iy, Iz=Iz, J=1e-3, nu=0.3,
                     ref_vector=REF_3D),
             Frame3D(2, _nodos(P1_3D, P2_3D), section=sec, ref_vector=REF_3D)),
        ]
        for con_material, con_seccion in pares:
            _cerca(con_seccion.compute_mass_matrix(lumping),
                   con_material.compute_mass_matrix(lumping), rtol=1e-12)


# ---------------------------------------------------------------------------
# Masa total y energía cinética de movimientos rígidos
# ---------------------------------------------------------------------------

def _velocidad_giro_2d(p1, p2, omega):
    """Giro rígido de velocidad angular ``ω`` alrededor del nodo 1:
    ``v = ω e_z × (x − x₁)`` y giro nodal ``ω`` en los dos nodos."""
    dx, dy = p2[0] - p1[0], p2[1] - p1[1]
    return np.array([0.0, 0.0, omega, -omega * dy, omega * dx, omega])


def _velocidad_giro_3d(p1, p2, omega):
    d = np.subtract(p2, p1)
    return np.concatenate([np.zeros(3), omega, np.cross(omega, d), omega])


def _elementos_2d():
    """(elemento, ρA, ρI o None) de los cuatro marcos 2D, con los dos caminos
    de los Euler-Bernoulli."""
    rhoA = RHO_S * A_S
    out = []
    for cls in MARCOS_EB_2D:
        out.append((cls(1, _nodos(P1_2D, P2_2D), _acero(), A=A_S, I=I_S), rhoA, None))
        out.append((cls(2, _nodos(P1_2D, P2_2D), section=_seccion_mixta()),
                    _rho_a_mixta(), None))
    out.append((Frame2DTimoshenko(3, _nodos(P1_2D, P2_2D), _acero(), A=A_S, I=I_S,
                                  As=5 / 6 * A_S, nu=0.3), rhoA, RHO_S * I_S))
    return out


@pytest.mark.parametrize("lumping", ["consistent", "lumped"])
def test_traslacion_rigida_2d(lumping):
    """``½u̇ᵀMu̇ = ρAL|v|²/2`` para una traslación ``v`` oblicua cualquiera."""
    L = _longitud(P1_2D, P2_2D)
    v = np.array([1.7, -0.9])
    for el, rhoA, _ in _elementos_2d():
        u = np.array([v[0], v[1], 0.0, v[0], v[1], 0.0])
        T_cin = 0.5 * u @ el.compute_mass_matrix(lumping) @ u
        assert T_cin == pytest.approx(rhoA * L * (v @ v) / 2, rel=1e-13), type(el).__name__


def test_giro_rigido_2d_consistente():
    """Giro rígido alrededor de un extremo: la hermítica reproduce exactamente
    ``½·ρA·L³/3·ω²`` (la interpolación cúbica contiene el campo lineal de
    velocidades). Sin inercia rotatoria en Euler-Bernoulli; Timoshenko suma
    ``½·ρI·L·ω²``."""
    L = _longitud(P1_2D, P2_2D)
    omega = 3.1
    u = _velocidad_giro_2d(P1_2D, P2_2D, omega)
    for el, rhoA, rhoI in _elementos_2d():
        esperada = 0.5 * rhoA * L**3 / 3 * omega**2
        if rhoI is not None:
            esperada += 0.5 * rhoI * L * omega**2
        T_cin = 0.5 * u @ el.compute_mass_matrix("consistent") @ u
        assert T_cin == pytest.approx(esperada, rel=1e-12), type(el).__name__


def _elementos_3d():
    return [
        (Frame3D(1, _nodos(P1_3D, P2_3D), _acero(), A=A_S, Iy=IY_S, Iz=IZ_S,
                 J=0.7 * (IY_S + IZ_S), nu=0.3, ref_vector=REF_3D),
         RHO_S * A_S, RHO_S * (IY_S + IZ_S)),
        (Frame3D(2, _nodos(P1_3D, P2_3D), section=_seccion_mixta(), ref_vector=REF_3D),
         _rho_a_mixta(), _rho_jp_mixta()),
    ]


@pytest.mark.parametrize("lumping", ["consistent", "lumped"])
def test_traslacion_rigida_3d(lumping):
    L = _longitud(P1_3D, P2_3D)
    v = np.array([1.7, -0.9, 0.6])
    u = np.concatenate([v, np.zeros(3), v, np.zeros(3)])
    for el, rhoA, _ in _elementos_3d():
        T_cin = 0.5 * u @ el.compute_mass_matrix(lumping) @ u
        assert T_cin == pytest.approx(rhoA * L * (v @ v) / 2, rel=1e-13)


def test_giro_rigido_3d_consistente():
    """Giro rígido de velocidad angular ``ω`` (oblicua) alrededor del nodo 1:
    ``½·ρA·L³/3·|ω_⊥|² + ½·ρJ_p·L·ω_∥²``, con ``ω_∥`` la componente sobre el
    eje de la barra (torsión) y ``ω_⊥`` la perpendicular (flexión en los dos
    planos). Un signo equivocado en el plano ``xz`` rompe la igualdad."""
    L = _longitud(P1_3D, P2_3D)
    e_x = np.subtract(P2_3D, P1_3D) / L
    for omega in (np.array([0.8, -1.3, 2.1]), np.array([0.0, 0.0, 1.0])):
        w_par = omega @ e_x
        w_perp2 = omega @ omega - w_par**2
        u = _velocidad_giro_3d(P1_3D, P2_3D, omega)
        for el, rhoA, rhoJp in _elementos_3d():
            esperada = 0.5 * rhoA * L**3 / 3 * w_perp2 + 0.5 * rhoJp * L * w_par**2
            T_cin = 0.5 * u @ el.compute_mass_matrix("consistent") @ u
            assert T_cin == pytest.approx(esperada, rel=1e-12)


def test_giro_rigido_por_plano_local_3d():
    """Giro alrededor de ``y_local`` y de ``z_local`` por separado: cada plano
    de flexión da ``½·ρA·L³/3·ω²``, sin inercia rotatoria."""
    el, rhoA, _ = _elementos_3d()[0]
    L = _longitud(P1_3D, P2_3D)
    for eje in (1, 2):                          # y_local, z_local
        omega = 2.3 * el.lam[eje]
        u = _velocidad_giro_3d(P1_3D, P2_3D, omega)
        T_cin = 0.5 * u @ el.compute_mass_matrix("consistent") @ u
        assert T_cin == pytest.approx(0.5 * rhoA * L**3 / 3 * 2.3**2, rel=1e-12)


# ---------------------------------------------------------------------------
# Frame2DTimoshenko conserva su masa con inercia rotatoria
# ---------------------------------------------------------------------------

class TestTimoshenkoSinCambio:

    def _elemento(self):
        return Frame2DTimoshenko(1, _nodos(P1_2D, P2_2D), _acero(), A=A_S, I=I_S,
                                 As=5 / 6 * A_S, nu=0.3)

    def test_consistente_con_inercia_rotatoria(self):
        el = self._elemento()
        L = _longitud(P1_2D, P2_2D)
        _cerca(_local(el, el.compute_mass_matrix("consistent")),
               _consistente_2d(RHO_S * A_S, L, rhoI=RHO_S * I_S))

    def test_concentrada_con_inercia_rotatoria(self):
        el = self._elemento()
        L = _longitud(P1_2D, P2_2D)
        _cerca(_local(el, el.compute_mass_matrix("lumped")),
               _concentrada_2d(RHO_S * A_S, L, rhoI=RHO_S * I_S))

    def test_diferencia_con_euler_bernoulli_es_la_inercia_rotatoria(self):
        """Mismo elemento: Timoshenko − Euler-Bernoulli = ``ρIL/6·[[2,1],[1,2]]``
        en ``(rz_i, rz_j)`` y nada más."""
        t = self._elemento()
        e = Frame2DEuler(2, _nodos(P1_2D, P2_2D), _acero(), A=A_S, I=I_S)
        L = _longitud(P1_2D, P2_2D)
        D = _local(t, t.compute_mass_matrix("consistent") - e.compute_mass_matrix("consistent"))
        esperada = np.zeros((6, 6))
        esperada[np.ix_([2, 5], [2, 5])] = _barra(RHO_S * I_S * L)
        np.testing.assert_allclose(D, esperada, rtol=0,
                                   atol=1e-12 * np.abs(t.compute_mass_matrix()).max())


# ---------------------------------------------------------------------------
# Modal: la masa consistente sin inercia rotatoria da cotas superiores de las
# frecuencias de Euler-Bernoulli
# ---------------------------------------------------------------------------

def test_frecuencias_biapoyada_cota_superior_de_euler_bernoulli():
    """Viga biapoyada robusta, ``r/L = 0.05``: con la masa consistente (la de
    Rayleigh-Ritz con las funciones de Hermite, conforme) las frecuencias
    discretas son cotas superiores de ``ω_n = (nπ/L)²√(EI/ρA)``.

    Escalas: el error de discretización de los cúbicos de Hermite con masa
    consistente es ``≈ (nπ/N)⁴/1440`` en ``ω`` (3.4e-5 en el modo 3 con
    ``N = 20``); la inercia rotatoria (viga de Rayleigh) bajaría ``ω_n`` en
    ``≈ ½(nπr/L)²`` = 1.2 %, 4.9 % y 11 %. La ventana
    ``[ω_n, ω_n(1 + 1e-3)]`` separa las dos cosas. ``ux = 0`` en todos los
    nodos aísla los modos de flexión."""
    L_tot, n_el = 4.0, 20
    r = 0.05 * L_tot
    A = 0.02
    I = A * r**2
    dom = Domain()
    nodos = [dom.add_node(i + 1, [i * L_tot / n_el, 0.0]) for i in range(n_el + 1)]
    mat = _acero()
    for i in range(n_el):
        dom.add_element(Frame2DEuler(i + 1, [nodos[i], nodos[i + 1]], mat, A=A, I=I))
    for nd in nodos:
        nd.fix_dof("ux", 0.0)
    nodos[0].fix_dof("uy", 0.0)
    nodos[-1].fix_dof("uy", 0.0)
    dom.generate_equation_numbers(verbose=False)
    w = ModalSolver(Assembler(dom), n_modes=3).solve().frequencies_rad
    w_eb = np.array([(n * math.pi / L_tot)**2 * math.sqrt(E_S * I / (RHO_S * A))
                     for n in (1, 2, 3)])
    ratio = w / w_eb
    assert np.all(ratio >= 1.0 - 1e-12), ratio
    assert np.all(ratio - 1.0 < 1e-3), ratio
