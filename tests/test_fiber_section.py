"""Tests unitarios de la sección por fibras y del camino con sección de los
marcos (ADR 0021, spec ``docs/specs/FiberSection.md``).

Cubren la implementación: construcción y validaciones de la sección y de sus
primitivas, cinemática y signos a nivel de sección, tangente frente a
diferencias finitas en régimen inelástico (sección y elementos), guardas de
los marcos, estado por punto de Gauss, masa y peso propio, simetría agregada,
YAML de extremo a extremo, exportación VTK y barras que desplazan el material
del parche (``displaces``, spec §3).

Los ocho casos de validación del contrato de la spec (``acceptance``) viven
aparte, en ``tests/validation/test_fiber_section_validation.py``.

Coeficientes físicos distintos de 1 a propósito (acero y concreto en
unidades SI): un error dimensional no se esconde tras ``E = A = 1``.
"""
import math
import re
import textwrap

import numpy as np
import pytest

import solidum  # noqa: F401 — autodescubrimiento
from solidum.core.domain import Domain
from solidum.core.material import Material
from solidum.core.node import Node
from solidum.core.section import Section
from solidum.elements.frame import Frame2DEuler, Frame2DEulerCorot, Frame2DTimoshenko
from solidum.elements.frame3d import Frame3D
from solidum.entry import run_yaml
from solidum.materials.damage_1d import IsotropicDamage1D
from solidum.materials.elastic import Elastic1D
from solidum.materials.elastic_2d import Elastic2D
from solidum.materials.plastic_1d import Elastoplastic1D
from solidum.math.assembly import Assembler
from solidum.math.solvers._shared import domain_is_symmetric
from solidum.registry import Registry, SectionRegistry
from solidum.sections import BarLayer, CircularPatch, Fiber, FiberSection, RectangularPatch
from solidum.utils.yaml_parser import YamlParser, YamlValidationError

E_S = 2.0e11          # acero [Pa]
E_C = 3.0e10          # concreto [Pa]
RHO_S = 7850.0
RHO_C = 2400.0
B, H = 0.3, 0.5       # sección rectangular [m]


def _acero(**kw):
    return Elastic1D(E=E_S, density=RHO_S, **kw)


def _rect(mat, ny=20, nz=1, y=(-H / 2, H / 2), z=(-B / 2, B / 2)):
    return RectangularPatch(mat, y=y, z=z, ny=ny, nz=nz)


def _nodos_2d(p1=(0.0, 0.0), p2=(3.0, 4.0)):
    return [Node(1, list(p1)), Node(2, list(p2))]


def _nodos_3d(p1=(0.0, 0.0, 0.0), p2=(1.0, 2.0, 3.0)):
    return [Node(1, list(p1)), Node(2, list(p2))]


def _fd_jacobian(fun, x, d):
    """Jacobiano por diferencias centradas de ``fun`` en ``x``."""
    x = np.asarray(x, dtype=float)
    f0 = fun(x)
    J = np.zeros((f0.size, x.size))
    for j in range(x.size):
        xp, xm = x.copy(), x.copy()
        xp[j] += d[j] if np.ndim(d) else d
        xm[j] -= d[j] if np.ndim(d) else d
        J[:, j] = (fun(xp) - fun(xm)) / (xp[j] - xm[j])
    return J


# ---------------------------------------------------------------------------
# Construcción y validaciones
# ---------------------------------------------------------------------------

class TestConstruccion:

    def test_seccion_solo_de_fibras_unilaterales_se_rechaza(self):
        """Decisión del usuario (2026-09-30): sin ninguna fibra bilateral la
        sección pierde toda la rigidez al comprimirse y K queda singular."""
        from solidum.materials.cable_1d import CableMaterial1D
        cable = CableMaterial1D(E=E_S, density=RHO_S)
        with pytest.raises(ValueError, match="unilateral"):
            FiberSection(patches=[_rect(cable, ny=4)])
        with pytest.raises(ValueError, match="unilateral"):
            FiberSection(fibers=[Fiber(0.1, 0.0, 1e-4, cable), Fiber(-0.1, 0.0, 1e-4, cable)])

    def test_fibras_unilaterales_mezcladas_se_admiten(self):
        """Tirantes o tendones dentro de una sección de concreto: uso legítimo."""
        from solidum.materials.cable_1d import CableMaterial1D
        cable = CableMaterial1D(E=E_S, density=RHO_S)
        concreto = Elastic1D(E=E_C, density=RHO_C)
        sec = FiberSection(patches=[_rect(concreto, ny=4)],
                           fibers=[Fiber(-0.2, 0.0, 1e-4, cable)])
        assert sec.n_fibers == 5
        # Con la sección entera comprimida, el concreto sigue aportando rigidez.
        _, k_s, _ = sec.response_2d(np.array([-1e-4, 0.0]), sec.initial_states())
        assert k_s[0, 0] > 0.0

    def test_registro_y_contrato(self):
        assert SectionRegistry.get("FiberSection") is FiberSection
        assert SectionRegistry.YAML_SECTION == "sections"
        assert Registry.for_spec_kind("section") is SectionRegistry
        assert issubclass(FiberSection, Section)
        # Definida después de la familia de materiales que referencia.
        secciones = [f.YAML_SECTION for f in Registry.families()]
        assert secciones.index("sections") > secciones.index("materials")

    def test_api_basica(self):
        concreto = Elastic1D(E=E_C, density=RHO_C)
        acero = _acero()
        sec = FiberSection(
            patches=[_rect(concreto, ny=4, nz=2),
                     BarLayer(acero, start=(-0.2, -0.1), end=(-0.2, 0.1), n=3, area=3.1e-4)],
            fibers=[Fiber(0.2, 0.0, 5.0e-4, acero)],
            GJ=1.0e6,
        )
        assert sec.n_fibers == 8 + 3 + 1
        for arr in (sec.y, sec.z, sec.area):
            assert isinstance(arr, np.ndarray) and arr.shape == (12,)
        assert sec.materials() == [concreto, acero]
        assert sec.GJ == 1.0e6
        assert sec.total_area() == pytest.approx(B * H + 3 * 3.1e-4 + 5.0e-4, rel=1e-14)
        # Orden de declaración: parches y después fibras sueltas.
        assert sec.y[-1] == 0.2 and sec.area[-1] == 5.0e-4

    def test_materiales_distintos_por_identidad(self):
        acero = _acero()
        sec = FiberSection(patches=[_rect(acero), _rect(acero, y=(0.3, 0.4))])
        assert sec.materials() == [acero]
        # La lista devuelta es una copia.
        sec.materials().append("x")
        assert sec.materials() == [acero]

    def test_initial_states(self):
        sec = FiberSection(patches=[_rect(_acero(), ny=5)])
        s1, s2 = sec.initial_states(), sec.initial_states()
        assert s1 == [None] * 5 and s1 is not s2

    @pytest.mark.parametrize("kwargs", [{}, {"patches": []}, {"patches": [], "fibers": []}])
    def test_sin_fibras(self, kwargs):
        with pytest.raises(ValueError, match="no tiene fibras"):
            FiberSection(**kwargs)

    def test_material_no_uniaxial(self):
        with pytest.raises(ValueError, match="STRAIN_DIM"):
            _rect(Elastic2D(E=E_S, nu=0.3))
        with pytest.raises(ValueError, match="compute_state"):
            Fiber(0.0, 0.0, 1e-4, object())

    @pytest.mark.parametrize("fabrica, patron", [
        (lambda m: RectangularPatch(m, y=(0.1, -0.1), z=(-0.1, 0.1), ny=2, nz=1), "y"),
        (lambda m: RectangularPatch(m, y=(-0.1, 0.1), z=(-0.1, 0.1), ny=0, nz=1), "ny"),
        (lambda m: RectangularPatch(m, y=(-0.1, 0.1), z=(-0.1, 0.1), ny=2.5, nz=1), "ny"),
        (lambda m: RectangularPatch(m, y=(-0.1, 0.1), z=(-0.1, 0.1), ny=2, nz=True), "nz"),
        (lambda m: RectangularPatch(m, y=(-0.1,), z=(-0.1, 0.1), ny=2, nz=1), "y"),
        (lambda m: CircularPatch(m, r_outer=0.1, r_inner=0.1, n_rad=2, n_circ=8), "r_outer > r_inner"),
        (lambda m: CircularPatch(m, r_outer=0.1, r_inner=-0.01, n_rad=2, n_circ=8), "r_inner"),
        (lambda m: CircularPatch(m, r_outer=-0.1, n_rad=2, n_circ=8), "r_outer"),
        (lambda m: CircularPatch(m, r_outer=0.1, n_rad=2, n_circ=0), "n_circ"),
        (lambda m: CircularPatch(m, r_outer=0.1, n_circ=8), "n_rad"),
        (lambda m: BarLayer(m, start=(0, 0), end=(0, 1), n=3, area=0.0), "area"),
        (lambda m: BarLayer(m, start=(0, 0), end=(0, 1), n=0, area=1e-4), "n="),
        (lambda m: BarLayer(m, start=(0, 0, 0), end=(0, 1), n=2, area=1e-4), "start"),
        (lambda m: Fiber(0.0, 0.0, -1e-4, m), "area"),
        (lambda m: Fiber(float("nan"), 0.0, 1e-4, m), "finito"),
    ])
    def test_primitiva_invalida(self, fabrica, patron):
        with pytest.raises(ValueError, match=patron):
            fabrica(_acero())

    @pytest.mark.parametrize("GJ", [0.0, -1.0, "abc"])
    def test_GJ_invalido(self, GJ):
        with pytest.raises(ValueError, match="GJ"):
            FiberSection(patches=[_rect(_acero())], GJ=GJ)

    def test_listas_mal_formadas(self):
        with pytest.raises(ValueError, match="debe ser una lista"):
            FiberSection(patches=_rect(_acero()))
        with pytest.raises(ValueError, match=r"patches\[0\] no es una primitiva"):
            FiberSection(patches=[object()])


# ---------------------------------------------------------------------------
# Primitivas de discretización
# ---------------------------------------------------------------------------

class TestPrimitivas:

    @pytest.mark.parametrize("n", [1, 2, 7, 50])
    def test_rectangulo_por_capas(self, n):
        """``EI`` de ``n`` capas = ``E·b·h³/12·(1 − 1/n²)`` (punto medio)."""
        sec = FiberSection(patches=[_rect(_acero(), ny=n, nz=1)])
        np.testing.assert_allclose(sec.area.sum(), B * H, rtol=1e-14)
        assert abs((sec.area * sec.y).sum()) < 1e-15 * B * H * H
        s, k, _ = sec.response_2d([0.0, 0.0], sec.initial_states())
        np.testing.assert_array_equal(s, [0.0, 0.0])
        EI = E_S * B * H**3 / 12.0 * (1.0 - 1.0 / n**2)
        np.testing.assert_allclose(k[0, 0], E_S * B * H, rtol=1e-13)
        np.testing.assert_allclose(k[1, 1], EI, rtol=1e-13)
        assert abs(k[0, 1]) < 1e-12 * E_S * B * H * H

    def test_rejilla_3d(self):
        ny, nz = 6, 4
        sec = FiberSection(patches=[_rect(_acero(), ny=ny, nz=nz)])
        assert sec.n_fibers == ny * nz
        np.testing.assert_allclose(sec.area, B * H / (ny * nz), rtol=1e-14)
        # Centros de celda (rejilla ordenada por y y después por z).
        np.testing.assert_allclose(sec.y.reshape(ny, nz)[:, 0],
                                   -H / 2 + (np.arange(ny) + 0.5) * H / ny, rtol=1e-12)
        np.testing.assert_allclose(sec.z.reshape(ny, nz)[0, :],
                                   -B / 2 + (np.arange(nz) + 0.5) * B / nz, rtol=1e-12)
        _, k, _ = sec.response_3d([0.0, 0.0, 0.0], None)
        EIz = E_S * B * H**3 / 12.0 * (1.0 - 1.0 / ny**2)
        EIy = E_S * H * B**3 / 12.0 * (1.0 - 1.0 / nz**2)
        np.testing.assert_allclose(np.diag(k), [E_S * B * H, EIy, EIz], rtol=1e-13)
        assert np.abs(k - np.diag(np.diag(k))).max() < 1e-12 * np.abs(k).max()

    def test_seccion_excentrica_acopla_N_y_M(self):
        """Coordenadas desde el eje del elemento, no desde el centroide."""
        y0, y1 = 0.1, 0.6
        sec = FiberSection(patches=[_rect(_acero(), ny=10, y=(y0, y1))])
        A = B * (y1 - y0)
        y_c = 0.5 * (y0 + y1)
        _, k, _ = sec.response_2d([0.0, 0.0], None)
        np.testing.assert_allclose(k[0, 1], -E_S * A * y_c, rtol=1e-12)
        eps0 = 1.0e-4
        s, _, _ = sec.response_2d([eps0, 0.0], None)
        np.testing.assert_allclose(s, [E_S * eps0 * A, -E_S * eps0 * A * y_c], rtol=1e-12)

    @pytest.mark.parametrize("n_rad, n_circ", [(1, 1), (1, 2), (3, 8), (5, 17)])
    def test_circulo_area_exacta(self, n_rad, n_circ):
        ro, ri = 0.2, 0.05
        sec = FiberSection(patches=[CircularPatch(_acero(), r_outer=ro, r_inner=ri,
                                                  n_rad=n_rad, n_circ=n_circ)])
        assert sec.n_fibers == n_rad * n_circ
        np.testing.assert_allclose(sec.area.sum(), math.pi * (ro**2 - ri**2), rtol=1e-14)

    def test_circulo_centroides_de_sector(self):
        r = 0.2
        # Medio disco (n_circ = 2): centroide a 4r/(3π) sobre la bisectriz;
        # el primer sector empieza en el ángulo 0 (+y) y gira hacia +z.
        sec = FiberSection(patches=[CircularPatch(_acero(), r_outer=r, n_rad=1, n_circ=2)])
        np.testing.assert_allclose([sec.y[0], sec.z[0]], [0.0, 4 * r / (3 * math.pi)],
                                   atol=1e-15)
        np.testing.assert_allclose([sec.y[1], sec.z[1]], [0.0, -4 * r / (3 * math.pi)],
                                   atol=1e-15)
        # Cuartos (n_circ = 4): primer sector centrado en 45°.
        sec = FiberSection(patches=[CircularPatch(_acero(), r_outer=r, center=(0.1, -0.3),
                                                  n_rad=1, n_circ=4)])
        alpha = math.pi / 4
        r_c = (2.0 / 3.0) * r * math.sin(alpha) / alpha
        np.testing.assert_allclose([sec.y[0], sec.z[0]],
                                   [0.1 + r_c * math.cos(alpha), -0.3 + r_c * math.sin(alpha)],
                                   rtol=1e-14)
        # Corona entera en un sector: centroide en el centro.
        sec = FiberSection(patches=[CircularPatch(_acero(), r_outer=r, center=(0.1, -0.3),
                                                  n_rad=1, n_circ=1)])
        np.testing.assert_allclose([sec.y[0], sec.z[0]], [0.1, -0.3], atol=1e-15)

    def test_circulo_hueco_inercia_converge_como_1_sobre_n2(self):
        ro, ri = 0.2, 0.15
        I_ex = math.pi * (ro**4 - ri**4) / 4
        errores = []
        for n in (2, 4, 8):
            sec = FiberSection(patches=[CircularPatch(_acero(), r_outer=ro, r_inner=ri,
                                                      n_rad=n, n_circ=8 * n)])
            I_y = (sec.area * sec.z**2).sum()
            I_z = (sec.area * sec.y**2).sum()
            np.testing.assert_allclose(I_y, I_z, rtol=1e-12)
            errores.append(abs(I_z - I_ex) / I_ex)
        for e1, e2 in zip(errores, errores[1:]):
            assert 3.5 < e1 / e2 < 4.5

    def test_capa_de_barras(self):
        sec = FiberSection(patches=[BarLayer(_acero(), start=(-0.2, -0.1), end=(0.2, 0.3),
                                             n=3, area=2.0e-4)])
        np.testing.assert_allclose(sec.y, [-0.2, 0.0, 0.2], atol=1e-15)
        np.testing.assert_allclose(sec.z, [-0.1, 0.1, 0.3], atol=1e-15)
        np.testing.assert_allclose(sec.area, 2.0e-4)
        una = FiberSection(patches=[BarLayer(_acero(), start=(0.1, 0.2), end=(0.5, 0.5),
                                             n=1, area=2.0e-4)])
        np.testing.assert_allclose([una.y[0], una.z[0]], [0.1, 0.2])

    def test_fibra_suelta(self):
        sec = FiberSection(fibers=[Fiber(y=0.1, z=-0.2, area=3.0e-4, material=_acero())])
        assert (sec.y[0], sec.z[0], sec.area[0]) == (0.1, -0.2, 3.0e-4)


# ---------------------------------------------------------------------------
# Cinemática, signos y fuerzas internas de la sección
# ---------------------------------------------------------------------------

class TestCinematica:

    def _sec(self):
        return FiberSection(patches=[_rect(_acero(), ny=4, nz=3)])

    def test_deformaciones_de_fibra(self):
        sec = self._sec()
        e2 = [1.0e-4, 2.0e-3]
        np.testing.assert_allclose(sec.strains_2d(e2), e2[0] - sec.y * e2[1], rtol=1e-14)
        e3 = [1.0e-4, -3.0e-3, 2.0e-3]
        np.testing.assert_allclose(sec.strains_3d(e3), e3[0] + sec.z * e3[1] - sec.y * e3[2],
                                   rtol=1e-14)

    def test_signos_reglas_5(self):
        """``M_z > 0`` tracciona ``y < 0``; ``M_y > 0`` tracciona ``z > 0``."""
        d = 0.2
        sec = FiberSection(fibers=[Fiber(-d, 0.0, 1e-3, _acero()), Fiber(d, 0.0, 1e-3, _acero()),
                                   Fiber(0.0, -d, 1e-3, _acero()), Fiber(0.0, d, 1e-3, _acero())])
        kz = 1.0e-3
        eps = sec.strains_2d([0.0, kz])
        assert eps[0] > 0 > eps[1]            # y < 0 en tracción
        s, _, _ = sec.response_2d([0.0, kz], None)
        assert s[1] > 0
        ky = 1.0e-3
        eps = sec.strains_3d([0.0, ky, 0.0])
        assert eps[3] > 0 > eps[2]            # z > 0 en tracción
        s, _, _ = sec.response_3d([0.0, ky, 0.0], None)
        assert s[1] > 0 and abs(s[2]) < 1e-9 * abs(s[1])

    def test_no_muta_los_estados_recibidos(self):
        mat = Elastoplastic1D(E=E_S, sigma_y=250e6, H=2e9)
        sec = FiberSection(patches=[_rect(mat, ny=6)])
        states = sec.initial_states()
        _, _, nuevos = sec.response_2d([0.0, 0.05], states)
        assert states == [None] * 6
        assert all(isinstance(st, dict) for st in nuevos)

    def test_dimensiones_incorrectas(self):
        sec = self._sec()
        with pytest.raises(ValueError, match="deformaciones generalizadas"):
            sec.response_2d([0.0, 0.0, 0.0], None)
        with pytest.raises(ValueError, match="estados de fibra"):
            sec.response_3d([0.0, 0.0, 0.0], [None] * 3)


# ---------------------------------------------------------------------------
# Tangente frente a diferencias finitas (régimen inelástico)
# ---------------------------------------------------------------------------

def _seccion_mixta():
    """Concreto con daño (ablandamiento) + acero elastoplástico con
    endurecimiento, excéntrica: ejercita E_t negativas, nulas y positivas."""
    concreto = IsotropicDamage1D(E=E_C, kappa_0=1.0e-4, alpha=0.99, density=RHO_C)
    acero = Elastoplastic1D(E=E_S, sigma_y=420e6, H=0.01 * E_S, density=RHO_S)
    return FiberSection(
        patches=[RectangularPatch(concreto, y=(-0.25, 0.27), z=(-0.15, 0.15), ny=13, nz=5),
                 BarLayer(acero, start=(-0.21, -0.11), end=(-0.21, 0.11), n=3, area=5.1e-4),
                 BarLayer(acero, start=(0.22, -0.11), end=(0.22, 0.11), n=2, area=3.1e-4)],
        GJ=2.0e7,
    )


class TestTangente:

    @pytest.mark.parametrize("dim", [2, 3])
    def test_k_s_frente_a_diferencias_finitas(self, dim):
        sec = _seccion_mixta()
        resp = sec.response_2d if dim == 2 else sec.response_3d
        # Historia: carga hasta e1 (daño y fluencia) y commit de los estados.
        e1 = np.array([2.0e-4, 1.1e-2]) if dim == 2 else np.array([2.0e-4, -4.3e-3, 1.1e-2])
        _, _, states = resp(e1, sec.initial_states())
        e2 = 1.07 * e1
        s, k, _ = resp(e2, states)
        J = _fd_jacobian(lambda e: resp(e, states)[0], e2, 1e-7 * np.abs(e2).max())
        np.testing.assert_allclose(k, J, rtol=0, atol=1e-6 * np.abs(k).max())
        # Régimen realmente inelástico: hay tangentes negativas y reducidas.
        _, _, st2 = resp(e2, states)
        danos = [st["damage"] for st in st2 if st and "damage" in st]
        alfas = [st["alpha"] for st in st2 if st and "alpha" in st]
        assert max(danos) > 0.5 and max(alfas) > 0.0
        np.testing.assert_array_equal(k, k.T)

    def _fd_elemento(self, el, u, d=1e-8):
        K, _ = el.compute_element_state(u)
        J = _fd_jacobian(lambda x: el.compute_element_state(x)[1], u, d)
        return K, J

    def test_frame2d_euler(self):
        mat = Elastoplastic1D(E=E_S, sigma_y=250e6, H=0.02 * E_S)
        sec = FiberSection(patches=[_rect(mat, ny=16)])
        el = Frame2DEuler(1, _nodos_2d(), section=sec, n_ip=3)
        u1 = np.array([0.0, 0.0, 0.0, 1.0e-3, -2.0e-3, 1.5e-2])
        el.compute_element_state(u1)
        el.commit_state()
        K, J = self._fd_elemento(el, 1.1 * u1)
        np.testing.assert_allclose(K, J, rtol=0, atol=1e-6 * np.abs(K).max())

    def test_frame2d_euler_corot_con_giro_grande(self):
        mat = Elastoplastic1D(E=E_S, sigma_y=250e6, H=0.02 * E_S)
        sec = FiberSection(patches=[_rect(mat, ny=16)])
        el = Frame2DEulerCorot(1, _nodos_2d((0.0, 0.0), (3.0, 1.0)), section=sec, n_ip=4)
        # Giro de sólido rígido grande + flexión que plastifica.
        u = np.array([0.0, 0.0, 0.9 + 0.02, -0.97, 1.87, 0.9 + 0.03])
        K, J = self._fd_elemento(el, u)
        np.testing.assert_allclose(K, J, rtol=0, atol=1e-6 * np.abs(K).max())

    def test_frame3d_con_torsion(self):
        sec = _seccion_mixta()
        el = Frame3D(1, _nodos_3d(), section=sec, n_ip=3)
        rng = np.random.default_rng(3)
        u = rng.normal(size=12) * 4.0e-3
        K, J = self._fd_elemento(el, u)
        np.testing.assert_allclose(K, J, rtol=0, atol=1e-6 * np.abs(K).max())


# ---------------------------------------------------------------------------
# Revisión adversaria 2026-09-29: términos cruzados y objetividad
# ---------------------------------------------------------------------------

def _k_hermite(L):
    """∫₀ᴸ dᵀd dx con d = [N₁'', N₂'', N₃'', N₄''] (spec §5): rigidez de
    flexión de Hermite para EI = 1 en el orden [v₁, θ₁, v₂, θ₂]."""
    return np.array([
        [12 / L**3, 6 / L**2, -12 / L**3, 6 / L**2],
        [6 / L**2, 4 / L, -6 / L**2, 2 / L],
        [-12 / L**3, -6 / L**2, 12 / L**3, -6 / L**2],
        [6 / L**2, 2 / L, -6 / L**2, 4 / L],
    ])


@pytest.mark.parametrize("n_ip", [2, 3])
def test_frame3d_seccion_excentrica_y_asimetrica_forma_cerrada(n_ip):
    """K elástica de ``Frame3D`` con una sección excéntrica y asimétrica
    (primeros momentos ``Q_y``, ``Q_z`` y producto ``I_yz`` no nulos) frente a
    la forma cerrada escrita a mano en los 12 DOF, incluidos los
    acoplamientos axial–flexión y flexión y–flexión z con sus signos.

    Con ``k_s = [[EA, Q_z, −Q_y], [Q_z, I_yy, −I_yz], [−Q_y, −I_yz, I_zz]]``
    (spec §4, momentos ponderados por ``E``), ``κ_z = d·[v₁, θz₁, v₂, θz₂]`` y
    ``κ_y = S·d·[w₁, θy₁, w₂, θy₂]`` con ``S = diag(−1, 1, −1, 1)`` (``w' = −θ_y``)
    y ``∫d dx = g = [0, −1, 0, 1]``::

        K_xy,xy = I_zz·k_H      K_xz,xz = I_yy·S k_H S     K_xz,xy = −I_yz·S k_H
        K_ax,xy = −Q_y·a⊗g      K_ax,xz = Q_z·a⊗(S g)      a = [−1/L, 1/L]

    Los tests de reproducción elástica usan secciones simétricas, donde estos
    bloques son nulos y no ven un error de signo en ellos.
    """
    c = Elastic1D(E=E_C)
    s = _acero()
    sec = FiberSection(
        patches=[RectangularPatch(c, y=(-0.18, 0.31), z=(-0.07, 0.22), ny=7, nz=5),
                 BarLayer(s, start=(-0.15, -0.05), end=(0.26, 0.19), n=4, area=4.9e-4)],
        fibers=[Fiber(0.29, -0.06, 8.0e-4, s)],
        GJ=3.7e7,
    )
    Ei = np.array([E_C] * 35 + [E_S] * 5)
    EA = (Ei * sec.area).sum()
    Qy, Qz = (Ei * sec.area * sec.y).sum(), (Ei * sec.area * sec.z).sum()
    Iyy, Izz = (Ei * sec.area * sec.z**2).sum(), (Ei * sec.area * sec.y**2).sum()
    Iyz = (Ei * sec.area * sec.y * sec.z).sum()
    assert min(abs(Qy) / math.sqrt(EA * Izz), abs(Qz) / math.sqrt(EA * Iyy),
               abs(Iyz) / math.sqrt(Iyy * Izz)) > 0.2     # acoplamientos relevantes

    L = 2.7
    p1 = np.array([0.5, -0.3, 0.2])
    el = Frame3D(1, _nodos_3d(tuple(p1), tuple(p1 + L * np.array([2.0, 3.0, 6.0]) / 7.0)),
                 section=sec, n_ip=n_ip)
    K_local = el.T @ el.compute_global_stiffness() @ el.T.T

    kH, S = _k_hermite(L), np.diag([-1.0, 1.0, -1.0, 1.0])
    g, a = np.array([0.0, -1.0, 0.0, 1.0]), np.array([-1.0 / L, 1.0 / L])
    ax, xy, xz, tx = [0, 6], [1, 5, 7, 11], [2, 4, 8, 10], [3, 9]
    Kc = np.zeros((12, 12))
    Kc[np.ix_(ax, ax)] = EA / L * np.array([[1.0, -1.0], [-1.0, 1.0]])
    Kc[np.ix_(tx, tx)] = sec.GJ / L * np.array([[1.0, -1.0], [-1.0, 1.0]])
    Kc[np.ix_(xy, xy)] = Izz * kH
    Kc[np.ix_(xz, xz)] = Iyy * S @ kH @ S
    Kc[np.ix_(xz, xy)] = -Iyz * S @ kH
    Kc[np.ix_(xy, xz)] = Kc[np.ix_(xz, xy)].T
    Kc[np.ix_(ax, xy)] = -Qy * np.outer(a, g)
    Kc[np.ix_(ax, xz)] = Qz * np.outer(a, S @ g)
    Kc[np.ix_(xy, ax)] = Kc[np.ix_(ax, xy)].T
    Kc[np.ix_(xz, ax)] = Kc[np.ix_(ax, xz)].T
    np.testing.assert_allclose(K_local, Kc, rtol=0, atol=1e-12 * np.abs(Kc).max())
    # Signos concretos que la forma cerrada fija (primer nodo):
    assert K_local[0, 4] == pytest.approx(Qz / L, rel=1e-12)       # u₁–θy₁
    assert K_local[0, 5] == pytest.approx(-Qy / L, rel=1e-12)      # u₁–θz₁
    assert K_local[4, 5] == pytest.approx(-4.0 * Iyz / L, rel=1e-12)
    assert K_local[2, 1] == pytest.approx(12.0 * Iyz / L**3, rel=1e-12)


@pytest.mark.parametrize("nombre", ["Frame2DEuler", "Frame2DEulerCorot", "Frame3D"])
def test_fibra_excentrica_traccionada_da_el_momento_de_la_excentricidad(nombre):
    """Blinda el signo de la fila de curvatura de ``B`` en cada marco por
    estática pura. Una fibra única en ``(e_y, e_z)`` con extensión axial
    ``δ`` (giros nulos) lleva la fuerza ``N = E A δ/L`` a una distancia
    ``r = e_y ŷ + e_z ẑ`` del eje; en el nodo j actúan ``+N x̂`` y el momento
    ``r × N x̂ = −e_y N ẑ + e_z N ŷ`` (en el i, los opuestos).

    Con una sección simétrica, invertir el signo de toda la fila de ``κ_z``
    (o de ``κ_y``) de ``B`` no cambia ``K`` ni ``F``: sólo qué fibras se
    traccionan. Antes de este test, esa mutación pasaba todos los tests del
    corotacional (revisión adversaria 2026-09-29).
    """
    e_y, e_z, A_f, L, delta = 0.12, -0.07, 1.0e-3, 2.5, 1.0e-4
    sec = FiberSection(fibers=[Fiber(e_y, e_z, A_f, _acero())], GJ=1.0e6)
    N = E_S * A_f * delta / L
    if nombre == "Frame3D":
        el = Frame3D(1, _nodos_3d((0.0, 0.0, 0.0), (L, 0.0, 0.0)), section=sec,
                     ref_vector=[0.0, 1.0, 0.0])
        u = np.zeros(12)
        u[6] = delta
        _, F = el.compute_element_state(u)
        esperado = np.zeros(12)
        esperado[[0, 6]] = [-N, N]
        esperado[[4, 10]] = [-e_z * N, e_z * N]       # My
        esperado[[5, 11]] = [e_y * N, -e_y * N]       # Mz
    else:
        cls = Frame2DEuler if nombre == "Frame2DEuler" else Frame2DEulerCorot
        el = cls(1, _nodos_2d((0.0, 0.0), (L, 0.0)), section=sec)
        _, F = el.compute_element_state(np.array([0.0, 0.0, 0.0, delta, 0.0, 0.0]))
        # El corotacional mide Δl = l − L₀ exacto: N difiere en O(δ²/L²).
        esperado = np.array([-N, 0.0, e_y * N, N, 0.0, -e_y * N])
    np.testing.assert_allclose(F, esperado, rtol=1e-9, atol=1e-9 * N)


@pytest.mark.parametrize("phi", [0.7, 2.8, -3.1, 5.0])
def test_corotacional_objetivo_con_historia_plastica(phi):
    """Un giro de sólido rígido ``Q(φ)`` (también ``|φ| > π``) aplicado a una
    configuración deformada con historia plástica comprometida rota las
    fuerzas y la tangente (``F(Qu) = Q·F(u)``, ``K(Qu) = Q K Qᵀ``) sin cambiar
    la plasticidad de las fibras; en particular, las fuerzas residuales de la
    historia rotan con el elemento y no aparecen fuerzas espurias."""
    mat = Elastoplastic1D(E=E_S, sigma_y=250e6, H=0.01 * E_S)
    sec = FiberSection(patches=[_rect(mat, ny=20)])
    p1, p2 = np.array([0.4, -0.2]), np.array([2.2, 2.2])
    el = Frame2DEulerCorot(1, _nodos_2d(tuple(p1), tuple(p2)), section=sec, n_ip=3)
    L0 = np.linalg.norm(p2 - p1)
    kap = 5 * 250e6 / E_S / (H / 2)
    u_def = np.array([0.0, 0.0, -kap * L0 / 2, 6e-4, 8e-4, 1.3 * kap * L0 / 2])
    el.compute_element_state(0.7 * u_def)
    el.commit_state()

    def rotado(u):
        Q = np.array([[math.cos(phi), -math.sin(phi)], [math.sin(phi), math.cos(phi)]])
        x1, x2 = p1 + u[0:2], p2 + u[3:5]
        x2n = x1 + Q @ (x2 - x1)
        return np.array([*(x1 - p1), u[2] + phi, *(x2n - p2), u[5] + phi]), Q

    def alfas():
        return np.array([d["alpha"] for st in el.state.vars_trial for d in st])

    K0, F0 = el.compute_element_state(u_def)
    a0 = alfas()
    u_r, Q = rotado(u_def)
    K1, F1 = el.compute_element_state(u_r)
    R = np.eye(6)
    R[0:2, 0:2] = R[3:5, 3:5] = Q
    np.testing.assert_allclose(F1, R @ F0, rtol=0, atol=1e-10 * np.abs(F0).max())
    np.testing.assert_allclose(K1, R @ K0 @ R.T, rtol=0, atol=1e-10 * np.abs(K0).max())
    np.testing.assert_allclose(alfas(), a0, rtol=1e-10, atol=0)
    assert a0.max() > 0.0
    # Residual de la historia a deformación corrotada nula: rota con el elemento.
    _, Fz = el.compute_element_state(np.zeros(6))
    _, Fzr = el.compute_element_state(rotado(np.zeros(6))[0])
    assert np.abs(Fz).max() > 0.0
    np.testing.assert_allclose(Fzr, R @ Fz, rtol=0, atol=1e-10 * np.abs(Fz).max())


# ---------------------------------------------------------------------------
# Marcos con sección: contrato del elemento y estado
# ---------------------------------------------------------------------------

_FABRICAS_CON_SECCION = {
    "Frame2DEuler": lambda sec, n_ip=3: Frame2DEuler(1, _nodos_2d(), section=sec, n_ip=n_ip),
    "Frame2DEulerCorot": lambda sec, n_ip=3: Frame2DEulerCorot(1, _nodos_2d(), section=sec,
                                                               n_ip=n_ip),
    "Frame3D": lambda sec, n_ip=3: Frame3D(1, _nodos_3d(), section=sec, n_ip=n_ip),
}


def _seccion_elastica(GJ=2.0e7):
    return FiberSection(patches=[_rect(_acero(), ny=10, nz=6)], GJ=GJ)


class TestMarcosConSeccion:

    def test_los_marcos_declaran_la_referencia_a_secciones(self):
        for cls in (Frame2DEuler, Frame2DEulerCorot, Frame3D):
            assert cls.REFERENCE_KWARGS["section"] is SectionRegistry
        assert "section" not in Frame2DTimoshenko.REFERENCE_KWARGS

    @pytest.mark.parametrize("nombre", sorted(_FABRICAS_CON_SECCION))
    @pytest.mark.parametrize("n_ip", [2, 3, 5])
    def test_estado_por_punto_de_gauss(self, nombre, n_ip):
        sec = _seccion_elastica()
        el = _FABRICAS_CON_SECCION[nombre](sec, n_ip)
        assert el.N_INTEGRATION_POINTS == n_ip and el.n_ip == n_ip
        assert type(el).N_INTEGRATION_POINTS == 1     # la clase no cambia
        assert len(el.state.vars) == n_ip
        assert all(v == [None] * sec.n_fibers for v in el.state.vars)
        dim = 3 if nombre == "Frame3D" else 2
        assert all(np.array_equal(s, np.zeros(dim)) for s in el.state.stresses)
        assert el.materials() == sec.materials()
        assert el.material is None

    @pytest.mark.parametrize("nombre", sorted(_FABRICAS_CON_SECCION))
    def test_n_ip_por_defecto_es_3(self, nombre):
        cls = {"Frame2DEuler": Frame2DEuler, "Frame2DEulerCorot": Frame2DEulerCorot,
               "Frame3D": Frame3D}[nombre]
        nodos = _nodos_3d() if nombre == "Frame3D" else _nodos_2d()
        el = cls(1, nodos, section=_seccion_elastica())
        assert el.N_INTEGRATION_POINTS == 3

    @pytest.mark.parametrize("nombre", sorted(_FABRICAS_CON_SECCION))
    def test_rigidez_simetrica_semidefinida_y_sin_fuerzas_por_traslacion(self, nombre):
        el = _FABRICAS_CON_SECCION[nombre](_seccion_elastica())
        K = el.compute_global_stiffness()
        escala = np.abs(K).max()
        # La local es simétrica exacta; la rotación a globales redondea.
        np.testing.assert_allclose(K, K.T, rtol=0, atol=1e-14 * escala)
        assert np.linalg.eigvalsh(K).min() > -1e-9 * escala
        n = len(el.DOF_NAMES)
        for eje, dof in enumerate(el.DOF_NAMES):
            if dof.startswith("u"):
                t = np.zeros(K.shape[0])
                t[eje::n] = 1.0
                assert np.abs(K @ t).max() < 1e-9 * escala

    @pytest.mark.parametrize("nombre", sorted(_FABRICAS_CON_SECCION))
    def test_evaluaciones_auxiliares_no_tocan_el_trial(self, nombre):
        mat = Elastoplastic1D(E=E_S, sigma_y=250e6, H=2e9)
        el = _FABRICAS_CON_SECCION[nombre](FiberSection(patches=[_rect(mat, ny=8, nz=2)],
                                                        GJ=2e7))
        n = len(el.nodes) * len(el.DOF_NAMES)
        u = np.linspace(-1.0, 1.0, n) * 5.0e-3
        el.compute_element_state(u)
        antes = el.state.snapshot_trial()
        dom = Domain()
        for nd in el.nodes:
            dom.nodes[nd.id] = nd
        dom.elements[el.id] = el
        dom.generate_equation_numbers(verbose=False)
        U = np.zeros(dom.total_dofs)
        el.compute_global_stiffness()
        el.internal_forces(U + 1e-3)
        el.compute_internal_forces(U + 1e-3)
        despues = el.state.snapshot_trial()
        assert antes[0] == despues[0]
        for a, b in zip(antes[1], despues[1]):
            np.testing.assert_array_equal(a, b)

    def test_commit_guarda_estados_y_fuerzas_de_seccion(self):
        mat = Elastoplastic1D(E=E_S, sigma_y=250e6, H=0.0)
        sec = FiberSection(patches=[_rect(mat, ny=20)])
        el = Frame2DEuler(1, _nodos_2d((0.0, 0.0), (2.0, 0.0)), section=sec, n_ip=3)
        # Curvatura uniforme (giros opuestos) muy por encima de la de fluencia.
        kappa = 4.0 * 2.0 * 250e6 / (E_S * H)
        L = 2.0
        u = np.array([0.0, 0.0, -kappa * L / 2, 0.0, 0.0, kappa * L / 2])
        el.compute_element_state(u)
        el.commit_state()
        M_p = 250e6 * B * H * H / 4
        for s in el.state.stresses:
            assert M_p * 0.95 < s[1] < M_p
        assert any(st["alpha"] > 0 for st in el.state.vars[0])
        # Descarga a desplazamiento nulo: queda momento residual (historia).
        _, F0 = el.compute_element_state(np.zeros(6))
        assert np.abs(F0).max() > 0.01 * M_p
        # Un elemento virgen no tiene fuerzas en u = 0.
        virgen = Frame2DEuler(2, _nodos_2d((0.0, 0.0), (2.0, 0.0)), section=sec, n_ip=3)
        assert np.abs(virgen.compute_element_state(np.zeros(6))[1]).max() == 0.0

    def test_diccionario_legado_con_seccion(self):
        el = Frame2DEuler(1, _nodos_2d(), section=_seccion_elastica(), n_ip=4)
        dom = Domain()
        for nd in el.nodes:
            dom.nodes[nd.id] = nd
        dom.elements[1] = el
        dom.generate_equation_numbers(verbose=False)
        U = np.zeros(dom.total_dofs)
        U[el.nodes[1].dofs["rz"]] = 1.0e-3
        out = el.compute_internal_forces(U)
        assert out["section_forces"].shape == (4, 2)
        assert out["section_deformations"].shape == (4, 2)
        assert "stress" not in out

    def test_torsion_elastica_desacoplada(self):
        GJ = 3.3e7
        el = Frame3D(1, _nodos_3d((0.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
                     section=_seccion_elastica(GJ=GJ))
        u = np.zeros(12)
        u[9] = 1.0e-3        # θx en el nodo j
        K, F = el.compute_element_state(u)
        np.testing.assert_allclose(F[[3, 9]], [-GJ * 1e-3 / 2.0, GJ * 1e-3 / 2.0], rtol=1e-14)
        assert np.abs(np.delete(F, [3, 9])).max() == 0.0
        np.testing.assert_allclose(K[3, 3], GJ / 2.0, rtol=1e-14)


# ---------------------------------------------------------------------------
# Guardas de los marcos
# ---------------------------------------------------------------------------

class TestGuardas:

    @pytest.mark.parametrize("fabrica", [
        lambda m: Frame2DEuler(1, _nodos_2d(), m, A=0.15, I=3e-3),
        lambda m: Frame2DEulerCorot(1, _nodos_2d(), m, A=0.15, I=3e-3),
        lambda m: Frame3D(1, _nodos_3d(), m, A=0.15, Iy=1e-3, Iz=3e-3, J=2e-3, nu=0.3),
    ])
    @pytest.mark.parametrize("material", [
        Elastoplastic1D(E=E_S, sigma_y=250e6),
        IsotropicDamage1D(E=E_C, kappa_0=1e-4, alpha=0.99),
    ])
    def test_material_inelastico_remite_a_la_seccion(self, fabrica, material):
        with pytest.raises(ValueError) as exc:
            fabrica(material)
        mensaje = str(exc.value)
        assert "IS_LINEAR_ELASTIC" in mensaje
        assert "FiberSection" in mensaje and "'section'" in mensaje
        assert type(material).__name__ in mensaje

    def test_material_sin_declaracion_se_rechaza(self):
        """Un material duck-typed sin ``IS_LINEAR_ELASTIC`` no pasa: el
        default es ``False`` (hay que declararlo)."""
        class Lineal:
            STRAIN_DIM = 1

            def compute_state(self, strain, state=None):
                return E_S * strain, E_S, state

        with pytest.raises(ValueError, match="IS_LINEAR_ELASTIC"):
            Frame2DEuler(1, _nodos_2d(), Lineal(), A=0.15, I=3e-3)
        Lineal.IS_LINEAR_ELASTIC = True
        Frame2DEuler(1, _nodos_2d(), Lineal(), A=0.15, I=3e-3)

    def test_timoshenko_solo_elastico(self):
        Frame2DTimoshenko(1, _nodos_2d(), _acero(), A=0.15, I=3e-3, As=0.125)
        with pytest.raises(ValueError, match="todavía no admite sección por fibras"):
            Frame2DTimoshenko(1, _nodos_2d(), Elastoplastic1D(E=E_S, sigma_y=250e6),
                              A=0.15, I=3e-3, As=0.125, nu=0.3)

    def test_dimension_antes_que_linealidad(self):
        """Un material 2D se rechaza por su dimensión (validación de la base)."""
        with pytest.raises(ValueError, match="STRAIN_DIM"):
            Frame2DEuler(1, _nodos_2d(), Elastic2D(E=E_S, nu=0.3), A=0.15, I=3e-3)

    def test_exactamente_un_camino(self):
        sec = _seccion_elastica()
        with pytest.raises(ValueError, match="no ambos"):
            Frame2DEuler(1, _nodos_2d(), _acero(), section=sec)
        with pytest.raises(ValueError, match="falta la constitutiva"):
            Frame2DEuler(1, _nodos_2d())
        with pytest.raises(ValueError, match=r"no se declaran \['A'\]"):
            Frame2DEulerCorot(1, _nodos_2d(), A=0.15, section=sec)
        with pytest.raises(ValueError, match=r"no se declaran \['nu'\]"):
            Frame3D(1, _nodos_3d(), nu=0.3, section=sec)
        with pytest.raises(ValueError, match=r"obligatorias las propiedades geométricas \['I'\]"):
            Frame2DEuler(1, _nodos_2d(), _acero(), A=0.15)
        with pytest.raises(ValueError, match=r"\['Iy', 'J'\]"):
            Frame3D(1, _nodos_3d(), _acero(), A=0.15, Iz=3e-3, nu=0.3)

    @pytest.mark.parametrize("n_ip", [1, 0, 2.0, True, "3"])
    def test_n_ip_invalido(self, n_ip):
        with pytest.raises(ValueError, match="n_ip"):
            Frame2DEuler(1, _nodos_2d(), section=_seccion_elastica(), n_ip=n_ip)

    def test_n_ip_solo_con_seccion(self):
        with pytest.raises(ValueError, match="n_ip=5 sólo aplica con 'section'"):
            Frame2DEuler(1, _nodos_2d(), _acero(), A=0.15, I=3e-3, n_ip=5)

    def test_section_debe_ser_una_seccion(self):
        with pytest.raises(ValueError, match="sección de marco"):
            Frame2DEuler(1, _nodos_2d(), section=_rect(_acero()))

    def test_frame3d_exige_GJ(self):
        sec = FiberSection(patches=[_rect(_acero(), ny=4, nz=4)])
        with pytest.raises(ValueError, match="no declara GJ"):
            Frame3D(1, _nodos_3d(), section=sec)
        # Los marcos planos no la necesitan.
        Frame2DEuler(1, _nodos_2d(), section=sec)


# ---------------------------------------------------------------------------
# Masa y peso propio
# ---------------------------------------------------------------------------

def _dominio(*elementos):
    dom = Domain()
    for el in elementos:
        for nd in el.nodes:
            dom.nodes.setdefault(nd.id, nd)
        dom.elements[el.id] = el
    dom.generate_equation_numbers(verbose=False)
    return dom


class TestMasa:

    def test_propiedades_masicas_de_seccion_mixta(self):
        concreto = Elastic1D(E=E_C, density=RHO_C)
        acero = _acero()
        a_b, y_b = 5.0e-4, -0.2
        sec = FiberSection(patches=[_rect(concreto, ny=10, nz=4),
                                    BarLayer(acero, start=(y_b, -0.1), end=(y_b, 0.1),
                                             n=2, area=a_b)])
        mp = sec.mass_properties()
        np.testing.assert_allclose(mp["rhoA"], RHO_C * B * H + RHO_S * 2 * a_b, rtol=1e-13)
        np.testing.assert_allclose(
            mp["rhoIz"], RHO_C * B * H**3 / 12 * (1 - 1 / 100) + RHO_S * 2 * a_b * y_b**2,
            rtol=1e-13)
        np.testing.assert_allclose(
            mp["rhoIy"], RHO_C * H * B**3 / 12 * (1 - 1 / 16) + RHO_S * 2 * a_b * 0.1**2,
            rtol=1e-13)

    def test_densidad_ausente(self):
        sec = FiberSection(patches=[_rect(Elastic1D(E=E_C)), _rect(_acero(), y=(0.3, 0.4))])
        with pytest.raises(ValueError, match=r"\['Elastic1D'\]"):
            sec.mass_properties()
        el = Frame2DEuler(1, _nodos_2d(), section=sec)
        with pytest.raises(ValueError, match="Materiales sin `density`"):
            Assembler(_dominio(el)).assemble_mass_matrix()

    @pytest.mark.parametrize("lumping", ["consistent", "lumped"])
    def test_misma_masa_que_el_camino_con_material(self, lumping):
        """Con densidad uniforme, la masa con sección es la del marco elástico
        con ``A = ΣA_i`` e ``I = ΣA_i y_i²`` (spec §6)."""
        mat = _acero()
        sec = FiberSection(patches=[_rect(mat, ny=10, nz=6)], GJ=2e7)
        A = sec.area.sum()
        Iz = (sec.area * sec.y**2).sum()
        Iy = (sec.area * sec.z**2).sum()
        pares = [
            (Frame2DEuler(1, _nodos_2d(), mat, A=A, I=Iz),
             Frame2DEuler(2, _nodos_2d(), section=sec)),
            (Frame2DEulerCorot(1, _nodos_2d(), mat, A=A, I=Iz),
             Frame2DEulerCorot(2, _nodos_2d(), section=sec)),
            (Frame3D(1, _nodos_3d(), mat, A=A, Iy=Iy, Iz=Iz, J=1e-3, nu=0.3),
             Frame3D(2, _nodos_3d(), section=sec)),
        ]
        for con_material, con_seccion in pares:
            M1 = con_material.compute_mass_matrix(lumping)
            M2 = con_seccion.compute_mass_matrix(lumping)
            np.testing.assert_allclose(M2, M1, rtol=0, atol=1e-12 * np.abs(M1).max())

    def test_masa_total(self):
        concreto = Elastic1D(E=E_C, density=RHO_C)
        sec = FiberSection(patches=[_rect(concreto, ny=5)],
                           fibers=[Fiber(-0.2, 0.0, 1e-3, _acero())], GJ=1e7)
        rhoA = RHO_C * B * H + RHO_S * 1e-3
        for el in (Frame2DEuler(1, _nodos_2d(), section=sec),
                   Frame3D(2, _nodos_3d(), section=sec)):
            L = el.L0
            for lumping in ("consistent", "lumped"):
                M = el.compute_mass_matrix(lumping)
                n = len(el.DOF_NAMES)
                # Traslación rígida en x global: masa total ρA·L.
                t = np.zeros(M.shape[0])
                t[0::n] = 1.0
                np.testing.assert_allclose(t @ M @ t, rhoA * L, rtol=1e-12)


class TestPesoPropio:

    def _seccion(self):
        concreto = Elastic1D(E=E_C, density=RHO_C)
        return FiberSection(patches=[_rect(concreto, ny=6),
                                     BarLayer(_acero(), start=(-0.2, -0.1), end=(-0.2, 0.1),
                                              n=3, area=4e-4)], GJ=1e7)

    def test_peso_propio_2d_usa_rhoA(self):
        sec = self._seccion()
        rhoA = RHO_C * B * H + RHO_S * 3 * 4e-4
        L = 4.0
        for cls in (Frame2DEuler, Frame2DEulerCorot):
            el = cls(1, [Node(1, [0.0, 0.0]), Node(2, [L, 0.0])], section=sec)
            asm = Assembler(_dominio(el))
            F = asm.assemble_self_weight([0.0, -9.81])
            q = -9.81 * rhoA
            esperado = np.array([0.0, q * L / 2, q * L**2 / 12, 0.0, q * L / 2, -q * L**2 / 12])
            np.testing.assert_allclose(F, esperado, rtol=1e-13)
            np.testing.assert_allclose(asm.equivalent_loads[1], esperado, rtol=1e-13)

    def test_peso_propio_3d_como_el_camino_con_material_equivalente(self):
        sec = self._seccion()
        rhoA = RHO_C * B * H + RHO_S * 3 * 4e-4
        A = sec.total_area()
        mat_eq = Elastic1D(E=E_S, density=rhoA / A)
        con_seccion = Frame3D(1, _nodos_3d(), section=sec)
        con_material = Frame3D(1, _nodos_3d(), mat_eq, A=A, Iy=1e-3, Iz=1e-3, J=1e-3, nu=0.3)
        g = [0.3, -1.2, -9.81]
        F1 = Assembler(_dominio(con_seccion)).assemble_self_weight(g)
        F2 = Assembler(_dominio(con_material)).assemble_self_weight(g)
        np.testing.assert_allclose(F1, F2, rtol=1e-13, atol=1e-13 * np.abs(F2).max())

    def test_fuerza_de_cuerpo_usa_el_area_bruta(self):
        sec = self._seccion()
        A = sec.total_area()
        b = np.array([0.0, -5.0e3])
        el = Frame2DEuler(1, _nodos_2d(), section=sec)
        ref = Frame2DEuler(1, _nodos_2d(), _acero(), A=A, I=1e-3)
        np.testing.assert_allclose(el.compute_body_load(b), ref.compute_body_load(b), rtol=1e-14)

    def test_densidad_ausente_en_una_fibra(self):
        sec = FiberSection(patches=[_rect(_acero())],
                           fibers=[Fiber(0.0, 0.0, 1e-4, Elastoplastic1D(E=E_S, sigma_y=250e6))])
        el = Frame2DEuler(1, _nodos_2d(), section=sec)
        with pytest.raises(ValueError, match=r"\['Elastoplastic1D'\]"):
            Assembler(_dominio(el)).assemble_self_weight([0.0, -9.81])


# ---------------------------------------------------------------------------
# Simetría agregada (ADR 0003) sobre los materiales de las fibras
# ---------------------------------------------------------------------------

class _NoSimetrico1D(Material):
    STRAIN_DIM = 1
    STATE_SCHEMA = {}
    IS_SYMMETRIC = False

    def __init__(self, E=E_S, density=None):
        self.E = E
        self.density = density

    def compute_state(self, strain, state_vars=None):
        return self.E * strain, self.E, state_vars


def test_simetria_agregada_ve_los_materiales_de_las_fibras():
    simetrica = FiberSection(patches=[_rect(_acero())])
    mixta = FiberSection(patches=[_rect(_acero())], fibers=[Fiber(0.0, 0.0, 1e-4, _NoSimetrico1D())])
    assert domain_is_symmetric(_dominio(Frame2DEuler(1, _nodos_2d(), section=simetrica)))
    assert not domain_is_symmetric(_dominio(Frame2DEuler(1, _nodos_2d(), section=mixta)))


def test_element_materials_por_defecto():
    mat = _acero()
    assert Frame2DEuler(1, _nodos_2d(), mat, A=0.1, I=1e-3).materials() == [mat]
    assert Frame2DTimoshenko(1, _nodos_2d(), mat, A=0.1, I=1e-3, As=0.08).materials() == [mat]


# ---------------------------------------------------------------------------
# YAML de extremo a extremo
# ---------------------------------------------------------------------------

def _escribir(tmp_path, texto):
    ruta = tmp_path / "modelo.yaml"
    ruta.write_text(textwrap.dedent(texto), encoding="utf-8")
    return ruta


def _resolver_yaml(ruta):
    """``(dominio, SolveResult)`` del modelo YAML: ``run_yaml`` y el dominio
    que arma, para leer los DOFs por nodo."""
    captura = {}
    original = YamlParser.parse

    def parse(self):
        captura["dominio"] = original(self)
        return captura["dominio"]

    YamlParser.parse = parse
    try:
        res = run_yaml(ruta)
    finally:
        YamlParser.parse = original
    return captura["dominio"], res


_VOLADIZO_2D = """\
nodes:
  - {id: 1, coords: [0.0, 0.0]}
  - {id: 2, coords: [1.0, 0.0]}
  - {id: 3, coords: [2.0, 0.0]}
  - {id: 4, coords: [3.0, 0.0]}
  - {id: 5, coords: [4.0, 0.0]}
materials:
  - {id: 1, type: Elastic1D, E: 2.0e11, density: 7850.0}
sections:
  - id: 7
    type: FiberSection
    patches:
      - {shape: rectangle, material: 1, y: [-0.25, 0.25], z: [-0.15, 0.15], ny: 50, nz: 1}
elements:
  - {id: 1, type: Frame2DEuler, nodes: [1, 2], section: 7, n_ip: 2}
  - {id: 2, type: Frame2DEuler, nodes: [2, 3], section: 7}
  - {id: 3, type: Frame2DEuler, nodes: [3, 4], section: 7, n_ip: 4}
  - {id: 4, type: Frame2DEuler, nodes: [4, 5], section: 7}
boundary_conditions_by_node:
  - {node_id: 1, ux: 0.0, uy: 0.0, rz: 0.0}
point_loads:
  - {node_id: 5, uy: -1.0e5}
solver:
  type: LinearSolver
"""


class TestYaml:

    def test_voladizo_elastico(self, tmp_path):
        """Flecha de punta ``PL³/(3EI)`` con ``EI`` de 50 capas: los
        elementos cúbicos son nodalmente exactos."""
        dom, res = _resolver_yaml(_escribir(tmp_path, _VOLADIZO_2D))
        P, L = -1.0e5, 4.0
        EI = E_S * B * H**3 / 12 * (1 - 1 / 50**2)
        uy = res.U[dom.nodes[5].dofs["uy"]]
        np.testing.assert_allclose(uy, P * L**3 / (3 * EI), rtol=1e-10)
        # Fuerzas internas de extremo: M(0) = P·L (hogging, M < 0).
        np.testing.assert_allclose(res.element_forces[1].components["M"][0], P * L, rtol=1e-10)
        np.testing.assert_allclose(res.element_forces[4].components["M"][1], 0.0,
                                   atol=1e-8 * abs(P * L))

    def test_seccion_con_todas_las_formas(self, tmp_path):
        texto = _VOLADIZO_2D.replace(
            "      - {shape: rectangle, material: 1, y: [-0.25, 0.25], z: [-0.15, 0.15], ny: 50, nz: 1}\n",
            "      - {shape: rectangle, material: 1, y: [-0.25, 0.25], z: [-0.15, 0.15], ny: 10, nz: 2}\n"
            "      - {shape: circle, material: 1, r_outer: 0.2, r_inner: 0.1, center: [0.0, 0.5], n_rad: 2, n_circ: 8}\n"
            "      - {shape: layer, material: 2, start: [-0.2, -0.1], end: [-0.2, 0.1], n: 3, area: 3.1e-4}\n"
            "    fibers:\n"
            "      - {y: 0.1, z: 0.0, area: 1.0e-4, material: 2}\n"
            "    GJ: 1.0e6\n",
        ).replace(
            "  - {id: 1, type: Elastic1D, E: 2.0e11, density: 7850.0}\n",
            "  - {id: 1, type: Elastic1D, E: 3.0e10, density: 2400.0}\n"
            "  - {id: 2, type: Elastoplastic1D, E: 2.0e11, sigma_y: 4.2e8, density: 7850.0}\n",
        )
        parser = YamlParser(str(_escribir(tmp_path, texto)))
        dom = parser.parse()
        sec = parser.sections[7]
        assert isinstance(sec, FiberSection)
        assert sec.n_fibers == 20 + 16 + 3 + 1
        assert sec.GJ == 1.0e6
        assert [type(m).__name__ for m in sec.materials()] == ["Elastic1D", "Elastoplastic1D"]
        assert all(el.section is sec for el in dom.elements.values())
        assert [el.n_ip for el in dom.elements.values()] == [2, 3, 4, 3]

    def test_no_lineal_con_plasticidad(self, tmp_path):
        """Voladizo elastoplástico con ``M(0) = 0.9·M_p > M_y``: converge y la
        flecha supera a la elástica."""
        M_p = 2.5e8 * B * H * H / 4
        P = -0.9 * M_p / 4.0
        texto = _VOLADIZO_2D.replace(
            "  - {id: 1, type: Elastic1D, E: 2.0e11, density: 7850.0}\n",
            "  - {id: 1, type: Elastoplastic1D, E: 2.0e11, sigma_y: 2.5e8, H: 2.0e9}\n",
        ).replace("ny: 50", "ny: 20").replace(
            "  - {node_id: 5, uy: -1.0e5}", f"  - {{node_id: 5, uy: {P!r}}}"
        ).replace("  type: LinearSolver", "  type: NonlinearSolver\n  num_steps: 5")
        dom, res = _resolver_yaml(_escribir(tmp_path, texto))
        assert res.converged
        EI = E_S * B * H**3 / 12 * (1 - 1 / 20**2)
        elastica = P * 4.0**3 / (3 * EI)
        # Más flexible que la elástica.
        assert res.U[dom.nodes[5].dofs["uy"]] < 1.05 * elastica
        # M(0) = P·L por equilibrio, también con plasticidad.
        np.testing.assert_allclose(res.element_forces[1].components["M"][0], P * 4.0,
                                   rtol=1e-6)

    def test_frame3d_torsion_con_GJ_del_yaml(self, tmp_path):
        texto = """\
        nodes:
          - {id: 1, coords: [0.0, 0.0, 0.0]}
          - {id: 2, coords: [2.5, 0.0, 0.0]}
        materials:
          - {id: 1, type: Elastic1D, E: 2.0e11}
        sections:
          - id: 1
            type: FiberSection
            GJ: 4.0e6
            patches:
              - {shape: circle, material: 1, r_outer: 0.1, n_rad: 4, n_circ: 16}
        elements:
          - {id: 1, type: Frame3D, nodes: [1, 2], section: 1}
        boundary_conditions_by_node:
          - {node_id: 1, ux: 0.0, uy: 0.0, uz: 0.0, rx: 0.0, ry: 0.0, rz: 0.0}
        point_loads:
          - {node_id: 2, rx: 1.0e3}
        solver:
          type: LinearSolver
        """
        dom, res = _resolver_yaml(_escribir(tmp_path, texto))
        rx = res.U[dom.nodes[2].dofs["rx"]]
        np.testing.assert_allclose(rx, 1.0e3 * 2.5 / 4.0e6, rtol=1e-12)
        np.testing.assert_allclose(res.element_forces[1].components["T"], [1.0e3, 1.0e3],
                                   rtol=1e-12)

    def test_material_inexistente_en_la_seccion(self, tmp_path):
        texto = _VOLADIZO_2D.replace("{shape: rectangle, material: 1,",
                                     "{shape: rectangle, material: 9,")
        with pytest.raises(YamlValidationError) as exc:
            YamlParser(str(_escribir(tmp_path, texto))).parse()
        mensaje = str(exc.value)
        assert "sections[0] (id=7)" in mensaje
        assert "patches[0].material" in mensaje
        assert "referencia a material inexistente (9)" in mensaje

    def test_seccion_inexistente_en_el_elemento(self, tmp_path):
        texto = _VOLADIZO_2D.replace("nodes: [2, 3], section: 7", "nodes: [2, 3], section: 8")
        with pytest.raises(YamlValidationError,
                           match=r"referencia a sección inexistente \(section=8\)"):
            YamlParser(str(_escribir(tmp_path, texto))).parse()

    @pytest.mark.parametrize("cambio, patron", [
        (("shape: rectangle", "shape: hexagon"), "shape='hexagon' desconocida"),
        (("ny: 50, nz: 1}", "ny: 50, nz: 1, nx: 3}"), r"parámetros \['nx'\] no aceptados"),
        (("ny: 50", "ny: 0"), r"sections\[0\] \(id=7\): patches\[0\]: RectangularPatch: ny=0"),
        (("{shape: rectangle, material: 1, ", "{shape: rectangle, "), "falta el campo obligatorio 'material'"),
    ])
    def test_errores_de_la_seccion_con_contexto(self, tmp_path, cambio, patron):
        texto = _VOLADIZO_2D.replace(*cambio)
        with pytest.raises(YamlValidationError, match=patron):
            YamlParser(str(_escribir(tmp_path, texto))).parse()

    @pytest.mark.parametrize("familia, bloque", [
        ("sección", "sections:\n"
                    "  - id: 7\n"
                    "    type: FiberSection\n"
                    "    fibers:\n"
                    "      - {y: 0.0, z: 0.0, area: 1.0, material: 1}\n"),
        ("material", "materials:\n"
                     "  - {id: 1, type: Elastic1D, E: 1.0}\n"),
    ])
    def test_id_de_familia_duplicado(self, tmp_path, familia, bloque):
        """Un id repetido en una familia reemplazaba en silencio al objeto
        anterior: los marcos con ``section: 7`` tomaban otra sección (o las
        fibras con ``material: 1`` otro material) distinta de la declarada
        primero. Revisión adversaria 2026-09-29."""
        cabecera = bloque.split("\n", 1)[0] + "\n"
        texto = _VOLADIZO_2D.replace(cabecera, bloque, 1)
        with pytest.raises(YamlValidationError,
                           match=rf"id de {familia} duplicado \({7 if familia == 'sección' else 1}\)"):
            YamlParser(str(_escribir(tmp_path, texto))).parse()

    def test_ids_de_familias_distintas_son_independientes(self, tmp_path):
        """Sección 1 y material 1 conviven: cada familia tiene su espacio de ids."""
        texto = _VOLADIZO_2D.replace("  - id: 7\n", "  - id: 1\n").replace("section: 7", "section: 1")
        dom, res = _resolver_yaml(_escribir(tmp_path, texto))
        assert res.converged

    def test_parametro_desconocido_de_la_seccion(self, tmp_path):
        texto = _VOLADIZO_2D.replace("    type: FiberSection\n", "    type: FiberSection\n    Gj: 1.0\n")
        with pytest.raises(YamlValidationError, match="parámetro 'Gj' no aceptado por 'FiberSection'"):
            YamlParser(str(_escribir(tmp_path, texto))).parse()

    def test_material_y_seccion_a_la_vez(self, tmp_path):
        texto = _VOLADIZO_2D.replace("nodes: [2, 3], section: 7",
                                     "nodes: [2, 3], section: 7, material: 1, A: 0.15, I: 3.0e-3")
        with pytest.raises(ValueError, match="no ambos"):
            YamlParser(str(_escribir(tmp_path, texto))).parse()

    def test_marco_inelastico_con_material_se_rechaza(self, tmp_path):
        texto = _VOLADIZO_2D.replace(
            "  - {id: 1, type: Elastic1D, E: 2.0e11, density: 7850.0}\n",
            "  - {id: 1, type: Elastic1D, E: 2.0e11, density: 7850.0}\n"
            "  - {id: 2, type: Elastoplastic1D, E: 2.0e11, sigma_y: 2.5e8}\n",
        ).replace("nodes: [2, 3], section: 7", "nodes: [2, 3], material: 2, A: 0.15, I: 3.0e-3")
        with pytest.raises(ValueError, match="FiberSection"):
            YamlParser(str(_escribir(tmp_path, texto))).parse()


# ---------------------------------------------------------------------------
# Exportación VTK
# ---------------------------------------------------------------------------

class TestVtk:

    def test_marcos_con_seccion_exportan_fuerzas_internas(self, tmp_path):
        meshio = pytest.importorskip("meshio")
        from solidum.entry import run
        from solidum.math.solvers import NonlinearSolver
        from solidum.utils.vtk_exporter import VtkExporter

        mat = Elastoplastic1D(E=E_S, sigma_y=250e6, H=2e9)
        sec = FiberSection(patches=[_rect(mat, ny=20)])
        dom = Domain()
        nodos = [dom.add_node(i + 1, [float(i), 0.0]) for i in range(5)]
        for k in range(4):
            dom.add_element(Frame2DEuler(k + 1, [nodos[k], nodos[k + 1]], section=sec))
        for d in ("ux", "uy", "rz"):
            nodos[0].fix_dof(d, 0.0)
        dom.generate_equation_numbers(verbose=False)
        F = np.zeros(dom.total_dofs)
        M_p = 250e6 * B * H * H / 4
        F[nodos[-1].dofs["uy"]] = -0.8 * M_p / 4.0
        asm = Assembler(dom)
        res = run(dom, assembler=asm, solver=NonlinearSolver(asm, num_steps=4), F_applied=F)

        ruta = str(tmp_path / "marcos.vtu")
        VtkExporter(dom).export(ruta, U=res.U, element_forces=res.element_forces)
        malla = meshio.read(ruta)
        assert {b.type for b in malla.cells} == {"line"}
        for comp in ("N", "V", "M"):
            for extremo, k in (("i", 0), ("j", 1)):
                valores = malla.cell_data[f"{comp}_{extremo}"][0]
                esperados = [res.element_forces[e].components[comp][k] for e in range(1, 5)]
                np.testing.assert_allclose(valores, esperados, rtol=1e-12, atol=1e-6)
        # M_i del primer elemento = P·L (equilibrio); los plastificados
        # tienen estado interno.
        np.testing.assert_allclose(malla.cell_data["M_i"][0][0], -0.8 * M_p, rtol=1e-6)
        assert malla.cell_data["Internal_State"][0][0] > 0.0

        # Sin element_forces, internal_forces(U) (sin carga distribuida es lo mismo).
        ruta2 = str(tmp_path / "marcos2.vtu")
        VtkExporter(dom).export(ruta2, U=res.U)
        np.testing.assert_allclose(meshio.read(ruta2).cell_data["M_i"][0],
                                   malla.cell_data["M_i"][0], rtol=1e-12)

    def test_celdas_sin_la_componente_valen_nan(self, tmp_path):
        meshio = pytest.importorskip("meshio")
        from solidum.elements.truss import Truss2D
        from solidum.utils.vtk_exporter import VtkExporter

        n1, n2, n3 = Node(1, [0.0, 0.0]), Node(2, [1.0, 0.0]), Node(3, [1.0, 1.0])
        marco = Frame2DEuler(1, [n1, n2], section=_seccion_elastica())
        barra = Truss2D(2, [n2, n3], _acero(), A=1e-3)
        dom = _dominio(marco, barra)
        U = np.zeros(dom.total_dofs)
        U[n2.dofs["ux"]] = 1e-4
        ruta = str(tmp_path / "mixto.vtu")
        VtkExporter(dom).export(ruta, U=U)
        m_i = meshio.read(ruta).cell_data["M_i"][0]
        assert np.isfinite(m_i[0]) and np.isnan(m_i[1])


# ---------------------------------------------------------------------------
# Barras que desplazan el material del parche: displaces (spec §3, decisión
# del usuario 2026-09-30)
# ---------------------------------------------------------------------------
#
# Sin displaces, las fibras del parche cubren también el sitio de las barras
# y en la posición de cada barra el área cuenta dos veces:
# N_b = (σ_s(ε_b) + σ_c(ε_b))·A_b. Con displaces, la sección añade en esa
# posición una fibra de compensación de área −A_b con el material desplazado,
# con su propia ley y su propia historia.
#
# Tolerancias: las comparaciones frente a sumas escritas a mano difieren sólo
# en el orden de las sumas y en el redondeo de cada producto. Para n
# sumandos, |fl(Σx_i) − Σx_i| ≤ γ_(n+3)·Σ|x_i| con γ_m ≈ m·ε (Higham, §4.2;
# los +3 son los productos σ·A·w); con n ≤ 400 y ε = 2.2e-16, γ < 1e-13.
# De ahí _GAMMA, que multiplica siempre la suma de los valores absolutos de
# los términos, entrada a entrada, nunca el resultado.

_GAMMA = 1.0e-13

_NY, _NZ = 10, 4
# Capas de barras fuera de los centros de celda (caso general): (y, z0, z1, n, área).
_CAPAS = [(-0.21, -0.11, 0.11, 3, 5.1e-4),
          (0.22, -0.11, 0.11, 2, 3.1e-4)]
_SUELTA = (0.03, 0.07, 2.0e-4)          # fibra suelta (y, z, área)

# Historia de carga, descarga, recarga más allá e inversión, e = [ε₀, κ_y, κ_z]:
# en el primer estado el concreto se daña y el acero fluye en tracción
# (y = −0.21) y en compresión (y = 0.22).
_TRAYECTORIA_3D = [np.array(e) for e in ([2.0e-4, -3.0e-3, 1.2e-2],
                                         [1.0e-4, -1.5e-3, 5.0e-3],
                                         [3.0e-4, -5.0e-3, 1.8e-2],
                                         [-2.0e-4, 4.0e-3, -1.4e-2])]


def _materiales_inelasticos():
    """Concreto con daño (ablandamiento) y acero elastoplástico."""
    concreto = IsotropicDamage1D(E=E_C, kappa_0=1.0e-4, alpha=0.99, density=RHO_C)
    acero = Elastoplastic1D(E=E_S, sigma_y=420e6, H=0.01 * E_S, density=RHO_S)
    return concreto, acero


def _seccion_armada(concreto, acero, displaces=True, GJ=2.0e7):
    """Parche de concreto ``_NY × _NZ`` con dos capas de barras y una fibra
    suelta de acero; con ``displaces`` todas desplazan el concreto."""
    d = concreto if displaces else None
    return FiberSection(
        patches=[_rect(concreto, ny=_NY, nz=_NZ)]
        + [BarLayer(acero, start=(y, z0), end=(y, z1), n=n, area=a, displaces=d)
           for y, z0, z1, n, a in _CAPAS],
        fibers=[Fiber(*_SUELTA, acero, displaces=d)],
        GJ=GJ,
    )


def _barras() -> np.ndarray:
    """``(y, z, A)`` de cada barra de :func:`_seccion_armada`, en el orden de
    declaración (misma fórmula de posición que ``BarLayer``)."""
    filas = []
    for y, z0, z1, n, a in _CAPAS:
        for t in np.linspace(0.0, 1.0, n):
            filas.append((y, z0 + t * (z1 - z0), a))
    filas.append(_SUELTA)
    return np.array(filas)


class _SeccionSinDobleConteoAMano:
    """La sección de :func:`_seccion_armada` sin doble conteo, evaluada a mano
    fibra a fibra con los materiales::

        s = Σ_celdas σ_c(ε_k) A_k w_k + Σ_barras [σ_s(ε_b) − σ_c(ε_b)] A_b w_b

    con ``w = [1, z, −y]`` y ``ε = ε₀ + z κ_y − y κ_z``. El concreto restado en
    cada barra lleva su propia historia (la de un punto de concreto en esa
    posición). ``evaluar`` guarda los estados como committed (la trayectoria
    se recorre paso a paso) y devuelve también ``Σ|términos|`` para la cota de
    redondeo."""

    def __init__(self, concreto, acero):
        yc, zc, ac = _rect(concreto, ny=_NY, nz=_NZ).discretize()
        barras = _barras()
        self.celdas = [(y, z, a, concreto) for y, z, a in zip(yc, zc, ac)]
        self.barras = [(y, z, a, acero) for y, z, a in barras]
        self.desplazado = [(y, z, -a, concreto) for y, z, a in barras]
        self.puntos = self.celdas + self.barras + self.desplazado
        self.estados = [None] * len(self.puntos)

    def evaluar(self, e):
        s, k = np.zeros(3), np.zeros((3, 3))
        cota_s, cota_k = np.zeros(3), np.zeros((3, 3))
        s_desplazado = np.zeros(3)
        nuevos = []
        for j, ((y, z, a, mat), st) in enumerate(zip(self.puntos, self.estados)):
            w = np.array([1.0, z, -y])
            sig, et, st_n = mat.compute_state(e[0] + z * e[1] - y * e[2], st)
            s += sig * a * w
            k += et * a * np.outer(w, w)
            cota_s += abs(sig * a) * np.abs(w)
            cota_k += abs(et * a) * np.outer(np.abs(w), np.abs(w))
            if j >= len(self.celdas) + len(self.barras):
                s_desplazado -= sig * a * w       # σ_c(ε_b)·A_b·w_b
            nuevos.append(st_n)
        self.estados = nuevos
        return s, k, cota_s, cota_k, s_desplazado


def _bits_iguales(a, b) -> bool:
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    return a.shape == b.shape and a.tobytes() == b.tobytes()


class TestDisplaces:

    def test_area_desplazada_mayor_que_la_declarada_se_rechaza(self):
        """El área desplazada de un material no puede superar su área
        declarada: la sección quedaría con área neta negativa (hallazgo de
        la revisión del 2026-09-30: una inclusión de 0.2 m² en 0.15 m² de
        concreto daba EA < 0 sin aviso)."""
        concreto = Elastic1D(E=E_C, density=RHO_C)
        blando = Elastic1D(E=1.0e9, density=1000.0)
        with pytest.raises(ValueError, match="área neta negativa"):
            FiberSection(patches=[_rect(concreto, ny=4)],
                         fibers=[Fiber(0.0, 0.0, 0.2, blando, displaces=concreto)])
        # Justo por debajo del área declarada (B·H = 0.15 m²) se admite.
        FiberSection(patches=[_rect(concreto, ny=4)],
                     fibers=[Fiber(0.0, 0.0, 0.149, blando, displaces=concreto)])

    def test_fibras_de_compensacion(self):
        concreto, acero = Elastic1D(E=E_C, density=RHO_C), _acero()
        sec = _seccion_armada(concreto, acero)
        barras = _barras()
        nc, nb = _NY * _NZ, len(barras)
        assert sec.n_fibers == nc + 2 * nb
        assert sec.materials() == [concreto, acero]
        assert len(sec.initial_states()) == sec.n_fibers
        # Al final, tras las fibras declaradas; misma posición, área opuesta.
        np.testing.assert_array_equal(sec.is_compensation, [False] * (nc + nb) + [True] * nb)
        np.testing.assert_array_equal(sec.compensated_index, np.arange(nc, nc + nb))
        comp = sec.is_compensation
        np.testing.assert_array_equal(sec.y[comp], sec.y[nc:nc + nb])
        np.testing.assert_array_equal(sec.z[comp], sec.z[nc:nc + nb])
        np.testing.assert_array_equal(sec.area[comp], -barras[:, 2])
        np.testing.assert_allclose(sec.y[nc:nc + nb], barras[:, 0], rtol=0, atol=1e-16)
        np.testing.assert_allclose(sec.z[nc:nc + nb], barras[:, 1], rtol=0, atol=1e-16)
        # Una fibra de compensación tiene la deformación de la barra.
        eps = sec.strains_3d([1.0e-4, -2.0e-3, 3.0e-3])
        np.testing.assert_allclose(eps[comp], eps[sec.compensated_index], rtol=1e-15)
        # Sólo las primitivas con displaces generan compensación.
        mixta = FiberSection(
            patches=[_rect(concreto),
                     BarLayer(acero, start=(-0.2, -0.1), end=(-0.2, 0.1), n=2, area=1e-4,
                              displaces=concreto)],
            fibers=[Fiber(0.2, 0.0, 1e-4, acero)])
        assert mixta.is_compensation.sum() == 2 and mixta.n_fibers == 20 + 3 + 2
        assert "compensation=2" in repr(mixta)
        assert "compensation" not in repr(FiberSection(patches=[_rect(concreto)]))

    def test_rigidez_elastica_sin_doble_conteo(self):
        """``k_s`` con displaces = la forma a mano con el concreto neto: el
        parche de concreto más ``(E_s − E_c)·A_b`` en cada barra. En particular
        ``EA = E_c·(A_c − ΣA_b) + E_s·ΣA_b`` y
        ``EI_z = E_c·b h³/12·(1 − 1/n_y²) + (E_s − E_c)·ΣA_b y_b²``."""
        concreto, acero = Elastic1D(E=E_C, density=RHO_C), _acero()
        sec = _seccion_armada(concreto, acero)
        yb, zb, ab = _barras().T
        dE = E_S - E_C
        EA = E_C * (B * H - ab.sum()) + E_S * ab.sum()
        EIz = E_C * B * H**3 / 12 * (1 - 1 / _NY**2) + dE * (ab * yb**2).sum()
        EIy = E_C * H * B**3 / 12 * (1 - 1 / _NZ**2) + dE * (ab * zb**2).sum()
        # El parche es simétrico respecto del eje: los acoplamientos los dan
        # sólo las barras (k₀₁ = ΣE A z, k₀₂ = −ΣE A y, k₁₂ = −ΣE A y z).
        Sz, Sy, Syz = dE * (ab * zb).sum(), dE * (ab * yb).sum(), dE * (ab * yb * zb).sum()
        k3_ref = np.array([[EA, Sz, -Sy], [Sz, EIy, -Syz], [-Sy, -Syz, EIz]])
        _, k3, _ = sec.response_3d([0.0, 0.0, 0.0], sec.initial_states())
        _, k2, _ = sec.response_2d([0.0, 0.0], sec.initial_states())
        # Cota de redondeo: Σ|E_i A_i| |w_i||w_i|ᵀ (E por fibra: concreto,
        # acero de las barras y concreto de las compensaciones).
        E_i = np.array([E_C] * (_NY * _NZ) + [E_S] * len(ab) + [E_C] * len(ab))
        W_abs = np.abs(np.column_stack([np.ones(sec.n_fibers), sec.z, -sec.y]))
        cota = W_abs.T @ (W_abs * np.abs(E_i * sec.area)[:, None])
        assert np.all(np.abs(k3 - k3_ref) <= _GAMMA * cota)
        idx = np.ix_([0, 2], [0, 2])
        assert np.all(np.abs(k2 - k3_ref[idx]) <= _GAMMA * cota[idx])
        # Sin displaces, la diferencia es exactamente el concreto contado dos
        # veces: E_c·ΣA_b w_b w_bᵀ.
        sin = _seccion_armada(concreto, acero, displaces=False)
        _, k3_sin, _ = sin.response_3d([0.0, 0.0, 0.0], None)
        Wb = np.column_stack([np.ones(len(ab)), zb, -yb])
        doble = E_C * Wb.T @ (Wb * ab[:, None])
        assert np.all(np.abs(k3_sin - k3 - doble) <= 2 * _GAMMA * cota)

    def test_inelastica_igual_a_la_suma_a_mano_fibra_a_fibra(self):
        """Carga, descarga, recarga e inversión con concreto dañado y acero
        plastificado: ``s`` y ``k_s`` con displaces (2D y 3D) son los de la
        sección sin doble conteo evaluada a mano, y la diferencia con la
        sección sin displaces es exactamente el concreto contado dos veces,
        ``Σ σ_c(ε_b)·A_b·w_b`` con la historia propia de cada barra."""
        for dim in (2, 3):
            concreto, acero = _materiales_inelasticos()
            sec = _seccion_armada(concreto, acero)
            sin = _seccion_armada(concreto, acero, displaces=False)
            mano = _SeccionSinDobleConteoAMano(concreto, acero)
            resp = sec.response_3d if dim == 3 else sec.response_2d
            resp_sin = sin.response_3d if dim == 3 else sin.response_2d
            sel = [0, 1, 2] if dim == 3 else [0, 2]
            st, st_sin = sec.initial_states(), sin.initial_states()
            for paso, e3 in enumerate(_TRAYECTORIA_3D):
                e = e3 if dim == 3 else e3[[0, 2]]
                if dim == 2:
                    e3 = np.array([e3[0], 0.0, e3[2]])
                s, k, st = resp(e, st)
                s_sin, _, st_sin = resp_sin(e, st_sin)
                s_m, k_m, cota_s, cota_k, s_desp = mano.evaluar(e3)
                s_m, k_m, s_desp = s_m[sel], k_m[np.ix_(sel, sel)], s_desp[sel]
                cota_s, cota_k = cota_s[sel], cota_k[np.ix_(sel, sel)]
                assert np.all(np.abs(s - s_m) <= _GAMMA * cota_s), (dim, paso)
                assert np.all(np.abs(k - k_m) <= _GAMMA * cota_k), (dim, paso)
                assert np.all(np.abs(s_sin - s - s_desp) <= 2 * _GAMMA * cota_s), (dim, paso)
                # Régimen realmente inelástico, también en el concreto restado.
                danos = [st[i]["damage"] for i in np.flatnonzero(sec.is_compensation)]
                assert min(danos) > 0.0 and max(danos) > 0.5
                assert max(st[i]["alpha"] for i in range(_NY * _NZ, sec.n_fibers)
                           if "alpha" in st[i]) > 0.0

    def test_inelastica_igual_a_la_seccion_sin_la_zona_de_las_barras(self):
        """2D, barras en el centro de una capa del parche: la sección con
        displaces es la que discretiza el concreto sin la zona de las barras
        (esa capa con área ``A_capa − n·A_b``), porque el concreto de la capa
        y las fibras de compensación comparten deformación e historia."""
        concreto, acero = _materiales_inelasticos()
        ny = 20
        yc, _, ac = _rect(concreto, ny=ny).discretize()
        capas = {1: (3, 5.1e-4), ny - 2: (2, 3.1e-4)}         # capa: (n, A_b)
        sec = FiberSection(patches=[_rect(concreto, ny=ny)] + [
            BarLayer(acero, start=(yc[k], -0.11), end=(yc[k], 0.11), n=n, area=a,
                     displaces=concreto) for k, (n, a) in capas.items()])
        fibras = [Fiber(yc[k], 0.0, ac[k] - capas.get(k, (0, 0.0))[0] * capas.get(k, (0, 0.0))[1],
                        concreto) for k in range(ny)]
        for k, (n, a) in capas.items():
            fibras += [Fiber(yc[k], 0.0, a, acero) for _ in range(n)]
        neta = FiberSection(fibers=fibras)
        # Cota de redondeo con |σ| ≤ σ_y + H·α (acero, α < 5e-2 aquí) y |E_t| ≤ E_s,
        # que acotan también al concreto (|σ_c| ≤ E_c·κ₀, |E_t,c| ≤ E_c).
        sig_max, E_max = 420e6 + 0.01 * E_S * 5e-2, E_S
        W_abs = np.abs(np.column_stack([np.ones(sec.n_fibers), -sec.y]))
        cota_s = sig_max * (W_abs.T @ np.abs(sec.area))
        cota_k = E_max * (W_abs.T @ (W_abs * np.abs(sec.area)[:, None]))
        st, st_n = sec.initial_states(), neta.initial_states()
        for e3 in _TRAYECTORIA_3D:
            e = e3[[0, 2]]
            s, k, st = sec.response_2d(e, st)
            s_n, k_n, st_n = neta.response_2d(e, st_n)
            assert np.all(np.abs(s - s_n) <= 2 * _GAMMA * cota_s)
            assert np.all(np.abs(k - k_n) <= 2 * _GAMMA * cota_k)
        assert max(d["alpha"] for d in st if "alpha" in d) < 5e-2

    def test_tangente_con_compensacion_frente_a_diferencias_finitas(self):
        concreto, acero = _materiales_inelasticos()
        sec = _seccion_armada(concreto, acero)
        e1 = _TRAYECTORIA_3D[0]
        _, _, st = sec.response_3d(e1, sec.initial_states())
        e2 = 1.07 * e1
        _, k, _ = sec.response_3d(e2, st)
        J = _fd_jacobian(lambda e: sec.response_3d(e, st)[0], e2, 1e-7 * np.abs(e2).max())
        np.testing.assert_allclose(k, J, rtol=0, atol=1e-6 * np.abs(k).max())
        np.testing.assert_array_equal(k, k.T)

    def test_masa_y_area_bruta(self):
        concreto, acero = Elastic1D(E=E_C, density=RHO_C), _acero()
        sec = _seccion_armada(concreto, acero)
        sin = _seccion_armada(concreto, acero, displaces=False)
        yb, zb, ab = _barras().T
        # Área bruta: la del rectángulo; sin displaces las barras se suman.
        # Redondeo: ~50 sumandos positivos, γ·Σ|A_i|/A < 1e-14.
        np.testing.assert_allclose(sec.total_area(), B * H, rtol=1e-14)
        np.testing.assert_allclose(sin.total_area(), B * H + ab.sum(), rtol=1e-14)
        mp = sec.mass_properties()
        dr = RHO_S - RHO_C
        np.testing.assert_allclose(mp["rhoA"], RHO_C * (B * H - ab.sum()) + RHO_S * ab.sum(),
                                   rtol=1e-13)
        np.testing.assert_allclose(
            mp["rhoIz"], RHO_C * B * H**3 / 12 * (1 - 1 / _NY**2) + dr * (ab * yb**2).sum(),
            rtol=1e-13)
        np.testing.assert_allclose(
            mp["rhoIy"], RHO_C * H * B**3 / 12 * (1 - 1 / _NZ**2) + dr * (ab * zb**2).sum(),
            rtol=1e-13)
        # En el marco: masa total ρA·L y fuerza de cuerpo con el área bruta.
        el = Frame2DEuler(1, _nodos_2d(), section=sec)
        for lumping in ("consistent", "lumped"):
            M = el.compute_mass_matrix(lumping)
            t = np.zeros(6)
            t[0::3] = 1.0
            np.testing.assert_allclose(t @ M @ t, mp["rhoA"] * el.L0, rtol=1e-12)
        b = np.array([0.0, -5.0e3])
        ref = Frame2DEuler(1, _nodos_2d(), _acero(), A=B * H, I=1e-3)
        f_ref = ref.compute_body_load(b)
        # La componente axial es nula: queda el redondeo, del orden de ε·max|f|,
        # que cambia con la BLAS de cada plataforma (en el CI de Linux 5.6e-14
        # frente a 3.3e-14): tolerancia absoluta escalada, no sólo relativa.
        np.testing.assert_allclose(el.compute_body_load(b), f_ref, rtol=1e-14,
                                   atol=1e-12 * np.abs(f_ref).max())

    def test_sin_displaces_nada_cambia_bit_a_bit(self):
        """Sin displaces (omitido o ``None``) la sección es la de antes: las
        fibras declaradas tal cual y la misma aritmética de spec §4
        (``s = Wᵀ(σ∘A)``, ``k_s = Wᵀ diag(E_t∘A) W`` simetrizada), bit a bit."""
        concreto, acero = _materiales_inelasticos()
        prims = ([_rect(concreto, ny=_NY, nz=_NZ)]
                 + [BarLayer(acero, start=(y, z0), end=(y, z1), n=n, area=a)
                    for y, z0, z1, n, a in _CAPAS])
        suelta = Fiber(*_SUELTA, acero)
        omitido = FiberSection(patches=prims, fibers=[suelta], GJ=2.0e7)
        con_none = _seccion_armada(concreto, acero, displaces=False)
        partes = [p.discretize() for p in prims + [suelta]]
        y, z, A = (np.concatenate([p[j] for p in partes]) for j in range(3))
        mats = [p.material for p, d in zip(prims + [suelta], partes) for _ in d[0]]
        for sec in (omitido, con_none):
            assert _bits_iguales(sec.y, y) and _bits_iguales(sec.z, z) and _bits_iguales(sec.area, A)
            assert not sec.is_compensation.any() and sec.compensated_index.size == 0
            assert sec.materials() == [concreto, acero]
            assert sec.total_area() == float(A.sum())
            rho_a = np.array([m.density for m in mats]) * A
            assert sec.mass_properties() == {"rhoA": float(rho_a.sum()),
                                             "rhoIy": float((rho_a * z * z).sum()),
                                             "rhoIz": float((rho_a * y * y).sum())}
            for dim in (2, 3):
                W = (np.column_stack([np.ones(len(A)), -y]) if dim == 2
                     else np.column_stack([np.ones(len(A)), z, -y]))
                resp = sec.response_2d if dim == 2 else sec.response_3d
                st = st_ref = sec.initial_states()
                for e3 in _TRAYECTORIA_3D:
                    e = e3[[0, 2]] if dim == 2 else e3
                    s, k, st = resp(e, st)
                    eps = (W @ e).tolist()
                    out = [m.compute_state(eps[i], st_ref[i]) for i, m in enumerate(mats)]
                    sig = np.array([o[0] for o in out])
                    tan = np.array([o[1] for o in out])
                    st_ref = [o[2] for o in out]
                    s_ref = W.T @ (sig * A)
                    k_ref = W.T @ (W * (tan * A)[:, None])
                    k_ref = 0.5 * (k_ref + k_ref.T)
                    assert _bits_iguales(s, s_ref) and _bits_iguales(k, k_ref)
                    assert st == st_ref

    def test_guarda_de_unilaterales_no_cuenta_las_fibras_de_compensacion(self):
        """Las fibras de compensación restan rigidez, no la aportan: una
        sección de fibras declaradas unilaterales se rechaza aunque desplacen
        un material bilateral."""
        from solidum.materials.cable_1d import CableMaterial1D
        cable = CableMaterial1D(E=E_S, density=RHO_S)
        otro_cable = CableMaterial1D(E=1.5e11, density=RHO_S)
        concreto = Elastic1D(E=E_C, density=RHO_C)
        with pytest.raises(ValueError, match="unilateral"):
            FiberSection(fibers=[Fiber(0.1, 0.0, 1e-4, cable, displaces=concreto),
                                 Fiber(-0.1, 0.0, 1e-4, cable)])
        with pytest.raises(ValueError, match="unilateral"):
            FiberSection(patches=[_rect(cable, ny=4),
                                  BarLayer(otro_cable, start=(-0.2, 0.0), end=(0.2, 0.0),
                                           n=2, area=1e-4, displaces=cable)])
        # Tendones unilaterales dentro de concreto, con descuento: se admiten.
        # Comprimida, el tendón no aporta y su hueco no tiene concreto; en
        # tracción el tendón sustituye al concreto. Redondeo: 8 sumandos.
        a_t = 4.0e-4
        sec = FiberSection(patches=[_rect(concreto, ny=4),
                                    BarLayer(cable, start=(-0.2, -0.1), end=(-0.2, 0.1), n=2,
                                             area=a_t, displaces=concreto)])
        _, k_c, _ = sec.response_2d([-1e-4, 0.0], sec.initial_states())
        np.testing.assert_allclose(k_c[0, 0], E_C * (B * H - 2 * a_t), rtol=1e-14)
        _, k_t, _ = sec.response_2d([1e-4, 0.0], sec.initial_states())
        np.testing.assert_allclose(k_t[0, 0], E_C * (B * H - 2 * a_t) + E_S * 2 * a_t, rtol=1e-14)

    @pytest.mark.parametrize("fabrica, patron", [
        (lambda s: BarLayer(s, start=(0, 0), end=(0, 1), n=2, area=1e-4, displaces=s),
         "mismo material"),
        (lambda s: Fiber(0.0, 0.0, 1e-4, s, displaces=s), "mismo material"),
        (lambda s: BarLayer(s, start=(0, 0), end=(0, 1), n=2, area=1e-4,
                            displaces=Elastic2D(E=E_C, nu=0.2)), r"\('displaces'\) declara STRAIN_DIM=3"),
        (lambda s: Fiber(0.0, 0.0, 1e-4, s, displaces=object()),
         "'displaces' debe ser un material con compute_state"),
    ])
    def test_displaces_invalido(self, fabrica, patron):
        with pytest.raises(ValueError, match=patron):
            fabrica(_acero())

    def test_material_desplazado_ausente_de_la_seccion(self):
        """``displaces`` debe ser el material (el mismo objeto) de alguna fibra
        declarada: si no, restaría un área que no existe."""
        concreto = Elastic1D(E=E_C, density=RHO_C)
        otro = Elastic1D(E=E_C, density=RHO_C)          # igual, pero otro objeto
        with pytest.raises(ValueError, match="ninguna fibra declarada"):
            FiberSection(patches=[_rect(concreto),
                                  BarLayer(_acero(), start=(-0.2, -0.1), end=(-0.2, 0.1), n=2,
                                           area=1e-4, displaces=otro)])

    def test_vtk_internal_state_sin_las_fibras_de_compensacion(self):
        """``Internal_State`` promedia sólo las fibras declaradas: una fibra de
        compensación no es un punto material."""
        from solidum.utils.vtk_exporter import _state_var_avg
        concreto, acero = _materiales_inelasticos()
        sec = _seccion_armada(concreto, acero)
        L, kappa = 2.0, 1.2e-2
        el = Frame2DEuler(1, _nodos_2d((0.0, 0.0), (L, 0.0)), section=sec, n_ip=3)
        el.compute_element_state(np.array([0.0, 0.0, -kappa * L / 2, 2e-4 * L, 0.0, kappa * L / 2]))
        el.commit_state()
        n_decl = int((~sec.is_compensation).sum())

        def escalar(st):
            return st["damage"] if "damage" in st else st["alpha"]

        declaradas = [escalar(st) for sv in el.state.vars for st in sv[:n_decl] if st is not None]
        todas = [escalar(st) for sv in el.state.vars for st in sv if st is not None]
        assert _state_var_avg(el) == pytest.approx(np.mean(declaradas), rel=1e-14)
        assert np.mean(todas) != pytest.approx(np.mean(declaradas), rel=1e-6)
        assert sec.material_point_states(el.state.vars[0]) == el.state.vars[0][:n_decl]


_ARMADA_YAML = """\
nodes:
  - {id: 1, coords: [0.0, 0.0]}
  - {id: 2, coords: [2.0, 0.0]}
  - {id: 3, coords: [4.0, 0.0]}
materials:
  - {id: 1, type: Elastic1D, E: 3.0e10, density: 2400.0}
  - {id: 2, type: Elastic1D, E: 2.0e11, density: 7850.0}
sections:
  - id: 7
    type: FiberSection
    patches:
      - {shape: rectangle, material: 1, y: [-0.25, 0.25], z: [-0.15, 0.15], ny: 20, nz: 1}
      - {shape: layer, material: 2, start: [-0.2, -0.1], end: [-0.2, 0.1], n: 3, area: 5.1e-4, displaces: 1}
      - {shape: layer, material: 2, start: [0.2, -0.1], end: [0.2, 0.1], n: 3, area: 5.1e-4}
    fibers:
      - {y: 0.0, z: 0.0, area: 2.0e-4, material: 2, displaces: 1}
elements:
  - {id: 1, type: Frame2DEuler, nodes: [1, 2], section: 7}
  - {id: 2, type: Frame2DEuler, nodes: [2, 3], section: 7}
boundary_conditions_by_node:
  - {node_id: 1, ux: 0.0, uy: 0.0, rz: 0.0}
point_loads:
  - {node_id: 3, ux: 1.0e6, uy: -1.0e5}
solver:
  type: LinearSolver
"""


class TestDisplacesYaml:

    def _simetrico(self):
        """Las dos capas desplazan el concreto: sección simétrica, sin
        acoplamiento axial–flexión."""
        return _ARMADA_YAML.replace("n: 3, area: 5.1e-4}", "n: 3, area: 5.1e-4, displaces: 1}")

    def test_voladizo_con_displaces(self, tmp_path):
        """Flecha y alargamiento de punta con ``EI`` y ``EA`` del concreto neto
        (elementos cúbicos: nodalmente exactos para carga en la punta)."""
        dom, res = _resolver_yaml(_escribir(tmp_path, self._simetrico()))
        sec = next(iter(dom.elements.values())).section
        assert sec.is_compensation.sum() == 7
        assert [m.E for m in sec.materials()] == [E_C, E_S]
        A_s = 6 * 5.1e-4 + 2.0e-4
        EA = E_C * (B * H - A_s) + E_S * A_s
        EI = E_C * B * H**3 / 12 * (1 - 1 / 20**2) + (E_S - E_C) * 6 * 5.1e-4 * 0.2**2
        L = 4.0
        np.testing.assert_allclose(res.U[dom.nodes[3].dofs["ux"]], 1.0e6 * L / EA, rtol=1e-10)
        np.testing.assert_allclose(res.U[dom.nodes[3].dofs["uy"]], -1.0e5 * L**3 / (3 * EI),
                                   rtol=1e-10)

    def test_displaces_es_opcional_por_primitiva(self, tmp_path):
        parser = YamlParser(str(_escribir(tmp_path, _ARMADA_YAML)))
        parser.parse()
        sec = parser.sections[7]
        # La primera capa y la fibra suelta desplazan; la segunda capa no.
        assert sec.is_compensation.sum() == 4
        np.testing.assert_array_equal(sec.compensated_index, [20, 21, 22, 26])

    @pytest.mark.parametrize("cambio, contexto", [
        (("area: 5.1e-4, displaces: 1}", "area: 5.1e-4, displaces: true}"),
         r"patches\[1\]\.displaces"),
        (("material: 2, displaces: 1}", "material: 2, displaces: true}"), r"fibers\[0\]\.displaces"),
        (("{shape: rectangle, material: 1,", "{shape: rectangle, material: true,"),
         r"patches\[0\]\.material"),
    ])
    def test_un_booleano_no_es_un_id(self, tmp_path, cambio, contexto):
        """``displaces: true`` (leído como bandera) no puede resolver en silencio
        al material de id 1 por ``True == 1``: es un error que nombra el
        contexto. Revisión adversaria 2026-09-30: sin la guarda del lector la
        sección se construía con 4 fibras de compensación del material 1."""
        texto = _ARMADA_YAML.replace(*cambio, 1)
        assert texto != _ARMADA_YAML
        with pytest.raises(YamlValidationError) as exc:
            YamlParser(str(_escribir(tmp_path, texto))).parse()
        msg = str(exc.value)
        assert re.search(contexto, msg), msg
        assert "True no es un id de material" in msg, msg

    @pytest.mark.parametrize("cambio, patrones", [
        (("area: 5.1e-4, displaces: 1}", "area: 5.1e-4, displaces: 9}"),
         [r"sections\[0\] \(id=7\), patches\[1\]\.displaces",
          r"referencia a material inexistente \(9\)"]),
        (("material: 2, displaces: 1}", "material: 2, displaces: 9}"),
         [r"fibers\[0\]\.displaces", r"referencia a material inexistente \(9\)"]),
        (("area: 5.1e-4, displaces: 1}", "area: 5.1e-4, displaces: 2}"),
         [r"patches\[1\]: BarLayer: 'displaces' es el mismo material"]),
        (("area: 5.1e-4, displaces: 1}", "area: 5.1e-4, displaces: 3}"),
         [r"patches\[1\]: BarLayer: el material Elastic2D \('displaces'\) declara STRAIN_DIM=3"]),
        (("area: 5.1e-4, displaces: 1}", "area: 5.1e-4, displaces: 4}"),
         ["no es el de ninguna fibra declarada"]),
        (("ny: 20, nz: 1}", "ny: 20, nz: 1, displaces: 1}"),
         [r"parámetros \['displaces'\] no aceptados por RectangularPatch"]),
    ])
    def test_errores_de_displaces(self, tmp_path, cambio, patrones):
        texto = _ARMADA_YAML.replace(*cambio, 1).replace(
            "  - {id: 2, type: Elastic1D, E: 2.0e11, density: 7850.0}\n",
            "  - {id: 2, type: Elastic1D, E: 2.0e11, density: 7850.0}\n"
            "  - {id: 3, type: Elastic2D, E: 3.0e10, nu: 0.2}\n"
            "  - {id: 4, type: Elastic1D, E: 3.0e10, density: 2400.0}\n")
        with pytest.raises(YamlValidationError) as exc:
            YamlParser(str(_escribir(tmp_path, texto))).parse()
        for patron in patrones:
            assert re.search(patron, str(exc.value)), (patron, str(exc.value))
