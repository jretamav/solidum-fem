"""Tests del elemento ``Hex8Thermal`` (Etapa 8).

Cubre el bloque ``acceptance`` de ``docs/specs/Hex8Thermal.md``.

Al ser spec de extensión, estos tests se concentran en lo que **cambia**
respecto al ``Quad4Thermal``: dimensión 3, ausencia de espesor, flujo por
cara y los 4 modos hourglass de la cuadratura reducida. Lo compartido ya
está blindado en ``test_quad4_thermal.py`` sobre la misma base
``_ThermalSolid``.

El test central es ``TestCrossCheck2D3D``: el análogo térmico del
cross-check 3D↔2D `plane_strain` que en A.bis blindó los materiales.
"""
import numpy as np
import pytest

from solidum.core.domain import Domain
from solidum.core.node import Node
from solidum.elements.thermal.hex8_thermal import Hex8Thermal
from solidum.elements.thermal.quad4_thermal import Quad4Thermal
from solidum.materials.elastic_3d import Elastic3D
from solidum.materials.thermal_conduction import ThermalConduction
from solidum.math.assembly import Assembler
from solidum.math.solvers.linear import LinearSolver

K_ACERO = 45.0
C_ACERO = 460.0
RHO_ACERO = 7850.0

_CUBO = [
    [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0],
]


def _material(**kw):
    return ThermalConduction(k=K_ACERO, dim=3, **kw)


def _cubo(material=None, quadrature=None, escala=1.0):
    nodes = [Node(i + 1, [c * escala for c in xyz]) for i, xyz in enumerate(_CUBO)]
    mat = material if material is not None else _material()
    return Hex8Thermal(1, nodes, mat, quadrature=quadrature)


def _distorsionado():
    coords = [
        [0.0, 0.0, 0.0], [1.2, 0.1, 0.05], [1.15, 1.3, 0.0], [0.1, 1.1, 0.1],
        [0.05, 0.1, 1.1], [1.25, 0.0, 1.2], [1.1, 1.2, 1.15], [0.0, 1.15, 1.05],
    ]
    nodes = [Node(i + 1, c) for i, c in enumerate(coords)]
    return Hex8Thermal(1, nodes, _material())


class TestContrato:
    def test_declara_dimension_3(self):
        el = _cubo()
        assert el.DOF_NAMES == ["T"]
        assert el.FLUX_DIM == 3
        assert el.N_NODES == 8
        assert el.N_INTEGRATION_POINTS == 8

    def test_registrado(self):
        from solidum.registry import ElementRegistry

        assert "Hex8Thermal" in ElementRegistry.names()

    def test_sin_thickness(self):
        """El volumen sale de la geometría, como en los sólidos 3D."""
        el = _cubo()
        assert el.thickness == 1.0
        # Escalar la geometría ×2 multiplica el volumen por 8.
        Q = 1000.0
        assert _cubo(escala=2.0).compute_body_source(Q).sum() == pytest.approx(
            Q * 8.0, rel=1e-11
        )

    def test_caras_paritarias_con_hex8_mecanico(self):
        from solidum.elements.solid_3d.hex8 import Hex8

        assert Hex8Thermal.FACE_NODES == Hex8.FACE_NODES


class TestMatrizDeConductividad:
    def test_simetria(self):
        for el in (_cubo(), _distorsionado()):
            Ke = el.compute_conductivity_matrix()
            np.testing.assert_allclose(Ke, Ke.T, atol=1e-11)

    def test_modo_nulo_y_rango(self):
        """rango 7 de 8: un único modo nulo, la temperatura uniforme."""
        for el in (_cubo(), _distorsionado()):
            Ke = el.compute_conductivity_matrix()
            np.testing.assert_allclose(Ke @ np.ones(8), np.zeros(8), atol=1e-9)
            assert np.linalg.matrix_rank(Ke, tol=1e-9) == 7

            eig = np.linalg.eigvalsh(Ke)
            assert np.min(eig) > -1e-9
            assert int(np.sum(np.abs(eig) < 1e-9)) == 1

    def test_hourglass_con_integracion_reducida(self):
        """4 modos hourglass, frente al único del Quad4Thermal reducido.

        Documenta la limitación; no la corrige (sin estabilización).
        """
        el = _cubo(quadrature="hex_1x1x1")
        Ke = el.compute_conductivity_matrix()
        assert np.linalg.matrix_rank(Ke, tol=1e-9) == 3
        n_nulos = int(np.sum(np.abs(np.linalg.eigvalsh(Ke)) < 1e-9))
        assert n_nulos == 5  # 1 físico + 4 hourglass

    def test_escala_con_la_conductividad(self):
        """K_e es lineal en k: triplicar k triplica la matriz.

        El `atol` escalado es necesario porque las entradas fuera de la
        diagonal principal del cubo regular son ceros numéricos (~1e-15):
        compararlas por tolerancia relativa no tiene sentido.
        """
        e1 = _cubo(material=ThermalConduction(k=K_ACERO, dim=3))
        e2 = _cubo(material=ThermalConduction(k=3 * K_ACERO, dim=3))
        K1 = e1.compute_conductivity_matrix()
        np.testing.assert_allclose(
            e2.compute_conductivity_matrix(),
            3.0 * K1,
            rtol=1e-12,
            atol=1e-12 * np.max(np.abs(K1)),
        )


class TestPatchTest:
    def test_campo_lineal_3d_exacto(self):
        """T = a + b·x + c·y + d·z ⇒ ∇T = (b, c, d) en todos los Gauss."""
        a, b, c, d = 2.0, 3.0, 5.0, -7.0
        for el in (_cubo(), _distorsionado()):
            coords = np.array([n.coordinates[:3] for n in el.nodes])
            T_e = (
                a + b * coords[:, 0] + c * coords[:, 1] + d * coords[:, 2]
            )
            for p in el.points:
                B, _ = el._gradient(p, coords)
                np.testing.assert_allclose(B @ T_e, [b, c, d], rtol=1e-10)

    def test_flujo_anisotropo_3d(self):
        k = np.diag([50.0, 20.0, 5.0])
        el = _cubo(material=ThermalConduction(k=k))
        coords = np.array([n.coordinates[:3] for n in el.nodes])
        T_e = coords[:, 0] + coords[:, 1] + coords[:, 2]

        for idx, n in enumerate(el.nodes):
            n.dofs["T"] = idx
        gs = el.compute_gauss_state(T_e)

        for g in range(len(el.points)):
            np.testing.assert_allclose(
                gs["flux"][g], [-50.0, -20.0, -5.0], rtol=1e-10
            )


class TestCapacidad:
    def _mat(self):
        return ThermalConduction(
            k=K_ACERO, c=C_ACERO, density=RHO_ACERO, dim=3
        )

    def test_capacidad_total_conservada(self):
        el = _cubo(material=self._mat())
        esperado = RHO_ACERO * C_ACERO * 1.0  # V_e = 1
        for modo in ("consistent", "lumped"):
            assert el.compute_capacity_matrix(modo).sum() == pytest.approx(
                esperado, rel=1e-11
            )

    def test_lumped_diagonal_positiva(self):
        C = _cubo(material=self._mat()).compute_capacity_matrix("lumped")
        np.testing.assert_allclose(C, np.diag(np.diag(C)), atol=1e-10)
        assert np.all(np.diag(C) > 0.0)

    def test_default_es_lumped(self):
        el = _cubo(material=self._mat())
        np.testing.assert_allclose(
            el.compute_capacity_matrix(),
            el.compute_capacity_matrix("lumped"),
            rtol=1e-14,
        )


class TestFlujoPorCara:
    def test_las_seis_caras(self):
        """Σf = -q̄·A con reparto A/4 y cero en los nodos fuera de la cara."""
        q_bar = 100.0
        el = _cubo()
        for face in range(6):
            f = el.compute_face_flux(face, q_bar)
            assert f.sum() == pytest.approx(-q_bar * 1.0, rel=1e-11)

            activos = [i for i in range(8) if abs(f[i]) > 1e-12]
            assert sorted(activos) == sorted(el.FACE_NODES[face])
            for i in activos:
                assert f[i] == pytest.approx(-q_bar * 1.0 / 4.0, rel=1e-11)

    def test_signo_saliente_extrae_energia(self):
        el = _cubo()
        assert el.compute_face_flux(0, 100.0).sum() < 0.0
        assert el.compute_face_flux(0, -100.0).sum() > 0.0

    def test_area_de_cara_escalada(self):
        el = _cubo(escala=3.0)  # caras de 3×3 = 9
        assert el.compute_face_flux(0, 10.0).sum() == pytest.approx(
            -10.0 * 9.0, rel=1e-11
        )

    def test_cara_invalida(self):
        with pytest.raises(ValueError, match="face=6"):
            _cubo().compute_face_flux(6, 1.0)


class TestFuenteVolumetrica:
    def test_uniforme_exacta(self):
        Q = 1000.0
        assert _cubo().compute_body_source(Q).sum() == pytest.approx(
            Q * 1.0, rel=1e-11
        )

    def test_invariante_ante_distorsion(self):
        """Σf = Q·V_e con V_e integrado por la propia cuadratura."""
        Q = 500.0
        el = _distorsionado()
        coords = np.array([n.coordinates[:3] for n in el.nodes])
        volumen = sum(
            el._gradient(p, coords)[1] * w
            for p, w in zip(el.points, el.weights)
        )
        assert el.compute_body_source(Q).sum() == pytest.approx(
            Q * volumen, rel=1e-11
        )


class TestValidaciones:
    def test_rechaza_material_mecanico(self):
        nodes = [Node(i + 1, c) for i, c in enumerate(_CUBO)]
        with pytest.raises(TypeError, match="familia térmica"):
            Hex8Thermal(1, nodes, Elastic3D(E=210e9, nu=0.3))

    def test_rechaza_material_2d(self):
        nodes = [Node(i + 1, c) for i, c in enumerate(_CUBO)]
        with pytest.raises(ValueError, match="FLUX_DIM"):
            Hex8Thermal(1, nodes, ThermalConduction(k=K_ACERO, dim=2))

    def test_rechaza_numero_de_nodos(self):
        nodes = [Node(i + 1, c) for i, c in enumerate(_CUBO[:6])]
        with pytest.raises(ValueError, match="requiere 8 nodos"):
            Hex8Thermal(1, nodes, _material())

    def test_jacobiano_invertido(self):
        """Numeración que produce volumen negativo."""
        coords = list(_CUBO)
        coords[0], coords[1] = coords[1], coords[0]
        nodes = [Node(i + 1, c) for i, c in enumerate(coords)]
        el = Hex8Thermal(1, nodes, _material())
        with pytest.raises(ValueError, match="Jacobiano"):
            el.compute_conductivity_matrix()


class TestCrossCheck2D3D:
    """El test clave de esta spec.

    La misma pared plana resuelta en 2D y con una capa de Hex8Thermal con
    las caras z adiabáticas debe dar campos idénticos nodo a nodo. Es el
    análogo del cross-check 3D↔2D `plane_strain` de A.bis, y cubre errores
    en la tercera dimensión sin necesitar mallador — razón por la que el
    benchmark de la esfera hueca pudo diferirse.
    """

    L, H, W = 2.0, 0.5, 0.7
    NX, NY = 6, 2
    T1, T2 = 100.0, 20.0

    def _resolver_2d(self):
        mat = ThermalConduction(k=K_ACERO, dim=2)
        dom = Domain()
        nid = {}
        c = 1
        for j in range(self.NY + 1):
            for i in range(self.NX + 1):
                dom.add_node(c, [self.L * i / self.NX, self.H * j / self.NY])
                nid[(i, j)] = c
                c += 1
        eid = 1
        for j in range(self.NY):
            for i in range(self.NX):
                ns = [
                    dom.nodes[nid[(i, j)]], dom.nodes[nid[(i + 1, j)]],
                    dom.nodes[nid[(i + 1, j + 1)]], dom.nodes[nid[(i, j + 1)]],
                ]
                dom.add_element(Quad4Thermal(eid, ns, mat, thickness=self.W))
                eid += 1
        for j in range(self.NY + 1):
            dom.nodes[nid[(0, j)]].fix_dof("T", self.T1)
            dom.nodes[nid[(self.NX, j)]].fix_dof("T", self.T2)

        asm = Assembler(dom)
        ndof = sum(len(n.dofs) for n in dom.nodes.values())
        T = LinearSolver(asm).solve(np.zeros(ndof))
        return {
            (i, j): T[dom.nodes[nid[(i, j)]].dofs["T"]]
            for j in range(self.NY + 1)
            for i in range(self.NX + 1)
        }

    def _resolver_3d(self):
        mat = ThermalConduction(k=K_ACERO, dim=3)
        dom = Domain()
        nid = {}
        c = 1
        for kk in range(2):
            for j in range(self.NY + 1):
                for i in range(self.NX + 1):
                    dom.add_node(
                        c,
                        [self.L * i / self.NX, self.H * j / self.NY, self.W * kk],
                    )
                    nid[(i, j, kk)] = c
                    c += 1
        eid = 1
        for j in range(self.NY):
            for i in range(self.NX):
                ns = [
                    dom.nodes[nid[(i, j, 0)]], dom.nodes[nid[(i + 1, j, 0)]],
                    dom.nodes[nid[(i + 1, j + 1, 0)]], dom.nodes[nid[(i, j + 1, 0)]],
                    dom.nodes[nid[(i, j, 1)]], dom.nodes[nid[(i + 1, j, 1)]],
                    dom.nodes[nid[(i + 1, j + 1, 1)]], dom.nodes[nid[(i, j + 1, 1)]],
                ]
                dom.add_element(Hex8Thermal(eid, ns, mat))
                eid += 1
        for kk in range(2):
            for j in range(self.NY + 1):
                dom.nodes[nid[(0, j, kk)]].fix_dof("T", self.T1)
                dom.nodes[nid[(self.NX, j, kk)]].fix_dof("T", self.T2)

        asm = Assembler(dom)
        ndof = sum(len(n.dofs) for n in dom.nodes.values())
        T = LinearSolver(asm).solve(np.zeros(ndof))
        return {
            (i, j): T[dom.nodes[nid[(i, j, 0)]].dofs["T"]]
            for j in range(self.NY + 1)
            for i in range(self.NX + 1)
        }

    def test_campos_identicos_nodo_a_nodo(self):
        t2 = self._resolver_2d()
        t3 = self._resolver_3d()
        for key in t2:
            assert t3[key] == pytest.approx(t2[key], abs=1e-10)

    def test_ambos_reproducen_el_analitico(self):
        t3 = self._resolver_3d()
        for j in range(self.NY + 1):
            for i in range(self.NX + 1):
                x = self.L * i / self.NX
                exacto = self.T1 + (self.T2 - self.T1) * x / self.L
                assert t3[(i, j)] == pytest.approx(exacto, abs=1e-10)

    def test_campo_independiente_de_z(self):
        """Caras z adiabáticas ⇒ el problema es efectivamente 2D."""
        mat = ThermalConduction(k=K_ACERO, dim=3)
        dom = Domain()
        nid = {}
        c = 1
        for kk in range(2):
            for j in range(2):
                for i in range(4):
                    dom.add_node(c, [i * 0.5, j * 0.5, kk * 0.5])
                    nid[(i, j, kk)] = c
                    c += 1
        eid = 1
        for i in range(3):
            ns = [
                dom.nodes[nid[(i, 0, 0)]], dom.nodes[nid[(i + 1, 0, 0)]],
                dom.nodes[nid[(i + 1, 1, 0)]], dom.nodes[nid[(i, 1, 0)]],
                dom.nodes[nid[(i, 0, 1)]], dom.nodes[nid[(i + 1, 0, 1)]],
                dom.nodes[nid[(i + 1, 1, 1)]], dom.nodes[nid[(i, 1, 1)]],
            ]
            dom.add_element(Hex8Thermal(eid, ns, mat))
            eid += 1
        for kk in range(2):
            for j in range(2):
                dom.nodes[nid[(0, j, kk)]].fix_dof("T", 80.0)
                dom.nodes[nid[(3, j, kk)]].fix_dof("T", 10.0)

        asm = Assembler(dom)
        ndof = sum(len(n.dofs) for n in dom.nodes.values())
        T = LinearSolver(asm).solve(np.zeros(ndof))

        for i in range(4):
            for j in range(2):
                t0 = T[dom.nodes[nid[(i, j, 0)]].dofs["T"]]
                t1 = T[dom.nodes[nid[(i, j, 1)]].dofs["T"]]
                assert t1 == pytest.approx(t0, abs=1e-10)
