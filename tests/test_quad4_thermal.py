"""Tests del elemento ``Quad4Thermal`` (Etapa 8).

Cubre el bloque ``acceptance`` de ``docs/specs/Quad4Thermal.md``.

Coeficientes físicos no unitarios a propósito (acero: k=45, c=460, ρ=7850;
espesor 0.3 m): un factor perdido o un error dimensional quedaría
enmascarado si todo valiera 1.
"""
import numpy as np
import pytest

from solidum.core.domain import Domain
from solidum.core.node import Node
from solidum.elements.thermal.quad4_thermal import Quad4Thermal
from solidum.materials.elastic_2d import Elastic2D
from solidum.materials.thermal_conduction import ThermalConduction
from solidum.math.assembly import Assembler
from solidum.math.solvers.linear import LinearSolver

K_ACERO = 45.0
C_ACERO = 460.0
RHO_ACERO = 7850.0
THK = 0.3


def _material(dim=2, **kw):
    return ThermalConduction(k=K_ACERO, dim=dim, **kw)


def _unit_square(thickness=THK, material=None, quadrature=None):
    nodes = [
        Node(1, [0.0, 0.0]),
        Node(2, [1.0, 0.0]),
        Node(3, [1.0, 1.0]),
        Node(4, [0.0, 1.0]),
    ]
    mat = material if material is not None else _material()
    return Quad4Thermal(1, nodes, mat, thickness=thickness, quadrature=quadrature)


def _distorted(thickness=THK):
    nodes = [
        Node(1, [0.0, 0.0]),
        Node(2, [2.3, 0.15]),
        Node(3, [2.05, 1.7]),
        Node(4, [0.2, 1.35]),
    ]
    return Quad4Thermal(1, nodes, _material(), thickness=thickness)


class TestContrato:
    def test_declara_un_dof_escalar(self):
        el = _unit_square()
        assert el.DOF_NAMES == ["T"]
        assert el.FLUX_DIM == 2
        assert el.N_NODES == 4
        assert el.N_INTEGRATION_POINTS == 4

    def test_registrado_como_elemento(self):
        from solidum.registry import ElementRegistry

        assert "Quad4Thermal" in ElementRegistry.names()

    def test_es_clase_distinta_del_quad4_mecanico(self):
        from solidum.elements.solid_2d.quad4 import Quad4

        assert Quad4Thermal is not Quad4
        assert Quad4.DOF_NAMES == ["ux", "uy"]
        assert Quad4Thermal.DOF_NAMES == ["T"]

    def test_sin_estado_interno(self):
        """La conducción de Fourier no tiene variables internas."""
        el = _unit_square()
        assert el.state is None
        el.commit_state()  # no debe fallar


class TestMatrizDeConductividad:
    def test_simetria_exacta(self):
        for el in (_unit_square(), _distorted()):
            Ke = el.compute_conductivity_matrix()
            np.testing.assert_allclose(Ke, Ke.T, atol=1e-12)

    def test_modo_nulo_temperatura_uniforme(self):
        """Un campo uniforme no produce gradiente ⇒ tampoco flujo.

        Análogo térmico de los modos de sólido rígido: por eso el sistema
        global sólo se vuelve resoluble al imponer un Dirichlet.
        """
        for el in (_unit_square(), _distorted()):
            Ke = el.compute_conductivity_matrix()
            np.testing.assert_allclose(Ke @ np.ones(4), np.zeros(4), atol=1e-10)

    def test_rango_y_autovalores(self):
        el = _unit_square()
        Ke = el.compute_conductivity_matrix()
        eig = np.linalg.eigvalsh(Ke)

        assert np.min(eig) > -1e-10                     # semidefinida positiva
        assert np.sum(np.abs(eig) < 1e-9) == 1          # exactamente un cero
        assert np.linalg.matrix_rank(Ke, tol=1e-9) == 3

    def test_escala_con_la_conductividad(self):
        """K_e es lineal en k: doblar k dobla la matriz."""
        el1 = _unit_square(material=ThermalConduction(k=K_ACERO, dim=2))
        el2 = _unit_square(material=ThermalConduction(k=2 * K_ACERO, dim=2))
        np.testing.assert_allclose(
            el2.compute_conductivity_matrix(),
            2.0 * el1.compute_conductivity_matrix(),
            rtol=1e-13,
        )

    def test_escala_con_el_espesor(self):
        e1 = _unit_square(thickness=0.3)
        e2 = _unit_square(thickness=0.9)
        np.testing.assert_allclose(
            e2.compute_conductivity_matrix(),
            3.0 * e1.compute_conductivity_matrix(),
            rtol=1e-13,
        )

    def test_hourglass_con_integracion_reducida(self):
        """Documenta la limitación de la cuadratura 1×1, no la corrige."""
        el = _unit_square(quadrature="1x1")
        Ke = el.compute_conductivity_matrix()
        assert np.linalg.matrix_rank(Ke, tol=1e-9) == 2  # frente a los 3 debidos


class TestPatchTest:
    def test_gradiente_de_campo_lineal_exacto(self):
        """T = a + b·x + c·y ⇒ ∇T = (b, c) en TODOS los Gauss."""
        a, b, c = 2.0, 3.0, 5.0
        for el in (_unit_square(), _distorted()):
            coords = np.array([n.coordinates[:2] for n in el.nodes])
            T_e = a + b * coords[:, 0] + c * coords[:, 1]
            for p in el.points:
                B, _ = el._gradient(p, coords)
                np.testing.assert_allclose(B @ T_e, [b, c], rtol=1e-11)

    def test_flujo_de_campo_lineal(self):
        """q = -k·∇T con el gradiente del campo lineal impuesto."""
        el = _unit_square()
        coords = np.array([n.coordinates[:2] for n in el.nodes])
        T_e = 3.0 * coords[:, 0] + 5.0 * coords[:, 1]

        for idx, n in enumerate(el.nodes):
            n.dofs["T"] = idx
        gs = el.compute_gauss_state(T_e)

        for g in range(len(el.points)):
            np.testing.assert_allclose(gs["grad_T"][g], [3.0, 5.0], rtol=1e-11)
            np.testing.assert_allclose(
                gs["flux"][g], [-K_ACERO * 3.0, -K_ACERO * 5.0], rtol=1e-11
            )


class TestMatrizDeCapacidad:
    def _mat_con_capacidad(self):
        return ThermalConduction(
            k=K_ACERO, c=C_ACERO, density=RHO_ACERO, dim=2
        )

    def test_capacidad_total_conservada_en_ambos_modos(self):
        """Σ C = ρ·c·A·t: la energía almacenable no depende del reparto."""
        el = _unit_square(material=self._mat_con_capacidad())
        esperado = RHO_ACERO * C_ACERO * 1.0 * THK

        for modo in ("consistent", "lumped"):
            C = el.compute_capacity_matrix(modo)
            assert C.sum() == pytest.approx(esperado, rel=1e-12)

    def test_lumped_es_diagonal_y_positiva(self):
        el = _unit_square(material=self._mat_con_capacidad())
        C = el.compute_capacity_matrix("lumped")
        np.testing.assert_allclose(C, np.diag(np.diag(C)), atol=1e-12)
        assert np.all(np.diag(C) > 0.0)

    def test_consistente_es_simetrica_y_definida_positiva(self):
        el = _unit_square(material=self._mat_con_capacidad())
        C = el.compute_capacity_matrix("consistent")
        np.testing.assert_allclose(C, C.T, atol=1e-12)
        assert np.min(np.linalg.eigvalsh(C)) > 0.0

    def test_default_es_lumped(self):
        """Al revés que en dinámica estructural; razón física en la spec."""
        el = _unit_square(material=self._mat_con_capacidad())
        np.testing.assert_allclose(
            el.compute_capacity_matrix(),
            el.compute_capacity_matrix("lumped"),
            rtol=1e-14,
        )

    def test_capacidad_sin_c_ni_density_falla_accionable(self):
        el = _unit_square()  # material sin c ni density
        with pytest.raises(ValueError) as exc:
            el.compute_capacity_matrix()
        msg = str(exc.value)
        assert "capacidad" in msg
        assert "LinearSolver" in msg or "estacionario" in msg

    def test_conductividad_no_exige_capacidad(self):
        """El camino estacionario no consulta c ni density."""
        el = _unit_square()
        Ke = el.compute_conductivity_matrix()
        assert np.isfinite(Ke).all()


class TestCargas:
    def test_fuente_volumetrica_uniforme(self):
        """Σf = Q·A·t exacto, invariante ante distorsión."""
        Q = 1000.0
        el = _unit_square()
        assert el.compute_body_source(Q).sum() == pytest.approx(
            Q * 1.0 * THK, rel=1e-12
        )

        dis = _distorted()
        # Área del cuadrilátero distorsionado por la fórmula del cordón.
        c = np.array([n.coordinates[:2] for n in dis.nodes])
        area = 0.5 * abs(
            np.dot(c[:, 0], np.roll(c[:, 1], -1))
            - np.dot(c[:, 1], np.roll(c[:, 0], -1))
        )
        assert dis.compute_body_source(Q).sum() == pytest.approx(
            Q * area * THK, rel=1e-12
        )

    def test_flujo_en_borde_reparto_mitad(self):
        q_bar = 100.0
        el = _unit_square()
        f = el.compute_edge_flux(0, q_bar)  # borde n0–n1, longitud 1

        esperado = -0.5 * 1.0 * q_bar * THK
        assert f[0] == pytest.approx(esperado, rel=1e-13)
        assert f[1] == pytest.approx(esperado, rel=1e-13)
        assert f[2] == 0.0 and f[3] == 0.0
        assert f.sum() == pytest.approx(-q_bar * 1.0 * THK, rel=1e-13)

    def test_signo_saliente_extrae_energia(self):
        """q̄ > 0 es flujo saliente ⇒ contribución negativa (enfría)."""
        el = _unit_square()
        assert el.compute_edge_flux(0, 100.0).sum() < 0.0
        assert el.compute_edge_flux(0, -100.0).sum() > 0.0

    def test_todos_los_bordes(self):
        el = _unit_square()
        for edge in range(4):
            f = el.compute_edge_flux(edge, 50.0)
            a, c = el.EDGE_NODES[edge]
            assert f[a] != 0.0 and f[c] != 0.0
            otros = [i for i in range(4) if i not in (a, c)]
            assert all(f[i] == 0.0 for i in otros)

    def test_borde_invalido_rechazado(self):
        el = _unit_square()
        with pytest.raises(ValueError, match="edge=4"):
            el.compute_edge_flux(4, 1.0)


class TestValidaciones:
    def test_rechaza_material_mecanico(self):
        """Familias paralelas: un material mecánico no cuela."""
        nodes = [Node(i + 1, c) for i, c in enumerate(
            [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
        )]
        with pytest.raises(TypeError, match="familia térmica"):
            Quad4Thermal(1, nodes, Elastic2D(E=210e9, nu=0.3))

    def test_rechaza_material_3d(self):
        nodes = [Node(i + 1, c) for i, c in enumerate(
            [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
        )]
        with pytest.raises(ValueError, match="FLUX_DIM"):
            Quad4Thermal(1, nodes, ThermalConduction(k=K_ACERO, dim=3))

    def test_rechaza_numero_de_nodos_incorrecto(self):
        nodes = [Node(i + 1, c) for i, c in enumerate(
            [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]
        )]
        with pytest.raises(ValueError, match="requiere 4 nodos"):
            Quad4Thermal(1, nodes, _material())

    def test_rechaza_thickness_no_positivo(self):
        nodes = [Node(i + 1, c) for i, c in enumerate(
            [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
        )]
        with pytest.raises(ValueError, match="thickness"):
            Quad4Thermal(1, nodes, _material(), thickness=0.0)

    def test_jacobiano_degenerado_abortado(self):
        """Nodos en orden horario ⇒ det J < 0."""
        nodes = [
            Node(1, [0.0, 0.0]),
            Node(2, [0.0, 1.0]),
            Node(3, [1.0, 1.0]),
            Node(4, [1.0, 0.0]),
        ]
        el = Quad4Thermal(1, nodes, _material())
        with pytest.raises(ValueError, match="Jacobiano"):
            el.compute_conductivity_matrix()


class TestParedPlana:
    """Benchmark con solución analítica cerrada: perfil lineal exacto."""

    @staticmethod
    def _resolver(nx=8, ny=3, L=2.0, H=0.5, T1=100.0, T2=20.0):
        mat = ThermalConduction(k=K_ACERO, dim=2)
        dom = Domain()
        nid = {}
        c = 1
        for j in range(ny + 1):
            for i in range(nx + 1):
                dom.add_node(c, [L * i / nx, H * j / ny])
                nid[(i, j)] = c
                c += 1

        eid = 1
        for j in range(ny):
            for i in range(nx):
                ns = [
                    dom.nodes[nid[(i, j)]],
                    dom.nodes[nid[(i + 1, j)]],
                    dom.nodes[nid[(i + 1, j + 1)]],
                    dom.nodes[nid[(i, j + 1)]],
                ]
                dom.add_element(Quad4Thermal(eid, ns, mat, thickness=1.0))
                eid += 1

        for j in range(ny + 1):
            dom.nodes[nid[(0, j)]].fix_dof("T", T1)
            dom.nodes[nid[(nx, j)]].fix_dof("T", T2)

        asm = Assembler(dom)
        ndof = sum(len(n.dofs) for n in dom.nodes.values())
        T = LinearSolver(asm).solve(np.zeros(ndof))
        return dom, nid, T

    def test_perfil_lineal_exacto(self):
        nx, ny, L, T1, T2 = 8, 3, 2.0, 100.0, 20.0
        dom, nid, T = self._resolver(nx=nx, ny=ny, L=L, T1=T1, T2=T2)

        for j in range(ny + 1):
            for i in range(nx + 1):
                x = L * i / nx
                exacto = T1 + (T2 - T1) * x / L
                got = T[dom.nodes[nid[(i, j)]].dofs["T"]]
                assert got == pytest.approx(exacto, abs=1e-10)

    def test_campo_independiente_de_y(self):
        """Bordes superior/inferior adiabáticos ⇒ conducción 1D pura."""
        nx, ny = 8, 3
        dom, nid, T = self._resolver(nx=nx, ny=ny)
        for i in range(nx + 1):
            col = [T[dom.nodes[nid[(i, j)]].dofs["T"]] for j in range(ny + 1)]
            assert np.ptp(col) < 1e-10

    def test_flujo_uniforme_en_los_gauss(self):
        """q_x = k·(T1−T2)/L constante; q_y nulo."""
        nx, L, T1, T2 = 8, 2.0, 100.0, 20.0
        dom, _, T = self._resolver(nx=nx, L=L, T1=T1, T2=T2)
        q_esperado = -K_ACERO * (T2 - T1) / L

        for el in dom.elements.values():
            gs = el.compute_gauss_state(T)
            for g in range(len(el.points)):
                assert gs["flux"][g, 0] == pytest.approx(q_esperado, rel=1e-9)
                assert gs["flux"][g, 1] == pytest.approx(0.0, abs=1e-8)
