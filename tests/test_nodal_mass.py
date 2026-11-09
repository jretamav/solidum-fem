"""Masas concentradas en nudos (ADR 0024).

Criterios de aceptación (spec de ``ModalSolver``, bloque ``acceptance``):

1. Una masa ``m`` en el ``ux`` de un nudo aumenta ``M_uxux`` exactamente en
   ``m`` y no cambia ninguna otra entrada, con masa consistente y concentrada.
2. Oscilador de un grado de libertad: columna en voladizo sin masa con masa
   ``m`` en la cabeza, ``T = 2π·√(m/k)`` con ``k = 3EI/L³``, diferencia
   relativa ≤ 1e-6. Variante rotacional: inercia ``J`` en el giro de la
   cabeza, ``k_θ = EI/L``.
3. Edificio de cortante de tres pisos (vigas rígidas, columnas sin masa,
   masas de piso en los nudos): periodos y modos iguales a los del problema
   3×3 resuelto a mano, diferencia relativa ≤ 1e-6.

Además: la API y el bloque YAML validan sus datos, la caché de ``M`` ve las
masas añadidas después, la regla del ADR 0008 no cambia, y los consumidores
de ``M`` (espectral, armónico) usan la masa de nudo.

Las columnas no tienen masa propia (``density = 0.0``), así que casi todos
los grados de libertad quedan sin masa: estas pruebas dependen también de la
condensación estática del ``ModalSolver`` (``test_modal_massless_dofs.py``).
Los valores físicos son deliberadamente distintos de 1.
"""
import itertools
import math
import os
import tempfile
import unittest

import numpy as np
import scipy.linalg as sla

from solidum import Assembler, Domain, run_yaml
from solidum.elements.frame import Frame2DEuler
from solidum.elements.truss import Truss2D
from solidum.materials.elastic import Elastic1D
from solidum.math.solvers import HarmonicSolver, ModalSolver, ResponseSpectrumSolver
from solidum.utils.yaml_parser import YamlParser, YamlValidationError

E = 2.0e11          # Pa
L = 3.5             # m
A = 6.0e-3          # m²
I = 8.5e-5          # m⁴
M_TIP = 2.5e4       # kg
J_TIP = 3.2e3       # kg·m²
RTOL_ACCEPT = 1.0e-6


def cantilever(n_el: int, density: float | None = 0.0):
    """Columna vertical en voladizo (empotrada en y = 0), Frame2DEuler."""
    dom = Domain()
    mat = Elastic1D(E=E, density=density) if density is not None else Elastic1D(E=E)
    nodes = [dom.add_node(i + 1, [0.0, i * L / n_el]) for i in range(n_el + 1)]
    for k in range(n_el):
        dom.add_element(Frame2DEuler(k + 1, [nodes[k], nodes[k + 1]], material=mat, A=A, I=I))
    for dof in ("ux", "uy", "rz"):
        nodes[0].fix_dof(dof, 0.0)
    return dom, nodes


def shear_building(heights, masses, inertias, *, span=6.0, n_el=2, area=0.09):
    """Marco de cortante: dos columnas por nivel, sin masa, y vigas rígidas.

    La viga rígida se impone por cinemática, que es la hipótesis exacta del
    edificio de cortante: en cada nivel ``uy = 0`` y ``rz = 0`` en los dos
    nudos (columnas axialmente rígidas, sin giro de nudo) y ``ux`` del nudo
    derecho igual al del izquierdo (restricción lineal). La masa del piso se
    reparte por mitades entre los dos nudos; la del nudo esclavo llega al
    maestro por ``T``. Cada columna es biempotrada: ``k = 12EI/h³``.
    """
    dom = Domain()
    ids, eids = itertools.count(1), itertools.count(1)
    mat = Elastic1D(E=E, density=0.0)

    def node(x, y):
        return dom.add_node(next(ids), [x, y])

    def column(a, b, inertia):
        x, y0, y1 = a.coordinates[0], a.coordinates[1], b.coordinates[1]
        pts = [a] + [node(x, y0 + (y1 - y0) * j / n_el) for j in range(1, n_el)] + [b]
        for p, q in zip(pts[:-1], pts[1:]):
            dom.add_element(Frame2DEuler(next(eids), [p, q], material=mat, A=area, I=inertia))

    base = (node(0.0, 0.0), node(span, 0.0))
    prev, y, floors = base, 0.0, []
    for h, m, inertia in zip(heights, masses, inertias):
        y += h
        left, right = node(0.0, y), node(span, y)
        column(prev[0], left, inertia)
        column(prev[1], right, inertia)
        floors.append((left, right, m))
        prev = (left, right)
    for n in base:
        for dof in ("ux", "uy", "rz"):
            n.fix_dof(dof, 0.0)
    for left, right, m in floors:
        for n in (left, right):
            n.fix_dof("uy", 0.0)
            n.fix_dof("rz", 0.0)
        dom.add_linear_constraint(slave=(right.id, "ux"), masters=[(left.id, "ux")],
                                  coefficients=[1.0])
        dom.add_nodal_mass(left.id, ux=0.5 * m)
        dom.add_nodal_mass(right.id, ux=0.5 * m)
    return dom, [f[0] for f in floors]


def shear_building_3x3(heights, masses, inertias):
    """El problema 3×3 a mano: rigidez de entrepiso 2·12EI/h³."""
    k = [2.0 * 12.0 * E * inertia / h ** 3 for h, inertia in zip(heights, inertias)]
    K = np.array([[k[0] + k[1], -k[1], 0.0],
                  [-k[1], k[1] + k[2], -k[2]],
                  [0.0, -k[2], k[2]]])
    return K, np.diag(masses)


class TestEnsamblajeDeLaMasaDeNudo(unittest.TestCase):
    """Criterio 1 y el contrato del ensamblaje."""

    def test_masa_en_ux_suma_exactamente_m_y_nada_mas(self):
        for lumping in ("consistent", "lumped"):
            with self.subTest(lumping=lumping):
                dom, nodes = cantilever(3, density=7850.0)
                dom.generate_equation_numbers()
                M0 = Assembler(dom).assemble_mass_matrix(lumping).toarray()
                dom.add_nodal_mass(nodes[-1].id, ux=M_TIP)
                M1 = Assembler(dom).assemble_mass_matrix(lumping).toarray()
                i = nodes[-1].dofs["ux"]
                expected = M0.copy()
                expected[i, i] = M0[i, i] + M_TIP
                np.testing.assert_array_equal(M1, expected)

    def test_inercia_rotacional_en_rz(self):
        dom, nodes = cantilever(3, density=7850.0)
        dom.generate_equation_numbers()
        M0 = Assembler(dom).assemble_mass_matrix().toarray()
        dom.add_nodal_mass(nodes[2].id, rz=J_TIP)
        M1 = Assembler(dom).assemble_mass_matrix().toarray()
        i = nodes[2].dofs["rz"]
        expected = M0.copy()
        expected[i, i] = M0[i, i] + J_TIP
        np.testing.assert_array_equal(M1, expected)

    def test_masa_conserva_el_patron_de_K(self):
        dom, nodes = cantilever(3, density=0.0)
        dom.add_nodal_mass(nodes[-1].id, ux=M_TIP, uy=M_TIP)
        asm = Assembler(dom)
        asm.assemble_system()
        M = asm.assemble_mass_matrix()
        np.testing.assert_array_equal(M.indptr, asm.K_global.indptr)
        np.testing.assert_array_equal(M.indices, asm.K_global.indices)

    def test_varias_declaraciones_se_suman(self):
        dom, nodes = cantilever(2)
        dom.add_nodal_mass(nodes[-1].id, ux=1.2e4)
        dom.add_nodal_mass(nodes[-1].id, ux=0.3e4, uy=0.7e4)
        self.assertEqual(dom.nodal_masses[(nodes[-1].id, "ux")], 1.5e4)
        self.assertEqual(dom.nodal_masses[(nodes[-1].id, "uy")], 0.7e4)

    def test_masa_anadida_despues_del_primer_ensamblaje_no_se_ignora(self):
        dom, nodes = cantilever(2, density=7850.0)
        asm = Assembler(dom)
        M_a = asm.assemble_mass_matrix()
        M_a_copy = M_a.toarray()
        dom.add_nodal_mass(nodes[-1].id, ux=M_TIP)
        M_b = asm.assemble_mass_matrix()
        i = nodes[-1].dofs["ux"]
        self.assertEqual(M_b[i, i], M_a_copy[i, i] + M_TIP)
        np.testing.assert_array_equal(M_a.toarray(), M_a_copy)   # la anterior no cambia
        self.assertIs(asm.assemble_mass_matrix(), M_b)          # y la nueva se cachea

    def test_material_sin_density_sigue_siendo_error(self):
        """ADR 0008: la masa de nudo no exime de declarar la densidad."""
        dom, nodes = cantilever(2, density=None)
        dom.add_nodal_mass(nodes[-1].id, ux=M_TIP)
        with self.assertRaisesRegex(ValueError, "density"):
            Assembler(dom).assemble_mass_matrix()

    def test_grado_de_libertad_que_el_nodo_no_tiene(self):
        dom = Domain()
        n1, n2 = dom.add_node(1, [0.0, 0.0]), dom.add_node(2, [2.0, 0.0])
        dom.add_element(Truss2D(1, [n1, n2], Elastic1D(E=E, density=7850.0), A=A))
        dom.add_nodal_mass(2, rz=J_TIP)
        with self.assertRaisesRegex(ValueError, r"\(nodo 2, 'rz'\)"):
            Assembler(dom).assemble_mass_matrix()


class TestValidacionDeLaApi(unittest.TestCase):

    def setUp(self):
        self.dom, self.nodes = cantilever(2)
        self.tip = self.nodes[-1].id

    def test_nodo_inexistente(self):
        with self.assertRaisesRegex(ValueError, "ningún nodo 99"):
            self.dom.add_nodal_mass(99, ux=M_TIP)

    def test_sin_grados_de_libertad(self):
        with self.assertRaisesRegex(ValueError, "al menos un grado de libertad"):
            self.dom.add_nodal_mass(self.tip)

    def test_grado_de_libertad_no_mecanico(self):
        with self.assertRaisesRegex(ValueError, "capacidad calorífica"):
            self.dom.add_nodal_mass(self.tip, T=1.0)

    def test_valores_invalidos_y_nada_se_guarda(self):
        for bad in (-1.0, math.nan, math.inf, True, "1000"):
            with self.subTest(value=bad):
                with self.assertRaisesRegex(ValueError, "finito no negativo"):
                    self.dom.add_nodal_mass(self.tip, ux=M_TIP, uy=bad)
        self.assertEqual(self.dom.nodal_masses, {})

    def test_acepta_tipos_numericos_de_numpy(self):
        self.dom.add_nodal_mass(self.tip, ux=np.float64(M_TIP), uy=np.int64(7))
        self.assertEqual(self.dom.nodal_masses[(self.tip, "uy")], 7.0)


class TestOsciladorDeUnGradoDeLibertad(unittest.TestCase):
    """Criterio 2: columna sin masa con masa concentrada en la cabeza."""

    def test_periodo_traslacional(self):
        k = 3.0 * E * I / L ** 3
        T_exact = 2.0 * math.pi * math.sqrt(M_TIP / k)
        for n_el in (1, 4):
            for lumping in ("consistent", "lumped"):
                with self.subTest(n_el=n_el, lumping=lumping):
                    dom, nodes = cantilever(n_el)
                    dom.add_nodal_mass(nodes[-1].id, ux=M_TIP)
                    res = ModalSolver(Assembler(dom), n_modes=1, lumping=lumping).solve()
                    self.assertLessEqual(abs(res.periods[0] / T_exact - 1.0), RTOL_ACCEPT)

    def test_forma_modal_en_los_grados_de_libertad_sin_masa(self):
        """Los grados de libertad sin masa se recuperan de forma exacta: el
        modo es la elástica de la carga en la punta,
        ``v(y)/v(L) = y²(3L − y)/(2L³)`` y ``θ(L) = −3/(2L)·v(L)`` (giro
        horario de una columna que se inclina hacia +x), sin ``uy``."""
        dom, nodes = cantilever(4)
        dom.add_nodal_mass(nodes[-1].id, ux=M_TIP)
        res = ModalSolver(Assembler(dom), n_modes=1).solve()
        phi = res.modes[:, 0]
        tip = phi[nodes[-1].dofs["ux"]]
        self.assertAlmostEqual(abs(tip), 1.0 / math.sqrt(M_TIP), delta=1e-12 / math.sqrt(M_TIP))
        for n in nodes:
            y = n.coordinates[1]
            self.assertAlmostEqual(phi[n.dofs["ux"]] / tip, y ** 2 * (3 * L - y) / (2 * L ** 3),
                                   delta=1e-9)
            self.assertAlmostEqual(phi[n.dofs["uy"]] / tip, 0.0, delta=1e-9)
        self.assertAlmostEqual(phi[nodes[-1].dofs["rz"]] / tip, -1.5 / L, delta=1e-9)

    def test_periodo_rotacional(self):
        """Inercia J en el giro de la cabeza, ux libre y sin masa: la rigidez
        al giro con la punta libre de trasladarse es ``EI/L``."""
        T_exact = 2.0 * math.pi * math.sqrt(J_TIP / (E * I / L))
        dom, nodes = cantilever(4)
        dom.add_nodal_mass(nodes[-1].id, rz=J_TIP)
        res = ModalSolver(Assembler(dom), n_modes=1).solve()
        self.assertLessEqual(abs(res.periods[0] / T_exact - 1.0), RTOL_ACCEPT)


class TestEdificioDeCortante(unittest.TestCase):
    """Criterio 3: tres pisos contra el problema 3×3 a mano."""

    def test_pisos_iguales_contra_la_forma_cerrada(self):
        """Cadena empotrada-libre de tres masas y resortes iguales:
        ``ω_j = 2·√(k/m)·sin((2j−1)π/14)`` y ``φ_j(i) ∝ sin(i(2j−1)π/7)``."""
        h, m, inertia = 3.2, 4.0e4, 3.0e-4
        k = 2.0 * 12.0 * E * inertia / h ** 3
        for lumping in ("consistent", "lumped"):
            with self.subTest(lumping=lumping):
                dom, floors = shear_building([h] * 3, [m] * 3, [inertia] * 3)
                res = ModalSolver(Assembler(dom), n_modes=3, lumping=lumping).solve()
                for j in range(1, 4):
                    omega = 2.0 * math.sqrt(k / m) * math.sin((2 * j - 1) * math.pi / 14)
                    self.assertLessEqual(abs(res.frequencies_rad[j - 1] / omega - 1.0), RTOL_ACCEPT)
                    shape = np.array([res.modes[f.dofs["ux"], j - 1] for f in floors])
                    exact = np.array([math.sin(i * (2 * j - 1) * math.pi / 7) for i in (1, 2, 3)])
                    np.testing.assert_allclose(shape / shape[-1], exact / exact[-1],
                                               rtol=0.0, atol=RTOL_ACCEPT)

    def test_pisos_distintos_contra_el_problema_3x3(self):
        heights, masses = [4.0, 3.2, 3.0], [5.2e4, 4.6e4, 3.1e4]
        inertias = [4.2e-4, 3.0e-4, 2.2e-4]
        K3, M3 = shear_building_3x3(heights, masses, inertias)
        lam, vec = sla.eigh(K3, M3)                 # M3-normalizados: Σ mᵢφᵢ² = 1
        dom, floors = shear_building(heights, masses, inertias)
        res = ModalSolver(Assembler(dom), n_modes=3).solve()
        np.testing.assert_allclose(res.periods, 2.0 * math.pi / np.sqrt(lam), rtol=RTOL_ACCEPT)
        for j in range(3):
            shape = np.array([res.modes[f.dofs["ux"], j] for f in floors])
            ref = vec[:, j] * np.sign(vec[-1, j] * shape[-1])
            # Misma normalización: la masa del piso, repartida entre los dos
            # nudos de un nivel que se mueven juntos, suma mᵢ.
            np.testing.assert_allclose(shape, ref, rtol=0.0,
                                       atol=RTOL_ACCEPT * np.max(np.abs(ref)))


class TestConsumidoresDeLaMatrizDeMasa(unittest.TestCase):
    """El espectral y el armónico usan la masa de nudo sin código propio."""

    def test_espectral_masa_total_y_respuesta_de_azotea(self):
        heights, masses = [4.0, 3.2, 3.0], [5.2e4, 4.6e4, 3.1e4]
        inertias = [4.2e-4, 3.0e-4, 2.2e-4]
        Sa = 2.45                                    # m/s²
        dom, floors = shear_building(heights, masses, inertias)
        res = ResponseSpectrumSolver(
            Assembler(dom), n_modes=3, direction={"dof_name": "ux"},
            spectrum={"type": "constant_acceleration", "Sa": Sa},
        ).solve()
        self.assertAlmostEqual(res.total_mass / sum(masses), 1.0, delta=1e-12)
        self.assertAlmostEqual(np.sum(res.effective_masses) / sum(masses), 1.0, delta=1e-10)
        # Azotea por SRSS, a mano: uₙ = Γₙ·φₙ·Sa/ωₙ², Γₙ = φₙᵀ·M·1.
        K3, M3 = shear_building_3x3(heights, masses, inertias)
        lam, vec = sla.eigh(K3, M3)
        gamma = vec.T @ M3 @ np.ones(3)
        roof = math.sqrt(sum((gamma[n] * vec[2, n] * Sa / lam[n]) ** 2 for n in range(3)))
        self.assertLessEqual(abs(res.u_combined[floors[-1].dofs["ux"]] / roof - 1.0), RTOL_ACCEPT)

    def test_armonico_oscilador_de_un_grado_de_libertad(self):
        """Respuesta estacionaria sin amortiguamiento: ``u = F/(k − Ω²m)``."""
        k = 3.0 * E * I / L ** 3
        Omega, F0 = 0.6 * math.sqrt(k / M_TIP), 1.5e4
        dom, nodes = cantilever(4)
        dom.add_nodal_mass(nodes[-1].id, ux=M_TIP)
        dom.generate_equation_numbers()
        F = np.zeros(dom.total_dofs)
        F[nodes[-1].dofs["ux"]] = F0
        res = HarmonicSolver(Assembler(dom), omega=np.array([Omega]), F_amplitude=F).solve()
        u = res.u_complex[nodes[-1].dofs["ux"], 0]
        self.assertAlmostEqual(u.real / (F0 / (k - Omega ** 2 * M_TIP)), 1.0, delta=1e-9)
        self.assertAlmostEqual(u.imag, 0.0, delta=1e-15)


def _cantilever_yaml(n_el: int, extra: str) -> str:
    nodes = "\n".join(f"  - {{id: {i + 1}, coords: [0.0, {i * L / n_el!r}]}}" for i in range(n_el + 1))
    elems = "\n".join(
        f"  - {{id: {k + 1}, type: Frame2DEuler, material: 1, nodes: [{k + 1}, {k + 2}], "
        f"A: {A!r}, I: {I!r}}}" for k in range(n_el))
    return (f"nodes:\n{nodes}\nmaterials:\n  - {{id: 1, type: Elastic1D, E: {E!r}, density: 0.0}}\n"
            f"elements:\n{elems}\nboundary_conditions:\n  - {{node_id: 1, ux: 0.0, uy: 0.0, rz: 0.0}}\n"
            f"{extra}")


class TestYaml(unittest.TestCase):

    def _write(self, text: str) -> str:
        fd, path = tempfile.mkstemp(suffix=".yaml")
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            f.write(text)
        self.addCleanup(os.remove, path)
        return path

    def test_oscilador_desde_yaml(self):
        path = self._write(_cantilever_yaml(3, (
            f"nodal_masses:\n  - {{node_id: 4, ux: {M_TIP!r}}}\n"
            "solver:\n  type: ModalSolver\n  n_modes: 1\n")))
        res = run_yaml(path)
        T_exact = 2.0 * math.pi * math.sqrt(M_TIP / (3.0 * E * I / L ** 3))
        self.assertLessEqual(abs(res.periods[0] / T_exact - 1.0), RTOL_ACCEPT)

    def test_el_bloque_llega_al_dominio(self):
        path = self._write(_cantilever_yaml(3, (
            "nodal_masses:\n"
            "  - {node_id: 4, ux: 1.0e4, uy: 1.0e4}\n"
            "  - {node_id: 4, ux: 2.0e3, rz: 5.0e2}\n")))
        dom = YamlParser(path).parse()
        self.assertEqual(dom.nodal_masses, {(4, "ux"): 1.2e4, (4, "uy"): 1.0e4, (4, "rz"): 5.0e2})

    def test_errores_de_forma_juntos(self):
        path = self._write(_cantilever_yaml(3, (
            "nodal_masses:\n"
            "  - {node_id: 4, ux: -1.0}\n"
            "  - {node_id: 4, T: 3.0}\n"
            "  - {ux: 1.0}\n"
            "  - {node_id: 77, ux: 1.0}\n"
            "  - {node_id: 4}\n")))
        with self.assertRaises(YamlValidationError) as ctx:
            YamlParser(path).parse()
        text = str(ctx.exception)
        for fragment in ("nodal_masses[0].ux", "nodal_masses[1]: 'T'", "nodal_masses[2]: falta",
                         "nodal_masses[3]: 'node_id=77'", "nodal_masses[4]: no declara"):
            self.assertIn(fragment, text)

    def test_grado_de_libertad_ausente_en_el_nodo(self):
        text = "\n".join([
            "nodes:", "  - {id: 1, coords: [0.0, 0.0]}", "  - {id: 2, coords: [2.0, 0.0]}",
            "materials:", "  - {id: 1, type: Elastic1D, E: 2.0e11, density: 7850.0}",
            "elements:", "  - {id: 1, type: Truss2D, material: 1, A: 0.01, nodes: [1, 2]}",
            "nodal_masses:", "  - {node_id: 2, ux: 10.0, rz: 3.0}", ""])
        with self.assertRaisesRegex(YamlValidationError, "no tiene el grado de libertad 'rz'"):
            YamlParser(self._write(text)).parse()

    def test_lista_requerida(self):
        path = self._write(_cantilever_yaml(1, "nodal_masses: {node_id: 2, ux: 1.0}\n"))
        with self.assertRaisesRegex(YamlValidationError, "debe ser una lista"):
            YamlParser(path).parse()


if __name__ == "__main__":
    unittest.main()
