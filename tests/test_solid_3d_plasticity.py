"""Tests de integración del pipeline Hex8/Tet4 + VonMises3D + NonlinearSolver.

Validan el flujo completo con sólidos 3D y plasticidad J2 sobre Voigt 6D
del proyecto (ADR 0012). Cierra la cobertura de integración 3D + material
no lineal, equivalente al test de integración Quad4 + VonMises2D del 2D.

Acceptance de ``docs/specs/VonMises3D.md`` sección ``acceptance.integration``.
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from solidum.core.domain import Domain
from solidum.elements.solid_3d import Hex8, Tet4
from solidum.materials.damage_3d import IsotropicDamage3D
from solidum.materials.drucker_prager_3d import DruckerPrager3D
from solidum.materials.von_mises_3d import VonMises3D
from solidum.math.assembly import Assembler
from solidum.math.convergence import ConvergenceCriterion
from solidum.math.solvers import NonlinearSolver


HEX8_UNIT_COORDS = [
    (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
]
# Nodos en cada cara del cubo unitario Hex8 (0-indexados en la lista de nodos)
HEX8_FACE_XMIN = (0, 3, 4, 7)     # x = 0
HEX8_FACE_XMAX = (1, 2, 5, 6)     # x = 1 (cara +ξ del elemento)
HEX8_FACE_YMIN = (0, 1, 4, 5)     # y = 0
HEX8_FACE_ZMIN = (0, 1, 2, 3)     # z = 0


def _build_hex8_unit_cube(material):
    """Cubo unitario Hex8 con coordenadas estándar."""
    dom = Domain()
    nodes = [dom.add_node(i + 1, list(c)) for i, c in enumerate(HEX8_UNIT_COORDS)]
    elem = Hex8(1, nodes, material)
    dom.add_element(elem)
    return dom, elem, nodes


class TestHex8VonMises3DUniaxialFree(unittest.TestCase):
    """Hex8 + VonMises3D bajo tracción uniaxial libre con Poisson activo.

    Setup: cubo unitario con tres planos de simetría (x=0, y=0, z=0)
    como rodillos (BC tipo roller). Las caras x=1, y=1, z=1 quedan
    libres salvo por la tracción aplicada en x=1. Bajo esta BC y carga
    horizontal, el estado tensorial es σ = diag(σ_xx, 0, 0) — uniaxial
    puro con contracción lateral por Poisson (elástico) o flujo plástico
    (post-yield).
    """

    def setUp(self):
        self.E = 2.0e5
        self.nu = 0.3
        self.sigma_y = 250.0
        self.H = 1.0e4
        self.mat = VonMises3D(
            E=self.E, nu=self.nu, sigma_y=self.sigma_y, H=self.H,
        )
        self.dom, self.elem, self.nodes = _build_hex8_unit_cube(self.mat)

        # Rodillos en los tres planos de simetría
        for i in HEX8_FACE_XMIN:
            self.nodes[i].fix_dof('ux', 0.0)
        for i in HEX8_FACE_YMIN:
            self.nodes[i].fix_dof('uy', 0.0)
        for i in HEX8_FACE_ZMIN:
            self.nodes[i].fix_dof('uz', 0.0)

        self.dom.generate_equation_numbers(verbose=False)
        self.assembler = Assembler(self.dom)

    def _apply_traction_xmax(self, F_total: float) -> np.ndarray:
        """Reparte ``F_total`` uniformemente en los 4 nodos de la cara x=1."""
        F_ext = np.zeros(self.dom.total_dofs)
        for i in HEX8_FACE_XMAX:
            F_ext[self.nodes[i].dofs['ux']] = F_total / 4.0
        return F_ext

    def test_elastic_uniaxial_response(self):
        """Carga 80% del yield: σ_xx = F/A, σ_yy ≈ σ_zz ≈ 0, ε_xx = σ/E."""
        sigma_target = 0.8 * self.sigma_y
        F_ext = self._apply_traction_xmax(sigma_target)  # área = 1

        conv = ConvergenceCriterion(rtol_force=1.0e-9, rtol_disp=1.0e-9)
        solver = NonlinearSolver(self.assembler, convergence=conv, num_steps=2)
        U = solver.solve(F_ext)

        # u_x del nodo 1 = ε_xx (cubo de longitud 1)
        u_x = U[self.nodes[1].dofs['ux']]
        self.assertAlmostEqual(u_x, sigma_target / self.E, places=7)

        # Estado en puntos de Gauss
        gs = self.elem.compute_gauss_state(U)
        sigma_xx_avg = gs['stress'][:, 0].mean()
        sigma_yy_avg = gs['stress'][:, 1].mean()
        sigma_zz_avg = gs['stress'][:, 2].mean()
        self.assertAlmostEqual(sigma_xx_avg, sigma_target, delta=sigma_target * 1.0e-6)
        self.assertAlmostEqual(sigma_yy_avg, 0.0, delta=sigma_target * 1.0e-6)
        self.assertAlmostEqual(sigma_zz_avg, 0.0, delta=sigma_target * 1.0e-6)
        # α = 0 en régimen elástico
        for state in self.elem.state.vars:
            self.assertAlmostEqual(state['alpha'], 0.0, places=12)

    def test_plastic_uniaxial_yield_at_sigma_y(self):
        """Carga 30% sobre el yield uniaxial: σ_xx ≈ σ_y, plasticidad activa."""
        # En uniaxial puro 3D el yield es σ_xx = σ_y.
        sigma_target = 1.3 * self.sigma_y
        F_ext = self._apply_traction_xmax(sigma_target)

        conv = ConvergenceCriterion(rtol_force=1.0e-7, rtol_disp=1.0e-7)
        solver = NonlinearSolver(self.assembler, convergence=conv, num_steps=10)
        U = solver.solve(F_ext)

        u_x = U[self.nodes[1].dofs['ux']]
        # u_x > σ_y/E porque entramos en régimen plástico
        self.assertGreater(u_x, self.sigma_y / self.E)

        gs = self.elem.compute_gauss_state(U)
        sigma_xx_avg = gs['stress'][:, 0].mean()
        # σ_xx tras return mapping ≈ σ_target (equilibrio con F externa)
        self.assertAlmostEqual(sigma_xx_avg, sigma_target, delta=sigma_target * 1.0e-3)
        # σ_yy, σ_zz se mantienen ≈ 0 (estado uniaxial preservado por las BCs)
        self.assertAlmostEqual(gs['stress'][:, 1].mean(), 0.0, delta=sigma_target * 1.0e-3)
        self.assertAlmostEqual(gs['stress'][:, 2].mean(), 0.0, delta=sigma_target * 1.0e-3)

        # α > 0 en todos los puntos de Gauss
        for state in self.elem.state.vars:
            self.assertGreater(state['alpha'], 0.0)
            # Incompresibilidad plástica preservada
            eps_p = state['eps_p']
            self.assertAlmostEqual(
                eps_p[0] + eps_p[1] + eps_p[2], 0.0, places=12,
                msg=f"tr(ε^p) = {eps_p[0] + eps_p[1] + eps_p[2]:.3e} ≠ 0"
            )


class TestHex8VonMises3DConfinedVsStandalone(unittest.TestCase):
    """Hex8 + VonMises3D con confinamiento transversal vs material standalone.

    Setup: cubo unitario con todas las caras laterales confinadas:
    - x=0: ux=0;  x=1: ux libre, tracción aplicada
    - y=0, y=1: uy=0 en TODOS los nodos (confinamiento total en y)
    - z=0, z=1: uz=0 en TODOS los nodos (confinamiento total en z)

    Bajo esta BC, el campo es uniforme con ε = (ε_xx, 0, 0, 0, 0, 0) en
    todos los puntos de Gauss — análogo 3D al test plane strain confinado
    del 2D. Permite comparación directa con el material ejecutado standalone
    sobre la misma trayectoria monotónica de ε_xx.
    """

    def setUp(self):
        self.E = 2.0e5
        self.nu = 0.3
        self.sigma_y = 250.0
        self.H = 1.0e4
        self.mat = VonMises3D(
            E=self.E, nu=self.nu, sigma_y=self.sigma_y, H=self.H,
        )
        self.dom, self.elem, self.nodes = _build_hex8_unit_cube(self.mat)

        # Confinamiento total transverso: uy=0 y uz=0 en TODOS los nodos
        for n in self.nodes:
            n.fix_dof('uy', 0.0)
            n.fix_dof('uz', 0.0)
        # ux=0 en la cara x=0
        for i in HEX8_FACE_XMIN:
            self.nodes[i].fix_dof('ux', 0.0)

        self.dom.generate_equation_numbers(verbose=False)
        self.assembler = Assembler(self.dom)

    def _apply_traction_xmax(self, F_total: float) -> np.ndarray:
        F_ext = np.zeros(self.dom.total_dofs)
        for i in HEX8_FACE_XMAX:
            F_ext[self.nodes[i].dofs['ux']] = F_total / 4.0
        return F_ext

    def test_plastic_regime_matches_standalone_material(self):
        """σ_xx FEM (avg gauss) coincide con material standalone bajo el mismo path ε_xx."""
        # Yield elástico confinado: σ_xx = (K + 4G/3)·ε_xx, con yield cuando
        # ||s|| = sqrt(2/3)·σ_y → ε_xx = σ_y/(2G).
        G = self.E / (2.0 * (1.0 + self.nu))
        K = self.E / (3.0 * (1.0 - 2.0 * self.nu))
        eps_yield = self.sigma_y / (2.0 * G)
        sigma_yield_elastic = (K + 4.0 * G / 3.0) * eps_yield
        # Cargar 50% por encima del yield para entrar bien en plástico
        F_total = 1.5 * sigma_yield_elastic
        F_ext = self._apply_traction_xmax(F_total)

        conv = ConvergenceCriterion(rtol_force=1.0e-7, rtol_disp=1.0e-7)
        solver = NonlinearSolver(self.assembler, convergence=conv, num_steps=10)
        U = solver.solve(F_ext)

        u_x_final = U[self.nodes[1].dofs['ux']]
        self.assertGreater(u_x_final, eps_yield)

        gs = self.elem.compute_gauss_state(U)
        sigma_xx_fem = gs['stress'][:, 0].mean()
        alpha_fem = self.elem.state.vars[0]['alpha']

        # Material standalone bajo trayectoria monotónica ε_xx con ε_yy=ε_zz=0
        # (J2 + endurecimiento lineal monotónico es path-independent en ε).
        mat_ref = VonMises3D(
            E=self.E, nu=self.nu, sigma_y=self.sigma_y, H=self.H,
        )
        state_ref = None
        sigma_ref_last = None
        for eps_xx in np.linspace(0.0, u_x_final, 50):
            eps = np.array([eps_xx, 0.0, 0.0, 0.0, 0.0, 0.0])
            sigma_ref_last, _, state_ref = mat_ref.compute_state(eps, state_vars=state_ref)

        # σ_xx FEM vs σ_xx material standalone
        self.assertAlmostEqual(
            sigma_xx_fem, sigma_ref_last[0],
            delta=abs(sigma_ref_last[0]) * 1.0e-3
        )
        # α coincide
        self.assertAlmostEqual(
            alpha_fem, state_ref['alpha'],
            delta=abs(state_ref['alpha']) * 1.0e-3 + 1.0e-10
        )
        self.assertGreater(alpha_fem, 0.0)


class TestHex8VonMises3DConvergence(unittest.TestCase):
    """Newton global converge en pocas iteraciones por la tangente algorítmica consistente."""

    def setUp(self):
        self.E = 2.0e5
        self.nu = 0.3
        self.sigma_y = 250.0
        self.H = 1.0e4
        self.mat = VonMises3D(
            E=self.E, nu=self.nu, sigma_y=self.sigma_y, H=self.H,
        )
        self.dom, self.elem, self.nodes = _build_hex8_unit_cube(self.mat)
        # Confinamiento total transverso
        for n in self.nodes:
            n.fix_dof('uy', 0.0)
            n.fix_dof('uz', 0.0)
        for i in HEX8_FACE_XMIN:
            self.nodes[i].fix_dof('ux', 0.0)
        self.dom.generate_equation_numbers(verbose=False)
        self.assembler = Assembler(self.dom)

    def test_plastic_converges_in_few_iter(self):
        """Carga plástica significativa converge sin agotar max_iter."""
        # σ_xx = 1.5 · σ_yield_elastico_confinado
        G = self.E / (2.0 * (1.0 + self.nu))
        K = self.E / (3.0 * (1.0 - 2.0 * self.nu))
        eps_yield = self.sigma_y / (2.0 * G)
        sigma_yield_elastic = (K + 4.0 * G / 3.0) * eps_yield
        F_total = 1.5 * sigma_yield_elastic

        F_ext = np.zeros(self.dom.total_dofs)
        for i in HEX8_FACE_XMAX:
            F_ext[self.nodes[i].dofs['ux']] = F_total / 4.0

        conv = ConvergenceCriterion(rtol_force=1.0e-8, rtol_disp=1.0e-8)
        solver = NonlinearSolver(
            self.assembler, convergence=conv, num_steps=4, max_iter=8
        )
        U = solver.solve(F_ext)
        # Si llegamos aquí, convergió. Verificar que sí hubo plasticidad
        # (sanity check del setup).
        self.assertGreater(self.elem.state.vars[0]['alpha'], 0.0)


class TestTet4VonMises3DBasic(unittest.TestCase):
    """Smoke test Tet4 + VonMises3D.

    Tet4 con CST 3D (1 punto de Gauss baricéntrico) bajo tracción uniforme
    cruzando el yield. Verifica que el pipeline funciona end-to-end con
    el otro elemento sólido 3D del catálogo y que tr(ε^p) = 0 se preserva.
    """

    def test_pipeline_basico(self):
        E, nu, sigma_y, H = 2.0e5, 0.3, 250.0, 1.0e4
        mat = VonMises3D(E=E, nu=nu, sigma_y=sigma_y, H=H)

        dom = Domain()
        # Tet4 con vértice opuesto suficientemente alejado para evitar
        # jacobiano degenerado en el patch test.
        coords = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
        nodes = [dom.add_node(i + 1, list(c)) for i, c in enumerate(coords)]
        elem = Tet4(1, nodes, mat)
        dom.add_element(elem)

        # Confinamiento transverso total + tracción en x del nodo 1.
        for n in nodes:
            n.fix_dof('uy', 0.0)
            n.fix_dof('uz', 0.0)
        nodes[0].fix_dof('ux', 0.0)
        nodes[2].fix_dof('ux', 0.0)
        nodes[3].fix_dof('ux', 0.0)

        dom.generate_equation_numbers(verbose=False)
        assembler = Assembler(dom)

        # Carga del orden del yield confinado
        G = E / (2.0 * (1.0 + nu))
        K = E / (3.0 * (1.0 - 2.0 * nu))
        eps_yield = sigma_y / (2.0 * G)
        F_total = 1.3 * (K + 4.0 * G / 3.0) * eps_yield * 0.5  # área tet face = 0.5

        F_ext = np.zeros(dom.total_dofs)
        F_ext[nodes[1].dofs['ux']] = F_total

        conv = ConvergenceCriterion(rtol_force=1.0e-7, rtol_disp=1.0e-7)
        solver = NonlinearSolver(assembler, convergence=conv, num_steps=8)
        U = solver.solve(F_ext)

        # Plasticidad activa
        state = elem.state.vars[0]
        self.assertGreater(state['alpha'], 0.0)
        # Incompresibilidad plástica preservada
        eps_p = state['eps_p']
        self.assertAlmostEqual(
            eps_p[0] + eps_p[1] + eps_p[2], 0.0, places=12,
            msg=f"Tet4: tr(ε^p) = {eps_p[0] + eps_p[1] + eps_p[2]:.3e} ≠ 0"
        )


class TestHex8DruckerPrager3DElastic(unittest.TestCase):
    """Hex8 + DruckerPrager3D bajo carga por debajo del yield → respuesta elástica."""

    def test_elastic_response(self):
        E, nu = 2.0e4, 0.3
        mat = DruckerPrager3D(
            E=E, nu=nu, cohesion=10.0, phi_deg=30.0, psi_deg=10.0, H=100.0,
            variant='outer_cone',
        )
        dom, elem, nodes = _build_hex8_unit_cube(mat)

        # Rodillos en tres planos de simetría
        for i in HEX8_FACE_XMIN:
            nodes[i].fix_dof('ux', 0.0)
        for i in HEX8_FACE_YMIN:
            nodes[i].fix_dof('uy', 0.0)
        for i in HEX8_FACE_ZMIN:
            nodes[i].fix_dof('uz', 0.0)
        dom.generate_equation_numbers(verbose=False)
        assembler = Assembler(dom)

        # Carga muy pequeña → bien debajo del yield
        sigma_target = 0.1   # << cualquier umbral de fluencia
        F_ext = np.zeros(dom.total_dofs)
        for i in HEX8_FACE_XMAX:
            F_ext[nodes[i].dofs['ux']] = sigma_target / 4.0

        conv = ConvergenceCriterion(rtol_force=1.0e-9, rtol_disp=1.0e-9)
        solver = NonlinearSolver(assembler, convergence=conv, num_steps=2)
        U = solver.solve(F_ext)

        # α = 0 en todos los Gauss points (régimen elástico)
        for state in elem.state.vars:
            self.assertAlmostEqual(state['alpha'], 0.0, places=12)


class TestHex8DruckerPrager3DRegularPlastic(unittest.TestCase):
    """Hex8 + DruckerPrager3D no asociado en rama regular.

    Setup: cubo con confinamiento transversal total + cortante aplicado vía
    desplazamiento prescrito en la cara superior. Genera estado de cortante
    desviador dominante → rama regular (no apex).
    """

    def test_regular_branch_converges_with_asymmetric_tangent(self):
        E, nu = 2.0e4, 0.3
        # No asociada para activar tangente asimétrica
        mat = DruckerPrager3D(
            E=E, nu=nu, cohesion=10.0, phi_deg=30.0, psi_deg=10.0, H=500.0,
            variant='outer_cone',
        )
        # IS_SYMMETRIC=False (no asociada) → despachador usará LU
        self.assertFalse(mat.IS_SYMMETRIC)

        dom, elem, nodes = _build_hex8_unit_cube(mat)
        # BC: cara inferior empotrada; caras laterales con uy = uz = 0;
        # cara superior con desplazamiento ux > 0 prescrito (induce cortante xz)
        for i in HEX8_FACE_ZMIN:
            nodes[i].fix_dof('ux', 0.0)
            nodes[i].fix_dof('uy', 0.0)
            nodes[i].fix_dof('uz', 0.0)
        # Cortante puro: cara superior solo en ux ≠ 0
        u_top = 5.0e-3
        for i in (4, 5, 6, 7):   # cara z=1
            nodes[i].fix_dof('ux', u_top)
            nodes[i].fix_dof('uy', 0.0)
            nodes[i].fix_dof('uz', 0.0)
        # Confinamiento lateral suave en caras y=0, y=1
        for i in HEX8_FACE_YMIN:
            nodes[i].fix_dof('uy', 0.0)
        for i in (3, 2, 7, 6):    # cara y=1
            nodes[i].fix_dof('uy', 0.0)

        dom.generate_equation_numbers(verbose=False)
        assembler = Assembler(dom)

        F_ext = np.zeros(dom.total_dofs)   # solo BC Dirichlet
        conv = ConvergenceCriterion(rtol_force=1.0e-7, rtol_disp=1.0e-7)
        solver = NonlinearSolver(assembler, convergence=conv, num_steps=8)
        solver.solve(F_ext)

        # Plasticidad activa en todos los Gauss points
        for state in elem.state.vars:
            self.assertGreater(state['alpha'], 0.0)
            # Invariante cinemática: tr(ε^p) = 3·η_g·α
            eps_p = state['eps_p']
            tr_eps_p = eps_p[0] + eps_p[1] + eps_p[2]
            expected = 3.0 * mat.eta_g * state['alpha']
            self.assertAlmostEqual(
                tr_eps_p, expected, delta=abs(expected) * 1e-8 + 1e-12,
                msg=f"tr(ε^p)={tr_eps_p:.6e} ≠ 3·η_g·α={expected:.6e}"
            )


class TestHex8DruckerPrager3DApex(unittest.TestCase):
    """Hex8 + DruckerPrager3D bajo carga hidrostática traccionante → ápice.

    Bajo expansión radial uniforme del cubo (todos los nodos desplazados
    desde el centro), el estado en cada Gauss point es puramente hidrostático
    traccionante. Si la expansión es suficientemente grande, el predictor
    cae más allá del cono → return al ápice.
    """

    def test_apex_branch_under_hydrostatic_expansion(self):
        E, nu = 2.0e4, 0.3
        mat = DruckerPrager3D(
            E=E, nu=nu, cohesion=10.0, phi_deg=30.0, psi_deg=10.0, H=100.0,
            variant='outer_cone',
        )
        dom, elem, nodes = _build_hex8_unit_cube(mat)

        # Expansión hidrostática: cada nodo se desplaza desde el centro del cubo
        # (0.5, 0.5, 0.5) por un factor ε_iso. Todos los DOFs son Dirichlet.
        eps_iso = 1.5e-2   # suficientemente grande para forzar apex
        center = np.array([0.5, 0.5, 0.5])
        for n in nodes:
            pos = np.array(n.coordinates[:3], dtype=float)
            u = eps_iso * (pos - center)
            n.fix_dof('ux', u[0])
            n.fix_dof('uy', u[1])
            n.fix_dof('uz', u[2])

        dom.generate_equation_numbers(verbose=False)
        assembler = Assembler(dom)

        F_ext = np.zeros(dom.total_dofs)
        conv = ConvergenceCriterion(rtol_force=1.0e-7, rtol_disp=1.0e-7)
        solver = NonlinearSolver(assembler, convergence=conv, num_steps=4)
        U = solver.solve(F_ext)

        # En todos los Gauss points: estado hidrostático ⇒ σ_xx ≈ σ_yy ≈ σ_zz,
        # cortantes ≈ 0, √J₂ ≈ 0 (rama de ápice).
        gs = elem.compute_gauss_state(U)
        for k in range(8):
            sigma_gp = gs['stress'][k]
            self.assertAlmostEqual(sigma_gp[0], sigma_gp[1],
                                   delta=abs(sigma_gp[0]) * 1e-6)
            self.assertAlmostEqual(sigma_gp[0], sigma_gp[2],
                                   delta=abs(sigma_gp[0]) * 1e-6)
            for j in range(3, 6):
                self.assertAlmostEqual(
                    sigma_gp[j], 0.0, delta=abs(sigma_gp[0]) * 1e-6,
                    msg=f"GP {k}: cortante {j} = {sigma_gp[j]:.3e} ≠ 0"
                )
        # α > 0 en todos los Gauss points
        for state in elem.state.vars:
            self.assertGreater(state['alpha'], 0.0)


class TestTet4DruckerPrager3DBasic(unittest.TestCase):
    """Smoke test Tet4 + DruckerPrager3D."""

    def test_pipeline_basico(self):
        E, nu = 2.0e4, 0.3
        mat = DruckerPrager3D(
            E=E, nu=nu, cohesion=10.0, phi_deg=30.0, psi_deg=10.0, H=200.0,
            variant='outer_cone',
        )

        dom = Domain()
        coords = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
        nodes = [dom.add_node(i + 1, list(c)) for i, c in enumerate(coords)]
        elem = Tet4(1, nodes, mat)
        dom.add_element(elem)

        # Cortante puro entre nodo 0 y nodo 3 (vértice opuesto)
        for n in nodes:
            n.fix_dof('uy', 0.0)
            n.fix_dof('uz', 0.0)
        nodes[0].fix_dof('ux', 0.0)
        nodes[1].fix_dof('ux', 0.0)
        nodes[2].fix_dof('ux', 0.0)
        # Solo nodo 3 libre en ux (en z=1) → cortante xz
        nodes[3].fix_dof('ux', 8.0e-3)

        dom.generate_equation_numbers(verbose=False)
        assembler = Assembler(dom)

        F_ext = np.zeros(dom.total_dofs)
        conv = ConvergenceCriterion(rtol_force=1.0e-7, rtol_disp=1.0e-7)
        solver = NonlinearSolver(assembler, convergence=conv, num_steps=6)
        solver.solve(F_ext)

        state = elem.state.vars[0]
        self.assertGreater(state['alpha'], 0.0)
        # Invariante cinemática
        eps_p = state['eps_p']
        tr_eps_p = eps_p[0] + eps_p[1] + eps_p[2]
        expected = 3.0 * mat.eta_g * state['alpha']
        self.assertAlmostEqual(
            tr_eps_p, expected, delta=abs(expected) * 1e-8 + 1e-12,
            msg=f"Tet4 DP: tr(ε^p)={tr_eps_p:.6e} ≠ 3·η_g·α={expected:.6e}"
        )


class TestHex8IsotropicDamage3DElastic(unittest.TestCase):
    """Hex8 + IsotropicDamage3D bajo carga por debajo del umbral: respuesta elástica."""

    def test_elastic_response(self):
        E, nu, kappa_0, alpha = 2.0e4, 0.2, 1.0e-3, 500.0
        mat = IsotropicDamage3D(E=E, nu=nu, kappa_0=kappa_0, alpha=alpha)
        dom, elem, nodes = _build_hex8_unit_cube(mat)

        # Rodillos en tres planos de simetría
        for i in HEX8_FACE_XMIN:
            nodes[i].fix_dof('ux', 0.0)
        for i in HEX8_FACE_YMIN:
            nodes[i].fix_dof('uy', 0.0)
        for i in HEX8_FACE_ZMIN:
            nodes[i].fix_dof('uz', 0.0)
        dom.generate_equation_numbers(verbose=False)
        assembler = Assembler(dom)

        # Carga pequeña → ε_xx ~ 1e-5 << kappa_0
        sigma_target = E * 1.0e-5
        F_ext = np.zeros(dom.total_dofs)
        for i in HEX8_FACE_XMAX:
            F_ext[nodes[i].dofs['ux']] = sigma_target / 4.0

        conv = ConvergenceCriterion(rtol_force=1.0e-9, rtol_disp=1.0e-9)
        solver = NonlinearSolver(assembler, convergence=conv, num_steps=2)
        solver.solve(F_ext)

        # d = 0 en todos los Gauss points
        for state in elem.state.vars:
            self.assertAlmostEqual(state['damage'], 0.0, places=14)
            self.assertAlmostEqual(state['kappa'], kappa_0, places=14)


class TestHex8IsotropicDamage3DActiveLoad(unittest.TestCase):
    """Hex8 + IsotropicDamage3D bajo carga uniaxial activa con daño + tangente
    algorítmica consistente que activa LU.

    Aplicamos un desplazamiento prescrito ux en la cara +x de magnitud suficiente
    para superar κ_0; el resto de la cinemática se ajusta por equilibrio
    (confinamiento transversal libre con Poisson activo).
    """

    def test_loading_activates_damage_and_converges(self):
        E, nu, kappa_0, alpha = 2.0e4, 0.2, 1.0e-4, 500.0
        mat = IsotropicDamage3D(E=E, nu=nu, kappa_0=kappa_0, alpha=alpha)
        self.assertFalse(mat.IS_SYMMETRIC)   # tangente consistente asimétrica → LU

        dom, elem, nodes = _build_hex8_unit_cube(mat)
        # Rodillos en tres planos de simetría
        for i in HEX8_FACE_XMIN:
            nodes[i].fix_dof('ux', 0.0)
        for i in HEX8_FACE_YMIN:
            nodes[i].fix_dof('uy', 0.0)
        for i in HEX8_FACE_ZMIN:
            nodes[i].fix_dof('uz', 0.0)
        # Desplazamiento prescrito en cara +x para forzar daño activo
        u_x_prescribed = 8.0 * kappa_0   # ε_xx ~ 8·κ_0 ⇒ d significativo
        for i in HEX8_FACE_XMAX:
            nodes[i].fix_dof('ux', u_x_prescribed)
        dom.generate_equation_numbers(verbose=False)
        assembler = Assembler(dom)

        F_ext = np.zeros(dom.total_dofs)   # solo BC Dirichlet
        conv = ConvergenceCriterion(rtol_force=1.0e-7, rtol_disp=1.0e-7)
        # 8 pasos para entrar progresivamente en daño
        solver = NonlinearSolver(
            assembler, convergence=conv, num_steps=8, max_iter=10
        )
        solver.solve(F_ext)

        # d > 0 y κ > κ_0 en todos los Gauss points
        for k, state in enumerate(elem.state.vars):
            self.assertGreater(
                state['damage'], 0.0,
                msg=f"GP {k}: damage = {state['damage']:.3e} debería ser > 0"
            )
            self.assertGreater(
                state['kappa'], kappa_0,
                msg=f"GP {k}: kappa = {state['kappa']:.3e} debería ser > κ_0"
            )


class TestTet4IsotropicDamage3DBasic(unittest.TestCase):
    """Smoke test Tet4 + IsotropicDamage3D."""

    def test_pipeline_basico(self):
        E, nu, kappa_0, alpha = 2.0e4, 0.2, 1.0e-4, 500.0
        mat = IsotropicDamage3D(E=E, nu=nu, kappa_0=kappa_0, alpha=alpha)

        dom = Domain()
        coords = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
        nodes = [dom.add_node(i + 1, list(c)) for i, c in enumerate(coords)]
        elem = Tet4(1, nodes, mat)
        dom.add_element(elem)

        # Confinamiento transversal total + tracción en nodo 1
        for n in nodes:
            n.fix_dof('uy', 0.0)
            n.fix_dof('uz', 0.0)
        nodes[0].fix_dof('ux', 0.0)
        nodes[2].fix_dof('ux', 0.0)
        nodes[3].fix_dof('ux', 0.0)
        # Desplazamiento prescrito en nodo 1 que active daño
        nodes[1].fix_dof('ux', 5.0 * kappa_0)

        dom.generate_equation_numbers(verbose=False)
        assembler = Assembler(dom)

        F_ext = np.zeros(dom.total_dofs)
        conv = ConvergenceCriterion(rtol_force=1.0e-7, rtol_disp=1.0e-7)
        solver = NonlinearSolver(assembler, convergence=conv, num_steps=6)
        solver.solve(F_ext)

        state = elem.state.vars[0]
        self.assertGreater(state['damage'], 0.0)
        self.assertGreater(state['kappa'], kappa_0)


if __name__ == "__main__":
    unittest.main()
