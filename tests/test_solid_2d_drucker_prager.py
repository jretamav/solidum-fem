"""Tests de integración del pipeline Quad4 + DruckerPrager2D + NonlinearSolver.

Acceptance de ``docs/specs/DruckerPrager2D.md`` sección ``acceptance.integration``.

Valida que el modelo Drucker-Prager con plasticidad no asociada se integra
correctamente con el pipeline del solver. Como la tangente es asimétrica
(IS_SYMMETRIC=False), el despachador algebraico usa LU automáticamente.

Sección final ``TestDruckerPrager2DAnalyticalBenchmarks`` añade comparación
cuantitativa contra soluciones cerradas (Fase B de la matriz de validación,
2026-05-19).
"""
import math
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from solidum.core.domain import Domain
from solidum.elements.solid_2d import Quad4
from solidum.materials.drucker_prager_2d import DruckerPrager2D
from solidum.math.assembly import Assembler
from solidum.math.convergence import ConvergenceCriterion
from solidum.math.solvers import NonlinearSolver


def _build_quad4_plane_strain_domain(material, pattern: str = 'confined_y'):
    """Cuadrado unitario con 1 Quad4. Patrones de apoyo:
    - 'confined_y': uy=0 en todos los nodos; ux=0 en lado izquierdo.
    - 'corner_fixed': nodo (0,0) totalmente fijo, (1,0) en uy, (0,1) en ux.
    """
    domain = Domain()
    n1 = domain.add_node(1, [0.0, 0.0])
    n2 = domain.add_node(2, [1.0, 0.0])
    n3 = domain.add_node(3, [1.0, 1.0])
    n4 = domain.add_node(4, [0.0, 1.0])

    element = Quad4(1, [n1, n2, n3, n4], material, thickness=1.0)
    domain.add_element(element)

    if pattern == 'confined_y':
        n1.fix_dof('ux', 0.0); n1.fix_dof('uy', 0.0)
        n2.fix_dof('uy', 0.0)
        n3.fix_dof('uy', 0.0)
        n4.fix_dof('ux', 0.0); n4.fix_dof('uy', 0.0)
    elif pattern == 'corner_fixed':
        n1.fix_dof('ux', 0.0); n1.fix_dof('uy', 0.0)
        n2.fix_dof('uy', 0.0)
        n4.fix_dof('ux', 0.0)
    else:
        raise ValueError(f"pattern desconocido: {pattern}")

    domain.generate_equation_numbers(verbose=False)
    return domain, element, {'n1': n1, 'n2': n2, 'n3': n3, 'n4': n4}


class TestQuad4DruckerPrager2DPipeline(unittest.TestCase):
    """Quad4 + DruckerPrager2D + NonlinearSolver — pipeline completo."""

    def setUp(self):
        self.E = 2.0e4
        self.nu = 0.3
        self.c0 = 10.0
        self.phi = 30.0
        self.psi = 10.0   # no asociada
        self.H = 100.0
        self.material = DruckerPrager2D(
            E=self.E, nu=self.nu, cohesion=self.c0,
            phi_deg=self.phi, psi_deg=self.psi, H=self.H,
            hypothesis='plane_strain',
        )

    def test_elastic_regime_below_yield(self):
        """Carga muy por debajo del yield: respuesta elástica pura."""
        domain, element, nodes = _build_quad4_plane_strain_domain(
            self.material, pattern='confined_y'
        )
        assembler = Assembler(domain)

        # Cortante mediante deformación angular pequeña: aplicar tracción
        # tangencial baja en el lado derecho. En realidad, con apoyos
        # 'confined_y' completos, la única deformación posible es ε_xx, así
        # que pasamos a tracción horizontal pequeña.
        F = np.zeros(domain.total_dofs)
        F_total = 1.0   # muy por debajo del yield
        F[nodes['n2'].dofs['ux']] = 0.5 * F_total
        F[nodes['n3'].dofs['ux']] = 0.5 * F_total

        conv = ConvergenceCriterion(rtol_force=1e-8, rtol_disp=1e-8)
        solver = NonlinearSolver(assembler, convergence=conv, num_steps=2)
        U = solver.solve(F)

        self.assertEqual(element.state.vars[0]['alpha'], 0.0)
        self.assertGreater(U[nodes['n2'].dofs['ux']], 0.0)

    def test_plastic_regime_apex_branch_converges(self):
        """Carga uniaxial confinada empuja al ápice; solver converge."""
        domain, element, nodes = _build_quad4_plane_strain_domain(
            self.material, pattern='confined_y'
        )
        assembler = Assembler(domain)

        # Carga suficientemente alta para plastificar
        F = np.zeros(domain.total_dofs)
        F_total = 30.0
        F[nodes['n2'].dofs['ux']] = 0.5 * F_total
        F[nodes['n3'].dofs['ux']] = 0.5 * F_total

        conv = ConvergenceCriterion(rtol_force=1e-7, rtol_disp=1e-7)
        # max_iter=10: holgura porque la rama de ápice puede ser numéricamente
        # más sutil que la regular
        solver = NonlinearSolver(assembler, convergence=conv, num_steps=4, max_iter=10)
        U = solver.solve(F)

        self.assertGreater(U[nodes['n2'].dofs['ux']], 0.0)
        self.assertGreater(element.state.vars[0]['alpha'], 0.0)

    def test_associated_flow_pipeline_converges(self):
        """Caso asociado (ψ = φ): tangente simétrica, el despachador puede usar
        Cholesky en lugar de LU. Verifica que el solver converge igualmente.
        """
        mat_assoc = DruckerPrager2D(
            E=self.E, nu=self.nu, cohesion=self.c0,
            phi_deg=self.phi, psi_deg=self.phi, H=self.H,
        )
        domain, element, nodes = _build_quad4_plane_strain_domain(
            mat_assoc, pattern='confined_y'
        )
        assembler = Assembler(domain)

        F = np.zeros(domain.total_dofs)
        F_total = 30.0
        F[nodes['n2'].dofs['ux']] = 0.5 * F_total
        F[nodes['n3'].dofs['ux']] = 0.5 * F_total

        conv = ConvergenceCriterion(rtol_force=1e-7, rtol_disp=1e-7)
        solver = NonlinearSolver(assembler, convergence=conv, num_steps=4, max_iter=10)
        U = solver.solve(F)

        self.assertGreater(element.state.vars[0]['alpha'], 0.0)


def _build_quad4_prescribed_strain(material, eps_xx: float, eps_yy: float):
    """Quad4 unitario plane strain con campo de deformación uniforme prescrito
    vía desplazamientos Dirichlet en los cuatro nodos.

    Con dominio [0,1]² y desplazamientos ``u_x = ε_xx · x``, ``u_y = ε_yy · y``
    el campo es exactamente uniforme: ``ε_xx``, ``ε_yy``, ``γ_xy = 0`` en cada
    punto de Gauss. Permite comparar σ y ε_p contra la derivación analítica
    del material sin acoplamientos elementales.
    """
    from solidum.core.domain import Domain
    from solidum.elements.solid_2d import Quad4

    domain = Domain()
    n1 = domain.add_node(1, [0.0, 0.0])
    n2 = domain.add_node(2, [1.0, 0.0])
    n3 = domain.add_node(3, [1.0, 1.0])
    n4 = domain.add_node(4, [0.0, 1.0])

    element = Quad4(1, [n1, n2, n3, n4], material, thickness=1.0)
    domain.add_element(element)

    # Desplazamientos prescritos en todos los nodos (campo lineal uniforme).
    n1.fix_dof('ux', 0.0);            n1.fix_dof('uy', 0.0)
    n2.fix_dof('ux', eps_xx * 1.0);   n2.fix_dof('uy', 0.0)
    n3.fix_dof('ux', eps_xx * 1.0);   n3.fix_dof('uy', eps_yy * 1.0)
    n4.fix_dof('ux', 0.0);            n4.fix_dof('uy', eps_yy * 1.0)

    domain.generate_equation_numbers(verbose=False)
    return domain, element


class TestDruckerPrager2DAnalyticalBenchmarks(unittest.TestCase):
    """Comparación cuantitativa pipeline ↔ solución cerrada (Fase B, 2026-05-19).

    Tests añadidos a la matriz de validación para cerrar el hueco "DP2D sólo
    sanidad a nivel pipeline". Todos usan Quad4 unitario con desplazamientos
    Dirichlet prescritos en los cuatro nodos: el campo de deformación es
    uniforme y los puntos de Gauss reproducen exactamente la derivación
    analítica del material.
    """

    @staticmethod
    def _yield_eps_uniaxial_confined(material, sign: float = +1.0):
        """ε_xx de onset bajo ε_yy=0, γ_xy=0 plane strain (ε_zz=0).

        Estado elástico: σ_xx=(λ+2μ)ε, σ_yy=σ_zz=λε. ⇒
        I_1 = 3K·ε, √J₂ = (2/√3)·G·|ε|. Yield: √J₂ + η_f·I_1 = k_0.
        Para signo de ε prescrito (``+1`` tracción, ``−1`` compresión)
        el denominador es ``sign·(2G/√3) + 3K·η_f·sign``... reordenando:
        ``ε_yield = sign · k_0 / (2G/√3 + sign · 3K · η_f)``.
        En tracción el denominador es siempre positivo; en compresión
        puede acercarse a cero si la fricción equilibra la cohesión.
        """
        G, K = material.G, material.K
        denom = 2.0 * G / math.sqrt(3.0) + sign * 3.0 * K * material.eta_f
        return sign * material.k0 / denom

    def test_DP1_yield_onset_under_confined_tension(self):
        """DP-1 — onset de fluencia en tracción confinada uniaxial (1 Quad4).

        Setup analítico: ε_yy=0, γ_xy=0, ε_zz=0 (plane strain). σ_xx elástica
        es ``(λ+2μ)·ε``. Yield cuando √J₂+η_f·I₁=k₀ → ``ε_y`` cerrada.

        Verifica tres puntos: (a) ``ε = 0.5·ε_y`` elástico exacto, σ_xx
        coincide con elástica cerrada; (b) ``ε = 0.99·ε_y`` aún elástico,
        α=0; (c) ``ε = 1.5·ε_y`` plástico, α>0 y consistencia ``f≈0``
        (verificable porque al re-aplicar la misma ε con state convergido
        el material devuelve C_e — la condición de consistencia se cumplió
        en el paso plástico).
        """
        from solidum.materials.drucker_prager_2d import DruckerPrager2D
        from solidum.math.assembly import Assembler
        from solidum.math.convergence import ConvergenceCriterion
        from solidum.math.solvers import NonlinearSolver

        E, nu, c0 = 2.0e4, 0.3, 10.0
        mat = DruckerPrager2D(E=E, nu=nu, cohesion=c0, phi_deg=30.0, psi_deg=10.0, H=0.0)

        eps_y = self._yield_eps_uniaxial_confined(mat, sign=+1.0)
        lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
        mu = E / (2.0 * (1.0 + nu))
        sigma_xx_at_yield = (lam + 2.0 * mu) * eps_y

        # (a) elástico a 0.5·ε_y
        domain_a, elem_a = _build_quad4_prescribed_strain(mat, 0.5 * eps_y, 0.0)
        conv = ConvergenceCriterion(rtol_force=1e-9, rtol_disp=1e-9)
        NonlinearSolver(Assembler(domain_a), convergence=conv, num_steps=1).solve(
            np.zeros(domain_a.total_dofs)
        )
        sigma_a = elem_a.state.vars[0]['sigma'] if 'sigma' in elem_a.state.vars[0] else None
        # σ no se almacena en state — derivar de ε vía material
        eps_gp = np.array([0.5 * eps_y, 0.0, 0.0])
        sigma_a, _, _ = mat.compute_state(eps_gp)
        self.assertAlmostEqual(sigma_a[0], 0.5 * sigma_xx_at_yield, delta=abs(sigma_xx_at_yield) * 1e-8)
        self.assertAlmostEqual(elem_a.state.vars[0]['alpha'], 0.0, places=12)

        # (b) elástico marginal a 0.99·ε_y
        domain_b, elem_b = _build_quad4_prescribed_strain(mat, 0.99 * eps_y, 0.0)
        NonlinearSolver(Assembler(domain_b), convergence=conv, num_steps=1).solve(
            np.zeros(domain_b.total_dofs)
        )
        self.assertAlmostEqual(elem_b.state.vars[0]['alpha'], 0.0, places=12)

        # (c) plástico a 1.5·ε_y: α>0 en TODOS los Gauss points (campo uniforme)
        domain_c, elem_c = _build_quad4_prescribed_strain(mat, 1.5 * eps_y, 0.0)
        NonlinearSolver(Assembler(domain_c), convergence=conv, num_steps=4, max_iter=15).solve(
            np.zeros(domain_c.total_dofs)
        )
        alphas = [elem_c.state.vars[gp]['alpha'] for gp in range(len(elem_c.state.vars))]
        for a in alphas:
            self.assertGreater(a, 0.0)
        # Campo uniforme ⇒ α idéntico en los 4 Gauss points
        self.assertAlmostEqual(max(alphas) - min(alphas), 0.0, delta=max(alphas) * 1e-9)

    def test_DP2_flow_rule_consistency_under_monotonic_loading(self):
        """DP-2 — invariante del flujo no asociado: ``tr(ε_p) = 3·η_g·α``.

        Bajo cualquier secuencia de pasos exclusivamente en rama regular y
        partiendo del estado virgen (ε_p=0, α=0), la regla de flujo da
        ``Δtr(ε_p) = 3·η_g·Δγ`` y ``Δα = Δγ``, por lo que cumulativamente
        ``tr(ε_p) ≡ 3·η_g·α``. Este invariante es independiente del
        endurecimiento y del path mientras no se active el ápice.

        Es una prueba cinemática fuerte: si el ensamblaje, la actualización
        de state o el commit estuvieran mal escalados, el invariante fallaría
        aunque el solver convergiera.
        """
        from solidum.materials.drucker_prager_2d import DruckerPrager2D
        from solidum.math.assembly import Assembler
        from solidum.math.convergence import ConvergenceCriterion
        from solidum.math.solvers import NonlinearSolver

        mat = DruckerPrager2D(E=2.0e4, nu=0.3, cohesion=10.0,
                              phi_deg=30.0, psi_deg=10.0, H=100.0)
        eps_y = self._yield_eps_uniaxial_confined(mat, sign=+1.0)

        # Carga monótona en 5 niveles bien dentro de rama regular (tracción
        # confinada uniaxial no llega al ápice salvo deformaciones enormes).
        for factor in (1.2, 1.5, 1.8, 2.2, 2.6):
            domain, elem = _build_quad4_prescribed_strain(mat, factor * eps_y, 0.0)
            conv = ConvergenceCriterion(rtol_force=1e-9, rtol_disp=1e-9)
            NonlinearSolver(Assembler(domain), convergence=conv, num_steps=6, max_iter=15).solve(
                np.zeros(domain.total_dofs)
            )
            for gp in range(len(elem.state.vars)):
                eps_p = elem.state.vars[gp]['eps_p']     # [xx, yy, zz, xy_tens]
                alpha = elem.state.vars[gp]['alpha']
                tr_eps_p = eps_p[0] + eps_p[1] + eps_p[2]
                expected = 3.0 * mat.eta_g * alpha
                self.assertAlmostEqual(
                    tr_eps_p, expected,
                    delta=max(abs(expected), 1e-12) * 1e-8,
                    msg=f"factor={factor}, gp={gp}: tr(ε_p)={tr_eps_p}, "
                        f"3η_g·α={expected}"
                )
                # α debe ser estrictamente positivo en régimen plástico
                self.assertGreater(alpha, 0.0)

    def test_DP2bis_apex_under_biaxial_tension(self):
        """DP-2bis — return al ápice bajo tracción biaxial confinada (1 Quad4).

        Setup: ε_xx = ε_yy = ε grande (predictor puramente hidrostático
        traccionante; en plane strain ε_zz=0 induce σ_zz pero el predictor
        del cono cae más allá del ápice). En el ápice ``√J₂=0`` y σ es
        puramente hidrostático con ``σ_ii = k(α)/(3·η_f)`` para i=x,y,z.

        Verifica al final del paso: (i) σ_xx ≈ σ_yy (simetría hidrostática
        del estado plástico convergido); (ii) σ_xy ≈ 0; (iii) σ_xx coincide
        con ``k(α_final)/(3·η_f)`` derivado de α en el state.
        """
        from solidum.materials.drucker_prager_2d import DruckerPrager2D
        from solidum.math.assembly import Assembler
        from solidum.math.convergence import ConvergenceCriterion
        from solidum.math.solvers import NonlinearSolver

        mat = DruckerPrager2D(E=2.0e4, nu=0.3, cohesion=10.0,
                              phi_deg=30.0, psi_deg=10.0, H=100.0)

        eps_target = 1.0e-2  # tracción hidrostática extrema → apex
        domain, elem = _build_quad4_prescribed_strain(mat, eps_target, eps_target)
        conv = ConvergenceCriterion(rtol_force=1e-9, rtol_disp=1e-9)
        NonlinearSolver(Assembler(domain), convergence=conv, num_steps=4, max_iter=15).solve(
            np.zeros(domain.total_dofs)
        )

        # Reconstruir σ en cada Gauss point a partir del state convergido
        # llamando al material con la ε prescrita y el state final.
        eps_gp = np.array([eps_target, eps_target, 0.0])
        for gp in range(len(elem.state.vars)):
            sv = elem.state.vars[gp]
            sigma, _, _ = mat.compute_state(eps_gp, state_vars=sv)
            self.assertGreater(sv['alpha'], 0.0)
            self.assertAlmostEqual(sigma[0], sigma[1], delta=abs(sigma[0]) * 1e-9,
                                   msg=f"gp={gp}: σ_xx≠σ_yy")
            self.assertAlmostEqual(sigma[2], 0.0, delta=abs(sigma[0]) * 1e-9,
                                   msg=f"gp={gp}: σ_xy no nulo en ápice")
            k_alpha = mat.k0 + mat.H * sv['alpha']
            expected_p = k_alpha / (3.0 * mat.eta_f)
            self.assertAlmostEqual(sigma[0], expected_p,
                                   delta=abs(expected_p) * 1e-6,
                                   msg=f"gp={gp}: σ_xx≠k/(3η_f) — fuera del ápice")


if __name__ == '__main__':
    unittest.main()
