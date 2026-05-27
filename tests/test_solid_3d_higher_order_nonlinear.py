"""Cross-check: sólidos 3D cuadráticos × materiales 3D no lineales.

Sub-fase 5 de A.ter (2026-05-27). Cierra la matriz elemento × material
3D verificando que cada combinación de elemento cuadrático
({Hex20, Hex27, Tet10}) con cada material no lineal 3D
({VonMises3D, DruckerPrager3D, IsotropicDamage3D}) pasa por el pipeline
``Element + NonlinearSolver`` sin bugs de cableado y produce estados
internos sensatos (variables internas no-cero en el régimen post-yield
o post-pico).

Patrón paritario con los smoke tests existentes para ``Tet4`` en
``tests/test_solid_3d_plasticity.py``: un elemento, BCs de
confinamiento simple, multi-step Newton, verificación de que el
material entra en régimen no lineal.

Cierra la matriz §1 de MATRIZ.md para las 9 celdas ○ → ✓.
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from solidum.core.domain import Domain
from solidum.elements.solid_3d import Hex20, Hex27, Tet10
from solidum.materials.damage_3d import IsotropicDamage3D
from solidum.materials.drucker_prager_3d import DruckerPrager3D
from solidum.materials.von_mises_3d import VonMises3D
from solidum.math.assembly import Assembler
from solidum.math.convergence import ConvergenceCriterion
from solidum.math.solvers import NonlinearSolver


# =============================================================================
# Coordenadas de referencia (cubo unitario para hexaedros, tet ref para Tet10)
# =============================================================================

HEX20_UNIT = [
    (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
    (0.5, 0.0, 0.0), (1.0, 0.5, 0.0), (0.5, 1.0, 0.0), (0.0, 0.5, 0.0),
    (0.5, 0.0, 1.0), (1.0, 0.5, 1.0), (0.5, 1.0, 1.0), (0.0, 0.5, 1.0),
    (0.0, 0.0, 0.5), (1.0, 0.0, 0.5), (1.0, 1.0, 0.5), (0.0, 1.0, 0.5),
]

HEX27_UNIT = HEX20_UNIT + [
    (0.0, 0.5, 0.5), (1.0, 0.5, 0.5),
    (0.5, 0.0, 0.5), (0.5, 1.0, 0.5),
    (0.5, 0.5, 0.0), (0.5, 0.5, 1.0),
    (0.5, 0.5, 0.5),
]

TET10_REF = [
    (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0),
    (0.5, 0.0, 0.0), (0.5, 0.5, 0.0), (0.0, 0.5, 0.0),
    (0.0, 0.0, 0.5), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5),
]


def _build_single(element_cls, coords, material):
    """Construye un dominio con un único elemento del tipo dado."""
    dom = Domain()
    nodes = [dom.add_node(i + 1, list(c)) for i, c in enumerate(coords)]
    elem = element_cls(1, nodes, material)
    dom.add_element(elem)
    return dom, elem, nodes


def _apply_confined_uniaxial_bcs(dom, prescribed_ux: float):
    """BCs de confinamiento total transverso + tracción/desplazamiento uniaxial.

    - Todos los nodos: ``u_y = u_z = 0`` (confinamiento total transverso).
    - Nodos en ``x = 0``: ``u_x = 0`` (cara empotrada).
    - Nodos en ``x = 1``: ``u_x = prescribed_ux`` (desplazamiento prescrito).
    - Nodos con ``0 < x < 1``: ``u_x`` libre (solver lo resuelve).

    Este confinamiento es independiente del número y tipo de nodos del
    elemento; funciona para cualquier elemento sólido 3D del catálogo.
    """
    for n in dom.nodes.values():
        n.fix_dof('uy', 0.0)
        n.fix_dof('uz', 0.0)
        if abs(n.coordinates[0]) < 1.0e-12:
            n.fix_dof('ux', 0.0)
        elif abs(n.coordinates[0] - 1.0) < 1.0e-12:
            n.fix_dof('ux', prescribed_ux)


def _solve_steps(dom, num_steps: int = 6):
    dom.generate_equation_numbers(verbose=False)
    assembler = Assembler(dom)
    F_ext = np.zeros(dom.total_dofs)
    conv = ConvergenceCriterion(rtol_force=1.0e-7, rtol_disp=1.0e-7)
    solver = NonlinearSolver(assembler, convergence=conv, num_steps=num_steps)
    return solver.solve(F_ext)


# =============================================================================
# VonMises3D × {Hex20, Hex27, Tet10}
# =============================================================================

class _VonMises3DSmokeMixin:
    """Verifica que el material plastifica (α > 0 en algún punto de Gauss)
    bajo tracción confinada que rebasa el yield uniaxial."""

    ELEMENT_CLS = None
    COORDS = ()

    def test_pipeline_basico(self):
        E, nu, sigma_y, H = 2.0e5, 0.3, 250.0, 1.0e4
        mat = VonMises3D(E=E, nu=nu, sigma_y=sigma_y, H=H)
        dom, elem, _ = _build_single(self.ELEMENT_CLS, self.COORDS, mat)

        # Tracción uniaxial confinada con desplazamiento prescrito
        # > yield strain por un factor cómodo.
        eps_yield_unconfined = sigma_y / E
        _apply_confined_uniaxial_bcs(dom, prescribed_ux=2.5 * eps_yield_unconfined)
        _solve_steps(dom, num_steps=6)

        alphas = [s['alpha'] for s in elem.state.vars]
        max_alpha = max(alphas)
        self.assertGreater(
            max_alpha, 0.0,
            msg=f"{self.ELEMENT_CLS.__name__} + VonMises3D no plastifica: "
                f"max α = {max_alpha:.3e} (esperado > 0)",
        )


class TestHex20VonMises3DSmoke(_VonMises3DSmokeMixin, unittest.TestCase):
    ELEMENT_CLS = Hex20
    COORDS = HEX20_UNIT


class TestHex27VonMises3DSmoke(_VonMises3DSmokeMixin, unittest.TestCase):
    ELEMENT_CLS = Hex27
    COORDS = HEX27_UNIT


class TestTet10VonMises3DSmoke(_VonMises3DSmokeMixin, unittest.TestCase):
    ELEMENT_CLS = Tet10
    COORDS = TET10_REF


# =============================================================================
# DruckerPrager3D × {Hex20, Hex27, Tet10}
# =============================================================================

class _DruckerPrager3DSmokeMixin:
    """Verifica plasticidad friccional bajo cortante puro confinado."""

    ELEMENT_CLS = None
    COORDS = ()
    PRESCRIBED_UZ = 0.0  # subclase puede sobreescribir

    def test_pipeline_basico(self):
        E, nu = 2.0e4, 0.3
        mat = DruckerPrager3D(
            E=E, nu=nu, cohesion=10.0, phi_deg=30.0, psi_deg=10.0, H=200.0,
            variant='outer_cone',
        )
        dom, elem, _ = _build_single(self.ELEMENT_CLS, self.COORDS, mat)

        # Cortante: desplazamiento prescrito en x sobre cara x=1, confinamiento.
        _apply_confined_uniaxial_bcs(dom, prescribed_ux=8.0e-3)
        _solve_steps(dom, num_steps=8)

        alphas = [s['alpha'] for s in elem.state.vars]
        max_alpha = max(alphas)
        self.assertGreater(
            max_alpha, 0.0,
            msg=f"{self.ELEMENT_CLS.__name__} + DruckerPrager3D no plastifica: "
                f"max α = {max_alpha:.3e}",
        )


class TestHex20DruckerPrager3DSmoke(_DruckerPrager3DSmokeMixin, unittest.TestCase):
    ELEMENT_CLS = Hex20
    COORDS = HEX20_UNIT


class TestHex27DruckerPrager3DSmoke(_DruckerPrager3DSmokeMixin, unittest.TestCase):
    ELEMENT_CLS = Hex27
    COORDS = HEX27_UNIT


class TestTet10DruckerPrager3DSmoke(_DruckerPrager3DSmokeMixin, unittest.TestCase):
    ELEMENT_CLS = Tet10
    COORDS = TET10_REF


# =============================================================================
# IsotropicDamage3D × {Hex20, Hex27, Tet10}
# =============================================================================

class _IsotropicDamage3DSmokeMixin:
    """Verifica que el daño se activa (ω > 0, κ > κ_0) bajo tracción
    confinada que rebasa el umbral de daño."""

    ELEMENT_CLS = None
    COORDS = ()

    def test_pipeline_basico(self):
        E, nu, kappa_0, alpha_dam = 2.0e4, 0.2, 1.0e-4, 500.0
        mat = IsotropicDamage3D(E=E, nu=nu, kappa_0=kappa_0, alpha=alpha_dam)
        dom, elem, _ = _build_single(self.ELEMENT_CLS, self.COORDS, mat)

        # Desplazamiento prescrito en x rebasando κ_0 con margen cómodo
        _apply_confined_uniaxial_bcs(dom, prescribed_ux=5.0 * kappa_0)
        _solve_steps(dom, num_steps=6)

        damages = [s['damage'] for s in elem.state.vars]
        kappas = [s['kappa'] for s in elem.state.vars]
        max_dam = max(damages)
        max_kappa = max(kappas)
        self.assertGreater(
            max_dam, 0.0,
            msg=f"{self.ELEMENT_CLS.__name__} + IsotropicDamage3D no daña: "
                f"max ω = {max_dam:.3e}",
        )
        self.assertGreater(
            max_kappa, kappa_0,
            msg=f"{self.ELEMENT_CLS.__name__} + IsotropicDamage3D no crece κ: "
                f"max κ = {max_kappa:.3e} (esperaba > κ_0 = {kappa_0})",
        )


class TestHex20IsotropicDamage3DSmoke(_IsotropicDamage3DSmokeMixin, unittest.TestCase):
    ELEMENT_CLS = Hex20
    COORDS = HEX20_UNIT


class TestHex27IsotropicDamage3DSmoke(_IsotropicDamage3DSmokeMixin, unittest.TestCase):
    ELEMENT_CLS = Hex27
    COORDS = HEX27_UNIT


class TestTet10IsotropicDamage3DSmoke(_IsotropicDamage3DSmokeMixin, unittest.TestCase):
    ELEMENT_CLS = Tet10
    COORDS = TET10_REF


if __name__ == "__main__":
    unittest.main()
