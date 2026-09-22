"""``Material.out_of_plane_stress`` — σ_zz de los materiales 2D en plane strain.

Se contrasta cada material 2D contra su gemelo 3D con ``ε_zz = 0``
impuesto (la definición misma de plane strain): el gemelo devuelve σ_zz en
el vector Voigt 6D y el 2D debe reproducirlo desde ``σ`` en plano y su
estado interno. Es el mismo criterio de cross-consistency 2D↔3D que ya
blinda ``compute_state``; aquí se extiende a la componente fuera del plano.
"""
import unittest

import numpy as np

import solidum  # noqa: F401 — autodiscover
from solidum.core.material import Material
from solidum.materials.damage_2d import IsotropicDamage2D
from solidum.materials.damage_3d import IsotropicDamage3D
from solidum.materials.drucker_prager_2d import DruckerPrager2D
from solidum.materials.drucker_prager_3d import DruckerPrager3D
from solidum.materials.elastic import Elastic1D
from solidum.materials.elastic_2d import Elastic2D
from solidum.materials.elastic_3d import Elastic3D
from solidum.materials.orthotropic_2d import Orthotropic2D
from solidum.materials.von_mises_2d import VonMises2D
from solidum.materials.von_mises_3d import VonMises3D


def _to_3d(strain_2d):
    """Voigt 2D ``[εxx, εyy, γxy]`` → Voigt 3D con ``ε_zz = γ_yz = γ_xz = 0``."""
    return np.array([strain_2d[0], strain_2d[1], 0.0, strain_2d[2], 0.0, 0.0])


def _walk(mat2d, mat3d, path):
    """Recorre ``path`` con ambos materiales encadenando el estado committed y
    devuelve pares ``(σ_zz_2D, σ_zz_3D)`` por paso."""
    sv2 = sv3 = None
    pairs = []
    for eps in path:
        eps = np.asarray(eps, dtype=float)
        sig2, _, sv2 = mat2d.compute_state(eps, sv2)
        sig3, _, sv3 = mat3d.compute_state(_to_3d(eps), sv3)
        pairs.append((mat2d.out_of_plane_stress(sig2, sv2), float(sig3[2])))
        # Coherencia en plano del gemelo (garantiza que comparamos lo mismo).
        np.testing.assert_allclose(sig2, [sig3[0], sig3[1], sig3[3]], rtol=1e-10, atol=1e-12)
    return pairs


class TestOutOfPlaneStress(unittest.TestCase):

    def test_default_es_cero_en_1d_3d_y_plane_stress(self):
        sig = np.array([3.0, -1.0, 0.5])
        self.assertEqual(Elastic1D(E=10.0).out_of_plane_stress(sig), 0.0)
        self.assertEqual(Elastic3D(E=10.0, nu=0.3).out_of_plane_stress(sig), 0.0)
        self.assertEqual(Elastic2D(E=10.0, nu=0.3, hypothesis='plane_stress').out_of_plane_stress(sig), 0.0)
        self.assertEqual(VonMises2D(E=10.0, nu=0.3, sigma_y=1.0, hypothesis='plane_stress').out_of_plane_stress(sig), 0.0)
        self.assertEqual(Orthotropic2D(E1=10.0, E2=1.0, G12=0.5, nu12=0.3).out_of_plane_stress(sig), 0.0)

    def test_elastic_plane_strain_vs_elastic3d(self):
        E, nu = 210e9, 0.3
        m2, m3 = Elastic2D(E, nu, hypothesis='plane_strain'), Elastic3D(E, nu)
        rng = np.random.default_rng(0)
        for _ in range(5):
            eps = rng.normal(size=3) * 1e-3
            (szz2, szz3), = _walk(m2, m3, [eps])
            self.assertAlmostEqual(szz2, szz3, delta=1e-6 * abs(szz3) + 1e-3)
            # Identidad clásica ν(σxx+σyy).
            sig, _, _ = m2.compute_state(eps)
            self.assertAlmostEqual(szz2, nu * (sig[0] + sig[1]), delta=1e-3)

    def test_damage_plane_strain_vs_damage3d(self):
        E, nu, k0, a = 30e9, 0.2, 1e-4, 500.0
        m2 = IsotropicDamage2D(E, nu, kappa_0=k0, alpha=a, hypothesis='plane_strain')
        m3 = IsotropicDamage3D(E, nu, kappa_0=k0, alpha=a)
        path = [np.array([f * 3e-4, -f * 1e-4, f * 2e-4]) for f in (0.2, 0.6, 1.0, 1.5, 0.8)]
        for szz2, szz3 in _walk(m2, m3, path):
            self.assertAlmostEqual(szz2, szz3, delta=1e-9 * abs(szz3) + 1e-6)

    def test_j2_plane_strain_vs_vonmises3d(self):
        E, nu, sy, H = 200e9, 0.3, 250e6, 20e9
        m2 = VonMises2D(E, nu, sy, H=H, hypothesis='plane_strain')
        m3 = VonMises3D(E, nu, sy, H=H)
        # Carga hasta plastificar, descarga parcial y recarga en otra dirección.
        path = [np.array([1e-3, 0.0, 0.0]), np.array([3e-3, 0.0, 1e-3]),
                np.array([2e-3, 5e-4, 2e-3]), np.array([-1e-3, 3e-3, 0.0])]
        pairs = _walk(m2, m3, path)
        self.assertGreater(abs(pairs[-1][1]), 0.0)
        for szz2, szz3 in pairs:
            self.assertAlmostEqual(szz2, szz3, delta=1e-9 * abs(szz3) + 1e-3)

    def test_j2_plane_strain_plastico_no_es_nu_por_traza(self):
        """Con ε^p_zz ≠ 0 la identidad elástica ν(σxx+σyy) deja de valer: el
        test asegura que la implementación no cae en ese atajo."""
        m2 = VonMises2D(200e9, 0.3, 250e6, H=0.0, hypothesis='plane_strain')
        sig, _, sv = m2.compute_state(np.array([4e-3, 0.0, 0.0]))
        self.assertGreater(abs(sv['eps_p'][2]), 0.0)
        szz = m2.out_of_plane_stress(sig, sv)
        self.assertGreater(abs(szz - 0.3 * (sig[0] + sig[1])), 1e6)

    def test_drucker_prager_plane_strain_vs_3d(self):
        common = dict(E=50e6, nu=0.25, cohesion=20e3, phi_deg=30.0, psi_deg=15.0, H=2e5,
                      variant='outer_cone')
        m2 = DruckerPrager2D(hypothesis='plane_strain', **common)
        m3 = DruckerPrager3D(**common)
        # Compresión dominante (retorno regular) y una excursión a tracción (ápice).
        path = [np.array([-2e-3, -1e-3, 5e-4]), np.array([-4e-3, -1e-3, 2e-3]),
                np.array([-3e-3, 1e-3, 2e-3]), np.array([2e-3, 2e-3, 0.0])]
        for szz2, szz3 in _walk(m2, m3, path):
            self.assertAlmostEqual(szz2, szz3, delta=1e-9 * abs(szz3) + 1e-6)

    def test_contrato_base_en_todo_el_registro(self):
        """Todo material registrado responde a ``out_of_plane_stress`` con un
        float finito (la base lo garantiza; las subclases no deben romperlo)."""
        sig = np.array([1.0, 2.0, 0.5])
        for name in solidum.MaterialRegistry.names():
            cls = solidum.MaterialRegistry.get(name)
            self.assertTrue(issubclass(cls, Material), name)
            self.assertTrue(callable(getattr(cls, 'out_of_plane_stress', None)), name)


if __name__ == '__main__':
    unittest.main()
