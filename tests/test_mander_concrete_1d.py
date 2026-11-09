"""``ManderConcrete1D``: verificación del material (spec ``ManderConcrete1D``).

Casos de ``acceptance.verification``: puntos notables de la envolvente,
deformación plástica de Mander (ecs. 31–34, evaluadas aquí de forma
independiente), descarga y recarga lineales, tracción con historia (factor de
Mander, descarga secante, cierre de la grieta), regularización por energía de
fractura en tracción y en compresión, función de confinamiento (ec. 29 de
Mander, aproximación de Chang y Mander 1994, ejemplo numérico de Mander),
tangente consistente por diferencias finitas y disipación no negativa en
trayectorias aleatorias.
"""
import logging
import math
import unittest

import numpy as np

import solidum  # noqa: F401  — autodiscover
from solidum.core.node import Node
from solidum.elements.truss import Truss2D
from solidum.materials.concrete_1d import ManderConcrete1D
from solidum.materials.concrete_confinement import (
    confined_strength_ratio,
    confinement_effectiveness_circular,
    confinement_effectiveness_rectangular,
    fib_fracture_energy,
    mander_confined_concrete,
    mander_elastic_modulus,
)

# Núcleo y recubrimiento de la spec (unidades SI).
CORE = dict(fc=34.5e6, eps_c=0.0045, eps_cu=0.015, E=25.0e9)
FT_COVER = 2.5e6
E_COVER = mander_elastic_modulus(25.0e6)
COVER = dict(fc=25.0e6, eps_c=0.002, eps_cu=0.004, eps_sp=0.006, E=E_COVER,
             ft=FT_COVER, eps_tu=10.0 * FT_COVER / E_COVER)


# ---------------------------------------------------------------------------
# Fórmulas independientes (no usan el código del material)
# ---------------------------------------------------------------------------

def popovics(e, fc, eps_c, E):
    r = E / (E - fc / eps_c)
    x = e / eps_c
    return fc * x * r / (r - 1.0 + x ** r)


def mander_plastic(e_un, f_un, eps_c, E):
    """Mander et al. (1988), ecs. 31–34."""
    a = max(eps_c / (eps_c + e_un), 0.09 * e_un / eps_c)
    e_a = a * math.sqrt(e_un * eps_c)
    return e_un - (e_un + e_a) * f_un / (f_un + E * e_a)


def chang_mander_K(q1, q2):
    """Chang y Mander (1994), ecs. 3-70a–e."""
    r = min(q1, q2) / max(q1, q2)
    x = 0.5 * (q1 + q2)
    A = 6.8886 - (0.6069 + 17.275 * r) * math.exp(-4.989 * r)
    B = 4.5 / ((5.0 / A) * (0.9849 - 0.6306 * math.exp(-3.8939 * r)) - 0.1) - 5.0
    return 1.0 + A * x * (0.1 + 0.9 / (1.0 + B * x))


def mander_K29(q):
    return -1.254 + 2.254 * math.sqrt(1.0 + 7.94 * q) - 2.0 * q


def walk(mat, path, state=None):
    """Recorre ``path`` comprometiendo cada punto; devuelve (σ, E_t, estados)."""
    sig, tan, states = [], [], []
    for eps in path:
        s, t, state = mat.compute_state(eps, state)
        sig.append(s)
        tan.append(t)
        states.append(state)
    return np.array(sig), np.array(tan), states


def stored_energy(mat, eps, sigma, state):
    """Energía libre del punto (spec §3–§4): ½·s·(e − e_p) en compresión y
    ½·σ·t en tracción, con t = ε + e_p."""
    em, h = state["eps_cmax"], state["h"]
    ep = mat.plastic_strain(em, mat.compression_envelope(em, h)[0]) if em > 0 else 0.0
    t = eps + ep
    return 0.5 * sigma * t


class TestEnvolvente(unittest.TestCase):

    def test_puntos_notables(self):
        core, cover = ManderConcrete1D(**CORE), ManderConcrete1D(**COVER)
        for mat, p in ((core, CORE), (cover, COVER)):
            with self.subTest(fc=p["fc"]):
                s, t, _ = mat.compute_state(-p["eps_c"])
                self.assertAlmostEqual(s / -p["fc"], 1.0, delta=1e-12)
                s, _, _ = mat.compute_state(-p["eps_cu"])
                self.assertAlmostEqual(
                    s / -popovics(p["eps_cu"], p["fc"], p["eps_c"], p["E"]), 1.0, delta=1e-12)
                _, t0, _ = mat.compute_state(-1e-12)
                self.assertAlmostEqual(t0 / p["E"], 1.0, delta=1e-6)
        s, _, _ = cover.compute_state(-COVER["eps_sp"])
        self.assertEqual(s, 0.0)
        eps_t = FT_COVER / E_COVER
        s, _, _ = cover.compute_state(eps_t)
        self.assertAlmostEqual(s / FT_COVER, 1.0, delta=1e-12)
        s, _, _ = cover.compute_state(COVER["eps_tu"])
        self.assertAlmostEqual(s, 0.0, delta=1e-9 * FT_COVER)
        # Núcleo: caída a cero al pasar eps_cu (eps_sp = eps_cu por omisión).
        s, _, _ = core.compute_state(-1.0001 * CORE["eps_cu"])
        self.assertEqual(s, 0.0)

    def test_validaciones(self):
        with self.assertRaisesRegex(ValueError, "secante al pico"):
            ManderConcrete1D(fc=30e6, eps_c=0.002, E=10e9, eps_cu=0.004)
        with self.assertRaisesRegex(ValueError, "uno de los dos"):
            ManderConcrete1D(fc=30e6, eps_c=0.002, E=27e9, eps_cu=0.004, ft=3e6)
        with self.assertRaisesRegex(ValueError, "uno de los dos"):
            ManderConcrete1D(fc=30e6, eps_c=0.002, E=27e9, eps_cu=0.004, ft=3e6,
                             eps_tu=1e-3, G_f=140.0)
        with self.assertRaisesRegex(ValueError, "sin tracción"):
            ManderConcrete1D(fc=30e6, eps_c=0.002, E=27e9, eps_cu=0.004, eps_tu=1e-3)
        with self.assertRaisesRegex(ValueError, "eps_sp"):
            ManderConcrete1D(fc=30e6, eps_c=0.002, E=27e9, eps_cu=0.004, eps_sp=0.003)
        with self.assertRaisesRegex(ValueError, "faltan"):
            ManderConcrete1D(fc=30e6, eps_c=0.002, E=27e9)


class TestDescargaEnCompresion(unittest.TestCase):

    def test_deformacion_plastica_de_mander(self):
        for p in (CORE, COVER):
            mat = ManderConcrete1D(**p)
            for eta in (0.25, 0.5, 1.0, 2.0, 3.0):
                em = eta * p["eps_c"]
                if em > mat.eps_sp:
                    continue
                with self.subTest(fc=p["fc"], eta=eta):
                    s_m = -mat.compute_state(-em)[0]
                    ep = mat.plastic_strain(em, s_m)
                    self.assertAlmostEqual(
                        ep, mander_plastic(em, s_m, p["eps_c"], p["E"]),
                        delta=1e-12 * em)
                    self.assertGreaterEqual(ep, 0.0)
                    self.assertLessEqual(ep, em)
                    if s_m > 0:
                        self.assertLessEqual(s_m / (em - ep), p["E"] * (1 + 1e-12))

    def test_descarga_y_recarga_por_la_misma_recta(self):
        mat = ManderConcrete1D(**CORE)
        for eta in (0.6, 1.8):
            em = eta * CORE["eps_c"]
            with self.subTest(eta=eta):
                s_m, _, state = mat.compute_state(-em)
                s_m = -s_m
                ep = mat.plastic_strain(em, s_m)
                E_r = s_m / (em - ep)
                down = -np.linspace(em, ep, 9)
                sig, tan, states = walk(mat, down[1:], state)
                np.testing.assert_allclose(-sig, E_r * (-down[1:] - ep), rtol=1e-12,
                                           atol=1e-9 * s_m)
                np.testing.assert_allclose(tan, E_r, rtol=1e-12)
                self.assertTrue(all(st["eps_cmax"] == em for st in states))
                up = -np.linspace(ep, em, 9)
                sig2, _, states2 = walk(mat, up[1:], states[-1])
                np.testing.assert_allclose(-sig2, E_r * (-up[1:] - ep), rtol=1e-12,
                                           atol=1e-9 * s_m)
                # Pasado e_m, de vuelta en la envolvente.
                s3, _, st3 = mat.compute_state(-1.05 * em, states2[-1])
                self.assertAlmostEqual(-s3 / popovics(1.05 * em, CORE["fc"], CORE["eps_c"],
                                                      CORE["E"]), 1.0, delta=1e-12)
                self.assertAlmostEqual(st3["eps_cmax"], 1.05 * em, delta=1e-18)

    def test_aplastada_no_resiste(self):
        mat = ManderConcrete1D(**COVER)
        _, _, st = mat.compute_state(-1.2 * COVER["eps_sp"])
        for eps in (-1.5 * COVER["eps_sp"], -0.5 * COVER["eps_sp"], 0.0, 1e-3):
            s, t, _ = mat.compute_state(eps, st)
            self.assertEqual((s, t), (0.0, 0.0))


class TestTraccionConHistoria(unittest.TestCase):

    def test_factor_de_mander(self):
        p = dict(fc=30e6, eps_c=0.002, E=mander_elastic_modulus(30e6), eps_cu=0.006,
                 ft=3e6, eps_tu=5 * 3e6 / mander_elastic_modulus(30e6))
        mat = ManderConcrete1D(**p)
        expected = {0.5: 0.96, 1.0: 0.78, 1.5: 0.50, 2.0: 0.16, 2.5: 0.0}
        for eta, phi_tab in expected.items():
            with self.subTest(eta=eta):
                em = eta * p["eps_c"]
                s_m, _, st = mat.compute_state(-em)
                ep = mander_plastic(em, -s_m, p["eps_c"], p["E"])
                phi = max(0.0, 1.0 - ep / p["eps_c"])
                self.assertAlmostEqual(phi, phi_tab, delta=0.006)   # tabla de la spec
                eps_t = p["ft"] / p["E"]
                s, t, _ = mat.compute_state(-ep + 0.5 * eps_t, st)
                self.assertAlmostEqual(s, phi * p["E"] * 0.5 * eps_t, delta=1e-9 * p["ft"])
                self.assertAlmostEqual(t, phi * p["E"], delta=1e-9 * p["E"])
                s, _, _ = mat.compute_state(-ep + eps_t, st)
                self.assertAlmostEqual(s, phi * p["ft"], delta=1e-9 * p["ft"])

    def test_descarga_secante_y_cierre_de_grieta(self):
        mat = ManderConcrete1D(**COVER)
        # Compresión previa (deja e_p y daño), después tracción hasta agrietar.
        s_m, _, st = mat.compute_state(-1.5 * COVER["eps_c"])
        em = 1.5 * COVER["eps_c"]
        ep = mat.plastic_strain(em, -s_m)
        E_r = -s_m / (em - ep)
        phi = 1.0 - ep / COVER["eps_c"]
        eps_t = FT_COVER / E_COVER
        t_m = 3.0 * eps_t
        s_crack, _, st = mat.compute_state(-ep + t_m, st)
        sig_env = phi * (FT_COVER - FT_COVER / (COVER["eps_tu"] - eps_t) * (t_m - eps_t))
        self.assertAlmostEqual(s_crack, sig_env, delta=1e-9 * FT_COVER)
        # Descarga secante hacia el origen desplazado −e_p.
        for frac in (0.75, 0.3):
            s, t, _ = mat.compute_state(-ep + frac * t_m, st)
            self.assertAlmostEqual(s, sig_env * frac, delta=1e-9 * FT_COVER)
            self.assertAlmostEqual(t, sig_env / t_m, delta=1e-9 * E_COVER)
        # Cierre: en compresión la rigidez es la de antes de agrietar.
        s, t, _ = mat.compute_state(-ep - 1e-5, st)
        self.assertAlmostEqual(t / E_r, 1.0, delta=1e-12)
        self.assertAlmostEqual(-s / (E_r * 1e-5), 1.0, delta=1e-9)

    def test_sin_traccion(self):
        mat = ManderConcrete1D(fc=30e6, eps_c=0.002, E=27e9, eps_cu=0.004)
        s, t, _ = mat.compute_state(1e-4)
        self.assertEqual((s, t), (0.0, 0.0))
        _, t0, _ = mat.compute_state(0.0)
        self.assertEqual(t0, 27e9)


class TestRegularizacion(unittest.TestCase):

    def test_energia_de_traccion_igual_a_Gf_entre_h(self):
        G_f = 140.0
        mat = ManderConcrete1D(fc=30e6, eps_c=0.002, E=27.4e9, eps_cu=0.004, ft=3e6, G_f=G_f)
        self.assertTrue(mat.uses_characteristic_length)
        for h in (0.05, 0.1, 0.5):
            with self.subTest(h=h):
                st = mat.initial_point_state(h)
                eps_t = 3e6 / 27.4e9
                eps_tu = mat.tension_softening_end(h)
                self.assertAlmostEqual(eps_tu, 2 * G_f / (3e6 * h), delta=1e-15)
                path = np.array([0.0, eps_t, eps_tu, 1.5 * eps_tu])
                sig, _, states = walk(mat, path, st)
                W = float(np.sum(0.5 * (sig[1:] + sig[:-1]) * np.diff(path)))
                self.assertAlmostEqual(W * h / G_f, 1.0, delta=1e-10)
                self.assertEqual(sig[-1], 0.0)

    def test_limite_de_h_al_construir(self):
        mat = ManderConcrete1D(fc=30e6, eps_c=0.002, E=27.4e9, eps_cu=0.004, ft=3e6, G_f=140.0)
        h_max = 2 * 140.0 * 27.4e9 / 3e6 ** 2
        mat.initial_point_state(0.99 * h_max)
        with self.assertRaisesRegex(ValueError, "retroceso"):
            mat.initial_point_state(1.01 * h_max)
        # El elemento lo rechaza al construirse, con su nombre.
        nodes = [Node(1, [0.0, 0.0]), Node(2, [1.01 * h_max, 0.0])]
        with self.assertRaisesRegex(ValueError, r"Truss2D\(id=7\).*retroceso"):
            Truss2D(7, nodes, mat, A=0.01)

    def test_sin_longitud_caracteristica(self):
        mat = ManderConcrete1D(fc=30e6, eps_c=0.002, E=27.4e9, eps_cu=0.004, ft=3e6, G_f=140.0)
        # Antes de agrietar no hace falta h; al agrietar, error claro.
        mat.compute_state(1e-5)
        with self.assertRaisesRegex(ValueError, "longitud característica"):
            mat.compute_state(1e-3)

    def test_area_posterior_al_pico_igual_a_Gfc_entre_h(self):
        G_fc = 2.0e4
        for p in (COVER, CORE):
            q = {k: v for k, v in p.items() if k not in ("ft", "eps_tu")}
            mat = ManderConcrete1D(G_fc=G_fc, **q)
            for h in (0.05, 0.1, 0.5):
                with self.subTest(fc=p["fc"], h=h):
                    lam = mat._stretch(h)
                    e_end = mat.eps_c + lam * (mat.eps_sp - mat.eps_c)
                    e_cu = mat.eps_c + lam * (mat.eps_cu - mat.eps_c)
                    area = 0.0
                    xi, w = np.polynomial.legendre.leggauss(80)
                    for a, b in ((mat.eps_c, e_cu), (e_cu, e_end)):
                        if b > a:
                            pts = 0.5 * (b - a) * xi + 0.5 * (a + b)
                            area += 0.5 * (b - a) * sum(
                                wi * mat.compression_envelope(x, h)[0] for x, wi in zip(pts, w))
                    self.assertAlmostEqual(area * h / G_fc, 1.0, delta=1e-10)
                    # Rama ascendente sin cambio.
                    self.assertEqual(mat.compression_envelope(0.5 * mat.eps_c, h)[0],
                                     mat.compression_envelope(0.5 * mat.eps_c)[0])


class TestConfinamiento(unittest.TestCase):

    def test_ecuacion_29_con_presiones_iguales(self):
        for q in (0.0, 0.05, 0.1, 0.2, 0.3):
            self.assertAlmostEqual(confined_strength_ratio(30e6, q * 30e6), mander_K29(q),
                                   delta=1e-12)

    def test_chang_mander(self):
        for q1, q2 in ((0.05, 0.10), (0.03, 0.12), (0.10, 0.20), (0.2, 0.05)):
            self.assertAlmostEqual(confined_strength_ratio(30e6, q1 * 30e6, q2 * 30e6),
                                   chang_mander_K(q1, q2), delta=1e-12)
        for q in np.linspace(0.01, 0.3, 30):
            self.assertLess(abs(confined_strength_ratio(30e6, q * 30e6, q * 30e6)
                                - mander_K29(q)), 0.003)

    def test_ejemplo_numerico_de_mander(self):
        """Mander et al. (1988), p. 1812: f'co = 30 MPa, f'ly = 2.7 MPa,
        f'lx = 5.1 MPa ⇒ f'cc/f'co = 1.65 leído en el ábaco (Fig. 4)."""
        K = confined_strength_ratio(30e6, 5.1e6, 2.7e6)
        self.assertAlmostEqual(K, 1.65, delta=0.01)

    def test_parametros_del_nucleo(self):
        k_e = confinement_effectiveness_rectangular(
            b_c=0.4, d_c=0.3, s_clear=0.09, w_clear=[0.15, 0.15, 0.1, 0.1], rho_cc=0.02)
        ref = ((1 - (2 * 0.15 ** 2 + 2 * 0.1 ** 2) / (6 * 0.4 * 0.3))
               * (1 - 0.09 / 0.8) * (1 - 0.09 / 0.6) / 0.98)
        self.assertAlmostEqual(k_e, ref, delta=1e-14)
        out = mander_confined_concrete(30e6, f_yh=420e6, eps_su=0.09, k_e=k_e,
                                       rho_x=0.008, rho_y=0.012)
        K = chang_mander_K(k_e * 0.008 * 420 / 30, k_e * 0.012 * 420 / 30)
        self.assertAlmostEqual(out["fc"] / (K * 30e6), 1.0, delta=1e-12)
        self.assertAlmostEqual(out["eps_c"], 0.002 * (1 + 5 * (K - 1)), delta=1e-15)
        self.assertAlmostEqual(out["eps_cu"], 0.004 + 1.4 * 0.02 * 420e6 * 0.09 / out["fc"],
                               delta=1e-15)
        self.assertAlmostEqual(out["E"], 5000e6 * math.sqrt(30.0), delta=1e-3)
        # Circular con espiral y la geometría en un dict.
        geo = dict(s_clear=0.05, d_s=0.5, rho_cc=0.015, spiral=True)
        out_c = mander_confined_concrete(30e6, f_yh=420e6, eps_su=0.09, k_e=geo, rho_s=0.01)
        k_c = (1 - 0.05 / 1.0) / 0.985
        self.assertAlmostEqual(confinement_effectiveness_circular(**geo), k_c, delta=1e-14)
        self.assertAlmostEqual(out_c["fc"] / 30e6, mander_K29(k_c * 0.5 * 0.01 * 420 / 30),
                               delta=1e-12)
        # El material desde el dict de confinamiento.
        mat = ManderConcrete1D(confinement=dict(fc0=30e6, f_yh=420e6, eps_su=0.09, k_e=k_e,
                                                rho_x=0.008, rho_y=0.012))
        self.assertEqual((mat.fc, mat.eps_c, mat.eps_cu, mat.E),
                         (out["fc"], out["eps_c"], out["eps_cu"], out["E"]))
        with self.assertRaisesRegex(ValueError, "no se declaran"):
            ManderConcrete1D(fc=30e6, confinement=dict(fc0=30e6, f_yh=420e6, eps_su=0.09,
                                                       k_e=0.8, rho_s=0.01))

    def test_energia_de_fractura_del_model_code(self):
        self.assertAlmostEqual(fib_fracture_energy(38e6), 73.0 * 38.0 ** 0.18, delta=1e-12)
        self.assertAlmostEqual(fib_fracture_energy(38e3, MPa=1e3, N_per_m=1e-3),
                               73.0 * 38.0 ** 0.18 * 1e-3, delta=1e-15)

    def test_fuera_de_validez_avisa(self):
        with self.assertLogs("solidum", level=logging.WARNING) as cm:
            confined_strength_ratio(30e6, 0.35 * 30e6, 0.2 * 30e6)
        self.assertTrue(any("validez" in m for m in cm.output))


class TestTangenteConsistente(unittest.TestCase):

    def _check(self, mat, state, eps_points):
        for eps in eps_points:
            d = 1e-9
            sp = mat.compute_state(eps + d, dict(state) if state else None)[0]
            sm = mat.compute_state(eps - d, dict(state) if state else None)[0]
            _, tan, _ = mat.compute_state(eps, dict(state) if state else None)
            fd = (sp - sm) / (2 * d)
            scale = max(abs(tan), 1e-3 * mat.E)
            self.assertAlmostEqual(fd, tan, delta=1e-6 * mat.E + 1e-5 * scale,
                                   msg=f"ε = {eps}")

    def test_todas_las_ramas(self):
        cover = ManderConcrete1D(**COVER)
        reg = ManderConcrete1D(fc=30e6, eps_c=0.002, E=27.4e9, eps_cu=0.004, eps_sp=0.006,
                               ft=3e6, G_f=140.0, G_fc=2e4)
        for mat, st0 in ((cover, None), (reg, reg.initial_point_state(0.1))):
            with self.subTest(mat=repr(mat)):
                eps_t = mat.ft / mat.E
                # Envolvente antes y después del pico, rama lineal, tracción.
                self._check(mat, st0, [-0.5e-3, -1.5e-3, -2.6e-3, -3.5e-3, -4.6e-3,
                                       0.5 * eps_t, 1.5 * eps_t])
                # Recta de daño y secante agrietada.
                _, _, st = mat.compute_state(-3.0e-3, st0)
                ep = mat.plastic_strain(3.0e-3, mat.compression_envelope(3.0e-3, st["h"])[0])
                self._check(mat, st, [-2.0e-3, -0.5 * (3e-3 + ep), -ep + 0.5 * eps_t])
                _, _, st2 = mat.compute_state(-ep + 3 * eps_t, st)
                self._check(mat, st2, [-ep + 1.5 * eps_t, -ep + 2.5 * eps_t])


class TestDisipacionNoNegativa(unittest.TestCase):

    def test_trayectorias_aleatorias(self):
        rng = np.random.default_rng(20261109)
        reg = ManderConcrete1D(fc=30e6, eps_c=0.002, E=27.4e9, eps_cu=0.004, eps_sp=0.006,
                               ft=3e6, G_f=140.0, G_fc=2e4)
        cases = [(ManderConcrete1D(**COVER), None), (ManderConcrete1D(**CORE), None),
                 (ManderConcrete1D(fc=30e6, eps_c=0.002, E=27.4e9, eps_cu=0.004,
                                   eps_sp=0.006), None),
                 (reg, reg.initial_point_state(0.05)), (reg, reg.initial_point_state(0.5))]
        for mat, st0 in cases:
            scale = mat.fc * mat.eps_c
            worst = 0.0
            for _ in range(40):
                targets = np.concatenate([[0.0], rng.uniform(-3.0 * mat.eps_c, 0.002, 10)])
                path = np.concatenate([np.linspace(a, b, 300)[1:]
                                       for a, b in zip(targets[:-1], targets[1:])])
                state, W, s_prev, e_prev = (dict(st0) if st0 else None), 0.0, 0.0, 0.0
                for eps in path:
                    s, _, state = mat.compute_state(eps, state)
                    W += 0.5 * (s + s_prev) * (eps - e_prev)
                    s_prev, e_prev = s, eps
                    worst = min(worst, (W - stored_energy(mat, eps, s, state)) / scale)
            with self.subTest(mat=repr(mat), h=st0 and st0["h"]):
                self.assertGreater(worst, -1e-6)


if __name__ == "__main__":
    unittest.main()
