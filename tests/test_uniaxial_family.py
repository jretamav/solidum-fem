"""Prueba común de la familia de materiales uniaxiales (spec ``ManderConcrete1D``,
caso ``familia_uniaxial``).

Recorre **todo** material 1D del programa principal con su muestra de
``MUESTRAS`` (``test_material_contract_sweep``, que ya obliga a cubrir cada
material nuevo): ningún modelo uniaxial —en particular, los histeréticos por
reglas, cuya construcción no garantiza la disipación— entra sin pasarla.

1. **Disipación no negativa.** (a) Desde el estado virgen, el trabajo
   realizado ``W(t) = ∫σ dε`` es la energía almacenada más la disipada, las
   dos ≥ 0: ``W(t) ≥ 0`` en todo instante de trayectorias cíclicas
   aleatorias. (b) En un ciclo cerrado de deformación repetido hasta que la
   respuesta se estabiliza (el ciclo repite exactamente los esfuerzos del
   anterior), el trabajo del ciclo no puede ser negativo: si lo fuera, al
   repetirlo ``W(t)`` decrecería sin límite, y ``W = ψ + D`` con ``ψ ≥ 0`` y
   ``D ≥ 0`` no puede hacerlo (postulado de Ilyushin).
2. **Tangente consistente.** Fuera de los quiebres (donde las derivadas por
   la izquierda y por la derecha difieren), la tangente coincide con la
   derivada por diferencias finitas centradas del esfuerzo, evaluada con el
   estado comprometido.
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.dirname(__file__)))

import solidum  # noqa: F401  — autodiscover
from solidum.registry import MaterialRegistry

from _programa_principal import nombres_programa_principal
from test_material_contract_sweep import MUESTRAS


def _uniaxiales():
    for nombre in nombres_programa_principal(MaterialRegistry):
        material = MaterialRegistry.create(nombre, **MUESTRAS[nombre])
        if material.STRAIN_DIM == 1:
            yield nombre, material


def _escala_de_deformacion(material):
    """Amplitud de los ciclos: varias veces la deformación característica
    del material (fluencia, umbral de daño, pico), o 1e-3 si es elástico."""
    E = float(np.asarray(material.compute_state(0.0, None)[1]).reshape(-1)[0])
    for attr in ("sigma_y", "fc"):
        if getattr(material, attr, None):
            return 4.0 * getattr(material, attr) / E
    if getattr(material, "kappa_0", None):
        return 20.0 * material.kappa_0
    return 1.0e-3


def _rigidez(material):
    """Rigidez de referencia: la mayor tangente inicial a uno y otro lado de
    ε = 0 (un material unilateral no tiene rigidez en uno de ellos)."""
    amp = 1e-9
    return max(abs(_sigma(material, d, None)[1]) for d in (amp, -amp))


def _ciclo(material, state, puntos, n=400):
    """Recorre los tramos rectos entre ``puntos``; devuelve (W, estado,
    esfuerzos del recorrido)."""
    W = 0.0
    s_prev = _sigma(material, puntos[0], state)[0]
    e_prev = puntos[0]
    sig = []
    for a, b in zip(puntos[:-1], puntos[1:]):
        for eps in np.linspace(a, b, n)[1:]:
            s, _, state = _sigma(material, eps, state)
            W += 0.5 * (s + s_prev) * (eps - e_prev)
            s_prev, e_prev = s, eps
            sig.append(s)
    return W, state, np.array(sig)


def _sigma(material, eps, state):
    s, t, st = material.compute_state(eps, None if state is None else dict(state))
    return float(np.asarray(s).reshape(-1)[0]), float(np.asarray(t).reshape(-1)[0]), st


class TestFamiliaUniaxial(unittest.TestCase):

    def test_hay_materiales_uniaxiales(self):
        self.assertGreaterEqual(len(list(_uniaxiales())), 5)

    def test_disipacion_no_negativa_desde_el_estado_virgen(self):
        rng = np.random.default_rng(20261109)
        for nombre, material in _uniaxiales():
            amp = _escala_de_deformacion(material)
            E = _rigidez(material)
            with self.subTest(material=nombre):
                worst = 0.0
                for _ in range(20):
                    targets = np.concatenate([[0.0], rng.uniform(-amp, amp, 8)])
                    path = np.concatenate([np.linspace(a, b, 400)[1:]
                                           for a, b in zip(targets[:-1], targets[1:])])
                    state, W, s_prev, e_prev = None, 0.0, 0.0, 0.0
                    for eps in path:
                        s, _, state = _sigma(material, eps, state)
                        W += 0.5 * (s + s_prev) * (eps - e_prev)
                        s_prev, e_prev = s, eps
                        worst = min(worst, W)
                # Tolerancia: el error de la regla del trapecio en los quiebres.
                self.assertGreater(worst, -1e-4 * E * amp ** 2,
                                   f"{nombre}: el trabajo realizado desde el estado virgen "
                                   f"se hizo negativo (crea energía)")

    def test_ciclos_cerrados_estabilizados(self):
        for nombre, material in _uniaxiales():
            amp = _escala_de_deformacion(material)
            E = _rigidez(material)
            n_estables = 0
            for lo, hi in ((-amp, amp), (-0.5 * amp, 0.25 * amp), (-amp, -0.3 * amp),
                           (0.0, 0.8 * amp)):
                with self.subTest(material=nombre, ciclo=(lo / amp, hi / amp)):
                    puntos = [0.0, hi, lo, 0.0]
                    state = None
                    for _ in range(4):                    # estabilizar la respuesta
                        _, state, sig_prev = _ciclo(material, state, puntos)
                    W, state, sig = _ciclo(material, state, puntos)
                    escala = max(np.abs(sig).max(), 1e-30)
                    if not np.allclose(sig, sig_prev, rtol=0.0, atol=1e-9 * escala):
                        continue          # la respuesta aún evoluciona
                    n_estables += 1
                    self.assertGreater(W, -1e-4 * E * amp ** 2,
                                       f"{nombre}: el ciclo cerrado crea energía")
            with self.subTest(material=nombre):
                self.assertGreater(n_estables, 0, f"{nombre}: ningún ciclo se estabilizó")

    def test_tangente_consistente_fuera_de_los_quiebres(self):
        rng = np.random.default_rng(7)
        for nombre, material in _uniaxiales():
            amp = _escala_de_deformacion(material)
            E = _rigidez(material)
            with self.subTest(material=nombre):
                n_checked = 0
                for _ in range(6):
                    # Historia previa aleatoria, comprometida punto a punto.
                    state = None
                    for eps in np.linspace(0.0, rng.uniform(-amp, amp), 50)[1:]:
                        state = _sigma(material, eps, state)[2]
                    for eps in rng.uniform(-amp, amp, 25):
                        d = 1e-7 * amp
                        s_p = _sigma(material, eps + d, state)[0]
                        s_0, tan, _ = _sigma(material, eps, state)
                        s_m = _sigma(material, eps - d, state)[0]
                        right, left = (s_p - s_0) / d, (s_0 - s_m) / d
                        if abs(right - left) > 1e-4 * E:
                            continue                      # quiebre: no derivable
                        n_checked += 1
                        self.assertAlmostEqual(
                            0.5 * (right + left), tan, delta=1e-5 * E,
                            msg=f"{nombre}: tangente inconsistente en ε = {eps:.6g}")
                self.assertGreater(n_checked, 50)


if __name__ == "__main__":
    unittest.main()
