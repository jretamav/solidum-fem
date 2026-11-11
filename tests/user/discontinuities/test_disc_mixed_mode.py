"""``CST_Embedded2D`` con la ley cohesiva de modo mixto ``CohesiveDamageMixedMode``.

Casos ``elemento_embebido_*`` y ``registro_de_interfaz`` (YAML = Python) del
bloque ``acceptance`` de docs/specs/CohesiveDamageMixedMode.md: la ley dentro
del elemento embebido con cinemática KOS (Retama 2010).
"""
from __future__ import annotations

import copy
import os
import sys
import tempfile

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(__file__))

from solidum.core.node import Node  # noqa: E402
from solidum.materials.elastic_2d import Elastic2D  # noqa: E402
from solidum.materials.interface import CohesiveDamageMixedMode  # noqa: E402
from solidum.user.discontinuities.embedded_cst import CST_Embedded2D  # noqa: E402
from test_disc_embedded_uniaxial_softening import G_F, K_E, SIGMA_T0, run  # noqa: E402

# Triángulos de 5 cm: con la bilineal (|T'| = 0.8·σ_t0²/G_F) y β = 1.5 el
# elemento de 10 cm de los ensayos de modo I viola la condición de tamaño en
# cortante (test_condicion_de_tamano_en_cortante).
GEOMETRIES = (
    (0.5 * np.array([[0.0, 0.0], [0.1, 0.0], [0.1, 0.05]]), np.array([1.0, 0.0])),
    (0.5 * np.array([[0.0, 0.0], [0.1, 0.01], [0.03, 0.06]]), np.array([0.95, 0.31])),
)


@pytest.mark.parametrize("shape", ["linear", "exponential", "bilinear"])
def test_tangente_condensada_con_deslizamiento(shape):
    """acceptance.elemento_embebido_tangente_condensada: grieta abierta y
    deslizando en carga activa; la tangente condensada del elemento frente a
    diferencias finitas centradas de su fuerza interna (con el Newton local
    del salto reconvergido en cada perturbación)."""
    w_ref = 2.0 * 100.0 / 3.0e6
    for coords, n in GEOMETRIES:
        law = CohesiveDamageMixedMode(sigma_t0=3.0e6, G_f=100.0, K_e=1.0e15, beta=1.5,
                                      softening=shape)
        emb = CST_Embedded2D(1, [Node(i, list(coords[i])) for i in range(3)],
                             Elastic2D(E=30.0e9, nu=0.2, hypothesis="plane_strain"),
                             law, thickness=0.02)
        emb._activate(n / np.linalg.norm(n), coords)
        ds = emb.discontinuity_state
        for frac, slip in ((0.3, 0.2), (0.6, 0.5)):
            jump = np.array([frac * w_ref, slip * w_ref])
            u = np.full(6, 1.0e-6)
            i = ds.solitary_node
            u[2 * i:2 * i + 2] += jump[0] * ds.normal + jump[1] * ds.tangent
            d_eff = np.hypot(jump[0], 1.5 * jump[1])
            ds.cohesive_state_committed = {"kappa": 0.5 * d_eff}
            ds.jump_trial = jump.copy()

            def f_int(v):
                return copy.deepcopy(emb).compute_element_state(v)

            K, _ = f_int(u)
            h = 1.0e-7 * frac * w_ref
            K_fd = np.zeros((6, 6))
            for j in range(6):
                up, um = u.copy(), u.copy()
                up[j] += h
                um[j] -= h
                K_fd[:, j] = (f_int(up)[1] - f_int(um)[1]) / (2.0 * h)
            err = np.linalg.norm(K - K_fd) / np.linalg.norm(K_fd)
            assert err < 1e-5, (shape, tuple(n), frac, err)
            # La ley transmite cortante: el salto convergido desliza y la
            # tracción tangencial no es nula.
            probe = copy.deepcopy(emb)
            probe.compute_element_state(u)
            t, _, _ = law.compute_traction(probe.discontinuity_state.jump_trial,
                                           ds.cohesive_state_committed)
            assert abs(t[1]) > 1e-3 * 3.0e6


def test_condicion_de_tamano_en_cortante():
    """Con el salto despejado a ``u`` fijo, ``K_jj`` es definida positiva sólo
    si el ablandamiento no supera la rigidez del volumen en ninguna dirección.
    Grieta paralela al lado opuesto al nodo solitario, a distancia ``h``: el
    volumen aporta ``l_d·t·E'/h`` en la normal y ``l_d·t·G/h`` en el
    deslizamiento; la ley, ``l_d·t·T'`` y ``l_d·t·β²·T'`` en carga. La
    condición en cortante, ``β²·|T'|·h < G``, es ``β²·2(1+ν)`` veces más
    estricta que la de tracción con ``E'`` del orden de ``E``."""
    from solidum.elements.solid_2d import compute_kinematics_tri3

    E, nu, th, beta = 30.0e9, 0.2, 0.02, 1.5
    G_mod = E / (2.0 * (1.0 + nu))
    law = CohesiveDamageMixedMode(sigma_t0=3.0e6, G_f=100.0, K_e=1.0e15, beta=beta,
                                  softening="bilinear")
    w_1 = law._curve.w_1
    slope = 0.8 * 3.0e6 / w_1                                  # |T'| en la primera rama
    for h in (0.1, 0.05):
        coords = h / 0.1 * np.array([[0.0, 0.0], [0.1, 0.0], [0.1, 0.05]])
        emb = CST_Embedded2D(1, [Node(i, list(coords[i])) for i in range(3)],
                             Elastic2D(E=E, nu=nu, hypothesis="plane_strain"), law,
                             thickness=th)
        emb._activate(np.array([1.0, 0.0]), coords)
        ds = emb.discontinuity_state
        i = ds.solitary_node
        B, detJ = compute_kinematics_tri3(coords)
        G = np.column_stack([ds.normal, ds.tangent])
        B_phi = B[:, 2 * i:2 * i + 2]
        K_bulk = G.T @ (B_phi.T @ emb.material.C @ B_phi) @ G * 0.5 * detJ * th
        assert K_bulk[1, 1] == pytest.approx(ds.l_d * th * G_mod / h, rel=1e-12)
        # Deslizamiento puro en carga sobre la primera rama de la bilineal.
        jump = np.array([0.0, 0.5 * w_1 / beta])
        _, D, _ = law.compute_traction(jump, {"kappa": 0.25 * w_1})
        assert D[1, 1] == pytest.approx(-beta ** 2 * slope, rel=1e-9)
        K_jj = K_bulk + ds.l_d * th * D
        positive = np.linalg.eigvalsh(0.5 * (K_jj + K_jj.T)).min() > 0.0
        assert positive == (beta ** 2 * slope * h < G_mod)
        assert positive == (h == 0.05)


@pytest.mark.parametrize("shape", ["linear", "exponential"])
def test_limite_modo_I_en_traccion_uniaxial(shape):
    """acceptance.elemento_embebido_limite_modo_I: con β = 0 la probeta
    uniaxial hasta la separación completa da la misma respuesta que con
    CohesiveDamageIsotropic (cuya solución exacta comprueba
    test_disc_embedded_uniaxial_softening)."""
    reference, _ = run(shape)
    mixed = CohesiveDamageMixedMode(sigma_t0=SIGMA_T0, G_f=G_F, K_e=K_E, beta=0.0,
                                    softening=shape)
    history, _ = run(shape, cohesive=mixed)
    assert len(history) == len(reference)
    F_peak = SIGMA_T0 * 0.05 * 0.02
    for (d_mix, F_mix, c_mix), (d_ref, F_ref, c_ref) in zip(history, reference):
        assert d_mix == d_ref and c_mix == c_ref
        assert F_mix == pytest.approx(F_ref, abs=1e-10 * F_peak)


def test_yaml_igual_a_python():
    """acceptance.registro_de_interfaz: la ley declarada en interface_materials
    llega al elemento con los mismos parámetros que en Python."""
    from solidum.utils.yaml_parser import YamlParser

    yaml_text = """
nodes:
  - {id: 1, coords: [0.0, 0.0]}
  - {id: 2, coords: [1.0, 0.0]}
  - {id: 3, coords: [0.0, 1.0]}

materials:
  - {id: 1, type: Elastic2D, E: 30.0e9, nu: 0.2, hypothesis: plane_strain}

user_modules: [discontinuities]
interface_materials:
  - {id: 7, type: CohesiveDamageMixedMode, sigma_t0: 3.0e6, G_f: 100.0, K_e: 1.0e15,
     beta: 1.5, softening: bilinear, knee_stress_ratio: 0.3333333333333333,
     final_opening_factor: 3.6}

elements:
  - {id: 1, type: CST_Embedded2D, nodes: [1, 2, 3], material: 1, cohesive_material: 7}
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False,
                                     encoding="utf-8") as f:
        f.write(yaml_text)
        path = f.name
    try:
        domain = YamlParser(path).parse()
    finally:
        os.unlink(path)
    law = next(iter(domain.elements.values())).cohesive_material
    ref = CohesiveDamageMixedMode(sigma_t0=3.0e6, G_f=100.0, K_e=1.0e15, beta=1.5,
                                  softening="bilinear", knee_stress_ratio=1.0 / 3.0,
                                  final_opening_factor=3.6)
    assert type(law) is CohesiveDamageMixedMode
    rng = np.random.default_rng(3)
    s1 = s2 = None
    for _ in range(30):
        jump = rng.uniform(-0.3, 1.0, 2) * 1.0e-4
        t1, D1, s1 = law.compute_traction(jump, s1)
        t2, D2, s2 = ref.compute_traction(jump, s2)
        np.testing.assert_array_equal(t1, t2)
        np.testing.assert_array_equal(D1, D2)
