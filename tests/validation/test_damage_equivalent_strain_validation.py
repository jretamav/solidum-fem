"""Validación ciega de la deformación equivalente configurable del daño isótropo.

Materiales: ``IsotropicDamage2D`` (esfuerzo plano ``plane_stress`` y
deformación plana ``plane_strain``) e ``IsotropicDamage3D``, con
``equivalent_strain ∈ {"energy", "mazars", "modified_von_mises"}`` y
``k = f_c/f_t`` obligatorio (y sólo admitido) con ``"modified_von_mises"``.

Este archivo se escribió sólo a partir de la física aprobada por el usuario
el 2026-09-29, **sin leer la implementación**: es una comprobación
independiente. Toda referencia numérica sale de fórmulas cerradas o de una
implementación de referencia propia, en numpy puro, de las tres medidas
(``_eps_eq_reference``), escrita desde las definiciones de abajo.

Modelo común a las tres medidas
-------------------------------
.. math::

   \\boldsymbol\\sigma = (1-d)\\,\\mathbf C_e\\,\\boldsymbol\\varepsilon, \\qquad
   \\kappa = \\max(\\kappa_0,\\ \\max_{hist}\\varepsilon_{eq}), \\qquad
   d(\\kappa) = 1 - \\frac{\\kappa_0}{\\kappa}\\,e^{-\\alpha(\\kappa-\\kappa_0)}
   \\ \\ (\\kappa > \\kappa_0)

con ``d = 0`` si ``κ ≤ κ₀`` y ``d ≤ DAMAGE_MAX``. Tangente de carga activa
con daño no saturado:
``C_alg = (1−d)·C_e − (∂d/∂κ)·(C_e·ε) ⊗ (∂ε_eq/∂ε)``,
``∂d/∂κ = (1−d)(1/κ + α)``. En descarga o saturación, ``C = (1−d)·C_e``.

Medidas de deformación equivalente
----------------------------------
1. ``energy`` (Simó y Ju 1987): ``ε_eq = √(εᵀ·C_e·ε / E)``, con ``C_e`` la
   matriz Voigt del propio material (en 2D la de esfuerzo plano o la de
   deformación plana).
2. ``mazars`` (Mazars 1986): ``ε_eq = √(Σ_i ⟨ε_i⟩₊²)``, con ``ε_i`` las
   deformaciones principales del tensor 3D completo.
3. ``modified_von_mises`` (de Vree, Brekelmans y van Gils 1995):
   ``ε_eq = (k−1)/(2k(1−2ν))·I1 + 1/(2k)·√( ((k−1)/(1−2ν))²·I1² + 12k/(1+ν)²·J2 )``,
   ``I1 = εxx+εyy+εzz``,
   ``J2 = (1/6)[(εxx−εyy)² + (εyy−εzz)² + (εzz−εxx)²] + (γxy²+γyz²+γxz²)/4``.

Deformación fuera del plano en 2D: esfuerzo plano
``ε_zz = −ν/(1−ν)·(εxx+εyy)`` (con daño escalar isótropo ``σ_zz = 0`` no
depende de ``d``); deformación plana ``ε_zz = 0``. En 2D ``z`` es dirección
principal (``γxz = γyz = 0``).

Voigt del proyecto: ``ε = [εxx, εyy, γxy]`` (2D) y
``[εxx, εyy, εzz, γxy, γyz, γxz]`` (3D), ``γ = 2ε_ij``.

Formas cerradas (E = 30 GPa, ν = 0.2, κ₀ = 1e-4, α = 1e3, k = 10)
------------------------------------------------------------------
- Tracción uniaxial de esfuerzo ``ε = [ε, −νε, −νε]``: las tres dan
  ``ε_eq = ε`` ⇒ el daño empieza en ``ε = κ₀ = 1e-4``.
- Compresión uniaxial de esfuerzo ``ε = [−e, νe, νe]``:
  energy ``ε_eq = e`` (umbral ``e = 1e-4``); mazars ``ε_eq = √2·ν·e``
  (umbral ``e = κ₀/(√2ν) = 3.5355e-4``); modified_von_mises ``ε_eq = e/k``
  (umbral ``e = k·κ₀ = 1e-3``).
- Cortante puro ``γxy = γ`` (principales ``±γ/2, 0``; ``I1 = 0``,
  ``J2 = γ²/4``): energy ``ε_eq = γ/√(2(1+ν))`` (umbral ``γ = 1.5492e-4``);
  mazars ``ε_eq = γ/2`` (umbral ``γ = 2e-4``); modified_von_mises
  ``ε_eq = √(3/k)·γ/(2(1+ν))`` (umbral ``γ = 4.3818e-4``).
- Deformación hidrostática ``ε = ±e·[1,1,1,0,0,0]``: energy
  ``ε_eq = √(3/(1−2ν))·e = √5·e`` con ambos signos; mazars ``√3·e`` en
  expansión y ``0`` en contracción; modified_von_mises
  ``3(k−1)/(k(1−2ν))·e = 4.5·e`` en expansión y ``0`` en contracción.
- Ley de daño: ``d(1.5κ₀) = 0.365847``, ``d(2κ₀) = 0.547581``.

Referencias
-----------
- Simó J.C., Ju J.W. (1987). Strain- and stress-based continuum damage
  models — I. Formulation. *Int. J. Solids Struct.* 23, 821-840.
- Mazars J. (1986). A description of micro- and macroscale damage of
  concrete structures. *Eng. Fract. Mech.* 25, 729-737.
- de Vree J.H.P., Brekelmans W.A.M., van Gils M.A.J. (1995). Comparison of
  nonlocal approaches in continuum damage mechanics. *Comput. Struct.* 55,
  581-588.
"""
import math
import os
import sys

import numpy as np
import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from solidum.constants import DAMAGE_MAX
from solidum.materials.damage_2d import IsotropicDamage2D
from solidum.materials.damage_3d import IsotropicDamage3D


# =============================================================================
# Datos físicos (valores distintos de 1 para que un factor olvidado se note).
# =============================================================================

E = 30.0e9
NU = 0.2
KAPPA_0 = 1.0e-4
ALPHA = 1.0e3
K_RATIO = 10.0              # k = f_c / f_t de modified_von_mises

MEASURES = ("energy", "mazars", "modified_von_mises")
HYPOTHESES_2D = ("plane_stress", "plane_strain")
CONFIGS = ("plane_stress", "plane_strain", "3d")

LAMBDA = E * NU / ((1.0 + NU) * (1.0 - 2.0 * NU))
MU = E / (2.0 * (1.0 + NU))
C_ZZ = -NU / (1.0 - NU)     # ε_zz = C_ZZ·(εxx + εyy) en esfuerzo plano

# Matrices elásticas Voigt del proyecto (cortante engineering: σ_xy = μ·γ_xy).
C_3D = np.array([
    [LAMBDA + 2 * MU, LAMBDA, LAMBDA, 0.0, 0.0, 0.0],
    [LAMBDA, LAMBDA + 2 * MU, LAMBDA, 0.0, 0.0, 0.0],
    [LAMBDA, LAMBDA, LAMBDA + 2 * MU, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, MU, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, MU, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, MU],
])
C_PS = E / (1.0 - NU ** 2) * np.array([
    [1.0, NU, 0.0],
    [NU, 1.0, 0.0],
    [0.0, 0.0, 0.5 * (1.0 - NU)],
])
C_PE = np.array([
    [LAMBDA + 2 * MU, LAMBDA, 0.0],
    [LAMBDA, LAMBDA + 2 * MU, 0.0],
    [0.0, 0.0, MU],
])
C_ELASTIC = {"plane_stress": C_PS, "plane_strain": C_PE, "3d": C_3D}

# Factores cerrados de ε_eq (ver docstring del módulo).
SHEAR_FACTOR = {                          # ε_eq = factor·|γ| en cortante puro
    "energy": 1.0 / math.sqrt(2.0 * (1.0 + NU)),
    "mazars": 0.5,
    "modified_von_mises": math.sqrt(3.0 / K_RATIO) / (2.0 * (1.0 + NU)),
}
COMPRESSION_FACTOR = {                    # ε_eq = factor·e en compresión uniaxial
    "energy": 1.0,
    "mazars": math.sqrt(2.0) * NU,
    "modified_von_mises": 1.0 / K_RATIO,
}
HYDRO_EXPANSION_FACTOR = {                # ε_eq = factor·e con ε = e·[1,1,1]
    "energy": math.sqrt(3.0 / (1.0 - 2.0 * NU)),
    "mazars": math.sqrt(3.0),
    "modified_von_mises": 3.0 * (K_RATIO - 1.0) / (K_RATIO * (1.0 - 2.0 * NU)),
}
HYDRO_CONTRACTION_FACTOR = {              # ε_eq = factor·e con ε = −e·[1,1,1]
    "energy": math.sqrt(3.0 / (1.0 - 2.0 * NU)),
    "mazars": 0.0,
    "modified_von_mises": 0.0,
}


# =============================================================================
# Referencia propia (numpy puro) — no usa ningún código de solidum.
# =============================================================================

def _damage_law(kappa):
    """Ley exponencial acotada por DAMAGE_MAX (0 si κ ≤ κ₀)."""
    if kappa <= KAPPA_0:
        return 0.0
    d = 1.0 - (KAPPA_0 / kappa) * math.exp(-ALPHA * (kappa - KAPPA_0))
    return min(d, DAMAGE_MAX)


def _embed_3d(eps2, hypothesis):
    """Estado 3D (Voigt 6) de un estado 2D: ε_zz según la hipótesis."""
    exx, eyy, gxy = (float(v) for v in eps2)
    ezz = C_ZZ * (exx + eyy) if hypothesis == "plane_stress" else 0.0
    return np.array([exx, eyy, ezz, gxy, 0.0, 0.0])


def _tensor_from_voigt(eps6):
    """Tensor 3×3 de deformaciones desde Voigt engineering (ε_ij = γ_ij/2)."""
    exx, eyy, ezz, gxy, gyz, gxz = (float(v) for v in eps6)
    return np.array([
        [exx, 0.5 * gxy, 0.5 * gxz],
        [0.5 * gxy, eyy, 0.5 * gyz],
        [0.5 * gxz, 0.5 * gyz, ezz],
    ])


def _voigt_from_tensor(t):
    """Voigt engineering desde el tensor 3×3 (γ_ij = 2ε_ij)."""
    return np.array([t[0, 0], t[1, 1], t[2, 2],
                     2.0 * t[0, 1], 2.0 * t[1, 2], 2.0 * t[0, 2]])


def _eps_eq_mazars_3d(eps6):
    principal = np.linalg.eigvalsh(_tensor_from_voigt(eps6))
    positive = np.maximum(principal, 0.0)
    return float(np.sqrt(np.sum(positive ** 2)))


def _eps_eq_mvm_3d(eps6, k=K_RATIO):
    exx, eyy, ezz, gxy, gyz, gxz = (float(v) for v in eps6)
    i1 = exx + eyy + ezz
    j2 = ((exx - eyy) ** 2 + (eyy - ezz) ** 2 + (ezz - exx) ** 2) / 6.0 \
        + (gxy ** 2 + gyz ** 2 + gxz ** 2) / 4.0
    a = (k - 1.0) / (2.0 * k * (1.0 - 2.0 * NU))
    b = (k - 1.0) / (1.0 - 2.0 * NU)
    c = 12.0 * k / (1.0 + NU) ** 2
    return a * i1 + math.sqrt(b ** 2 * i1 ** 2 + c * j2) / (2.0 * k)


def _eps_eq_reference(eps, measure, config):
    """ε_eq de referencia para un estado Voigt de la configuración dada.

    ``energy`` en 2D se evalúa con la matriz 2D (esfuerzo plano o
    deformación plana) directamente, no a través del 3D; ``mazars`` y
    ``modified_von_mises`` en 2D se evalúan sobre el tensor 3D completo con
    ``ε_zz`` de la hipótesis.
    """
    eps = np.asarray(eps, dtype=float)
    if measure == "energy":
        C = C_ELASTIC[config]
        return math.sqrt(max(float(eps @ C @ eps), 0.0) / E)
    eps6 = eps if config == "3d" else _embed_3d(eps, config)
    if measure == "mazars":
        return _eps_eq_mazars_3d(eps6)
    if measure == "modified_von_mises":
        return _eps_eq_mvm_3d(eps6)
    raise ValueError(measure)


def _normalized(direction, measure, config):
    """Escala ``direction`` para que su ε_eq de referencia valga κ₀.

    Las tres medidas son positivamente homogéneas de grado 1
    (ε_eq(λε) = λ·ε_eq(ε), λ > 0), así que ``s·_normalized(...)`` tiene
    ``ε_eq = s·κ₀``.
    """
    direction = np.asarray(direction, dtype=float)
    eq = _eps_eq_reference(direction, measure, config)
    assert eq > 0.0, "la dirección debe tener ε_eq > 0 con esta medida"
    return direction * (KAPPA_0 / eq)


# =============================================================================
# Utilidades sobre el material bajo prueba.
# =============================================================================

def _make(config, measure):
    kwargs = dict(E=E, nu=NU, kappa_0=KAPPA_0, alpha=ALPHA, equivalent_strain=measure)
    if measure == "modified_von_mises":
        kwargs["k"] = K_RATIO
    if config == "3d":
        return IsotropicDamage3D(**kwargs)
    return IsotropicDamage2D(hypothesis=config, **kwargs)


def _call(mat, eps, state):
    """compute_state sobre una copia del estado (el de entrada queda intacto)."""
    sigma, C_tan, new_state = mat.compute_state(np.asarray(eps, dtype=float), dict(state))
    return np.asarray(sigma, dtype=float), np.asarray(C_tan, dtype=float), new_state


def _fresh(mat, eps):
    """Evaluación desde el estado inicial (carga monótona desde cero)."""
    return _call(mat, eps, mat.initial_state())


def _uniaxial_stress_strain(eps_xx, config):
    """Deformación de un estado uniaxial de esfuerzo σ = σ_xx·e_x.

    3D: ``[ε, −νε, −νε, 0, 0, 0]``; 2D esfuerzo plano: ``[ε, −νε, 0]``
    (su ``ε_zz = −ν/(1−ν)·(1−ν)ε = −νε`` coincide con la del 3D).
    """
    if config == "3d":
        return np.array([eps_xx, -NU * eps_xx, -NU * eps_xx, 0.0, 0.0, 0.0])
    assert config == "plane_stress"
    return np.array([eps_xx, -NU * eps_xx, 0.0])


# =============================================================================
# 0. Autoconsistencia de la referencia del propio test (no usa el material).
# =============================================================================

@pytest.mark.parametrize("measure", MEASURES)
def test_reference_reproduces_closed_forms(measure):
    """La referencia numpy reproduce las formas cerradas del docstring.

    Blinda el propio test: si una fórmula cerrada o la referencia estuviese
    mal derivada, los tests de abajo compararían contra un valor erróneo.
    Tracción uniaxial ⇒ ε_eq = ε (las tres); compresión uniaxial ⇒
    COMPRESSION_FACTOR·e; cortante puro ⇒ SHEAR_FACTOR·|γ|; hidrostática ⇒
    HYDRO_*_FACTOR·e. En 2D, energy con C_PS debe coincidir con energy 3D
    sobre el estado embebido (εᵀC_ps ε = σ:ε del 3D con σ_zz = 0).
    """
    e = 3.0e-4
    for config in ("3d", "plane_stress"):
        tens = _uniaxial_stress_strain(e, config)
        comp = _uniaxial_stress_strain(-e, config)
        assert _eps_eq_reference(tens, measure, config) == pytest.approx(e, rel=1e-13)
        assert _eps_eq_reference(comp, measure, config) == pytest.approx(
            COMPRESSION_FACTOR[measure] * e, rel=1e-13)
    for gamma in (e, -e):
        assert _eps_eq_reference([0, 0, 0, gamma, 0, 0], measure, "3d") == pytest.approx(
            SHEAR_FACTOR[measure] * e, rel=1e-13)
    hydro = e * np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
    assert _eps_eq_reference(hydro, measure, "3d") == pytest.approx(
        HYDRO_EXPANSION_FACTOR[measure] * e, rel=1e-13)
    assert _eps_eq_reference(-hydro, measure, "3d") == pytest.approx(
        HYDRO_CONTRACTION_FACTOR[measure] * e, rel=1e-13, abs=1e-20)
    eps2 = np.array([1.3e-4, -0.7e-4, 0.9e-4])
    for hyp in HYPOTHESES_2D:
        e3 = _embed_3d(eps2, hyp)
        assert float(eps2 @ C_ELASTIC[hyp] @ eps2) == pytest.approx(float(e3 @ C_3D @ e3), rel=1e-13)


# =============================================================================
# (a) Tracción uniaxial de esfuerzo: umbral κ₀ y curva monótona.
# =============================================================================

@pytest.mark.parametrize("config", ["3d", "plane_stress"])
@pytest.mark.parametrize("measure", MEASURES)
def test_uniaxial_tension_onset_at_kappa0(measure, config):
    """Tracción uniaxial de esfuerzo: el daño empieza exactamente en ε = κ₀.

    Con ``ε = [ε, −νε, −νε]``:
    - energy: ``C_e·ε = [Eε, 0, 0]`` ⇒ ``εᵀC_eε/E = ε²`` ⇒ ``ε_eq = ε``.
    - mazars: la única principal positiva es ``ε`` ⇒ ``ε_eq = ε``.
    - modified_von_mises: ``I1 = (1−2ν)ε``, ``J2 = (1+ν)²ε²/3`` ⇒ la raíz vale
      ``(k+1)ε`` y ``ε_eq = [(k−1) + (k+1)]ε/(2k) = ε``.
    Por tanto ``d = 0`` en ``0.999κ₀`` y ``d > 0`` en ``1.001κ₀``
    (``d(1.001κ₀) = 1 − e^{−1e-4}/1.001 ≈ 1.10e-3``). El margen del 0.1 %
    supera en muchos órdenes el redondeo de ε_eq (~1e-15 relativo).
    """
    mat = _make(config, measure)
    _, _, st = _fresh(mat, _uniaxial_stress_strain(0.999 * KAPPA_0, config))
    assert st["damage"] == 0.0
    _, _, st = _fresh(mat, _uniaxial_stress_strain(1.001 * KAPPA_0, config))
    assert st["damage"] > 0.0
    assert st["kappa"] == pytest.approx(1.001 * KAPPA_0, rel=1e-12)
    assert st["damage"] == pytest.approx(_damage_law(1.001 * KAPPA_0), rel=1e-9)


@pytest.mark.parametrize("config", ["3d", "plane_stress"])
@pytest.mark.parametrize("measure", MEASURES)
def test_uniaxial_tension_monotonic_curve(measure, config):
    """Curva monótona en tracción uniaxial de esfuerzo con el estado arrastrado.

    Como ``ε_eq = ε`` (ver test anterior) y la carga es monótona, ``κ = ε`` y
    ``d(ε) = 1 − (κ₀/ε)·e^{−α(ε−κ₀)}`` (acotado a DAMAGE_MAX); además
    ``C_e·ε = [Eε, 0, 0, …]`` ⇒ ``σ_xx = (1−d)·E·ε`` y el resto de
    componentes nulas. El barrido llega a la saturación: en 45κ₀,
    ``(1/45)·e^{−4.4} = 2.7e-4 < 1 − DAMAGE_MAX``.

    Tolerancias: ``κ`` y ``d`` a 1e-12 relativo (redondeo de ε_eq ~1e-15,
    amplificado por ``(1−d)(1+ακ)/d`` = O(1) lejos del umbral); ``σ_xx`` a
    1e-10 relativo porque en la saturación ``1−d = 1e-3`` amplifica 1e3 el
    redondeo de ``d``; componentes nulas a 1e-12·E·ε (redondeo de
    ``λε − (λ+2μ)νε − λνε`` ~ ulp·λ·ε).
    """
    mat = _make(config, measure)
    st = mat.initial_state()
    sweep = np.concatenate([np.linspace(1.05, 20.0, 40), [30.0, 45.0, 60.0]]) * KAPPA_0
    saturated = False
    for eps_xx in sweep:
        eps = _uniaxial_stress_strain(float(eps_xx), config)
        sigma, _, st = _call(mat, eps, st)
        d_ref = _damage_law(float(eps_xx))
        assert st["kappa"] == pytest.approx(eps_xx, rel=1e-12)
        assert st["damage"] == pytest.approx(d_ref, rel=1e-12)
        assert sigma[0] == pytest.approx((1.0 - d_ref) * E * eps_xx, rel=1e-10)
        np.testing.assert_allclose(sigma[1:], 0.0, atol=1e-12 * E * eps_xx)
        saturated |= d_ref == DAMAGE_MAX
    assert saturated, "el barrido debe alcanzar DAMAGE_MAX"


# =============================================================================
# (b) Compresión uniaxial de esfuerzo: umbral propio de cada medida.
# =============================================================================

@pytest.mark.parametrize("config", ["3d", "plane_stress"])
@pytest.mark.parametrize("measure", MEASURES)
def test_uniaxial_compression_onset(measure, config):
    """Compresión uniaxial de esfuerzo ``ε = [−e, νe, νe]``.

    - energy: ``εᵀC_eε = σ:ε = E·e²`` ⇒ ``ε_eq = e``; umbral ``e = κ₀ = 1e-4``.
    - mazars: principales positivas ``νe, νe`` ⇒ ``ε_eq = √2·ν·e``; umbral
      ``e = κ₀/(√2ν) = 3.5355e-4``.
    - modified_von_mises: ``I1 = −(1−2ν)e``, raíz ``(k+1)e`` ⇒
      ``ε_eq = [−(k−1) + (k+1)]e/(2k) = e/k``; umbral ``e = kκ₀ = 1e-3``.

    Se comprueba el umbral por ambos lados (0.999 y 1.001 veces) y, en
    ``e = 1.5·umbral``: ``κ = 1.5κ₀``, ``d = d(1.5κ₀) = 0.365847`` y
    ``σ_xx = −(1−d)·E·e`` con el resto de componentes nulas.
    """
    mat = _make(config, measure)
    e_crit = KAPPA_0 / COMPRESSION_FACTOR[measure]
    _, _, st = _fresh(mat, _uniaxial_stress_strain(-0.999 * e_crit, config))
    assert st["damage"] == 0.0
    _, _, st = _fresh(mat, _uniaxial_stress_strain(-1.001 * e_crit, config))
    assert st["damage"] > 0.0

    e = 1.5 * e_crit
    sigma, _, st = _fresh(mat, _uniaxial_stress_strain(-e, config))
    d_ref = _damage_law(1.5 * KAPPA_0)
    assert d_ref == pytest.approx(0.365847, abs=1e-6)
    assert st["kappa"] == pytest.approx(1.5 * KAPPA_0, rel=1e-12)
    assert st["damage"] == pytest.approx(d_ref, rel=1e-12)
    assert sigma[0] == pytest.approx(-(1.0 - d_ref) * E * e, rel=1e-10)
    np.testing.assert_allclose(sigma[1:], 0.0, atol=1e-12 * E * e)


# =============================================================================
# (c) Cortante puro.
# =============================================================================

def _pure_shear_cases():
    """(config, índice Voigt del cortante): los tres cortantes del 3D y γxy en 2D."""
    return [("3d", 3), ("3d", 4), ("3d", 5), ("plane_stress", 2), ("plane_strain", 2)]


@pytest.mark.parametrize("config,slot", _pure_shear_cases())
@pytest.mark.parametrize("measure", MEASURES)
def test_pure_shear(measure, config, slot):
    """Cortante puro ``γ`` en un único componente cortante.

    Principales ``+γ/2, −γ/2, 0``; ``I1 = 0``; ``J2 = γ²/4``; en 2D
    ``εxx + εyy = 0`` ⇒ ``ε_zz = 0`` en ambas hipótesis, así que 2D y 3D dan
    lo mismo.
    - energy: ``εᵀC_eε = μγ²`` ⇒ ``ε_eq = √(μ/E)·|γ| = |γ|/√(2(1+ν))``
      = 0.645497·|γ|; umbral ``γ = 1.549193e-4``.
    - mazars: ``ε_eq = |γ|/2``; umbral ``γ = 2e-4``.
    - modified_von_mises: ``ε_eq = (1/2k)·√(12k/(1+ν)²·γ²/4)
      = √(3/k)·|γ|/(2(1+ν))`` = 0.228218·|γ|; umbral ``γ = 4.381780e-4``.

    Ambos signos de γ (las medidas son pares en el cortante) y, en 3D, los
    tres componentes cortantes (isotropía; detecta un factor 2 o un orden
    Voigt equivocado en un solo componente).
    """
    mat = _make(config, measure)
    n = 6 if config == "3d" else 3
    gamma_crit = KAPPA_0 / SHEAR_FACTOR[measure]
    for sign in (1.0, -1.0):
        def shear(g):
            eps = np.zeros(n)
            eps[slot] = sign * g
            return eps
        _, _, st = _fresh(mat, shear(0.999 * gamma_crit))
        assert st["damage"] == 0.0
        _, _, st = _fresh(mat, shear(1.001 * gamma_crit))
        assert st["damage"] > 0.0
        sigma, _, st = _fresh(mat, shear(2.0 * gamma_crit))
        d_ref = _damage_law(2.0 * KAPPA_0)
        assert st["kappa"] == pytest.approx(2.0 * KAPPA_0, rel=1e-12)
        assert st["damage"] == pytest.approx(d_ref, rel=1e-12)
        expected = np.zeros(n)
        expected[slot] = (1.0 - d_ref) * MU * sign * 2.0 * gamma_crit
        np.testing.assert_allclose(sigma, expected, rtol=1e-10, atol=1e-12 * E * gamma_crit)


# =============================================================================
# Hidrostática: la distinción tracción/compresión de mazars y modified_von_mises.
# =============================================================================

@pytest.mark.parametrize("measure", MEASURES)
def test_hydrostatic_expansion_and_contraction_3d(measure):
    """Deformación hidrostática ``ε = ±e·[1,1,1,0,0,0]`` en 3D.

    - energy: ``εᵀC_eε = 9K·e²`` con ``K = E/(3(1−2ν))`` ⇒
      ``ε_eq = √(3/(1−2ν))·e = √5·e`` con ambos signos (medida simétrica).
    - mazars: expansión ``√3·e``; contracción: ninguna principal positiva ⇒
      ``ε_eq = 0`` ⇒ no hay daño por grande que sea ``e``.
    - modified_von_mises: ``J2 = 0`` ⇒ ``ε_eq = a·I1 + a·|I1|`` con
      ``a = (k−1)/(2k(1−2ν))``: expansión ``6a·e = 4.5·e``; contracción
      ``ε_eq = 0`` (para ``k ≥ 1``).
    Umbral de expansión en ``e = κ₀/factor``; en contracción, energy daña en
    ``e = κ₀/√5`` y las otras dos no dañan ni con ``e = 100κ₀``.
    """
    mat = _make("3d", measure)
    ones = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])

    e_crit = KAPPA_0 / HYDRO_EXPANSION_FACTOR[measure]
    _, _, st = _fresh(mat, 0.999 * e_crit * ones)
    assert st["damage"] == 0.0
    _, _, st = _fresh(mat, 1.001 * e_crit * ones)
    assert st["damage"] > 0.0
    _, _, st = _fresh(mat, 2.0 * e_crit * ones)
    assert st["kappa"] == pytest.approx(2.0 * KAPPA_0, rel=1e-12)
    assert st["damage"] == pytest.approx(_damage_law(2.0 * KAPPA_0), rel=1e-12)

    if HYDRO_CONTRACTION_FACTOR[measure] > 0.0:
        e_crit_c = KAPPA_0 / HYDRO_CONTRACTION_FACTOR[measure]
        _, _, st = _fresh(mat, -0.999 * e_crit_c * ones)
        assert st["damage"] == 0.0
        _, _, st = _fresh(mat, -1.001 * e_crit_c * ones)
        assert st["damage"] > 0.0
    else:
        sigma, C_tan, st = _fresh(mat, -100.0 * KAPPA_0 * ones)
        assert st["damage"] == 0.0
        np.testing.assert_allclose(sigma, C_3D @ (-100.0 * KAPPA_0 * ones), rtol=1e-12)
        np.testing.assert_allclose(C_tan, C_3D, rtol=1e-12, atol=1e-12 * E)


def test_mazars_plane_strain_biaxial_contraction_never_damages():
    """Mazars en deformación plana con las dos principales del plano negativas.

    ``ε = [−1.0, −0.4, 0.5]·e``: centro ``−0.7e``, radio ``0.3905e`` ⇒
    principales ``−0.3095e, −1.0905e`` y ``ε_zz = 0`` ⇒ ``ε_eq = 0``: no hay
    daño aunque ``e = 100κ₀``. En esfuerzo plano, en cambio,
    ``ε_zz = −ν/(1−ν)·(−1.4e) = 0.35e > 0`` es la única principal positiva ⇒
    ``ε_eq = 0.35e``, que daña a partir de ``e = κ₀/0.35``.
    """
    direction = np.array([-1.0, -0.4, 0.5])
    mat_pe = _make("plane_strain", "mazars")
    _, _, st = _fresh(mat_pe, 100.0 * KAPPA_0 * direction)
    assert st["damage"] == 0.0

    mat_ps = _make("plane_stress", "mazars")
    e_crit = KAPPA_0 / (-C_ZZ * 1.4)
    _, _, st = _fresh(mat_ps, 0.999 * e_crit * direction)
    assert st["damage"] == 0.0
    _, _, st = _fresh(mat_ps, 2.0 * e_crit * direction)
    assert st["kappa"] == pytest.approx(2.0 * KAPPA_0, rel=1e-12)


# =============================================================================
# (d) Consistencia 2D ↔ 3D en carga, descarga y recarga.
# =============================================================================

# Direcciones de deformación en el plano (unidades arbitrarias; se escalan por
# medida para que la historia recorra ε_eq/κ₀ = 0 → 4 → 2 → 6).
DIRECTIONS_2D = {
    "uniaxial": (1.0, -NU, 0.0),          # uniaxial de esfuerzo en esfuerzo plano
    "biaxial": (1.0, 0.6, 0.0),
    "shear": (0.0, 0.0, 2.0),
    "general": (0.8, -0.3, 1.1),
    "compression": (-1.0, 0.3, 0.2),      # dominada por compresión
    "biaxial_contraction": (-1.0, -0.4, 0.5),
}


def _history_factors():
    """ε_eq/κ₀ a lo largo de la historia: carga 0→4, descarga 4→2, recarga 2→6."""
    return np.concatenate([np.linspace(0.0, 4.0, 41)[1:],
                           np.linspace(4.0, 2.0, 11)[1:],
                           np.linspace(2.0, 6.0, 21)[1:]])


@pytest.mark.parametrize("name", list(DIRECTIONS_2D))
@pytest.mark.parametrize("hypothesis", HYPOTHESES_2D)
@pytest.mark.parametrize("measure", MEASURES)
def test_2d_is_exact_restriction_of_3d(measure, hypothesis, name):
    """El 2D reproduce al 3D alimentado con su estado 3D, paso a paso.

    Estado 3D: ``[εxx, εyy, ε_zz, γxy, 0, 0]`` con
    ``ε_zz = −ν/(1−ν)(εxx+εyy)`` (esfuerzo plano) o ``0`` (deformación plana).
    - energy: con esa ``ε_zz``, ``σ_zz = 0`` ⇒ ``ε3ᵀC_3Dε3 = σ:ε =
      ε2ᵀC_psε2``; con ``ε_zz = 0``, ``ε3ᵀC_3Dε3 = ε2ᵀC_peε2``.
    - mazars y modified_von_mises se definen sobre el tensor 3D completo.
    Luego ``κ`` y ``d`` coinciden (1e-12 relativo; ``d`` además con 1e-14
    absoluto cerca del umbral) y ``σ`` en el plano también: con ``σ_zz = 0``
    las filas xx, yy, xy de ``C_3D·ε3`` son ``C_ps·ε2``. En esfuerzo plano el
    3D debe dar ``σ_zz = 0`` (redondeo de ``λ·tr ε + 2μ·ε_zz``).

    Mazars en deformación plana con ``biaxial_contraction`` tiene
    ``ε_eq ≡ 0`` (ver test anterior): la historia se escala entonces por
    la mayor componente y se exige que no haya daño en ninguno de los dos.
    """
    direction = np.asarray(DIRECTIONS_2D[name], dtype=float)
    eq_unit = _eps_eq_reference(direction, measure, hypothesis)
    # Umbral relativo: un autovalor nulo puede salir de eigvalsh como ±1e-17.
    expect_damage = eq_unit > 1e-12 * np.abs(direction).max()
    if expect_damage:
        unit = direction * (KAPPA_0 / eq_unit)
    else:
        unit = direction * (KAPPA_0 / np.abs(direction).max())

    m2 = _make(hypothesis, measure)
    m3 = _make("3d", measure)
    st2, st3 = m2.initial_state(), m3.initial_state()
    max_d = 0.0
    for s in _history_factors():
        eps2 = s * unit
        eps3 = _embed_3d(eps2, hypothesis)
        s2, _, st2 = _call(m2, eps2, st2)
        s3, _, st3 = _call(m3, eps3, st3)
        scale = E * np.abs(eps3).max()
        assert st2["kappa"] == pytest.approx(st3["kappa"], rel=1e-12)
        assert st2["damage"] == pytest.approx(st3["damage"], rel=1e-12, abs=1e-14)
        np.testing.assert_allclose(s2, [s3[0], s3[1], s3[3]], rtol=0, atol=1e-11 * scale)
        if hypothesis == "plane_stress":
            assert abs(s3[2]) <= 1e-11 * scale
        max_d = max(max_d, st2["damage"])
    if expect_damage:
        # En ε_eq = 6κ₀: d = 1 − e^{−0.5}/6 = 0.899.
        assert max_d == pytest.approx(_damage_law(6.0 * KAPPA_0), rel=1e-10)
    else:
        assert max_d == 0.0


# =============================================================================
# (e) Invariancia ante rotación de las direcciones principales.
# =============================================================================

def _rotation_zyz(a, b, c):
    def rz(t):
        return np.array([[math.cos(t), -math.sin(t), 0.0],
                         [math.sin(t), math.cos(t), 0.0],
                         [0.0, 0.0, 1.0]])

    def ry(t):
        return np.array([[math.cos(t), 0.0, math.sin(t)],
                         [0.0, 1.0, 0.0],
                         [-math.sin(t), 0.0, math.cos(t)]])
    return rz(a) @ ry(b) @ rz(c)


PRINCIPAL_3D = np.array([3.0, -1.0, 1.5]) * KAPPA_0
EULER_ANGLES = [(0.3, 1.1, -0.7), (2.0, 0.4, 1.3), (-1.2, 2.5, 0.2)]
PRINCIPAL_2D = np.array([2.5, -0.8]) * KAPPA_0
PLANE_ANGLES = [0.2, 0.9, 2.1, -1.3]


@pytest.mark.parametrize("measure", MEASURES)
def test_rotation_invariance_3d(measure):
    """Girar las direcciones principales con valores principales fijos no
    cambia ``κ`` ni ``d``.

    Las tres medidas son invariantes isótropos: energy depende de
    ``λ(tr ε)² + 2μ ε:ε``; mazars de las principales; modified_von_mises de
    ``I1`` y ``J2``. Valores principales ``(3.0, −1.0, 1.5)·κ₀``:
    energy ``ε_eq = 3.6893κ₀`` (``d = 0.792863``), mazars ``3.3541κ₀``
    (``d = 0.764394``), modified_von_mises ``5.4073κ₀`` (``d = 0.880983``):
    todas con daño no saturado. ``ε = R·diag(ε_i)·Rᵀ`` con
    ``R`` de Euler ZYZ, pasado a Voigt con ``γ = 2ε_ij``.

    Tolerancia 1e-12 relativa: la rotación y la descomposición espectral
    introducen errores de pocos ulp en las componentes (~1e-15 relativo a
    ``‖ε‖``), y ε_eq y d son lipschitzianos con constante O(1) aquí.
    """
    mat = _make("3d", measure)
    eps_principal = np.array([*PRINCIPAL_3D, 0.0, 0.0, 0.0])
    _, _, st0 = _fresh(mat, eps_principal)
    kappa_ref = _eps_eq_reference(eps_principal, measure, "3d")
    assert st0["kappa"] == pytest.approx(kappa_ref, rel=1e-12)
    assert 0.0 < st0["damage"] < DAMAGE_MAX
    for angles in EULER_ANGLES:
        R = _rotation_zyz(*angles)
        eps_rot = _voigt_from_tensor(R @ np.diag(PRINCIPAL_3D) @ R.T)
        assert np.abs(eps_rot[3:]).max() > 0.1 * KAPPA_0   # el giro es genuino
        _, _, st = _fresh(mat, eps_rot)
        assert st["kappa"] == pytest.approx(st0["kappa"], rel=1e-12)
        assert st["damage"] == pytest.approx(st0["damage"], rel=1e-12)


@pytest.mark.parametrize("hypothesis", HYPOTHESES_2D)
@pytest.mark.parametrize("measure", MEASURES)
def test_rotation_invariance_2d(measure, hypothesis):
    """Giro en el plano con principales ``(2.5, −0.8)·κ₀``.

    ``ε(θ) = [p1c² + p2s², p1s² + p2c², 2(p1−p2)sc]`` con ``c = cos θ``,
    ``s = sin θ``. ``εxx + εyy = p1 + p2`` no cambia, así que ``ε_zz`` de
    esfuerzo plano tampoco: el tensor 3D gira alrededor de ``z`` con sus
    principales fijas y ``κ``, ``d`` no deben cambiar (1e-12, como en 3D).
    """
    mat = _make(hypothesis, measure)
    p1, p2 = PRINCIPAL_2D
    eps0 = np.array([p1, p2, 0.0])
    _, _, st0 = _fresh(mat, eps0)
    assert st0["kappa"] == pytest.approx(_eps_eq_reference(eps0, measure, hypothesis), rel=1e-12)
    assert 0.0 < st0["damage"] < DAMAGE_MAX
    for theta in PLANE_ANGLES:
        c, s = math.cos(theta), math.sin(theta)
        eps = np.array([p1 * c * c + p2 * s * s,
                        p1 * s * s + p2 * c * c,
                        2.0 * (p1 - p2) * s * c])
        _, _, st = _fresh(mat, eps)
        assert st["kappa"] == pytest.approx(st0["kappa"], rel=1e-12)
        assert st["damage"] == pytest.approx(st0["damage"], rel=1e-12)


# =============================================================================
# (f) Tangente consistente frente a diferencias finitas centradas.
# (g) Simetría de la tangente de carga activa.
# =============================================================================

# Direcciones sin principales repetidas ni próximas a cero (ver
# resultados de la referencia en el docstring de _active_loading_point).
_FD_DIRECTIONS_2D = {
    "mixed": ([2.0, -0.5, 1.2], [1.0, 0.4, -0.7]),
    "two_positive": ([2.0, 0.8, 1.2], [0.6, 1.0, 0.5]),
    # εxx + εyy < 0 ⇒ en esfuerzo plano ε_zz > 0 es principal positiva: con
    # mazars activa la regla de la cadena de ε_zz en la tangente (en las otras
    # dos direcciones ε_zz < 0 y ⟨ε_zz⟩₊ = 0 la anula).
    "z_positive": ([-2.0, 0.5, 1.2], [-0.8, 0.3, 0.6]),
}
FD_DIRECTIONS = {
    "3d": {
        "mixed": ([2.0, -0.5, 0.9, 1.2, -0.4, 0.6], [1.0, 0.3, -0.6, 0.4, 0.8, -0.2]),
        "two_positive": ([2.0, 1.4, -0.8, 0.6, -0.3, 0.5], [0.9, 0.2, 0.5, -0.4, 0.3, 0.6]),
        "one_positive": ([-1.5, -0.8, 1.6, 0.7, -0.5, 0.4], [-0.6, 0.9, 0.4, 0.5, 0.2, -0.7]),
    },
    "plane_stress": _FD_DIRECTIONS_2D,
    "plane_strain": _FD_DIRECTIONS_2D,
}
FD_CASES = [(cfg, name) for cfg in CONFIGS for name in FD_DIRECTIONS[cfg]]
FD_STEP = 1.0e-9            # h = 1e-5·κ₀ (ver docstring del test)


def _active_loading_point(mat, measure, config, name):
    """Estado previo con ``κ = 1.5κ₀`` y punto de carga activa con ``ε_eq = 2κ₀``.

    Ambas direcciones se escalan con la referencia (homogeneidad de grado 1):
    ``prev = 1.5·N(dir_b)``, ``eps = 2·N(dir_a)``. El margen
    ``ε_eq(eps) − κ_prev = 0.5κ₀ = 5e-5`` es 5e4 veces el paso ``h``: las
    perturbaciones ``eps ± h·e_j`` siguen en carga activa.

    Principales de ``dir_a`` (sin repetidas ni nulas, salvo la ``ε_zz = 0``
    exacta de deformación plana, que no se perturba):
    3D mixed ≈ (2.1854, 0.8951, −0.6804); 3D two_positive ≈
    (2.1348, 1.3014, −0.8362); 3D one_positive ≈ (1.6344, −0.6652, −1.6692);
    2D mixed: (2.1365, −0.6365) y ``ε_zz`` = −0.375 (esfuerzo plano) o 0;
    2D two_positive: (2.2485, 0.5515) y ``ε_zz`` = −0.7 o 0;
    2D z_positive: (0.6365, −2.1365) y ``ε_zz`` = +0.375 o 0.
    """
    dir_a, dir_b = (np.asarray(v, dtype=float) for v in FD_DIRECTIONS[config][name])
    prev = 1.5 * _normalized(dir_b, measure, config)
    eps = 2.0 * _normalized(dir_a, measure, config)
    _, _, st_prev = _fresh(mat, prev)
    assert st_prev["kappa"] == pytest.approx(1.5 * KAPPA_0, rel=1e-12)
    return eps, st_prev


@pytest.mark.parametrize("config,name", FD_CASES)
@pytest.mark.parametrize("measure", MEASURES)
def test_consistent_tangent_vs_central_differences(measure, config, name):
    """``C_tan`` de carga activa = derivada exacta de ``σ(ε)`` (desde el mismo
    estado previo).

    ``C_fd[:, j] = [σ(ε + h e_j) − σ(ε − h e_j)]/(2h)``, ``h = 1e-9``.
    Error de las diferencias centradas, relativo a ``max|C| ~ E``:
    truncamiento ``(h/ε_s)²/6·O(1)`` con la escala de variación
    ``ε_s ~ κ = 2e-4`` ⇒ ~1e-11; redondeo ``ε_mach·|σ|/h ~
    2e-16·(3e-4 E)/1e-9`` ⇒ ~1e-10 E (hasta ~1e-9 E si el redondeo de
    ε_eq, p. ej. vía autovalores, es de algunos ulp); con la referencia
    propia del test el error medido es de 2e-11 a 7e-10·max|C| (el mayor en
    modified_von_mises con z_positive, donde ``ε_eq = 2κ₀`` exige
    ``|ε|`` ~ 1e-3 y el redondeo ``ε_mach·|σ|/h`` crece en proporción). El término de
    ablandamiento ``(∂d/∂κ)(C_eε)⊗(∂ε_eq/∂ε)`` es O(E)
    (``(1−d)(1+ακ)`` ≈ 0.54), y cualquier término olvidado (p. ej. la regla de
    la cadena de ``ε_zz`` en esfuerzo plano, factor ``ν/(1−ν) = 0.25``) es
    ≥ 1e-2·max|C|. La tolerancia ``1e-6·max|C|`` deja tres órdenes sobre el
    error de las diferencias y cuatro bajo el menor error físico.
    """
    mat = _make(config, measure)
    eps, st_prev = _active_loading_point(mat, measure, config, name)
    _, C_tan, st = _call(mat, eps, st_prev)
    assert st["kappa"] == pytest.approx(2.0 * KAPPA_0, rel=1e-12)
    assert st["damage"] == pytest.approx(_damage_law(2.0 * KAPPA_0), rel=1e-12)

    n = eps.size
    C_fd = np.zeros((n, n))
    for j in range(n):
        step = np.zeros(n)
        step[j] = FD_STEP
        sp, _, stp = _call(mat, eps + step, st_prev)
        sm, _, stm = _call(mat, eps - step, st_prev)
        assert stp["kappa"] > st_prev["kappa"] and stm["kappa"] > st_prev["kappa"]
        C_fd[:, j] = (sp - sm) / (2.0 * FD_STEP)
    scale = np.abs(C_fd).max()
    np.testing.assert_allclose(C_tan, C_fd, rtol=0, atol=1e-6 * scale)
    # La tangente de carga activa no es la secante: el término de
    # ablandamiento es O(E).
    assert np.abs(C_tan - (1.0 - st["damage"]) * C_ELASTIC[config]).max() > 1e-2 * scale


@pytest.mark.parametrize("config,name", FD_CASES)
@pytest.mark.parametrize("measure", MEASURES)
def test_active_tangent_symmetry(measure, config, name):
    """Simetría de la tangente de carga activa.

    - energy: ``∂ε_eq/∂ε = C_eε/(E·ε_eq)`` ⇒ ``C_alg = (1−d)C_e −
      (∂d/∂κ)/(E ε_eq)·(C_eε)⊗(C_eε)``, simétrica. Tolerancia 1e-12·max|C|:
      ``(a·u_i)(b·u_j)`` y ``(a·u_j)(b·u_i)`` difieren en pocos ulp.
    - mazars y modified_von_mises: ``∂ε_eq/∂ε`` no es paralelo a ``C_eε`` en
      general, así que ``(C_eε)⊗(∂ε_eq/∂ε)`` no es simétrico. En estos
      estados la referencia propia da ``max|C−Cᵀ|/max|C|`` entre 0.05 y
      0.93 (mazars) y entre 0.53 y 1.22 (modified_von_mises); se exige
      > 1e-3, cincuenta veces por debajo del menor.
    """
    mat = _make(config, measure)
    eps, st_prev = _active_loading_point(mat, measure, config, name)
    _, C_tan, _ = _call(mat, eps, st_prev)
    scale = np.abs(C_tan).max()
    asym = np.abs(C_tan - C_tan.T).max()
    if measure == "energy":
        assert asym <= 1e-12 * scale
    else:
        assert asym > 1e-3 * scale


# =============================================================================
# Estados generales frente a la referencia; descarga secante.
# =============================================================================

@pytest.mark.parametrize("config,name", FD_CASES)
@pytest.mark.parametrize("measure", MEASURES)
def test_general_state_matches_reference(measure, config, name):
    """En un estado general desde cero: ``κ = ε_eq`` de referencia,
    ``d = d(κ)``, ``σ = (1−d)·C_e·ε``.

    Los estados son las direcciones ``dir_a`` de la tangente escaladas con la
    referencia a ``ε_eq = 2.3κ₀`` (``d = 0.618219``): el material calcula su
    propia ``ε_eq`` de una combinación general de componentes y debe dar
    ``κ = 2.3κ₀``. Comprueba la definición de cada medida fuera de los
    estados especiales (en 2D, energy se referencia con la matriz 2D y
    mazars / modified_von_mises con el tensor 3D con ``ε_zz``). Tolerancias
    1e-12 en ``κ`` y ``d`` y 1e-10 en ``σ`` (redondeo de pocos ulp).
    """
    mat = _make(config, measure)
    dir_a = np.asarray(FD_DIRECTIONS[config][name][0], dtype=float)
    eps = 2.3 * _normalized(dir_a, measure, config)
    sigma, _, st = _fresh(mat, eps)
    kappa_ref = 2.3 * KAPPA_0
    d_ref = _damage_law(kappa_ref)
    assert d_ref == pytest.approx(0.618219, abs=1e-6)
    assert st["kappa"] == pytest.approx(kappa_ref, rel=1e-12)
    assert st["damage"] == pytest.approx(d_ref, rel=1e-12)
    C_e = C_ELASTIC[config]
    np.testing.assert_allclose(sigma, (1.0 - d_ref) * (C_e @ eps),
                               rtol=1e-10, atol=1e-12 * E * np.abs(eps).max())


@pytest.mark.parametrize("config", CONFIGS)
@pytest.mark.parametrize("measure", MEASURES)
def test_elastic_and_unloading_secant(measure, config):
    """Rama elástica y descarga: la tangente es la secante ``(1−d)·C_e``.

    1. Desde cero con ``ε_eq = 0.5κ₀``: ``d = 0``, ``σ = C_e·ε``, ``C = C_e``.
    2. Carga a ``ε_eq = 3κ₀`` (``d = 1 − e^{−0.2}/3 = 0.727090``) y descarga a
       la mitad: ``κ`` y ``d`` no cambian, ``σ = (1−d)·C_e·ε`` y
       ``C = (1−d)·C_e`` (simétrica con cualquier medida).
    """
    mat = _make(config, measure)
    C_e = C_ELASTIC[config]
    dir_a = np.asarray(FD_DIRECTIONS[config]["mixed"][0], dtype=float)
    unit = _normalized(dir_a, measure, config)

    eps_el = 0.5 * unit
    sigma, C_tan, st = _fresh(mat, eps_el)
    assert st["damage"] == 0.0
    np.testing.assert_allclose(sigma, C_e @ eps_el, rtol=1e-12, atol=1e-12 * E * KAPPA_0)
    np.testing.assert_allclose(C_tan, C_e, rtol=1e-12, atol=1e-12 * E)

    _, _, st_loaded = _fresh(mat, 3.0 * unit)
    d = _damage_law(3.0 * KAPPA_0)
    assert st_loaded["damage"] == pytest.approx(d, rel=1e-12)
    eps_un = 1.5 * unit
    sigma, C_tan, st = _call(mat, eps_un, st_loaded)
    assert st["kappa"] == pytest.approx(st_loaded["kappa"], rel=1e-14)
    assert st["damage"] == pytest.approx(st_loaded["damage"], rel=1e-14)
    np.testing.assert_allclose(sigma, (1.0 - d) * (C_e @ eps_un),
                               rtol=1e-10, atol=1e-12 * E * np.abs(eps_un).max())
    np.testing.assert_allclose(C_tan, (1.0 - d) * C_e, rtol=1e-10, atol=1e-12 * E)


# =============================================================================
# Contrato de parámetros.
# =============================================================================

@pytest.mark.parametrize("config", CONFIGS)
def test_default_measure_is_energy(config):
    """Sin ``equivalent_strain`` la medida es energy: en cortante puro el
    umbral es ``γ = κ₀·√(2(1+ν)) = 1.549193e-4`` (mazars daría 2e-4 y la
    norma del tensor retirada, ``ε_eq = |γ|/√2``, 1.414214e-4)."""
    kwargs = dict(E=E, nu=NU, kappa_0=KAPPA_0, alpha=ALPHA)
    if config == "3d":
        mat = IsotropicDamage3D(**kwargs)
        n, slot = 6, 3
    else:
        mat = IsotropicDamage2D(hypothesis=config, **kwargs)
        n, slot = 3, 2
    gamma_crit = KAPPA_0 / SHEAR_FACTOR["energy"]
    eps = np.zeros(n)
    eps[slot] = 0.999 * gamma_crit
    _, _, st = _fresh(mat, eps)
    assert st["damage"] == 0.0
    eps[slot] = 2.0 * gamma_crit
    _, _, st = _fresh(mat, eps)
    assert st["kappa"] == pytest.approx(2.0 * KAPPA_0, rel=1e-12)


def _build(config, **extra):
    kwargs = dict(E=E, nu=NU, kappa_0=KAPPA_0, alpha=ALPHA, **extra)
    if config == "3d":
        return IsotropicDamage3D(**kwargs)
    return IsotropicDamage2D(hypothesis=config, **kwargs)


@pytest.mark.parametrize("config", CONFIGS)
def test_parameter_contract(config):
    """``k`` es obligatorio y positivo con modified_von_mises y no se admite
    con energy ni mazars; un nombre de medida desconocido se rechaza (la
    norma del tensor ya no existe como opción)."""
    with pytest.raises(ValueError):
        _build(config, equivalent_strain="modified_von_mises")
    for bad_k in (0.0, -2.0):
        with pytest.raises(ValueError):
            _build(config, equivalent_strain="modified_von_mises", k=bad_k)
    for measure in ("energy", "mazars"):
        with pytest.raises(ValueError):
            _build(config, equivalent_strain=measure, k=K_RATIO)
    with pytest.raises(ValueError):
        _build(config, equivalent_strain="tensor_norm")
    # Las combinaciones válidas se construyen sin error.
    _build(config, equivalent_strain="energy")
    _build(config, equivalent_strain="mazars")
    _build(config, equivalent_strain="modified_von_mises", k=K_RATIO)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
