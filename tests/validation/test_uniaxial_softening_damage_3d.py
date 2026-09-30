"""Uniaxial confinado softening — IsotropicDamage3D vs curva σ-ε analítica.

Referencias
-----------
- Lemaitre J., Chaboche J.-L. (1990). *Mechanics of Solid Materials*,
  Cambridge UP, Cap. 7 (continuum damage mechanics).
- Simó J.C., Ju J.W. (1987). Strain- and stress-based continuum damage
  models — I. Formulation. *Int. J. Solids Struct.* 23, 821-840.

Concepto
--------
1 elemento ``Hex8`` con confinamiento transversal total (``u_y = u_z = 0``
en todos los nodos, ``u_x = 0`` en la cara ``x=0``) y desplazamiento ``u_x``
prescrito en la cara ``x=1``. Bajo esta restricción cinemática el campo es
uniforme con ``ε = (ε_xx, 0, 0, 0, 0, 0)`` en todos los Gauss points. La
deformación equivalente de cada medida (``equivalent_strain``, 2026-09-29)
en este estado, con ``ε_xx > 0``:

- ``energy``: ``εᵀ·C_e·ε/E = C_e^{xx,xx}·ε_xx²/E``, es decir
  ``ε_eq = ε_xx·√((1−ν)/[(1+ν)(1−2ν)])`` (``≈ 1,054·ε_xx`` con ``ν = 0,2``);
- ``mazars``: única deformación principal positiva ``ε_xx`` ⇒ ``ε_eq = ε_xx``;
- ``modified_von_mises``: ``I₁ = ε_xx``, ``J₂ = ε_xx²/3`` ⇒
  ``ε_eq = [b·ε_xx + √(b²·ε_xx² + 4k·ε_xx²/(1+ν)²)]/(2k)`` con
  ``b = (k−1)/(1−2ν)``.

Hasta el 2026-09-29 la medida era la norma del tensor, ``ε_eq = |ε_xx|`` en
este estado.

Bajo carga monotónica desde cero, la ley de daño exponencial centralizada
en ``solidum.materials._softening`` produce, para ``ε_eq > κ_0``:

.. math::

   d(\\varepsilon_{eq}) = 1 - \\frac{\\kappa_0}{\\varepsilon_{eq}}\\,
                        \\exp\\!\\left(-\\alpha\\,(\\varepsilon_{eq} - \\kappa_0)\\right)

(saturado a ``DAMAGE_MAX``). El esfuerzo nominal es ``σ = (1-d)·C_e·ε``, así
que en confinamiento triaxial uniaxial el componente activo es

.. math::

   \\sigma_{xx}(\\varepsilon_{xx}) = (1 - d) \\cdot C_e^{xx,xx} \\cdot \\varepsilon_{xx}

con ``C_e^{xx,xx} = E(1-\\nu)/[(1+\\nu)(1-2\\nu)]`` (módulo confinado axial,
no ``E``). Por simetría ``σ_yy = σ_zz = (1-d)·C_e^{yy,xx}·ε_xx`` con
``C_e^{yy,xx} = E·\\nu/[(1+\\nu)(1-2\\nu)]``.

Diseño del test
---------------
Múltiples solves independientes barriendo ``ε_xx`` desde por debajo del
umbral ``κ_0`` hasta varias veces ``κ_0`` (régimen plenamente dañado).
Cada solve parte de cero (sin estado heredado), lo cual es válido porque
para carga monotónica el daño es *path-independent*: el estado final
depende sólo del ``ε_eq`` máximo histórico, no de cómo se llegó.

Validaciones (con cada una de las tres medidas)
-----------------------------------------------
1. **Régimen elástico exacto** (``ε_eq < κ_0``): ``σ_xx = C_e^{xx,xx}·ε_xx``
   a precisión máquina, ``d = 0`` exacto.
2. **Curva σ-ε analítica** (``ε_xx ∈ [1,2·κ_0, 5·κ_0]``, 20 puntos): el
   ``σ_xx`` FEM coincide con la fórmula analítica a precisión < 1e-10 (el
   integrador 8-Gauss promedia perfectamente bajo campo uniforme).
3. **Saturación a DAMAGE_MAX** (``ε_xx >> κ_0``): ``d ≈ DAMAGE_MAX``,
   ``σ_xx ≈ (1-DAMAGE_MAX)·C_e^{xx,xx}·ε_xx``.
"""
import math
import os
import sys

import numpy as np
import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from solidum.constants import DAMAGE_MAX
from solidum.core.domain import Domain
from solidum.elements.solid_3d import Hex8
from solidum.materials.damage_3d import IsotropicDamage3D
from solidum.math.assembly import Assembler
from solidum.math.convergence import ConvergenceCriterion
from solidum.math.solvers import NonlinearSolver


# Parámetros físicos del benchmark.
E_YOUNG = 2.0e4
NU = 0.2
KAPPA_0 = 1.0e-4
ALPHA = 500.0
K_FC_FT = 10.0     # k = f_c/f_t de von Mises modificado

# Medidas de deformación equivalente y sus parámetros.
MEDIDAS = {"energy": {}, "mazars": {}, "modified_von_mises": {"k": K_FC_FT}}

# Módulos confinados de Elastic3D (filas 0 y 1 de C_e):
#   C_xx_xx = E(1-ν) / [(1+ν)(1-2ν)]   (módulo axial confinado)
#   C_yy_xx = E·ν     / [(1+ν)(1-2ν)]   (acoplamiento Poisson en confinamiento)
_COEF = E_YOUNG / ((1.0 + NU) * (1.0 - 2.0 * NU))
C_AXIAL_CONFINED = _COEF * (1.0 - NU)
C_LATERAL_CONFINED = _COEF * NU

# Cubo unitario Hex8.
HEX8_UNIT_COORDS = [
    (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
]
HEX8_FACE_XMIN = (0, 3, 4, 7)
HEX8_FACE_XMAX = (1, 2, 5, 6)


def _eps_eq_analytical(eps_xx: float, medida: str) -> float:
    """``ε_eq`` de la deformación uniaxial confinada ``ε = (ε_xx, 0, …, 0)``
    con ``ε_xx > 0`` (ver docstring del módulo)."""
    if medida == "energy":
        return eps_xx * math.sqrt((1.0 - NU) / ((1.0 + NU) * (1.0 - 2.0 * NU)))
    if medida == "mazars":
        return eps_xx
    b = (K_FC_FT - 1.0) / (1.0 - 2.0 * NU)
    return (b * eps_xx + math.sqrt(b * b * eps_xx ** 2
                                   + 4.0 * K_FC_FT * eps_xx ** 2 / (1.0 + NU) ** 2)) / (2.0 * K_FC_FT)


def _damage_analytical(eps_xx: float, medida: str) -> float:
    """``d(ε_eq)`` con la ley exponencial y el ``ε_eq`` de la medida."""
    eps_eq = _eps_eq_analytical(eps_xx, medida)
    if eps_eq <= KAPPA_0:
        return 0.0
    d = 1.0 - (KAPPA_0 / eps_eq) * math.exp(-ALPHA * (eps_eq - KAPPA_0))
    return min(d, DAMAGE_MAX)


def _sigma_xx_analytical(eps_xx: float, medida: str) -> float:
    """σ_xx(ε_xx) analítico bajo confinamiento triaxial uniaxial."""
    d = _damage_analytical(eps_xx, medida)
    return (1.0 - d) * C_AXIAL_CONFINED * eps_xx


def _solve_confined_uniaxial(eps_xx_target: float, medida: str, num_steps: int = 4):
    """Cubo confinado con u_x prescrito en cara +x.

    Confinamiento total transversal: u_y = u_z = 0 en todos los nodos.
    u_x = 0 en cara x=0; u_x = eps_xx_target en cara x=1. Cubo de longitud 1
    ⇒ ε_xx = eps_xx_target uniforme.
    """
    mat = IsotropicDamage3D(E=E_YOUNG, nu=NU, kappa_0=KAPPA_0, alpha=ALPHA,
                            equivalent_strain=medida, **MEDIDAS[medida])
    dom = Domain()
    nodes = [dom.add_node(i + 1, list(c)) for i, c in enumerate(HEX8_UNIT_COORDS)]
    elem = Hex8(1, nodes, mat)
    dom.add_element(elem)

    # Confinamiento transversal total: u_y = u_z = 0 en cada nodo
    for n in nodes:
        n.fix_dof('uy', 0.0)
        n.fix_dof('uz', 0.0)
    # u_x = 0 en x=0, u_x = eps_xx_target en x=1
    for i in HEX8_FACE_XMIN:
        nodes[i].fix_dof('ux', 0.0)
    for i in HEX8_FACE_XMAX:
        nodes[i].fix_dof('ux', eps_xx_target)

    dom.generate_equation_numbers(verbose=False)
    assembler = Assembler(dom)

    F_ext = np.zeros(dom.total_dofs)
    conv = ConvergenceCriterion(rtol_force=1.0e-9, rtol_disp=1.0e-9)
    solver = NonlinearSolver(assembler, convergence=conv, num_steps=num_steps)
    U = solver.solve(F_ext)
    return elem, U


# =============================================================================
# Test 1 — Régimen elástico (ε_eq < κ_0).
# =============================================================================

@pytest.mark.parametrize("medida", list(MEDIDAS))
def test_elastic_regime_exact(medida):
    """ε_eq < κ_0: σ_xx = C_e^{xx,xx}·ε_xx exacto, d = 0."""
    eps_xx = 0.5 * KAPPA_0
    assert _eps_eq_analytical(eps_xx, medida) < KAPPA_0     # sanidad del caso
    elem, U = _solve_confined_uniaxial(eps_xx, medida, num_steps=2)

    gs = elem.compute_gauss_state(U)
    # σ_xx promedio en los 8 Gauss points (campo uniforme ⇒ todos iguales)
    sigma_xx_fem = gs['stress'][:, 0].mean()
    sigma_yy_fem = gs['stress'][:, 1].mean()
    sigma_zz_fem = gs['stress'][:, 2].mean()

    sigma_xx_expected = C_AXIAL_CONFINED * eps_xx
    sigma_yy_expected = C_LATERAL_CONFINED * eps_xx

    np.testing.assert_allclose(
        sigma_xx_fem, sigma_xx_expected, rtol=1.0e-12,
        err_msg=f"σ_xx FEM={sigma_xx_fem:.6e} vs C·ε={sigma_xx_expected:.6e}"
    )
    np.testing.assert_allclose(
        sigma_yy_fem, sigma_yy_expected, rtol=1.0e-12,
        err_msg=f"σ_yy FEM={sigma_yy_fem:.6e} vs C·ε={sigma_yy_expected:.6e}"
    )
    np.testing.assert_allclose(
        sigma_zz_fem, sigma_yy_expected, rtol=1.0e-12,
        err_msg=f"σ_zz FEM={sigma_zz_fem:.6e} vs C·ε={sigma_yy_expected:.6e}"
    )
    # d = 0 en todos los Gauss points
    for k, state in enumerate(elem.state.vars):
        assert state['damage'] == 0.0, (
            f"GP {k}: d = {state['damage']:.3e} debe ser 0 con ε_xx={eps_xx:.3e}"
        )


# =============================================================================
# Test 2 — Curva σ-ε analítica (ε_eq > κ_0, barrido).
# =============================================================================

@pytest.mark.parametrize("medida", list(MEDIDAS))
def test_sigma_epsilon_curve_matches_analytical(medida):
    """Barrido de ε_xx en [1.2·κ_0, 5·κ_0]: σ_xx FEM = σ_xx analítico a precisión < 1e-10.

    El integrador 8-Gauss reproduce exactamente el campo uniforme bajo
    confinamiento total, así que la única fuente de error es el cómputo
    del kernel del material — que debe coincidir con la fórmula analítica
    a precisión máquina.
    """
    eps_values = np.linspace(1.2 * KAPPA_0, 5.0 * KAPPA_0, 20)
    for eps_xx in eps_values:
        elem, U = _solve_confined_uniaxial(float(eps_xx), medida)
        gs = elem.compute_gauss_state(U)

        sigma_xx_fem = gs['stress'][:, 0].mean()
        sigma_xx_expected = _sigma_xx_analytical(float(eps_xx), medida)
        d_expected = _damage_analytical(float(eps_xx), medida)
        assert d_expected > 0.0

        # σ_xx coincide a precisión muy alta
        rel_err = abs(sigma_xx_fem - sigma_xx_expected) / abs(sigma_xx_expected)
        assert rel_err < 1.0e-10, (
            f"ε_xx={eps_xx:.4e}: σ_xx FEM={sigma_xx_fem:.6e} vs "
            f"analítico={sigma_xx_expected:.6e} (err_rel={rel_err:.3e})"
        )

        # Variabilidad inter-Gauss point < 1e-10 (campo uniforme)
        sigma_std = gs['stress'][:, 0].std()
        sigma_avg_abs = abs(sigma_xx_fem)
        assert sigma_std / max(sigma_avg_abs, 1.0) < 1.0e-10, (
            f"ε_xx={eps_xx:.4e}: variabilidad inter-Gauss σ_xx = {sigma_std:.3e} "
            f"demasiado alta (campo debe ser uniforme)"
        )

        # d en cada Gauss point coincide con d analítico
        for k, state in enumerate(elem.state.vars):
            d_fem = state['damage']
            rel_err_d = abs(d_fem - d_expected) / max(d_expected, 1.0e-15)
            assert rel_err_d < 1.0e-12, (
                f"ε_xx={eps_xx:.4e}, GP {k}: d FEM={d_fem:.6e} vs "
                f"analítico={d_expected:.6e} (err_rel={rel_err_d:.3e})"
            )


# =============================================================================
# Test 3 — Saturación a DAMAGE_MAX.
# =============================================================================

@pytest.mark.parametrize("medida", list(MEDIDAS))
def test_saturation_at_damage_max(medida):
    """ε_xx muy grande: d alcanza DAMAGE_MAX, σ_xx = (1-DAMAGE_MAX)·C·ε_xx."""
    # ε_xx tal que el d analítico (sin cap) sería > DAMAGE_MAX por amplio margen.
    # Con α=500, κ_0=1e-4 y ε_eq ≥ ε_xx = 200·κ_0: d_uncapped ≥ 1 − 2.4e-7.
    eps_xx = 200.0 * KAPPA_0
    eps_eq = _eps_eq_analytical(eps_xx, medida)
    d_uncapped = 1.0 - (KAPPA_0 / eps_eq) * math.exp(-ALPHA * (eps_eq - KAPPA_0))
    assert d_uncapped > DAMAGE_MAX, (
        f"Sanity check del setup: d sin cap = {d_uncapped:.6e}, "
        f"debe ser > DAMAGE_MAX = {DAMAGE_MAX}"
    )

    elem, U = _solve_confined_uniaxial(eps_xx, medida, num_steps=8)
    gs = elem.compute_gauss_state(U)

    # d capeado en DAMAGE_MAX
    for k, state in enumerate(elem.state.vars):
        assert state['damage'] == DAMAGE_MAX, (
            f"GP {k}: d = {state['damage']:.6e} no coincide con DAMAGE_MAX = {DAMAGE_MAX}"
        )

    # σ_xx = (1-DAMAGE_MAX)·C·ε_xx
    sigma_xx_fem = gs['stress'][:, 0].mean()
    sigma_xx_expected = (1.0 - DAMAGE_MAX) * C_AXIAL_CONFINED * eps_xx
    np.testing.assert_allclose(
        sigma_xx_fem, sigma_xx_expected, rtol=1.0e-12,
        err_msg=f"σ_xx FEM={sigma_xx_fem:.6e} vs (1-DAMAGE_MAX)·C·ε={sigma_xx_expected:.6e}"
    )


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
