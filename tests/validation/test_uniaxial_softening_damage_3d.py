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
uniforme con ``ε = (ε_xx, 0, 0, 0, 0, 0)`` en todos los Gauss points, y la
norma de deformación equivalente del modelo es ``ε_eq = √(ε_xx²) = |ε_xx|``.

Bajo carga monotónica desde cero, la ley de daño exponencial centralizada
en ``solidum.materials._softening`` produce, para ``ε_xx > κ_0``:

.. math::

   d(\\varepsilon_{xx}) = 1 - \\frac{\\kappa_0}{\\varepsilon_{xx}}\\,
                        \\exp\\!\\left(-\\alpha\\,(\\varepsilon_{xx} - \\kappa_0)\\right)

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

Validaciones
------------
1. **Régimen elástico exacto** (``ε_xx < κ_0``): ``σ_xx = C_e^{xx,xx}·ε_xx``
   a precisión máquina, ``d = 0`` exacto.
2. **Curva σ-ε analítica** (``ε_xx > κ_0``, 20 puntos): el ``σ_xx`` FEM
   coincide con la fórmula analítica a precisión < 1e-10 (el integrador
   8-Gauss promedia perfectamente bajo campo uniforme).
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


def _damage_analytical(eps_xx: float) -> float:
    """``d(ε_xx)`` con ley exponencial; ``ε_eq = |ε_xx|`` bajo confinamiento."""
    eps_eq = abs(eps_xx)
    if eps_eq <= KAPPA_0:
        return 0.0
    d = 1.0 - (KAPPA_0 / eps_eq) * math.exp(-ALPHA * (eps_eq - KAPPA_0))
    return min(d, DAMAGE_MAX)


def _sigma_xx_analytical(eps_xx: float) -> float:
    """σ_xx(ε_xx) analítico bajo confinamiento triaxial uniaxial."""
    d = _damage_analytical(eps_xx)
    return (1.0 - d) * C_AXIAL_CONFINED * eps_xx


def _sigma_yy_analytical(eps_xx: float) -> float:
    """σ_yy = σ_zz analítico (acoplamiento Poisson confinado)."""
    d = _damage_analytical(eps_xx)
    return (1.0 - d) * C_LATERAL_CONFINED * eps_xx


def _solve_confined_uniaxial(eps_xx_target: float, num_steps: int = 4):
    """Cubo confinado con u_x prescrito en cara +x.

    Confinamiento total transversal: u_y = u_z = 0 en todos los nodos.
    u_x = 0 en cara x=0; u_x = eps_xx_target en cara x=1. Cubo de longitud 1
    ⇒ ε_xx = eps_xx_target uniforme.
    """
    mat = IsotropicDamage3D(E=E_YOUNG, nu=NU, kappa_0=KAPPA_0, alpha=ALPHA)
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
# Test 1 — Régimen elástico (ε_xx < κ_0).
# =============================================================================

def test_elastic_regime_exact():
    """ε_xx < κ_0: σ_xx = C_e^{xx,xx}·ε_xx exacto, d = 0."""
    eps_xx = 0.5 * KAPPA_0
    elem, U = _solve_confined_uniaxial(eps_xx, num_steps=2)

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
            f"GP {k}: d = {state['damage']:.3e} debe ser 0 con ε_xx={eps_xx:.3e} < κ_0={KAPPA_0:.3e}"
        )


# =============================================================================
# Test 2 — Curva σ-ε analítica (ε_xx > κ_0, barrido).
# =============================================================================

def test_sigma_epsilon_curve_matches_analytical():
    """Barrido de ε_xx en [1.2·κ_0, 5·κ_0]: σ_xx FEM = σ_xx analítico a precisión < 1e-10.

    El integrador 8-Gauss reproduce exactamente el campo uniforme bajo
    confinamiento total, así que la única fuente de error es el cómputo
    del kernel del material — que debe coincidir con la fórmula analítica
    a precisión máquina.
    """
    eps_values = np.linspace(1.2 * KAPPA_0, 5.0 * KAPPA_0, 20)
    for eps_xx in eps_values:
        elem, U = _solve_confined_uniaxial(float(eps_xx))
        gs = elem.compute_gauss_state(U)

        sigma_xx_fem = gs['stress'][:, 0].mean()
        sigma_xx_expected = _sigma_xx_analytical(float(eps_xx))
        d_expected = _damage_analytical(float(eps_xx))

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

def test_saturation_at_damage_max():
    """ε_xx muy grande: d alcanza DAMAGE_MAX, σ_xx = (1-DAMAGE_MAX)·C·ε_xx."""
    # ε_xx tal que el d analítico (sin cap) sería > DAMAGE_MAX por amplio margen.
    # Con α=500, κ_0=1e-4: ε_xx = 200·κ_0 da d_uncapped ≈ 1 - 2.4e-7 (saturado lejos).
    eps_xx = 200.0 * KAPPA_0
    d_uncapped = 1.0 - (KAPPA_0 / eps_xx) * math.exp(-ALPHA * (eps_xx - KAPPA_0))
    assert d_uncapped > DAMAGE_MAX, (
        f"Sanity check del setup: d sin cap = {d_uncapped:.6e}, "
        f"debe ser > DAMAGE_MAX = {DAMAGE_MAX}"
    )

    elem, U = _solve_confined_uniaxial(eps_xx, num_steps=8)
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
