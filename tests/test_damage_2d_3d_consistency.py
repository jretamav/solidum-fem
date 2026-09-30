"""``IsotropicDamage2D`` es la restricción exacta de ``IsotropicDamage3D``.

El daño 2D y el 3D miden la deformación equivalente con la norma del tensor
de deformaciones, ``ε_eq = √(ε:ε)``. En tensión plana la deformación fuera del
plano no es nula: con daño escalar isótropo, ``σ_zz = (1−d)·(λ·tr ε + 2μ·ε_zz)
= 0`` da ``ε_zz = −ν/(1−ν)·(ε_xx + ε_yy)`` para cualquier ``d``. Hasta el
2026-09-29 el 2D la omitía y, con el mismo estado de esfuerzo, retrasaba el
inicio del daño un factor ``√(1+2ν²)/√(1+ν²)`` respecto del 3D (hallazgo
clas-3 de la auditoría; corrección aprobada por el usuario).

Estos tests comparan, punto a punto de una historia de carga, el 2D con el 3D
alimentado con el estado de deformación 3D que corresponde: en tensión plana
con ``ε_zz`` de arriba y en deformación plana con ``ε_zz = 0``. Deben coincidir
el daño, la variable histórica y las componentes en el plano del esfuerzo, y
en tensión plana el 3D debe dar ``σ_zz = 0``.
"""
import numpy as np
import pytest

from solidum.materials.damage_2d import IsotropicDamage2D
from solidum.materials.damage_3d import IsotropicDamage3D

E, NU, KAPPA_0, ALPHA = 30.0e9, 0.2, 1.0e-4, 1.0e3


def _estado_3d(eps_2d, hipotesis):
    exx, eyy, gxy = eps_2d
    ezz = -NU / (1.0 - NU) * (exx + eyy) if hipotesis == "plane_stress" else 0.0
    return np.array([exx, eyy, ezz, gxy, 0.0, 0.0])


def _historia(direccion):
    """Carga monótona hasta bien entrado el ablandamiento, descarga parcial y
    recarga: recorre la rama elástica, el inicio del daño, el ablandamiento y
    la descarga secante."""
    s = np.concatenate([np.linspace(0.0, 4.0, 41)[1:],
                        np.linspace(4.0, 2.0, 11)[1:],
                        np.linspace(2.0, 6.0, 21)[1:]])
    return [f * KAPPA_0 * np.asarray(direccion) for f in s]


# Direcciones de deformación en el plano: uniaxial de esfuerzo en tensión
# plana (ε_yy = −ν ε_xx), biaxial, cortante puro y un estado general.
DIRECCIONES = {
    "uniaxial": (1.0, -NU, 0.0),
    "biaxial": (1.0, 0.6, 0.0),
    "cortante": (0.0, 0.0, 2.0),
    "general": (0.8, -0.3, 1.1),
}


@pytest.mark.parametrize("hipotesis", ["plane_stress", "plane_strain"])
@pytest.mark.parametrize("nombre", list(DIRECCIONES))
def test_el_2d_reproduce_al_3d_en_toda_la_historia(hipotesis, nombre):
    m2 = IsotropicDamage2D(E, NU, KAPPA_0, ALPHA, hypothesis=hipotesis)
    m3 = IsotropicDamage3D(E, NU, KAPPA_0, ALPHA)
    st2, st3 = m2.initial_state(), m3.initial_state()
    hubo_dano = False
    for eps in _historia(DIRECCIONES[nombre]):
        s2, _, st2 = m2.compute_state(eps, st2)
        s3, _, st3 = m3.compute_state(_estado_3d(eps, hipotesis), st3)
        assert st2["kappa"] == pytest.approx(st3["kappa"], rel=1e-12)
        assert st2["damage"] == pytest.approx(st3["damage"], rel=1e-12, abs=1e-14)
        escala = E * KAPPA_0
        np.testing.assert_allclose(s2, [s3[0], s3[1], s3[3]], rtol=0, atol=1e-10 * escala)
        if hipotesis == "plane_stress":
            assert abs(s3[2]) <= 1e-10 * escala
        hubo_dano |= st2["damage"] > 0.1
    assert hubo_dano, "la historia debe entrar en el ablandamiento"


def test_inicio_del_dano_en_traccion_uniaxial_de_esfuerzo():
    """Con σ = σ_xx·e_x, el daño empieza cuando la norma completa del tensor
    de deformaciones, ε_xx·√(1 + 2ν²), alcanza κ₀ (tensión plana)."""
    m2 = IsotropicDamage2D(E, NU, KAPPA_0, ALPHA, hypothesis="plane_stress")
    umbral = KAPPA_0 / np.sqrt(1.0 + 2.0 * NU ** 2)
    _, _, st = m2.compute_state(np.array([0.999 * umbral, -NU * 0.999 * umbral, 0.0]),
                                m2.initial_state())
    assert st["damage"] == 0.0
    _, _, st = m2.compute_state(np.array([1.01 * umbral, -NU * 1.01 * umbral, 0.0]),
                                m2.initial_state())
    assert st["damage"] > 0.0


@pytest.mark.parametrize("hipotesis", ["plane_stress", "plane_strain"])
def test_tangente_consistente_frente_a_diferencias_finitas(hipotesis):
    """La derivada de ε_eq incluye ahora el término de ε_zz: la tangente de
    carga activa debe seguir siendo la derivada exacta de σ(ε)."""
    m = IsotropicDamage2D(E, NU, KAPPA_0, ALPHA, hypothesis=hipotesis)
    estado = m.initial_state()
    _, _, estado = m.compute_state(np.array([1.5e-4, -0.4e-4, 0.9e-4]), estado)
    eps = np.array([2.0e-4, -0.5e-4, 1.2e-4])        # carga activa desde ese estado
    _, C, _ = m.compute_state(eps, estado)
    h = 1e-10
    C_fd = np.zeros((3, 3))
    for j in range(3):
        d = np.zeros(3)
        d[j] = h
        sp, _, _ = m.compute_state(eps + d, estado)
        sm, _, _ = m.compute_state(eps - d, estado)
        C_fd[:, j] = (sp - sm) / (2.0 * h)
    np.testing.assert_allclose(C, C_fd, rtol=1e-5, atol=1e-6 * np.abs(C_fd).max())
