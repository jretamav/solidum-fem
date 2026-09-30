"""Deformación equivalente configurable del daño isótropo 2D y 3D.

Cambio de formulación aprobado por el usuario el 2026-09-29:
``IsotropicDamage2D`` e ``IsotropicDamage3D`` retiran la norma del tensor
``√(ε:ε)`` (no estándar y mal calibrada: ``ε_eq = ε·√(1+2ν²)`` en tracción
uniaxial) y admiten tres medidas de la literatura con el parámetro
``equivalent_strain``:

- ``energy`` (Simó y Ju 1987, por defecto): ``ε_eq = √(εᵀ·C_e·ε/E)``;
- ``mazars`` (Mazars 1986): ``ε_eq = √(Σ⟨ε_i⟩₊²)``;
- ``modified_von_mises`` (de Vree et al. 1995), con ``k = f_c/f_t``.

Lo que se comprueba aquí:

1. el valor de ``ε_eq`` frente a una implementación de referencia escrita
   aparte sobre el tensor 3D (autovalores con ``numpy``, ``σ:ε`` con Lamé,
   ``J₂ = ½·e:e``), en 2D (tensión y deformación plana) y 3D;
2. la calibración ``ε_eq = ε`` en tracción uniaxial de esfuerzo y los
   valores de compresión uniaxial de cada medida;
3. la tangente consistente frente a diferencias finitas en carga activa,
   incluidos valores principales repetidos (Mazars);
4. el kernel por lotes idéntico bit a bit a ``compute_state``, y el
   ensamblaje por familias (serie y paralelo) igual al camino por elemento;
5. ``IS_SYMMETRIC`` por instancia y la simetría exacta con ``energy``;
6. la validación del constructor y la lectura desde YAML.
"""
import math
import os
import tempfile

import numpy as np
import pytest

import solidum
from solidum.materials.damage_2d import IsotropicDamage2D
from solidum.materials.damage_3d import IsotropicDamage3D
from solidum.math.batch import StateSchema
from solidum.registry import MaterialRegistry

E, NU, KAPPA_0, ALPHA, K = 30.0e9, 0.2, 1.0e-4, 1.0e3, 10.0

MEDIDAS = {"energy": {}, "mazars": {}, "modified_von_mises": {"k": K}}
CONFIGS = ["plane_stress", "plane_strain", "3d"]


def _material(config, medida, nu=NU):
    kw = dict(MEDIDAS[medida])
    if config == "3d":
        return IsotropicDamage3D(E, nu, KAPPA_0, ALPHA, equivalent_strain=medida, **kw)
    return IsotropicDamage2D(E, nu, KAPPA_0, ALPHA, hypothesis=config,
                             equivalent_strain=medida, **kw)


def _estado_3d(eps, config, nu=NU):
    """Estado 3D (Voigt 6D) asociado a ``eps`` según la configuración."""
    eps = np.asarray(eps, dtype=float)
    if config == "3d":
        return eps.copy()
    exx, eyy, gxy = eps
    ezz = -nu / (1.0 - nu) * (exx + eyy) if config == "plane_stress" else 0.0
    return np.array([exx, eyy, ezz, gxy, 0.0, 0.0])


def _tensor(e6):
    return np.array([[e6[0], 0.5 * e6[3], 0.5 * e6[5]],
                     [0.5 * e6[3], e6[1], 0.5 * e6[4]],
                     [0.5 * e6[5], 0.5 * e6[4], e6[2]]])


def _eps_eq_referencia(e6, medida, nu=NU, k=K):
    """ε_eq escrita sobre el tensor 3D, por un camino distinto del código."""
    t = _tensor(e6)
    tr = np.trace(t)
    if medida == "energy":
        lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
        mu = E / (2.0 * (1.0 + nu))
        sig = lam * tr * np.eye(3) + 2.0 * mu * t
        return math.sqrt(float(np.tensordot(sig, t)) / E)
    if medida == "mazars":
        w = np.linalg.eigvalsh(t)
        return math.sqrt(float(np.sum(np.maximum(w, 0.0) ** 2)))
    dev = t - tr / 3.0 * np.eye(3)
    J2 = 0.5 * float(np.tensordot(dev, dev))
    b = (k - 1.0) / (1.0 - 2.0 * nu)
    return (b * tr + math.sqrt(b * b * tr * tr + 12.0 * k / (1.0 + nu) ** 2 * J2)) / (2.0 * k)


def _rotacion(a, b, c):
    ca, sa, cb, sb, cc, sc = (math.cos(a), math.sin(a), math.cos(b), math.sin(b),
                              math.cos(c), math.sin(c))
    Rz = np.array([[ca, -sa, 0.0], [sa, ca, 0.0], [0.0, 0.0, 1.0]])
    Ry = np.array([[cb, 0.0, sb], [0.0, 1.0, 0.0], [-sb, 0.0, cb]])
    Rx = np.array([[1.0, 0.0, 0.0], [0.0, cc, -sc], [0.0, sc, cc]])
    return Rz @ Ry @ Rx


def _voigt(t):
    return np.array([t[0, 0], t[1, 1], t[2, 2], 2.0 * t[0, 1], 2.0 * t[1, 2], 2.0 * t[0, 2]])


def _kappa_de(material, eps):
    """``κ`` tras un paso desde ``κ = 0``: vale ``ε_eq(eps)``."""
    _, _, st = material.compute_state(np.asarray(eps, float), {"kappa": 0.0, "damage": 0.0})
    return st["kappa"]


# =============================================================================
# 1. Valor de ε_eq frente a la referencia independiente
# =============================================================================

@pytest.mark.parametrize("medida", list(MEDIDAS))
@pytest.mark.parametrize("config", CONFIGS)
def test_eps_eq_coincide_con_la_referencia(config, medida):
    m = _material(config, medida)
    rng = np.random.default_rng(20260929)
    dim = 6 if config == "3d" else 3
    for _ in range(25):
        eps = 3.0e-4 * rng.uniform(-1.0, 1.0, dim)
        esperado = _eps_eq_referencia(_estado_3d(eps, config), medida)
        assert _kappa_de(m, eps) == pytest.approx(esperado, rel=1e-12, abs=1e-20)


@pytest.mark.parametrize("medida", list(MEDIDAS))
@pytest.mark.parametrize("config", ["plane_stress", "3d"])
def test_traccion_uniaxial_de_esfuerzo_calibra_kappa0(config, medida):
    """Tracción uniaxial de esfuerzo (``ε_yy = ε_zz = −ν·ε``): ``ε_eq = ε`` con
    las tres medidas, así que el daño empieza en ``σ = E·κ₀ = f_t``. En
    compresión: energética ``|ε|``, Mazars ``√2·ν·|ε|``, von Mises
    modificado ``|ε|/k`` (sólo dañan las extensiones laterales en Mazars)."""
    m = _material(config, medida)
    compresion = {"energy": 1.0, "mazars": math.sqrt(2.0) * NU,
                  "modified_von_mises": 1.0 / K}[medida]
    for e in (0.5e-4, 2.0e-4):
        for signo, factor in ((1.0, 1.0), (-1.0, compresion)):
            x = signo * e
            eps = [x, -NU * x, 0.0] if config == "plane_stress" else [x, -NU * x, -NU * x, 0, 0, 0]
            assert _kappa_de(m, eps) == pytest.approx(factor * e, rel=1e-13)


@pytest.mark.parametrize("medida", list(MEDIDAS))
def test_umbral_del_dano_en_traccion_uniaxial(medida):
    """``d = 0`` justo por debajo de ``ε = κ₀`` y ``d > 0`` justo por encima."""
    m = _material("plane_stress", medida)
    for f, con_dano in ((0.999, False), (1.001, True)):
        x = f * KAPPA_0
        _, _, st = m.compute_state(np.array([x, -NU * x, 0.0]), m.initial_state())
        assert (st["damage"] > 0.0) is con_dano


def test_compresion_hidrostatica():
    """Sin extensiones no hay daño con Mazars ni con von Mises modificado
    (``k > 1``); la medida energética sí daña (es simétrica)."""
    eps = -3.0e-4 * np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
    assert _kappa_de(_material("3d", "mazars"), eps) == 0.0
    assert _kappa_de(_material("3d", "modified_von_mises"), eps) == pytest.approx(0.0, abs=1e-19)
    # ε = −e·I: σ_eff = −(3λ + 2μ)·e·I con 3λ + 2μ = E/(1 − 2ν), así que
    # εᵀ·C_e·ε/E = 3·e²/(1 − 2ν).
    esperado = 3.0e-4 * math.sqrt(3.0 / (1.0 - 2.0 * NU))
    assert _kappa_de(_material("3d", "energy"), eps) == pytest.approx(esperado, rel=1e-13)


@pytest.mark.parametrize("medida", list(MEDIDAS))
def test_invariancia_frente_a_rotaciones(medida):
    """Las tres medidas son funciones isótropas del tensor: ``ε_eq(R·ε·Rᵀ) =
    ε_eq(ε)`` en 3D, y en 2D bajo giros en el plano."""
    m3 = _material("3d", medida)
    t = _tensor([2.0e-4, -0.5e-4, 0.3e-4, 1.2e-4, 0.4e-4, -0.7e-4])
    ref = _kappa_de(m3, _voigt(t))
    for ang in ((0.3, -1.1, 2.0), (1.7, 0.4, -0.6)):
        R = _rotacion(*ang)
        assert _kappa_de(m3, _voigt(R @ t @ R.T)) == pytest.approx(ref, rel=1e-12)
    for config in ("plane_stress", "plane_strain"):
        m2 = _material(config, medida)
        t2 = np.array([[2.0e-4, 0.6e-4], [0.6e-4, -0.5e-4]])
        ref2 = _kappa_de(m2, [t2[0, 0], t2[1, 1], 2.0 * t2[0, 1]])
        for th in (0.4, 1.3, 2.9):
            c, s = math.cos(th), math.sin(th)
            R = np.array([[c, -s], [s, c]])
            r = R @ t2 @ R.T
            assert _kappa_de(m2, [r[0, 0], r[1, 1], 2.0 * r[0, 1]]) == pytest.approx(ref2, rel=1e-12)


# =============================================================================
# 2. Tangente consistente frente a diferencias finitas
# =============================================================================

_R = _rotacion(0.7, -0.4, 1.9)
ESTADOS = {
    "3d": {
        "general": 1.0e-4 * np.array([2.0, -0.5, 0.3, 1.2, 0.4, -0.7]),
        "todo_extension": 1.0e-4 * np.array([3.0, 2.0, 1.5, 0.4, 0.2, -0.3]),
        # Dos deformaciones principales iguales (3, 3, −1)·1e-4 giradas.
        "principales_repetidas": _voigt(_R @ np.diag([3.0e-4, 3.0e-4, -1.0e-4]) @ _R.T),
    },
    "2d": {
        "general": 1.0e-4 * np.array([2.0, -0.5, 1.2]),
        "todo_extension": 1.0e-4 * np.array([3.0, 2.0, 0.4]),
        # Principales del plano iguales (r = 0).
        "principales_repetidas": 1.0e-4 * np.array([2.5, 2.5, 0.0]),
    },
}


def _tangente_fd(m, eps, estado, h=1.0e-10):
    n = eps.size
    C = np.zeros((n, n))
    for j in range(n):
        d = np.zeros(n)
        d[j] = h
        sp, _, _ = m.compute_state(eps + d, estado)
        sm, _, _ = m.compute_state(eps - d, estado)
        C[:, j] = (sp - sm) / (2.0 * h)
    return C


@pytest.mark.parametrize("nombre", ["general", "todo_extension", "principales_repetidas"])
@pytest.mark.parametrize("medida", list(MEDIDAS))
@pytest.mark.parametrize("config", CONFIGS)
def test_tangente_consistente_frente_a_diferencias_finitas(config, medida, nombre):
    """Carga activa desde ``κ = κ₀``: la tangente cerrada es la derivada de
    ``σ(ε)``. Diferencia centrada con ``h = 1e-10`` sobre ``|ε| ~ 1e-4``:
    truncamiento ``O((h/ε)²) ~ 1e-12`` y redondeo ``~ u·|σ|/h ~ 1e-10``
    relativos, así que ``1e-6`` deja tres órdenes de margen."""
    m = _material(config, medida)
    eps = ESTADOS["3d" if config == "3d" else "2d"][nombre].copy()
    estado = {"kappa": KAPPA_0, "damage": 0.0}
    _, C, st = m.compute_state(eps, estado)
    assert 0.0 < st["damage"] < 0.99, "el estado debe estar en carga activa sin saturar"
    C_fd = _tangente_fd(m, eps, estado)
    np.testing.assert_allclose(C, C_fd, rtol=0.0, atol=1e-6 * np.abs(C_fd).max())


# =============================================================================
# 3. Simetría de la tangente e IS_SYMMETRIC por instancia
# =============================================================================

@pytest.mark.parametrize("medida", list(MEDIDAS))
@pytest.mark.parametrize("config", CONFIGS)
def test_is_symmetric_por_instancia(config, medida):
    m = _material(config, medida)
    assert type(m).IS_SYMMETRIC is False          # cota conservadora de clase
    assert m.IS_SYMMETRIC is (medida == "energy")
    eps = ESTADOS["3d" if config == "3d" else "2d"]["general"]
    _, C, st = m.compute_state(eps, {"kappa": KAPPA_0, "damage": 0.0})
    assert st["damage"] > 0.0
    if medida == "energy":
        # C_alg = (1−d)·C_e − (∂d/∂κ)/(E·ε_eq)·(C_e·ε)⊗(C_e·ε): simétrica bit a bit.
        np.testing.assert_array_equal(C, C.T)
    else:
        assert np.linalg.norm(C - C.T) / np.linalg.norm(C) > 1e-3


# =============================================================================
# 4. Camino por lotes ≡ camino por elemento
# =============================================================================

def _run_kernel(material, strain, S_old):
    schema = StateSchema(material.STATE_SCHEMA)
    S_new = np.zeros(schema.n_state)
    sigma = np.zeros(material.STRAIN_DIM)
    flag = np.zeros(1, dtype=np.int8)
    C = np.ascontiguousarray(material.batch_matrix(), dtype=np.float64)
    params = np.ascontiguousarray(material.batch_params(), dtype=np.float64)
    Ct = material.batch_kernel()(np.ascontiguousarray(strain, dtype=np.float64),
                                 S_old, S_new, params, C, sigma, np.zeros_like(C), flag)
    return sigma, np.asarray(Ct), S_new


@pytest.mark.parametrize("medida", list(MEDIDAS))
@pytest.mark.parametrize("config", CONFIGS)
def test_camino_por_lotes_identico_al_camino_por_elemento(config, medida):
    """Historia elástica, carga, descarga, recarga y compresión: esfuerzo,
    tangente y fila de estado idénticos bit a bit (misma función compilada)."""
    m = _material(config, medida)
    base = ESTADOS["3d" if config == "3d" else "2d"]["general"]
    schema = StateSchema(m.STATE_SCHEMA)
    S = np.zeros(schema.n_state)
    schema.pack(m.initial_state(), S)
    estado = m.initial_state()
    for f in (0.01, 1.0, 2.0, 1.0, 3.0, -2.0):
        eps = f * base
        s_ref, C_ref, estado = m.compute_state(eps, estado)
        s, C, S_new = _run_kernel(m, eps, S)
        np.testing.assert_array_equal(s, s_ref)
        np.testing.assert_array_equal(C, C_ref)
        fila = np.zeros(schema.n_state)
        schema.pack(estado, fila)
        np.testing.assert_array_equal(S_new, fila)
        S = S_new
    assert estado["damage"] > 0.0


def _malla(config, material):
    """Malla pequeña: 3×2 ``Quad4`` en 2D o 2×1×1 ``Hex8`` en 3D, lado 1."""
    from solidum.core.domain import Domain
    from solidum.elements.solid_2d import Quad4
    from solidum.elements.solid_3d import Hex8
    dom = Domain()
    if config == "3d":
        nx, ny, nz = 2, 1, 1
        idx = {}
        for k in range(nz + 1):
            for j in range(ny + 1):
                for i in range(nx + 1):
                    idx[i, j, k] = len(idx) + 1
                    dom.add_node(idx[i, j, k], [float(i), float(j), float(k)])
        for i in range(nx):
            c = [idx[i, 0, 0], idx[i + 1, 0, 0], idx[i + 1, 1, 0], idx[i, 1, 0],
                 idx[i, 0, 1], idx[i + 1, 0, 1], idx[i + 1, 1, 1], idx[i, 1, 1]]
            dom.add_element(Hex8(i + 1, [dom.get_node(n) for n in c], material))
    else:
        nx, ny = 3, 2
        idx = {}
        for j in range(ny + 1):
            for i in range(nx + 1):
                idx[i, j] = len(idx) + 1
                dom.add_node(idx[i, j], [float(i), float(j)])
        e = 1
        for j in range(ny):
            for i in range(nx):
                c = [idx[i, j], idx[i + 1, j], idx[i + 1, j + 1], idx[i, j + 1]]
                dom.add_element(Quad4(e, [dom.get_node(n) for n in c], material, thickness=0.1))
                e += 1
    dom.generate_equation_numbers(verbose=False)
    return dom


@pytest.mark.parametrize("parallel", [False, True])
@pytest.mark.parametrize("medida", list(MEDIDAS))
@pytest.mark.parametrize("config", CONFIGS)
def test_ensamblaje_por_lotes_coincide_con_el_camino_por_elemento(config, medida, parallel):
    """El kernel de familia (serie y paralelo) llama al núcleo del material
    como función de primera clase; con cada medida —también Mazars 3D, que
    llama a LAPACK (``eigh``) dentro del bucle— ``K`` y ``F_int`` coinciden
    con el camino por elemento salvo el orden de las sumas del ensamblaje
    (misma tolerancia que ``test_batch_assembly``)."""
    from solidum.math.assembly import Assembler
    dom = _malla(config, _material(config, medida))
    rng = np.random.default_rng(7)
    U = 6.0e-4 * rng.uniform(-1.0, 1.0, dom.total_dofs)
    K0, F0 = Assembler(dom, batch=False).assemble_non_linear_system(U)
    danos = [st["damage"] for el in dom.elements.values() for st in el.state.vars_trial]
    assert max(danos) > 0.1 and min(danos) < max(danos), "estados dañados y distintos"
    K1, F1 = Assembler(dom, batch=True, parallel=parallel).assemble_non_linear_system(U)
    A0, A1 = K0.toarray(), K1.toarray()
    escala = np.abs(A0).max()
    np.testing.assert_allclose(A1, A0, rtol=0.0, atol=1e-14 * escala)
    np.testing.assert_allclose(F1, F0, rtol=0.0, atol=1e-14 * np.abs(F0).max())


# =============================================================================
# 5. Validación del constructor
# =============================================================================

CLASES = {"IsotropicDamage2D": IsotropicDamage2D, "IsotropicDamage3D": IsotropicDamage3D}


@pytest.mark.parametrize("nombre", list(CLASES))
class TestValidacionDelConstructor:

    def _crear(self, nombre, **kw):
        return CLASES[nombre](E=E, nu=NU, kappa_0=KAPPA_0, alpha=ALPHA, **kw)

    def test_por_defecto_es_la_energetica(self, nombre):
        m = self._crear(nombre)
        assert m.equivalent_strain == "energy"
        assert m.k is None
        assert m.IS_SYMMETRIC is True

    @pytest.mark.parametrize("medida", ["norm", "Energy", "von_mises", "", None, 1])
    def test_medida_desconocida(self, nombre, medida):
        with pytest.raises(ValueError) as exc:
            self._crear(nombre, equivalent_strain=medida)
        assert nombre in str(exc.value)
        assert "equivalent_strain" in str(exc.value)

    def test_von_mises_modificado_exige_k(self, nombre):
        with pytest.raises(ValueError) as exc:
            self._crear(nombre, equivalent_strain="modified_von_mises")
        assert nombre in str(exc.value)
        assert "k" in str(exc.value)

    @pytest.mark.parametrize("medida", ["energy", "mazars"])
    def test_k_no_se_admite_con_otras_medidas(self, nombre, medida):
        with pytest.raises(ValueError) as exc:
            self._crear(nombre, equivalent_strain=medida, k=10.0)
        assert nombre in str(exc.value)

    @pytest.mark.parametrize("k", [0.0, -3.0, float("nan"), float("inf"), "diez", True])
    def test_k_debe_ser_positivo_y_finito(self, nombre, k):
        with pytest.raises(ValueError) as exc:
            self._crear(nombre, equivalent_strain="modified_von_mises", k=k)
        assert nombre in str(exc.value)

    def test_k_entero_se_admite(self, nombre):
        m = self._crear(nombre, equivalent_strain="modified_von_mises", k=12)
        assert m.k == 12.0 and isinstance(m.k, float)
        assert m.IS_SYMMETRIC is False

    def test_el_registro_construye_con_los_kwargs(self, nombre):
        m = MaterialRegistry.create(nombre, E=E, nu=NU, kappa_0=KAPPA_0, alpha=ALPHA,
                                    equivalent_strain="mazars")
        assert m.equivalent_strain == "mazars" and m.k is None


# =============================================================================
# 6. Lectura desde YAML
# =============================================================================

DELTA = 2.5 * KAPPA_0

_YAML_QUAD4 = """
nodes:
  - {{id: 1, coords: [0.0, 0.0]}}
  - {{id: 2, coords: [1.0, 0.0]}}
  - {{id: 3, coords: [1.0, 1.0]}}
  - {{id: 4, coords: [0.0, 1.0]}}
materials:
  - {{id: 1, type: IsotropicDamage2D, E: {E}, nu: {nu}, kappa_0: {k0}, alpha: {alpha},
     hypothesis: plane_stress{extra}}}
elements:
  - {{id: 1, type: Quad4, material: 1, thickness: 1.0, nodes: [1, 2, 3, 4]}}
boundary_conditions:
  - {{node_id: 1, ux: 0.0, uy: 0.0}}
  - {{node_id: 2, ux: {delta}, uy: 0.0}}
  - {{node_id: 3, ux: {delta}}}
  - {{node_id: 4, ux: 0.0}}
solver:
  type: NonlinearSolver
  num_steps: 5
"""

_YAML_HEX8 = """
nodes:
  - {{id: 1, coords: [0.0, 0.0, 0.0]}}
  - {{id: 2, coords: [1.0, 0.0, 0.0]}}
  - {{id: 3, coords: [1.0, 1.0, 0.0]}}
  - {{id: 4, coords: [0.0, 1.0, 0.0]}}
  - {{id: 5, coords: [0.0, 0.0, 1.0]}}
  - {{id: 6, coords: [1.0, 0.0, 1.0]}}
  - {{id: 7, coords: [1.0, 1.0, 1.0]}}
  - {{id: 8, coords: [0.0, 1.0, 1.0]}}
materials:
  - {{id: 1, type: IsotropicDamage3D, E: {E}, nu: {nu}, kappa_0: {k0}, alpha: {alpha}{extra}}}
elements:
  - {{id: 1, type: Hex8, material: 1, nodes: [1, 2, 3, 4, 5, 6, 7, 8]}}
boundary_conditions:
  - {{node_id: 1, ux: 0.0, uy: 0.0, uz: 0.0}}
  - {{node_id: 2, ux: {delta}, uy: 0.0, uz: 0.0}}
  - {{node_id: 3, ux: {delta}, uz: 0.0}}
  - {{node_id: 4, ux: 0.0, uz: 0.0}}
  - {{node_id: 5, ux: 0.0, uy: 0.0}}
  - {{node_id: 6, ux: {delta}, uy: 0.0}}
  - {{node_id: 7, ux: {delta}}}
  - {{node_id: 8, ux: 0.0}}
solver:
  type: NonlinearSolver
  num_steps: 5
"""

CASOS_YAML = {
    "quad4_mazars": (_YAML_QUAD4, ", equivalent_strain: mazars", "mazars", None, (2, 3)),
    "quad4_energia_por_defecto": (_YAML_QUAD4, "", "energy", None, (2, 3)),
    "hex8_von_mises_modificado": (_YAML_HEX8, ", equivalent_strain: modified_von_mises, k: 10.0",
                                  "modified_von_mises", 10.0, (2, 3, 6, 7)),
}


def _escribir(texto):
    fd, path = tempfile.mkstemp(suffix=".yaml", dir=os.path.dirname(os.path.abspath(__file__)))
    os.close(fd)
    with open(path, "w", encoding="utf-8") as f:
        f.write(texto)
    return path


@pytest.mark.parametrize("caso", list(CASOS_YAML))
def test_yaml_tension_uniaxial(caso):
    """Tracción uniaxial de esfuerzo por desplazamiento prescrito ``δ`` en un
    elemento: con las tres medidas ``ε_eq = δ``, así que ``d`` y la reacción
    ``(1−d)·E·δ`` son analíticos."""
    plantilla, extra, medida, k, nodos_derecha = CASOS_YAML[caso]
    path = _escribir(plantilla.format(E=E, nu=NU, k0=KAPPA_0, alpha=ALPHA,
                                      delta=DELTA, extra=extra))
    try:
        from solidum.utils.yaml_parser import YamlParser
        parser = YamlParser(path)
        parser.parse()
        mat = parser.materials[1]
        assert mat.equivalent_strain == medida
        assert mat.k == k
        assert mat.IS_SYMMETRIC is (medida == "energy")
        res = solidum.run_yaml(path)
    finally:
        os.remove(path)
    d = 1.0 - (KAPPA_0 / DELTA) * math.exp(-ALPHA * (DELTA - KAPPA_0))
    reaccion = sum(res.reactions_by_node[n]["ux"] for n in nodos_derecha)
    assert reaccion == pytest.approx((1.0 - d) * E * DELTA, rel=1e-8)


def test_yaml_von_mises_modificado_sin_k_se_rechaza():
    path = _escribir(_YAML_QUAD4.format(E=E, nu=NU, k0=KAPPA_0, alpha=ALPHA, delta=DELTA,
                                        extra=", equivalent_strain: modified_von_mises"))
    try:
        from solidum.utils.yaml_parser import YamlParser
        with pytest.raises(ValueError, match="exige el parámetro"):
            YamlParser(path).parse()
    finally:
        os.remove(path)
