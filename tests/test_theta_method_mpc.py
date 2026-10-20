"""Restricciones multipunto (MPC) en el problema térmico (auditoría 2026-09-30).

El ``ThetaMethodSolver`` y el estacionario térmico eliminan las MPC con el
modelo afín del ADR 0004, ``T = T_op·T_libre + g``. Dos familias de
comprobaciones:

- **Oráculo del nodo fundido**: atar con ``T_s = T_m`` dos nodos coincidentes
  equivale a fundirlos en uno, así que el modelo atado y el de malla conforme
  dan la misma historia de temperaturas hasta el redondeo. Se comprueba en 2D
  (``Quad4Thermal``) y 3D (``Hex8Thermal``), en el transitorio con
  ``θ ∈ {1, 1/2, 0}`` y capacidad concentrada o consistente, y en el
  estacionario (``LinearSolver``).
- **Dirichlet frente a MPC en el transitorio**. Hasta el 2026-09-30 el solver
  trataba como Dirichlet todo DOF eliminado del sistema reducido, esclavos de
  MPC incluidos: ``dirichlet_func`` sustituía el ``g`` completo y borraba el
  de los esclavos (la MPC se incumplía sin aviso), el aviso de ``T_initial``
  citaba los esclavos como Dirichlet contradichos, la guarda «sin Dirichlet»
  dejaba pasar un modelo con sólo MPC (la guarda se retiró del transitorio el
  2026-10-20), y el instante inicial almacenado usaba
  el valor declarado en vez de ``dirichlet_func(0)``. Ahora el solver usa el
  Dirichlet declarado (``ConstraintSet.dirichlet_dofs``) y reconstruye las
  filas esclavas con ``ConstraintSet.build_g``.
"""
from __future__ import annotations

import contextlib
import logging

import numpy as np
import pytest

from solidum.bc.constraints import ConstraintSet
from solidum.core.domain import Domain
from solidum.elements.thermal.hex8_thermal import Hex8Thermal
from solidum.elements.thermal.quad4_thermal import Quad4Thermal
from solidum.materials.thermal_conduction import ThermalConduction
from solidum.math.assembly import Assembler
from solidum.math.solvers.linear import LinearSolver
from solidum.math.solvers.theta_method import ThetaMethodSolver

# Acero estructural: coeficientes no unitarios a propósito.
K_ACERO, C_ACERO, RHO_ACERO = 45.0, 460.0, 7850.0
DUP = 100_000        # desplazamiento de los ids de los nodos duplicados
T_IZQ = 100.0        # temperatura impuesta en el borde izquierdo
Q_DER = 50.0         # calor nodal aportado en el borde derecho [W]


# ----------------------------------------------------------------------
# Captura de avisos (el logger ``solidum`` no propaga, ADR 0005)
# ----------------------------------------------------------------------

class _CapturaLog(logging.Handler):
    def __init__(self):
        super().__init__(level=logging.WARNING)
        self.mensajes: list[str] = []

    def emit(self, record):
        self.mensajes.append(record.getMessage())

    def con(self, texto: str) -> list[str]:
        return [m for m in self.mensajes if texto in m]


@contextlib.contextmanager
def capturar_avisos():
    handler = _CapturaLog()
    logger = logging.getLogger("solidum")
    logger.addHandler(handler)
    try:
        yield handler
    finally:
        logger.removeHandler(handler)


# ----------------------------------------------------------------------
# Modelos: fundido (malla conforme) y atado por MPC en una columna interior
# ----------------------------------------------------------------------

def placa_2d(atado: bool, *, nx=6, ny=3, L=0.3, H=0.15, ic=3,
             T_izq: float | None = T_IZQ, dirichlet_interfaz: bool = False,
             mpc_g: float = 0.0):
    """Placa de ``Quad4Thermal``. Con ``atado`` la columna ``ic`` de nodos se
    duplica para los elementos a su derecha y se impone ``T_dup = T_orig + g``.

    ``dirichlet_interfaz`` fija además ``T = 50`` en el nodo inferior de la
    columna ``ic`` (maestro de una MPC en el modelo atado): su esclavo queda
    determinado por el cierre transitivo sin ser un Dirichlet declarado.
    """
    mat = ThermalConduction(k=K_ACERO, c=C_ACERO, density=RHO_ACERO, dim=2)
    dom = Domain()
    nid = {}
    for j in range(ny + 1):
        for i in range(nx + 1):
            nid[i, j] = j * (nx + 1) + i + 1
            dom.add_node(nid[i, j], [L * i / nx, H * j / ny])
    dup = {}
    if atado:
        for j in range(ny + 1):
            dup[ic, j] = DUP + nid[ic, j]
            dom.add_node(dup[ic, j], [L * ic / nx, H * j / ny])
    e = 1
    for j in range(ny):
        for i in range(nx):
            esquinas = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)]
            ids = [dup[c] if (atado and c[0] == ic and i >= ic) else nid[c]
                   for c in esquinas]
            dom.add_element(Quad4Thermal(
                e, [dom.nodes[k] for k in ids], mat, thickness=0.02))
            e += 1
    if T_izq is not None:
        for j in range(ny + 1):
            dom.nodes[nid[0, j]].fix_dof("T", T_izq)
    if dirichlet_interfaz:
        dom.nodes[nid[ic, 0]].fix_dof("T", 50.0)
    if atado:
        for key, n in dup.items():
            dom.add_linear_constraint((n, "T"), [(nid[key], "T")], [1.0], g=mpc_g)
    dom.generate_equation_numbers(verbose=False)
    return dom, nid, dup


def placa_3d(atado: bool, *, nx=4, ny=1, nz=1, L=0.2, H=0.05, W=0.05, ic=2,
             dirichlet_interfaz: bool = False):
    """Prisma de ``Hex8Thermal`` con la misma construcción que :func:`placa_2d`."""
    mat = ThermalConduction(k=K_ACERO, c=C_ACERO, density=RHO_ACERO, dim=3)
    dom = Domain()
    nid = {}
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                nid[i, j, k] = (k * (ny + 1) + j) * (nx + 1) + i + 1
                dom.add_node(nid[i, j, k], [L * i / nx, H * j / ny, W * k / nz])
    dup = {}
    if atado:
        for k in range(nz + 1):
            for j in range(ny + 1):
                dup[ic, j, k] = DUP + nid[ic, j, k]
                dom.add_node(dup[ic, j, k], [L * ic / nx, H * j / ny, W * k / nz])
    e = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                esquinas = [(i, j, k), (i + 1, j, k), (i + 1, j + 1, k), (i, j + 1, k),
                            (i, j, k + 1), (i + 1, j, k + 1), (i + 1, j + 1, k + 1),
                            (i, j + 1, k + 1)]
                ids = [dup[c] if (atado and c[0] == ic and i >= ic) else nid[c]
                       for c in esquinas]
                dom.add_element(Hex8Thermal(e, [dom.nodes[q] for q in ids], mat))
                e += 1
    for k in range(nz + 1):
        for j in range(ny + 1):
            dom.nodes[nid[0, j, k]].fix_dof("T", T_IZQ)
    if dirichlet_interfaz:
        dom.nodes[nid[ic, 0, 0]].fix_dof("T", 50.0)
    if atado:
        for key, n in dup.items():
            dom.add_linear_constraint((n, "T"), [(nid[key], "T")], [1.0])
    dom.generate_equation_numbers(verbose=False)
    return dom, nid, dup


def calor_borde_derecho(dom, nid, q=Q_DER) -> np.ndarray:
    F = np.zeros(dom.total_dofs)
    imax = max(key[0] for key in nid)
    for key, n in nid.items():
        if key[0] == imax:
            F[dom.nodes[n].dofs["T"]] += q
    return F


def valores_declarados(dom) -> np.ndarray:
    """Vector ``(ndof,)`` con el valor Dirichlet declarado en cada DOF prescrito
    y cero en el resto: el patrón natural de un ``dirichlet_func``."""
    base = np.zeros(dom.total_dofs)
    for node in dom.nodes.values():
        for d, v in node.boundary_conditions.items():
            base[node.dofs[d]] = v
    return base


def al_modelo_atado(dom_atado, dom_fundido, campo_fundido: np.ndarray) -> np.ndarray:
    """Lleva un campo del modelo fundido al atado: cada duplicado toma el valor
    de su original."""
    campo_fundido = np.asarray(campo_fundido)
    out = np.zeros((dom_atado.total_dofs,) + campo_fundido.shape[1:])
    for n, node in dom_atado.nodes.items():
        origen = n - DUP if n > DUP else n
        out[node.dofs["T"]] = campo_fundido[dom_fundido.nodes[origen].dofs["T"]]
    return out


def _transitorio(builder, theta, lumping, *, dt, n_steps=12, T_initial=20.0,
                 escala_dirichlet=None, **kw_modelo):
    """Corre el modelo fundido y el atado; devuelve ``(ref, T_atado, doms)``."""
    out = {}
    for atado in (False, True):
        dom, nid, dup = builder(atado, **kw_modelo)
        F = calor_borde_derecho(dom, nid)
        kw = {}
        if escala_dirichlet is not None:
            base = valores_declarados(dom)
            kw["dirichlet_func"] = lambda t, b=base: b * escala_dirichlet(t)
        r = ThetaMethodSolver(
            Assembler(dom), dt=dt, n_steps=n_steps, T_initial=T_initial,
            theta=theta, lumping=lumping, F_func=lambda t, F=F: F, **kw,
        ).solve()
        assert r.converged
        out[atado] = (dom, nid, dup, r)
    dom_f, _, _, r_f = out[False]
    dom_a, nid_a, dup_a, r_a = out[True]
    ref = al_modelo_atado(dom_a, dom_f, r_f.T_history)
    return ref, r_a.T_history, (dom_a, nid_a, dup_a)


# ----------------------------------------------------------------------
# ConstraintSet: Dirichlet declarado y g para valores de Dirichlet dados
# ----------------------------------------------------------------------

class TestConstraintSetDirichlet:

    def test_dirichlet_dofs_son_los_declarados(self):
        """``add_dirichlet`` y ``add_linear`` sin maestros cuentan; las MPC no."""
        cs = ConstraintSet()
        cs.add_dirichlet(4, 1.0)
        cs.add_linear(0, masters=[], coefficients=[], g=2.0)
        cs.add_linear(2, masters=[3], coefficients=[1.0])
        np.testing.assert_array_equal(cs.dirichlet_dofs, [0, 4])
        cs.build(6)
        np.testing.assert_array_equal(cs.dirichlet_dofs, [0, 4])

    def test_el_esclavo_atado_a_un_dirichlet_no_cuenta(self):
        """El cierre deja sin maestros al esclavo de un maestro con Dirichlet,
        pero ese DOF no es un Dirichlet declarado, antes ni después de ``build``."""
        cs = ConstraintSet()
        cs.add_dirichlet(1, 50.0)
        cs.add_linear(0, masters=[1], coefficients=[1.0])
        np.testing.assert_array_equal(cs.dirichlet_dofs, [1])
        _T, g = cs.build(3)
        assert g[0] == 50.0
        np.testing.assert_array_equal(cs.dirichlet_dofs, [1])

    def test_redeclarar_tras_build_es_idempotente(self):
        """La comparación es con lo declarado, no con la restricción cerrada."""
        cs = ConstraintSet()
        cs.add_linear(0, masters=[1], coefficients=[1.0])
        cs.add_linear(1, masters=[2], coefficients=[1.0])
        T1, g1 = cs.build(3)
        cs.add_linear(0, masters=[1], coefficients=[1.0])   # no debe lanzar
        T2, g2 = cs.build(3)
        np.testing.assert_array_equal(T1.toarray(), T2.toarray())
        np.testing.assert_array_equal(g1, g2)

    def test_build_g_con_los_valores_declarados_reproduce_build(self):
        cs = ConstraintSet()
        cs.add_dirichlet(5, 7.0)
        cs.add_linear(0, masters=[1], coefficients=[2.0], g=1.0)
        cs.add_linear(1, masters=[5, 3], coefficients=[3.0, -1.0], g=0.5)
        _T, g = cs.build(6)
        declarados = np.zeros(6)
        declarados[5] = 7.0
        np.testing.assert_allclose(cs.build_g(6, declarados), g, rtol=0, atol=1e-14)

    def test_build_g_compone_el_cierre(self):
        """``s0 = 1 + 2·s1``, ``s1 = 0.5 + 3·d − m`` con ``d`` Dirichlet:
        ``g_s0 = 1 + 2·(0.5 + 3·ḡ_d)`` para cualquier ``ḡ_d``."""
        cs = ConstraintSet()
        cs.add_dirichlet(5, 7.0)
        cs.add_linear(0, masters=[1], coefficients=[2.0], g=1.0)
        cs.add_linear(1, masters=[5, 3], coefficients=[3.0, -1.0], g=0.5)
        valores = np.zeros(6)
        valores[5] = -4.0
        g = cs.build_g(6, valores)
        assert g[5] == -4.0
        assert g[1] == pytest.approx(0.5 + 3.0 * -4.0, abs=1e-14)
        assert g[0] == pytest.approx(1.0 + 2.0 * (0.5 + 3.0 * -4.0), abs=1e-14)
        assert g[2] == g[3] == g[4] == 0.0

    def test_build_g_solo_lee_los_dirichlet(self):
        cs = ConstraintSet()
        cs.add_dirichlet(2, 1.0)
        cs.add_linear(0, masters=[2], coefficients=[1.0], g=3.0)
        valores = np.array([99.0, 99.0, 10.0, 99.0])
        np.testing.assert_array_equal(cs.build_g(4, valores), [13.0, 0.0, 10.0, 0.0])

    def test_build_g_sin_mpc_es_exacto(self):
        """Sin MPC ``build_g`` es la selección exacta de los valores Dirichlet."""
        cs = ConstraintSet()
        for d, v in ((0, 0.1), (3, 1.0 / 3.0)):
            cs.add_dirichlet(d, v)
        valores = np.array([np.pi, 5.0, 6.0, np.e, 8.0])
        g = cs.build_g(5, valores)
        assert g[0] == np.pi and g[3] == np.e
        assert g[1] == g[2] == g[4] == 0.0

    def test_build_g_rechaza_un_tamano_incorrecto(self):
        cs = ConstraintSet()
        cs.add_dirichlet(0, 1.0)
        with pytest.raises(ValueError, match="dirichlet_values"):
            cs.build_g(4, np.zeros(3))


# ----------------------------------------------------------------------
# Oráculo del nodo fundido: transitorio y estacionario, 2D y 3D
# ----------------------------------------------------------------------

class TestOraculoNodoFundido:

    @pytest.mark.parametrize("theta, lumping, dt", [
        (1.0, "lumped", 5.0), (1.0, "consistent", 5.0),
        (0.5, "lumped", 5.0), (0.5, "consistent", 5.0),
        (0.0, "lumped", 0.5),
    ])
    def test_transitorio_2d(self, theta, lumping, dt):
        ref, T_a, _ = _transitorio(placa_2d, theta, lumping, dt=dt)
        np.testing.assert_allclose(T_a, ref, rtol=0, atol=1e-11)

    @pytest.mark.parametrize("theta, lumping, dt", [
        (1.0, "lumped", 5.0), (0.5, "consistent", 5.0), (0.0, "lumped", 0.5),
    ])
    def test_transitorio_3d(self, theta, lumping, dt):
        ref, T_a, _ = _transitorio(placa_3d, theta, lumping, dt=dt)
        np.testing.assert_allclose(T_a, ref, rtol=0, atol=1e-11)

    @pytest.mark.parametrize("builder", [placa_2d, placa_3d], ids=["2d", "3d"])
    def test_estacionario(self, builder):
        out = {}
        for atado in (False, True):
            dom, nid, _ = builder(atado)
            U = LinearSolver(Assembler(dom)).solve(calor_borde_derecho(dom, nid))
            out[atado] = (dom, U)
        (dom_f, U_f), (dom_a, U_a) = out[False], out[True]
        np.testing.assert_allclose(U_a, al_modelo_atado(dom_a, dom_f, U_f),
                                   rtol=0, atol=1e-11)

    @pytest.mark.parametrize("builder", [placa_2d, placa_3d], ids=["2d", "3d"])
    def test_dirichlet_variable_sobre_el_maestro_de_la_interfaz(self, builder):
        """``dirichlet_func`` con sólo los valores Dirichlet (el patrón natural,
        desde ``boundary_conditions``) y un Dirichlet sobre un maestro de MPC:
        el esclavo atado sigue al maestro en todos los instantes, t = 0 incluido.
        Hasta el 2026-09-30 el esclavo quedaba a 0."""
        ref, T_a, (dom_a, nid_a, dup_a) = _transitorio(
            builder, 1.0, "lumped", dt=5.0, dirichlet_interfaz=True,
            escala_dirichlet=lambda t: 1.0 + 0.01 * t,
        )
        np.testing.assert_allclose(T_a, ref, rtol=0, atol=1e-11)
        clave = next(iter(k for k in dup_a if all(c == 0 for c in k[1:])))
        m = dom_a.nodes[nid_a[clave]].dofs["T"]
        s = dom_a.nodes[dup_a[clave]].dofs["T"]
        np.testing.assert_array_equal(T_a[s], T_a[m])
        assert T_a[m, -1] == pytest.approx(50.0 * (1.0 + 0.01 * 60.0))


# ----------------------------------------------------------------------
# dirichlet_func con MPC
# ----------------------------------------------------------------------

class TestDirichletFuncConMPC:

    def test_el_salto_g_de_la_mpc_se_conserva(self):
        """MPC ``T_s = T_m + 5`` y ``dirichlet_func`` constante construido desde
        ``boundary_conditions``: el resultado es el del camino sin
        ``dirichlet_func`` y el salto se mantiene. Hasta el 2026-09-30 el salto
        se perdía (``T_s − T_m = 0``)."""
        res = []
        for con_func in (False, True):
            dom, nid, dup = placa_2d(True, mpc_g=5.0)
            kw = {}
            if con_func:
                base = valores_declarados(dom)
                kw["dirichlet_func"] = lambda t, b=base: b
            r = ThetaMethodSolver(
                Assembler(dom), dt=5.0, n_steps=10, T_initial=20.0, **kw,
            ).solve()
            res.append((dom, nid, dup, r))
        (_, _, _, r_sin), (dom, nid, dup, r_con) = res
        np.testing.assert_allclose(r_con.T_history, r_sin.T_history, rtol=0, atol=1e-11)
        for j in range(4):
            s = dom.nodes[dup[3, j]].dofs["T"]
            m = dom.nodes[nid[3, j]].dofs["T"]
            np.testing.assert_allclose(r_con.T_history[s] - r_con.T_history[m],
                                       5.0, rtol=0, atol=1e-11)

    def test_un_valor_en_un_dof_libre_se_rechaza(self):
        """Prescribir un DOF exige declararlo: un valor no nulo en un DOF libre
        no se ignora en silencio."""
        dom, nid, _ = placa_2d(True)
        base = valores_declarados(dom)
        base[dom.nodes[nid[5, 1]].dofs["T"]] = 30.0
        s = ThetaMethodSolver(Assembler(dom), dt=5.0, n_steps=2, T_initial=20.0,
                              dirichlet_func=lambda t: base)
        with pytest.raises(ValueError, match="sin Dirichlet declarado"):
            s.solve()

    def test_la_columna_inicial_usa_dirichlet_func_en_t0(self):
        """Declarado 100, ``dirichlet_func ≡ 300``: el instante inicial
        almacenado vale 300 en el borde y en el esclavo atado al maestro con
        Dirichlet, y el aviso de ``T_initial`` se evalúa contra 300. Hasta el
        2026-09-30 la columna 0 mostraba el valor declarado."""
        dom, nid, dup = placa_2d(True, dirichlet_interfaz=True)
        base = valores_declarados(dom)
        base[base == T_IZQ] = 300.0
        with capturar_avisos() as log:
            r = ThetaMethodSolver(
                Assembler(dom), dt=5.0, n_steps=3, T_initial=20.0,
                dirichlet_func=lambda t: base,
            ).solve()
        izq = [dom.nodes[nid[0, j]].dofs["T"] for j in range(4)]
        np.testing.assert_array_equal(r.T_history[izq, 0], 300.0)
        s = dom.nodes[dup[3, 0]].dofs["T"]
        np.testing.assert_array_equal(r.T_history[s], 50.0)
        aviso = log.con("T_initial contradice")
        assert len(aviso) == 1 and "en 5 DOF(s)" in aviso[0]

        # Con T_initial igual a dirichlet_func(0) en el borde no hay choque.
        T0 = np.full(dom.total_dofs, 300.0)
        T0[dom.nodes[nid[3, 0]].dofs["T"]] = 50.0
        with capturar_avisos() as log:
            ThetaMethodSolver(Assembler(dom), dt=5.0, n_steps=3, T_initial=T0,
                              dirichlet_func=lambda t: base).solve()
        assert not log.con("T_initial contradice")


# ----------------------------------------------------------------------
# Aviso de T_initial y guarda «sin Dirichlet»
# ----------------------------------------------------------------------

class TestAvisosYGuarda:

    @pytest.mark.parametrize("dirichlet_interfaz", [False, True])
    def test_el_aviso_de_T_initial_cita_solo_los_dirichlet(self, dirichlet_interfaz):
        """El modelo atado y el fundido emiten el mismo aviso: los esclavos de
        MPC no son Dirichlet contradichos. Hasta el 2026-09-30 el atado citaba
        también los esclavos (8 o 9 DOF frente a 4 o 5)."""
        avisos = {}
        for atado in (False, True):
            dom, _, _ = placa_2d(atado, dirichlet_interfaz=dirichlet_interfaz)
            with capturar_avisos() as log:
                ThetaMethodSolver(Assembler(dom), dt=5.0, n_steps=2,
                                  T_initial=20.0).solve()
            avisos[atado] = log.con("T_initial contradice")
        n = 5 if dirichlet_interfaz else 4
        assert len(avisos[False]) == 1 and f"en {n} DOF(s)" in avisos[False][0]
        assert avisos[True] == avisos[False]

    @pytest.mark.parametrize("theta, lumping", [(1.0, "lumped"), (0.5, "consistent")])
    def test_el_modelo_atado_sin_dirichlet_coincide_con_el_fundido(self, theta, lumping):
        """Sin Dirichlet el transitorio está bien planteado (``C`` definida
        positiva) y el oráculo del nodo fundido se cumple: el calor que entra
        por el borde derecho calienta igual la placa atada y la conforme. Hasta
        el 2026-09-30 el atado pasaba la guarda «sin Dirichlet» y el fundido
        no; hasta el 2026-10-20 ambos se rechazaban."""
        ref, T_a, _ = _transitorio(placa_2d, theta, lumping, dt=5.0, T_izq=None)
        np.testing.assert_allclose(T_a, ref, rtol=0, atol=1e-11)

    def test_una_restriccion_sin_maestros_cuenta_como_dirichlet(self):
        """``u_s = g`` declarada como restricción lineal sin maestros es un
        Dirichlet: el modelo corre y coincide con el que usa ``fix_dof``."""
        res = []
        for via_mpc in (False, True):
            dom, nid, _ = placa_2d(False, T_izq=None if via_mpc else T_IZQ)
            if via_mpc:
                for j in range(4):
                    dom.add_linear_constraint((nid[0, j], "T"), [], [], g=T_IZQ)
            F = calor_borde_derecho(dom, nid)
            with capturar_avisos() as log:
                r = ThetaMethodSolver(Assembler(dom), dt=5.0, n_steps=5,
                                      T_initial=20.0, F_func=lambda t, F=F: F).solve()
            res.append((r, log.con("T_initial contradice")))
        (r_fix, av_fix), (r_mpc, av_mpc) = res
        np.testing.assert_array_equal(r_mpc.T_history, r_fix.T_history)
        assert av_mpc == av_fix and len(av_fix) == 1

    def test_con_theta_0_el_aviso_dice_estabilidad_condicional(self):
        """θ = 0 no es A-estable: el aviso de choque inicial no debe decirlo."""
        dom, _, _ = placa_2d(False)
        with capturar_avisos() as log:
            ThetaMethodSolver(Assembler(dom), dt=0.5, n_steps=2, T_initial=20.0,
                              theta=0.0).solve()
        aviso = log.con("T_initial contradice")
        assert len(aviso) == 1
        assert "condicionalmente estable" in aviso[0]
        assert "A-estable pero" not in aviso[0]

    def test_con_theta_medio_el_aviso_sigue_diciendo_A_estable(self):
        dom, _, _ = placa_2d(False)
        with capturar_avisos() as log:
            ThetaMethodSolver(Assembler(dom), dt=5.0, n_steps=2, T_initial=20.0,
                              theta=0.5).solve()
        aviso = log.con("T_initial contradice")
        assert len(aviso) == 1 and "A-estable pero no L-estable" in aviso[0]
