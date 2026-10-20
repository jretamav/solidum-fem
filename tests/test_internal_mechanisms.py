"""Capa 3 de la red de seguridad: mecanismos internos por la energía del modo
(ADR 0019, ampliación del 2026-10-20, decisión del usuario).

Hasta el 2026-10-20 un mecanismo interno se declaraba por el tamaño de un
pivote, ``|p| < 1e-13·max|K|``. Medido entonces: el pivote de un mecanismo
real es ruido de redondeo que crece con el tamaño del modelo (1e-16 con 100
DOF, 2.4e-11 con 264 000), así que los mecanismos de modelos medianos pasaban
sin aviso (desde unos 4 000 DOF con SuperLU y 30 000 con Pardiso), y el
veredicto dependía de las unidades (un pórtico con zonas rígidas se rechazaba
en mm y no en m). Ahora un pivote pequeño frente a la diagonal sólo dispara
la búsqueda del modo de menor energía de ``K·v = λ·diag(K)·v``, y decide la
energía de ese modo en cada elemento que mueve.

1. Mecanismos que el criterio anterior dejaba pasar, localizados y descritos
   igual con SuperLU y con Pardiso.
2. Descripción: giro alrededor de un punto, traslación más giro, rótula
   esférica, campo escalar sin valor prescrito, tangente inicial.
3. Contrastes de rigidez extremos que no son mecanismos: zonas rígidas hasta
   ×1e10 en m y en mm, ν = 0.4999999.
4. El filtro no depende de las unidades, y sin pivotes pequeños la
   comprobación no cuesta nada.
5. Pardiso lee sus pivotes (``pardiso_getdiag``): detecta un mecanismo de
   28 560 DOF que su conteo de pivotes perturbados no ve.
6. Sin factorización utilizable (pivote exactamente nulo, o tan pequeño que
   la resolución no es finita) el modo se calcula con ``K + σ·diag(K)``: el
   veredicto y la descripción no dependen de que el redondeo deje el pivote
   en cero o en 1e-17, que cambia con la versión de SciPy (el CI con Python
   3.11 lo destapó el 2026-10-20).
"""
import logging
import os
import sys
import unittest
from unittest import mock

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.dirname(__file__))

from _malla_estructurada import Malla, empotrar  # noqa: E402

from solidum.materials.elastic import Elastic1D  # noqa: E402
from solidum.materials.thermal_conduction import ThermalConduction  # noqa: E402
from solidum.math.assembly import Assembler  # noqa: E402
from solidum.math.linalg import LUSolver  # noqa: E402
from solidum.math.linalg.dispatcher import _HAS_PARDISO  # noqa: E402
from solidum.math.solvers import LinearSolver, NonlinearSolver, model_checks  # noqa: E402
from solidum.math.solvers.diagnostics import IllPosedSystemError  # noqa: E402
from solidum.math.solvers.model_checks import find_internal_mechanisms  # noqa: E402
from solidum.registry import MaterialRegistry  # noqa: E402

DIRECT = ("lu",) + (("pardiso",) if _HAS_PARDISO else ())
E_ACERO = 210e9
HORMIGON = dict(E=30e9, nu=0.2)


def setUpModule():
    logging.disable(logging.WARNING)


def tearDownModule():
    logging.disable(logging.NOTSET)


def backend(nombre):
    if nombre == "lu":
        return LUSolver()
    from solidum.math.linalg.pardiso import PardisoSolver
    return PardisoSolver()


def factorizar(dom, nombre):
    """``(assembler, K reducida, factor)`` de la rigidez en ``U = 0``."""
    a = Assembler(dom)
    a.assemble_system()
    K, _, _, _ = a.reduce(a.K_global.tocsr(), np.zeros(dom.total_dofs))
    return a, K, backend(nombre).factorize(K)


# --- mecanismos ---------------------------------------------------------------

def placas_unidas_por_un_nodo(n):
    """Placa A [0,1]² empotrada en x = 0; placa B [1,2]×[1,2] que sólo comparte
    con A el nodo (1, 1). Carga sobre A: no activa el giro de B."""
    m = Malla()
    mat = MaterialRegistry.get("Elastic2D")(**HORMIGON)
    m.caja("Quad4", (0, 0), (1, 1), (n, n), mat)
    m.caja("Quad4", (1, 1), (1, 1), (n, n), mat)
    empotrar(m.nodos(lambda x: abs(x[0]) < 1e-9), ("ux", "uy"))
    return m.cerrar([(m.nodo((1.0, 0.5)), "uy", -1e5)])


def bloques_unidos(n, origen_b):
    """Bloque A [0,1]³ empotrado en x = 0 y bloque B de lado 1 con origen
    ``origen_b``, que comparte con A una arista o un nodo."""
    m = Malla()
    mat = MaterialRegistry.get("Elastic3D")(**HORMIGON)
    m.caja("Hex8", (0, 0, 0), (1, 1, 1), (n, n, n), mat)
    m.caja("Hex8", origen_b, (1, 1, 1), (n, n, n), mat)
    empotrar(m.nodos(lambda x: abs(x[0]) < 1e-9), ("ux", "uy", "uz"))
    return m.cerrar([(m.nodo((0.5, 0.5, 1.0)), "uz", -1e5)])


def placa_atada_por_una_mpc():
    """Placa B atada a la A sólo por el ux de un par de nodos: puede deslizar
    en vertical y girar alrededor de cualquier punto de la recta y = 0.5."""
    m = Malla()
    mat = MaterialRegistry.get("Elastic2D")(**HORMIGON)
    m.caja("Quad4", (0, 0), (1, 1), (2, 2), mat)
    m.caja("Quad4", (1.5, 0), (1, 1), (2, 2), mat)
    empotrar(m.nodos(lambda x: abs(x[0]) < 1e-9), ("ux", "uy"))
    a, b = m.nodo((1.0, 0.5)), m.nodo((1.5, 0.5))
    m.dom.add_linear_constraint((b.id, "ux"), [(a.id, "ux")], [1.0])
    return m.cerrar([(a, "ux", 1e4)])


def placas_termicas_desconectadas():
    """Dos placas térmicas sin nodos comunes; temperatura impuesta sólo en la
    primera. La capa 1 no lo ve: el campo T tiene valores prescritos."""
    m = Malla()
    mat = ThermalConduction(k=45.0, dim=2)
    m.caja("Quad4Thermal", (0, 0), (1, 1), (2, 2), mat)
    m.caja("Quad4Thermal", (1.5, 0), (1, 1), (2, 2), mat)
    for nodo in m.nodos(lambda x: abs(x[0]) < 1e-9):
        nodo.fix_dof("T", 20.0)
    return m.cerrar()


# --- bien planteados con contrastes extremos ------------------------------------

def portico_con_zonas_rigidas(contraste, unidades="m"):
    """Pórtico empotrado con zonas rígidas de 0,3 m en los nudos viga-columna
    (Frame2DEuler de A·c e I·c), en N y m o en N y mm."""
    f = 1e3 if unidades == "mm" else 1.0
    E = E_ACERO / f**2
    A, I = 5e-3 * f**2, 8e-5 * f**4
    m = Malla()
    mat = Elastic1D(E=E)

    def tramo(a, b, c=1.0, n=1):
        a, b = f * np.asarray(a, float), f * np.asarray(b, float)
        for k in range(n):
            m.elemento("Frame2DEuler", [a + (b - a) * k / n, a + (b - a) * (k + 1) / n],
                       mat, A=A * c, I=I * c)
    tramo((0, 0), (0, 4), n=4)
    tramo((6, 0), (6, 4), n=4)
    tramo((0, 4), (0.3, 4), c=contraste)
    tramo((0.3, 4), (5.7, 4), n=4)
    tramo((5.7, 4), (6, 4), c=contraste)
    empotrar([m.nodo((0, 0)), m.nodo((6 * f, 0))], ("ux", "uy", "rz"))
    return m.cerrar()


def placa_casi_incompresible(nu):
    m = Malla()
    m.caja("Quad4", (0, 0), (2, 1), (16, 8),
           MaterialRegistry.get("Elastic2D")(E=30e9, nu=nu, hypothesis="plane_strain"))
    empotrar(m.nodos(lambda x: abs(x[0]) < 1e-9), ("ux", "uy"))
    borde = m.nodos(lambda x: abs(x[0] - 2) < 1e-9)
    return m.cerrar([(n, "uy", -1e4 / len(borde)) for n in borde])


# ----------------------------------------------------------------------
# 1-2. Mecanismos: detectados, localizados y descritos
# ----------------------------------------------------------------------

class TestMecanismosLocalizados(unittest.TestCase):

    def rechazo(self, modelo, la, solver=LinearSolver, **kw):
        dom, F = modelo()
        with self.assertRaises(IllPosedSystemError) as cm:
            solver(Assembler(dom), linear_algebra=la, **kw).solve(F)
        return cm.exception

    def test_quad4_de_4288_dof_unidas_por_un_nodo(self):
        """El criterio anterior lo dejaba pasar con SuperLU (pivote 3.5e-12·max|K|)."""
        for la in DIRECT:
            with self.subTest(backend=la):
                err = self.rechazo(lambda: placas_unidas_por_un_nodo(32), la)
                self.assertEqual(err.zero_pivots, 1)
                self.assertEqual(err.mechanisms, [
                    "la parte formada por los elementos 1025–2048 puede moverse sin "
                    "deformarse: giro alrededor del eje z que pasa por (1, 1)"])
                self.assertIn("LinearSolver", str(err))

    def test_hex8_de_4104_dof_unidos_por_una_arista(self):
        """El criterio anterior lo dejaba pasar con una LDLᵀ tipo CHOLMOD."""
        for la in DIRECT:
            with self.subTest(backend=la):
                err = self.rechazo(lambda: bloques_unidos(8, (1, 0, 1)), la)
                self.assertEqual(err.zero_pivots, 1)
                self.assertEqual(err.mechanisms, [
                    "la parte formada por los elementos 513–1024 puede moverse sin "
                    "deformarse: giro alrededor del eje y que pasa por (1, 0.5, 1)"])

    def test_bloque_unido_por_un_nodo_es_una_rotula_esferica(self):
        for la in DIRECT:
            with self.subTest(backend=la):
                err = self.rechazo(lambda: bloques_unidos(2, (1, 1, 1)), la)
                self.assertEqual(err.zero_pivots, 3)
                self.assertEqual(err.mechanisms, [
                    "la parte formada por los elementos 9–16 puede moverse sin "
                    "deformarse de 3 formas independientes: giro alrededor de "
                    "cualquier eje que pasa por (1, 1, 1) (como una rótula esférica)"])

    def test_placa_atada_por_una_mpc_desliza_y_gira(self):
        """Dos movimientos libres de la misma parte: base legible (traslación
        pura y giro sin traslación libre), la misma con los dos backends."""
        for la in DIRECT:
            with self.subTest(backend=la):
                err = self.rechazo(placa_atada_por_una_mpc, la)
                self.assertEqual(err.zero_pivots, 2)
                self.assertEqual(err.mechanisms, [
                    "la parte formada por los elementos 5–8 puede moverse sin "
                    "deformarse de 2 formas independientes: traslación en la "
                    "dirección y; giro alrededor del eje z que pasa por (2, 0.5)"])

    def test_parte_termica_sin_valor_prescrito(self):
        for la in DIRECT:
            with self.subTest(backend=la):
                err = self.rechazo(placas_termicas_desconectadas, la)
                self.assertEqual(err.mechanisms, [
                    "la parte formada por los elementos 5–8 no está unida a ningún "
                    "valor prescrito del campo 'T': su nivel queda indeterminado"])

    def test_tangente_inicial_del_no_lineal(self):
        for la in DIRECT:
            with self.subTest(backend=la):
                err = self.rechazo(lambda: placas_unidas_por_un_nodo(16), la,
                                   solver=NonlinearSolver, num_steps=1)
                self.assertTrue(err.initial_tangent)
                self.assertEqual(err.mechanisms, [
                    "la parte formada por los elementos 257–512 puede moverse sin "
                    "deformarse: giro alrededor del eje z que pasa por (1, 1)"])


# ----------------------------------------------------------------------
# 3. Contrastes de rigidez extremos: no son mecanismos
# ----------------------------------------------------------------------

class TestContrastesQueNoSonMecanismos(unittest.TestCase):

    def test_zonas_rigidas_en_m_y_en_mm(self):
        """El filtro salta (hay pivotes pequeños) y la búsqueda del modo se hace,
        pero el modo deforma las columnas y la viga: no es un mecanismo.
        Con el criterio anterior, c = 1e7 en mm daba dos pivotes nulos."""
        for c in (1e6, 1e8, 1e10):
            for unidades in ("m", "mm"):
                for la in DIRECT:
                    with self.subTest(contraste=c, unidades=unidades, backend=la):
                        dom, _ = portico_con_zonas_rigidas(c, unidades)
                        a, K, factor = factorizar(dom, la)
                        if c >= 1e8:
                            self.assertGreater(factor.n_small_pivots, 0)
                        self.assertEqual(find_internal_mechanisms(a, K, factor), [])

    def test_material_casi_incompresible(self):
        """``1 − 2ν = 2e-7``: el modo de menor energía es isócoro y deforma
        todos los elementos (energía relativa ~1e-7). Sólo la capa 3: la
        capa 2 (equilibrio, ``‖F − K·u‖ ≤ 1e-8·‖F‖``) rechaza hoy este modelo
        con un residuo de 3e-8, el falso positivo pendiente de decisión
        (ADR 0019, «Lo que queda sin cubrir»)."""
        for la in DIRECT:
            with self.subTest(backend=la):
                dom, _ = placa_casi_incompresible(0.4999999)
                a, K, factor = factorizar(dom, la)
                self.assertEqual(find_internal_mechanisms(a, K, factor), [])


# ----------------------------------------------------------------------
# 4. Filtro: unidades y coste
# ----------------------------------------------------------------------

class TestFiltro(unittest.TestCase):

    def test_el_filtro_salta_en_m_y_en_mm(self):
        """El cociente pivote/diagonal es invariante a las unidades para una
        secuencia de pivotes dada, pero el pivoteo parcial de SuperLU puede
        elegir otra al cambiar la escala: con c = 1e8 da 6 pivotes pequeños
        en m y 2 en mm. El filtro salta en las dos, y el veredicto, que es lo
        que decide, es el mismo (``test_zonas_rigidas_en_m_y_en_mm``). Con
        el criterio anterior (``max|K|``) en m no saltaba nada y en mm daba
        dos pivotes nulos y un rechazo."""
        for unidades in ("m", "mm"):
            with self.subTest(unidades=unidades):
                dom, _ = portico_con_zonas_rigidas(1e8, unidades)
                self.assertGreater(factorizar(dom, "lu")[2].n_small_pivots, 0)

    def test_sin_pivotes_pequenos_no_se_busca_el_modo(self):
        """Un modelo normal no paga la comprobación: ni una resolución extra ni
        la pasada por elementos."""
        for la in DIRECT:
            with self.subTest(backend=la):
                dom, F = placa_casi_incompresible(0.3)
                with mock.patch.object(model_checks, "_element_blocks",
                                       side_effect=AssertionError("no debía llamarse")):
                    LinearSolver(Assembler(dom), linear_algebra=la).solve(F)


# ----------------------------------------------------------------------
# 5. Pardiso lee sus pivotes
# ----------------------------------------------------------------------

@unittest.skipUnless(_HAS_PARDISO, "pypardiso no instalado")
class TestPivotesDePardiso(unittest.TestCase):

    def test_detecta_lo_que_su_conteo_de_perturbados_no_ve(self):
        """Hex8 de 28 560 DOF: Pardiso no perturba ningún pivote (su umbral es
        1e-13 en la matriz escalada; el del mecanismo es 5.9e-13), así que el
        criterio anterior lo dejaba pasar. Con la diagonal del factor, el
        filtro salta y la energía del modo lo confirma."""
        from solidum.math.linalg.pardiso import PardisoSolver
        dom, _ = bloques_unidos(16, (1, 0, 1))
        a, K, factor = factorizar(dom, "pardiso")
        self.assertIsNotNone(PardisoSolver._getdiag)
        self.assertEqual(PardisoSolver._handle().get_iparm(14), 0)
        self.assertGreaterEqual(factor.n_small_pivots, 1)
        encontrados = find_internal_mechanisms(a, K, factor)
        self.assertEqual([m.dimension for m in encontrados], [1])
        self.assertIn("giro alrededor del eje y que pasa por (1, 0.5, 1)",
                      encontrados[0].description)


# ----------------------------------------------------------------------
# 6. Sin factorización utilizable: K + σ·diag(K)
# ----------------------------------------------------------------------

class _FactorNoFinito:
    """Factorización con un pivote del orden de 1e-300: resuelve con valores
    no finitos."""
    n_small_pivots = 1

    def solve(self, b):
        return np.full_like(b, np.inf)


class TestSinFactorizacionUtilizable(unittest.TestCase):

    CASOS = {
        "térmica": (placas_termicas_desconectadas,
                    "la parte formada por los elementos 5–8 no está unida a ningún "
                    "valor prescrito del campo 'T': su nivel queda indeterminado"),
        "Quad4": (lambda: placas_unidas_por_un_nodo(8),
                  "la parte formada por los elementos 65–128 puede moverse sin "
                  "deformarse: giro alrededor del eje z que pasa por (1, 1)"),
    }

    def test_mismo_veredicto_sin_factor_y_con_factor_no_finito(self):
        for nombre, (modelo, descripcion) in self.CASOS.items():
            dom, _ = modelo()
            a, K, factor = factorizar(dom, "lu")
            for etiqueta, f in (("con factor", factor), ("sin factor", None),
                                ("factor no finito", _FactorNoFinito())):
                with self.subTest(modelo=nombre, caso=etiqueta):
                    encontrados = find_internal_mechanisms(a, K, f)
                    self.assertEqual([(m.dimension, m.description) for m in encontrados],
                                     [(1, descripcion)])


if __name__ == "__main__":
    unittest.main()
