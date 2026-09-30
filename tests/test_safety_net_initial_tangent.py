"""Red de seguridad del ADR 0019 en los solvers estáticos no lineales y con el
backend iterativo (cierres del 2026-09-30).

Hasta el 2026-09-30 un **mecanismo interno que la carga no activa** pasaba sin
aviso fuera del estático lineal con solver directo: el Newton resolvía el
sistema singular pero consistente (Pardiso perturba el pivote nulo) y
devolvía unos desplazamientos con una parte arbitraria del mecanismo. Estas
pruebas fijan el comportamiento nuevo, sobre la armadura del ejemplo 11 del
manual de ejemplos (``_armadura_pratt.py``):

1. **Tangente inicial** (decisión del usuario del 2026-09-30):
   ``NonlinearSolver``, ``ArcLengthSolver`` e ``IndirectDisplacementSolver``
   cuentan los pivotes nulos de la tangente en ``U = 0``.

   - Modelo geométricamente lineal (``Truss2D``): la no linealidad del
     material sólo reduce la rigidez, así que el mecanismo lo es en todo el
     análisis ⇒ ``IllPosedSystemError``, cargado o no.
   - Con elementos geométricamente no lineales (``Truss2DCorot``) o materiales
     unilaterales (fibras de cable), la rigidez puede crecer al deformarse ⇒
     sólo un aviso.
   - Modelo bien planteado, con material lineal o no lineal: nada.
   - La comprobación reutiliza la primera factorización del solver: el
     número de factorizaciones no cambia.

2. **Backend iterativo**: no tiene pivotes que inspeccionar; lo dice con un
   aviso, una vez por análisis, en ``LinearSolver`` y en los no lineales. Con
   un backend directo no hay aviso.
3. **Factorización imposible** (pivote exactamente nulo, un nudo colgante):
   ``IllPosedSystemError`` también en el estático lineal, en vez del
   ``RuntimeError`` del backend.
4. **Familias de modelos** (revisión adversaria del 2026-09-30): la
   comprobación rechaza los mecanismos internos no activados en sólidos 3D,
   marcos con rótula interna y modelos con MPC, y calla en modelos bien
   planteados con plasticidad perfecta, daño, fibras, dominios mixtos, MPC,
   apoyos oblicuos y zonas rígidas ×1e6 en unidades SI.
"""
import contextlib
import logging
import os
import sys
import unittest
from unittest import mock

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.dirname(__file__))

from _armadura_pratt import AREA, E, P, armadura, nodo_de_control  # noqa: E402

from solidum.core.domain import Domain  # noqa: E402
from solidum.elements.frame import Frame2DEuler  # noqa: E402
from solidum.elements.truss import Truss2D, Truss2DCorot  # noqa: E402
from solidum.materials.cable_1d import CableMaterial1D  # noqa: E402
from solidum.materials.elastic import Elastic1D  # noqa: E402
from solidum.materials.plastic_1d import Elastoplastic1D  # noqa: E402
from solidum.math.assembly import Assembler  # noqa: E402
from solidum.math.linalg import LUSolver  # noqa: E402
from solidum.math.linalg.dispatcher import _HAS_PARDISO  # noqa: E402
from solidum.math.solvers import (  # noqa: E402
    ArcLengthSolver,
    IndirectDisplacementSolver,
    LinearSolver,
    NonlinearSolver,
)
from solidum.math.solvers import nonlinear as nonlinear_module  # noqa: E402
from solidum.math.solvers.diagnostics import IllPosedSystemError  # noqa: E402
from solidum.sections.fiber import Fiber, FiberSection  # noqa: E402

DIRECT = ("lu",) + (("pardiso",) if _HAS_PARDISO else ())
AVISO_TANGENTE = "matriz tangente inicial (U = 0)"
AVISO_ITERATIVO = "no se comprueba si el modelo tiene un mecanismo interno"


def _no_lineales(la):
    """Los tres solvers estáticos no lineales del programa principal."""
    return {
        "NonlinearSolver": lambda a, dom: NonlinearSolver(a, num_steps=2, linear_algebra=la),
        "ArcLengthSolver": lambda a, dom: ArcLengthSolver(a, max_lambda=1.0, linear_algebra=la),
        "IndirectDisplacementSolver": lambda a, dom: IndirectDisplacementSolver(
            a, max_lambda=1.0, linear_algebra=la,
            control=[(nodo_de_control(dom), "ux", 1.0)]),
    }


@contextlib.contextmanager
def avisos():
    """Registros WARNING del logger de los solvers durante el bloque.

    Handler propio sobre ``solidum.solvers`` (el logger ``solidum`` no
    propaga) y ``logging.disable`` anulado mientras dura: otro módulo de
    pruebas puede haberlo dejado activo (memoria del proyecto, ADR 0018)."""
    registros = []

    class _Lista(logging.Handler):
        def emit(self, record):
            registros.append(record.getMessage())

    logger = logging.getLogger("solidum.solvers")
    handler = _Lista(level=logging.WARNING)
    previo = logging.root.manager.disable
    logging.disable(logging.NOTSET)
    logger.addHandler(handler)
    try:
        yield registros
    finally:
        logger.removeHandler(handler)
        logging.disable(previo)


def _con(registros, texto):
    return [r for r in registros if texto in r]


# ----------------------------------------------------------------------
# 1. Tangente inicial de los solvers estáticos no lineales
# ----------------------------------------------------------------------

class TestMecanismoInternoEnLaTangenteInicial(unittest.TestCase):

    def test_geometricamente_lineal_se_rechaza_aunque_la_carga_no_lo_active(self):
        """Armadura sin diagonal con carga horizontal: la estática tiene
        solución y el Newton convergía; ahora se rechaza al factorizar K(0)."""
        for la in DIRECT:
            for nombre, crear in _no_lineales(la).items():
                with self.subTest(solver=nombre, backend=la):
                    dom, F = armadura(sin_diagonal=True, carga="horizontal")
                    with self.assertRaises(IllPosedSystemError) as cm:
                        crear(Assembler(dom), dom).solve(F)
                    err = cm.exception
                    self.assertTrue(err.initial_tangent)
                    self.assertGreaterEqual(err.zero_pivots, 1)
                    msg = str(err)
                    self.assertIn(nombre, msg)
                    self.assertIn("mecanismo interno", msg)
                    self.assertIn("análisis modal", msg)

    def test_geometricamente_lineal_con_carga_que_lo_activa(self):
        """Antes: bisecciones hasta ``SingularTangentError``. Ahora el mismo
        error de planteamiento, antes de iterar."""
        for la in DIRECT:
            for nombre, crear in _no_lineales(la).items():
                if nombre == "IndirectDisplacementSolver":
                    continue      # control sobre ux del rodillo: la carga vertical también lo mueve
                with self.subTest(solver=nombre, backend=la):
                    dom, F = armadura(sin_diagonal=True, carga="vertical")
                    with self.assertRaises(IllPosedSystemError):
                        crear(Assembler(dom), dom).solve(F)

    def test_corrotacional_solo_avisa(self):
        """Con ``Truss2DCorot`` la rigidez geométrica puede estabilizar el
        mecanismo al deformarse: aviso, una vez, y el análisis sigue."""
        for la in DIRECT:
            for nombre, crear in _no_lineales(la).items():
                with self.subTest(solver=nombre, backend=la):
                    dom, F = armadura(sin_diagonal=True, carga="horizontal",
                                      elemento=Truss2DCorot)
                    with avisos() as registros:
                        crear(Assembler(dom), dom).solve(F)
                    aviso = _con(registros, AVISO_TANGENTE)
                    self.assertEqual(len(aviso), 1, registros)
                    self.assertIn(nombre, aviso[0])
                    self.assertIn("elementos geométricamente no lineales (Truss2DCorot)",
                                  aviso[0])
                    self.assertIn("análisis modal", aviso[0])

    def test_corrotacional_sin_carga_que_lo_active_llega_a_la_solucion(self):
        """El Newton sobre la armadura corrotacional sin diagonal converge a
        la solución de la cuerda inferior: ux del rodillo = P·L/(EA)."""
        dom, F = armadura(sin_diagonal=True, carga="horizontal", elemento=Truss2DCorot)
        U = NonlinearSolver(Assembler(dom), num_steps=2, linear_algebra="lu").solve(F)
        ux = U[nodo_de_control(dom).dofs["ux"]]
        self.assertAlmostEqual(ux / (P * 12.0 / (E * AREA)), 1.0, places=3)

    def test_bien_planteado_no_avisa(self):
        for elemento in (Truss2D, Truss2DCorot):
            for carga in ("horizontal", "vertical"):
                for la in DIRECT:
                    for nombre, crear in _no_lineales(la).items():
                        with self.subTest(elemento=elemento.__name__, carga=carga,
                                          solver=nombre, backend=la):
                            dom, F = armadura(sin_diagonal=False, carga=carga,
                                              elemento=elemento)
                            with avisos() as registros:
                                U = crear(Assembler(dom), dom).solve(F)
                            self.assertEqual(_con(registros, AVISO_TANGENTE), [])
                            self.assertTrue(np.all(np.isfinite(U)))

    def test_material_no_lineal_bien_apoyado_no_avisa(self):
        """Cuerda inferior en fluencia con endurecimiento: la tangente pierde
        rigidez, pero la inicial es regular. Ni aviso ni error, y la fluencia
        se produce de verdad (desplazamiento mayor que el elástico)."""
        sigma_y = 0.5 * P / AREA            # la cuerda inferior lleva N = P
        mat = Elastoplastic1D(E=E, sigma_y=sigma_y, H=0.1 * E)
        for la in DIRECT:
            with self.subTest(backend=la):
                dom, F = armadura(sin_diagonal=False, carga="horizontal", material=mat)
                with avisos() as registros:
                    U = NonlinearSolver(Assembler(dom), num_steps=4, linear_algebra=la).solve(F)
                self.assertEqual(_con(registros, AVISO_TANGENTE), [])
                ux = U[nodo_de_control(dom).dofs["ux"]]
                self.assertGreater(ux, 1.5 * P * 12.0 / (E * AREA))

    def test_veredicto_segun_lo_que_puede_ganar_rigidez(self):
        """Lógica del veredicto, por la API pública (``check_initial_tangent``
        con un corrector que sólo recoge el inspector) y con factorizaciones
        simuladas, sobre un voladizo ``Frame2DEuler`` (geométricamente lineal)
        con sección por fibras.

        Sin tendones, una tangente inicial singular es un error. Con tendones
        de cable (``CableMaterial1D``: rigidez nula sin deformar, E al
        traccionarse) la rigidez de la sección puede crecer al deformarse,
        como la geométrica de un corrotacional: sólo aviso. En un modelo real
        la tangente de este voladizo es regular (el acero ya da EA y EI): el
        pivote nulo se simula para aislar la decisión del caso que la
        provoca."""
        from solidum.math.solvers.model_checks import check_initial_tangent
        from solidum.sections import RectangularPatch

        def voladizo(con_tendones):
            acero = Elastic1D(E=E)
            fibras = ([Fiber(0.1, 0.0, 1e-4, CableMaterial1D(E=E)),
                       Fiber(-0.1, 0.0, 1e-4, CableMaterial1D(E=E))] if con_tendones else [])
            seccion = FiberSection(patches=[RectangularPatch(acero, y=(-0.15, 0.15),
                                                             z=(-0.1, 0.1), ny=6, nz=1)],
                                   fibers=fibras or None)
            dom = Domain()
            n1, n2 = dom.add_node(1, [0.0, 0.0]), dom.add_node(2, [2.0, 0.0])
            dom.add_element(Frame2DEuler(1, [n1, n2], section=seccion))
            for dof in ("ux", "uy", "rz"):
                n1.fix_dof(dof, 0.0)
            dom.generate_equation_numbers(verbose=False)
            return Assembler(dom)

        class _Corrector:
            def inspect_next_factorization(self, inspector):
                self.inspector = inspector

        class _Factor:
            def __init__(self, n):
                self.n_zero_pivots = n

        def veredicto(con_tendones, factor):
            corrector = _Corrector()
            check_initial_tangent(voladizo(con_tendones), corrector, "NonlinearSolver")
            with avisos() as registros:
                corrector.inspector(factor)
            return registros

        # Regular: nada, con o sin tendones.
        for con in (False, True):
            self.assertEqual(veredicto(con, _Factor(0)), [])
        # Singular sin nada que pueda ganar rigidez: error.
        with self.assertRaises(IllPosedSystemError) as cm:
            veredicto(False, _Factor(2))
        self.assertEqual(cm.exception.zero_pivots, 2)
        with self.assertRaises(IllPosedSystemError) as cm:
            veredicto(False, None)                  # pivote exactamente nulo
        self.assertIsNone(cm.exception.zero_pivots)
        self.assertIn("exactamente nulo", str(cm.exception))
        # Singular con tendones: aviso que nombra el material.
        for factor in (_Factor(2), None):
            aviso = _con(veredicto(True, factor), AVISO_TANGENTE)
            self.assertEqual(len(aviso), 1)
            self.assertIn("materiales unilaterales (CableMaterial1D)", aviso[0])
        # Sin pivotes que inspeccionar (iterativo): el aviso del iterativo.
        for con in (False, True):
            registros = veredicto(con, _Factor(None))
            self.assertEqual(len(_con(registros, AVISO_ITERATIVO)), 1)
            self.assertEqual(_con(registros, AVISO_TANGENTE), [])

    def test_reutiliza_la_primera_factorizacion(self):
        """La comprobación no factoriza dos veces: el número de
        factorizaciones del análisis es el mismo con ella y sin ella."""
        def contar(activa):
            original = LUSolver.factorize
            llamadas = []

            def factorize(self, K):
                llamadas.append(K.shape[0])
                return original(self, K)

            dom, F = armadura(sin_diagonal=False, carga="vertical")
            parches = [mock.patch.object(LUSolver, "factorize", factorize)]
            if not activa:
                parches.append(mock.patch.object(nonlinear_module, "check_initial_tangent",
                                                 lambda *a, **k: None))
            with contextlib.ExitStack() as pila:
                for p in parches:
                    pila.enter_context(p)
                NonlinearSolver(Assembler(dom), num_steps=3, linear_algebra="lu").solve(F)
            return len(llamadas)

        con, sin = contar(True), contar(False)
        self.assertGreater(sin, 0)
        self.assertEqual(con, sin)

    def test_la_dinamica_no_comprueba(self):
        """En dinámica un mecanismo es legítimo (la masa lo regulariza): el
        Newton-Newmark no inspecciona la tangente inicial."""
        from solidum.math.solvers import NewtonNewmarkSolver
        dom, F = armadura(sin_diagonal=True, carga="horizontal",
                          material=Elastic1D(E=E, density=7850.0))
        s = NewtonNewmarkSolver(Assembler(dom), t_end=2e-3, dt=1e-3,
                                F_func=lambda t: F, linear_algebra="lu")
        with avisos() as registros:
            s.solve()
        self.assertEqual(_con(registros, AVISO_TANGENTE), [])


# ----------------------------------------------------------------------
# 2. Backend iterativo: el hueco deja de ser silencioso
# ----------------------------------------------------------------------

class TestAvisoDelBackendIterativo(unittest.TestCase):

    def test_lineal_avisa_una_vez_por_analisis(self):
        dom, F = armadura(sin_diagonal=False, carga="vertical")
        solver = LinearSolver(Assembler(dom), linear_algebra="iterative")
        with avisos() as registros:
            U1 = solver.solve(F)
            U2 = solver.solve(2.0 * F)        # factorización cacheada: mismo análisis
        aviso = _con(registros, AVISO_ITERATIVO)
        self.assertEqual(len(aviso), 1, registros)
        self.assertIn("LinearSolver", aviso[0])
        self.assertIn("análisis modal", aviso[0])
        np.testing.assert_allclose(U2, 2.0 * U1, rtol=1e-8, atol=1e-14)

    def test_no_lineales_avisan_una_vez_por_analisis(self):
        for nombre, crear in _no_lineales("iterative").items():
            with self.subTest(solver=nombre):
                dom, F = armadura(sin_diagonal=False, carga="horizontal")
                with avisos() as registros:
                    crear(Assembler(dom), dom).solve(F)
                aviso = _con(registros, AVISO_ITERATIVO)
                self.assertEqual(len(aviso), 1, registros)
                self.assertIn(nombre, aviso[0])

    def test_el_mecanismo_no_activado_pasa_con_el_aviso(self):
        """El hueco declarado: con el iterativo el mecanismo que la carga no
        activa no se detecta, pero ya no pasa en silencio."""
        dom, F = armadura(sin_diagonal=True, carga="horizontal")
        with avisos() as registros:
            U = LinearSolver(Assembler(dom), linear_algebra="iterative").solve(F)
        self.assertTrue(np.all(np.isfinite(U)))
        self.assertEqual(len(_con(registros, AVISO_ITERATIVO)), 1)

    def test_con_directo_no_hay_aviso(self):
        for la in DIRECT:
            with self.subTest(backend=la):
                dom, F = armadura(sin_diagonal=False, carga="vertical")
                with avisos() as registros:
                    LinearSolver(Assembler(dom), linear_algebra=la).solve(F)
                    NonlinearSolver(Assembler(dom), num_steps=2, linear_algebra=la).solve(F)
                self.assertEqual(_con(registros, AVISO_ITERATIVO), [])


# ----------------------------------------------------------------------
# 3. Factorización imposible: pivote exactamente nulo
# ----------------------------------------------------------------------

def _armadura_con_barra_colgante():
    """Armadura bien apoyada con una barra horizontal más que cuelga del
    rodillo: su nudo libre no tiene ninguna rigidez vertical (fila vacía)."""
    dom, F = armadura(sin_diagonal=False, carga="vertical")
    libre = dom.add_node(99, [15.0, 0.0])
    dom.add_element(Truss2D(99, [nodo_de_control(dom), libre], Elastic1D(E=E), A=AREA))
    dom.generate_equation_numbers(verbose=False)
    F2 = np.zeros(dom.total_dofs)
    F2[:F.size] = F
    return dom, F2


class TestFactorizacionImposible(unittest.TestCase):

    def test_lineal(self):
        for la in DIRECT:
            with self.subTest(backend=la):
                dom, F = _armadura_con_barra_colgante()
                with self.assertRaises(IllPosedSystemError) as cm:
                    LinearSolver(Assembler(dom), linear_algebra=la).solve(F)
                self.assertIsNotNone(cm.exception.factorization_error)
                self.assertIn("no pudo factorizar", str(cm.exception))

    def test_no_lineal(self):
        for la in DIRECT:
            with self.subTest(backend=la):
                dom, F = _armadura_con_barra_colgante()
                with self.assertRaises(IllPosedSystemError) as cm:
                    NonlinearSolver(Assembler(dom), linear_algebra=la).solve(F)
                self.assertTrue(cm.exception.initial_tangent)
                self.assertIn("exactamente nulo", str(cm.exception))


# ----------------------------------------------------------------------
# 4. Familias de modelos (revisión adversaria del 2026-09-30)
# ----------------------------------------------------------------------
#
# Las pruebas anteriores usan la armadura del ejemplo 11. Éstas fijan que la
# comprobación de la tangente inicial se comporta igual en las demás
# familias: detecta los mecanismos internos que la carga no activa en
# sólidos 3D, en marcos con una rótula interna y en modelos con restricciones
# lineales (MPC), y no se dispara en modelos bien planteados con plasticidad
# perfecta, daño, secciones por fibras, dominios mixtos, MPC, apoyos oblicuos
# ni zonas rígidas de rigidez 1e6 veces mayor (en unidades SI).

from solidum.math.solvers.corrector import NewtonCorrector  # noqa: E402
from solidum.math.convergence import ConvergenceCriterion  # noqa: E402
from solidum.math.solvers.model_checks import check_initial_tangent  # noqa: E402
from solidum.registry import ElementRegistry, MaterialRegistry  # noqa: E402

_EL = ElementRegistry.get
_MAT = MaterialRegistry.get

# Orden local de los nodos en el cuadrado / cubo unidad de referencia.
_QUAD4 = ((0, 0), (1, 0), (1, 1), (0, 1))
_HEX8 = ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
         (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1))


class _Malla:
    """Nodos por coordenadas (los compartidos se funden) y mallas
    estructuradas de Quad4/Hex8 sobre cajas."""

    def __init__(self):
        self.dom = Domain()
        self._ids = {}

    def nodo(self, x):
        clave = tuple(round(float(v), 9) for v in x)
        if clave not in self._ids:
            nid = len(self.dom.nodes) + 1
            while nid in self.dom.nodes:
                nid += 1
            self._ids[clave] = nid
            self.dom.add_node(nid, [float(v) for v in x])
        return self.dom.nodes[self._ids[clave]]

    def elemento(self, tipo, coords, *args, **kwargs):
        el = _EL(tipo)(len(self.dom.elements) + 1, [self.nodo(c) for c in coords],
                       *args, **kwargs)
        self.dom.add_element(el)
        return el

    def caja(self, tipo, origen, tam, n, material):
        ref = _QUAD4 if len(n) == 2 else _HEX8
        h = np.asarray(tam, float) / np.asarray(n, float)
        for celda in np.ndindex(*n):
            o = np.asarray(origen, float) + np.asarray(celda) * h
            self.elemento(tipo, [o + np.asarray(r) * h for r in ref], material)

    def nodos(self, condicion):
        return [self.dom.nodes[i] for x, i in self._ids.items() if condicion(np.asarray(x))]

    def cerrar(self, cargas=()):
        self.dom.generate_equation_numbers(verbose=False)
        F = np.zeros(self.dom.total_dofs)
        for nodo, dof, valor in cargas:
            F[nodo.dofs[dof]] += valor
        return self.dom, F


def _empotrar(nodos, dofs):
    for nodo in nodos:
        for dof in dofs:
            nodo.fix_dof(dof, 0.0)


# --- mecanismos internos que la carga no activa -----------------------------

def _bloques_hex8_unidos_por_una_arista():
    """Bloque A [0,1]³ empotrado en x = 0; bloque B [1,2]×[0,1]×[1,2] que sólo
    comparte con A la arista x = 1, z = 1: B gira alrededor de ella. La carga
    actúa sobre A."""
    m = _Malla()
    mat = _MAT("Elastic3D")(E=30e9, nu=0.2)
    m.caja("Hex8", (0, 0, 0), (1, 1, 1), (2, 2, 2), mat)
    m.caja("Hex8", (1, 0, 1), (1, 1, 1), (2, 2, 2), mat)
    _empotrar(m.nodos(lambda x: abs(x[0]) < 1e-9), ("ux", "uy", "uz"))
    return m.cerrar([(m.nodo((0.5, 0.5, 1.0)), "uz", -1e5)])


def _viga_con_rotula_interna_por_mpc():
    """Viga biapoyada Frame2DEuler con una rótula en el centro (nudo duplicado
    y MPC en ux y uy): mecanismo vertical. La carga horizontal en el rodillo
    sólo produce axil y no lo activa."""
    m = _Malla()
    mat = Elastic1D(E=E)
    for i in range(4):
        m.elemento("Frame2DEuler", [(0.75 * i, 0.0), (0.75 * (i + 1), 0.0)], mat, A=1e-2, I=1e-4)
    derecha = [m.dom.add_node(100 + i, [3.0 + 0.75 * i, 0.0]) for i in range(5)]
    for i in range(4):
        m.dom.add_element(Frame2DEuler(100 + i, [derecha[i], derecha[i + 1]], mat,
                                       A=1e-2, I=1e-4))
    _empotrar([m.nodo((0.0, 0.0))], ("ux", "uy"))
    derecha[-1].fix_dof("uy", 0.0)
    for dof in ("ux", "uy"):
        m.dom.add_linear_constraint((derecha[0].id, dof), [(m.nodo((3.0, 0.0)).id, dof)], [1.0])
    dom, _ = m.cerrar()
    F = np.zeros(dom.total_dofs)
    F[derecha[-1].dofs["ux"]] = 1e4
    return dom, F


def _placas_unidas_por_una_sola_mpc():
    """Dos placas Quad4; la segunda sólo está atada a la primera por el ux de
    un par de nodos (MPC): puede deslizar en vertical y girar. La carga
    horizontal en el nodo atado no lo activa."""
    m = _Malla()
    mat = _MAT("Elastic2D")(E=30e9, nu=0.2)
    m.caja("Quad4", (0, 0), (1, 1), (2, 2), mat)
    m.caja("Quad4", (1.5, 0), (1, 1), (2, 2), mat)
    _empotrar(m.nodos(lambda x: abs(x[0]) < 1e-9), ("ux", "uy"))
    a, b = m.nodo((1.0, 0.5)), m.nodo((1.5, 0.5))
    m.dom.add_linear_constraint((b.id, "ux"), [(a.id, "ux")], [1.0])
    return m.cerrar([(a, "ux", 1e4)])


# --- modelos bien planteados ------------------------------------------------

def _portico_con_fibras_elastoplasticas():
    """Pórtico empotrado con sección por fibras de acero elastoplástico
    perfecto (H = 0), cargado más allá del primer punto de fluencia (la
    carga lateral de colapso por desplazamiento lateral es 4·Mp/h ≈ 1,1 MN)."""
    from solidum.sections import RectangularPatch
    acero = Elastoplastic1D(E=E, sigma_y=250e6, H=0.0)
    seccion = FiberSection(patches=[RectangularPatch(acero, y=(-0.15, 0.15),
                                                     z=(-0.1, 0.1), ny=10, nz=1)])
    m = _Malla()
    vertices = [(0, 0), (0, 4), (6, 4), (6, 0)]
    for (x0, y0), (x1, y1) in zip(vertices[:-1], vertices[1:]):
        for k in range(4):
            m.elemento("Frame2DEuler", [(x0 + (x1 - x0) * k / 4, y0 + (y1 - y0) * k / 4),
                                        (x0 + (x1 - x0) * (k + 1) / 4, y0 + (y1 - y0) * (k + 1) / 4)],
                       section=seccion)
    _empotrar([m.nodo((0, 0)), m.nodo((6, 0))], ("ux", "uy", "rz"))
    return m.cerrar([(m.nodo((0, 4)), "ux", 8e5), (m.nodo((0, 4)), "uy", -1e5)])


def _placa_von_mises_perfecta():
    m = _Malla()
    m.caja("Quad4", (0, 0), (3, 1), (12, 4),
           _MAT("VonMises2D")(E=30e9, nu=0.3, sigma_y=20e6, H=0.0, hypothesis="plane_strain"))
    _empotrar(m.nodos(lambda x: abs(x[0]) < 1e-9), ("ux", "uy"))
    borde = m.nodos(lambda x: abs(x[0] - 3) < 1e-9)
    return m.cerrar([(n, "uy", -2e5 / len(borde)) for n in borde])


def _bloque_con_dano():
    m = _Malla()
    m.caja("Hex8", (0, 0, 0), (3, 1, 1), (6, 2, 2),
           _MAT("IsotropicDamage3D")(E=30e9, nu=0.2, kappa_0=1e-4, alpha=500.0))
    _empotrar(m.nodos(lambda x: abs(x[0]) < 1e-9), ("ux", "uy", "uz"))
    borde = m.nodos(lambda x: abs(x[0] - 3) < 1e-9)
    return m.cerrar([(n, "uz", -5e5 / len(borde)) for n in borde])


def _portico_con_zonas_rigidas(contraste=1e6):
    """Pórtico empotrado cuyos nudos viga-columna se modelan con zonas
    rígidas de 0,3 m: Frame2DEuler de A·c e I·c, práctica habitual. En
    unidades SI el menor pivote relativo es ~3e-9 con c = 1e6."""
    m = _Malla()
    mat = Elastic1D(E=E)

    def tramo(a, b, c=1.0, n=1):
        a, b = np.asarray(a, float), np.asarray(b, float)
        for k in range(n):
            m.elemento("Frame2DEuler", [a + (b - a) * k / n, a + (b - a) * (k + 1) / n], mat,
                       A=5e-3 * c, I=8e-5 * c)
    tramo((0, 0), (0, 4), n=4)
    tramo((6, 0), (6, 4), n=4)
    tramo((0, 4), (0.3, 4), c=contraste)
    tramo((0.3, 4), (5.7, 4), n=4)
    tramo((5.7, 4), (6, 4), c=contraste)
    _empotrar([m.nodo((0, 0)), m.nodo((6, 0))], ("ux", "uy", "rz"))
    return m.cerrar([(m.nodo((0, 4)), "ux", 5e4), (m.nodo((0, 4)), "uy", -1e5)])


def _muro_con_columna_embebida():
    """Muro Quad4 empotrado en la base con una columna Frame2DEuler que
    comparte con él cuatro nudos (su giro lo restringe el muro) y sobresale
    un metro por encima. La columna aporta el grado de libertad rz."""
    m = _Malla()
    m.caja("Quad4", (0, 0), (1, 2), (4, 8), _MAT("Elastic2D")(E=30e9, nu=0.2))
    _empotrar(m.nodos(lambda x: abs(x[1]) < 1e-9), ("ux", "uy"))
    mat = Elastic1D(E=E)
    for k in range(8):
        m.elemento("Frame2DEuler", [(0.5, 1.0 + 0.25 * k), (0.5, 1.0 + 0.25 * (k + 1))],
                   mat, A=1e-2, I=1e-4)
    return m.cerrar([(m.nodo((0.5, 3.0)), "ux", 1e4)])


def _marco_3d():
    m = _Malla()
    mat = Elastic1D(E=E)
    kw = dict(A=5e-3, Iy=5e-5, Iz=8e-5, J=1e-6, nu=0.3)
    columnas = [(0, 0), (5, 0), (5, 4), (0, 4)]
    for x, y in columnas:
        m.elemento("Frame3D", [(x, y, 0), (x, y, 3)], mat, ref_vector=[1, 0, 0], **kw)
    for (x0, y0), (x1, y1) in zip(columnas, columnas[1:] + columnas[:1]):
        m.elemento("Frame3D", [(x0, y0, 3), (x1, y1, 3)], mat, ref_vector=[0, 0, 1], **kw)
    _empotrar([m.nodo((x, y, 0)) for x, y in columnas], ("ux", "uy", "uz", "rx", "ry", "rz"))
    return m.cerrar([(m.nodo((0, 0, 3)), "ux", 1e4), (m.nodo((0, 0, 3)), "uz", -1e4)])


def _placas_atadas_con_apoyo_oblicuo():
    """Dos placas Quad4 no conformes atadas por MPC en todo el borde común y
    un apoyo deslizante oblicuo (30°) por MPC."""
    m = _Malla()
    mat = _MAT("Elastic2D")(E=30e9, nu=0.2)
    m.caja("Quad4", (0, 0), (1, 1), (4, 4), mat)
    m.caja("Quad4", (1.001, 0), (1, 1), (4, 4), mat)
    _empotrar([m.nodo((0.0, 0.0))], ("ux", "uy"))
    for nodo in m.nodos(lambda x: abs(x[0]) < 1e-9 and x[1] > 0):
        nodo.fix_dof("ux", 0.0)
    for nodo in m.nodos(lambda x: abs(x[0] - 1.0) < 1e-9):
        otro = m.nodo((1.001, nodo.coordinates[1]))
        for dof in ("ux", "uy"):
            m.dom.add_linear_constraint((otro.id, dof), [(nodo.id, dof)], [1.0])
    esquina = m.nodo((2.001, 0.0))
    m.dom.add_linear_constraint((esquina.id, "uy"), [(esquina.id, "ux")],
                                [float(np.tan(np.radians(30.0)))])
    return m.cerrar([(m.nodo((2.001, 1.0)), "uy", -1e4)])


class TestFamiliasDeModelos(unittest.TestCase):
    """La tangente inicial fuera de la armadura del ejemplo 11."""

    def test_mecanismos_internos_no_activados_se_rechazan(self):
        casos = {
            "Hex8, bloque unido por una arista": _bloques_hex8_unidos_por_una_arista,
            "Frame2DEuler, rótula interna por MPC": _viga_con_rotula_interna_por_mpc,
            "Quad4, placa atada por una sola MPC": _placas_unidas_por_una_sola_mpc,
        }
        for nombre, modelo in casos.items():
            for la in DIRECT:
                for solver in ("NonlinearSolver", "ArcLengthSolver"):
                    with self.subTest(modelo=nombre, backend=la, solver=solver):
                        dom, F = modelo()
                        a = Assembler(dom)
                        s = (NonlinearSolver(a, num_steps=2, linear_algebra=la)
                             if solver == "NonlinearSolver"
                             else ArcLengthSolver(a, max_lambda=1.0, linear_algebra=la))
                        with self.assertRaises(IllPosedSystemError) as cm:
                            s.solve(F)
                        self.assertTrue(cm.exception.initial_tangent)
                        self.assertGreaterEqual(cm.exception.zero_pivots, 1)

    def test_bien_planteados_no_avisan(self):
        casos = {
            "pórtico con fibras elastoplásticas perfectas": _portico_con_fibras_elastoplasticas,
            "placa von Mises perfecta": _placa_von_mises_perfecta,
            "bloque Hex8 con daño": _bloque_con_dano,
            "pórtico con zonas rígidas ×1e6 (SI)": _portico_con_zonas_rigidas,
            "muro Quad4 con columna embebida": _muro_con_columna_embebida,
            "marco 3D": _marco_3d,
        }
        for nombre, modelo in casos.items():
            for la in DIRECT:
                with self.subTest(modelo=nombre, backend=la):
                    dom, F = modelo()
                    with avisos() as registros:
                        U = NonlinearSolver(Assembler(dom), num_steps=2,
                                            linear_algebra=la).solve(F)
                    self.assertEqual(_con(registros, AVISO_TANGENTE), [])
                    self.assertTrue(np.all(np.isfinite(U)))

    def test_mpc_y_apoyo_oblicuo_no_disparan_la_comprobacion(self):
        """Restricciones lineales bien planteadas: la tangente inicial
        reducida (``u = T·u_libre + g``) es regular. Se comprueba con las
        piezas públicas —``check_initial_tangent`` sobre la primera
        factorización del corrector— porque es lo que decide el veredicto."""
        for la in DIRECT:
            with self.subTest(backend=la):
                dom, F = _placas_atadas_con_apoyo_oblicuo()
                a = Assembler(dom)
                K, _ = a.assemble_non_linear_system(np.zeros(dom.total_dofs))
                K_red, F_red, _, _ = a.reduce(K, F)
                corrector = NewtonCorrector(ConvergenceCriterion(), max_iter=5,
                                            is_symmetric=True, linear_algebra=la)
                check_initial_tangent(a, corrector, "NonlinearSolver")
                with avisos() as registros:
                    x = corrector.solve(K_red, F_red)
                self.assertEqual(_con(registros, AVISO_TANGENTE), [])
                self.assertTrue(np.all(np.isfinite(x)))


if __name__ == "__main__":
    unittest.main()
