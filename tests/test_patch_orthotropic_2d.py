"""Patch test de MacNeal-Harder con material **ortótropo** sobre los cinco
elementos sólidos 2D del catálogo.

Referencia
----------
MacNeal, R.H.; Harder, R.L. (1985). *A proposed standard set of problems to
test finite element accuracy*. Finite Elements in Analysis and Design, 1(1).

Qué añade sobre los patch tests existentes
------------------------------------------
`tests/test_patch_solid_2d.py` y `tests/test_patch_higher_order.py` ya cubren
el patch test para Quad4, Tri3, Quad8, Quad9 y Tri6 — pero **los cinco usan un
material isótropo** definido dentro del propio test. En un isótropo el bloque
de acoplamiento de la constitutiva es idénticamente nulo:

    C[0,2] = C[1,2] = 0

Ese cero **enmascara** toda una familia de errores. Un elemento puede reproducir
un campo lineal correctamente con material isótropo y fallar con uno anisótropo
si, por ejemplo, ensambla `B^T C B` con una C transpuesta, confunde el orden de
las componentes de Voigt, o pierde el factor 2 del cortante *engineering* al
recuperar esfuerzos — errores todos invisibles cuando `C[0,2] = C[1,2] = 0`.

Con `Orthotropic2D` a `theta` fuera de los ejes principales la constitutiva se
llena por completo (ADR 0013, spec §2.4) y esos caminos quedan ejercitados.

Qué se verifica
---------------
Sobre malla distorsionada con nodos interiores libres, imponiendo en la frontera
el campo lineal ``u = a0 + a1·x + a2·y``, ``v = b0 + b1·x + b2·y``:

1. Los nodos interiores adoptan exactamente el campo lineal.
2. La deformación en **cada punto de Gauss** es constante e igual a la
   analítica.
3. El **esfuerzo** en cada punto de Gauss es uniforme e igual a
   ``C_glob @ eps``, con `C_glob` la constitutiva rotada.

El punto 3 es el que no existía: se lee a través de `compute_gauss_state()`,
la API pública de recuperación de esfuerzos, de modo que se valida el camino
completo material -> elemento -> post-proceso y no sólo la matriz constitutiva.

Un patch test con σ verificado es, además, condición necesaria de convergencia:
si el esfuerzo no es uniforme bajo deformación uniforme, el elemento no
converge por mucho que los desplazamientos parezcan correctos.

Barrido de orientaciones
------------------------
`theta = 0°` es el caso degenerado (sin acoplamiento) y sirve de control:
verifica que el test pasa por la razón correcta y no por una tolerancia laxa.
`theta ∈ {30°, 45°, 90°}` activan el acoplamiento; 90° además intercambia los
papeles de E1 y E2, que es donde se detecta una convención de Poisson invertida.

Datos: bambú, órdenes de literatura general (E1/E2 = 18.75) — los mismos que
usa `tests/validation/test_off_axis_orthotropic.py`.
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from solidum.core.domain import Domain
from solidum.elements.solid_2d import Quad4, Quad8, Quad9, Tri3, Tri6
from solidum.materials.orthotropic_2d import Orthotropic2D
from solidum.math.assembly import Assembler
from solidum.math.solvers import LinearSolver

# Bambú, literatura general. No corresponde a una especie concreta —
# ver §Diálogo de docs/specs/Orthotropic2D.md.
BAMBU = dict(E1=15.0e9, E2=0.8e9, G12=0.7e9, nu12=0.35)

# theta = 0 es el control sin acoplamiento; el resto lo activan.
ANGULOS = (0.0, 30.0, 45.0, 90.0)

# Campo lineal de MacNeal-Harder: genera eps_xx = eps_yy = gamma_xy = 1e-3.
A1, A2 = 1.0e-3, 0.5e-3
B1, B2 = 0.5e-3, 1.0e-3
EXPECTED_STRAIN = np.array([A1, B2, A2 + B1])


def _linear_field(x, y):
    return A1 * x + A2 * y, B1 * x + B2 * y


def _verificar_gauss(caso, domain, U, material):
    """ε y σ constantes en todos los puntos de Gauss de todos los elementos.

    Lee el estado a través de `compute_gauss_state()` — la API pública de
    post-proceso— en vez de recalcular B en el test. Así se valida el camino
    que realmente recorre un usuario al pedir esfuerzos, con la constitutiva
    anisótropa dentro del elemento.

    Tolerancias relativas (ADR 0006): σ es del orden de 1e7 Pa, comparar contra
    cero absoluto convertiría el test en una medida de la magnitud de E1.
    """
    expected_stress = material.C @ EXPECTED_STRAIN
    escala_sig = np.abs(expected_stress).max()
    escala_eps = np.abs(EXPECTED_STRAIN).max()

    for elem in domain.elements.values():
        gs = elem.compute_gauss_state(U)
        n_g = gs['strain'].shape[0]
        caso.assertGreater(n_g, 0, f"elem {elem.id} sin puntos de Gauss")

        for idx in range(n_g):
            np.testing.assert_allclose(
                gs['strain'][idx], EXPECTED_STRAIN,
                rtol=0.0, atol=1e-9 * escala_eps,
                err_msg=f"ε no constante en elem {elem.id}, Gauss {idx}")
            np.testing.assert_allclose(
                gs['stress'][idx], expected_stress,
                rtol=0.0, atol=1e-9 * escala_sig,
                err_msg=f"σ no constante en elem {elem.id}, Gauss {idx}")


# =============================================================================
# Quad4 — malla de cinco quads distorsionados (MacNeal-Harder estándar)
# =============================================================================

class TestPatchOrtotropoQuad4(unittest.TestCase):
    """Patch test estándar de MacNeal-Harder, ahora con fibra fuera de eje."""

    EXT_NODES = {1: (0.00, 0.00), 2: (0.24, 0.00),
                 3: (0.24, 0.12), 4: (0.00, 0.12)}
    INT_NODES = {5: (0.04, 0.02), 6: (0.18, 0.03),
                 7: (0.16, 0.08), 8: (0.08, 0.08)}
    CONNECTIVITY = {1: (1, 2, 6, 5), 2: (2, 3, 7, 6), 3: (3, 4, 8, 7),
                    4: (4, 1, 5, 8), 5: (5, 6, 7, 8)}

    def _montar(self, theta):
        dom = Domain()
        for nid, xy in {**self.EXT_NODES, **self.INT_NODES}.items():
            dom.add_node(nid, list(xy))

        material = Orthotropic2D(**BAMBU, theta=theta)
        for eid, conn in self.CONNECTIVITY.items():
            nodes = [dom.get_node(nid) for nid in conn]
            dom.add_element(Quad4(eid, nodes, material, thickness=0.001))

        for nid, (x, y) in self.EXT_NODES.items():
            u_x, u_y = _linear_field(x, y)
            node = dom.get_node(nid)
            node.fix_dof('ux', u_x)
            node.fix_dof('uy', u_y)

        dom.generate_equation_numbers(verbose=False)
        return dom, material

    def test_nodos_interiores_siguen_el_campo_lineal(self):
        for theta in ANGULOS:
            with self.subTest(theta=theta):
                dom, _ = self._montar(theta)
                U = LinearSolver(Assembler(dom)).solve(np.zeros(dom.total_dofs))

                for nid, (x, y) in self.INT_NODES.items():
                    u_esp, v_esp = _linear_field(x, y)
                    node = dom.get_node(nid)
                    self.assertAlmostEqual(U[node.dofs['ux']], u_esp, places=12,
                                           msg=f"ux nodo {nid}, θ={theta}")
                    self.assertAlmostEqual(U[node.dofs['uy']], v_esp, places=12,
                                           msg=f"uy nodo {nid}, θ={theta}")

    def test_estado_uniforme_en_puntos_de_gauss(self):
        for theta in ANGULOS:
            with self.subTest(theta=theta):
                dom, material = self._montar(theta)
                U = LinearSolver(Assembler(dom)).solve(np.zeros(dom.total_dofs))
                _verificar_gauss(self, dom, U, material)

    def test_el_acoplamiento_esta_realmente_activo(self):
        """Guardia del propio test: a 45° la constitutiva debe estar llena.

        Sin esta comprobación, un error que anulase C[0,2] y C[1,2] dejaría
        pasar el patch test silenciosamente — el test seguiría verde pero ya no
        estaría probando nada anisótropo.
        """
        material = Orthotropic2D(**BAMBU, theta=45.0)
        escala = np.abs(material.C).max()
        self.assertGreater(abs(material.C[0, 2]) / escala, 1e-3,
                           "C[0,2] nulo a 45°: el test no ejercita acoplamiento")
        self.assertGreater(abs(material.C[1, 2]) / escala, 1e-3,
                           "C[1,2] nulo a 45°: el test no ejercita acoplamiento")

        # Y en ejes principales debe seguir siendo cero: el acoplamiento
        # aparece por la rotación, no por construir mal C_mat.
        en_ejes = Orthotropic2D(**BAMBU, theta=0.0)
        escala_0 = np.abs(en_ejes.C).max()
        self.assertLess(abs(en_ejes.C[0, 2]) / escala_0, 1e-14)
        self.assertLess(abs(en_ejes.C[1, 2]) / escala_0, 1e-14)


# =============================================================================
# Tri3 — cuadrado con nodo interior descentrado
# =============================================================================

class TestPatchOrtotropoTri3(unittest.TestCase):
    """CST reproduce campos lineales exactamente; con C llena también debe."""

    EXT_NODES = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (1.0, 1.0), 4: (0.0, 1.0)}
    INT_NODES = {5: (0.37, 0.42)}
    CONNECTIVITY = {1: (1, 2, 5), 2: (2, 3, 5), 3: (3, 4, 5), 4: (4, 1, 5)}

    def _montar(self, theta):
        dom = Domain()
        for nid, xy in {**self.EXT_NODES, **self.INT_NODES}.items():
            dom.add_node(nid, list(xy))

        material = Orthotropic2D(**BAMBU, theta=theta)
        for eid, conn in self.CONNECTIVITY.items():
            nodes = [dom.get_node(nid) for nid in conn]
            dom.add_element(Tri3(eid, nodes, material, thickness=1.0))

        for nid, (x, y) in self.EXT_NODES.items():
            u_x, u_y = _linear_field(x, y)
            node = dom.get_node(nid)
            node.fix_dof('ux', u_x)
            node.fix_dof('uy', u_y)

        dom.generate_equation_numbers(verbose=False)
        return dom, material

    def test_nodos_interiores_siguen_el_campo_lineal(self):
        for theta in ANGULOS:
            with self.subTest(theta=theta):
                dom, _ = self._montar(theta)
                U = LinearSolver(Assembler(dom)).solve(np.zeros(dom.total_dofs))

                for nid, (x, y) in self.INT_NODES.items():
                    u_esp, v_esp = _linear_field(x, y)
                    node = dom.get_node(nid)
                    self.assertAlmostEqual(U[node.dofs['ux']], u_esp, places=12,
                                           msg=f"ux nodo {nid}, θ={theta}")
                    self.assertAlmostEqual(U[node.dofs['uy']], v_esp, places=12,
                                           msg=f"uy nodo {nid}, θ={theta}")

    def test_estado_uniforme_en_puntos_de_gauss(self):
        for theta in ANGULOS:
            with self.subTest(theta=theta):
                dom, material = self._montar(theta)
                U = LinearSolver(Assembler(dom)).solve(np.zeros(dom.total_dofs))
                _verificar_gauss(self, dom, U, material)


# =============================================================================
# Quad8 / Quad9 — malla 2×2 con corner central desplazado (mapeo no afín)
# =============================================================================

def _malla_cuadratica(elem_cls, material, include_center,
                      x_shift=0.15, y_shift=0.10):
    """Malla 2×2 sobre [0,2]² con el corner central desplazado.

    Réplica de la geometría de `tests/test_patch_higher_order.py`: el mapeo
    isoparamétrico deja de ser afín, que es la condición bajo la cual el patch
    test tiene poder de detección real.
    """
    dom = Domain()
    grid = {}
    nid = 0

    def base_x(I):
        return I * 0.5

    def base_y(J):
        return J * 0.5

    for J in range(5):
        for I in range(5):
            es_centro_param = (I % 2 == 1) and (J % 2 == 1)
            if es_centro_param and not include_center:
                continue
            nid += 1
            x, y = base_x(I), base_y(J)
            if I == 2 and J == 2:
                x += x_shift
                y += y_shift
            # Midnodes vecinos al corner desplazado: al midpoint geométrico.
            elif (I, J) == (1, 2):
                x = 0.5 * (base_x(0) + (base_x(2) + x_shift))
                y = 0.5 * (base_y(2) + (base_y(2) + y_shift))
            elif (I, J) == (3, 2):
                x = 0.5 * ((base_x(2) + x_shift) + base_x(4))
                y = 0.5 * ((base_y(2) + y_shift) + base_y(2))
            elif (I, J) == (2, 1):
                x = 0.5 * ((base_x(2) + x_shift) + base_x(2))
                y = 0.5 * ((base_y(2) + y_shift) + base_y(0))
            elif (I, J) == (2, 3):
                x = 0.5 * ((base_x(2) + x_shift) + base_x(2))
                y = 0.5 * ((base_y(2) + y_shift) + base_y(4))
            grid[(I, J)] = dom.add_node(nid, [x, y])

    eid = 0
    for jy in range(2):
        for ix in range(2):
            I0, J0 = 2 * ix, 2 * jy
            conn = [grid[(I0, J0)], grid[(I0 + 2, J0)],
                    grid[(I0 + 2, J0 + 2)], grid[(I0, J0 + 2)],
                    grid[(I0 + 1, J0)], grid[(I0 + 2, J0 + 1)],
                    grid[(I0 + 1, J0 + 2)], grid[(I0, J0 + 1)]]
            if include_center:
                conn.append(grid[(I0 + 1, J0 + 1)])
            eid += 1
            dom.add_element(elem_cls(eid, conn, material, thickness=1.0))

    ext_ids = {(I, J) for (I, J) in grid if I in (0, 4) or J in (0, 4)}
    int_ids = set(grid.keys()) - ext_ids
    return dom, ext_ids, int_ids, grid


class _MixinPatchOrtotropoCuadratico:
    """Verificación común para Quad8 y Quad9 con material ortótropo."""

    ELEM_CLS = None
    INCLUDE_CENTER = False

    def _montar(self, theta):
        material = Orthotropic2D(**BAMBU, theta=theta)
        dom, ext_ids, int_ids, grid = _malla_cuadratica(
            self.ELEM_CLS, material, include_center=self.INCLUDE_CENTER)

        for ij in ext_ids:
            node = grid[ij]
            x, y = node.coordinates
            u_x, u_y = _linear_field(x, y)
            node.fix_dof('ux', u_x)
            node.fix_dof('uy', u_y)

        dom.generate_equation_numbers(verbose=False)
        return dom, material, int_ids, grid

    def test_nodos_interiores_siguen_el_campo_lineal(self):
        for theta in ANGULOS:
            with self.subTest(theta=theta):
                dom, _, int_ids, grid = self._montar(theta)
                U = LinearSolver(Assembler(dom)).solve(np.zeros(dom.total_dofs))

                for ij in int_ids:
                    node = grid[ij]
                    x, y = node.coordinates
                    u_esp, v_esp = _linear_field(x, y)
                    self.assertAlmostEqual(U[node.dofs['ux']], u_esp, places=11,
                                           msg=f"ux nodo {ij}, θ={theta}")
                    self.assertAlmostEqual(U[node.dofs['uy']], v_esp, places=11,
                                           msg=f"uy nodo {ij}, θ={theta}")

    def test_estado_uniforme_en_puntos_de_gauss(self):
        for theta in ANGULOS:
            with self.subTest(theta=theta):
                dom, material, _, _ = self._montar(theta)
                U = LinearSolver(Assembler(dom)).solve(np.zeros(dom.total_dofs))
                _verificar_gauss(self, dom, U, material)


class TestPatchOrtotropoQuad8(_MixinPatchOrtotropoCuadratico, unittest.TestCase):
    ELEM_CLS = Quad8
    INCLUDE_CENTER = False


class TestPatchOrtotropoQuad9(_MixinPatchOrtotropoCuadratico, unittest.TestCase):
    ELEM_CLS = Quad9
    INCLUDE_CENTER = True


# =============================================================================
# Tri6 — triangulación irregular con corner central descentrado
# =============================================================================

class TestPatchOrtotropoTri6(unittest.TestCase):
    """Cuadrado [0,1]² en 4 Tri6 con vértice interior descentrado."""

    def _montar(self, theta):
        material = Orthotropic2D(**BAMBU, theta=theta)
        dom = Domain()

        corners = {'sw': (0.0, 0.0), 'se': (1.0, 0.0), 'ne': (1.0, 1.0),
                   'nw': (0.0, 1.0), 'c': (0.37, 0.42)}
        ids = {}
        nid = 0
        for name, xy in corners.items():
            nid += 1
            ids[name] = dom.add_node(nid, list(xy))

        ext_mid = {'se_m': (0.5, 0.0), 'en_m': (1.0, 0.5),
                   'nw_m': (0.5, 1.0), 'ws_m': (0.0, 0.5)}
        for name, xy in ext_mid.items():
            nid += 1
            ids[name] = dom.add_node(nid, list(xy))

        cx, cy = corners['c']
        int_mid = {
            'sw_c': (0.5 * (0.0 + cx), 0.5 * (0.0 + cy)),
            'se_c': (0.5 * (1.0 + cx), 0.5 * (0.0 + cy)),
            'ne_c': (0.5 * (1.0 + cx), 0.5 * (1.0 + cy)),
            'nw_c': (0.5 * (0.0 + cx), 0.5 * (1.0 + cy)),
        }
        for name, xy in int_mid.items():
            nid += 1
            ids[name] = dom.add_node(nid, list(xy))

        # Numeración Solidum del Tri6: [v0, v1, v2, m01, m12, m20].
        triangulos = [
            ('sw', 'se', 'c', 'se_m', 'se_c', 'sw_c'),
            ('se', 'ne', 'c', 'en_m', 'ne_c', 'se_c'),
            ('ne', 'nw', 'c', 'nw_m', 'nw_c', 'ne_c'),
            ('nw', 'sw', 'c', 'ws_m', 'sw_c', 'nw_c'),
        ]
        eid = 0
        for tri in triangulos:
            eid += 1
            dom.add_element(Tri6(eid, [ids[n] for n in tri],
                                 material, thickness=1.0))

        for name in ('sw', 'se', 'ne', 'nw', 'se_m', 'en_m', 'nw_m', 'ws_m'):
            node = ids[name]
            x, y = node.coordinates
            u_x, u_y = _linear_field(x, y)
            node.fix_dof('ux', u_x)
            node.fix_dof('uy', u_y)

        dom.generate_equation_numbers(verbose=False)
        return dom, material, ids

    INTERIORES = ('c', 'sw_c', 'se_c', 'ne_c', 'nw_c')

    def test_nodos_interiores_siguen_el_campo_lineal(self):
        for theta in ANGULOS:
            with self.subTest(theta=theta):
                dom, _, ids = self._montar(theta)
                U = LinearSolver(Assembler(dom)).solve(np.zeros(dom.total_dofs))

                for name in self.INTERIORES:
                    node = ids[name]
                    x, y = node.coordinates
                    u_esp, v_esp = _linear_field(x, y)
                    self.assertAlmostEqual(U[node.dofs['ux']], u_esp, places=11,
                                           msg=f"ux nodo {name}, θ={theta}")
                    self.assertAlmostEqual(U[node.dofs['uy']], v_esp, places=11,
                                           msg=f"uy nodo {name}, θ={theta}")

    def test_estado_uniforme_en_puntos_de_gauss(self):
        for theta in ANGULOS:
            with self.subTest(theta=theta):
                dom, material, _ = self._montar(theta)
                U = LinearSolver(Assembler(dom)).solve(np.zeros(dom.total_dofs))
                _verificar_gauss(self, dom, U, material)


if __name__ == '__main__':
    unittest.main()
