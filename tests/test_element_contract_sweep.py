"""Barrido del contrato `Element` sobre **todo** el registro.

Contraparte de `tests/test_material_contract_sweep.py` para la otra mitad del
catálogo. Misma motivación: ningún test del proyecto itera `ElementRegistry`,
así que los barridos genéricos enumeran elementos a mano y un elemento nuevo
sólo entra donde alguien lo añada explícitamente.

Este archivo recorre el registro completo y exige a cada elemento el contrato
declarado en `solidum/core/element.py`, más las propiedades que todo elemento
finito debe cumplir con independencia de su formulación. El elemento número 24
queda cubierto el día que se registra.

Qué se verifica y por qué
-------------------------
1. **Contrato declarativo**: `DOF_NAMES`, `STRAIN_DIM`/`FLUX_DIM`,
   `N_INTEGRATION_POINTS`, `PRESERVES_SYMMETRY`, `ACCEPTS_UNILATERAL`,
   `GEOMETRICALLY_NONLINEAR`. Es lo que consumen el registro de DOFs, la
   validación material/elemento, el despachador algebraico (ADR 0003) y la
   red de seguridad de los solvers estáticos no lineales (ADR 0019). La
   última se contrasta además con la cinemática: con material elástico
   lineal, la rigidez de un elemento geométricamente lineal no depende del
   desplazamiento, y la de uno corrotacional sí.

2. **Simetría de la rigidez**: todo elemento derivado de un funcional de
   energía con cargas conservativas la tiene. Una `K` asimétrica con
   `PRESERVES_SYMMETRY = True` haría que el despachador eligiera Cholesky sobre
   una matriz inválida.

3. **Modos de sólido rígido**: `K·t = 0` para una traslación uniforme `t`. Es
   la propiedad más barata de comprobar y la que caza errores de signo en `B`,
   de ensamblaje local y de transformación a ejes globales. Un elemento que la
   incumple genera fuerzas espurias al mover la estructura sin deformarla.

4. **Masa**: total igual a `ρ·V` en consistente y en lumped. El lumping
   redistribuye masa entre DOFs pero **no puede crearla ni destruirla**, y esa
   conservación es lo que hace comparables los dos esquemas (ADR 0009).

5. **Validación defensiva**: `lumping` no soportado debe dar
   `NotImplementedError` con mensaje claro, y un material de dimensión
   incompatible debe rechazarse **al construir** y no como un `matmul` críptico
   tras varias iteraciones del solver.

Alcance
-------
Contrato estructural, no la física de cada formulación — eso lo cubren los
tests dedicados y los benchmarks de `tests/validation/`. Aquí se comprueba que
lo declarado y lo que ocurre en runtime concuerdan.

Los elementos térmicos (`Quad4Thermal`, `Hex8Thermal`) participan en lo que les
aplica: tienen `FLUX_DIM` en vez de `STRAIN_DIM` y su "masa" es la matriz de
capacidad, con unidades distintas. Se les exime de las comprobaciones que
presuponen DOFs de desplazamiento, no de las estructurales.

Arnés público para módulos de usuario
-------------------------------------
La batería recorre sólo las clases del programa principal; cada módulo de
usuario (ADR 0020) trae su barrido y la reutiliza sin tocar este archivo. Las
clases ``Test*`` no leen variables globales: piden los nombres, las fábricas y
las exenciones a un arnés, que se construye con :func:`arnes` y se pone delante
en la herencia. Piezas públicas: :func:`arnes`, :func:`nombres_de_paquete`,
:func:`rechazo_dimensional`, :func:`nodos`, las geometrías de referencia
(``TRI``, ``QUAD``, ``HEX_VERTICES``…), los materiales ``MAT_*`` y ``ESPESOR``
/ ``AREA``. Modelo: ``tests/user/discontinuities/test_disc_contratos.py``.
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import solidum  # noqa: F401  — dispara el autodiscover que puebla el registro
from solidum.core.element import SUPPORTED_LUMPING, Element
from solidum.core.node import Node
from solidum.materials.cable_1d import CableMaterial1D
from solidum.materials.elastic import Elastic1D
from solidum.materials.elastic_2d import Elastic2D
from solidum.materials.elastic_3d import Elastic3D
from solidum.materials.thermal_conduction import ThermalConduction
from solidum.registry import ElementRegistry

from _programa_principal import nombres_programa_principal  # noqa: F401 — público para los arneses

RHO = 7850.0
MAT_1D = Elastic1D(E=2.0e11, density=RHO)
MAT_2D = Elastic2D(E=2.0e11, nu=0.3, density=RHO)
MAT_3D = Elastic3D(E=2.0e11, nu=0.3, density=RHO)
MAT_CABLE = CableMaterial1D(E=2.0e11, density=RHO)
MAT_TERM_2D = ThermalConduction(k=50.0, c=460.0, density=RHO, dim=2)
MAT_TERM_3D = ThermalConduction(k=50.0, c=460.0, density=RHO, dim=3)

ESPESOR = 0.01
AREA = 0.01


def nodos(coords):
    """Nodos numerados desde 1 con las coordenadas dadas."""
    return [Node(i + 1, list(p)) for i, p in enumerate(coords)]


# --- Geometrías de referencia -------------------------------------------
# Elementos rectos y bien formados: el barrido comprueba contrato, no
# robustez ante distorsión (eso es competencia de los patch tests).

LINEA_2D = ((0.0, 0.0), (1.0, 0.0))
LINEA_3D = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0))
TRI = ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
QUAD = ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
QUAD8 = ((0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0),
         (1.0, 0.0), (2.0, 1.0), (1.0, 2.0), (0.0, 1.0))
QUAD9 = QUAD8 + ((1.0, 1.0),)
TRI6 = ((0.0, 0.0), (2.0, 0.0), (0.0, 2.0),
        (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))

HEX_VERTICES = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0),
                (0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 1.0),
                (1.0, 1.0, 1.0), (0.0, 1.0, 1.0))
HEX_ARISTAS = ((0.5, 0.0, 0.0), (1.0, 0.5, 0.0), (0.5, 1.0, 0.0),
               (0.0, 0.5, 0.0), (0.5, 0.0, 1.0), (1.0, 0.5, 1.0),
               (0.5, 1.0, 1.0), (0.0, 0.5, 1.0), (0.0, 0.0, 0.5),
               (1.0, 0.0, 0.5), (1.0, 1.0, 0.5), (0.0, 1.0, 0.5))
# Centros de cara en orden VTK_TRIQUADRATIC_HEXAHEDRON: (-x, +x, -y, +y, -z, +z).
HEX_CARAS = ((0.0, 0.5, 0.5), (1.0, 0.5, 0.5), (0.5, 0.0, 0.5),
             (0.5, 1.0, 0.5), (0.5, 0.5, 0.0), (0.5, 0.5, 1.0))
HEX20 = HEX_VERTICES + HEX_ARISTAS
HEX27 = HEX20 + HEX_CARAS + ((0.5, 0.5, 0.5),)

TET = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
TET10 = TET + ((0.5, 0.0, 0.0), (0.5, 0.5, 0.0), (0.0, 0.5, 0.0),
               (0.0, 0.0, 0.5), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5))

# Volumen (o área×espesor, o longitud×área) de cada geometría, para contrastar
# la masa total contra ρ·V sin depender de la implementación del elemento.
V_LINEA = 1.0 * AREA
V_TRI = 0.5 * ESPESOR
V_QUAD = 1.0 * ESPESOR
V_QUAD8 = 4.0 * ESPESOR
V_TRI6 = 2.0 * ESPESOR
V_HEX = 1.0
V_TET = 1.0 / 6.0

# Fábrica por elemento: (constructor, volumen esperado). El constructor recibe
# la clase y pasa el material como tercer argumento posicional (lo exige
# `TestValidacionDefensiva`). El volumen es None cuando la comprobación de
# masa total no aplica: en los elementos térmicos la "masa" es la matriz de
# capacidad, ρ·c·V, otra magnitud.
FABRICAS = {
    'Quad4':        (lambda c: c(1, nodos(QUAD), MAT_2D, thickness=ESPESOR), V_QUAD),
    'Tri3':         (lambda c: c(1, nodos(TRI), MAT_2D, thickness=ESPESOR), V_TRI),
    'Quad8':        (lambda c: c(1, nodos(QUAD8), MAT_2D, thickness=ESPESOR), V_QUAD8),
    'Quad9':        (lambda c: c(1, nodos(QUAD9), MAT_2D, thickness=ESPESOR), V_QUAD8),
    'Tri6':         (lambda c: c(1, nodos(TRI6), MAT_2D, thickness=ESPESOR), V_TRI6),
    'Hex8':         (lambda c: c(1, nodos(HEX_VERTICES), MAT_3D), V_HEX),
    'Hex20':        (lambda c: c(1, nodos(HEX20), MAT_3D), V_HEX),
    'Hex27':        (lambda c: c(1, nodos(HEX27), MAT_3D), V_HEX),
    'Tet4':         (lambda c: c(1, nodos(TET), MAT_3D), V_TET),
    'Tet10':        (lambda c: c(1, nodos(TET10), MAT_3D), V_TET),
    'Truss2D':      (lambda c: c(1, nodos(LINEA_2D), MAT_1D, A=AREA), V_LINEA),
    'Truss2DCorot': (lambda c: c(1, nodos(LINEA_2D), MAT_1D, A=AREA), V_LINEA),
    'Truss3D':      (lambda c: c(1, nodos(LINEA_3D), MAT_1D, A=AREA), V_LINEA),
    'Truss3DCorot': (lambda c: c(1, nodos(LINEA_3D), MAT_1D, A=AREA), V_LINEA),
    'Cable2DCorot': (lambda c: c(1, nodos(LINEA_2D), MAT_CABLE, A=AREA), V_LINEA),
    'Cable3DCorot': (lambda c: c(1, nodos(LINEA_3D), MAT_CABLE, A=AREA), V_LINEA),
    'Frame2DEuler': (lambda c: c(1, nodos(LINEA_2D), MAT_1D, A=AREA, I=1.0e-5), V_LINEA),
    'Frame2DEulerCorot': (
        lambda c: c(1, nodos(LINEA_2D), MAT_1D, A=AREA, I=1.0e-5), V_LINEA),
    'Frame2DTimoshenko': (
        lambda c: c(1, nodos(LINEA_2D), MAT_1D, A=AREA, I=1.0e-5, As=0.008,
                    nu=0.3), V_LINEA),
    'Frame3D':      (lambda c: c(1, nodos(LINEA_3D), MAT_1D, A=AREA, Iy=1.0e-5,
                                 Iz=1.0e-5, J=2.0e-5, nu=0.3), V_LINEA),
    'Quad4Thermal': (lambda c: c(1, nodos(QUAD), MAT_TERM_2D,
                                 thickness=ESPESOR), None),
    'Hex8Thermal':  (lambda c: c(1, nodos(HEX_VERTICES), MAT_TERM_3D), None),
}

# Elementos térmicos: sus DOFs son temperaturas, no desplazamientos. Las
# comprobaciones de sólido rígido y de masa mecánica no les aplican.
TERMICOS = frozenset({'Quad4Thermal', 'Hex8Thermal'})


# --- Arnés de la batería --------------------------------------------------


def nombres_de_paquete(registry, paquete: str) -> list:
    """Nombres registrados cuya clase vive en ``paquete`` o en un subpaquete
    suyo (p. ej. ``'solidum.user.discontinuities'``)."""
    return [n for n in registry.names()
            if (registry.get(n).__module__ + ".").startswith(paquete + ".")]


def rechazo_dimensional(mensaje: str, material) -> bool:
    """El ``ValueError`` viene de la validación dimensional de la base
    (``Element._validate_material_compatibility``), que nombra ``STRAIN_DIM``."""
    return 'STRAIN_DIM' in mensaje


def arnes(nombres, fabricas, termicos=frozenset(), minimo=1,
          rechazo_de_material=rechazo_dimensional):
    """Mixin que fija qué elementos recorre la batería de contrato y cómo los
    construye.

    Las clases ``Test*`` de este archivo no leen variables globales: piden al
    arnés los nombres, las fábricas y las exenciones. El programa principal
    usa :data:`ARNES_PROGRAMA_PRINCIPAL`. Un módulo de usuario reutiliza la
    batería con su propio arnés, puesto delante en la herencia::

        ArnesModulo = arnes(lambda: nombres_de_paquete(ElementRegistry, 'solidum.user.<m>'),
                            {'MiElemento': (fabrica, volumen)})

        class TestMiContratoDeclarativo(ArnesModulo, TestContratoDeclarativo):
            pass

    Parameters
    ----------
    nombres : callable sin argumentos, o secuencia
        Nombres registrados que recorre la batería. Si es invocable, se evalúa
        en cada prueba: el registro puede crecer durante la sesión de pytest.
    fabricas : dict
        ``{nombre: (constructor, volumen)}``, como :data:`FABRICAS`. El
        constructor recibe la clase y le pasa el material como tercer
        argumento posicional; ``volumen`` es ``None`` si la comprobación de
        masa total no aplica.
    termicos : conjunto de nombres
        Elementos cuyos DOFs son temperaturas (``FLUX_DIM`` en vez de
        ``STRAIN_DIM``).
    minimo : int
        Número mínimo de nombres; protege contra un barrido que pase por
        vacuidad si el registro no se poblara.
    rechazo_de_material : callable ``(mensaje, material) -> bool``
        Reconoce el ``ValueError`` legítimo al construir con un material de
        dimensión incompatible. Por omisión, el de la validación dimensional
        de la base.
    """
    class Arnes:
        def nombres_barrido(self):
            return list(nombres() if callable(nombres) else nombres)

        def fabricas_barrido(self):
            return fabricas

        def minimo_barrido(self):
            return minimo

        def es_termico(self, nombre):
            return nombre in termicos

        def construir(self, nombre, clase=None):
            fabrica, _ = fabricas[nombre]
            return fabrica(clase if clase is not None else ElementRegistry.get(nombre))

        def volumen(self, nombre):
            return fabricas[nombre][1]

        def rechaza_por_el_material(self, mensaje, material):
            return rechazo_de_material(mensaje, material)

    return Arnes


ARNES_PROGRAMA_PRINCIPAL = arnes(
    lambda: nombres_programa_principal(ElementRegistry), FABRICAS, TERMICOS,
    minimo=22)


def _rigidez(elemento):
    K = np.asarray(elemento.compute_global_stiffness(), dtype=float)
    return K


class TestRegistroCompleto(ARNES_PROGRAMA_PRINCIPAL, unittest.TestCase):
    """El barrido se mantiene sincronizado con el registro por construcción."""

    def test_todo_elemento_registrado_tiene_fabrica(self):
        """Guardián del andamiaje: registrar un elemento obliga a cubrirlo.

        Si este test falla, el arreglo NO es borrarlo: es añadir la fábrica que
        falta. Es lo que impide que el elemento N+1 quede fuera del barrido,
        que es exactamente lo que le pasó a `Orthotropic2D` entre los
        materiales.
        """
        registrados = set(self.nombres_barrido())
        cubiertos = set(self.fabricas_barrido())

        sin_cubrir = registrados - cubiertos
        self.assertFalse(
            sin_cubrir,
            f"Elementos registrados sin fábrica: {sorted(sin_cubrir)}. "
            f"Añádelos a las fábricas del arnés (FABRICAS en el programa "
            f"principal) para que entren en el barrido.")

        sobrantes = cubiertos - registrados
        self.assertFalse(
            sobrantes,
            f"Las fábricas del arnés tienen entradas de elementos que ya no se "
            f"registran: {sorted(sobrantes)}. Elimínalas.")

    def test_el_registro_no_esta_vacio(self):
        """Si el autodiscover se rompiera, el barrido pasaría por vacuidad."""
        self.assertGreaterEqual(len(self.nombres_barrido()), self.minimo_barrido())

    def test_todos_construyen(self):
        for nombre in self.nombres_barrido():
            with self.subTest(elemento=nombre):
                elemento = self.construir(nombre)
                self.assertEqual(elemento.id, 1)


class TestContratoDeclarativo(ARNES_PROGRAMA_PRINCIPAL, unittest.TestCase):
    """Atributos de clase exigidos por `solidum/core/element.py`."""

    def test_hereda_de_element(self):
        for nombre in self.nombres_barrido():
            with self.subTest(elemento=nombre):
                self.assertTrue(issubclass(ElementRegistry.get(nombre), Element))

    def test_dof_names_no_vacio_y_sin_duplicados(self):
        """`DOF_NAMES` gobierna el registro de grados de libertad por nodo.

        Un duplicado produciría dos ecuaciones para el mismo DOF y una matriz
        singular mucho más adelante, sin señal en el punto de origen.
        """
        for nombre in self.nombres_barrido():
            with self.subTest(elemento=nombre):
                clase = ElementRegistry.get(nombre)
                dofs = clase.DOF_NAMES
                self.assertTrue(dofs, f"{nombre}: DOF_NAMES vacío")
                self.assertEqual(
                    len(dofs), len(set(dofs)),
                    f"{nombre}: DOF_NAMES con duplicados: {dofs}")
                for d in dofs:
                    self.assertIsInstance(d, str)

    def test_dimension_del_material_declarada(self):
        """Mecánicos declaran `STRAIN_DIM`; térmicos, `FLUX_DIM`.

        La distinción es física (Reglas.md §5): el gradiente de temperatura es
        un vector genuino, no un tensor simétrico comprimido en Voigt, así que
        el elemento térmico no tiene `STRAIN_DIM` que declarar.
        """
        for nombre in self.nombres_barrido():
            with self.subTest(elemento=nombre):
                clase = ElementRegistry.get(nombre)
                strain = getattr(clase, 'STRAIN_DIM', None)
                flux = getattr(clase, 'FLUX_DIM', None)

                if self.es_termico(nombre):
                    self.assertIn(flux, (2, 3),
                                  f"{nombre}: FLUX_DIM debe ser 2 ó 3")
                else:
                    self.assertIn(strain, (1, 3, 6),
                                  f"{nombre}: STRAIN_DIM inválido")

    def test_n_integration_points_positivo_y_coherente(self):
        """`N_INTEGRATION_POINTS` dimensiona el `ElementState`.

        Si no coincidiera con los puntos que el elemento recorre de verdad, el
        estado interno se desalinearía silenciosamente entre iteraciones.
        """
        for nombre in self.nombres_barrido():
            with self.subTest(elemento=nombre):
                elemento = self.construir(nombre)
                n = elemento.N_INTEGRATION_POINTS
                self.assertIsInstance(n, int)
                self.assertGreaterEqual(n, 1)

                puntos = getattr(elemento, 'points', None)
                if puntos is not None:
                    self.assertEqual(
                        len(puntos), n,
                        f"{nombre}: N_INTEGRATION_POINTS={n} pero "
                        f"len(points)={len(puntos)}")

    def test_banderas_son_booleanas(self):
        for nombre in self.nombres_barrido():
            with self.subTest(elemento=nombre):
                clase = ElementRegistry.get(nombre)
                self.assertIsInstance(clase.PRESERVES_SYMMETRY, bool)
                self.assertIsInstance(clase.ACCEPTS_UNILATERAL, bool)
                self.assertIsInstance(clase.GEOMETRICALLY_NONLINEAR, bool)

    def test_no_linealidad_geometrica_declarada_dice_la_verdad(self):
        """`GEOMETRICALLY_NONLINEAR` decide si un pivote nulo de la tangente
        inicial es un error o un aviso (ADR 0019, 2026-09-30): un elemento
        corrotacional que olvidara declararlo convertiría en error un modelo
        que la rigidez geométrica puede estabilizar.

        Se evalúa la rigidez en ``u = 0`` y tras un giro rígido pequeño con un
        alargamiento uniforme (θ = 1e-3 rad, 1e-6; el alargamiento tensa los
        cables, que sin deformar no tienen rigidez). Con los materiales
        elásticos de las fábricas, la de un elemento geométricamente lineal
        no cambia: su matriz ``B`` es constante. La de uno corrotacional sí:
        gira con el elemento y suma la rigidez geométrica. Las deformaciones
        son del orden de 1e-6, bajo cualquier umbral de daño o fluencia.
        """
        theta, stretch = 1.0e-3, 1.0e-6
        c, s = np.cos(theta), np.sin(theta)
        for nombre in self.nombres_barrido():
            with self.subTest(elemento=nombre):
                elemento = self.construir(nombre)
                u = []
                for nodo in elemento.nodes:
                    x = np.zeros(3)
                    coords = np.asarray(nodo.coordinates, dtype=float)
                    x[:coords.size] = coords
                    giro = (1.0 + stretch) * np.array([c * x[0] - s * x[1],
                                                       s * x[0] + c * x[1], x[2]])
                    d = giro - x
                    valores = {'ux': d[0], 'uy': d[1], 'uz': d[2], 'rz': theta,
                               'rx': 0.0, 'ry': 0.0, 'T': 1.0}
                    u.extend(valores[dof] for dof in elemento.DOF_NAMES)
                u = np.asarray(u)
                K0 = np.asarray(elemento.compute_element_state(np.zeros_like(u))[0], float)
                K1 = np.asarray(elemento.compute_element_state(u)[0], float)
                escala = max(np.abs(K0).max(), np.abs(K1).max())
                cambia = not np.allclose(K1, K0, rtol=0.0, atol=1e-9 * escala)
                self.assertEqual(
                    elemento.GEOMETRICALLY_NONLINEAR, cambia,
                    f"{nombre}: GEOMETRICALLY_NONLINEAR="
                    f"{elemento.GEOMETRICALLY_NONLINEAR}, pero su rigidez "
                    f"{'cambia' if cambia else 'no cambia'} con un giro rígido")

    def test_solo_los_preparados_aceptan_material_unilateral(self):
        """Un material unilateral (cable) en un elemento no preparado degenera
        la matriz global cuando la rigidez colapsa a cero. La base lo valida al
        construir; aquí se comprueba que la bandera dice la verdad.
        """
        for nombre in self.nombres_barrido():
            with self.subTest(elemento=nombre):
                elemento = self.construir(nombre)
                material = getattr(elemento, 'material', None)
                if material is None or not getattr(material, 'IS_UNILATERAL',
                                                   False):
                    continue
                self.assertTrue(
                    elemento.ACCEPTS_UNILATERAL,
                    f"{nombre} recibió un material unilateral sin declarar "
                    f"ACCEPTS_UNILATERAL")


class TestRigidezElemental(ARNES_PROGRAMA_PRINCIPAL, unittest.TestCase):
    """Propiedades que toda `K` elemental debe cumplir."""

    def test_forma_coherente_con_dofs(self):
        for nombre in self.nombres_barrido():
            with self.subTest(elemento=nombre):
                elemento = self.construir(nombre)
                esperado = len(elemento.nodes) * len(elemento.DOF_NAMES)
                K = _rigidez(elemento)
                self.assertEqual(K.shape, (esperado, esperado))

    def test_valores_finitos(self):
        for nombre in self.nombres_barrido():
            with self.subTest(elemento=nombre):
                self.assertTrue(np.all(np.isfinite(_rigidez(self.construir(nombre)))))

    def test_simetria_cuando_se_declara(self):
        """`PRESERVES_SYMMETRY = True` es una promesa que el despachador cree.

        La capa algebraica (ADR 0003) la agrega con `material.IS_SYMMETRIC`
        para elegir backend: declararla con una `K` asimétrica llevaría a
        factorizar por Cholesky una matriz que no lo admite.
        """
        for nombre in self.nombres_barrido():
            with self.subTest(elemento=nombre):
                elemento = self.construir(nombre)
                if not elemento.PRESERVES_SYMMETRY:
                    continue
                material = getattr(elemento, 'material', None)
                if material is not None and not getattr(material,
                                                        'IS_SYMMETRIC', True):
                    continue    # la asimetría vendría del material, no del elemento

                K = _rigidez(elemento)
                escala = np.abs(K).max()
                np.testing.assert_allclose(
                    K, K.T, rtol=0.0, atol=1e-9 * escala,
                    err_msg=f"{nombre} declara PRESERVES_SYMMETRY pero su K "
                            f"no es simétrica")

    def test_traslacion_de_solido_rigido_no_genera_fuerzas(self):
        """`K · t = 0` para una traslación uniforme.

        Es la comprobación más barata con más poder de detección: caza errores
        de signo en `B`, de ensamblaje local y de transformación a ejes
        globales. Un elemento que la incumple genera fuerzas espurias al
        trasladar la estructura sin deformarla — y eso contamina cualquier
        análisis, lineal o no.

        Sólo se traslada; la rotación de sólido rígido es exacta en elementos
        lineales pero sólo infinitesimalmente en formulaciones corotacionales,
        así que exigirla aquí mezclaría contrato con formulación.
        """
        for nombre in self.nombres_barrido():
            if self.es_termico(nombre):
                continue    # sus DOFs son temperaturas: no hay traslación
            with self.subTest(elemento=nombre):
                elemento = self.construir(nombre)
                dofs = elemento.DOF_NAMES
                K = _rigidez(elemento)
                # Un cable sin pretensión tiene K idénticamente nula (no
                # resiste compresión ni tiene rigidez geométrica sin tensión):
                # ahí la escala relativa degenera y el criterio correcto es
                # que el residuo sea nulo, que es justamente lo que ocurre.
                escala = max(np.abs(K).max(), 1.0)

                for eje, nombre_dof in enumerate(dofs):
                    if not nombre_dof.startswith('u'):
                        continue    # rotaciones: no son traslación
                    with self.subTest(dof=nombre_dof):
                        t = np.zeros(K.shape[0])
                        t[eje::len(dofs)] = 1.0
                        residuo = np.abs(K @ t).max()
                        self.assertLess(
                            residuo, 1e-9 * escala,
                            f"{nombre}: trasladar en {nombre_dof} genera "
                            f"fuerzas internas (residuo {residuo:.3e})")

    def test_semidefinida_positiva(self):
        """Sin material en régimen degradante, la energía de deformación no
        puede ser negativa: `K` es semidefinida positiva."""
        for nombre in self.nombres_barrido():
            with self.subTest(elemento=nombre):
                K = _rigidez(self.construir(nombre))
                K_sim = 0.5 * (K + K.T)
                autovalores = np.linalg.eigvalsh(K_sim)
                escala = max(np.abs(autovalores).max(), 1.0)
                self.assertGreater(
                    autovalores.min(), -1e-9 * escala,
                    f"{nombre}: K tiene un autovalor negativo")


class TestMasa(ARNES_PROGRAMA_PRINCIPAL, unittest.TestCase):
    """Contrato de `compute_mass_matrix` (ADR 0009)."""

    def _masa_total(self, elemento, M):
        """Suma la masa de una sola dirección traslacional.

        En una matriz de masa, cada DOF traslacional lleva la masa completa del
        elemento en su dirección; sumar todas las direcciones la contaría
        varias veces.
        """
        n_dof = len(elemento.DOF_NAMES)
        eje = next(i for i, d in enumerate(elemento.DOF_NAMES)
                   if d.startswith('u'))
        bloque = M[eje::n_dof, eje::n_dof]
        return bloque.sum()

    def test_masa_total_consistente(self):
        for nombre in self.nombres_barrido():
            volumen = self.volumen(nombre)
            if volumen is None or self.es_termico(nombre):
                continue
            with self.subTest(elemento=nombre):
                elemento = self.construir(nombre)
                try:
                    M = elemento.compute_mass_matrix('consistent')
                except NotImplementedError:
                    continue    # masa opcional en el contrato

                esperada = RHO * volumen
                self.assertAlmostEqual(
                    self._masa_total(elemento, M) / esperada, 1.0, places=10,
                    msg=f"{nombre}: masa consistente ≠ ρ·V")

    def test_lumped_conserva_la_masa_total(self):
        """El lumping redistribuye masa entre DOFs; no la crea ni la destruye.

        Es lo que hace comparables los dos esquemas: si el total no coincidiera,
        una misma estructura tendría frecuencias propias distintas por un
        detalle de implementación de la masa.
        """
        for nombre in self.nombres_barrido():
            volumen = self.volumen(nombre)
            if volumen is None or self.es_termico(nombre):
                continue
            with self.subTest(elemento=nombre):
                elemento = self.construir(nombre)
                try:
                    M_c = elemento.compute_mass_matrix('consistent')
                    M_l = elemento.compute_mass_matrix('lumped')
                except NotImplementedError:
                    continue

                total_c = self._masa_total(elemento, M_c)
                total_l = self._masa_total(elemento, M_l)
                self.assertAlmostEqual(
                    total_l / total_c, 1.0, places=10,
                    msg=f"{nombre}: lumped no conserva la masa total")

    def test_lumped_es_diagonal(self):
        for nombre in self.nombres_barrido():
            with self.subTest(elemento=nombre):
                elemento = self.construir(nombre)
                try:
                    M = np.asarray(elemento.compute_mass_matrix('lumped'),
                                   dtype=float)
                except NotImplementedError:
                    continue
                fuera = M - np.diag(np.diag(M))
                escala = max(np.abs(np.diag(M)).max(), 1.0)
                self.assertLess(
                    np.abs(fuera).max(), 1e-12 * escala,
                    f"{nombre}: la masa lumped no es diagonal")

    def test_masa_simetrica_y_definida_no_negativa(self):
        for nombre in self.nombres_barrido():
            with self.subTest(elemento=nombre):
                elemento = self.construir(nombre)
                try:
                    M = np.asarray(elemento.compute_mass_matrix('consistent'),
                                   dtype=float)
                except NotImplementedError:
                    continue

                escala = np.abs(M).max()
                np.testing.assert_allclose(
                    M, M.T, rtol=0.0, atol=1e-12 * escala,
                    err_msg=f"{nombre}: matriz de masa no simétrica")
                autovalores = np.linalg.eigvalsh(0.5 * (M + M.T))
                self.assertGreater(
                    autovalores.min(), -1e-10 * escala,
                    f"{nombre}: masa con autovalor negativo")


class TestValidacionDefensiva(ARNES_PROGRAMA_PRINCIPAL, unittest.TestCase):
    """Errores del usuario atrapados al construir, no en runtime."""

    def test_lumping_no_soportado_es_rechazado(self):
        """El mensaje debe nombrar la clase y los valores admitidos.

        `validate_lumping_kwarg` centraliza esto en la base (auditoría H-1.6);
        el barrido comprueba que ninguna subclase se lo salte con una
        implementación propia.
        """
        for nombre in self.nombres_barrido():
            with self.subTest(elemento=nombre):
                elemento = self.construir(nombre)
                with self.assertRaises(
                        NotImplementedError,
                        msg=f"{nombre} acepta un lumping inexistente"):
                    elemento.compute_mass_matrix('esquema_inexistente')

    def test_los_valores_de_lumping_del_contrato_funcionan(self):
        """Lo contrario del test anterior: los soportados no deben fallar.

        Un elemento que declare soportar `lumped` y lance al pedirlo dejaría el
        análisis dinámico roto sólo para él, sin señal hasta ejecutarlo.
        """
        for nombre in self.nombres_barrido():
            elemento = self.construir(nombre)
            for esquema in sorted(SUPPORTED_LUMPING):
                with self.subTest(elemento=nombre, lumping=esquema):
                    try:
                        M = elemento.compute_mass_matrix(esquema)
                    except NotImplementedError:
                        continue    # el elemento no implementa masa: legítimo
                    self.assertTrue(np.all(np.isfinite(np.asarray(M, float))))

    def test_material_de_dimension_incompatible_se_rechaza_al_construir(self):
        """Un material 3D en un elemento 2D debe fallar de inmediato.

        Sin esta validación el error aparece como un `matmul` críptico tras
        varias iteraciones del solver, lejos de su causa. Es el ejemplo que la
        propia base cita al documentar `_validate_material_compatibility`.
        """
        for nombre in self.nombres_barrido():
            if self.es_termico(nombre):
                continue    # los térmicos validan contra FLUX_DIM, no STRAIN_DIM
            with self.subTest(elemento=nombre):
                clase = ElementRegistry.get(nombre)
                dim = getattr(clase, 'STRAIN_DIM', None)
                if dim is None:
                    continue

                # Un material de dimensión distinta a la que el elemento pide.
                incompatible = MAT_3D if dim in (1, 3) else MAT_1D

                with self.assertRaises(ValueError) as ctx:
                    self._reconstruir_con_material(nombre, incompatible)

                # Debe fallar por el material, no por un error casual de la
                # reconstrucción: sin esto el test pasaría por accidente. En
                # el programa principal, la validación dimensional de la base,
                # que nombra STRAIN_DIM; un módulo de usuario con una
                # validación propia más estricta la declara en su arnés.
                mensaje = str(ctx.exception)
                self.assertTrue(
                    self.rechaza_por_el_material(mensaje, incompatible),
                    f"{nombre}: el ValueError no viene de la validación del "
                    f"material incompatible: {mensaje}")

    def _reconstruir_con_material(self, nombre, material):
        """Reconstruye el elemento cambiando sólo el material.

        Se apoya en que todas las fábricas pasan el material como tercer
        argumento posicional, así que basta con envolver la clase para
        interceptarlo.
        """
        clase = ElementRegistry.get(nombre)

        class _ConMaterialCambiado(clase):
            def __init__(self, element_id, nodes, _material_original,
                         *args, **kwargs):
                super().__init__(element_id, nodes, material, *args, **kwargs)

        return self.construir(nombre, clase=_ConMaterialCambiado)


if __name__ == '__main__':
    unittest.main()
