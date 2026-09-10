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
   `N_INTEGRATION_POINTS`, `PRESERVES_SYMMETRY`, `ACCEPTS_UNILATERAL`. Es lo
   que consumen el registro de DOFs, la validación material/elemento y el
   despachador algebraico (ADR 0003).

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
from solidum.registry import CohesiveMaterialRegistry, ElementRegistry

RHO = 7850.0
MAT_1D = Elastic1D(E=2.0e11, density=RHO)
MAT_2D = Elastic2D(E=2.0e11, nu=0.3, density=RHO)
MAT_3D = Elastic3D(E=2.0e11, nu=0.3, density=RHO)
MAT_CABLE = CableMaterial1D(E=2.0e11, density=RHO)
MAT_TERM_2D = ThermalConduction(k=50.0, c=460.0, density=RHO, dim=2)
MAT_TERM_3D = ThermalConduction(k=50.0, c=460.0, density=RHO, dim=3)
MAT_COHESIVO = CohesiveMaterialRegistry.create(
    'CohesiveDamageIsotropic', sigma_t0=3.0e6, G_f=100.0, K_e=1.0e13,
    softening='exponential')

ESPESOR = 0.01
AREA = 0.01


def _nodos(coords):
    return [Node(i + 1, list(p)) for i, p in enumerate(coords)]


# --- Geometrías de referencia -------------------------------------------
# Elementos rectos y bien formados: el barrido comprueba contrato, no
# robustez ante distorsión (eso es competencia de los patch tests).

_LINEA_2D = ((0.0, 0.0), (1.0, 0.0))
_LINEA_3D = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0))
_TRI = ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
_QUAD = ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
_QUAD8 = ((0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0),
          (1.0, 0.0), (2.0, 1.0), (1.0, 2.0), (0.0, 1.0))
_QUAD9 = _QUAD8 + ((1.0, 1.0),)
_TRI6 = ((0.0, 0.0), (2.0, 0.0), (0.0, 2.0),
         (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))

_HEX_VERTICES = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0),
                 (0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 1.0),
                 (1.0, 1.0, 1.0), (0.0, 1.0, 1.0))
_HEX_ARISTAS = ((0.5, 0.0, 0.0), (1.0, 0.5, 0.0), (0.5, 1.0, 0.0),
                (0.0, 0.5, 0.0), (0.5, 0.0, 1.0), (1.0, 0.5, 1.0),
                (0.5, 1.0, 1.0), (0.0, 0.5, 1.0), (0.0, 0.0, 0.5),
                (1.0, 0.0, 0.5), (1.0, 1.0, 0.5), (0.0, 1.0, 0.5))
# Centros de cara en orden VTK_TRIQUADRATIC_HEXAHEDRON: (-x, +x, -y, +y, -z, +z).
_HEX_CARAS = ((0.0, 0.5, 0.5), (1.0, 0.5, 0.5), (0.5, 0.0, 0.5),
              (0.5, 1.0, 0.5), (0.5, 0.5, 0.0), (0.5, 0.5, 1.0))
_HEX20 = _HEX_VERTICES + _HEX_ARISTAS
_HEX27 = _HEX20 + _HEX_CARAS + ((0.5, 0.5, 0.5),)

_TET = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
_TET10 = _TET + ((0.5, 0.0, 0.0), (0.5, 0.5, 0.0), (0.0, 0.5, 0.0),
                 (0.0, 0.0, 0.5), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5))

# Volumen (o área×espesor, o longitud×área) de cada geometría, para contrastar
# la masa total contra ρ·V sin depender de la implementación del elemento.
_V_LINEA = 1.0 * AREA
_V_TRI = 0.5 * ESPESOR
_V_QUAD = 1.0 * ESPESOR
_V_QUAD8 = 4.0 * ESPESOR
_V_TRI6 = 2.0 * ESPESOR
_V_HEX = 1.0
_V_TET = 1.0 / 6.0

# Fábrica por elemento: (constructor, volumen esperado). El volumen es None
# cuando la comprobación de masa no aplica (elementos térmicos: su matriz de
# capacidad es ρ·c·V, otra magnitud; embedded: sin masa implementada).
FABRICAS = {
    'Quad4':        (lambda c: c(1, _nodos(_QUAD), MAT_2D, thickness=ESPESOR), _V_QUAD),
    'Tri3':         (lambda c: c(1, _nodos(_TRI), MAT_2D, thickness=ESPESOR), _V_TRI),
    'Quad8':        (lambda c: c(1, _nodos(_QUAD8), MAT_2D, thickness=ESPESOR), _V_QUAD8),
    'Quad9':        (lambda c: c(1, _nodos(_QUAD9), MAT_2D, thickness=ESPESOR), _V_QUAD8),
    'Tri6':         (lambda c: c(1, _nodos(_TRI6), MAT_2D, thickness=ESPESOR), _V_TRI6),
    'Hex8':         (lambda c: c(1, _nodos(_HEX_VERTICES), MAT_3D), _V_HEX),
    'Hex20':        (lambda c: c(1, _nodos(_HEX20), MAT_3D), _V_HEX),
    'Hex27':        (lambda c: c(1, _nodos(_HEX27), MAT_3D), _V_HEX),
    'Tet4':         (lambda c: c(1, _nodos(_TET), MAT_3D), _V_TET),
    'Tet10':        (lambda c: c(1, _nodos(_TET10), MAT_3D), _V_TET),
    'Truss2D':      (lambda c: c(1, _nodos(_LINEA_2D), MAT_1D, A=AREA), _V_LINEA),
    'Truss2DCorot': (lambda c: c(1, _nodos(_LINEA_2D), MAT_1D, A=AREA), _V_LINEA),
    'Truss3D':      (lambda c: c(1, _nodos(_LINEA_3D), MAT_1D, A=AREA), _V_LINEA),
    'Truss3DCorot': (lambda c: c(1, _nodos(_LINEA_3D), MAT_1D, A=AREA), _V_LINEA),
    'Cable2DCorot': (lambda c: c(1, _nodos(_LINEA_2D), MAT_CABLE, A=AREA), _V_LINEA),
    'Cable3DCorot': (lambda c: c(1, _nodos(_LINEA_3D), MAT_CABLE, A=AREA), _V_LINEA),
    'Frame2DEuler': (lambda c: c(1, _nodos(_LINEA_2D), MAT_1D, A=AREA, I=1.0e-5), _V_LINEA),
    'Frame2DEulerCorot': (
        lambda c: c(1, _nodos(_LINEA_2D), MAT_1D, A=AREA, I=1.0e-5), _V_LINEA),
    'Frame2DTimoshenko': (
        lambda c: c(1, _nodos(_LINEA_2D), MAT_1D, A=AREA, I=1.0e-5, As=0.008,
                    nu=0.3), _V_LINEA),
    'Frame3D':      (lambda c: c(1, _nodos(_LINEA_3D), MAT_1D, A=AREA, Iy=1.0e-5,
                                 Iz=1.0e-5, J=2.0e-5, nu=0.3), _V_LINEA),
    'Quad4Thermal': (lambda c: c(1, _nodos(_QUAD), MAT_TERM_2D,
                                 thickness=ESPESOR), None),
    'Hex8Thermal':  (lambda c: c(1, _nodos(_HEX_VERTICES), MAT_TERM_3D), None),
    'CST_Embedded2D': (
        lambda c: c(1, _nodos(_TRI), MAT_2D, cohesive_material=MAT_COHESIVO,
                    thickness=ESPESOR), None),
}

# Elementos térmicos: sus DOFs son temperaturas, no desplazamientos. Las
# comprobaciones de sólido rígido y de masa mecánica no les aplican.
TERMICOS = frozenset({'Quad4Thermal', 'Hex8Thermal'})


def _construir(nombre):
    fabrica, _ = FABRICAS[nombre]
    return fabrica(ElementRegistry.get(nombre))


def _volumen(nombre):
    return FABRICAS[nombre][1]


def _rigidez(elemento):
    K = np.asarray(elemento.compute_global_stiffness(), dtype=float)
    return K


class TestRegistroCompleto(unittest.TestCase):
    """El barrido se mantiene sincronizado con el registro por construcción."""

    def test_todo_elemento_registrado_tiene_fabrica(self):
        """Guardián del andamiaje: registrar un elemento obliga a cubrirlo.

        Si este test falla, el arreglo NO es borrarlo: es añadir la fábrica que
        falta. Es lo que impide que el elemento N+1 quede fuera del barrido,
        que es exactamente lo que le pasó a `Orthotropic2D` entre los
        materiales.
        """
        registrados = set(ElementRegistry.names())
        cubiertos = set(FABRICAS)

        sin_cubrir = registrados - cubiertos
        self.assertFalse(
            sin_cubrir,
            f"Elementos registrados sin fábrica: {sorted(sin_cubrir)}. "
            f"Añádelos a FABRICAS para que entren en el barrido.")

        sobrantes = cubiertos - registrados
        self.assertFalse(
            sobrantes,
            f"FABRICAS tiene entradas de elementos que ya no se registran: "
            f"{sorted(sobrantes)}. Elimínalas.")

    def test_el_registro_no_esta_vacio(self):
        """Si el autodiscover se rompiera, el barrido pasaría por vacuidad."""
        self.assertGreaterEqual(len(ElementRegistry.names()), 23)

    def test_todos_construyen(self):
        for nombre in ElementRegistry.names():
            with self.subTest(elemento=nombre):
                elemento = _construir(nombre)
                self.assertEqual(elemento.id, 1)


class TestContratoDeclarativo(unittest.TestCase):
    """Atributos de clase exigidos por `solidum/core/element.py`."""

    def test_hereda_de_element(self):
        for nombre in ElementRegistry.names():
            with self.subTest(elemento=nombre):
                self.assertTrue(issubclass(ElementRegistry.get(nombre), Element))

    def test_dof_names_no_vacio_y_sin_duplicados(self):
        """`DOF_NAMES` gobierna el registro de grados de libertad por nodo.

        Un duplicado produciría dos ecuaciones para el mismo DOF y una matriz
        singular mucho más adelante, sin señal en el punto de origen.
        """
        for nombre in ElementRegistry.names():
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
        for nombre in ElementRegistry.names():
            with self.subTest(elemento=nombre):
                clase = ElementRegistry.get(nombre)
                strain = getattr(clase, 'STRAIN_DIM', None)
                flux = getattr(clase, 'FLUX_DIM', None)

                if nombre in TERMICOS:
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
        for nombre in ElementRegistry.names():
            with self.subTest(elemento=nombre):
                elemento = _construir(nombre)
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
        for nombre in ElementRegistry.names():
            with self.subTest(elemento=nombre):
                clase = ElementRegistry.get(nombre)
                self.assertIsInstance(clase.PRESERVES_SYMMETRY, bool)
                self.assertIsInstance(clase.ACCEPTS_UNILATERAL, bool)

    def test_solo_los_preparados_aceptan_material_unilateral(self):
        """Un material unilateral (cable) en un elemento no preparado degenera
        la matriz global cuando la rigidez colapsa a cero. La base lo valida al
        construir; aquí se comprueba que la bandera dice la verdad.
        """
        for nombre in ElementRegistry.names():
            with self.subTest(elemento=nombre):
                elemento = _construir(nombre)
                material = getattr(elemento, 'material', None)
                if material is None or not getattr(material, 'IS_UNILATERAL',
                                                   False):
                    continue
                self.assertTrue(
                    elemento.ACCEPTS_UNILATERAL,
                    f"{nombre} recibió un material unilateral sin declarar "
                    f"ACCEPTS_UNILATERAL")


class TestRigidezElemental(unittest.TestCase):
    """Propiedades que toda `K` elemental debe cumplir."""

    def test_forma_coherente_con_dofs(self):
        for nombre in ElementRegistry.names():
            with self.subTest(elemento=nombre):
                elemento = _construir(nombre)
                esperado = len(elemento.nodes) * len(elemento.DOF_NAMES)
                K = _rigidez(elemento)
                self.assertEqual(K.shape, (esperado, esperado))

    def test_valores_finitos(self):
        for nombre in ElementRegistry.names():
            with self.subTest(elemento=nombre):
                self.assertTrue(np.all(np.isfinite(_rigidez(_construir(nombre)))))

    def test_simetria_cuando_se_declara(self):
        """`PRESERVES_SYMMETRY = True` es una promesa que el despachador cree.

        La capa algebraica (ADR 0003) la agrega con `material.IS_SYMMETRIC`
        para elegir backend: declararla con una `K` asimétrica llevaría a
        factorizar por Cholesky una matriz que no lo admite.
        """
        for nombre in ElementRegistry.names():
            with self.subTest(elemento=nombre):
                elemento = _construir(nombre)
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
        for nombre in ElementRegistry.names():
            if nombre in TERMICOS:
                continue    # sus DOFs son temperaturas: no hay traslación
            with self.subTest(elemento=nombre):
                elemento = _construir(nombre)
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
        for nombre in ElementRegistry.names():
            with self.subTest(elemento=nombre):
                K = _rigidez(_construir(nombre))
                K_sim = 0.5 * (K + K.T)
                autovalores = np.linalg.eigvalsh(K_sim)
                escala = max(np.abs(autovalores).max(), 1.0)
                self.assertGreater(
                    autovalores.min(), -1e-9 * escala,
                    f"{nombre}: K tiene un autovalor negativo")


class TestMasa(unittest.TestCase):
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
        for nombre in ElementRegistry.names():
            volumen = _volumen(nombre)
            if volumen is None or nombre in TERMICOS:
                continue
            with self.subTest(elemento=nombre):
                elemento = _construir(nombre)
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
        for nombre in ElementRegistry.names():
            volumen = _volumen(nombre)
            if volumen is None or nombre in TERMICOS:
                continue
            with self.subTest(elemento=nombre):
                elemento = _construir(nombre)
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
        for nombre in ElementRegistry.names():
            with self.subTest(elemento=nombre):
                elemento = _construir(nombre)
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
        for nombre in ElementRegistry.names():
            with self.subTest(elemento=nombre):
                elemento = _construir(nombre)
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


class TestValidacionDefensiva(unittest.TestCase):
    """Errores del usuario atrapados al construir, no en runtime."""

    def test_lumping_no_soportado_es_rechazado(self):
        """El mensaje debe nombrar la clase y los valores admitidos.

        `validate_lumping_kwarg` centraliza esto en la base (auditoría H-1.6);
        el barrido comprueba que ninguna subclase se lo salte con una
        implementación propia.
        """
        for nombre in ElementRegistry.names():
            with self.subTest(elemento=nombre):
                elemento = _construir(nombre)
                with self.assertRaises(
                        NotImplementedError,
                        msg=f"{nombre} acepta un lumping inexistente"):
                    elemento.compute_mass_matrix('esquema_inexistente')

    def test_los_valores_de_lumping_del_contrato_funcionan(self):
        """Lo contrario del test anterior: los soportados no deben fallar.

        Un elemento que declare soportar `lumped` y lance al pedirlo dejaría el
        análisis dinámico roto sólo para él, sin señal hasta ejecutarlo.
        """
        for nombre in ElementRegistry.names():
            elemento = _construir(nombre)
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
        for nombre in ElementRegistry.names():
            if nombre in TERMICOS:
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
                # reconstrucción: sin esto el test pasaría por accidente.
                # Se acepta cualquiera de las dos validaciones legítimas: la
                # dimensional de la base, o una whitelist más estricta como la
                # de `CST_Embedded2D`, que sólo admite bulks elásticos porque
                # la discrete approach concentra la disipación en Γ_d (ADR
                # 0010). Ambas atrapan el error al construir, que es lo que
                # este test protege.
                mensaje = str(ctx.exception)
                self.assertTrue(
                    'STRAIN_DIM' in mensaje or type(incompatible).__name__ in mensaje,
                    f"{nombre}: el ValueError no identifica el material "
                    f"incompatible: {mensaje}")

    @staticmethod
    def _reconstruir_con_material(nombre, material):
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

        fabrica, _ = FABRICAS[nombre]
        return fabrica(_ConMaterialCambiado)


if __name__ == '__main__':
    unittest.main()
