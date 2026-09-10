"""Barrido del contrato `Material` sobre **todo** el registro.

Por qué existe
--------------
Los tests genéricos del proyecto enumeran materiales **a mano**: listas de
imports hardcoded (`test_materials_unit.py`, `test_material_deep_coverage.py`,
`test_defensive_validations.py`) o materiales *ad hoc* definidos dentro del
propio test (`test_patch_solid_2d.py`, `test_rigid_body_modes.py`). Ninguno
itera `MaterialRegistry`.

La consecuencia es que **un material nuevo no entra solo en ningún barrido**.
Entra sólo si alguien recuerda añadirlo a cada lista, y con 13 materiales
registrados eso ya dejó de ocurrir: el docstring de `test_materials_unit.py`
todavía habla de "los seis materiales registrados". `Orthotropic2D` estuvo
ausente de doce barridos genéricos hasta que se auditó explícitamente.

Este archivo invierte la carga de la prueba. Recorre `MaterialRegistry.names()`
y exige a cada material cumplir el contrato de `solidum/core/material.py`. El
material número 14 queda cubierto el día que se registra, sin que nadie tenga
que acordarse — que es exactamente el criterio de `Reglas.md §1`: cada decisión
se evalúa contra el coste de añadir el componente N+1.

El mecanismo autodefensivo
--------------------------
`MUESTRAS` mapea cada material a un juego de parámetros válidos, necesario
porque las firmas de constructor son legítimamente heterogéneas (un `Elastic1D`
pide `E`; un `DruckerPrager2D` pide cuatro parámetros). Esa tabla es el único
punto de mantenimiento manual, y `test_toda_clase_registrada_tiene_muestra`
falla en cuanto se registra un material que no figura en ella. Así el olvido es
imposible: o se añade la entrada, o la suite se pone roja.

Alcance
-------
Se verifica el contrato **estructural y declarativo**, no la física de cada
material — eso es competencia de sus tests dedicados y de los benchmarks de
`tests/validation/`. Aquí se comprueba que lo declarado y lo que ocurre en
runtime concuerdan, que es donde se cuelan las incoherencias silenciosas.

Registros paralelos
-------------------
`CohesiveMaterialRegistry` y `ThermalMaterialRegistry` son deliberadamente
distintos (ADR 0010 y Etapa 8): un cohesivo relaciona tracción con salto de
desplazamiento y un térmico flujo con gradiente, ninguno σ con ε en Voigt. No
comparten este contrato y por eso no se barren aquí; sí se comprueba que no
haya materiales traspapelados entre registros.
"""
import inspect
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import solidum  # noqa: F401  — dispara el autodiscover que puebla el registro
from solidum.core.material import Material
from solidum.registry import (CohesiveMaterialRegistry, MaterialRegistry,
                              ThermalMaterialRegistry)

# Parámetros válidos y físicamente razonables para construir cada material.
# Sólo los obligatorios: los opcionales se dejan en su default a propósito,
# porque el contrato debe cumplirse también en la construcción mínima.
MUESTRAS = {
    'CableMaterial1D':   dict(E=2.0e11),
    'DruckerPrager2D':   dict(E=2.0e10, nu=0.2, cohesion=1.0e6, phi_deg=30.0),
    'DruckerPrager3D':   dict(E=2.0e10, nu=0.2, cohesion=1.0e6, phi_deg=30.0),
    'Elastic1D':         dict(E=2.0e11),
    'Elastic2D':         dict(E=2.0e11, nu=0.3),
    'Elastic3D':         dict(E=2.0e11, nu=0.3),
    'Elastoplastic1D':   dict(E=2.0e11, sigma_y=2.5e8),
    'IsotropicDamage1D': dict(E=2.0e10, kappa_0=1.0e-4, alpha=0.99),
    'IsotropicDamage2D': dict(E=2.0e10, nu=0.2, kappa_0=1.0e-4, alpha=0.99),
    'IsotropicDamage3D': dict(E=2.0e10, nu=0.2, kappa_0=1.0e-4, alpha=0.99),
    'Orthotropic2D':     dict(E1=15.0e9, E2=0.8e9, G12=0.7e9, nu12=0.35),
    'VonMises2D':        dict(E=2.0e11, nu=0.3, sigma_y=2.5e8),
    'VonMises3D':        dict(E=2.0e11, nu=0.3, sigma_y=2.5e8),
}

# Dimensiones de deformación admitidas por el contrato (docstring de Material).
STRAIN_DIMS_VALIDAS = (1, 3, 6)


def _construir(nombre):
    return MaterialRegistry.create(nombre, **MUESTRAS[nombre])


def _deformacion_pequena(dim, escala=1.0e-8):
    """Deformación en régimen indudablemente elástico.

    Se mantiene muy por debajo de cualquier umbral de fluencia o daño de las
    muestras, de modo que las comprobaciones estructurales no dependan del
    régimen constitutivo en que caiga cada material.
    """
    eps = np.zeros(dim)
    eps[0] = escala
    return eps


class TestRegistroCompleto(unittest.TestCase):
    """El barrido se mantiene sincronizado con el registro por construcción."""

    def test_toda_clase_registrada_tiene_muestra(self):
        """Guardián del andamiaje: registrar un material obliga a cubrirlo.

        Si este test falla, el arreglo NO es borrarlo: es añadir la entrada
        que falta en `MUESTRAS`. Ese es todo el mecanismo que impide que el
        material N+1 vuelva a quedarse fuera de los barridos genéricos.
        """
        registrados = set(MaterialRegistry.names())
        cubiertos = set(MUESTRAS)

        sin_cubrir = registrados - cubiertos
        self.assertFalse(
            sin_cubrir,
            f"Materiales registrados sin entrada en MUESTRAS: "
            f"{sorted(sin_cubrir)}. Añádelos a la tabla para que entren en el "
            f"barrido del contrato.")

        sobrantes = cubiertos - registrados
        self.assertFalse(
            sobrantes,
            f"MUESTRAS tiene entradas de materiales que ya no se registran: "
            f"{sorted(sobrantes)}. Elimínalas.")

    def test_el_registro_no_esta_vacio(self):
        """Si el autodiscover se rompiera, el barrido pasaría por vacuidad."""
        self.assertGreaterEqual(len(MaterialRegistry.names()), 13)

    def test_los_registros_paralelos_no_se_mezclan(self):
        """Cohesivos y térmicos viven en registros propios (ADR 0010, Etapa 8).

        Un material traspapelado al registro principal llegaría al parser YAML
        como si fuera un material de bulk y fallaría de forma críptica al
        construir un elemento sólido.
        """
        principal = set(MaterialRegistry.names())
        for otro_registro, etiqueta in ((CohesiveMaterialRegistry, 'cohesivo'),
                                        (ThermalMaterialRegistry, 'térmico')):
            solapamiento = principal & set(otro_registro.names())
            self.assertFalse(
                solapamiento,
                f"Materiales en el registro principal y en el {etiqueta} a la "
                f"vez: {sorted(solapamiento)}")


class TestContratoDeclarativo(unittest.TestCase):
    """Atributos de clase exigidos por `solidum/core/material.py`."""

    def test_hereda_de_material(self):
        for nombre in MaterialRegistry.names():
            with self.subTest(material=nombre):
                clase = MaterialRegistry.get(nombre)
                self.assertTrue(
                    issubclass(clase, Material),
                    f"{nombre} no hereda de Material: los elementos validan "
                    f"STRAIN_DIM asumiendo el contrato de la base.")

    def test_strain_dim_declarado_y_valido(self):
        """`STRAIN_DIM` es lo que los elementos validan al construirse.

        No tiene default en la base a propósito: declararlo es obligatorio, y
        omitirlo debe romper aquí y no en un `matmul` críptico en runtime.
        """
        for nombre in MaterialRegistry.names():
            with self.subTest(material=nombre):
                clase = MaterialRegistry.get(nombre)
                self.assertTrue(
                    hasattr(clase, 'STRAIN_DIM'),
                    f"{nombre} no declara STRAIN_DIM")
                self.assertIn(clase.STRAIN_DIM, STRAIN_DIMS_VALIDAS)

    def test_primary_state_var_coherente_con_el_estado(self):
        """`PRIMARY_STATE_VAR` nombra una clave real del dict de estado.

        Lo consume `VtkExporter` para exportar la variable interna sin
        hardcodear nombres: un nombre que no exista produce una exportación
        vacía sin error visible.
        """
        for nombre in MaterialRegistry.names():
            with self.subTest(material=nombre):
                material = _construir(nombre)
                clave = material.PRIMARY_STATE_VAR

                if clave is None:
                    continue    # material sin variables internas

                self.assertIsInstance(clave, str)
                _, _, estado = material.compute_state(
                    _deformacion_pequena(material.STRAIN_DIM), None)
                self.assertIsInstance(
                    estado, dict,
                    f"{nombre} declara PRIMARY_STATE_VAR='{clave}' pero no "
                    f"devuelve un dict de estado")
                self.assertIn(
                    clave, estado,
                    f"{nombre} declara PRIMARY_STATE_VAR='{clave}', ausente "
                    f"del estado devuelto: {sorted(estado)}")

    def test_banderas_son_booleanas(self):
        for nombre in MaterialRegistry.names():
            with self.subTest(material=nombre):
                clase = MaterialRegistry.get(nombre)
                self.assertIsInstance(clase.IS_SYMMETRIC, bool)
                self.assertIsInstance(clase.IS_UNILATERAL, bool)

    def test_density_es_opcional_en_el_constructor(self):
        """ADR 0008: sólo peso propio y masa la exigen.

        Un material que la hiciera obligatoria rompería todo análisis estático
        que no la declara.
        """
        for nombre in MaterialRegistry.names():
            with self.subTest(material=nombre):
                material = _construir(nombre)
                self.assertIsNone(
                    material.density,
                    f"{nombre} asigna densidad sin que se la pidan; el default "
                    f"debe ser None para que el fallo por masa ausente sea "
                    f"explícito (ADR 0008)")

    def test_density_se_acepta_cuando_se_declara(self):
        for nombre in MaterialRegistry.names():
            with self.subTest(material=nombre):
                firma = inspect.signature(MaterialRegistry.get(nombre).__init__)
                self.assertIn(
                    'density', firma.parameters,
                    f"{nombre} no admite `density`: quedaría excluido de peso "
                    f"propio y análisis dinámico")
                material = MaterialRegistry.create(
                    nombre, **MUESTRAS[nombre], density=2400.0)
                self.assertEqual(material.density, 2400.0)


class TestContratoDeComputeState(unittest.TestCase):
    """Firma y valores de retorno de `compute_state`, el único método del
    contrato que los elementos invocan."""

    def test_firma_acepta_estado_opcional(self):
        """`compute_state(strain, state_vars=None)`.

        ADR 0013 §2 planea ampliarla con `orientation=None` para migrar la
        orientación material al elemento. Que el segundo parámetro tenga
        default es lo que hace esa migración aditiva y no ruptura.
        """
        for nombre in MaterialRegistry.names():
            with self.subTest(material=nombre):
                firma = inspect.signature(
                    MaterialRegistry.get(nombre).compute_state)
                params = [p for p in firma.parameters if p != 'self']
                self.assertGreaterEqual(len(params), 2)
                segundo = firma.parameters[params[1]]
                self.assertIsNot(
                    segundo.default, inspect.Parameter.empty,
                    f"{nombre}: el estado debe ser opcional")

    def test_devuelve_terna_con_formas_correctas(self):
        for nombre in MaterialRegistry.names():
            with self.subTest(material=nombre):
                material = _construir(nombre)
                dim = material.STRAIN_DIM
                resultado = material.compute_state(
                    _deformacion_pequena(dim), None)

                self.assertEqual(len(resultado), 3,
                                 f"{nombre}: se esperan (sigma, C, state_vars)")
                sigma, C, _ = resultado
                sigma = np.atleast_1d(np.asarray(sigma, dtype=float))
                C = np.atleast_2d(np.asarray(C, dtype=float))

                self.assertEqual(sigma.shape, (dim,),
                                 f"{nombre}: σ incompatible con STRAIN_DIM")
                self.assertEqual(C.shape, (dim, dim),
                                 f"{nombre}: C incompatible con STRAIN_DIM")

    def test_valores_finitos(self):
        """Ni NaN ni inf: un valor no finito envenena el ensamblaje entero y
        aparece mucho después como una factorización fallida."""
        for nombre in MaterialRegistry.names():
            with self.subTest(material=nombre):
                material = _construir(nombre)
                sigma, C, _ = material.compute_state(
                    _deformacion_pequena(material.STRAIN_DIM), None)
                self.assertTrue(np.all(np.isfinite(np.asarray(sigma, float))))
                self.assertTrue(np.all(np.isfinite(np.asarray(C, float))))

    def test_deformacion_nula_no_produce_esfuerzo(self):
        """Estado natural sin esfuerzos iniciales: σ(0) = 0.

        Todos los materiales del catálogo parten de un estado libre de
        esfuerzos. Un pretensado o una deformación térmica inicial romperían
        esto legítimamente — y entonces este test es el sitio donde documentar
        la excepción, no un obstáculo que sortear.
        """
        for nombre in MaterialRegistry.names():
            with self.subTest(material=nombre):
                material = _construir(nombre)
                dim = material.STRAIN_DIM
                sigma, C, _ = material.compute_state(np.zeros(dim), None)
                sigma = np.atleast_1d(np.asarray(sigma, dtype=float))
                escala = np.abs(np.atleast_2d(np.asarray(C, float))).max()
                np.testing.assert_allclose(
                    sigma, np.zeros(dim), rtol=0.0, atol=1e-12 * escala,
                    err_msg=f"{nombre}: σ(ε=0) ≠ 0")

    def test_no_muta_la_deformacion_recibida(self):
        """El elemento reutiliza el array de deformación entre puntos de Gauss;
        mutarlo in situ corrompería los siguientes de forma difícil de rastrear.
        """
        for nombre in MaterialRegistry.names():
            with self.subTest(material=nombre):
                material = _construir(nombre)
                eps = _deformacion_pequena(material.STRAIN_DIM)
                copia = eps.copy()
                material.compute_state(eps, None)
                np.testing.assert_array_equal(
                    eps, copia, err_msg=f"{nombre} mutó el strain recibido")

    def test_tangente_elastica_es_simetrica_y_definida_positiva(self):
        """En régimen elástico toda tangente del catálogo es SPD.

        Los materiales con `IS_SYMMETRIC = False` (daño isótropo, Drucker-Prager
        no asociado) lo son por su tangente **en carga activa**: la bandera es
        una cota superior conservadora que informa al despachador algebraico
        (ADR 0003) de que no puede asumir Cholesky, no una promesa de asimetría
        punto a punto. Aquí se sondea el régimen elástico, donde todas
        coinciden, y por eso la comprobación es incondicional.
        """
        for nombre in MaterialRegistry.names():
            with self.subTest(material=nombre):
                material = _construir(nombre)
                _, C, _ = material.compute_state(
                    _deformacion_pequena(material.STRAIN_DIM), None)
                C = np.atleast_2d(np.asarray(C, dtype=float))

                escala = np.abs(C).max()
                np.testing.assert_allclose(
                    C, C.T, rtol=0.0, atol=1e-12 * escala,
                    err_msg=f"{nombre}: tangente elástica no simétrica")
                self.assertGreater(
                    np.linalg.eigvalsh(0.5 * (C + C.T)).min(), 0.0,
                    f"{nombre}: tangente elástica no definida positiva")

    def test_is_symmetric_no_subdeclara(self):
        """`IS_SYMMETRIC = True` sí es una promesa firme, y se verifica.

        La bandera gobierna la elección de Cholesky en el despachador (ADR
        0003): declararla `True` con una tangente asimétrica produciría una
        factorización inválida. Al revés no es defecto —declarar `False` y ser
        simétrico sólo renuncia a una optimización—, así que sólo se exige la
        dirección que puede romper algo.
        """
        for nombre in MaterialRegistry.names():
            with self.subTest(material=nombre):
                material = _construir(nombre)
                if not material.IS_SYMMETRIC:
                    continue

                _, C, _ = material.compute_state(
                    _deformacion_pequena(material.STRAIN_DIM), None)
                C = np.atleast_2d(np.asarray(C, dtype=float))
                escala = np.abs(C).max()
                np.testing.assert_allclose(
                    C, C.T, rtol=0.0, atol=1e-12 * escala,
                    err_msg=f"{nombre} declara IS_SYMMETRIC=True pero su "
                            f"tangente no es simétrica: el despachador "
                            f"elegiría Cholesky sobre una matriz inválida")

    def test_estado_none_es_aceptable_en_el_primer_paso(self):
        """El primer `compute_state` de un análisis llega con `state_vars=None`
        porque aún no hay historia. Ningún material puede exigir estado previo.
        """
        for nombre in MaterialRegistry.names():
            with self.subTest(material=nombre):
                material = _construir(nombre)
                material.compute_state(
                    _deformacion_pequena(material.STRAIN_DIM), None)

    def test_estado_devuelto_se_readmite(self):
        """El estado que un material produce debe poder realimentarlo.

        Es el ciclo que el solver ejecuta en cada iteración de Newton: el
        `state_vars` trial de una iteración entra como estado de la siguiente.
        """
        for nombre in MaterialRegistry.names():
            with self.subTest(material=nombre):
                material = _construir(nombre)
                eps = _deformacion_pequena(material.STRAIN_DIM)
                _, _, estado = material.compute_state(eps, None)
                sigma, C, _ = material.compute_state(eps, estado)
                self.assertTrue(np.all(np.isfinite(np.asarray(sigma, float))))
                self.assertTrue(np.all(np.isfinite(np.asarray(C, float))))


class TestLinealidadDeLosElasticos(unittest.TestCase):
    """Los materiales sin variables internas son lineales: σ(k·ε) = k·σ(ε).

    Se identifican por `PRIMARY_STATE_VAR is None`, que es precisamente la
    declaración de "sin historia". Comprobarlo cierra el círculo entre lo
    declarado y el comportamiento: un material con historia real que declarase
    `None` fallaría aquí.

    `CableMaterial1D` queda excluido por ser unilateral —sin rigidez en
    compresión—, lo que rompe la homogeneidad para escalares negativos por
    diseño, no por defecto. La bandera `IS_UNILATERAL` es justamente lo que lo
    hace detectable sin nombrarlo a mano.
    """

    def test_homogeneidad(self):
        for nombre in MaterialRegistry.names():
            material = _construir(nombre)
            if material.PRIMARY_STATE_VAR is not None:
                continue
            if material.IS_UNILATERAL:
                continue

            with self.subTest(material=nombre):
                eps = _deformacion_pequena(material.STRAIN_DIM)
                sigma_1, _, _ = material.compute_state(eps, None)
                for k in (2.0, -3.5, 100.0):
                    sigma_k, _, _ = material.compute_state(k * eps, None)
                    escala = np.abs(np.asarray(sigma_1, float)).max()
                    np.testing.assert_allclose(
                        np.asarray(sigma_k, float),
                        k * np.asarray(sigma_1, float),
                        rtol=1e-12, atol=1e-12 * escala,
                        err_msg=f"{nombre} no es lineal pese a declararse sin "
                                f"variables internas (factor {k})")

    def test_tangente_constante(self):
        for nombre in MaterialRegistry.names():
            material = _construir(nombre)
            if material.PRIMARY_STATE_VAR is not None:
                continue
            if material.IS_UNILATERAL:
                continue

            with self.subTest(material=nombre):
                dim = material.STRAIN_DIM
                _, C_ref, _ = material.compute_state(np.zeros(dim), None)
                C_ref = np.atleast_2d(np.asarray(C_ref, dtype=float))
                for factor in (1.0e-8, 1.0e-4, 1.0e-2):
                    _, C, _ = material.compute_state(
                        _deformacion_pequena(dim, factor), None)
                    np.testing.assert_allclose(
                        np.atleast_2d(np.asarray(C, float)), C_ref,
                        rtol=1e-12,
                        err_msg=f"{nombre}: tangente no constante pese a "
                                f"declararse sin variables internas")


# Deuda de validación detectada por este barrido el 2026-09-10, la primera vez
# que un test recorrió el registro completo. NO son excepciones legítimas: son
# materiales que no validan lo que sus pares sí validan, y la inconsistencia es
# histórica (cada material se escribió en su momento con el criterio de
# entonces), no una decisión de diseño.
#
# Se listan aquí en vez de silenciar el test para que la deuda quede contada y
# visible. Cuando un material se corrija, se retira de la lista y el test lo
# exige a partir de entonces — la lista sólo puede encogerse.
#
# Pendiente de decisión del usuario: tocar seis materiales validados es un
# cambio transversal sobre código con física ya verificada, y no corresponde a
# la IA emprenderlo por iniciativa propia (Reglas.md §2 y §3).

SIN_VALIDAR_MODULO = frozenset({
    'Elastoplastic1D',      # plastic_1d.py: asigna E sin comprobar signo
})

SIN_VALIDAR_DENSIDAD = frozenset({
    'CableMaterial1D',
    'DruckerPrager2D',      # su par 3D sí valida — inconsistencia entre gemelos
    'Elastoplastic1D',
    'IsotropicDamage1D',
    'IsotropicDamage2D',    # su par 3D sí valida — inconsistencia entre gemelos
    'VonMises2D',           # su par 3D sí valida — inconsistencia entre gemelos
})


class TestValidacionDefensiva(unittest.TestCase):
    """Validación mínima común a todo material del catálogo.

    Un módulo elástico no positivo produce rigidez singular o negativa; una
    densidad negativa, masa negativa. Sin validación al construir, ambos
    aparecen mucho más tarde como una factorización imposible o un modo de
    vibración imaginario, lejos de su causa.

    `Reglas.md §1` lo pide explícitamente: validación temprana al construir con
    mensajes claros, sobre fallos crípticos en runtime.
    """

    def test_rechaza_modulo_no_positivo(self):
        for nombre in MaterialRegistry.names():
            if nombre in SIN_VALIDAR_MODULO:
                continue
            with self.subTest(material=nombre):
                muestra = dict(MUESTRAS[nombre])
                clave = 'E' if 'E' in muestra else 'E1'
                muestra[clave] = -1.0
                with self.assertRaises(
                        ValueError,
                        msg=f"{nombre} acepta {clave} negativo sin protestar"):
                    MaterialRegistry.create(nombre, **muestra)

    def test_rechaza_densidad_negativa(self):
        for nombre in MaterialRegistry.names():
            if nombre in SIN_VALIDAR_DENSIDAD:
                continue
            with self.subTest(material=nombre):
                with self.assertRaises(
                        ValueError,
                        msg=f"{nombre} acepta densidad negativa"):
                    MaterialRegistry.create(nombre, **MUESTRAS[nombre],
                                            density=-1.0)

    def test_las_listas_de_deuda_solo_encogen(self):
        """Trinquete: si un material ya validaba, no puede dejar de hacerlo.

        Comprueba que ningún nombre de las listas de deuda esté ya corregido
        —en cuyo caso hay que retirarlo de la lista— y que no se hayan colado
        nombres de materiales inexistentes.
        """
        registrados = set(MaterialRegistry.names())

        for lista, etiqueta in ((SIN_VALIDAR_MODULO, 'módulo'),
                                (SIN_VALIDAR_DENSIDAD, 'densidad')):
            fantasmas = lista - registrados
            self.assertFalse(
                fantasmas,
                f"La lista de deuda de {etiqueta} nombra materiales que ya no "
                f"se registran: {sorted(fantasmas)}")

        for nombre in sorted(SIN_VALIDAR_MODULO):
            with self.subTest(material=nombre, criterio='módulo'):
                muestra = dict(MUESTRAS[nombre])
                clave = 'E' if 'E' in muestra else 'E1'
                muestra[clave] = -1.0
                try:
                    MaterialRegistry.create(nombre, **muestra)
                except ValueError:
                    self.fail(
                        f"{nombre} ya valida el módulo: retíralo de "
                        f"SIN_VALIDAR_MODULO para que el test lo exija.")

        for nombre in sorted(SIN_VALIDAR_DENSIDAD):
            with self.subTest(material=nombre, criterio='densidad'):
                try:
                    MaterialRegistry.create(nombre, **MUESTRAS[nombre],
                                            density=-1.0)
                except ValueError:
                    self.fail(
                        f"{nombre} ya valida la densidad: retíralo de "
                        f"SIN_VALIDAR_DENSIDAD para que el test lo exija.")


if __name__ == '__main__':
    unittest.main()
