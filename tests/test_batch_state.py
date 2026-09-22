"""Estado interno por arreglos (ADR 0014, fase 1).

Verifica la pieza de infraestructura que hace posible el ensamblaje por
lotes sin tocar todavía el ensamblador:

1. ``STATE_SCHEMA`` es obligatorio en todo material registrado (mecánico y
   térmico) y **coincide** con lo que ``compute_state`` devuelve — claves y
   formas. Un esquema que no describa el estado real produciría filas
   corruptas de forma silenciosa.
2. ``StateSchema`` empaqueta y desempaqueta sin pérdida.
3. ``BatchedElementState`` se comporta como ``ElementState``: las vistas
   devuelven lo mismo que las listas, asignar escribe la fila y ``commit``
   promueve trial → committed sólo para ese elemento.
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import solidum  # noqa: F401
from solidum.core.element_state import ElementState
from solidum.math.batch import BatchedElementState, FamilyState, StateSchema
from solidum.registry import MaterialRegistry, ThermalMaterialRegistry

from test_material_contract_sweep import MUESTRAS


def _strain_that_evolves(material):
    """Deformación grande para que las variables internas cambien."""
    dim = material.STRAIN_DIM
    if dim == 1:
        return 5.0e-2
    eps = np.zeros(dim)
    eps[0] = 5.0e-2
    eps[-1] = 2.0e-2
    return eps


class TestEsquemaObligatorio(unittest.TestCase):

    def test_todo_material_declara_state_schema(self):
        for nombre in MaterialRegistry.names():
            with self.subTest(material=nombre):
                clase = MaterialRegistry.get(nombre)
                self.assertIsInstance(
                    clase.STATE_SCHEMA, dict,
                    f"{nombre}: STATE_SCHEMA debe ser un dict ({{}} si no hay historia)")

    def test_materiales_termicos_declaran_state_schema(self):
        for nombre in ThermalMaterialRegistry.names():
            with self.subTest(material=nombre):
                clase = ThermalMaterialRegistry.get(nombre)
                self.assertIsInstance(clase.STATE_SCHEMA, dict)

    def test_esquema_coincide_con_el_estado_devuelto(self):
        """Claves y formas del dict devuelto == esquema, en régimen elástico
        y tras una deformación grande (rama inelástica)."""
        for nombre in MaterialRegistry.names():
            with self.subTest(material=nombre):
                material = MaterialRegistry.create(nombre, **MUESTRAS[nombre])
                schema = StateSchema(material.STATE_SCHEMA)
                eps = _strain_that_evolves(material)
                _, _, estado = material.compute_state(eps, None)
                schema.validate(estado)
                _, _, estado2 = material.compute_state(eps, estado)
                schema.validate(estado2)

    def test_esquema_coincide_con_primary_state_var(self):
        for nombre in MaterialRegistry.names():
            with self.subTest(material=nombre):
                clase = MaterialRegistry.get(nombre)
                if clase.PRIMARY_STATE_VAR is None:
                    self.assertEqual(clase.STATE_SCHEMA, {},
                                     f"{nombre}: sin PRIMARY_STATE_VAR pero con esquema")
                else:
                    self.assertIn(clase.PRIMARY_STATE_VAR, clase.STATE_SCHEMA)

    def test_estado_inicial_respeta_el_esquema(self):
        """``initial_state`` es lo que el kernel ve en el primer paso; debe
        ser equivalente a ``state_vars=None`` para ``compute_state``."""
        for nombre in MaterialRegistry.names():
            with self.subTest(material=nombre):
                material = MaterialRegistry.create(nombre, **MUESTRAS[nombre])
                schema = StateSchema(material.STATE_SCHEMA)
                inicial = material.initial_state()
                schema.validate(inicial)
                eps = _strain_that_evolves(material)
                s_none, C_none, st_none = material.compute_state(eps, None)
                s_ini, C_ini, st_ini = material.compute_state(eps, inicial)
                np.testing.assert_array_equal(np.asarray(s_none, float), np.asarray(s_ini, float))
                np.testing.assert_array_equal(np.asarray(C_none, float), np.asarray(C_ini, float))
                if st_none is not None:
                    for k in st_none:
                        np.testing.assert_array_equal(
                            np.asarray(st_none[k], float), np.asarray(st_ini[k], float))


class TestStateSchema(unittest.TestCase):

    def test_pack_unpack_sin_perdida(self):
        schema = StateSchema({'eps_p': (4,), 'alpha': (), 'M': (2, 3)})
        self.assertEqual(schema.n_state, 4 + 1 + 6)
        self.assertEqual(schema.offsets, (0, 4, 5))
        estado = {'eps_p': np.arange(4.0), 'alpha': 7.5, 'M': np.arange(6.0).reshape(2, 3)}
        fila = np.zeros(schema.n_state)
        schema.pack(estado, fila)
        vuelta = schema.unpack(fila)
        np.testing.assert_array_equal(vuelta['eps_p'], estado['eps_p'])
        self.assertEqual(vuelta['alpha'], 7.5)
        self.assertIsInstance(vuelta['alpha'], float)
        np.testing.assert_array_equal(vuelta['M'], estado['M'])
        # El dict desempaquetado es independiente de la fila.
        vuelta['eps_p'][0] = -1.0
        self.assertEqual(fila[0], 0.0)

    def test_esquema_vacio_devuelve_none(self):
        schema = StateSchema({})
        self.assertEqual(schema.n_state, 0)
        self.assertIsNone(schema.unpack(np.zeros(0)))
        schema.validate(None)

    def test_pack_none_restaura_el_estado_inicial(self):
        schema = StateSchema({'kappa': (), 'damage': ()})
        fila = np.array([3.0, 0.5])
        schema.pack(None, fila, initial=np.array([1.0e-4, 0.0]))
        np.testing.assert_array_equal(fila, [1.0e-4, 0.0])

    def test_pack_rechaza_estado_incompatible(self):
        schema = StateSchema({'eps_p': (4,), 'alpha': ()})
        fila = np.zeros(5)
        with self.assertRaises(ValueError):
            schema.pack({'eps_p': np.zeros(3), 'alpha': 0.0}, fila)
        with self.assertRaises(ValueError):
            schema.pack({'alpha': 0.0}, fila)
        with self.assertRaises(ValueError):
            schema.validate({'eps_p': np.zeros(4), 'alpha': 0.0, 'extra': 1.0})


class TestBatchedElementState(unittest.TestCase):

    def setUp(self):
        self.schema = StateSchema({'eps_p': (4,), 'alpha': ()})
        self.fs = FamilyState(n_elem=3, n_gp=4, schema=self.schema, n_sigma=3)
        self.estado = BatchedElementState(self.fs, elem_index=1)

    def test_se_comporta_como_element_state(self):
        self.assertIsInstance(self.estado, ElementState)
        self.assertEqual(self.estado.num_ip, 4)
        self.assertEqual(len(self.estado.vars), 4)
        self.assertEqual(len(self.estado.stresses_trial), 4)
        # Estado inicial (ceros) → dict con claves del esquema.
        v = self.estado.vars[0]
        self.assertEqual(set(v), {'eps_p', 'alpha'})
        np.testing.assert_array_equal(v['eps_p'], np.zeros(4))
        self.assertEqual(v['alpha'], 0.0)

    def test_asignar_escribe_la_fila_correcta(self):
        self.estado.vars_trial[2] = {'eps_p': np.array([1., 2., 3., 4.]), 'alpha': 0.25}
        self.estado.stresses_trial[2] = np.array([10., 20., 30.])
        fila = 1 * 4 + 2
        np.testing.assert_array_equal(self.fs.S_trial[fila], [1., 2., 3., 4., 0.25])
        np.testing.assert_array_equal(self.fs.sig_trial[fila], [10., 20., 30.])
        # El resto de filas no se toca.
        self.assertEqual(np.abs(self.fs.S_trial).sum(), 10.25)
        # Lectura de vuelta (índice negativo también).
        self.assertEqual(self.estado.vars_trial[-2]['alpha'], 0.25)
        np.testing.assert_array_equal(self.estado.stresses_trial[2], [10., 20., 30.])

    def test_commit_de_elemento_solo_afecta_a_ese_elemento(self):
        otro = BatchedElementState(self.fs, elem_index=0)
        otro.vars_trial[0] = {'eps_p': np.ones(4), 'alpha': 9.0}
        self.estado.vars_trial[1] = {'eps_p': 2 * np.ones(4), 'alpha': 1.0}
        self.estado.commit()
        self.assertEqual(self.estado.vars[1]['alpha'], 1.0)
        self.assertEqual(otro.vars[0]['alpha'], 0.0)      # sin commit
        self.fs.commit()
        self.assertEqual(otro.vars[0]['alpha'], 9.0)

    def test_commit_de_familia_deja_trial_igual_a_committed(self):
        self.estado.vars_trial[3] = {'eps_p': np.zeros(4), 'alpha': 4.0}
        self.fs.commit()
        self.assertEqual(self.estado.vars[3]['alpha'], 4.0)
        self.assertEqual(self.estado.vars_trial[3]['alpha'], 4.0)

    def test_adopt_y_to_plain_son_inversos(self):
        plano = ElementState(4, init_stress=np.zeros(3))
        plano.vars_trial[1] = {'eps_p': np.arange(4.0), 'alpha': 0.5}
        plano.stresses_trial[1] = np.array([1., 2., 3.])
        plano.commit()
        self.estado.adopt(plano)
        vuelta = self.estado.to_plain()
        self.assertEqual(vuelta.vars[1]['alpha'], 0.5)
        np.testing.assert_array_equal(vuelta.vars[1]['eps_p'], np.arange(4.0))
        np.testing.assert_array_equal(vuelta.stresses[1], [1., 2., 3.])
        self.assertEqual(vuelta.vars[0]['alpha'], 0.0)

    def test_esquema_vacio_devuelve_none_como_hoy(self):
        fs = FamilyState(2, 1, StateSchema({}), 3)
        st = BatchedElementState(fs, 0)
        self.assertIsNone(st.vars[0])
        st.vars_trial[0] = None       # lo que devuelve un material sin historia
        self.assertIsNone(st.vars_trial[0])

    def test_memoria_por_punto_de_gauss(self):
        """80 B/gp para J2 2D (4+1 estado ×2 + 3 σ ×2, en dobles) + 1 B de flag."""
        fs = FamilyState(1000, 4, self.schema, 3)
        self.assertEqual(fs.nbytes, 4000 * (2 * 5 * 8 + 2 * 3 * 8 + 1))


if __name__ == '__main__':
    unittest.main()
