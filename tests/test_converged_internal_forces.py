"""Fuerzas internas del paso confirmado (ADR 0040, punto 4).

``Assembler.internal_forces_at(U)`` devuelve, sin ensamblar, la ``F_int`` con
la que el criterio aceptó el paso si ``U`` es el estado confirmado y nada lo
cambió después; si no, ensambla. Invariantes:

1. La historia no ensambla de más: con columnas de reacción y energía, los
   ensamblajes de un análisis son los mismos que sin historia (más el de
   validación de las energías en U = 0).
2. Con daño (implícito o IMPL-EX) la ``F_int`` registrada es idéntica a la de
   un reensamblaje tras el commit; con *return mapping* (J2) difiere sólo en
   el redondeo: el arquitecto midió ≤ 5·10⁻¹⁶ relativo; la cota de la prueba,
   10⁻¹³ del máximo de ``|F_int|``, deja margen para otras plataformas sin
   dejar pasar un error de formulación (que sería del orden de la
   tolerancia de convergencia, 10⁻⁶ o más).
3. El registro lo invalidan una U distinta en un bit, una evaluación nueva,
   ``prepare_all_steps``, ``set_pseudo_time_increment`` e ``invalidate()``; y
   es de sólo lectura.
4. Los resultados no dependen de la salida: U idéntica con y sin columnas de
   historia.
5. La reacción de ``SolveResult`` sale de la misma ``F_int`` que la última
   fila de la historia.
"""
import csv
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

import solidum
from solidum.math.assembly import Assembler

HISTORIA = """\
  history:
    columns:
      - {name: u, displacement: {group: xmax, dof: ux}}
      - {name: P, reaction: {group: xmax, dof: ux}}
"""
# Energías: sólo con materiales que descargan por la secante (ADR 0038).
ENERGIAS = """\
      - {name: W, energy: external_work}
      - {name: D, energy: dissipated}
"""


def _texto(material, pasos=30, u=0.04, historia=True):
    salida = "output:\n  file_name: m\n  log_level: ERROR\n"
    if historia:
        salida += HISTORIA + (ENERGIAS if "Damage" in material else "")
    return f"""\
materials:
  - {material}
mesh_element: Quad4
mesh_thickness: 10.0
mesh:
  structured:
    blocks:
      - {{material: 1, x: {{length: 60, n: 6}}, y: {{length: 20, n: 2}}}}
    groups:
      Origen: "near(x, 0) and near(y, 0)"
boundary_conditions_by_group:
  xmin: {{ux: 0.0}}
  Origen: {{uy: 0.0}}
  xmax: {{ux: {u}}}
solver: {{type: NonlinearSolver, num_steps: {pasos}}}
""" + salida


def _dano(integracion):
    return ("{id: 1, type: IsotropicDamageCrackBand2D, E: 30000.0, nu: 0.0, f_t: 3.0, "
            f"G_f: 0.1, integration: {integracion}}}")


J2 = "{id: 1, type: VonMises2D, E: 70.0, nu: 0.3, sigma_y: 0.243, H: 0.135, hypothesis: plane_strain}"


def _correr(texto):
    """Corre el YAML escribiendo su salida; devuelve el modelo, el resultado
    y la carpeta."""
    d = Path(tempfile.mkdtemp())
    ruta = d / "m.yaml"
    ruta.write_text(texto, encoding="utf-8")
    modelo = solidum.load_yaml(ruta)
    return modelo, modelo.solve(write_output=True), d


class _ContarEnsamblajes:
    def __enter__(self):
        original, self.n = Assembler.assemble_non_linear_system, 0

        def contar(asm, U):
            self.n += 1
            return original(asm, U)
        self._p = mock.patch.object(Assembler, "assemble_non_linear_system", contar)
        self._p.start()
        return self

    def __exit__(self, *exc):
        self._p.stop()


class TestSinEnsamblajesDeMas(unittest.TestCase):

    def test_la_historia_no_ensambla_en_cada_paso(self):
        texto_con, texto_sin = _texto(_dano("impl-ex")), _texto(_dano("impl-ex"), historia=False)
        with _ContarEnsamblajes() as con:
            _correr(texto_con)
        with _ContarEnsamblajes() as sin:
            _correr(texto_sin)
        self.assertEqual(con.n, sin.n + 1)          # el de validar las energías en U = 0


class TestIgualAUnReensamblaje(unittest.TestCase):

    def _comparar(self, material):
        """En el último paso, justo después del commit (la historia y las
        reacciones lo consultan ahí): el registro sin ensamblar y un
        reensamblaje."""
        d = Path(tempfile.mkdtemp())
        ruta = d / "m.yaml"
        ruta.write_text(_texto(material, historia=False), encoding="utf-8")
        modelo = solidum.load_yaml(ruta)
        a, capt = modelo.assembler, {}

        def cb(step, U, lam):
            if lam >= 1.0 - 1e-12:
                with _ContarEnsamblajes() as k:
                    capt["registrada"] = np.array(a.internal_forces_at(U))
                capt["n"] = k.n
                capt["reensamblada"] = a.assemble_non_linear_system(U)[1].copy()
        modelo.solve(step_callback=cb)
        self.assertEqual(capt["n"], 0)               # vino del registro
        return capt["registrada"], capt["reensamblada"]

    def test_con_dano_identica_bit_a_bit(self):
        for integracion in ("implicit", "impl-ex"):
            with self.subTest(integracion=integracion):
                registrada, reensamblada = self._comparar(_dano(integracion))
                self.assertTrue(np.array_equal(registrada, reensamblada))

    def test_con_j2_igual_salvo_el_redondeo(self):
        registrada, reensamblada = self._comparar(J2)
        escala = float(np.max(np.abs(reensamblada)))
        self.assertGreater(escala, 0.0)
        np.testing.assert_allclose(registrada, reensamblada, rtol=0.0, atol=1e-13 * escala)


class TestInvalidacion(unittest.TestCase):

    def setUp(self):
        self.modelo, self.res, _ = _correr(_texto(_dano("impl-ex"), pasos=10, historia=False))
        self.a = self.modelo.assembler
        self.U = np.array(self.res.U)
        # Un paso convergido: evaluación en U y commit.
        self.a.assemble_non_linear_system(self.U)
        self.a.commit_all_states()

    def _ensambla(self, U):
        with _ContarEnsamblajes() as k:
            self.a.internal_forces_at(U)
        return k.n

    def test_vigente_tras_el_commit_y_de_solo_lectura(self):
        self.assertEqual(self._ensambla(self.U), 0)
        F = self.a.internal_forces_at(self.U)
        with self.assertRaises(ValueError):
            F[0] = 1.0

    def test_lo_invalidan_los_cambios_de_estado(self):
        otro = self.U.copy()
        i = int(np.argmax(np.abs(otro)))
        otro[i] = np.nextafter(otro[i], np.inf)
        self.assertEqual(self._ensambla(otro), 1)            # U distinta en un bit
        self.assertEqual(self._ensambla(self.U), 1)          # evaluación nueva en medio

        # prepare_all_steps al final: deja abierto un intento de paso.
        for cambio in ("invalidate", "set_pseudo_time_increment", "prepare_all_steps"):
            with self.subTest(cambio=cambio):
                self.a.assemble_non_linear_system(self.U)
                self.a.commit_all_states()
                self.assertEqual(self._ensambla(self.U), 0)
                if cambio == "prepare_all_steps":
                    self.a.prepare_all_steps(self.U)
                elif cambio == "set_pseudo_time_increment":
                    self.a.set_pseudo_time_increment(0.1)
                else:
                    self.a.invalidate()
                self.assertEqual(self._ensambla(self.U), 1)

    def test_un_commit_sin_evaluacion_no_deja_registro(self):
        self.a.commit_all_states()                           # ninguna evaluación desde el anterior
        self.assertEqual(self._ensambla(self.U), 1)


class TestResultadosIndependientesDeLaSalida(unittest.TestCase):

    def test_U_identica_con_y_sin_historia(self):
        for material in (_dano("impl-ex"), _dano("implicit"), J2):
            with self.subTest(material=material[:40]):
                _, con, _ = _correr(_texto(material))
                _, sin, _ = _correr(_texto(material, historia=False))
                self.assertTrue(np.array_equal(con.U, sin.U))

    def test_las_reacciones_de_solve_result_salen_del_registro(self):
        """``build_solve_result`` no ensambla: toma la misma ``F_int`` que la
        historia en el último paso, y la reacción coincide con la última
        fila."""
        from solidum import entry
        original, conteo = entry.build_solve_result, {}

        def envoltura(*args, **kwargs):
            with _ContarEnsamblajes() as k:
                resultado = original(*args, **kwargs)
            conteo["n"] = k.n
            return resultado
        with mock.patch.object(entry, "build_solve_result", envoltura):
            modelo, res, d = _correr(_texto(J2))
        self.assertEqual(conteo["n"], 0)
        with open(d / "m_historia.csv", newline="", encoding="utf-8") as f:
            P_csv = list(csv.reader(f))[-1][3]
        dofs = [n.dofs["ux"] for n in modelo.domain.nodes.values()
                if abs(n.coordinates[0] - 60.0) < 1e-9]
        self.assertEqual(f"{float(np.sum(res.R[dofs])):.10g}", P_csv)

if __name__ == "__main__":
    unittest.main()
