"""Barrido del contrato de solvers sobre **todo** el registro.

Tercera pieza, junto a `tests/test_material_contract_sweep.py` y
`tests/test_element_contract_sweep.py`, de la misma idea: ningún test del
proyecto iteraba los registros, así que cada componente nuevo entraba en los
barridos genéricos sólo si alguien lo añadía a mano.

Qué gobierna `PIPELINE_KIND`
----------------------------
Es el atributo declarativo que decide **por qué entrypoint** se despacha un
solver (regla C del ADR 0009). `solidum/entry.py` mantiene el conjunto canónico
`_KNOWN_PIPELINE_KINDS` y rechaza en `run_yaml` cualquier solver con un valor
fuera de él. Un solver nuevo que lo declare mal —o que no lo declare— no es
alcanzable desde YAML: no falla al importarse, simplemente queda inaccesible
para el usuario. Este barrido convierte eso en un fallo ruidoso al registrarlo.

Qué se verifica
---------------
Contrato estructural, no la convergencia ni la física de cada esquema temporal
— eso lo cubren los tests dedicados de cada solver y los benchmarks de
`tests/validation/`. Aquí: que el `PIPELINE_KIND` declarado sea despachable,
que la firma de construcción y de `solve` sea la que el entrypoint espera, y
que la clasificación estática/transitoria concuerde con los parámetros que el
solver realmente pide.

Por qué no se ejecuta cada solver
---------------------------------
Resolver un problema con los trece exigiría un modelo válido por familia
(estático, modal, transitorio, armónico, espectral, térmico), con condiciones
de contorno y parámetros de paso ajustados a cada esquema. Eso ya existe,
repartido en los tests dedicados, y duplicarlo aquí sería frágil sin añadir
cobertura. El valor de este archivo está en lo que ningún test dedicado mira:
el contrato **transversal** que hace despachable a un solver.
"""
import inspect
import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import solidum  # noqa: F401  — dispara el autodiscover que puebla el registro
from solidum.entry import _KNOWN_PIPELINE_KINDS
from solidum.registry import SolverRegistry

# Parámetros que cada familia de pipeline exige más allá del assembler.
# Es la clasificación que `entry.py` presupone al despachar: un solver
# transitorio sin `dt` no podría invocarse desde `run_transient`.
PARAMETROS_ESPERADOS = {
    'static': set(),
    'modal': {'n_modes'},
    'transient': {'t_end', 'dt'},
    'harmonic': set(),
    'spectrum': {'n_modes', 'direction', 'spectrum'},
    'thermal_transient': {'dt', 'n_steps', 'T_initial'},
}


def _firma_init(nombre):
    return inspect.signature(SolverRegistry.get(nombre).__init__)


def _parametros_requeridos(nombre):
    """Parámetros obligatorios del constructor, sin `self` ni `assembler`."""
    firma = _firma_init(nombre)
    return {
        p for p, v in firma.parameters.items()
        if p not in ('self', 'assembler')
        and v.default is inspect.Parameter.empty
        and v.kind not in (v.VAR_POSITIONAL, v.VAR_KEYWORD)
    }


class TestRegistroCompleto(unittest.TestCase):

    def test_el_registro_no_esta_vacio(self):
        """Si el autodiscover se rompiera, el barrido pasaría por vacuidad."""
        self.assertGreaterEqual(len(SolverRegistry.names()), 13)

    def test_toda_familia_de_pipeline_tiene_al_menos_un_solver(self):
        """Un `PIPELINE_KIND` canónico sin solver que lo implemente señala un
        entrypoint muerto en `entry.py`, o un solver que se dejó de registrar.
        """
        declarados = {
            getattr(SolverRegistry.get(n), 'PIPELINE_KIND', None)
            for n in SolverRegistry.names()
        }
        huerfanas = _KNOWN_PIPELINE_KINDS - declarados
        self.assertFalse(
            huerfanas,
            f"Familias de pipeline sin ningún solver registrado: "
            f"{sorted(huerfanas)}")


class TestContratoDeclarativo(unittest.TestCase):
    """`PIPELINE_KIND` y la firma que el entrypoint presupone."""

    def test_pipeline_kind_declarado(self):
        """Sin él, `run_yaml` no sabe por qué entrypoint despachar el solver."""
        for nombre in SolverRegistry.names():
            with self.subTest(solver=nombre):
                clase = SolverRegistry.get(nombre)
                self.assertTrue(
                    hasattr(clase, 'PIPELINE_KIND'),
                    f"{nombre} no declara PIPELINE_KIND: quedaría inalcanzable "
                    f"desde YAML sin error visible al importarlo")

    def test_pipeline_kind_es_despachable(self):
        """El valor debe estar en el conjunto canónico de `entry.py`.

        Un typo (`'estatic'`, `'transiente'`) no rompe nada al registrar el
        solver — lo deja inaccesible, que es un fallo mucho más difícil de
        diagnosticar que una excepción.
        """
        for nombre in SolverRegistry.names():
            with self.subTest(solver=nombre):
                kind = SolverRegistry.get(nombre).PIPELINE_KIND
                self.assertIn(
                    kind, _KNOWN_PIPELINE_KINDS,
                    f"{nombre}: PIPELINE_KIND={kind!r} no está entre los "
                    f"despachables {sorted(_KNOWN_PIPELINE_KINDS)}")

    def test_primer_parametro_es_el_assembler(self):
        """Todos los entrypoints construyen `Solver(assembler, ...)`."""
        for nombre in SolverRegistry.names():
            with self.subTest(solver=nombre):
                params = [p for p in _firma_init(nombre).parameters
                          if p != 'self']
                self.assertTrue(params, f"{nombre}: constructor sin parámetros")
                self.assertEqual(
                    params[0], 'assembler',
                    f"{nombre}: el primer parámetro debe ser `assembler`")

    def test_parametros_coherentes_con_la_familia(self):
        """Lo declarado y lo que el constructor pide deben concordar.

        Un solver marcado como `transient` que no acepte `dt` no es invocable
        por `run_transient`; uno marcado `static` que exija `t_end` tampoco lo
        es por `run`. La discrepancia sólo aparecería al ejecutar el YAML.
        """
        for nombre in SolverRegistry.names():
            with self.subTest(solver=nombre):
                clase = SolverRegistry.get(nombre)
                # Un PIPELINE_KIND fuera del canon ya lo caza el test anterior;
                # aquí se salta para que el fallo salga con su mensaje propio y
                # no como un KeyError crudo de esta tabla.
                if clase.PIPELINE_KIND not in PARAMETROS_ESPERADOS:
                    continue
                esperados = PARAMETROS_ESPERADOS[clase.PIPELINE_KIND]
                requeridos = _parametros_requeridos(nombre)

                faltantes = esperados - requeridos
                self.assertFalse(
                    faltantes,
                    f"{nombre} ({clase.PIPELINE_KIND}) no pide {sorted(faltantes)}, "
                    f"que su entrypoint le pasa")

    def test_los_estaticos_no_piden_parametros_temporales(self):
        """Un solver estático que exigiera `dt` revelaría una clasificación
        equivocada — y quedaría inconstruible desde el entrypoint estático."""
        temporales = {'dt', 't_end', 'n_steps'}
        for nombre in SolverRegistry.names():
            clase = SolverRegistry.get(nombre)
            if clase.PIPELINE_KIND != 'static':
                continue
            with self.subTest(solver=nombre):
                intrusos = _parametros_requeridos(nombre) & temporales
                self.assertFalse(
                    intrusos,
                    f"{nombre} se declara `static` pero exige {sorted(intrusos)}")


class TestContratoDeSolve(unittest.TestCase):
    """`solve` es el método que todos los entrypoints invocan."""

    def test_expone_solve(self):
        for nombre in SolverRegistry.names():
            with self.subTest(solver=nombre):
                clase = SolverRegistry.get(nombre)
                self.assertTrue(
                    callable(getattr(clase, 'solve', None)),
                    f"{nombre} no expone un `solve` invocable")

    def test_los_estaticos_reciben_las_cargas_en_solve(self):
        """`static` es la única familia que pasa el vector de fuerzas a `solve`.

        El resto lo recibe en el constructor —una historia temporal, un
        espectro, una excitación armónica no son un vector de cargas puntual—
        y `solve()` se invoca sin argumentos. La división es la que `entry.py`
        presupone al despachar; romperla dejaría el solver inalcanzable desde
        su entrypoint.
        """
        for nombre in SolverRegistry.names():
            clase = SolverRegistry.get(nombre)
            if clase.PIPELINE_KIND != 'static':
                continue
            with self.subTest(solver=nombre):
                params = [p for p in inspect.signature(clase.solve).parameters
                          if p != 'self']
                self.assertTrue(
                    params,
                    f"{nombre}.solve no acepta las cargas aplicadas")
                primero = inspect.signature(clase.solve).parameters[params[0]]
                self.assertIs(
                    primero.default, inspect.Parameter.empty,
                    f"{nombre}.solve: las cargas deben ser obligatorias")

    def test_los_no_estaticos_resuelven_sin_argumentos(self):
        """Transitorios, modales, armónicos y espectrales: `solve()` a secas.

        Sus datos de entrada llegan por el constructor, así que cualquier
        parámetro adicional de `solve` debe tener default; si no, el
        entrypoint —que lo llama sin argumentos— fallaría.
        """
        for nombre in SolverRegistry.names():
            clase = SolverRegistry.get(nombre)
            if clase.PIPELINE_KIND == 'static':
                continue
            with self.subTest(solver=nombre):
                firma = inspect.signature(clase.solve)
                obligatorios = [
                    p for p, v in firma.parameters.items()
                    if p != 'self'
                    and v.default is inspect.Parameter.empty
                    and v.kind not in (v.VAR_POSITIONAL, v.VAR_KEYWORD)
                ]
                self.assertFalse(
                    obligatorios,
                    f"{nombre}.solve exige {obligatorios}, pero su entrypoint "
                    f"lo invoca sin argumentos")

    def test_los_no_estaticos_devuelven_un_resultado_tipado(self):
        """Las dataclasses de `solidum.results` son la API pública del análisis.

        `run_yaml` valida el tipo devuelto; un solver que retornase un array
        crudo rompería a los consumidores del resultado sin fallar al
        registrarse.
        """
        import solidum.results as resultados

        tipos = tuple(
            obj for obj in vars(resultados).values()
            if isinstance(obj, type) and obj.__name__.endswith('Result')
        )
        self.assertTrue(tipos, "no se encontraron dataclasses de resultado")

        nombres_validos = {t.__name__ for t in tipos}
        for nombre in SolverRegistry.names():
            clase = SolverRegistry.get(nombre)
            if clase.PIPELINE_KIND == 'static':
                continue
            with self.subTest(solver=nombre):
                anotacion = inspect.signature(clase.solve).return_annotation
                self.assertIsNot(
                    anotacion, inspect.Signature.empty,
                    f"{nombre}.solve no anota su tipo de retorno")
                texto = (anotacion if isinstance(anotacion, str)
                         else getattr(anotacion, '__name__', str(anotacion)))
                self.assertTrue(
                    any(v in texto for v in nombres_validos),
                    f"{nombre}.solve devuelve {texto!r}, que no es una "
                    f"dataclass de resultado de `solidum.results`")

    def test_step_callback_es_opcional_donde_existe(self):
        """`_invoke_solve` en `entry.py` detecta por firma si pasar el callback.

        Un solver que lo declarase obligatorio rompería ese despacho: el
        entrypoint lo llama sin él cuando no hay callback que pasar.
        """
        for nombre in SolverRegistry.names():
            with self.subTest(solver=nombre):
                firma = inspect.signature(SolverRegistry.get(nombre).solve)
                if 'step_callback' not in firma.parameters:
                    continue
                self.assertIsNot(
                    firma.parameters['step_callback'].default,
                    inspect.Parameter.empty,
                    f"{nombre}.solve exige `step_callback`; debe ser opcional")


class TestCoberturaDeLaTabla(unittest.TestCase):
    """La tabla de familias se mantiene sincronizada con `entry.py`."""

    def test_la_tabla_cubre_todas_las_familias_canonicas(self):
        """Si `entry.py` añade una familia y esta tabla no, el test de
        coherencia de parámetros lanzaría `KeyError` en vez de fallar con un
        mensaje útil. Este test lo convierte en un fallo explícito.
        """
        sin_entrada = _KNOWN_PIPELINE_KINDS - set(PARAMETROS_ESPERADOS)
        self.assertFalse(
            sin_entrada,
            f"Familias de pipeline sin entrada en PARAMETROS_ESPERADOS: "
            f"{sorted(sin_entrada)}. Añádelas.")

        sobrantes = set(PARAMETROS_ESPERADOS) - _KNOWN_PIPELINE_KINDS
        self.assertFalse(
            sobrantes,
            f"PARAMETROS_ESPERADOS nombra familias que `entry.py` ya no "
            f"reconoce: {sorted(sobrantes)}")


if __name__ == '__main__':
    unittest.main()
