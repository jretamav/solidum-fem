"""Carga de módulos de usuario (ADR 0020, P4).

Un módulo de usuario es un subpaquete de ``solidum/user/`` que el programa
principal no importa: lo carga el modelo que lo pide (``user_modules`` en el
YAML, ``solidum.load_user_module`` en Python). Sólo puede añadir entradas a los
registros, nunca sustituir las existentes. Se prueba con módulos de juguete
escritos en un directorio temporal y añadidos a ``solidum.user.__path__`` en
un subproceso: sus registros no deben contaminar la sesión de pytest.
"""
import subprocess
import sys
import textwrap
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

_JUGUETE = {
    "__init__.py": '"""Módulo de usuario de juguete para las pruebas."""\n',
    "resortes.py": textwrap.dedent('''\
        from solidum.registry import Registry


        class ResorteRegistry(Registry):
            SPEC_KIND = "resorte_juguete"
            YAML_SECTION = "resortes_juguete"
            YAML_LABEL = "resorte"


        @ResorteRegistry.register
        class ResorteLineal:
            def __init__(self, k: float):
                self.k = k
        '''),
}

_BARRA = textwrap.dedent("""\
    nodes:
      - {id: 1, coords: [0.0, 0.0]}
      - {id: 2, coords: [1.0, 0.0]}
    materials:
      - {id: 1, type: Elastic1D, E: 200.0e9}
    elements:
      - {id: 1, type: Truss2D, nodes: [1, 2], material: 1, A: 1.0e-3}
    """)
_MODELO = _BARRA + "resortes_juguete:\n  - {id: 1, type: ResorteLineal, k: 2.0}\n"

_SCRIPT = textwrap.dedent('''\
    import sys
    from pathlib import Path

    import solidum
    import solidum.user

    solidum.user.__path__.append(sys.argv[1])
    carpeta = Path(sys.argv[2])
    from solidum.utils.yaml_parser import YamlParser, YamlValidationError

    def parsear(texto):
        ruta = carpeta / "modelo.yaml"
        ruta.write_text(texto, encoding="utf-8")
        parser = YamlParser(str(ruta))
        parser.parse()
        return parser

    modelo = Path(sys.argv[3]).read_text(encoding="utf-8")
    assert "juguete" in solidum.available_user_modules()
    assert not any(m.startswith("solidum.user.juguete") for m in sys.modules), \\
        "import solidum cargó el módulo de usuario"

    try:
        parsear(modelo)
        raise AssertionError("sin user_modules la sección del módulo debía fallar")
    except YamlValidationError as exc:
        assert "Secciones desconocidas: ['resortes_juguete']" in str(exc), exc
        assert "decláralo en 'user_modules'" in str(exc), exc
        assert "juguete" in str(exc), exc

    try:
        parsear(modelo + "user_modules: [inexistente]\\n")
        raise AssertionError("un módulo inexistente debía fallar")
    except YamlValidationError as exc:
        assert "Módulo de usuario 'inexistente' desconocido" in str(exc), exc

    parser = parsear(modelo + "user_modules: [juguete]\\n")
    assert parser.resortes_juguete[1].k == 2.0
    assert "solidum.user.juguete.resortes" in sys.modules

    solidum.load_user_module("juguete")        # dos veces: no hace nada
    parser = parsear(modelo + "user_modules: juguete\\n")
    assert parser.resortes_juguete[1].k == 2.0
    print("OK")
    ''')


def _ejecutar(tmp_path, nombre: str, archivos: dict, script: str, modelo: str) -> None:
    """Escribe el módulo de juguete ``nombre`` y ejecuta ``script`` en un
    subproceso con argumentos ``(carpeta de módulos, carpeta de trabajo,
    YAML base)``; el script termina imprimiendo ``OK``."""
    raiz = tmp_path / "modulos"
    paquete = raiz / nombre
    paquete.mkdir(parents=True)
    for archivo, texto in archivos.items():
        (paquete / archivo).write_text(texto, encoding="utf-8")
    ruta_script = tmp_path / "prueba.py"
    ruta_script.write_text(script, encoding="utf-8")
    ruta_modelo = tmp_path / "modelo_base.yaml"
    ruta_modelo.write_text(modelo, encoding="utf-8")

    r = subprocess.run(
        [sys.executable, str(ruta_script), str(raiz), str(tmp_path), str(ruta_modelo)],
        cwd=ROOT, capture_output=True, text=True, encoding="utf-8", timeout=300,
    )
    assert r.returncode == 0, r.stdout + r.stderr
    assert r.stdout.strip().endswith("OK")


def test_modulo_de_usuario_se_carga_solo_cuando_el_modelo_lo_pide(tmp_path):
    _ejecutar(tmp_path, "juguete", _JUGUETE, _SCRIPT, _MODELO)


# Módulo de juguete que intenta sustituir una clase del programa principal
# (además de añadir una familia y un material propios, para comprobar que la
# carga fallida se deshace entera).
_SOMBRA = {
    "__init__.py": '"""Módulo de usuario de juguete que reasigna Elastic2D."""\n',
    "sombra.py": textwrap.dedent('''\
        from solidum.materials.elastic_2d import Elastic2D
        from solidum.registry import MaterialRegistry, Registry


        class SombraRegistry(Registry):
            SPEC_KIND = "sombra_juguete"
            YAML_SECTION = "sombras_juguete"


        @MaterialRegistry.register
        class MaterialNuevoDeSombra(Elastic2D):
            pass


        @MaterialRegistry.register("Elastic2D")
        class ElasticoSombra(Elastic2D):
            pass
        '''),
}

_SCRIPT_SOMBRA = textwrap.dedent('''\
    import sys
    from pathlib import Path

    import solidum
    import solidum.user
    from solidum.materials.elastic_2d import Elastic2D
    from solidum.registry import MaterialRegistry, Registry

    solidum.user.__path__.append(sys.argv[1])
    carpeta = Path(sys.argv[2])
    from solidum.utils.yaml_parser import YamlParser, YamlValidationError

    def intacto():
        assert MaterialRegistry.get("Elastic2D") is Elastic2D
        assert "MaterialNuevoDeSombra" not in MaterialRegistry.names()
        assert "sombras_juguete" not in [f.YAML_SECTION for f in Registry.families()]
        assert Registry.for_spec_kind("sombra_juguete") is None
        assert not any(m.startswith("solidum.user.sombra") for m in sys.modules)

    for _ in range(2):             # el segundo intento vuelve a fallar
        try:
            solidum.load_user_module("sombra")
            raise AssertionError("reasignar Elastic2D debía fallar")
        except ValueError as exc:
            texto = str(exc)
            assert "Módulo de usuario 'sombra': reasigna entradas ya registradas" in texto, texto
            assert ("MaterialRegistry 'Elastic2D': solidum.materials.elastic_2d.Elastic2D "
                    "→ solidum.user.sombra.sombra.ElasticoSombra") in texto, texto
            assert "sólo puede añadir entradas" in texto, texto
        intacto()

    ruta = carpeta / "modelo.yaml"
    ruta.write_text(Path(sys.argv[3]).read_text(encoding="utf-8")
                    + "user_modules: [sombra]\\n", encoding="utf-8")
    try:
        YamlParser(str(ruta)).parse()
        raise AssertionError("el YAML con el módulo que reasigna debía fallar")
    except YamlValidationError as exc:
        assert "user_modules: Módulo de usuario 'sombra': reasigna" in str(exc), exc
    intacto()
    print("OK")
    ''')


def test_un_modulo_no_puede_reasignar_una_entrada_existente(tmp_path):
    """Un módulo de usuario sólo añade entradas a los registros. Si pudiera
    sustituir ``Elastic2D``, cambiaría en silencio el programa principal para
    todos los modelos posteriores de la sesión. La carga fallida se deshace
    entera, y un segundo intento vuelve a fallar."""
    _ejecutar(tmp_path, "sombra", _SOMBRA, _SCRIPT_SOMBRA, _BARRA)


def test_modulo_inexistente_desde_python():
    import pytest
    import solidum
    with pytest.raises(ValueError, match="Módulo de usuario 'no_existe' desconocido"):
        solidum.load_user_module("no_existe")
