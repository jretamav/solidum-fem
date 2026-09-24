"""Carga de módulos de usuario (ADR 0020, P4).

Un módulo de usuario es un subpaquete de ``solidum/user/`` que el programa
principal no importa: lo carga el modelo que lo pide (``user_modules`` en el
YAML, ``solidum.load_user_module`` en Python). Se prueba con un módulo de
juguete escrito en un directorio temporal y añadido a ``solidum.user.__path__``
en un subproceso: sus registros no deben contaminar la sesión de pytest.
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
            _kind = "Resorte"
            SPEC_KIND = "resorte_juguete"
            YAML_SECTION = "resortes_juguete"
            YAML_LABEL = "resorte"


        @ResorteRegistry.register
        class ResorteLineal:
            def __init__(self, k: float):
                self.k = k
        '''),
}

_MODELO = textwrap.dedent("""\
    nodes:
      - {id: 1, coords: [0.0, 0.0]}
      - {id: 2, coords: [1.0, 0.0]}
    materials:
      - {id: 1, type: Elastic1D, E: 200.0e9}
    elements:
      - {id: 1, type: Truss2D, nodes: [1, 2], material: 1, A: 1.0e-3}
    resortes_juguete:
      - {id: 1, type: ResorteLineal, k: 2.0}
    """)

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


def test_modulo_de_usuario_se_carga_solo_cuando_el_modelo_lo_pide(tmp_path):
    raiz = tmp_path / "modulos"
    paquete = raiz / "juguete"
    paquete.mkdir(parents=True)
    for nombre, texto in _JUGUETE.items():
        (paquete / nombre).write_text(texto, encoding="utf-8")
    script = tmp_path / "prueba.py"
    script.write_text(_SCRIPT, encoding="utf-8")
    modelo = tmp_path / "modelo_base.yaml"
    modelo.write_text(_MODELO, encoding="utf-8")

    r = subprocess.run(
        [sys.executable, str(script), str(raiz), str(tmp_path), str(modelo)],
        cwd=ROOT, capture_output=True, text=True, encoding="utf-8", timeout=300,
    )
    assert r.returncode == 0, r.stdout + r.stderr
    assert r.stdout.strip().endswith("OK")


def test_modulo_inexistente_desde_python():
    import pytest
    import solidum
    with pytest.raises(ValueError, match="Módulo de usuario 'no_existe' desconocido"):
        solidum.load_user_module("no_existe")
