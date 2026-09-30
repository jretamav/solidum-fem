"""Pureza del programa principal frente a los módulos de usuario (ADR 0020).

El programa principal es Solidum estándar: todo ``solidum/`` salvo los
módulos de usuario ``solidum/user/<nombre>/``. No los importa ni los nombra,
como el programa principal de FEAP no conoce los elementos de usuario que se
le enlazan. Estas pruebas lo vigilan con tres controles:

1. **Estructural** (``ast``): ningún archivo del programa principal importa un
   módulo de usuario, de ninguna forma: ``import solidum.user.<x>``,
   ``from solidum.user.<x> import …``, ``from solidum.user import <x>``, sus
   formas relativas (``from ..user import x``, ``from ...user.x import y``),
   dentro de una función o a nivel de módulo, por cadena
   (``importlib.import_module("solidum.user." + n)``, ``__import__``,
   ``sys.modules["solidum.user.x"]``) ni con ``load_user_module("<nombre>")``
   y nombre literal. Lo único que el programa principal importa de
   ``solidum.user`` es el cargador genérico (``load_user_module``,
   ``available_user_modules``), que sólo construye rutas a partir del nombre
   que le pasa el modelo.
2. **Léxico**: ningún archivo del programa principal menciona lo que aporta un
   módulo de usuario. Los nombres vigilados no se mantienen a mano: se derivan
   cargando cada módulo en un proceso limpio (lo que añade a los registros, las
   familias que declara con su sección YAML y su ``SPEC_KIND``, las
   referencias nuevas de sus elementos, su ``__all__`` y su ruta de
   importación).
3. **En ejecución**: en un proceso limpio, ``import solidum`` no carga ningún
   módulo de usuario y toda clase registrada pertenece al programa principal.

Un módulo de usuario nuevo no tiene que darse de alta aquí, y todas las
pruebas pasan también si no hay ninguno.
"""
import ast
import json
import re
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SOLIDUM = ROOT / "solidum"

# Complemento opcional, a mano, para nombres característicos que la derivación
# no ve (constantes de nivel de módulo). Sólo se aplica a los módulos que
# existen; ningún módulo está obligado a figurar aquí.
NOMBRES_ADICIONALES = {
    "discontinuities": ["LOCAL_JUMP_RTOL", "LOCAL_JUMP_MAX_ITER"],
}


def _archivos_del_programa_principal():
    for p in sorted(SOLIDUM.rglob("*.py")):
        partes = p.relative_to(SOLIDUM).parts
        # solidum/user/__init__.py (la carga genérica) es del programa
        # principal; los subpaquetes de solidum/user/ son los módulos.
        if partes[0] == "user" and len(partes) > 2:
            continue
        yield p


# ---------------------------------------------------------------------------
# 1. Control estructural
# ---------------------------------------------------------------------------

_PAQUETE = "solidum.user"
# Única API de ``solidum.user`` que el programa principal puede importar.
_CARGADOR = {"load_user_module", "available_user_modules"}
# Funciones que importan un módulo a partir de su nombre.
_IMPORTADORES = {"import_module", "__import__", "_discover_package",
                 "find_spec", "run_module"}


def _nombre_de_modulo(p: Path) -> tuple:
    """``(ruta de importación, es paquete)`` del archivo ``p`` de ``solidum/``."""
    partes = ("solidum",) + p.relative_to(SOLIDUM).with_suffix("").parts
    if partes[-1] == "__init__":
        return ".".join(partes[:-1]), True
    return ".".join(partes), False


def _resolver(modulo: str, es_paquete: bool, nivel: int, destino) -> str:
    """Ruta absoluta de ``from <nivel puntos><destino> import …`` en ``modulo``."""
    base = modulo.split(".") if es_paquete else modulo.split(".")[:-1]
    if nivel > 1:
        base = base[: len(base) - (nivel - 1)]
    return ".".join(base + ([destino] if destino else []))


def _violaciones_estructurales(fuente: str, modulo: str, es_paquete: bool = False) -> list:
    """Formas de importar un módulo de usuario en ``fuente``, el código del
    módulo ``modulo``. Devuelve ``[(línea, descripción)]``."""
    arbol = ast.parse(fuente)
    # Cadenas sueltas (docstrings): pueden citar rutas sin importar nada.
    sueltas = {id(n.value) for n in ast.walk(arbol)
               if isinstance(n, ast.Expr) and isinstance(n.value, ast.Constant)}
    hallazgos = []
    for n in ast.walk(arbol):
        linea = getattr(n, "lineno", "?")
        if isinstance(n, ast.Import):
            for alias in n.names:
                if alias.name.startswith(_PAQUETE + "."):
                    hallazgos.append((linea, f"import {alias.name}"))
        elif isinstance(n, ast.ImportFrom):
            ruta = (_resolver(modulo, es_paquete, n.level, n.module) if n.level
                    else (n.module or ""))
            if ruta.startswith(_PAQUETE + "."):
                hallazgos.append((linea, f"from {ruta} import …"))
            elif ruta == _PAQUETE:
                ajenos = [a.name for a in n.names if a.name not in _CARGADOR]
                if ajenos:
                    hallazgos.append((linea, f"from {_PAQUETE} import {ajenos}"))
        elif (isinstance(n, ast.Constant) and isinstance(n.value, str)
              and id(n) not in sueltas and n.value.startswith(_PAQUETE + ".")):
            hallazgos.append((linea, f"cadena {n.value!r}"))
        elif isinstance(n, ast.Call):
            f = n.func
            nombre = f.attr if isinstance(f, ast.Attribute) else getattr(f, "id", "")
            if nombre in _IMPORTADORES:
                cadenas = [c.value for c in ast.walk(n)
                           if isinstance(c, ast.Constant) and isinstance(c.value, str)]
                if any(c.startswith(_PAQUETE) for c in cadenas):
                    hallazgos.append((linea, f"{nombre}(… '{_PAQUETE}…')"))
            if nombre == "load_user_module":
                argumentos = [*n.args, *(k.value for k in n.keywords)]
                if any(isinstance(a, (ast.Constant, ast.JoinedStr)) for a in argumentos):
                    hallazgos.append((linea, "load_user_module('<nombre literal>')"))
    return hallazgos


def test_el_programa_principal_no_importa_modulos_de_usuario():
    hallazgos = []
    for p in _archivos_del_programa_principal():
        modulo, es_paquete = _nombre_de_modulo(p)
        for linea, que in _violaciones_estructurales(
                p.read_text(encoding="utf-8"), modulo, es_paquete):
            hallazgos.append(f"{p.relative_to(ROOT)}:{linea}: {que}")
    assert not hallazgos, (
        "El programa principal importa módulos de usuario:\n" + "\n".join(hallazgos))


# Cada forma de escape se prueba contra el detector, para que la prueba de
# arriba no quede en verde por un detector ciego.
_ESCAPES = {
    "absoluta": "import solidum.user.discontinuities\n",
    "absoluta_alias": "import solidum.user.discontinuities as d\n",
    "from_submodulo": "from solidum.user.discontinuities import registry\n",
    "from_paquete": "from solidum.user import discontinuities\n",
    "from_asterisco": "from solidum.user import *\n",
    "relativa_paquete": "from ..user import discontinuities\n",
    "relativa_submodulo": "from ..user.discontinuities.registry import X\n",
    "perezosa": "def f():\n    from ..user import discontinuities as d\n    return d\n",
    "cadena_concatenada": ("import importlib\n"
                           "m = importlib.import_module('solidum.user.' + nombre)\n"),
    "cadena_fstring": "import importlib\nm = importlib.import_module(f'solidum.user.{n}')\n",
    "cadena_por_partes": ("import importlib\n"
                          "m = importlib.import_module('solidum.user' + '.' + n)\n"),
    "dunder_import": "m = __import__('solidum.user.discontinuities')\n",
    "sys_modules": "import sys\nm = sys.modules['solidum.user.discontinuities']\n",
    "carga_literal": ("from solidum.user import load_user_module\n"
                      "load_user_module('discontinuities')\n"),
    "carga_literal_kw": "import solidum\nsolidum.load_user_module(name='discontinuities')\n",
}


@pytest.mark.parametrize("fuente", list(_ESCAPES.values()), ids=list(_ESCAPES))
def test_el_detector_ve_cada_forma_de_importar(fuente):
    assert _violaciones_estructurales(fuente, "solidum.utils.ayuda")


def test_el_detector_resuelve_la_profundidad_de_las_relativas():
    assert _violaciones_estructurales(
        "from ...user.discontinuities import registry\n",
        "solidum.elements.solid_2d.ayuda")
    assert _violaciones_estructurales(
        "from ..user import tracking\n", "solidum.utils", es_paquete=True)
    # En un módulo de solidum/elements/solid_2d/, ``..user`` es
    # solidum.elements.user: no es la carpeta de módulos de usuario.
    assert not _violaciones_estructurales(
        "from ..user import tracking\n", "solidum.elements.solid_2d.ayuda")


_LEGITIMOS = {
    "cargador": "from solidum.user import available_user_modules, load_user_module\n",
    "cargador_relativo": "from ..user import load_user_module\n",
    "paquete_del_cargador": "import solidum.user\n",
    "carga_por_variable": "def f(nombre):\n    return load_user_module(nombre)\n",
    "docstring": '"""Los módulos de usuario viven en solidum.user.<nombre>."""\n',
    "otro_modulo_por_cadena": ("import importlib\n"
                               "importlib.import_module('solidum.materials.elastic')\n"),
}


@pytest.mark.parametrize("fuente", list(_LEGITIMOS.values()), ids=list(_LEGITIMOS))
def test_el_detector_admite_el_cargador_generico(fuente):
    assert _violaciones_estructurales(fuente, "solidum.utils.ayuda") == []


# ---------------------------------------------------------------------------
# 2. Control léxico con nombres derivados
# ---------------------------------------------------------------------------

_DERIVAR = r"""
import json
import sys

import solidum
import solidum.user
from solidum.registry import ElementRegistry, Registry

solidum.user.__path__.extend(sys.argv[1:])


def registros():
    vistos, pendientes = [], list(Registry.__subclasses__())
    while pendientes:
        reg = pendientes.pop()
        if reg not in vistos:
            vistos.append(reg)
            pendientes.extend(reg.__subclasses__())
    return vistos


def entradas():
    return {(reg, nombre): cls for reg in registros() for nombre, cls in reg._items.items()}


# Referencias que ya usa algún elemento del programa principal ('material').
claves_principales = {k for n in ElementRegistry.names()
                      for k in ElementRegistry.get(n).REFERENCE_KWARGS}
derivados = {}
for modulo in solidum.available_user_modules():
    antes, registros_antes = entradas(), set(registros())
    paquete = solidum.load_user_module(modulo)
    nombres = {f"solidum.user.{modulo}"} | set(getattr(paquete, "__all__", ()))
    for (reg, nombre), cls in entradas().items():
        if (reg, nombre) in antes:
            continue
        nombres |= {nombre, cls.__name__}
        nombres |= set(getattr(cls, "REFERENCE_KWARGS", {})) - claves_principales
    for reg in set(registros()) - registros_antes:
        nombres.add(reg.__name__)
        nombres |= {reg.__dict__.get(a) for a in ("YAML_SECTION", "SPEC_KIND")} - {None}
    derivados[modulo] = sorted(nombres)
print(json.dumps(derivados))
"""

# Módulo de juguete para comprobar la derivación: se añade a
# ``solidum.user.__path__`` sólo en el subproceso.
_JUGUETE = "juguete_pureza"
_ARCHIVOS_JUGUETE = {
    "__init__.py": textwrap.dedent('''\
        """Módulo de usuario de juguete para derivar sus nombres vigilados."""
        from solidum.user.juguete_pureza.resortes import ResorteLineal

        AYUDA_DE_JUGUETE = 1
        __all__ = ["ResorteLineal", "AYUDA_DE_JUGUETE"]
        '''),
    "resortes.py": textwrap.dedent('''\
        from solidum.elements.truss import Truss2D
        from solidum.registry import ElementRegistry, MaterialRegistry, Registry


        class ResorteRegistry(Registry):
            SPEC_KIND = "resorte_pureza"
            YAML_SECTION = "resortes_pureza"


        @ResorteRegistry.register
        class ResorteLineal:
            def __init__(self, k: float):
                self.k = k


        @ElementRegistry.register("BarraConResorteDeJuguete")
        class BarraConResorte(Truss2D):
            REFERENCE_KWARGS = {"material": MaterialRegistry,
                                "resorte_de_juguete": ResorteRegistry}
        '''),
}


@pytest.fixture(scope="module")
def nombres_derivados(tmp_path_factory):
    """Nombres vigilados de cada módulo de usuario instalado, más el de juguete."""
    raiz = tmp_path_factory.mktemp("modulos_pureza")
    (raiz / _JUGUETE).mkdir()
    for nombre, texto in _ARCHIVOS_JUGUETE.items():
        (raiz / _JUGUETE / nombre).write_text(texto, encoding="utf-8")
    r = subprocess.run([sys.executable, "-c", _DERIVAR, str(raiz)], cwd=ROOT,
                       capture_output=True, text=True, encoding="utf-8", timeout=300)
    assert r.returncode == 0, r.stdout + r.stderr
    return json.loads(r.stdout.strip().splitlines()[-1])


def test_la_derivacion_recoge_lo_que_aporta_el_modulo(nombres_derivados):
    nombres = set(nombres_derivados[_JUGUETE])
    assert {
        "solidum.user.juguete_pureza",          # ruta de importación
        "ResorteRegistry", "resortes_pureza", "resorte_pureza",   # familia
        "ResorteLineal",                        # clase registrada en la familia
        "BarraConResorteDeJuguete", "BarraConResorte",   # alias y clase
        "resorte_de_juguete",                   # referencia nueva
        "AYUDA_DE_JUGUETE",                     # __all__
    } <= nombres
    # Lo que el módulo reutiliza del programa principal no se vigila.
    assert not {"material", "materials", "Truss2D", "MaterialRegistry"} & nombres


def test_se_derivan_los_nombres_de_todos_los_modulos(nombres_derivados):
    """Sustituye a la antigua lista a mano: ningún módulo queda sin vigilar."""
    import solidum
    assert set(nombres_derivados) - {_JUGUETE} == set(solidum.available_user_modules())


def test_el_programa_principal_no_nombra_a_los_modulos(nombres_derivados):
    hallazgos = []
    for modulo in sorted(set(nombres_derivados) - {_JUGUETE}):
        nombres = set(nombres_derivados[modulo]) | set(NOMBRES_ADICIONALES.get(modulo, ()))
        patron = re.compile("|".join(
            r"\b" + re.escape(x) + r"\b" for x in sorted(nombres, key=len, reverse=True)))
        for p in _archivos_del_programa_principal():
            for n, linea in enumerate(p.read_text(encoding="utf-8").splitlines(), 1):
                if patron.search(linea):
                    hallazgos.append(f"[{modulo}] {p.relative_to(ROOT)}:{n}: {linea.strip()}")
    assert not hallazgos, (
        "El programa principal nombra módulos de usuario:\n" + "\n".join(hallazgos))


# ---------------------------------------------------------------------------
# 3. Control en ejecución
# ---------------------------------------------------------------------------

_SCRIPT = """
import sys
import solidum
from solidum.registry import Registry
cargados = sorted(m for m in sys.modules if m.startswith("solidum.user."))
ajenos = []
registros = [r for r in Registry.__subclasses__()]
for reg in registros:
    for nombre in reg.names():
        mod = reg.get(nombre).__module__
        if mod.startswith("solidum.user."):
            ajenos.append(f"{reg.__name__}:{nombre}:{mod}")
secciones = [f.YAML_SECTION for f in Registry.families()]
print(repr((cargados, ajenos, secciones)))
"""


def test_import_solidum_no_carga_modulos_de_usuario():
    r = subprocess.run([sys.executable, "-c", _SCRIPT], cwd=ROOT,
                       capture_output=True, text=True, encoding="utf-8", timeout=300)
    assert r.returncode == 0, r.stderr
    cargados, ajenos, secciones = eval(r.stdout.strip().splitlines()[-1])
    assert cargados == [], f"import solidum cargó módulos de usuario: {cargados}"
    assert ajenos == [], f"Clases de módulos de usuario registradas: {ajenos}"
    assert secciones == ["materials", "thermal_materials"]
