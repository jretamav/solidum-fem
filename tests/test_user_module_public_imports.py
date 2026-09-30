"""Los módulos de usuario sólo usan el contrato público del programa principal.

Contraparte de ``tests/test_main_program_purity.py``, que vigila la otra
dirección (el programa principal no conoce los módulos). ADR 0020: un módulo
de usuario se conecta al programa principal por piezas públicas (P1-P8). Si se
apoya en un detalle interno —un módulo ``_shared``, un atributo ``_kind``—, el
programa principal ya no puede reorganizarlo sin romper el módulo, y aparece
la presión para congelar internos o para meter casos especiales.

Se comprueba, en cada módulo de ``solidum/user/<nombre>/``:

1. **Importaciones** (análisis sintáctico con ``ast``): ninguna importación del
   programa principal pasa por una ruta con un componente privado
   (``solidum.elements.solid_2d._shared``) ni trae un nombre con guion bajo
   (``from solidum.autodiscover import _discover_package``). Las importaciones
   entre archivos del propio módulo quedan libres. Lo mismo se exige a sus
   tests, ``tests/user/<nombre>/``.
2. **Registros**: una subclase de ``solidum.registry.Registry`` definida en el
   módulo declara sólo atributos públicos (``SPEC_KIND``, ``YAML_SECTION``,
   ``YAML_LABEL``). El almacén ``_items`` lo crea la base. Las subclases se
   localizan en tiempo de ejecución, tras cargar el módulo, no por el nombre
   de sus bases.

``solidum/user/__init__.py`` es el cargador (``load_user_module``) y pertenece
al programa principal: no se analiza.
"""
import ast
import inspect
import textwrap
import unittest
from pathlib import Path

import solidum
from solidum.registry import Registry
from solidum.user import available_user_modules

ROOT = Path(__file__).resolve().parents[1]


def es_programa_principal(modulo: str) -> bool:
    """``solidum`` y sus subpaquetes, salvo los módulos de usuario
    (``solidum.user.<nombre>``). El cargador ``solidum.user`` es del programa
    principal."""
    if modulo != "solidum" and not modulo.startswith("solidum."):
        return False
    return not modulo.startswith("solidum.user.")


def _absoluto(nodo: ast.ImportFrom, paquete: str) -> str:
    """Ruta absoluta de un ``from … import`` resuelto desde ``paquete``."""
    if nodo.level == 0:
        return nodo.module or ""
    partes = paquete.split(".")
    base = partes[:len(partes) - (nodo.level - 1)]
    return ".".join(base + ([nodo.module] if nodo.module else []))


def importaciones_privadas(fuente: str, paquete: str) -> list:
    """Importaciones del programa principal que no usan su contrato público.

    ``paquete`` es el paquete al que pertenece el archivo (para resolver las
    importaciones relativas). Devuelve ``[(línea, descripción)]``.
    """
    hallazgos = []
    for nodo in ast.walk(ast.parse(fuente)):
        if isinstance(nodo, ast.ImportFrom):
            modulo = _absoluto(nodo, paquete)
            if not es_programa_principal(modulo):
                continue
            ruta = [c for c in modulo.split(".") if c.startswith("_")]
            nombres = [a.name for a in nodo.names if a.name.startswith("_")]
            if ruta or nombres:
                hallazgos.append((nodo.lineno, f"from {modulo} import "
                                               f"{', '.join(a.name for a in nodo.names)}"))
        elif isinstance(nodo, ast.Import):
            for alias in nodo.names:
                if (es_programa_principal(alias.name)
                        and any(c.startswith("_") for c in alias.name.split("."))):
                    hallazgos.append((nodo.lineno, f"import {alias.name}"))
    return hallazgos


def atributos_privados_declarados(clase) -> list:
    """Atributos con guion bajo (no dunder) asignados en el cuerpo de
    ``clase``."""
    arbol = ast.parse(textwrap.dedent(inspect.getsource(clase)))
    definicion = next(n for n in arbol.body if isinstance(n, ast.ClassDef))
    privados = []
    for sentencia in definicion.body:
        if isinstance(sentencia, ast.Assign):
            objetivos = sentencia.targets
        elif isinstance(sentencia, ast.AnnAssign):
            objetivos = [sentencia.target]
        else:
            continue
        privados += [o.id for o in objetivos
                     if isinstance(o, ast.Name) and o.id.startswith("_")
                     and not o.id.startswith("__")]
    return privados


def _subclases(clase):
    for sub in clase.__subclasses__():
        yield sub
        yield from _subclases(sub)


def _archivos(carpeta: Path) -> list:
    return sorted(p for p in carpeta.rglob("*.py") if "__pycache__" not in p.parts)


def _paquete_de(archivo: Path) -> str:
    """Paquete al que pertenece un archivo de ``solidum/``."""
    partes = archivo.relative_to(ROOT).with_suffix("").parts
    return ".".join(partes if archivo.name == "__init__.py" else partes[:-1])


class TestDetector(unittest.TestCase):
    """El detector encuentra lo que dice encontrar: sin esto, los tests de
    abajo pasarían por vacuidad si el análisis dejara de ver algo."""

    def test_importaciones(self):
        fuente = textwrap.dedent("""\
            from solidum.elements.solid_2d._shared import compute_integrands
            from solidum.autodiscover import _discover_package
            import solidum.math.solvers._shared
            from ...elements.solid_2d._shared import compute_kinematics_tri3
            from solidum.elements.solid_2d import compute_integrands
            from solidum.registry import Registry
            from solidum.user.otro_modulo._interno import algo
            from ._interno import _propio
            """)
        lineas = [linea for linea, _ in
                  importaciones_privadas(fuente, "solidum.user.ejemplo")]
        self.assertEqual(lineas, [1, 2, 3, 4])

    def test_atributos_privados_en_un_registro(self):
        class RegistroMal(Registry):
            _kind = "Resorte"
            _items: dict = {}
            YAML_LABEL = "resorte"

        class RegistroBien(Registry):
            YAML_LABEL = "resorte"

        self.assertEqual(atributos_privados_declarados(RegistroMal), ["_kind", "_items"])
        self.assertEqual(atributos_privados_declarados(RegistroBien), [])


class TestModulosDeUsuario(unittest.TestCase):

    def setUp(self):
        self.modulos = available_user_modules()
        if not self.modulos:
            self.skipTest("no hay módulos de usuario instalados")

    def test_sus_importaciones_del_programa_principal_son_publicas(self):
        for nombre in self.modulos:
            carpetas = [ROOT / "solidum" / "user" / nombre, ROOT / "tests" / "user" / nombre]
            archivos = [a for c in carpetas if c.is_dir() for a in _archivos(c)]
            self.assertTrue(archivos, f"módulo '{nombre}' sin archivos que analizar")
            for archivo in archivos:
                # Los tests no son un paquete: sus importaciones relativas no
                # llegan al programa principal.
                paquete = (_paquete_de(archivo) if archivo.is_relative_to(ROOT / "solidum")
                           else "tests")
                hallazgos = importaciones_privadas(archivo.read_text(encoding="utf-8"), paquete)
                with self.subTest(archivo=str(archivo.relative_to(ROOT))):
                    self.assertFalse(
                        hallazgos,
                        "importa del programa principal fuera de su contrato público "
                        "(usar la ruta pública; si no existe, exponerla en el "
                        "programa principal como pieza genérica):\n"
                        + "\n".join(f"  línea {n}: {texto}" for n, texto in hallazgos))

    def test_sus_registros_declaran_solo_atributos_publicos(self):
        for nombre in self.modulos:
            solidum.load_user_module(nombre)
            paquete = f"solidum.user.{nombre}"
            registros = [c for c in _subclases(Registry)
                         if (c.__module__ + ".").startswith(paquete + ".")]
            for registro in registros:
                with self.subTest(registro=f"{registro.__module__}.{registro.__qualname__}"):
                    self.assertEqual(
                        atributos_privados_declarados(registro), [],
                        "una subclase de Registry declara sólo YAML_LABEL, "
                        "SPEC_KIND y YAML_SECTION; _items lo crea la base")


if __name__ == "__main__":
    unittest.main()
