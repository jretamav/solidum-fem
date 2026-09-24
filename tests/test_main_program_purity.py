"""Pureza del programa principal frente a los módulos de usuario (ADR 0020).

El programa principal es Solidum estándar: todo ``solidum/`` salvo los
módulos de usuario ``solidum/user/<nombre>/``. No los importa ni los nombra,
como el programa principal de FEAP no conoce los elementos de usuario que se
le enlazan. Estas pruebas lo vigilan:

1. Ningún archivo del programa principal menciona los nombres de un módulo de
   usuario (clases, constantes, claves YAML, ruta de importación).
2. En un proceso limpio, ``import solidum`` no carga ningún módulo de usuario
   y toda clase registrada pertenece al programa principal.

Un módulo de usuario nuevo añade aquí sus nombres característicos.
"""
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SOLIDUM = ROOT / "solidum"

# Nombres característicos de cada módulo de usuario.
NOMBRES_POR_MODULO = {
    "discontinuities": [
        r"CST_Embedded2D", r"CohesiveMaterial", r"CohesiveDamageIsotropic",
        r"DiscontinuityState", r"EMBEDDED_LOCAL_JUMP", r"cohesive_material",
        r"solidum\.user\.discontinuities", r"(?i:embebid|embedded)",
    ],
}


def _archivos_del_programa_principal():
    for p in sorted(SOLIDUM.rglob("*.py")):
        partes = p.relative_to(SOLIDUM).parts
        # solidum/user/__init__.py (la carga genérica) es del programa
        # principal; los subpaquetes de solidum/user/ son los módulos.
        if partes[0] == "user" and len(partes) > 2:
            continue
        yield p


@pytest.mark.parametrize("modulo", sorted(NOMBRES_POR_MODULO))
def test_el_programa_principal_no_nombra_al_modulo(modulo):
    patron = re.compile("|".join(NOMBRES_POR_MODULO[modulo]))
    hallazgos = []
    for p in _archivos_del_programa_principal():
        for n, linea in enumerate(p.read_text(encoding="utf-8").splitlines(), 1):
            if patron.search(linea):
                hallazgos.append(f"{p.relative_to(ROOT)}:{n}: {linea.strip()}")
    assert not hallazgos, (
        f"El programa principal nombra el módulo de usuario '{modulo}':\n"
        + "\n".join(hallazgos))


def test_los_modulos_de_usuario_existen():
    """Cada módulo vigilado sigue existiendo (si se retira, se quita de aquí)."""
    import solidum
    assert set(NOMBRES_POR_MODULO) <= set(solidum.available_user_modules())


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
