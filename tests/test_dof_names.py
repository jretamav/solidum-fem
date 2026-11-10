"""Nombres de los grados de libertad: convención única (``solidum.core.dof_names``).

Criterios: las constantes son las del programa (los elementos térmicos usan
``TEMPERATURE``); la dimensión de las traslaciones sale de ``uz``; y ningún
módulo del programa principal indexa los grados de libertad con el nombre
literal de la temperatura ni repite la regla de ``uz``, para que un cambio
de convención se haga en un solo sitio.
"""
import re
from pathlib import Path

from solidum.core.dof_names import (
    MECHANICAL,
    MECHANICAL_ORDER,
    ROTATIONS,
    TEMPERATURE,
    TRANSLATIONS,
    translation_dims,
)
from solidum.core.node import Node
from solidum.elements.thermal import Hex8Thermal, Quad4Thermal

PACKAGE = Path(__file__).resolve().parents[1] / "solidum"


def test_constantes():
    assert MECHANICAL_ORDER == ("ux", "uy", "uz", "rx", "ry", "rz")
    assert MECHANICAL == frozenset(TRANSLATIONS + ROTATIONS)
    assert TEMPERATURE == "T"
    assert Quad4Thermal.DOF_NAMES == [TEMPERATURE]
    assert Hex8Thermal.DOF_NAMES == [TEMPERATURE]


def test_dimension_de_las_traslaciones():
    plano, espacio = Node(1, [0.0, 0.0]), Node(2, [0.0, 0.0, 0.0])
    for name in ("ux", "uy"):
        plano.add_dof(name)
        espacio.add_dof(name)
    espacio.add_dof("uz")
    assert translation_dims([plano]) == 2
    assert translation_dims([plano, espacio]) == 3


def test_sin_literales_fuera_de_la_convencion():
    patterns = [re.compile(r"""\.dofs\[\s*["']T["']\s*\]"""),
                re.compile(r"""["']uz["']\s+in\s+\w+\.dofs""")]
    offenders = []
    for path in PACKAGE.rglob("*.py"):
        rel = path.relative_to(PACKAGE).as_posix()
        if rel == "core/dof_names.py" or rel.startswith("user/"):
            continue
        text = path.read_text(encoding="utf-8")
        for pattern in patterns:
            for match in pattern.finditer(text):
                line = text.count("\n", 0, match.start()) + 1
                offenders.append(f"{rel}:{line}: {match.group(0)}")
    assert not offenders, "Usar solidum.core.dof_names:\n" + "\n".join(offenders)
