"""Specs del módulo de usuario de discontinuidades (ADR 0020).

Las specs del módulo viven en ``docs/user/discontinuities/specs/``. Se validan
con la misma herramienta que las del programa principal —esquema y contrato
declarativo frente a la clase registrada—, con el módulo cargado, y además
con la interfaz propia de la familia cohesiva (``JUMP_DIM``,
``PRIMARY_STATE_VAR``, ``IS_SYMMETRIC``), que el programa principal no
conoce.
"""
from pathlib import Path

import pytest

import solidum
from solidum.tools.spec import (
    collect_specs,
    cross_check_with_registry,
    parse_spec,
    validate_schema,
)

solidum.load_user_module("discontinuities")
from solidum.user.discontinuities import CohesiveMaterialRegistry  # noqa: E402

SPECS_DIR = Path(__file__).resolve().parents[3] / "docs" / "user" / "discontinuities" / "specs"
SPEC_PATHS = collect_specs(SPECS_DIR)

# Interfaz propia de la familia cohesiva: atributo de clase ↔ clave de la spec.
_INTERFAZ_COHESIVA = (
    ("JUMP_DIM", "jump_dim"),
    ("PRIMARY_STATE_VAR", "primary_state_var"),
    ("IS_SYMMETRIC", "is_symmetric"),
)


def test_hay_specs():
    assert {p.stem for p in SPEC_PATHS} == {"CST_Embedded2D", "CohesiveDamageIsotropic"}


@pytest.mark.parametrize("spec_path", SPEC_PATHS, ids=[p.name for p in SPEC_PATHS])
def test_spec_schema_y_contrato(spec_path):
    spec = parse_spec(spec_path)
    validate_schema(spec)
    assert spec.name == spec_path.stem
    errors = cross_check_with_registry(spec)
    assert not errors, "\n".join(errors)


@pytest.mark.parametrize("spec_path", [p for p in SPEC_PATHS
                                       if parse_spec(p).kind == "cohesive_material"],
                         ids=lambda p: p.name)
def test_interfaz_cohesiva(spec_path):
    spec = parse_spec(spec_path)
    cls = CohesiveMaterialRegistry.get(spec.name)
    iface = spec.contract["interface"]
    for attr, key in _INTERFAZ_COHESIVA:
        if key not in iface:
            continue                  # opcional en la spec; lo declarado debe coincidir
        assert hasattr(cls, attr), f"{cls.__name__} no declara {attr}"
        assert getattr(cls, attr) == iface[key], (
            f"{attr}={getattr(cls, attr)!r} (código) vs {key}={iface[key]!r} (spec)")
