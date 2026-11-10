"""Specs del módulo de usuario de discontinuidades (ADR 0020).

Las specs del módulo viven en ``docs/user/discontinuities/specs/``. Se validan
con la misma herramienta que las del programa principal —esquema y contrato
declarativo frente a la clase registrada—, con el módulo cargado. La ley
cohesiva es un material de interfaz del programa principal (ADR 0033) y su
spec se valida con las de él.
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

SPECS_DIR = Path(__file__).resolve().parents[3] / "docs" / "user" / "discontinuities" / "specs"
SPEC_PATHS = collect_specs(SPECS_DIR)


def test_hay_specs():
    assert {p.stem for p in SPEC_PATHS} == {"CST_Embedded2D", "GlobalCrackTracking",
                                            "PrescribedCrackPath"}


@pytest.mark.parametrize("spec_path", SPEC_PATHS, ids=[p.name for p in SPEC_PATHS])
def test_spec_schema_y_contrato(spec_path):
    spec = parse_spec(spec_path)
    validate_schema(spec)
    assert spec.name == spec_path.stem
    errors = cross_check_with_registry(spec)
    assert not errors, "\n".join(errors)
