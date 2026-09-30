"""Specs del módulo de usuario ``dissipation_control`` (ADR 0020).

Las specs del módulo viven en ``docs/user/dissipation_control/specs/``. Se
validan con la misma herramienta que las del programa principal —esquema y
contrato declarativo frente a la clase registrada—, con el módulo cargado.
Además, como la spec del solver es una spec corta tipo extensión (Reglas.md
§4), se comprueba que su padre existe en el programa principal y que la
clase hereda de él.
"""
from pathlib import Path

import pytest

import solidum
from solidum.registry import Registry
from solidum.tools.spec import (
    collect_specs,
    cross_check_with_registry,
    parse_spec,
    validate_schema,
)

solidum.load_user_module("dissipation_control")

ROOT = Path(__file__).resolve().parents[3]
SPECS_DIR = ROOT / "docs" / "user" / "dissipation_control" / "specs"
MAIN_SPECS_DIR = ROOT / "docs" / "specs"
SPEC_PATHS = collect_specs(SPECS_DIR)


def test_hay_specs():
    assert {p.stem for p in SPEC_PATHS} == {"DissipationArcLengthSolver"}


@pytest.mark.parametrize("spec_path", SPEC_PATHS, ids=[p.name for p in SPEC_PATHS])
def test_spec_schema_y_contrato(spec_path):
    spec = parse_spec(spec_path)
    validate_schema(spec)
    assert spec.name == spec_path.stem
    errors = cross_check_with_registry(spec)
    assert not errors, "\n".join(errors)


@pytest.mark.parametrize("spec_path", SPEC_PATHS, ids=[p.name for p in SPEC_PATHS])
def test_spec_padre_en_el_programa_principal(spec_path):
    spec = parse_spec(spec_path)
    padre = spec.contract.get("parent_spec")
    if padre is None:
        pytest.skip("spec sin padre")
    assert (MAIN_SPECS_DIR / f"{padre}.md").is_file(), (
        f"{spec.name}: parent_spec={padre!r} no tiene spec en docs/specs/")
    registro = Registry.for_spec_kind(spec.kind)
    assert issubclass(registro.get(spec.name), registro.get(padre))
