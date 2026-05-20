"""Tests de coherencia de las specs de componentes.

Cada spec en docs/specs/ debe:
 1. Parsear y validar schema.
 2. Si status ∈ {implemented, validated}: la clase debe estar registrada y
    su contrato declarativo debe coincidir con el YAML de la spec.
"""
from pathlib import Path

import pytest

from solidum.tools.spec import (
    collect_specs,
    cross_check_with_registry,
    parse_spec,
    validate_schema,
)

SPECS_DIR = Path(__file__).resolve().parents[1] / "docs" / "specs"


def _spec_ids(paths):
    return [p.name for p in paths]


SPEC_PATHS = collect_specs(SPECS_DIR)


@pytest.mark.parametrize("spec_path", SPEC_PATHS, ids=_spec_ids(SPEC_PATHS))
def test_spec_schema(spec_path):
    spec = parse_spec(spec_path)
    validate_schema(spec)
    assert spec.name == spec_path.stem, (
        f"spec.name='{spec.name}' debe coincidir con nombre de archivo "
        f"'{spec_path.stem}'"
    )


@pytest.mark.parametrize("spec_path", SPEC_PATHS, ids=_spec_ids(SPEC_PATHS))
def test_spec_matches_code(spec_path):
    spec = parse_spec(spec_path)
    validate_schema(spec)
    errors = cross_check_with_registry(spec)
    assert not errors, "\n".join(errors)
