"""Tests de coherencia de las specs de componentes.

Cada spec en docs/specs/ debe:
 1. Parsear y validar schema, incluido el bloque ``admission``: el criterio de
    admisión al programa principal (Reglas.md §1) con sus referencias en las
    listas cerradas de textos y de códigos generales.
 2. Si status ∈ {implemented, validated}: la clase debe estar registrada y
    su contrato declarativo debe coincidir con el YAML de la spec.

Y en sentido inverso, del código a la spec: todo componente registrado del
programa principal tiene su spec en docs/specs/ con ``admission`` válido. Así
ninguna formulación entra al programa principal sin declarar, de forma
revisable, por qué es estándar.
"""
import re
from pathlib import Path

import pytest

from solidum.tools.spec import (
    ADMISSION_GENERAL_CODES,
    ADMISSION_TEXTBOOKS,
    SpecError,
    admission_entries,
    collect_specs,
    cross_check_with_registry,
    parse_spec,
    validate_schema,
)

REPO = Path(__file__).resolve().parents[1]
SPECS_DIR = REPO / "docs" / "specs"
REGLAS = REPO / "Reglas.md"


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


# ---------------------------------------------------------------------------
# Del código a la spec: todo componente registrado del programa principal
# declara su admisión.
# ---------------------------------------------------------------------------

def _main_program_components():
    """``(kind, nombre)`` de cada componente registrado del programa principal.

    Se filtra por el módulo que define la clase (y el registro), no por qué
    módulos están cargados: otro test de la misma sesión puede haber cargado
    un módulo de usuario, cuyas clases y registros propios viven en
    ``solidum.user.*`` (ADR 0020).
    """
    import solidum  # noqa: F401 — dispara autodiscover
    from solidum.registry import Registry

    components = []
    for kind in Registry.spec_kinds():
        registry = Registry.for_spec_kind(kind)
        if registry.__module__.startswith("solidum.user."):
            continue
        for name in registry.names():
            if registry.get(name).__module__.startswith("solidum.user."):
                continue
            components.append((kind, name))
    return components


MAIN_COMPONENTS = _main_program_components()


def test_hay_componentes_del_programa_principal():
    assert MAIN_COMPONENTS, "ningún componente registrado tras import solidum"


@pytest.mark.parametrize("kind,name", MAIN_COMPONENTS,
                         ids=[name for _, name in MAIN_COMPONENTS])
def test_componente_registrado_tiene_spec_con_admission(kind, name):
    path = SPECS_DIR / f"{name}.md"
    assert path.is_file(), (
        f"{name} está registrado en el programa principal pero no tiene spec en "
        f"docs/specs/{name}.md. Todo componente estándar entra con su spec y su "
        f"bloque 'admission' (Reglas.md §1 y §4); una formulación que no cumpla "
        f"el criterio va a un módulo de usuario (ADR 0020)."
    )
    spec = parse_spec(path)
    assert spec.kind == kind, (
        f"docs/specs/{name}.md declara kind={spec.kind!r}, pero la clase está "
        f"registrada como {kind!r}"
    )
    validate_schema(spec, main_program=True)


# ---------------------------------------------------------------------------
# Las listas cerradas del código son las de Reglas.md §1.
# ---------------------------------------------------------------------------

def _lista_de_reglas(pattern: str) -> list:
    text = REGLAS.read_text(encoding="utf-8")
    match = re.search(pattern, text)
    assert match, (
        f"no se encontró en Reglas.md §1 la enumeración que busca {pattern!r}; "
        "si cambió la redacción del criterio de admisión, actualiza este test"
    )
    return sorted(item.strip() for item in match.group(1).split(","))


def test_listas_cerradas_coinciden_con_reglas():
    textos = _lista_de_reglas(r"textos de referencia del área \(([^)]*)\)")
    codigos = _lista_de_reglas(r"código general establecido \(([^)]*)\)")
    assert textos == sorted(ADMISSION_TEXTBOOKS.values()), (
        "ADMISSION_TEXTBOOKS (solidum/tools/spec.py) no coincide con Reglas.md §1")
    assert codigos == sorted(ADMISSION_GENERAL_CODES.values()), (
        "ADMISSION_GENERAL_CODES (solidum/tools/spec.py) no coincide con Reglas.md §1")


TEMPLATE_PATHS = sorted(SPECS_DIR.glob("_template_*.md"))


@pytest.mark.parametrize("template", TEMPLATE_PATHS, ids=_spec_ids(TEMPLATE_PATHS))
def test_plantilla_declara_admission(template):
    admission = parse_spec(template).contract.get("admission") or {}
    assert set(admission) == {"textbooks", "general_codes"}


# ---------------------------------------------------------------------------
# Validador del bloque admission.
# ---------------------------------------------------------------------------

_VALID_ADMISSION = (
    "admission:\n"
    "  textbooks:\n"
    "    - hughes\n"
    '    - bathe: "§9.2"\n'
    "  general_codes:\n"
    "    - opensees: \"integrator Newmark\"\n"
)


def _write_spec(directory: Path, admission: str):
    directory.mkdir(parents=True)
    path = directory / "X.md"
    path.write_text(
        "# X\n\n```yaml\nname: X\nkind: solver\nstatus: draft\n"
        "interface: {}\n" + admission + "```\n",
        encoding="utf-8",
    )
    return parse_spec(path)


def test_admission_valida_normalizada(tmp_path):
    spec = _write_spec(tmp_path / "docs" / "specs", _VALID_ADMISSION)
    validate_schema(spec)
    assert admission_entries(spec) == {
        "textbooks": [("hughes", None), ("bathe", "§9.2")],
        "general_codes": [("opensees", "integrator Newmark")],
    }


@pytest.mark.parametrize("admission", ["", "admission:\n"], ids=["ausente", "vacio"])
def test_admission_obligatoria_en_el_programa_principal(tmp_path, admission):
    spec = _write_spec(tmp_path / "docs" / "specs", admission)
    with pytest.raises(SpecError, match="falta el bloque 'admission'"):
        validate_schema(spec)


def test_admission_opcional_en_modulo_de_usuario(tmp_path):
    spec = _write_spec(tmp_path / "docs" / "user" / "mod" / "specs", "")
    validate_schema(spec)
    with pytest.raises(SpecError, match="admission"):
        validate_schema(spec, main_program=True)


@pytest.mark.parametrize("admission,match", [
    ("admission:\n  textbooks: [hughes]\n  general_codes: [nastran]\n", "lista cerrada"),
    ("admission:\n  textbooks: [simo_hughes]\n  general_codes: [abaqus]\n", "lista cerrada"),
    ("admission:\n  textbooks: []\n  general_codes: [abaqus]\n", "al menos una"),
    ("admission:\n  textbooks: [hughes]\n", "general_codes"),
    ("admission:\n  textbooks: [hughes, hughes]\n  general_codes: [abaqus]\n", "repetido"),
    ("admission:\n  textbooks: [hughes]\n  general_codes: [abaqus]\n  otros: [x]\n",
     "desconocidas"),
    ("admission:\n  textbooks:\n    - {hughes: a, bathe: b}\n  general_codes: [abaqus]\n",
     "inválida"),
    ("admission: [hughes, abaqus]\n", "mapping"),
])
def test_admission_invalida(tmp_path, admission, match):
    spec = _write_spec(tmp_path / "docs" / "specs", admission)
    with pytest.raises(SpecError, match=match):
        validate_schema(spec)


def test_admission_invalida_tambien_en_modulo_de_usuario(tmp_path):
    """En un módulo de usuario el bloque es opcional, pero si se declara vale
    lo mismo que en el programa principal."""
    spec = _write_spec(tmp_path / "docs" / "user" / "mod" / "specs",
                       "admission:\n  textbooks: [hughes]\n  general_codes: [nastran]\n")
    with pytest.raises(SpecError, match="lista cerrada"):
        validate_schema(spec)
