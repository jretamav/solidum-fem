"""Cobertura del manual de referencia: ninguna spec puede quedarse fuera.

Motivación (saneamiento 2026-08-25): el builder del *Reference manual*
enumeraba a mano las specs que entraban al PDF. Esa lista se quedó atrás
respecto al repositorio — 23 de 46 specs no llegaban al manual (todo el
subsistema 3D, la fractura embebida y 8 solvers) — **sin producir ningún
error**, mientras la portada del PDF afirmaba estar generada automáticamente
desde ``docs/specs/``.

El defecto real no era la lista desactualizada sino que fallara en silencio.
Estos tests convierten el invariante en algo que la suite protege: si un
componente futuro entra al repositorio con spec pero no llega al manual, la
suite se pone roja en vez de emitir un PDF incompleto.

Desde el ADR 0020 hay dos procedencias. Las specs del programa principal
(``docs/specs/``) se reparten en los capítulos del cuerpo; las de los módulos
de usuario (``docs/user/<módulo>/specs/``) van al apéndice «Módulos de
usuario», agrupadas por módulo. El invariante vale para las dos.
"""
import shutil
import sys
from pathlib import Path

import pytest

from solidum.tools.spec import collect_specs, parse_spec

ROOT = Path(__file__).resolve().parents[1]
SPECS_DIR = ROOT / "docs" / "specs"
USER_DOCS_DIR = ROOT / "docs" / "user"

sys.path.insert(0, str(ROOT / "manuals"))
build_reference_manual = pytest.importorskip(
    "build_reference_manual",
    reason="builder del manual de referencia no importable",
)


def _grouped_names() -> list[str]:
    names: list[str] = []
    for _chapter, components in build_reference_manual.build_groups():
        names.extend(components)
    return names


def test_toda_spec_entra_en_algun_capitulo():
    """Cada spec de ``docs/specs/`` aparece exactamente una vez en el manual."""
    grouped = _grouped_names()
    expected = {parse_spec(p).name for p in collect_specs(SPECS_DIR)}

    faltantes = expected - set(grouped)
    assert not faltantes, (
        "Specs que no llegarían al Reference manual: "
        f"{sorted(faltantes)}. Revisar `_classify` en "
        "manuals/build_reference_manual.py."
    )

    duplicadas = [n for n in grouped if grouped.count(n) > 1]
    assert not duplicadas, f"Specs duplicadas entre capítulos: {sorted(set(duplicadas))}"


def test_capitulos_declarados_y_no_vacios():
    """Todo capítulo emitido está en ``CHAPTER_ORDER`` y lleva contenido."""
    groups = build_reference_manual.build_groups()
    orden = build_reference_manual.CHAPTER_ORDER

    for chapter, components in groups:
        assert chapter in orden, f"Capítulo '{chapter}' fuera de CHAPTER_ORDER"
        assert components, f"Capítulo '{chapter}' emitido vacío"

    emitidos = [ch for ch, _ in groups]
    assert emitidos == sorted(emitidos, key=orden.index), (
        "El orden de los capítulos no respeta CHAPTER_ORDER"
    )


def test_spec_no_clasificable_aborta_el_build(tmp_path, monkeypatch):
    """Una spec no clasificable detiene el build en vez de omitirse callando.

    Es la garantía central del saneamiento: el modo de fallo anterior era la
    omisión silenciosa.
    """
    specs = collect_specs(SPECS_DIR)
    fuente = next(p for p in specs if parse_spec(p).name == "Hex20")

    sandbox = tmp_path / "specs"
    sandbox.mkdir()
    texto = fuente.read_text(encoding="utf-8")
    (sandbox / "Hex20.md").write_text(
        texto.replace("strain_dim: 6", "strain_dim: 99", 1), encoding="utf-8"
    )

    monkeypatch.setattr(build_reference_manual, "SPECS_DIR", sandbox)
    with pytest.raises(SystemExit) as exc:
        build_reference_manual.build_groups()

    assert "Hex20" in str(exc.value), (
        "El aborto debe nombrar la spec culpable para que sea accionable"
    )


def _user_specs_on_disk() -> set[tuple[str, str]]:
    """(módulo, nombre) de toda spec de ``docs/user/*/specs/``."""
    return {(d.parent.name, parse_spec(p).name)
            for d in USER_DOCS_DIR.glob("*/specs") if d.is_dir()
            for p in collect_specs(d)}


def test_toda_spec_de_modulo_de_usuario_entra_en_el_apendice():
    """Cada spec de ``docs/user/*/specs/`` aparece exactamente una vez en el
    apéndice «Módulos de usuario», bajo su módulo, y ninguna en el cuerpo."""
    esperadas = _user_specs_on_disk()
    assert esperadas, "No hay specs de módulos de usuario: ¿se movió docs/user/?"

    grupos = build_reference_manual.build_user_groups()
    obtenidas = [(modulo, nombre) for modulo, nombres in grupos for nombre in nombres]

    faltantes = esperadas - set(obtenidas)
    assert not faltantes, (
        "Specs de módulos de usuario que no llegarían al Reference manual: "
        f"{sorted(faltantes)}. Revisar `build_user_groups` en "
        "manuals/build_reference_manual.py."
    )
    duplicadas = [x for x in obtenidas if obtenidas.count(x) > 1]
    assert not duplicadas, f"Specs duplicadas en el apéndice: {sorted(set(duplicadas))}"

    modulos = [m for m, _ in grupos]
    assert modulos == sorted(set(modulos)), "Un capítulo por módulo, en orden alfabético"

    # El cuerpo del manual es el programa principal: ninguna spec de un
    # módulo de usuario entra en sus capítulos (las etiquetas chocarían).
    cuerpo = set(_grouped_names())
    en_el_cuerpo = {n for _, n in esperadas} & cuerpo
    assert not en_el_cuerpo, f"Specs de módulos de usuario en el cuerpo: {sorted(en_el_cuerpo)}"


def test_specs_de_modulo_sin_paquete_abortan_el_build(tmp_path, monkeypatch):
    """Una carpeta ``docs/user/<módulo>/specs/`` sin módulo de usuario que la
    respalde detiene el build nombrándolo, en vez de omitir sus specs."""
    fuente = next(p for d in sorted(USER_DOCS_DIR.glob("*/specs")) for p in collect_specs(d))
    specs = tmp_path / "fantasma" / "specs"
    specs.mkdir(parents=True)
    shutil.copy(fuente, specs / fuente.name)

    monkeypatch.setattr(build_reference_manual, "USER_DOCS_DIR", tmp_path)
    with pytest.raises(SystemExit) as exc:
        build_reference_manual.build_user_groups()

    assert "fantasma" in str(exc.value)
