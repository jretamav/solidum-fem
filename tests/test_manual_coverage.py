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
usuario», agrupadas por módulo, junto con las secciones «Módulos de usuario —
<módulo>» de los catálogos. El invariante vale para las dos.

Los tests del apéndice no dependen de ningún módulo concreto: sin módulos de
usuario se saltan, y los que prueban los abortos del build usan módulos
fantasma en directorios temporales.
"""
import re
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


def _modulos_de_usuario() -> list[str]:
    """Módulos de usuario instalados en ``solidum/user/``."""
    import solidum
    return solidum.available_user_modules()


def _saltar_sin_modulos_de_usuario() -> None:
    """El apéndice sólo existe si hay algún módulo de usuario: sin ninguno, sus
    tests no tienen qué comprobar (el programa principal no los necesita)."""
    if not _modulos_de_usuario():
        pytest.skip("No hay módulos de usuario en solidum/user/: no hay apéndice.")


def test_documentacion_de_modulo_si_y_solo_si_existe_el_modulo():
    """``docs/user/<m>/`` existe si y sólo si existe el módulo de usuario
    ``solidum/user/<m>/``, y guarda sus specs en ``docs/user/<m>/specs/``.

    Sustituye a suponer que hay un módulo concreto: retirar un módulo exige
    retirar sus documentos, y un módulo nuevo exige los suyos (Reglas.md §4).
    """
    modulos = set(_modulos_de_usuario())
    carpetas = ({d.name for d in USER_DOCS_DIR.iterdir()
                 if d.is_dir() and not d.name.startswith(("_", "."))}
                if USER_DOCS_DIR.is_dir() else set())

    assert not carpetas - modulos, (
        f"docs/user/ documenta módulos que no existen en solidum/user/: "
        f"{sorted(carpetas - modulos)}")
    assert not modulos - carpetas, (
        f"Módulos de usuario sin documentación en docs/user/<módulo>/: "
        f"{sorted(modulos - carpetas)}")
    sin_specs = sorted(m for m in modulos
                       if not collect_specs(USER_DOCS_DIR / m / "specs"))
    assert not sin_specs, (
        f"Módulos de usuario sin specs en docs/user/<módulo>/specs/: {sin_specs}")


def test_toda_spec_de_modulo_de_usuario_entra_en_el_apendice():
    """Cada spec de ``docs/user/*/specs/`` aparece exactamente una vez en el
    apéndice «Módulos de usuario», bajo su módulo, y ninguna en el cuerpo."""
    _saltar_sin_modulos_de_usuario()
    esperadas = _user_specs_on_disk()

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
    respalde detiene el build nombrándolo, en vez de omitir sus specs.

    La spec copiada es del programa principal: el aborto ocurre antes de leerla,
    y así el test no depende de que exista algún módulo de usuario."""
    fuente = collect_specs(SPECS_DIR)[0]
    specs = tmp_path / "fantasma" / "specs"
    specs.mkdir(parents=True)
    shutil.copy(fuente, specs / fuente.name)

    monkeypatch.setattr(build_reference_manual, "USER_DOCS_DIR", tmp_path)
    with pytest.raises(SystemExit) as exc:
        build_reference_manual.build_user_groups()

    assert "fantasma" in str(exc.value)


# ---------------------------------------------------------------------------
# Secciones «Módulos de usuario — <módulo>» de los catálogos: van al capítulo
# del módulo en el apéndice, no a los capítulos de catálogo del cuerpo.
# ---------------------------------------------------------------------------

_CATALOGO_DE_PRUEBA = """# Catálogo de prueba

## Truss2D — barra

Texto del programa principal.

```python
# Módulos de usuario — otro   (comentario de un listado, no un encabezado)
x = 1
```

# Módulos de usuario — fantasma

> Nota del módulo.

## Resorte — elemento de juguete

# Elementos térmicos

## Quad4Thermal — cuadrilátero térmico
"""


def test_seccion_de_modulo_de_usuario_se_separa_del_catalogo():
    """La sección llega del H1 «Módulos de usuario — <m>» al H1 siguiente; los
    ``#`` dentro de bloques de código no son encabezados."""
    cuerpo, secciones = build_reference_manual.split_user_module_sections(
        _CATALOGO_DE_PRUEBA, "prueba.md")

    assert list(secciones) == ["fantasma"]
    seccion = secciones["fantasma"]
    assert "## Resorte" in seccion and "Nota del módulo" in seccion
    assert "Módulos de usuario" not in seccion, "El H1 lo sustituye el capítulo del apéndice"
    assert "Quad4Thermal" not in seccion, "La sección termina en el H1 siguiente"

    assert "Resorte" not in cuerpo and "Nota del módulo" not in cuerpo
    assert "fantasma" not in cuerpo
    for texto in ("## Truss2D", "# Módulos de usuario — otro", "# Elementos térmicos",
                  "## Quad4Thermal"):
        assert texto in cuerpo


def test_seccion_de_modulo_sin_nombre_aborta_el_build():
    """Un H1 «Módulos de usuario» que no nombra módulo detiene el build: sus
    entradas no pueden repartirse, y no deben caer en el cuerpo."""
    with pytest.raises(SystemExit) as exc:
        build_reference_manual.split_user_module_sections(
            "# Catálogo\n\n# Módulos de usuario\n\n## Resorte — x\n", "catalogo_prueba.md")
    assert "catalogo_prueba.md" in str(exc.value)


def test_seccion_de_catalogo_de_modulo_inexistente_aborta_el_build(tmp_path, monkeypatch):
    """Una sección de catálogo de un módulo que no existe en ``solidum/user/``
    detiene el build nombrándolo, en vez de omitirse o caer en el cuerpo."""
    catalogo = tmp_path / "catalogo_prueba.md"
    catalogo.write_text(_CATALOGO_DE_PRUEBA, encoding="utf-8")
    monkeypatch.setattr(build_reference_manual, "APPENDIX_CHAPTERS",
                        [("Prueba — catálogo", str(catalogo))])

    with pytest.raises(SystemExit) as exc:
        build_reference_manual.build_catalog_chapters()
    assert "fantasma" in str(exc.value)


def test_seccion_de_catalogo_va_al_capitulo_de_su_modulo(tmp_path, monkeypatch):
    """Con un módulo que sí existe, su sección sale del cuerpo y se asigna a
    él, con el título y el archivo del catálogo del que procede."""
    _saltar_sin_modulos_de_usuario()
    modulo = _modulos_de_usuario()[0]
    catalogo = tmp_path / "catalogo_prueba.md"
    catalogo.write_text(_CATALOGO_DE_PRUEBA.replace("fantasma", modulo), encoding="utf-8")
    monkeypatch.setattr(build_reference_manual, "APPENDIX_CHAPTERS",
                        [("Prueba — catálogo", str(catalogo))])

    cuerpo, por_modulo = build_reference_manual.build_catalog_chapters()

    assert [titulo for titulo, _, _ in cuerpo] == ["Prueba — catálogo"]
    assert "Resorte" not in cuerpo[0][2]
    assert list(por_modulo) == [modulo]
    [(titulo, fuente, seccion)] = por_modulo[modulo]
    assert titulo == "Prueba — catálogo" and fuente == catalogo
    assert "## Resorte" in seccion


_ENCABEZADO_TEX = re.compile(r"^\\(?:sub){0,2}section\*?\{(.*)\}[ \t]*$", re.MULTILINE)


def _primer_termino(titulo: str) -> str:
    """Primera palabra de un título LaTeX o Markdown, sin macros, llaves ni
    comillas de código (``CST\\_Embedded2D --- …`` → ``CST_Embedded2D``)."""
    texto = titulo.replace("\\allowbreak{}", "").replace("\\_", "_")
    texto = re.sub(r"\\[A-Za-z]+\*?", " ", texto)
    for ch in "{}`*":
        texto = texto.replace(ch, " ")
    partes = texto.split()
    return partes[0] if partes else ""


def _encabezados_md(md: str) -> list[str]:
    """Títulos de los encabezados ``##`` a ``####`` fuera de bloques de código."""
    titulos, en_bloque = [], False
    for linea in md.split("\n"):
        if linea.lstrip().startswith("```"):
            en_bloque = not en_bloque
        elif not en_bloque and (m := re.match(r"#{2,4}[ \t]+(.+)$", linea)):
            titulos.append(m.group(1).strip())
    return titulos


@pytest.fixture(scope="module")
def manual_tex() -> str:
    """El .tex completo del Reference manual, sin compilar."""
    _saltar_sin_modulos_de_usuario()
    return build_reference_manual.assemble()


def test_ningun_componente_de_modulo_de_usuario_en_el_cuerpo_del_manual(manual_tex):
    """Todo lo de un módulo de usuario va tras ``\\appendix``: ni su spec ni su
    entrada de catálogo aparecen como encabezado en los capítulos del cuerpo.

    Hasta 2026-09-29 las secciones «Módulos de usuario — <módulo>» de los
    catálogos entraban en los capítulos de catálogo del cuerpo como una
    \\section más, al mismo nivel que los componentes del programa principal:
    el builder descartaba su H1."""
    specs_de_usuario = {n for _, n in _user_specs_on_disk()}
    assert manual_tex.count("\\appendix\n") == 1, "El apéndice debe abrirse una vez"
    cuerpo, apendice = manual_tex.split("\\appendix\n", 1)

    en_el_cuerpo = sorted({t for t in _ENCABEZADO_TEX.findall(cuerpo)
                           if _primer_termino(t) in specs_de_usuario})
    assert not en_el_cuerpo, (
        "Encabezados de componentes de módulos de usuario en el cuerpo del "
        f"manual: {en_el_cuerpo}")

    for nombre in sorted(specs_de_usuario):
        etiqueta = f"\\label{{spec:{nombre}}}"
        assert etiqueta not in cuerpo, f"La spec {nombre} entra en el cuerpo"
        assert apendice.count(etiqueta) == 1, f"La spec {nombre} no está en el apéndice"

    # Y nada de las secciones de catálogo de los módulos se pierde por el
    # camino: cada uno de sus encabezados reaparece en el apéndice.
    en_el_apendice = {_primer_termino(t) for t in _ENCABEZADO_TEX.findall(apendice)}
    _, por_modulo = build_reference_manual.build_catalog_chapters()
    perdidos = sorted({t for secciones in por_modulo.values() for _, _, md in secciones
                       for t in _encabezados_md(md)
                       if _primer_termino(t) not in en_el_apendice})
    assert not perdidos, f"Entradas de catálogo de módulos de usuario perdidas: {perdidos}"
