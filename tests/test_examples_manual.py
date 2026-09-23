"""Verifica lo que afirma el manual de ejemplos (``manuals/Example_manual.pdf``).

Cada capítulo de ``manuals/sources/examples/`` está ligado a una carpeta de
``examples/`` cuyo ``run.py`` expone ``calcular()`` y ``TOLERANCIAS``. El
manual cita cifras de ``resultados.json`` sin escribirlas a mano; este módulo
cierra el lazo:

1. ``calcular()`` se ejecuta y cada error declarado en ``TOLERANCIAS`` debe
   quedar por debajo de su cota: la física que el capítulo afirma se sostiene.
2. El capítulo renderizado con los resultados recién calculados debe ser
   idéntico, carácter a carácter, al renderizado con el ``resultados.json``
   versionado: si un cambio del código mueve una cifra **que el manual
   imprime**, el manual quedó desactualizado y hay que regenerarlo
   (``python manuals/build_example_manual.py --recalcular``).

   Se compara el texto y no los números crudos a propósito: las cifras al
   nivel del redondeo varían entre ejecuciones (Pardiso reparte las sumas
   entre hilos) y entre plataformas, sin que cambie nada de lo que el lector
   ve. Medido 2026-09-23: un error de 2e-8 difería en la sexta cifra entre dos
   ejecuciones idénticas.
3. Todo marcador ``{{r:…}}`` del capítulo existe en ``resultados.json``.

No compila LaTeX ni dibuja figuras (``calcular()`` no usa matplotlib).
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "manuals"))
sys.path.insert(0, str(ROOT / "examples"))

from build_example_manual import (  # noqa: E402
    SRC_DIR,
    chapter_sources,
    expand_placeholders,
)
from ejemplos_comun import _a_json  # noqa: E402

CAPITULOS = chapter_sources()
IDS = [ex.name for _, ex in CAPITULOS]


def _cargar(ex_dir: Path):
    spec = importlib.util.spec_from_file_location(f"ejemplo_{ex_dir.name}", ex_dir / "run.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_CACHE: dict[str, tuple] = {}


def _calcular(ex_dir: Path):
    if ex_dir.name not in _CACHE:
        module = _cargar(ex_dir)
        # Ida y vuelta por JSON: lo que se compara es lo que el manual leería.
        _CACHE[ex_dir.name] = (module, json.loads(json.dumps(_a_json(module.calcular()))))
    return _CACHE[ex_dir.name]


def _valor(res: dict, ruta: str):
    for parte in ruta.split("."):
        res = res[parte]
    return res


@pytest.mark.parametrize("ex_dir", [ex for _, ex in CAPITULOS], ids=IDS)
def test_errores_dentro_de_tolerancia(ex_dir):
    module, res = _calcular(ex_dir)
    assert module.TOLERANCIAS, f"{ex_dir.name}/run.py no declara TOLERANCIAS"
    for ruta, tol in module.TOLERANCIAS.items():
        valor = _valor(res, ruta)
        if isinstance(tol, tuple):          # intervalo (mínimo, máximo)
            lo, hi = tol
            assert lo <= valor <= hi, f"{ex_dir.name}: {ruta} = {valor:.4e} fuera de [{lo}, {hi}]"
        else:                               # cota superior
            assert valor <= tol, f"{ex_dir.name}: {ruta} = {valor:.3e} > {tol:.1e}"


@pytest.mark.parametrize("src,ex_dir", CAPITULOS, ids=IDS)
def test_manual_al_dia_con_el_codigo(src, ex_dir, tmp_path):
    _, res = _calcular(ex_dir)
    md = src.read_text(encoding="utf-8")
    impreso = expand_placeholders(md, ex_dir, src.name)
    # Mismo capítulo, con los resultados recién calculados.
    (tmp_path / "resultados.json").write_text(json.dumps(res), encoding="utf-8")
    for f in ex_dir.iterdir():                    # listados {{yaml:…}} / {{py:…}}
        if f.suffix in (".yaml", ".py"):
            (tmp_path / f.name).write_bytes(f.read_bytes())
    actual = expand_placeholders(md, tmp_path, src.name)
    difs = [(a, b) for a, b in zip(impreso.splitlines(), actual.splitlines()) if a != b]
    assert not difs, (
        f"El capítulo {src.name} imprime cifras que el código ya no produce "
        f"({len(difs)} líneas; primera: {difs[0][0][:160]!r} → {difs[0][1][:160]!r}). "
        "Regenerar con 'python manuals/build_example_manual.py --recalcular'.")


@pytest.mark.parametrize("src,ex_dir", CAPITULOS, ids=IDS)
def test_marcadores_del_capitulo_resuelven(src, ex_dir):
    expand_placeholders(src.read_text(encoding="utf-8"), ex_dir, src.name)


def test_hay_capitulos():
    assert CAPITULOS, f"sin capítulos en {SRC_DIR}"
