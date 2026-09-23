"""Genera Example_manual.pdf — problemas resueltos con Solidum FEM.

Cada capítulo es un archivo ``manuals/sources/examples/NN_nombre.md`` asociado
a una carpeta de ``examples/`` mediante un comentario en su cabecera::

    <!-- ejemplo: voladizo_marco -->

La carpeta contiene el modelo y un ``run.py`` que resuelve el problema y
escribe ``resultados.json`` y las figuras. **El capítulo no lleva cifras
escritas a mano**: las toma de ``resultados.json`` al compilar, mediante
marcadores que este builder sustituye antes de convertir a LaTeX:

``{{r:clave.subclave:fmt}}``
    Valor de ``resultados.json`` con el formato de Python ``fmt`` (opcional).
    Los formatos exponenciales se escriben como potencia de diez y los
    negativos en modo matemático (signo menos tipográfico, no guion).
``{{yaml:archivo}}``, ``{{py:archivo}}``, ``{{py:archivo#funcion}}``
    Listado del archivo de la carpeta del ejemplo (o sólo de esa función),
    copiado del archivo real: el manual no puede mostrar un modelo distinto
    del que se ejecuta.

Las figuras se citan con la sintaxis común del conversor
(``![pie](ruta){#fig:nombre}`` y ``@fig:nombre``, ver ``md_to_latex``).

Uso::

    python manuals/build_example_manual.py               # usa los resultados guardados
    python manuals/build_example_manual.py --recalcular  # ejecuta antes cada run.py

Salida: ``manuals/Example_manual.tex`` y ``manuals/Example_manual.pdf``.
"""
from __future__ import annotations

import ast
import importlib.util
import json
import re
import sys
import textwrap
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT / "manuals" / "sources" / "examples"
EXAMPLES_DIR = ROOT / "examples"
OUT_DIR = ROOT / "manuals"
OUT_TEX = OUT_DIR / "Example_manual.tex"
OUT_PDF = OUT_DIR / "Example_manual.pdf"

sys.path.insert(0, str(OUT_DIR))
from build_reference_manual import (  # noqa: E402
    md_to_latex,
    report_broken_links,
    set_link_context,
)
from build_user_manual import (  # noqa: E402
    build_colofon,
    build_preamble,
    compile_pdf,
    reroute_advertencia_boxes,
)

SOBRE_EJEMPLOS = r"""\noindent Este manual reúne \textbf{problemas resueltos} con Solidum FEM. Cada capítulo plantea un problema de mecánica de sólidos, lo modela paso a paso ---con el archivo \texttt{.yaml} o con el API de Python---, compara los resultados con una solución analítica o de referencia publicada y comenta qué enseña el ejemplo sobre el método o sobre el programa.

\textbf{Ninguna cifra ni figura está escrita a mano.} Cada ejemplo tiene una carpeta en \texttt{examples/} con su modelo y un script \texttt{run.py} que resuelve el problema, guarda los resultados en \texttt{resultados.json} y dibuja las figuras; el manual los lee de ahí al compilarse. Un test de la suite (\texttt{tests/test\_examples\_manual.py}) ejecuta los mismos scripts y comprueba que los errores citados siguen dentro de su tolerancia: lo que este manual afirma se verifica en cada cambio del código.

Para reproducir un ejemplo basta con
\begin{center}\texttt{python examples/\textit{carpeta}/run.py}\end{center}

Para la sintaxis completa del archivo de entrada, consultar el \emph{Manual de Usuario} (\texttt{manuals/User\_manual.pdf}); para la formulación de cada elemento, material y solver, el \emph{Manual de Referencia} (\texttt{manuals/Reference\_manual.pdf}).

Este manual se regenera con:
\begin{center}\texttt{python manuals/build\_example\_manual.py}\end{center}

\noindent No editar este PDF manualmente; cualquier corrección debe hacerse sobre los capítulos en \texttt{manuals/sources/examples/} o sobre los scripts de \texttt{examples/}.

"""

PREAMBLE = build_preamble(
    titulo="Manual de Ejemplos",
    sigla="FF-ME",
    cabecera="Ejemplos",
    subtitulo="Problemas resueltos y verificados contra\\\\\n"
              "            soluciones analíticas y de referencia",
    sobre=SOBRE_EJEMPLOS,
    subject="Problemas resueltos con Solidum FEM, verificados contra soluciones analíticas",
    keywords="elementos finitos, mecánica de sólidos, ejemplos, verificación, Solidum FEM",
)

POSTAMBLE = r"""
\end{document}
"""

_EJEMPLO_RE = re.compile(r"<!--\s*ejemplo:\s*([\w-]+)\s*-->")
_RESULT_RE = re.compile(r"\{\{r:([\w.]+)(?::([^}]*))?\}\}")
_LISTING_RE = re.compile(r"\{\{(yaml|py):([\w./-]+)(?:#(\w+))?\}\}")


def _lookup(res: dict, key: str, source: str):
    value = res
    for part in key.split("."):
        if not isinstance(value, dict) or part not in value:
            raise KeyError(f"{source}: '{key}' no está en resultados.json")
        value = value[part]
    return value


def format_value(value, spec: str | None) -> str:
    """Formatea una cifra para el texto del manual.

    ``1.2e-15`` → ``$1.2\\times 10^{-15}$``; negativos entre ``$`` para que
    el signo sea un menos y no un guion."""
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        return str(value)
    s = format(value, spec) if spec else repr(value)
    m = re.fullmatch(r"(-?\d+(?:\.\d+)?)e([+-]\d+)", s)
    if m:
        mant, exp = m.group(1), int(m.group(2))
        power = f"10^{{{exp}}}"
        return f"${power}$" if mant in ("1", "1.0") else f"${mant}\\times {power}$"
    return f"${s}$" if s.startswith("-") else s


def _function_source(path: Path, name: str) -> str:
    text = path.read_text(encoding="utf-8")
    for node in ast.walk(ast.parse(text)):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name:
            start = min([node.lineno] + [d.lineno for d in node.decorator_list])
            lines = text.splitlines()[start - 1:node.end_lineno]
            return textwrap.dedent("\n".join(lines))
    raise KeyError(f"{path.name}: no hay función '{name}'")


def expand_placeholders(md: str, ex_dir: Path, source: str) -> str:
    """Sustituye los marcadores ``{{r:…}}``, ``{{yaml:…}}`` y ``{{py:…}}``."""
    res_path = ex_dir / "resultados.json"
    if not res_path.exists():
        raise FileNotFoundError(
            f"{source}: falta {res_path.relative_to(ROOT)}; ejecutar "
            f"'python examples/{ex_dir.name}/run.py' o compilar con --recalcular")
    res = json.loads(res_path.read_text(encoding="utf-8"))

    def _listing(m: re.Match) -> str:
        lang, fname, func = m.group(1), m.group(2), m.group(3)
        path = ex_dir / fname
        if not path.exists():
            raise FileNotFoundError(f"{source}: no existe {path.relative_to(ROOT)}")
        code = _function_source(path, func) if func else path.read_text(encoding="utf-8")
        fence = "python" if lang == "py" else "yaml"
        return f"```{fence}\n{code.rstrip()}\n```"

    md = _LISTING_RE.sub(_listing, md)
    md = _RESULT_RE.sub(lambda m: format_value(_lookup(res, m.group(1), source), m.group(2)), md)
    leftover = re.findall(r"\{\{[^}]*\}\}", md)
    if leftover:
        raise ValueError(f"{source}: marcadores no reconocidos {leftover[:3]}")
    return md


def chapter_sources() -> list[tuple[Path, Path]]:
    """``[(fuente .md, carpeta del ejemplo)]`` en orden de capítulo."""
    out = []
    for src in sorted(SRC_DIR.glob("*.md")):
        m = _EJEMPLO_RE.search(src.read_text(encoding="utf-8"))
        if not m:
            raise ValueError(f"{src.name}: falta el comentario '<!-- ejemplo: carpeta -->'")
        ex_dir = EXAMPLES_DIR / m.group(1)
        if not (ex_dir / "run.py").exists():
            raise FileNotFoundError(f"{src.name}: {ex_dir.relative_to(ROOT)}/run.py no existe")
        out.append((src, ex_dir))
    if not out:
        raise SystemExit(f"[!] No hay fuentes en {SRC_DIR}")
    return out


def run_example(ex_dir: Path) -> None:
    """Ejecuta ``main()`` del ``run.py`` del ejemplo (resultados + figuras)."""
    spec = importlib.util.spec_from_file_location(f"ejemplo_{ex_dir.name}", ex_dir / "run.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.main()


def assemble() -> str:
    parts = [PREAMBLE]
    set_link_context(manual="examples", internal_specs=set(), internal_adrs=set(),
                     heading_offset=0)
    for src, ex_dir in chapter_sources():
        md = src.read_text(encoding="utf-8")
        title = next((ln[2:].strip() for ln in md.splitlines() if ln.startswith("# ")), src.stem)
        md = expand_placeholders(md, ex_dir, src.name)
        md = re.sub(r"<!--.*?-->", "", md, flags=re.DOTALL)
        set_link_context(source=src, label_prefix=ex_dir.name)
        ltx = reroute_advertencia_boxes(md_to_latex(md))
        # La cabecera lleva sólo lo que precede a los dos puntos del título:
        # el título completo no cabe junto a "Solidum FEM --- Ejemplos".
        short = title.split(":")[0].strip()
        parts.append(f"\\chapter{{{title}}}\n\\chaptermark{{{short}}}\n"
                     f"\\label{{cap:{ex_dir.name}}}\n")
        parts.append(ltx)
        parts.append("\n\\newpage\n")
    parts.append(build_colofon(
        titulo="Manual de Ejemplos", sigla="FF-ME",
        fuentes="\\texttt{manuals/sources/examples/} y los scripts de \\texttt{examples/}",
        script="build\\_example\\_manual.py"))
    report_broken_links()
    return "\n".join(parts) + POSTAMBLE


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--recalcular" in argv:
        for _, ex_dir in chapter_sources():
            print(f"Ejecutando examples/{ex_dir.name}/run.py ...")
            run_example(ex_dir)
    print(f"Leyendo fuentes desde: {SRC_DIR}")
    tex_content = assemble()
    OUT_TEX.write_text(tex_content, encoding="utf-8")
    print(f"  -> .tex escrito en: {OUT_TEX} ({len(tex_content):,} chars)")
    print("Compilando con lualatex (2 pasadas)...")
    if compile_pdf(OUT_TEX):
        print(f"  -> PDF generado: {OUT_PDF}")
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
