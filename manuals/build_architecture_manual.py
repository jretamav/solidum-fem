"""Genera Architecture_manual.pdf — visión arquitectural sin entrar al detalle.

Reutiliza la conversión Markdown→LaTeX y las tablas Unicode de
`build_reference_manual.py`. Las fuentes son los archivos numerados en
`manuals/sources/architecture/`, que se ensamblan en orden alfabético como
capítulos del libro.

Uso:
    python manuals/build_architecture_manual.py

Salida:
    manuals/Architecture_manual.tex
    manuals/Architecture_manual.pdf
"""
from __future__ import annotations

import datetime
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT / "manuals" / "sources" / "architecture"
OUT_DIR = ROOT / "manuals"
OUT_TEX = OUT_DIR / "Architecture_manual.tex"
OUT_PDF = OUT_DIR / "Architecture_manual.pdf"

# Reutilizar md_to_latex del builder de referencia.
sys.path.insert(0, str(OUT_DIR))
from build_reference_manual import md_to_latex  # noqa: E402


PREAMBLE = r"""\documentclass[11pt,letterpaper,oneside]{report}

% Motor: lualatex (Unicode nativo, sin necesidad de inputenc/fontenc).
\usepackage{fontspec}
\usepackage[spanish,es-tabla]{babel}
\usepackage[margin=2.5cm, headheight=15pt]{geometry}
\usepackage{listings}
\usepackage{xcolor}
\usepackage{hyperref}
\usepackage{titlesec}
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{amsfonts}
\usepackage{booktabs}
\usepackage{tabularx}
\usepackage{array}
\usepackage{ragged2e}
\usepackage{fancyhdr}
\usepackage{microtype}
\usepackage{tikz}
\usetikzlibrary{positioning, arrows.meta, shapes.geometric}
\usepackage[most]{tcolorbox}

% Columna X de tabularx con justificación a la izquierda sin guionado agresivo.
\newcolumntype{Y}{>{\RaggedRight\arraybackslash\footnotesize}X}

\definecolor{yamlkey}{RGB}{0, 102, 204}
\definecolor{yamlval}{RGB}{204, 0, 0}
\definecolor{yamlcomment}{RGB}{120, 120, 120}
\definecolor{bglight}{RGB}{245, 245, 245}
\definecolor{darkgray}{RGB}{50, 50, 50}
\definecolor{pendienteorange}{RGB}{200, 100, 0}

% Comando para texto pendiente: etiqueta [Pendiente] + texto en cursiva,
% ambos en naranja oscuro. La etiqueta se sintetiza en el conversor MD->LaTeX
% a partir del marcador `[PENDIENTE: ...]` en las fuentes.
\newcommand{\pendiente}[1]{\textcolor{pendienteorange}{\textit{[Pendiente]}~\textcolor{pendienteorange}{\textit{#1}}}}

% Cuadro de informacion complementaria: borde discreto, fondo muy claro,
% un titulo en negrita en la parte superior. Se usa para apartados como
% Justificacion, Riesgo de arquitectura, Razones del rechazo, etc.
\newtcolorbox{infobox}[1]{
  enhanced,
  colback=bglight, colframe=darkgray, colbacktitle=bglight,
  coltitle=yamlkey, fonttitle=\bfseries\small,
  title=#1, attach boxed title to top left={yshift=-2mm, xshift=4mm},
  boxed title style={colback=white, colframe=white, sharp corners},
  boxrule=0.4pt, arc=2pt,
  left=4mm, right=4mm, top=2mm, bottom=2mm,
  before skip=8pt, after skip=8pt,
}

\hypersetup{
    colorlinks=true,
    linkcolor=yamlkey,
    urlcolor=yamlkey,
    pdftitle={Manual de Arquitectura - Solidum FEM},
    pdfauthor={Jaime Retama Velasco}
}

\titleformat{\chapter}[display]
  {\normalfont\huge\bfseries\color{yamlkey}}{\chaptertitlename\ \thechapter}{20pt}{\Huge}
\titleformat{\section}{\Large\bfseries\color{darkgray}}{\thesection}{1em}{}
\titleformat{\subsection}{\large\bfseries\color{darkgray}}{\thesubsection}{1em}{}
\titleformat{\subsubsection}{\normalsize\bfseries\color{darkgray}}{\thesubsubsection}{1em}{}

% Redefinicion de \chaptermark para que el encabezado muestre
% "Cap. N - Titulo" en lugar del texto en mayusculas por defecto.
\renewcommand{\chaptermark}[1]{\markboth{Cap.\ \thechapter\ ---\ #1}{}}

\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\textbf{\color{darkgray}Solidum FEM --- Arquitectura}}
\fancyhead[R]{\color{darkgray}\nouppercase{\leftmark}}
\fancyfoot[C]{\thepage}
\fancyfoot[R]{\footnotesize\color{darkgray}FF-MA}
\renewcommand{\headrulewidth}{0.4pt}
\renewcommand{\footrulewidth}{0.4pt}

\lstset{
  basicstyle=\ttfamily\footnotesize,
  backgroundcolor=\color{bglight},
  frame=single,
  rulecolor=\color{lightgray},
  breaklines=true,
  showstringspaces=false,
  extendedchars=true,
  captionpos=b,
  literate={á}{{\'a}}1 {é}{{\'e}}1 {í}{{\'i}}1 {ó}{{\'o}}1 {ú}{{\'u}}1
           {Á}{{\'A}}1 {É}{{\'E}}1 {Í}{{\'I}}1 {Ó}{{\'O}}1 {Ú}{{\'U}}1
           {ñ}{{\~n}}1 {Ñ}{{\~N}}1 {ü}{{\"u}}1 {Ü}{{\"U}}1
           {§}{{\S{}}}1
           {ε}{{$\varepsilon$}}1 {σ}{{$\sigma$}}1 {α}{{$\alpha$}}1
           {β}{{$\beta$}}1 {γ}{{$\gamma$}}1 {δ}{{$\delta$}}1
           {θ}{{$\theta$}}1 {κ}{{$\kappa$}}1 {λ}{{$\lambda$}}1
           {μ}{{$\mu$}}1 {ν}{{$\nu$}}1 {π}{{$\pi$}}1
           {ρ}{{$\rho$}}1 {τ}{{$\tau$}}1 {φ}{{$\varphi$}}1
           {ω}{{$\omega$}}1 {Δ}{{$\Delta$}}1 {Ω}{{$\Omega$}}1
           {Φ}{{$\Phi$}}1
           {⇔}{{$\Leftrightarrow$}}1 {⇒}{{$\Rightarrow$}}1 {⇐}{{$\Leftarrow$}}1
           {→}{{$\to$}}1 {←}{{$\leftarrow$}}1 {↔}{{$\leftrightarrow$}}1
           {≤}{{$\le$}}1 {≥}{{$\ge$}}1 {≠}{{$\neq$}}1
           {·}{{$\cdot$}}1 {×}{{$\times$}}1 {±}{{$\pm$}}1
           {°}{{$^{\circ}$}}1 {²}{{$^{2}$}}1 {³}{{$^{3}$}}1
           {₀}{{$_{0}$}}1 {₁}{{$_{1}$}}1 {₂}{{$_{2}$}}1
           {ᵀ}{{$^{\top}$}}1 {⁺}{{$^{+}$}}1 {⁻}{{$^{-}$}}1
           {—}{{---}}1 {–}{{--}}1 {−}{{-}}1
           {“}{{``}}1 {”}{{''}}1
}

\begin{document}

\begin{titlepage}
    \centering
    \vspace*{2cm}
    {\Huge\bfseries\color{yamlkey} Solidum FEM \par}
    \vspace{1cm}
    {\LARGE Manual de Arquitectura \par}
    \vspace{0.3cm}
    {\normalsize\color{darkgray} Sigla: \texttt{FF-MA} \par}
    \vspace{0.5cm}
    \rule{\linewidth}{0.5mm} \par
    \vspace{2cm}
    {\large Visión arquitectural del programa: capas, bloques funcionales,\\
            mecanismos transversales y dirección de evolución\par}
    \vspace{2cm}
    {\large Lectura del arquitecto, sin entrar al detalle de cada componente\par}
    \vspace{1cm}
    {\Large \textbf{Autor:} Jaime Retama Velasco \par}
    \vspace{0.5cm}
    {\large Facultad de Estudios Superiores Aragón \\ Universidad Nacional Autónoma de México \par}
    \vfill
    {\large \today \par}
\end{titlepage}

\pagenumbering{roman}
\setcounter{page}{1}
\tableofcontents
\newpage

\chapter*{Sobre este manual}
\addcontentsline{toc}{chapter}{Sobre este manual}
\noindent Este manual ofrece una \textbf{visión arquitectural} de Solidum FEM. Su objetivo es que el lector pueda reconstruir mentalmente cómo está organizado el programa, qué piezas tiene y cómo encajan, sin entrar al detalle de la formulación de cada elemento ni al algoritmo concreto de cada solver.

Para la \emph{referencia formal} de cada componente (ecuaciones, matrices $\mathbf B$, contratos YAML, criterios de aceptación), consultar \texttt{manuals/Reference\_manual.pdf}, generado automáticamente desde \texttt{docs/specs/}.

Para una guía orientada al uso del programa (sintaxis YAML, ejemplos, post-procesamiento), consultar \texttt{manuals/User\_manual.pdf}.

Este manual se regenera con:
\begin{center}\texttt{python manuals/build\_architecture\_manual.py}\end{center}

\noindent No editar este PDF manualmente; cualquier corrección debe hacerse sobre los archivos fuente en \texttt{manuals/sources/architecture/}.

\newpage
\pagenumbering{arabic}
\setcounter{page}{1}

"""

POSTAMBLE = r"""
\end{document}
"""


def chapter_title(stem: str) -> str:
    """Convierte '01_filosofia' → 'Filosofía' usando el primer H1 del archivo."""
    md = (SRC_DIR / f"{stem}.md").read_text(encoding="utf-8")
    for line in md.splitlines():
        if line.startswith("# "):
            return line[2:].strip()
    return stem


# Marcador en fuentes Markdown:
#   [PENDIENTE: texto del pendiente]
# Se traduce a \pendiente{...} en LaTeX y se acumula en una lista global
# para construir el apéndice consolidado.
PENDIENTE_RE = re.compile(r"\[PENDIENTE:\s*(.+?)\]", re.DOTALL)


def extract_pendientes(md: str, chapter_title: str) -> tuple[str, list[str]]:
    """Sustituye `[PENDIENTE: ...]` por placeholders y devuelve la lista.

    El reemplazo se hace ANTES de la conversión a LaTeX para evitar que el
    intérprete de Markdown trate los corchetes como sintaxis de enlace.
    Tras la conversión, los placeholders se sustituyen por `\\pendiente{...}`.
    """
    found: list[str] = []

    def _repl(m: re.Match) -> str:
        text = m.group(1).strip()
        # Colapsar saltos de línea internos para que el placeholder no se rompa.
        text = re.sub(r"\s+", " ", text)
        idx = len(found)
        found.append(text)
        return f"@@PEND{idx}@@"

    md_clean = PENDIENTE_RE.sub(_repl, md)
    return md_clean, found


def restore_pendientes(ltx: str, pendientes: list[str]) -> str:
    """Sustituye los placeholders @@PENDn@@ por \\pendiente{...} con el texto.

    El texto del pendiente es ya texto en español que debe pasar por el
    conversor para que sus signos LaTeX se escapen correctamente. Como el
    placeholder se inyectó antes de la conversión, el texto restaurado debe
    convertirse de forma aislada y sin la fase de placeholders, lo que se
    logra invocando md_to_latex sobre cada fragmento.
    """
    for idx, text in enumerate(pendientes):
        latex_text = md_to_latex(text).strip()
        # md_to_latex envuelve a veces en párrafos; quitar líneas vacías finales.
        ltx = ltx.replace(f"@@PEND{idx}@@", f"\\pendiente{{{latex_text}}}")
    return ltx


def assemble() -> str:
    parts = [PREAMBLE]
    sources = sorted(SRC_DIR.glob("*.md"))
    if not sources:
        raise SystemExit(f"[!] No hay fuentes en {SRC_DIR}")

    pendientes_globales: list[tuple[str, str]] = []  # (capítulo, texto)

    for idx, src in enumerate(sources, start=1):
        title = chapter_title(src.stem)
        md = src.read_text(encoding="utf-8")
        md_clean, pendientes = extract_pendientes(md, title)
        for p in pendientes:
            pendientes_globales.append((title, p))
        ltx = md_to_latex(md_clean)
        ltx = restore_pendientes(ltx, pendientes)
        # Inyectar etiquetas de capítulo (cap:N) y de ADR (adr:000N) para
        # los hipervínculos cruzados internos.
        ltx = inject_adr_labels(ltx)
        parts.append(f"\\chapter{{{title}}}\n\\label{{cap:{idx}}}\n")
        parts.append(ltx)
        parts.append("\n\\newpage\n")

    # Apéndice consolidado de pendientes.
    if pendientes_globales:
        parts.append(build_appendix_pendientes(pendientes_globales))

    # Página final de colofón.
    parts.append(build_colofon())

    body = "\n".join(parts)

    # Pasada final de hipervínculos cruzados sobre todo el cuerpo.
    body = link_cross_references(body)

    return body + POSTAMBLE


# Patrones que se transforman a hipervínculos en el LaTeX final. La pasada se
# realiza sobre el documento completo, después de la sustitucion de
# placeholders, para no interferir con el escape de caracteres ni con la
# preservacion de bloques verbatim.

ADR_LABEL_RE = re.compile(
    r"\\section\{ADR\s+(\d{4})\b([^}]*)\}"
)


def inject_adr_labels(ltx: str) -> str:
    """Añade `\\label{adr:NNNN}` después de cada sección "ADR NNNN".

    Solo se aplica al capítulo 8 (decisiones de arquitectura) en la práctica;
    si la cadena no aparece en otros capítulos, la regex no encuentra
    coincidencias y la función es transparente.
    """
    def _repl(m: re.Match) -> str:
        num = m.group(1)
        return m.group(0) + f"\n\\label{{adr:{num}}}"

    return ADR_LABEL_RE.sub(_repl, ltx)


# Mapa de capítulos por número (1..N) → label.
CHAPTER_REF_RE = re.compile(r"cap[íi]tulo\s+(\d{1,2})")
# ADR mencionado en el cuerpo del texto.
ADR_REF_RE = re.compile(r"\bADR\s+(\d{4})\b")


def link_cross_references(ltx: str) -> str:
    """Sustituye menciones a "capítulo N" y "ADR NNNN" por hipervínculos.

    No toca las apariciones dentro de comandos de etiquetado o de definición
    (\\label, \\chapter, \\section), ni los nombres en los títulos
    propiamente dichos del capítulo de ADR (que ya son los destinos).
    """
    # Excluir las menciones que aparecen dentro de \section{ADR NNNN ...}
    # u otros comandos: el truco es procesar línea a línea y saltar las que
    # comienzan por \section, \chapter o \label. Las menciones a hipervincular
    # estan en el cuerpo del texto, no en encabezados.
    out_lines: list[str] = []
    skip_prefixes = ("\\section", "\\chapter", "\\label", "\\subsection")
    for line in ltx.splitlines():
        stripped = line.lstrip()
        if any(stripped.startswith(p) for p in skip_prefixes):
            out_lines.append(line)
            continue

        def _repl_cap(m: re.Match) -> str:
            n = m.group(1)
            return f"\\hyperref[cap:{n}]{{{m.group(0)}}}"

        def _repl_adr(m: re.Match) -> str:
            n = m.group(1)
            return f"\\hyperref[adr:{n}]{{{m.group(0)}}}"

        line = CHAPTER_REF_RE.sub(_repl_cap, line)
        line = ADR_REF_RE.sub(_repl_adr, line)
        out_lines.append(line)
    return "\n".join(out_lines)


def get_git_commit() -> str:
    """Devuelve el hash corto del commit actual o 'no disponible' si no hay git."""
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=ROOT, capture_output=True, text=True, timeout=5,
        )
        if result.returncode == 0:
            return result.stdout.strip()
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass
    return "no disponible"


def get_git_dirty() -> bool:
    """Indica si hay modificaciones sin commit en el árbol de trabajo."""
    try:
        result = subprocess.run(
            ["git", "status", "--porcelain"],
            cwd=ROOT, capture_output=True, text=True, timeout=5,
        )
        if result.returncode == 0:
            return bool(result.stdout.strip())
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass
    return False


def build_colofon() -> str:
    """Página final con metadatos de la compilación: commit, fecha, sigla."""
    commit = get_git_commit()
    dirty = " (con cambios sin confirmar)" if get_git_dirty() else ""
    fecha = datetime.datetime.now().strftime("%d de %B de %Y, %H:%M")
    # Tradución manual de los meses al español (no asumimos locale del sistema).
    meses = {
        "January": "enero", "February": "febrero", "March": "marzo",
        "April": "abril", "May": "mayo", "June": "junio",
        "July": "julio", "August": "agosto", "September": "septiembre",
        "October": "octubre", "November": "noviembre", "December": "diciembre",
    }
    for en, es in meses.items():
        fecha = fecha.replace(en, es)

    return (
        "\\clearpage\n"
        "\\thispagestyle{empty}\n"
        "\\vspace*{\\fill}\n"
        "\\begin{center}\n"
        "\\rule{0.4\\textwidth}{0.4pt}\\\\[1em]\n"
        "{\\bfseries\\color{yamlkey} Solidum FEM --- Manual de Arquitectura}\\\\[0.5em]\n"
        f"Sigla del manual: \\texttt{{FF-MA}}\\\\[0.5em]\n"
        f"Compilado el {fecha}\\\\[0.5em]\n"
        f"Commit Git: \\texttt{{{commit}}}{dirty}\\\\[1em]\n"
        "\\rule{0.4\\textwidth}{0.4pt}\\\\[2em]\n"
        "\\parbox{0.75\\textwidth}{\\small\\itshape\\centering\n"
        "Documento generado automáticamente desde los archivos fuente en "
        "\\texttt{manuals/sources/architecture/} mediante "
        "\\texttt{python manuals/build\\_architecture\\_manual.py}. "
        "No editar directamente el PDF: cualquier corrección debe realizarse "
        "sobre el archivo fuente correspondiente y regenerar el manual.}\n"
        "\\end{center}\n"
        "\\vspace*{\\fill}\n"
    )


def build_appendix_pendientes(items: list[tuple[str, str]]) -> str:
    """Construye el apéndice "Pendientes de desarrollo" agrupando por capítulo."""
    out: list[str] = []
    out.append("\\appendix\n")
    out.append("\\chapter{Pendientes de desarrollo}\n")
    out.append(
        "\\noindent Este apéndice consolida, en una sola vista, todos los "
        "elementos marcados a lo largo del manual como pendientes de "
        "desarrollo. Cada entrada indica el capítulo de procedencia y el "
        "texto del pendiente. La lista se genera de forma automática a "
        "partir de los marcadores \\texttt{[PENDIENTE: ...]} insertados en "
        "las fuentes; cualquier modificación se realiza sobre el archivo "
        "fuente correspondiente, no sobre este apéndice.\n\n"
    )

    # Agrupar conservando el orden de aparición.
    grouped: dict[str, list[str]] = {}
    order: list[str] = []
    for chap, text in items:
        if chap not in grouped:
            grouped[chap] = []
            order.append(chap)
        grouped[chap].append(text)

    for chap in order:
        out.append(f"\\section*{{{chap}}}\n")
        out.append(
            f"\\addcontentsline{{toc}}{{section}}{{{chap}}}\n"
        )
        out.append("\\begin{itemize}\n")
        for text in grouped[chap]:
            latex_text = md_to_latex(text).strip()
            out.append(f"  \\item \\pendiente{{{latex_text}}}\n")
        out.append("\\end{itemize}\n\n")

    out.append("\\newpage\n")
    return "\n".join(out)


def compile_pdf(tex_path: Path) -> bool:
    lualatex = shutil.which("lualatex")
    if lualatex is None:
        print("[!] lualatex no encontrado en PATH. Genera el .tex pero no compila.")
        return False
    cwd = tex_path.parent
    for i in range(2):
        result = subprocess.run(
            [lualatex, "-interaction=nonstopmode", "-halt-on-error", tex_path.name],
            cwd=cwd,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        if result.returncode != 0:
            print(f"[!] lualatex (pasada {i+1}) fallo:")
            tail = result.stdout[-3000:]
            try:
                print(tail)
            except UnicodeEncodeError:
                # Console encoding (cp1252 en Windows) puede no aceptar
                # algunos glifos Unicode del log; degradamos a ASCII.
                print(tail.encode("ascii", errors="replace").decode("ascii"))
            return False
    return True


def main() -> int:
    print(f"Leyendo fuentes desde: {SRC_DIR}")
    tex_content = assemble()
    OUT_TEX.write_text(tex_content, encoding="utf-8")
    print(f"  -> .tex escrito en: {OUT_TEX} ({len(tex_content):,} chars)")

    print("Compilando con lualatex (2 pasadas)...")
    ok = compile_pdf(OUT_TEX)
    if ok:
        print(f"  -> PDF generado: {OUT_PDF}")
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
