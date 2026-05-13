"""Genera User_manual.pdf — guía dirigida al usuario final de Fenix FEM.

Reutiliza la conversión Markdown→LaTeX de `build_reference_manual.py`. Las
fuentes son los archivos numerados en `manuals/sources/user/`, ensamblados en
orden alfabético como capítulos del manual.

Uso:
    python manuals/build_user_manual.py

Salida:
    manuals/User_manual.tex
    manuals/User_manual.pdf
"""
from __future__ import annotations

import datetime
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT / "manuals" / "sources" / "user"
OUT_DIR = ROOT / "manuals"
OUT_TEX = OUT_DIR / "User_manual.tex"
OUT_PDF = OUT_DIR / "User_manual.pdf"

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
\usepackage[most]{tcolorbox}

\newcolumntype{Y}{>{\RaggedRight\arraybackslash\footnotesize}X}

\definecolor{yamlkey}{RGB}{0, 102, 204}
\definecolor{yamlval}{RGB}{204, 0, 0}
\definecolor{yamlcomment}{RGB}{120, 120, 120}
\definecolor{bglight}{RGB}{245, 245, 245}
\definecolor{darkgray}{RGB}{50, 50, 50}
\definecolor{notebackground}{RGB}{240, 248, 255}
\definecolor{noteborder}{RGB}{0, 102, 204}
\definecolor{warnbackground}{RGB}{255, 245, 230}
\definecolor{warnborder}{RGB}{204, 102, 0}

% Cuadro de informacion: usado por los bloques `callout Nota`.
\newtcolorbox{infobox}[1]{
  enhanced,
  colback=notebackground, colframe=noteborder, colbacktitle=notebackground,
  coltitle=noteborder, fonttitle=\bfseries,
  title=#1, attach boxed title to top left={yshift=-2mm, xshift=4mm},
  boxed title style={colback=white, colframe=white, sharp corners},
  boxrule=0.4pt, arc=2pt,
  left=4mm, right=4mm, top=2mm, bottom=2mm,
  before skip=8pt, after skip=8pt,
}

% Cuadro de advertencia: usado por los bloques `callout Advertencia`.
% Comparte el comando ``\begin{infobox}{Advertencia}'' del conversor; los
% bloques marcados con titulo "Advertencia" reciben el color naranja via
% redefinicion local de las claves de color (ver macro \warnbox abajo).
\newtcolorbox{warnbox}[1]{
  enhanced,
  colback=warnbackground, colframe=warnborder, colbacktitle=warnbackground,
  coltitle=warnborder, fonttitle=\bfseries,
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
    pdftitle={Manual de Usuario - Fenix FEM},
    pdfauthor={Jaime Retama Velasco}
}

\titleformat{\chapter}[display]
  {\normalfont\huge\bfseries\color{yamlkey}}{\chaptertitlename\ \thechapter}{20pt}{\Huge}
\titleformat{\section}{\Large\bfseries\color{darkgray}}{\thesection}{1em}{}
\titleformat{\subsection}{\large\bfseries\color{darkgray}}{\thesubsection}{1em}{}
\titleformat{\subsubsection}{\normalsize\bfseries\color{darkgray}}{\thesubsubsection}{1em}{}

\renewcommand{\chaptermark}[1]{\markboth{Cap.\ \thechapter\ ---\ #1}{}}

\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\textbf{\color{darkgray}Fenix FEM --- Usuario}}
\fancyhead[R]{\color{darkgray}\nouppercase{\leftmark}}
\fancyfoot[C]{\thepage}
\fancyfoot[R]{\footnotesize\color{darkgray}FF-MU}
\renewcommand{\headrulewidth}{0.4pt}
\renewcommand{\footrulewidth}{0.4pt}

\lstdefinelanguage{yaml}{
  keywords={nodes, materials, elements, mesh, mesh_material, mesh_thickness, mesh_quadrature, mesh_physical_groups, quadrature, boundary_conditions, boundary_conditions_by_node, boundary_conditions_by_coord, boundary_conditions_by_group, point_loads, point_loads_by_node, point_loads_by_coord, point_loads_by_group, body_force, gravity, density, linear_constraints, slave, masters, coefficients, g, dof, node, solver, output, file_name, pre_process, export, results, text_export, frequency, nodal_results, element_results, format, id, type, coords, ref_vector, A, I, Iy, Iz, J, As, E, nu, sigma_y, H, hypothesis, ux, uy, uz, rx, ry, rz, tol, coord, val, num_steps, max_iter, adaptive, tol_iter, max_lambda, initial_dl, max_steps, material, thickness, kappa_0, alpha, cohesion, phi_deg, psi_deg, variant, n_modes, sigma, which, tolerance, lumping, t_end, dt, beta, gamma, rayleigh, alpha_rayleigh, beta_rayleigh, xi1, omega1, xi2, omega2, u0, u0_dot, F_func, freeze_tangent_after_iter, convergence, rtol_force, rtol_disp, atol_force_factor, atol_disp_factor, caption, linear_algebra, min_delta_lambda},
  keywordstyle=\color{yamlkey}\bfseries,
  ndkeywords={true, false, null, all, all_steps, last_step, plane_stress, plane_strain, plane_strain_matched, outer_cone, inner_cone, consistent, lumped, LinearSolver, NonlinearSolver, ArcLengthSolver, ModalSolver, NewmarkSolver, NewtonNewmarkSolver, Elastic1D, Elastic2D, Elastoplastic1D, VonMises2D, IsotropicDamage1D, IsotropicDamage2D, DruckerPrager2D, CableMaterial1D, Truss2D, Truss2DCorot, Truss3D, Truss3DCorot, Cable2DCorot, Cable3DCorot, Frame2DEuler, Frame2DTimoshenko, Frame2DEulerCorot, Frame3D, Quad4, Quad8, Quad9, Tri3, Tri6},
  ndkeywordstyle=\color{yamlval}\bfseries,
  identifierstyle=\color{black},
  sensitive=false,
  comment=[l]{\#},
  commentstyle=\color{yamlcomment}\ttfamily,
  stringstyle=\color{yamlval}\ttfamily,
  morestring=[b]',
  morestring=[b]"
}

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
           {ε}{{$\varepsilon$}}1 {σ}{{$\sigma$}}1 {α}{{$\alpha$}}1
           {β}{{$\beta$}}1 {γ}{{$\gamma$}}1 {δ}{{$\delta$}}1
           {θ}{{$\theta$}}1 {κ}{{$\kappa$}}1 {λ}{{$\lambda$}}1
           {μ}{{$\mu$}}1 {ν}{{$\nu$}}1 {π}{{$\pi$}}1
           {ρ}{{$\rho$}}1 {τ}{{$\tau$}}1 {φ}{{$\varphi$}}1
           {ω}{{$\omega$}}1 {Δ}{{$\Delta$}}1 {Ω}{{$\Omega$}}1
           {⇔}{{$\Leftrightarrow$}}1 {⇒}{{$\Rightarrow$}}1
           {→}{{$\to$}}1 {←}{{$\leftarrow$}}1
           {≤}{{$\le$}}1 {≥}{{$\ge$}}1 {≠}{{$\neq$}}1
           {·}{{$\cdot$}}1 {×}{{$\times$}}1 {±}{{$\pm$}}1
           {°}{{$^{\circ}$}}1 {²}{{$^{2}$}}1 {³}{{$^{3}$}}1
           {—}{{---}}1 {–}{{--}}1
           {“}{{``}}1 {”}{{''}}1
}

\begin{document}

\begin{titlepage}
    \centering
    \vspace*{2cm}
    {\Huge\bfseries\color{yamlkey} Fenix FEM \par}
    \vspace{1cm}
    {\LARGE Plataforma de Análisis por Elementos Finitos en Python \par}
    \vspace{0.3cm}
    {\normalsize\color{darkgray} Sigla: \texttt{FF-MU} \par}
    \vspace{0.5cm}
    \rule{\linewidth}{0.5mm} \par
    \vspace{2cm}
    {\huge\bfseries Manual de Usuario \par}
    \vspace{2cm}
    {\large Sintaxis del archivo YAML, catálogos de componentes,\\
            ejemplos completos y workflow de ejecución\par}
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
\noindent Este manual está dirigido al \textbf{usuario final} de Fenix FEM: cubre la sintaxis del archivo de entrada \texttt{.yaml}, el catálogo completo de elementos, materiales y solvers, ejemplos de uso y el workflow de post-procesamiento. No entra al detalle de la implementación interna ni a la formulación matemática completa de cada componente.

Para la \emph{referencia formal} de cada componente (ecuaciones, matrices $\mathbf B$, criterios de aceptación), consultar \texttt{manuals/Reference\_manual.pdf}, generado automáticamente desde \texttt{docs/specs/}.

Para la \emph{visión arquitectural} del programa (capas, bloques funcionales, mecanismos transversales y dirección de evolución), consultar \texttt{manuals/Architecture\_manual.pdf}.

Este manual se regenera con:
\begin{center}\texttt{python manuals/build\_user\_manual.py}\end{center}

\noindent No editar este PDF manualmente; cualquier corrección debe hacerse sobre los archivos fuente en \texttt{manuals/sources/user/}.

\newpage
\pagenumbering{arabic}
\setcounter{page}{1}

"""

POSTAMBLE = r"""
\end{document}
"""


def chapter_title(stem: str) -> str:
    """Convierte '01_introduccion' → 'Introducción' usando el primer H1 del archivo."""
    md = (SRC_DIR / f"{stem}.md").read_text(encoding="utf-8")
    for line in md.splitlines():
        if line.startswith("# "):
            return line[2:].strip()
    return stem


# Patron del bloque `callout Advertencia` tras la conversion del md_to_latex:
# el conversor envuelve siempre en `infobox`. Aqui re-etiquetamos las cajas
# con titulo "Advertencia" para que usen el entorno `warnbox` (naranja).
WARN_BOX_RE = re.compile(
    r"\\begin\{infobox\}\{Advertencia\}(.*?)\\end\{infobox\}",
    re.DOTALL,
)


def reroute_advertencia_boxes(ltx: str) -> str:
    """Convierte las cajas con titulo `Advertencia` al entorno `warnbox`."""
    return WARN_BOX_RE.sub(
        lambda m: f"\\begin{{warnbox}}{{Advertencia}}{m.group(1)}\\end{{warnbox}}",
        ltx,
    )


def assemble() -> str:
    parts = [PREAMBLE]
    sources = sorted(SRC_DIR.glob("*.md"))
    if not sources:
        raise SystemExit(f"[!] No hay fuentes en {SRC_DIR}")

    for idx, src in enumerate(sources, start=1):
        title = chapter_title(src.stem)
        md = src.read_text(encoding="utf-8")
        ltx = md_to_latex(md)
        ltx = reroute_advertencia_boxes(ltx)
        parts.append(f"\\chapter{{{title}}}\n\\label{{cap:{idx}}}\n")
        parts.append(ltx)
        parts.append("\n\\newpage\n")

    parts.append(build_colofon())
    return "\n".join(parts) + POSTAMBLE


def get_git_commit() -> str:
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
    """Página final con metadatos: commit, fecha, sigla."""
    commit = get_git_commit()
    dirty = " (con cambios sin confirmar)" if get_git_dirty() else ""
    fecha = datetime.datetime.now().strftime("%d de %B de %Y, %H:%M")
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
        "{\\bfseries\\color{yamlkey} Fenix FEM --- Manual de Usuario}\\\\[0.5em]\n"
        f"Sigla del manual: \\texttt{{FF-MU}}\\\\[0.5em]\n"
        f"Compilado el {fecha}\\\\[0.5em]\n"
        f"Commit Git: \\texttt{{{commit}}}{dirty}\\\\[1em]\n"
        "\\rule{0.4\\textwidth}{0.4pt}\\\\[2em]\n"
        "\\parbox{0.75\\textwidth}{\\small\\itshape\\centering\n"
        "Documento generado automáticamente desde los archivos fuente en "
        "\\texttt{manuals/sources/user/} mediante "
        "\\texttt{python manuals/build\\_user\\_manual.py}. "
        "No editar directamente el PDF: cualquier corrección debe realizarse "
        "sobre el archivo fuente correspondiente y regenerar el manual.}\n"
        "\\end{center}\n"
        "\\vspace*{\\fill}\n"
    )


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
