"""Genera Reference_manual.pdf a partir de docs/specs/*.md.

Prototipo autónomo (sin dependencias externas más allá de la stdlib).
Convierte un subconjunto controlado de Markdown a LaTeX y compila con
lualatex (motor con soporte nativo de Unicode). Si en el futuro pandoc
está disponible en el sistema, sustituir la función `md_to_latex` por una
invocación a pandoc — el resto (ensamblaje, preamble, agrupación) se
mantiene.

Uso:
    python manuals/build_reference_manual.py

Salida:
    manuals/Reference_manual.tex
    manuals/Reference_manual.pdf
"""
from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPECS_DIR = ROOT / "docs" / "specs"
SOURCES_DIR = ROOT / "manuals" / "sources"
OUT_DIR = ROOT / "manuals"
OUT_TEX = OUT_DIR / "Reference_manual.tex"
OUT_PDF = OUT_DIR / "Reference_manual.pdf"

GROUPS: list[tuple[str, list[str]]] = [
    ("Elementos 1D — Armaduras", ["Truss2D", "Truss2DCorot", "Truss3D", "Truss3DCorot"]),
    ("Elementos 1D — Cables", ["Cable2DCorot", "Cable3DCorot"]),
    ("Elementos 1D — Marcos / Vigas", ["Frame2DEuler", "Frame2DTimoshenko", "Frame2DEulerCorot", "Frame3D"]),
    ("Modelos Constitutivos (Materiales)", ["CableMaterial1D"]),
]

# Capítulos finales que NO derivan de specs (referencia técnica de plumbing
# arquitectural). Cada entrada: (título_capítulo, archivo_fuente).
APPENDIX_CHAPTERS: list[tuple[str, str]] = [
    ("Anexos técnicos", "anexo_capa_algebraica.md"),
]

UNICODE_MAP: dict[str, str] = {
    "σ": r"\ensuremath{\sigma}",
    "Σ": r"\ensuremath{\Sigma}",
    "ε": r"\ensuremath{\varepsilon}",
    "α": r"\ensuremath{\alpha}",
    "β": r"\ensuremath{\beta}",
    "γ": r"\ensuremath{\gamma}",
    "δ": r"\ensuremath{\delta}",
    "Δ": r"\ensuremath{\Delta}",
    "θ": r"\ensuremath{\theta}",
    "Θ": r"\ensuremath{\Theta}",
    "λ": r"\ensuremath{\lambda}",
    "Λ": r"\ensuremath{\Lambda}",
    "μ": r"\ensuremath{\mu}",
    "ν": r"\ensuremath{\nu}",
    "π": r"\ensuremath{\pi}",
    "ρ": r"\ensuremath{\rho}",
    "τ": r"\ensuremath{\tau}",
    "φ": r"\ensuremath{\varphi}",
    "ϕ": r"\ensuremath{\phi}",
    "Φ": r"\ensuremath{\Phi}",
    "χ": r"\ensuremath{\chi}",
    "ψ": r"\ensuremath{\psi}",
    "Ψ": r"\ensuremath{\Psi}",
    "ω": r"\ensuremath{\omega}",
    "Ω": r"\ensuremath{\Omega}",
    "κ": r"\ensuremath{\kappa}",
    "η": r"\ensuremath{\eta}",
    "ζ": r"\ensuremath{\zeta}",
    "ξ": r"\ensuremath{\xi}",
    "Ξ": r"\ensuremath{\Xi}",
    "→": r"\ensuremath{\to}",
    "←": r"\ensuremath{\leftarrow}",
    "↔": r"\ensuremath{\leftrightarrow}",
    "⇒": r"\ensuremath{\Rightarrow}",
    "⇐": r"\ensuremath{\Leftarrow}",
    "⇔": r"\ensuremath{\Leftrightarrow}",
    "≤": r"\ensuremath{\le}",
    "≥": r"\ensuremath{\ge}",
    "≠": r"\ensuremath{\neq}",
    "≈": r"\ensuremath{\approx}",
    "≲": r"\ensuremath{\lesssim}",
    "≳": r"\ensuremath{\gtrsim}",
    "∈": r"\ensuremath{\in}",
    "∂": r"\ensuremath{\partial}",
    "∇": r"\ensuremath{\nabla}",
    "∞": r"\ensuremath{\infty}",
    "∫": r"\ensuremath{\int}",
    "√": r"\ensuremath{\surd}",
    "·": r"\ensuremath{\cdot}",
    "×": r"\ensuremath{\times}",
    "±": r"\ensuremath{\pm}",
    "°": r"\ensuremath{^{\circ}}",
    "—": "---",
    "–": "--",
    "“": "``",
    "”": "''",
    "‘": "`",
    "’": "'",
    "…": r"\ldots{}",
    "ℝ": r"\ensuremath{\mathbb{R}}",
    "⟨": r"\ensuremath{\langle}",
    "⟩": r"\ensuremath{\rangle}",
    "‖": r"\ensuremath{\|}",
    "ʹ": "'",
    "′": r"\ensuremath{'}",
    "²": r"\ensuremath{^{2}}",
    "³": r"\ensuremath{^{3}}",
    "₀": r"\ensuremath{_{0}}",
    "₁": r"\ensuremath{_{1}}",
    "₂": r"\ensuremath{_{2}}",
    "₃": r"\ensuremath{_{3}}",
    "ₙ": r"\ensuremath{_{n}}",
    "ₜ": r"\ensuremath{_{t}}",
    "ᵀ": r"\ensuremath{^{\top}}",
    "⁺": r"\ensuremath{^{+}}",
    "⁻": r"\ensuremath{^{-}}",
}

LATEX_ESCAPE_TEXT = {
    "&": r"\&",
    "%": r"\%",
    "#": r"\#",
    "_": r"\_",
    "$": r"\$",
}


def _save(content: str, prefix: str, store: dict, counter: list[int]) -> str:
    key = f"@@{prefix}{counter[0]}@@"
    counter[0] += 1
    store[key] = content
    return key


def md_to_latex(md: str) -> str:
    """Convierte el subconjunto de Markdown usado en specs a LaTeX.

    Estrategia:
      1. Extraer bloques verbatim (código, math) a placeholders.
      2. Aplicar transformaciones MD → LaTeX al texto restante.
      3. Escapar caracteres LaTeX especiales en el texto.
      4. Restaurar los placeholders.
    """
    placeholders: dict[str, str] = {}
    counter = [0]

    # 1a. Bloques de código ```lang ... ```
    def _code_block(m: re.Match) -> str:
        lang = (m.group(1) or "").strip().lower()
        code = m.group(2)
        if lang in ("yaml", "yml"):
            ltx = f"\\begin{{lstlisting}}[language=yaml]\n{code}\n\\end{{lstlisting}}"
        elif lang in ("python", "py"):
            ltx = f"\\begin{{lstlisting}}[language=Python]\n{code}\n\\end{{lstlisting}}"
        elif lang == "bash":
            ltx = f"\\begin{{lstlisting}}[language=bash]\n{code}\n\\end{{lstlisting}}"
        else:
            ltx = f"\\begin{{lstlisting}}\n{code}\n\\end{{lstlisting}}"
        return _save(ltx, "CODE", placeholders, counter)

    md = re.sub(r"```([\w]*)\n(.*?)```", _code_block, md, flags=re.DOTALL)

    # 1b. Math display $$...$$
    def _math_display(m: re.Match) -> str:
        body = m.group(1).strip()
        return _save(f"\\[\n{body}\n\\]", "MDISP", placeholders, counter)

    md = re.sub(r"\$\$(.+?)\$\$", _math_display, md, flags=re.DOTALL)

    # 1c. Math inline $...$ (no salto de línea, no doble $)
    def _math_inline(m: re.Match) -> str:
        return _save(f"${m.group(1)}$", "MINL", placeholders, counter)

    md = re.sub(r"\$(?!\$)([^\n$]+?)\$", _math_inline, md)

    # 1d. Inline code `...`
    def _inline_code(m: re.Match) -> str:
        body = m.group(1)
        # Escapar dentro de \texttt{}
        body = body.replace("\\", r"\textbackslash{}")
        body = body.replace("{", r"\{").replace("}", r"\}")
        body = body.replace("_", r"\_").replace("&", r"\&")
        body = body.replace("#", r"\#").replace("%", r"\%")
        body = body.replace("$", r"\$")
        # Sustituir Unicode dentro del inline code (no llega la fase 4)
        for ch, cmd in UNICODE_MAP.items():
            body = body.replace(ch, cmd)
        return _save(f"\\texttt{{{body}}}", "ICODE", placeholders, counter)

    md = re.sub(r"`([^`\n]+?)`", _inline_code, md)

    # 2. Escapar caracteres LaTeX especiales en el texto restante (los
    #    placeholders @@...@@ no contienen ninguno, así que están a salvo).
    for ch, esc in LATEX_ESCAPE_TEXT.items():
        md = re.sub(
            rf"(?<!\\){re.escape(ch)}",
            esc.replace("\\", "\\\\"),
            md,
        )

    # 3. Transformaciones MD → LaTeX

    # Encabezados (# H1 lo descartamos: el título lo provee el ensamblador)
    # IMPORTANTE: el escape previo convirtió '#' → '\#'. Restauramos la sintaxis
    # MD para los encabezados: una línea que empiece por uno o más '\#' es un header.
    md = re.sub(r"^\\#\\#\\#\\#\s+(.+)$", r"\\subsubsection{\1}", md, flags=re.MULTILINE)
    md = re.sub(r"^\\#\\#\\#\s+(.+)$", r"\\subsection{\1}", md, flags=re.MULTILINE)
    md = re.sub(r"^\\#\\#\s+(.+)$", r"\\section{\1}", md, flags=re.MULTILINE)
    md = re.sub(r"^\\#\s+.+$", "", md, flags=re.MULTILINE)  # descartar H1

    # Negritas y cursivas
    md = re.sub(r"\*\*([^*\n]+?)\*\*", r"\\textbf{\1}", md)
    md = re.sub(r"(?<![*\\])\*([^*\n]+?)\*(?!\*)", r"\\textit{\1}", md)

    # Enlaces [texto](url) — preservamos solo el texto
    md = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r"\1", md)

    # Reglas horizontales
    md = re.sub(r"^---+\s*$", "", md, flags=re.MULTILINE)

    # Blockquotes simples (> texto)
    md = re.sub(r"^>\s*(.+)$", r"\\textit{\1}\n", md, flags=re.MULTILINE)

    # Listas con viñeta — bloque de líneas consecutivas iniciadas por `- ` o `* `
    def _itemize(match: re.Match) -> str:
        items = [ln for ln in match.group(0).split("\n") if ln.strip()]
        out = ["\\begin{itemize}"]
        for it in items:
            it = re.sub(r"^\s*[*-]\s+", "", it)
            out.append(f"  \\item {it}")
        out.append("\\end{itemize}\n")
        return "\n".join(out)

    md = re.sub(r"(?:^[ \t]*[*-]\s+.+\n?)+", _itemize, md, flags=re.MULTILINE)

    # Listas numeradas
    def _enumerate_list(match: re.Match) -> str:
        items = [ln for ln in match.group(0).split("\n") if ln.strip()]
        out = ["\\begin{enumerate}"]
        for it in items:
            it = re.sub(r"^\s*\d+\.\s+", "", it)
            out.append(f"  \\item {it}")
        out.append("\\end{enumerate}\n")
        return "\n".join(out)

    md = re.sub(r"(?:^[ \t]*\d+\.\s+.+\n?)+", _enumerate_list, md, flags=re.MULTILINE)

    # 4. Sustitución de Unicode matemático/Griego por comandos LaTeX
    #    (se hace ANTES de restaurar placeholders para no clobbear el contenido
    #    de bloques de código; los Unicode dentro de lstlisting los maneja la
    #    opción `literate` del preamble).
    for ch, cmd in UNICODE_MAP.items():
        md = md.replace(ch, cmd)

    # 5. Restaurar placeholders verbatim
    for key, val in placeholders.items():
        md = md.replace(key, val)

    return md


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
\usepackage{fancyhdr}
\usepackage{microtype}

\definecolor{yamlkey}{RGB}{0, 102, 204}
\definecolor{yamlval}{RGB}{204, 0, 0}
\definecolor{yamlcomment}{RGB}{120, 120, 120}
\definecolor{bglight}{RGB}{245, 245, 245}
\definecolor{darkgray}{RGB}{50, 50, 50}

\hypersetup{
    colorlinks=true,
    linkcolor=yamlkey,
    urlcolor=yamlkey,
    pdftitle={Manual de Referencia - Fenix FEM},
    pdfauthor={Jaime Retama Velasco}
}

\titleformat{\chapter}[display]
  {\normalfont\huge\bfseries\color{yamlkey}}{\chaptertitlename\ \thechapter}{20pt}{\Huge}
\titleformat{\section}{\Large\bfseries\color{darkgray}}{\thesection}{1em}{}
\titleformat{\subsection}{\large\bfseries\color{darkgray}}{\thesubsection}{1em}{}
\titleformat{\subsubsection}{\normalsize\bfseries\color{darkgray}}{\thesubsubsection}{1em}{}

\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\textbf{\color{darkgray}Fenix FEM --- Referencia}}
\fancyhead[R]{\color{darkgray}\leftmark}
\fancyfoot[C]{\thepage}
\renewcommand{\headrulewidth}{0.4pt}
\renewcommand{\footrulewidth}{0.4pt}

\lstdefinelanguage{yaml}{
  keywords={nodes, materials, elements, mesh, name, kind, status, interface, parameters, signature, conventions, validity, out_of_scope, acceptance, references, material_contract, dof_names, n_nodes, strain_dim, n_integration_points, type, required, desc, sign, voigt, node_orientation, configuration, expected_behaviour, strain_kind, primary_state_var, state_passthrough, numerical_caveats, setup, expect, tol_rel, tol_abs},
  keywordstyle=\color{yamlkey}\bfseries,
  ndkeywords={true, false, null},
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
           {≲}{{$\lesssim$}}1 {≳}{{$\gtrsim$}}1 {≈}{{$\approx$}}1
           {·}{{$\cdot$}}1 {×}{{$\times$}}1 {±}{{$\pm$}}1
           {°}{{$^{\circ}$}}1 {²}{{$^{2}$}}1 {³}{{$^{3}$}}1
           {₀}{{$_{0}$}}1 {₁}{{$_{1}$}}1 {₂}{{$_{2}$}}1
           {ᵀ}{{$^{\top}$}}1 {⁺}{{$^{+}$}}1 {⁻}{{$^{-}$}}1
           {⟨}{{$\langle$}}1 {⟩}{{$\rangle$}}1 {‖}{{$\|$}}1
           {—}{{---}}1 {–}{{--}}1 {−}{{-}}1
           {“}{{``}}1 {”}{{''}}1
           {ʹ}{{'}}1 {′}{{$'$}}1
           {…}{{\ldots{}}}1
           {ϕ}{{$\phi$}}1 {Σ}{{$\Sigma$}}1 {Λ}{{$\Lambda$}}1
           {η}{{$\eta$}}1 {ζ}{{$\zeta$}}1 {ξ}{{$\xi$}}1 {Ξ}{{$\Xi$}}1
           {ψ}{{$\psi$}}1 {Ψ}{{$\Psi$}}1 {Θ}{{$\Theta$}}1 {χ}{{$\chi$}}1
           {ℝ}{{$\mathbb{R}$}}1
           {∈}{{$\in$}}1 {∂}{{$\partial$}}1 {∇}{{$\nabla$}}1
           {∞}{{$\infty$}}1 {∫}{{$\int$}}1 {√}{{$\surd$}}1
           {₃}{{$_{3}$}}1 {ₙ}{{$_{n}$}}1 {ₜ}{{$_{t}$}}1
}

\begin{document}

\begin{titlepage}
    \centering
    \vspace*{2cm}
    {\Huge\bfseries\color{yamlkey} Fenix FEM \par}
    \vspace{1cm}
    {\LARGE Manual de Referencia \par}
    \vspace{0.5cm}
    \rule{\linewidth}{0.5mm} \par
    \vspace{2cm}
    {\large Especificaciones físicas, formulaciones numéricas y contratos\\
            de los componentes del programa\par}
    \vspace{2cm}
    {\large Generado automáticamente desde \texttt{docs/specs/}\par}
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

\chapter*{Sobre este Manual}
\addcontentsline{toc}{chapter}{Sobre este Manual}
\noindent Este manual contiene la \textbf{referencia formal} de los componentes de Fenix FEM: cada elemento finito y cada modelo constitutivo se documenta con su especificación física (ecuaciones), su formulación numérica (matrices $\mathbf B$, rigidez tangente, integración) y su contrato YAML.

El contenido se genera automáticamente desde los archivos \texttt{docs/specs/*.md} del repositorio, que constituyen la \emph{fuente única de verdad} de cada componente. No editar este PDF manualmente; cualquier corrección debe hacerse sobre la spec correspondiente y regenerar el manual mediante:
\begin{center}\texttt{python manuals/build\_reference\_manual.py}\end{center}

Para una guía orientada al uso del programa (sintaxis YAML, ejemplos, post-procesamiento), consulte el \textbf{Manual de Usuario} en \texttt{manuals/User\_manual.pdf}.

\newpage
\pagenumbering{arabic}
\setcounter{page}{1}

"""

POSTAMBLE = r"""
\end{document}
"""


def assemble() -> str:
    parts = [PREAMBLE]
    for chapter_name, components in GROUPS:
        parts.append(f"\\chapter{{{chapter_name}}}\n")
        for comp in components:
            spec_path = SPECS_DIR / f"{comp}.md"
            if not spec_path.exists():
                print(f"  [!] Spec no encontrada: {spec_path}")
                continue
            md = spec_path.read_text(encoding="utf-8")
            ltx = md_to_latex(md)
            parts.append(f"\\section{{{comp}}}\n")
            parts.append(ltx)
            parts.append("\n\\newpage\n")

    # Anexos técnicos (no derivan de specs; viven en manuals/sources/).
    for chapter_name, source_file in APPENDIX_CHAPTERS:
        source_path = SOURCES_DIR / source_file
        if not source_path.exists():
            print(f"  [!] Anexo no encontrado: {source_path}")
            continue
        md = source_path.read_text(encoding="utf-8")
        ltx = md_to_latex(md)
        parts.append(f"\\chapter{{{chapter_name}}}\n")
        parts.append(ltx)
        parts.append("\n\\newpage\n")

    parts.append(POSTAMBLE)
    return "\n".join(parts)


def compile_pdf(tex_path: Path) -> bool:
    lualatex = shutil.which("lualatex")
    if lualatex is None:
        print("[!] lualatex no encontrado en PATH. Genera el .tex pero no compila.")
        return False
    cwd = tex_path.parent
    for i in range(2):  # dos pasadas para resolver TOC
        result = subprocess.run(
            [lualatex, "-interaction=nonstopmode", "-halt-on-error", tex_path.name],
            cwd=cwd,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        if result.returncode != 0:
            print(f"[!] lualatex (pasada {i+1}) falló:")
            print(result.stdout[-3000:])
            return False
    return True


def main() -> int:
    print(f"Leyendo specs desde: {SPECS_DIR}")
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
