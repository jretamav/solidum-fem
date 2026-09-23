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

import os
import re
import shutil
import subprocess
import sys
import unicodedata
from pathlib import Path
from urllib.parse import quote

ROOT = Path(__file__).resolve().parents[1]

# El builder se invoca como script (`python manuals/build_reference_manual.py`),
# por lo que la raíz del repositorio no está en sys.path y `solidum.tools.spec`
# —del que ahora se deriva la agrupación de capítulos— no sería importable.
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

SPECS_DIR = ROOT / "docs" / "specs"
SOURCES_DIR = ROOT / "manuals" / "sources"
OUT_DIR = ROOT / "manuals"
OUT_TEX = OUT_DIR / "Reference_manual.tex"
OUT_PDF = OUT_DIR / "Reference_manual.pdf"

# ---------------------------------------------------------------------------
# Agrupación de specs en capítulos — DERIVADA, no enumerada.
#
# Historia: hasta 2026-08 este módulo llevaba una lista literal con los nombres
# de las specs que entraban al manual. Esa lista se quedó atrás respecto al
# repositorio (23 de 46 specs no llegaban al PDF: todo el subsistema 3D, la
# fractura embebida y 8 solvers) sin producir ningún error — el manual se
# anunciaba como "generado automáticamente desde docs/specs/" mientras omitía
# la mitad del catálogo en silencio. Reglas.md §1 prescribe justamente lo
# contrario: descubrimiento automático sobre listas manuales.
#
# Ahora la PERTENENCIA de cada spec a un capítulo se infiere de su propio
# contrato (`kind`, `interface.strain_dim`, `interface.dof_names`), que ya
# está validado por tests/test_specs.py. Lo único que queda declarado a mano
# es el ORDEN de presentación de los capítulos, que es criterio editorial y
# no un dato del dominio.
#
# Una spec que no encaje en ningún capítulo ABORTA el build nombrándola. El
# fallo silencioso era el defecto real; añadir un componente en el futuro debe
# propagarse solo o detenerse ruidosamente.
# ---------------------------------------------------------------------------

# Orden editorial de los capítulos. Cada entrada es la etiqueta que aparece
# como \chapter en el PDF; el orden de esta lista es el orden del manual.
CHAPTER_ORDER: list[str] = [
    "Elementos 1D — Armaduras",
    "Elementos 1D — Cables",
    "Elementos 1D — Marcos / Vigas",
    "Elementos 2D — Sólidos",
    "Elementos 2D — Discontinuidades embebidas",
    "Elementos 3D — Sólidos",
    "Elementos Térmicos — 2D",
    "Elementos Térmicos — 3D",
    "Modelos Constitutivos — 1D",
    "Modelos Constitutivos — 2D",
    "Modelos Constitutivos — 3D",
    "Modelos Constitutivos — Cohesivos",
    "Modelos Constitutivos — Térmicos",
    "Esquemas de Solución — Estáticos",
    "Esquemas de Solución — Modal y dinámicos",
]

# Solvers: la partición estático/dinámico no es inferible del contrato (no hay
# campo que la declare), así que se declara aquí por nombre. Es la única
# clasificación que permanece enumerada; cualquier solver no listado cae al
# capítulo de estáticos sólo si su spec lo permite — si no, el build aborta.
_STATIC_SOLVERS = {
    "LinearSolver",
    "NonlinearSolver",
    "ArcLengthSolver",
    "DissipationArcLengthSolver",
}


def _classify(spec) -> str:
    """Devuelve el capítulo al que pertenece una spec, o lanza SpecError.

    La clasificación usa exclusivamente datos ya presentes y validados en el
    contrato de la spec, de modo que un componente nuevo se ubica solo.
    """
    from solidum.tools.spec import SpecError

    kind = spec.kind
    iface = spec.contract.get("interface") or {}
    strain_dim = iface.get("strain_dim")
    dof_names = iface.get("dof_names") or []
    name = spec.name

    if kind == "cohesive_material":
        return "Modelos Constitutivos — Cohesivos"

    if kind == "thermal_material":
        return "Modelos Constitutivos — Térmicos"

    if kind == "material":
        by_dim = {
            1: "Modelos Constitutivos — 1D",
            3: "Modelos Constitutivos — 2D",
            6: "Modelos Constitutivos — 3D",
        }
        if strain_dim in by_dim:
            return by_dim[strain_dim]
        raise SpecError(
            f"{name}: material con strain_dim={strain_dim!r} no clasificable. "
            f"Esperado uno de {sorted(by_dim)}."
        )

    if kind == "solver":
        if name in _STATIC_SOLVERS:
            return "Esquemas de Solución — Estáticos"
        return "Esquemas de Solución — Modal y dinámicos"

    if kind == "element":
        # Los elementos térmicos se agrupan por su campo, no por strain_dim
        # (no tienen deformación); la dimensión la da `flux_dim`.
        if iface.get("field") == "temperature":
            flux_dim = iface.get("flux_dim")
            by_flux = {
                2: "Elementos Térmicos — 2D",
                3: "Elementos Térmicos — 3D",
            }
            if flux_dim in by_flux:
                return by_flux[flux_dim]
            raise SpecError(
                f"{name}: elemento térmico con flux_dim={flux_dim!r} no "
                f"clasificable. Esperado uno de {sorted(by_flux)}."
            )

        # Sólidos: la dimensión del tensor de deformaciones separa 2D de 3D.
        if strain_dim == 6:
            return "Elementos 3D — Sólidos"
        if strain_dim == 3:
            # El embebido se distingue por exigir DOS materiales: el bulk del
            # continuo y uno cohesivo que gobierna el salto en Gamma_d. Ese
            # segundo contrato sólo lo declaran los elementos con
            # discontinuidad interior, y ya está en el YAML de la spec.
            if "cohesive" in (spec.contract.get("material_contract") or {}):
                return "Elementos 2D — Discontinuidades embebidas"
            return "Elementos 2D — Sólidos"
        if strain_dim == 1:
            # Dentro de los 1D, las rotaciones nodales separan vigas de barras,
            # y el material unilateral separa el cable de la armadura.
            if any(d.startswith("r") for d in dof_names):
                return "Elementos 1D — Marcos / Vigas"
            if "cable" in name.lower():
                return "Elementos 1D — Cables"
            return "Elementos 1D — Armaduras"
        raise SpecError(
            f"{name}: elemento con strain_dim={strain_dim!r} no clasificable."
        )

    raise SpecError(f"{name}: kind={kind!r} sin capítulo asignado.")


def build_groups() -> list[tuple[str, list[str]]]:
    """Agrupa todas las specs del repositorio en capítulos ordenados.

    Recorre `docs/specs/` completo — no una lista literal — de modo que un
    componente nuevo entra al manual por el mero hecho de tener spec.
    """
    from solidum.tools.spec import SpecError, collect_specs, parse_spec

    buckets: dict[str, list[str]] = {ch: [] for ch in CHAPTER_ORDER}
    unclassified: list[str] = []

    for path in collect_specs(SPECS_DIR):
        spec = parse_spec(path)
        try:
            chapter = _classify(spec)
        except SpecError as exc:
            unclassified.append(str(exc))
            continue
        if chapter not in buckets:
            unclassified.append(
                f"{spec.name}: capítulo '{chapter}' no está en CHAPTER_ORDER."
            )
            continue
        buckets[chapter].append(spec.name)

    if unclassified:
        raise SystemExit(
            "[!] Specs no clasificables — el manual estaría incompleto:\n  "
            + "\n  ".join(unclassified)
        )

    # Orden estable dentro de cada capítulo: por número de nodos cuando el
    # contrato lo declara (Quad4 antes que Quad8), alfabético en el resto.
    def _sort_key(name: str):
        spec = parse_spec(SPECS_DIR / f"{name}.md")
        n_nodes = (spec.contract.get("interface") or {}).get("n_nodes")
        return (n_nodes if isinstance(n_nodes, int) else 0, name)

    return [(ch, sorted(buckets[ch], key=_sort_key)) for ch in CHAPTER_ORDER if buckets[ch]]


# Capítulos finales que NO derivan de specs (referencia técnica de plumbing
# arquitectural o catálogos transversales). Cada entrada: (título_capítulo,
# ruta del archivo fuente relativa a la raíz del repositorio).
APPENDIX_CHAPTERS: list[tuple[str, str]] = [
    ("Elementos — catálogo", "docs/catalogo_elementos.md"),
    ("Modelos constitutivos — catálogo", "docs/catalogo_materiales.md"),
    ("Solvers — catálogo", "docs/catalogo_solvers.md"),
    ("Anexos técnicos", "manuals/sources/anexo_capa_algebraica.md"),
]

UNICODE_MAP: dict[str, str] = {
    "✓": r"\ensuremath{\checkmark}",
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
    # `^` y `~` son activos también fuera de \texttt{}: en texto corriente
    # inician superíndice y espacio irrompible. Aparecen en encabezados que
    # citan símbolos sin delimitadores matemáticos (p. ej. "Matrices B
    # estándar + B^φ" en la spec del CST_Embedded2D), donde abortaban la
    # compilación con "Missing $ inserted".
    "^": r"\textasciicircum{}",
    "~": r"\textasciitilde{}",
    # `"`, `<` y `>` son atajos activos de babel en español ("a, "e, <<, >>).
    # Medido 2026-09-23: `frequency: "all_steps"` salía `.all_steps"` y
    # `"Elastic1D"` salía `.Elastic1D"` en los PDF, sin error de compilación;
    # dentro de \texttt, `{"ux": …}` abortaba ("Bad character code") y
    # `<carpeta>` pasaba el resto del párrafo a versalitas.
    '"': r"\textquotedbl{}",
    "<": r"\textless{}",
    ">": r"\textgreater{}",
}



# ----------------------------------------------------------------------
# Enlaces vivos (2026-09-23). Los manuales se leen SÓLO en pantalla, así que
# las referencias de las fuentes Markdown deben poder pulsarse. Hasta ahora el
# conversor descartaba todos los enlaces y dejaba sólo el texto (288 enlaces
# perdidos, más 431 menciones "ADR NNNN" sin vínculo en el Reference).
#
# Reglas:
#   - URL web              -> \href a la URL.
#   - Spec incluida en el manual que se compila -> enlace interno a su sección.
#   - ADR con sección en el manual (el de arquitectura) -> enlace interno.
#   - Cualquier otro archivo del repositorio -> su página en GitHub (rama
#     `main`, que muestra siempre la versión vigente).
#   - Ruta que no existe -> texto sin enlace y aviso al terminar el build.
# ----------------------------------------------------------------------
GITHUB_BLOB = "https://github.com/jretamav/solidum-fem/blob/main/"
GITHUB_TREE = "https://github.com/jretamav/solidum-fem/tree/main/"

_LINK_CTX: dict = {"source_dir": ROOT, "manual": None,
                   "internal_specs": set(), "internal_adrs": set()}
BROKEN_LINKS: list[tuple[str, str]] = []  # (archivo fuente, destino)
_ADR_FILES: dict[str, str] = {}


def set_link_context(*, source: Path | None = None, manual: str | None = None,
                     internal_specs=None, internal_adrs=None,
                     heading_offset: int | None = None,
                     label_prefix: str | None = None) -> None:
    """Fija el contexto con el que ``md_to_latex`` resuelve los enlaces.

    ``source`` es el archivo Markdown que se convierte (para resolver las rutas
    relativas); ``manual`` es ``"reference"``, ``"user"``, ``"architecture"`` o
    ``"examples"``; ``internal_specs`` / ``internal_adrs`` son los destinos que
    existen dentro del manual que se compila (el resto va a GitHub).
    ``label_prefix`` separa las etiquetas de figura de cada capítulo: dos
    ejemplos pueden llamar ``diagramas`` a su figura sin chocar."""
    if source is not None:
        _LINK_CTX["source_dir"] = Path(source).resolve().parent
        _LINK_CTX["source_name"] = Path(source).name
    if manual is not None:
        _LINK_CTX["manual"] = manual
    if internal_specs is not None:
        _LINK_CTX["internal_specs"] = set(internal_specs)
    if internal_adrs is not None:
        _LINK_CTX["internal_adrs"] = set(internal_adrs)
    if heading_offset is not None:
        _LINK_CTX["heading_offset"] = int(heading_offset)
    if label_prefix is not None:
        _LINK_CTX["label_prefix"] = label_prefix


def _adr_file(num: str) -> str | None:
    if not _ADR_FILES:
        for f in (ROOT / "docs" / "adr").glob("[0-9][0-9][0-9][0-9]-*.md"):
            _ADR_FILES[f.name[:4]] = f"docs/adr/{f.name}"
    return _ADR_FILES.get(num)


def _href(url: str, text: str) -> str:
    return "\\href{" + url.replace("%", "\\%").replace("#", "\\#") + "}{" + text + "}"


def _unescape_url(url: str) -> str:
    """Deshace el escape LaTeX que el conversor ya aplicó a la URL."""
    for esc, ch in (("\\textasciicircum{}", "^"), ("\\textasciitilde{}", "~"),
                    ("\\textquotedbl{}", '"'), ("\\textless{}", "<"), ("\\textgreater{}", ">"),
                    ("\\_", "_"), ("\\#", "#"), ("\\%", "%"), ("\\&", "&"), ("\\$", "$")):
        url = url.replace(esc, ch)
    return url.strip()


def _link_latex(text: str, url: str) -> str:
    raw = _unescape_url(url)
    if raw.startswith(("http://", "https://", "mailto:")):
        return _href(raw, text)
    path, _, frag = raw.partition("#")
    if not path:  # ancla dentro del mismo documento: sin destino fiable
        return text
    target = (_LINK_CTX["source_dir"] / path).resolve()
    try:
        rel = target.relative_to(ROOT)
    except ValueError:
        BROKEN_LINKS.append((_LINK_CTX.get("source_name", "?"), raw))
        return text
    if not target.exists():
        BROKEN_LINKS.append((_LINK_CTX.get("source_name", "?"), raw))
        return text
    rel_posix = rel.as_posix()
    manual = _LINK_CTX["manual"]
    if (rel_posix.startswith("docs/specs/") and target.suffix == ".md"
            and target.stem in _LINK_CTX["internal_specs"]):
        return f"\\hyperref[spec:{target.stem}]{{{text}}}"
    if rel_posix.startswith("docs/adr/") and target.name[:4] in _LINK_CTX["internal_adrs"]:
        return f"\\hyperref[adr:{target.name[:4]}]{{{text}}}"
    base = GITHUB_TREE if target.is_dir() else GITHUB_BLOB
    url_out = base + quote(rel_posix, safe="/-_.~")
    if frag:
        url_out += "#" + frag
    return _href(url_out, text)


_ADR_MENTION = re.compile(r"\bADR (\d{4})\b")
_NO_LINK_LINE = ("\\section", "\\subsection", "\\subsubsection", "\\chapter",
                 "\\label", "\\caption", "\\paragraph")


def _link_adr_mentions(md: str) -> str:
    """Convierte cada "ADR NNNN" del texto corriente en enlace: interno si el
    manual tiene la sección del ADR, a GitHub si no. No toca los títulos (un
    enlace dentro de un título rompe marcadores e índice) ni el código, las
    matemáticas o las tablas, que en esta fase ya son marcadores de posición."""
    internal = _LINK_CTX["internal_adrs"]

    def _repl(m: re.Match) -> str:
        num = m.group(1)
        if num in internal:
            return f"\\hyperref[adr:{num}]{{{m.group(0)}}}"
        f = _adr_file(num)
        return _href(GITHUB_BLOB + f, m.group(0)) if f else m.group(0)

    out = []
    for line in md.split("\n"):
        if line.lstrip().startswith(_NO_LINK_LINE):
            out.append(line)
        else:
            out.append(_ADR_MENTION.sub(_repl, line))
    return "\n".join(out)


def report_broken_links() -> None:
    if BROKEN_LINKS:
        print(f"  [!] {len(BROKEN_LINKS)} enlace(s) a rutas inexistentes (quedan como texto):")
        for src, dst in BROKEN_LINKS:
            print(f"      {src}: {dst}")


def _fig_label(name: str) -> str:
    prefix = _LINK_CTX.get("label_prefix") or ""
    return f"fig:{prefix}:{name}" if prefix else f"fig:{name}"


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
    #     Casos especiales:
    #       lang == "latex" o "tex": el contenido se inserta verbatim en el
    #         documento LaTeX, sin envoltorio `lstlisting`. Sirve para incluir
    #         diagramas TikZ, formulas extensas o cualquier fragmento LaTeX
    #         que deba renderizarse como tal y no como codigo fuente.
    def _code_block(m: re.Match) -> str:
        lang_raw = (m.group(1) or "").strip()
        lang = lang_raw.lower()
        code = m.group(2)
        # Bloques `callout Titulo`: el contenido se procesa recursivamente
        # como Markdown y se envuelve en el entorno `infobox` con el titulo
        # indicado tras la palabra clave. Solo disponible en manuales que
        # carguen el entorno (preambles de los builders de arquitectura y
        # de usuario; el de referencia no lo carga, por lo que los callouts
        # solo deben usarse en fuentes consumidas por aquellos).
        if lang.startswith("callout"):
            title = lang_raw[len("callout"):].strip() or "Nota"
            inner_md = code.strip()
            for key, val in placeholders.items():
                inner_md = inner_md.replace(key, val)
            inner_latex = md_to_latex(inner_md).strip()
            ltx = (
                f"\\begin{{infobox}}{{{title}}}\n"
                f"{inner_latex}\n"
                f"\\end{{infobox}}"
            )
            return _save(ltx, "CALLOUT", placeholders, counter)
        if lang in ("latex", "tex"):
            ltx = code  # se inserta tal cual al documento
        elif lang in ("yaml", "yml"):
            ltx = f"\\begin{{lstlisting}}[language=yaml]\n{code}\n\\end{{lstlisting}}"
        elif lang in ("python", "py"):
            ltx = f"\\begin{{lstlisting}}[language=Python]\n{code}\n\\end{{lstlisting}}"
        elif lang == "bash":
            ltx = f"\\begin{{lstlisting}}[language=bash]\n{code}\n\\end{{lstlisting}}"
        else:
            ltx = f"\\begin{{lstlisting}}\n{code}\n\\end{{lstlisting}}"
        return _save(ltx, "CODE", placeholders, counter)

    # Nota: la cabecera del fence acepta cualquier caracter no-newline para
    # soportar titulos con espacios (e.g. ``` callout Riesgo de arquitectura ```).
    # El uso original `[\w]*` solo capturaba la primera palabra y dejaba el
    # resto del titulo fuera del match, haciendo que la regex no encontrase
    # el bloque y el callout cayese al fallback de `lstlisting` generico.
    md = re.sub(r"```([^\n]*)\n(.*?)```", _code_block, md, flags=re.DOTALL)

    # 1a'. Figuras (2026-09-23, manual de ejemplos). Una línea sola
    #      ``![pie](ruta){#fig:nombre width=80%}`` — sintaxis de
    #      pandoc-crossref; el bloque {…} es opcional. La ruta es relativa a la
    #      fuente Markdown, y si apunta a un .png con un .pdf hermano se usa el
    #      PDF (vectorial, nítido a cualquier zoom en pantalla); así la misma
    #      fuente se ve en GitHub (PNG) y en el manual (PDF). Una figura que no
    #      existe aborta aquí: LaTeX abortaría igual, con un mensaje peor.
    def _figure(m: re.Match) -> str:
        caption, path, attrs = m.group(1).strip(), m.group(2), m.group(3) or ""
        target = (_LINK_CTX["source_dir"] / path).resolve()
        if target.suffix.lower() == ".png" and target.with_suffix(".pdf").exists():
            target = target.with_suffix(".pdf")
        if not target.exists():
            raise FileNotFoundError(
                f"md_to_latex: figura inexistente {path!r} en {_LINK_CTX.get('source_name', '?')}")
        rel = os.path.relpath(target, OUT_DIR).replace(os.sep, "/")
        width = re.search(r"width=(\d+(?:\.\d+)?)%", attrs)
        frac = float(width.group(1)) / 100 if width else 0.9
        label = re.search(r"#fig:([\w-]+)", attrs)
        # [H]: la figura va donde la cita el texto. En pantalla se lee en
        # scroll continuo, así que un hueco al pie de página molesta menos que
        # una figura que flota a la página siguiente.
        out = ["\\begin{figure}[H]", "\\centering",
               f"\\includegraphics[width={frac:g}\\linewidth]{{{rel}}}"]
        if caption:
            out.append(f"\\caption{{{md_to_latex(caption).strip()}}}")
        if label:
            out.append(f"\\label{{{_fig_label(label.group(1))}}}")
        out.append("\\end{figure}")
        return _save("\n".join(out), "FIG", placeholders, counter)

    md = re.sub(r"^!\[([^\]]*)\]\(([^)\s]+)\)(\{[^}\n]*\})?[ \t]*$", _figure, md,
                flags=re.MULTILINE)

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
        # Caracteres LaTeX activos que aparecen en código (e.g. B^T, ~user).
        # Dentro de \texttt{} `^` inicia superscript y `~` non-breaking space;
        # ambos rompen la compilación si llegan literales. Se traducen a las
        # macros \textasciicircum y \textasciitilde para que aparezcan como
        # glifos ASCII normales sin entrar en modo matemático.
        body = body.replace("^", r"\textasciicircum{}")
        body = body.replace("~", r"\textasciitilde{}")
        # Atajos activos de babel en español: ver LATEX_ESCAPE_TEXT.
        body = body.replace('"', r"\textquotedbl{}")
        body = body.replace("<", r"\textless{}").replace(">", r"\textgreater{}")
        # Sustituir Unicode dentro del inline code (no llega la fase 4)
        for ch, cmd in UNICODE_MAP.items():
            body = body.replace(ch, cmd)
        # Rutas y nombres largos (``solidum/math/linalg/iterative.py``,
        # ``ITERATIVE_MAX_RESTARTS``) no tienen dónde partirse y desbordaban
        # el margen. Se permite cortar tras '/', '.' y '_' sólo en los largos.
        # Umbral 16 (antes 22): varias rutas medianas seguidas
        # ("placa_agujero.msh", "resultados.json") también desbordaban.
        if len(m.group(1)) > 16:
            body = (body.replace("/", "/\\allowbreak{}")
                        .replace(".", ".\\allowbreak{}")
                        .replace("\\_", "\\_\\allowbreak{}"))
        return _save(f"\\texttt{{{body}}}", "ICODE", placeholders, counter)

    md = re.sub(r"`([^`\n]+?)`", _inline_code, md)

    # 1d'. Referencias a figura: ``@fig:nombre`` → "figura N" enlazada (con
    #      ``@Fig:`` → "Figura N", para principio de frase).
    def _fig_ref(m: re.Match) -> str:
        word = "Figura" if m.group(1) == "Fig" else "figura"
        label = _fig_label(m.group(2))
        return _save(f"\\hyperref[{label}]{{{word}~\\ref*{{{label}}}}}", "FREF",
                     placeholders, counter)

    md = re.sub(r"(?<![@\w])@(fig|Fig):([\w-]+)", _fig_ref, md)

    # 1e. Tablas estilo pipe Markdown.
    #     Patron: una linea de cabecera "| a | b |", separador "|---|---|"
    #     y una o mas filas de datos. Se convierten a `tabular` LaTeX y se
    #     guardan como placeholder para sobrevivir al escape posterior.
    def _table_block(m: re.Match) -> str:
        block = m.group(0).strip()
        lines = [ln.strip() for ln in block.split("\n") if ln.strip()]
        if len(lines) < 2:
            return block

        # Si la primera linea es un caption con formato "[TABLA: texto]",
        # se extrae para envolver la tabla en `\\begin{table}\\caption{...}`
        # y obtener numeracion automatica. Si no, la tabla se emite sin
        # numerar (mismo comportamiento previo).
        caption: str | None = None
        cap_match = re.match(r"^\[TABLA:\s*(.+?)\]$", lines[0])
        if cap_match:
            caption = cap_match.group(1).strip()
            lines = lines[1:]
            if len(lines) < 2:
                return block

        def _split_row(row: str) -> list[str]:
            row = row.strip()
            if row.startswith("|"):
                row = row[1:]
            if row.endswith("|"):
                row = row[:-1]
            return [c.strip() for c in row.split("|")]

        header = _split_row(lines[0])
        ncols = len(header)
        # Cada columna usa el tipo `Y` (X de tabularx con RaggedRight y
        # \footnotesize). Reparte \textwidth y permite el ajuste de linea
        # automatico, evitando que celdas con texto largo desborden el
        # ancho de la pagina.
        align_spec = "Y" * ncols

        # Convertir contenido de cada celda recursivamente (md -> latex) para
        # cubrir negritas, código inline y símbolos Unicode dentro de la tabla.
        # Antes de la llamada recursiva, restauramos los placeholders del
        # nivel exterior (@@ICODEn@@, @@MINLn@@, etc.) que ya pudieran
        # haberse insertado en el texto de la celda — de lo contrario, la
        # llamada interior no los reconoce y los emite literalmente al PDF.
        def _cell(text: str) -> str:
            for key, val in placeholders.items():
                text = text.replace(key, val)
            return md_to_latex(text).strip()

        body_rows = [_split_row(ln) for ln in lines[2:]]

        if caption is not None:
            # Tabla con caption: envoltorio `table` + `\caption` para
            # numeracion automatica y entrada en la lista de tablas.
            cap_latex = md_to_latex(caption).strip()
            out = [
                "\\begin{table}[h]",
                "\\centering",
                f"\\caption{{{cap_latex}}}",
                f"\\begin{{tabularx}}{{\\textwidth}}{{{align_spec}}}",
            ]
        else:
            out = [
                f"\\begin{{center}}\n"
                f"\\begin{{tabularx}}{{\\textwidth}}{{{align_spec}}}"
            ]
        out.append("\\hline")
        out.append(" & ".join(f"\\textbf{{{_cell(c)}}}" for c in header) + " \\\\")
        out.append("\\hline")
        for row in body_rows:
            cells = list(row) + [""] * (ncols - len(row))
            out.append(" & ".join(_cell(c) for c in cells[:ncols]) + " \\\\")
        out.append("\\hline")
        out.append("\\end{tabularx}")
        if caption is not None:
            out.append("\\end{table}")
        else:
            out.append("\\end{center}")
        latex = "\n".join(out)
        return _save(latex, "TABLE", placeholders, counter)

    md = re.sub(
        r"(?:^\[TABLA:[^\]]+\][ \t]*\n)?"    # caption opcional
        r"(?:^\|[^\n]+\|[ \t]*\n)"           # cabecera
        r"(?:\|[ \t:|\-]+\|[ \t]*\n)"         # separador --- | ---
        r"(?:\|[^\n]+\|[ \t]*\n?)+",           # filas
        _table_block,
        md,
        flags=re.MULTILINE,
    )

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
    # Nivel LaTeX de cada encabezado Markdown. `heading_offset` lo baja un
    # nivel cuando el documento ya vive dentro de una \section (las specs del
    # Reference manual): así sus partes cuelgan de la spec en el índice y en
    # los marcadores, en vez de quedar al mismo nivel que ella.
    #
    # Un título que ya trae su propio número ("1. Problema físico", "4.2 …")
    # se emite sin numerar por LaTeX, con entrada manual en el índice: LaTeX
    # le añadía otro número delante ("15.19.1" + "0." se leía "15.19.10.") y
    # las referencias del texto (§7) dejaban de coincidir con lo visible.
    _levels = ("section", "subsection", "subsubsection", "paragraph", "subparagraph")
    _offset = _LINK_CTX.get("heading_offset", 0)

    def _heading(m: re.Match) -> str:
        depth = len(m.group(1)) // 2 - 2          # '\#\#' -> 0, '\#\#\#' -> 1 ...
        cmd = _levels[min(depth + _offset, len(_levels) - 1)]
        title = m.group(2).strip()
        if re.match(r"^\d+(\.\d+)*\.?\s", title):
            return (f"\\phantomsection\n\\{cmd}*{{{title}}}\n"
                    f"\\addcontentsline{{toc}}{{{cmd}}}{{{title}}}")
        return f"\\{cmd}{{{title}}}"

    md = re.sub(r"^((?:\\#){2,4})\s+(.+)$", _heading, md, flags=re.MULTILINE)
    md = re.sub(r"^\\#\s+.+$", "", md, flags=re.MULTILINE)  # descartar H1

    # Negritas y cursivas
    md = re.sub(r"\*\*([^*\n]+?)\*\*", r"\\textbf{\1}", md)
    md = re.sub(r"(?<![*\\])\*([^*\n]+?)\*(?!\*)", r"\\textit{\1}", md)

    # Enlaces [texto](url): vivos (ver _link_latex). Se guardan como
    # marcadores para que el enlazado de "ADR NNNN" no los anide.
    def _md_link(m: re.Match) -> str:
        text = m.group(1)
        for ch, cmd in UNICODE_MAP.items():
            text = text.replace(ch, cmd)
        return _save(_link_latex(text, m.group(2)), "LINK", placeholders, counter)

    md = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", _md_link, md)
    md = _link_adr_mentions(md)

    # Reglas horizontales
    md = re.sub(r"^---+\s*$", "", md, flags=re.MULTILINE)

    # Blockquotes simples (> texto); el `>` ya llega escapado.
    md = re.sub(r"^\\textgreater\{\}\s*(.+)$", r"\\textit{\1}\n", md, flags=re.MULTILINE)

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

    # 5. Restaurar placeholders verbatim. Hay marcadores anidados (un enlace
    #    cuyo texto es código en línea: [`Tri3`](Tri3.md)), así que se repite
    #    hasta que no quede ninguno; un marcador superviviente aparecería
    #    literalmente en el PDF ("@@ICODE24@@"), de modo que se aborta.
    for _ in range(8):
        before = md
        for key, val in placeholders.items():
            md = md.replace(key, val)
        if md == before:
            break
    leftover = re.findall(r"@@[A-Z]+\d+@@", md)
    leftover = [k for k in leftover if k in placeholders]
    if leftover:
        raise RuntimeError(f"md_to_latex: marcadores sin restaurar {leftover[:5]}")

    return md


# Fuentes de los tres manuales (2026-09-23). Latin Modern sigue siendo la
# fuente principal, pero su variante monoespaciada no tiene griego ni la
# mayoría de símbolos matemáticos, y lualatex descartaba esos glifos EN
# SILENCIO: el Reference manual perdía ~6 200 caracteres (σ, ε, α, ≤, ⇒, ∈…),
# casi todos en los contratos YAML de las specs, donde el PDF mostraba un
# hueco en lugar de la letra. Con la cadena de respaldo de luaotfload, cada
# glifo que Latin Modern no tiene se toma de DejaVu (instalada con la
# distribución TeX) o, para dingbats como ✅, de Segoe UI Symbol/Emoji. El
# resto del texto no cambia de aspecto. Definido aquí una vez y reutilizado
# por los builders de User y Architecture.
FONT_SETUP = r"""\usepackage{fontspec}
% Respaldo de glifos: ver FONT_SETUP en build_reference_manual.py.
\directlua{
  luaotfload.add_fallback("solidumserif", {"DejaVuSerif:mode=node;", "DejaVuSans:mode=node;", "SegoeUISymbol:mode=node;", "SegoeUIEmoji:mode=node;"})
  luaotfload.add_fallback("solidummono", {"DejaVuSansMono:mode=node;", "DejaVuSans:mode=node;", "SegoeUISymbol:mode=node;", "SegoeUIEmoji:mode=node;"})
}
\setmainfont{Latin Modern Roman}[RawFeature={fallback=solidumserif}]
\setsansfont{Latin Modern Sans}[RawFeature={fallback=solidumserif}]
\setmonofont{Latin Modern Mono}[RawFeature={fallback=solidummono}]
"""



# Maquetación para pantalla (2026-09-23). El usuario fijó que los manuales se
# leen SÓLO en formato digital. Página más estrecha que carta con márgenes de
# pantalla y el MISMO ancho de texto (16,6 cm, que necesitan tablas y código):
# al ajustar al ancho de la ventana el texto se ve ~12 % mayor sin alargar las
# líneas. microtype (protrusión + expansión, completo sólo en LuaLaTeX y
# pdfLaTeX) y \emergencystretch evitan que el texto desborde el margen, que
# con márgenes estrechos llegaría al borde. El PDF abre con el panel de
# marcadores, ajustado al ancho, y cada página enlaza al índice.
SCREEN_GEOMETRY = (r"\usepackage[paperwidth=19.2cm, paperheight=25.6cm, hmargin=1.3cm, "
                   r"top=2.1cm, bottom=1.7cm, headheight=15pt, headsep=0.45cm, "
                   r"footskip=0.75cm]{geometry}")


def code_block_characters() -> list[str]:
    """Caracteres no ASCII que aparecen dentro de bloques ``` de todas las
    fuentes de los manuales (specs, catálogos, anexos, capítulos), más los
    archivos de ``examples/`` que el manual de ejemplos incrusta enteros como
    listado (``{{yaml:…}}``, ``{{py:…}}``)."""
    files = (list((ROOT / "docs" / "specs").glob("*.md"))
             + list((ROOT / "docs").glob("catalogo_*.md"))
             + list((ROOT / "manuals" / "sources").rglob("*.md")))
    embedded = (list((ROOT / "examples").glob("*/*.yaml"))
                + list((ROOT / "examples").glob("*/*.py")))
    chars: set[str] = set()
    blocks = [b for f in files
              for b in re.findall(r"```[^\n]*\n(.*?)```", f.read_text(encoding="utf-8"),
                                  flags=re.DOTALL)]
    blocks += [f.read_text(encoding="utf-8") for f in embedded]
    for block in blocks:
        # Sin marcas combinantes (categoría Mn): declaradas como un carácter
        # de una columna se separarían de la letra a la que acentúan.
        chars.update(c for c in block if ord(c) > 127 and not c.isspace()
                     and unicodedata.category(c) != "Mn")
    return sorted(chars)


def with_screen_setup(preamble: str, *, subject: str, keywords: str) -> str:
    """Adapta un preámbulo pensado para papel a lectura en pantalla."""
    for old, new in (
        ("\\documentclass[11pt,letterpaper,oneside]{report}", "\\documentclass[11pt,oneside]{report}"),
        ("\\usepackage[margin=2.5cm, headheight=15pt]{geometry}", SCREEN_GEOMETRY),
        # es-nodecimaldot: babel en español cambia el punto decimal por coma
        # sólo dentro de las fórmulas, y el texto, las tablas y el código usan
        # punto. Medido 2026-09-23: "0.78 %" en el texto y "1 + 0,78 (h/L)²"
        # en la ecuación del mismo párrafo.
        ("\\usepackage[spanish,es-tabla]{babel}",
         "\\usepackage[spanish,es-tabla,es-nodecimaldot]{babel}"),
    ):
        assert preamble.count(old) == 1, old
        preamble = preamble.replace(old, new, 1)
    mono = "\\setmonofont{Latin Modern Mono}[RawFeature={fallback=solidummono}]\n"
    assert preamble.count(mono) == 1
    # graphicx + float: las figuras que emite md_to_latex (`![…](…)`) usan
    # \includegraphics y la colocación [H]; cargados en todos los manuales
    # para que el conversor no dependa de qué preámbulo lo invoca.
    preamble = preamble.replace(
        mono, mono + "\\usepackage[protrusion=true,expansion=true]{microtype}\n"
              "\\usepackage{graphicx}\n\\usepackage{float}\n", 1)
    # Caracteres no ASCII en bloques de código: sin declararlos en `literate`,
    # listings no sabe que ocupan una columna y los recoloca (medido:
    # "(|ε| ≲ 1e-2)" salía "ε(|| ≲ 1e-2)"). Se declaran todos los que
    # aparecen en las fuentes, salvo los que el preámbulo ya declara.
    lit_start = preamble.find("literate=")
    assert lit_start >= 0
    declared = set(re.findall(r"\{(.)\}\{\{", preamble[lit_start:]))
    extra = [c for c in code_block_characters() if c not in declared]
    if extra:
        entries = " ".join("{" + c + "}{{" + c + "}}1" for c in extra)
        preamble = preamble.replace("literate=", "literate=" + entries + "\n           ", 1)
    screen = (
        "% --- Lectura en pantalla: ver with_screen_setup en build_reference_manual.py\n"
        "% Índice: columnas de número anchas (numeraciones de 2-3 niveles se\n"
        "% montaban sobre el título), índice compacto y marcadores profundos.\n"
        "\\makeatletter\n"
        "\\renewcommand*\\l@section{\\@dottedtocline{1}{1.5em}{3.3em}}\n"
        "\\renewcommand*\\l@subsection{\\@dottedtocline{2}{4.8em}{4.2em}}\n"
        "\\renewcommand*\\l@subsubsection{\\@dottedtocline{3}{9.0em}{5.0em}}\n"
        "\\makeatother\n"
        "\\setcounter{tocdepth}{2}\n"
        "\\definecolor{enlaceinterno}{RGB}{0, 84, 166}\n"
        "\\definecolor{enlaceexterno}{RGB}{0, 128, 110}\n"
        "\\hypersetup{\n"
        "    colorlinks=true,\n"
        "    linkcolor=enlaceinterno, citecolor=enlaceinterno,\n"
        "    urlcolor=enlaceexterno, filecolor=enlaceexterno,\n"
        "    bookmarksopen=true, bookmarksopenlevel=0, bookmarksnumbered=true, bookmarksdepth=3,\n"
        "    pdfpagemode=UseOutlines, pdfstartview=FitH, pdfdisplaydoctitle=true,\n"
        "    pdflang=es-ES,\n"
        f"    pdfsubject={{{subject}}},\n"
        f"    pdfkeywords={{{keywords}}},\n"
        "}\n"
        "\\setlength{\\emergencystretch}{3em}\n"
        "\\fancyfoot[L]{\\hyperlink{indice}{\\footnotesize Índice}}\n"
    )
    begin = "\\begin{document}\n"
    assert preamble.count(begin) == 1
    # El ancla del índice va como primera línea del .toc, así queda al
    # principio de la lista (bajo el título) y no en la página anterior.
    preamble = preamble.replace(
        begin, screen + begin + "\\addtocontents{toc}{\\protect\\hypertarget{indice}{}}\n", 1)
    return preamble


def with_font_setup(preamble: str) -> str:
    """Sustituye el ``\\usepackage{fontspec}`` del preámbulo por FONT_SETUP."""
    assert preamble.count("\\usepackage{fontspec}\n") == 1
    return preamble.replace("\\usepackage{fontspec}\n", FONT_SETUP, 1)


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
% Corchetes dobles \llbracket ⟦u⟧ \rrbracket: notación del salto de
% desplazamientos en las specs de discontinuidad embebida y material cohesivo
% (ADR 0010). No están en amssymb; sin stmaryrd la compilación aborta con
% "Undefined control sequence".
\usepackage{stmaryrd}
\usepackage{booktabs}
\usepackage{tabularx}
\usepackage{array}
\usepackage{ragged2e}
\usepackage{fancyhdr}
\usepackage{microtype}

\newcolumntype{Y}{>{\RaggedRight\arraybackslash\footnotesize}X}

\definecolor{yamlkey}{RGB}{0, 102, 204}
\definecolor{yamlval}{RGB}{204, 0, 0}
\definecolor{yamlcomment}{RGB}{120, 120, 120}
\definecolor{bglight}{RGB}{245, 245, 245}
\definecolor{darkgray}{RGB}{50, 50, 50}

\hypersetup{
    colorlinks=true,
    linkcolor=yamlkey,
    urlcolor=yamlkey,
    pdftitle={Manual de Referencia - Solidum FEM},
    pdfauthor={Jaime Retama Velasco}
}

\titleformat{\chapter}[display]
  {\normalfont\huge\bfseries\color{yamlkey}}{\chaptertitlename\ \thechapter}{20pt}{\Huge}
\titleformat{\section}{\Large\bfseries\color{darkgray}}{\thesection}{1em}{}
\titleformat{\subsection}{\large\bfseries\color{darkgray}}{\thesubsection}{1em}{}
\titleformat{\subsubsection}{\normalsize\bfseries\color{darkgray}}{\thesubsubsection}{1em}{}

\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\textbf{\color{darkgray}Solidum FEM --- Referencia}}
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
    {\Huge\bfseries\color{yamlkey} Solidum FEM \par}
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
\noindent Este manual contiene la \textbf{referencia formal} de los componentes de Solidum FEM: cada elemento finito y cada modelo constitutivo se documenta con su especificación física (ecuaciones), su formulación numérica (matrices $\mathbf B$, rigidez tangente, integración) y su contrato YAML.

El contenido se genera automáticamente desde los archivos \texttt{docs/specs/*.md} del repositorio, que constituyen la \emph{fuente única de verdad} de cada componente. No editar este PDF manualmente; cualquier corrección debe hacerse sobre la spec correspondiente y regenerar el manual mediante:
\begin{center}\texttt{python manuals/build\_reference\_manual.py}\end{center}

Para una guía orientada al uso del programa (sintaxis YAML, ejemplos, post-procesamiento), consulte el \textbf{Manual de Usuario} en \texttt{manuals/User\_manual.pdf}.

\newpage
\pagenumbering{arabic}
\setcounter{page}{1}

"""
PREAMBLE = with_font_setup(PREAMBLE)
PREAMBLE = with_screen_setup(
    PREAMBLE,
    subject="Especificación de cada elemento, material y solver de Solidum FEM",
    keywords="elementos finitos, mecánica de sólidos, especificaciones, Solidum FEM",
)

POSTAMBLE = r"""
\end{document}
"""


def _escape_title(text: str) -> str:
    """Escapa un nombre de componente para usarlo como título LaTeX.

    Los nombres de spec son identificadores de código (`CST_Embedded2D`), y el
    guión bajo es carácter activo en LaTeX: sin escapar aborta la compilación
    con "Missing $ inserted". Sólo aparece desde que el manual dejó de
    enumerar specs a mano y empezó a recorrer `docs/specs/` completo.
    """
    return text.replace("\\", r"\textbackslash{}").replace("_", r"\_")


def assemble() -> str:
    parts = [PREAMBLE]
    groups = build_groups()
    set_link_context(manual="reference", internal_adrs=set(),
                     internal_specs={c for _, comps in groups for c in comps
                                     if (SPECS_DIR / f"{c}.md").exists()})
    for chapter_name, components in groups:
        parts.append(f"\\chapter{{{_escape_title(chapter_name)}}}\n")
        for comp in components:
            spec_path = SPECS_DIR / f"{comp}.md"
            if not spec_path.exists():
                print(f"  [!] Spec no encontrada: {spec_path}")
                continue
            md = spec_path.read_text(encoding="utf-8")
            set_link_context(source=spec_path, heading_offset=1)
            ltx = md_to_latex(md)
            parts.append(f"\\section{{{_escape_title(comp)}}}\n\\label{{spec:{comp}}}\n")
            parts.append(ltx)
            parts.append("\n\\newpage\n")

    # Anexos (no derivan de specs; rutas relativas a la raíz del repo).
    for chapter_name, source_rel in APPENDIX_CHAPTERS:
        source_path = ROOT / source_rel
        if not source_path.exists():
            print(f"  [!] Anexo no encontrado: {source_path}")
            continue
        md = source_path.read_text(encoding="utf-8")
        set_link_context(source=source_path, heading_offset=0)
        ltx = md_to_latex(md)
        parts.append(f"\\chapter{{{_escape_title(chapter_name)}}}\n")
        parts.append(ltx)
        parts.append("\n\\newpage\n")

    parts.append(POSTAMBLE)
    report_broken_links()
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
