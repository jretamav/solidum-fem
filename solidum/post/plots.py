"""Figuras de los resultados de un modelo YAML (ADR 0039).

Se dibujan a partir de lo que Solidum escribe (ADR 0036 y 0038): los VTK por
instante con su ``.pvd``, el ``<nombre>_final.vtu`` y la historia en CSV. Así
las mismas figuras salen al terminar ``python -m solidum modelo.yaml`` o, sin
volver a calcular, con ``python -m solidum.post modelo.yaml``.

``output.plots`` es una lista de figuras:

- ``{type: mesh}``: la malla sin deformar, con los nodos con grados de
  libertad restringidos y las cargas nodales;
- ``{type: deformed, at: [f1, ...], field: Internal_State, scale: auto}``:
  la deformada en los instantes escritos más cercanos a cada factor de carga
  (``final`` para el último), coloreada por un campo del VTK (por celda, o
  ``displacement``: la magnitud del desplazamiento), un panel por instante con
  la misma escala de color;
- ``{type: curve, x: columna, y: columna | [columnas], compare: [otros.yaml]}``:
  columnas de la historia, también de otros modelos (estudios de malla).

Mallas 2D (sólidos planos y elementos de línea en el plano); las 3D, todavía
no. Colores (guía de visualización del proyecto): curvas con la paleta
categórica en orden fijo; campos positivos con una rampa azul de un solo
tono, campos con signo con la divergente azul–gris–rojo centrada en cero.
"""
from __future__ import annotations

import csv
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

from solidum.logging import get_logger

_log = get_logger("post")

PLOT_TYPES = ("mesh", "deformed", "curve")
_KEYS = {
    "mesh": {"type", "file", "title", "zoom", "markers"},
    "deformed": {"type", "file", "title", "zoom", "at", "field", "field_label", "scale", "limits",
                 "edges", "label_column", "show_undeformed"},
    "curve": {"type", "file", "title", "x", "y", "compare", "labels", "xlabel", "ylabel"},
}
_FORMATS = (".png", ".pdf", ".svg")

# Paleta de referencia de la guía de visualización (modo claro, impresión).
_SERIES = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
_BLUE_RAMP = ("#f7f9fc", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95",
              "#0d366b")
# Con más polígonos que esto, los bordes de los elementos tapan el campo.
_EDGES_AUTO_MAX = 3000
_DIVERGING = ("#184f95", "#3987e5", "#9ec5f4", "#f0efec", "#f3b0ad", "#e66767", "#b8302f")
_INK, _INK_2, _GRID = "#0b0b0b", "#52514e", "#d9d8d4"
_EDGE = "#3b3a37"

# Contorno de cada celda VTK de 2D en el orden de sus nodos (vértices y nodos
# intermedios alternados).
_OUTLINE = {
    "triangle": (0, 1, 2), "triangle6": (0, 3, 1, 4, 2, 5),
    "quad": (0, 1, 2, 3), "quad8": (0, 4, 1, 5, 2, 6, 3, 7), "quad9": (0, 4, 1, 5, 2, 6, 3, 7),
}
_LINES = {"line": (0, 1), "line3": (0, 2, 1)}
_CELLS_3D = {"tetra", "tetra10", "hexahedron", "hexahedron20", "hexahedron27", "wedge", "pyramid"}


# ----------------------------------------------------------------------------
# Validación
# ----------------------------------------------------------------------------

def validate_plots_block(plots) -> list[str]:
    """Errores de forma de ``output.plots`` (todos juntos)."""
    if plots is None:
        return []
    if not isinstance(plots, list) or not plots:
        return ["output.plots: se esperaba una lista de figuras ({type: mesh | deformed | curve, ...})."]
    errors = []
    files = set()
    for i, p in enumerate(plots):
        ctx = f"output.plots[{i}]"
        if not isinstance(p, dict) or p.get("type") not in PLOT_TYPES:
            errors.append(f"{ctx}: se esperaba un diccionario con type: {' | '.join(PLOT_TYPES)}.")
            continue
        kind = p["type"]
        extra = set(p) - _KEYS[kind]
        if extra:
            errors.append(f"{ctx}: claves no admitidas {sorted(extra)} en una figura '{kind}'; se "
                          f"admiten {sorted(_KEYS[kind])}.")
        f = p.get("file")
        if f is not None:
            if not isinstance(f, str) or not f:
                errors.append(f"{ctx}.file: se esperaba un nombre de archivo.")
            elif Path(f).suffix and Path(f).suffix.lower() not in _FORMATS:
                errors.append(f"{ctx}.file: formato {Path(f).suffix!r} no admitido; {list(_FORMATS)}.")
            elif f in files:
                errors.append(f"{ctx}.file: {f!r} repetido.")
            files.add(f)
        if "zoom" in p and not _numbers(p["zoom"], 4):
            errors.append(f"{ctx}.zoom: se esperaba [xmin, xmax, ymin, ymax].")
        if kind == "deformed":
            at = p.get("at", ["final"])
            if not isinstance(at, list) or not at or not all(
                    v == "final" or _is_number(v) for v in at):
                errors.append(f"{ctx}.at: se esperaba una lista de factores de carga (o 'final').")
            if p.get("edges", "auto") not in ("auto", True, False):
                errors.append(f"{ctx}.edges: se esperaba 'auto', true o false.")
            sc = p.get("scale", "auto")
            if not (sc == "auto" or (_is_number(sc) and sc >= 0)):
                errors.append(f"{ctx}.scale: se esperaba 'auto' o un número ≥ 0.")
            if "limits" in p and not _numbers(p["limits"], 2):
                errors.append(f"{ctx}.limits: se esperaba [mínimo, máximo].")
            if "field" in p and not isinstance(p["field"], str):
                errors.append(f"{ctx}.field: se esperaba el nombre de un campo.")
        if kind == "curve":
            if not isinstance(p.get("x"), str):
                errors.append(f"{ctx}.x: falta la columna del eje x.")
            ys = p.get("y")
            if not (isinstance(ys, str) or (isinstance(ys, list) and ys
                                           and all(isinstance(v, str) for v in ys))):
                errors.append(f"{ctx}.y: se esperaba una columna o una lista de columnas.")
            comp = p.get("compare", [])
            if not isinstance(comp, list) or not all(isinstance(v, str) for v in comp):
                errors.append(f"{ctx}.compare: se esperaba una lista de archivos YAML.")
            labels = p.get("labels")
            if labels is not None and not (isinstance(labels, list)
                                           and all(isinstance(v, str) for v in labels)):
                errors.append(f"{ctx}.labels: se esperaba una lista de nombres.")
    return errors


def _is_number(v) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def _numbers(v, n) -> bool:
    return isinstance(v, list) and len(v) == n and all(_is_number(x) for x in v)


# ----------------------------------------------------------------------------
# Lectura de la salida de Solidum
# ----------------------------------------------------------------------------

class ResultFiles:
    """Archivos de resultados de un modelo: ``<carpeta>/<nombre>*``."""

    def __init__(self, base_dir, base_name, history_file=None):
        self.base_dir = Path(base_dir)
        self.base_name = base_name
        self.history_path = self.base_dir / (history_file or f"{base_name}_historia.csv")

    @classmethod
    def of_yaml(cls, yaml_path) -> "ResultFiles":
        import yaml
        yaml_path = Path(yaml_path)
        with open(yaml_path, encoding="utf-8") as f:
            data = yaml.safe_load(f) or {}
        out = data.get("output") or {}
        hist = out.get("history") or {}
        return cls(yaml_path.parent, out.get("file_name") or yaml_path.stem, hist.get("file"))

    def instants(self) -> list[tuple[float, Path]]:
        """``[(factor, archivo)]`` de los VTK por instante, según el ``.pvd``."""
        pvd = self.base_dir / f"{self.base_name}.pvd"
        if not pvd.is_file():
            return []
        root = ET.parse(pvd).getroot()
        return [(float(ds.get("timestep")), self.base_dir / ds.get("file"))
                for ds in root.iter("DataSet")]

    def final(self) -> Path:
        return self.base_dir / f"{self.base_name}_final.vtu"

    def history(self) -> dict[str, np.ndarray]:
        if not self.history_path.is_file():
            raise FileNotFoundError(f"No existe la historia {self.history_path}: declare "
                                    "output.history y resuelva el modelo.")
        with open(self.history_path, newline="", encoding="utf-8") as f:
            rows = list(csv.reader(f))
        data = np.array(rows[1:], dtype=float) if len(rows) > 1 else np.zeros((0, len(rows[0])))
        return {name: data[:, k] for k, name in enumerate(rows[0])}


def _read_vtu(path: Path):
    import meshio
    if not path.is_file():
        raise FileNotFoundError(f"No existe {path}: resuelva el modelo con write_output.")
    m = meshio.read(str(path))
    kinds = {b.type for b in m.cells}
    if kinds & _CELLS_3D:
        raise ValueError(f"{path.name}: las figuras de mallas 3D todavía no están disponibles "
                         "(use los VTK).")
    return m


def _polygons(m, U=None, scale=0.0):
    """Polígonos de las celdas 2D, segmentos de las de línea y, por bloque, a
    cuál pertenece cada uno (para los campos por celda)."""
    pts = np.asarray(m.points[:, :2], dtype=float)
    if U is not None:
        pts = pts + scale * U[:, :2]
    polys, poly_src, segs, seg_src = [], [], [], []
    for b, block in enumerate(m.cells):
        if block.type in _OUTLINE:
            order = list(_OUTLINE[block.type])
            polys.extend(pts[block.data[:, order]])
            poly_src.extend((b, k) for k in range(len(block.data)))
        elif block.type in _LINES:
            segs.extend(pts[block.data[:, list(_LINES[block.type])]])
            seg_src.extend((b, k) for k in range(len(block.data)))
    return pts, polys, poly_src, segs, seg_src


def _cell_values(m, field: str, src) -> np.ndarray:
    if field == "displacement":
        mag = np.linalg.norm(np.asarray(m.point_data["Displacements"], dtype=float), axis=1)
        return np.array([mag[m.cells[b].data[k]].mean() for b, k in src])
    if field in m.cell_data:
        return np.array([float(m.cell_data[field][b][k]) for b, k in src])
    if field in m.point_data:
        v = np.asarray(m.point_data[field], dtype=float)
        v = np.linalg.norm(v, axis=1) if v.ndim == 2 else v
        return np.array([v[m.cells[b].data[k]].mean() for b, k in src])
    avail = sorted(set(m.cell_data) | set(m.point_data) | {"displacement"})
    raise ValueError(f"Campo {field!r} inexistente en el VTK; disponibles: {avail}.")


# ----------------------------------------------------------------------------
# Dibujo
# ----------------------------------------------------------------------------

def _pyplot():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - depende del entorno
        raise ImportError("Las figuras necesitan matplotlib (pip install matplotlib, o el grupo "
                          "opcional 'plot' de solidum).") from exc
    plt.rcParams.update({
        "font.size": 9, "axes.edgecolor": _INK_2, "axes.labelcolor": _INK, "xtick.color": _INK_2,
        "ytick.color": _INK_2, "axes.titlesize": 9, "axes.titlecolor": _INK,
        "figure.dpi": 100, "savefig.dpi": 200, "axes.spines.top": False,
        "axes.spines.right": False,
    })
    return plt


def _frame(ax, pts, zoom):
    ax.set_aspect("equal")
    if zoom:
        ax.set_xlim(zoom[0], zoom[1])
        ax.set_ylim(zoom[2], zoom[3])
    else:
        lo, hi = pts.min(axis=0), pts.max(axis=0)
        pad = 0.03 * float(np.max(hi - lo))
        ax.set_xlim(lo[0] - pad, hi[0] + pad)
        ax.set_ylim(lo[1] - pad, hi[1] + pad)
    ax.tick_params(length=2)


def _panel_size(pts, zoom, n_panels=1, width=7.0):
    if zoom:
        w, h = zoom[1] - zoom[0], zoom[3] - zoom[2]
    else:
        span = pts.max(axis=0) - pts.min(axis=0)
        w, h = float(span[0]) or 1.0, float(span[1]) or 1.0
    panel = float(np.clip(width * h / w, 1.1, 6.0)) + 0.5
    return width, min(panel * n_panels, 11.0)


def _output_path(files: ResultFiles, spec: dict, index: int) -> Path:
    name = spec.get("file") or f"{files.base_name}_{spec['type']}_{index + 1}"
    p = Path(name)
    if not p.suffix:
        p = p.with_suffix(".png")
    return p if p.is_absolute() else files.base_dir / p


def plot_mesh(files: ResultFiles, spec: dict, out: Path) -> Path:
    plt = _pyplot()
    from matplotlib.collections import LineCollection, PolyCollection
    pre = files.base_dir / f"{files.base_name}_pre_proceso.vtu"
    final = _read_vtu(files.final())
    m = _read_vtu(pre) if pre.is_file() else final
    pts, polys, _, segs, _ = _polygons(m)
    fig, ax = plt.subplots(figsize=_panel_size(pts, spec.get("zoom")))
    if polys:
        ax.add_collection(PolyCollection(polys, facecolors="#f6f5f2", edgecolors=_EDGE,
                                         linewidths=0.3))
    if segs:
        ax.add_collection(LineCollection(segs, colors=_SERIES[0], linewidths=1.5))
    if spec.get("markers", True):
        handles = []
        sup = np.asarray(final.point_data.get("Supports", np.zeros((len(pts), 3))),
                         dtype=float)[:, :2] != 0.0
        u_fin = np.asarray(final.point_data.get("Displacements", np.zeros((len(pts), 3))),
                           dtype=float)[:, :2]
        # Un grado de libertad restringido que se movió tiene desplazamiento impuesto.
        imposed = np.any(sup & (u_fin != 0.0), axis=1)
        fixed = np.where(np.any(sup, axis=1) & ~imposed)[0]
        if fixed.size:
            handles.append(ax.scatter(pts[fixed, 0], pts[fixed, 1], marker="^", s=36,
                                      color=_INK_2, zorder=3, label="apoyos"))
        moved = np.where(imposed)[0]
        if moved.size:
            handles.append(ax.scatter(pts[moved, 0], pts[moved, 1], marker="v", s=40,
                                      color=_SERIES[1], zorder=3,
                                      label="desplazamiento impuesto"))
        F = np.asarray(m.point_data.get("External_Forces", np.zeros((len(pts), 3))), dtype=float)[:, :2]
        loaded = np.where(np.linalg.norm(F, axis=1) > 0.0)[0]
        if loaded.size:
            size = float(np.max(pts.max(axis=0) - pts.min(axis=0)))
            d = F[loaded] / np.max(np.linalg.norm(F[loaded], axis=1)) * 0.06 * size
            ax.quiver(pts[loaded, 0] - d[:, 0], pts[loaded, 1] - d[:, 1], d[:, 0], d[:, 1],
                      angles="xy", scale_units="xy", scale=1.0, color=_SERIES[1], width=0.004,
                      zorder=4)
            handles.append(plt.Line2D([], [], color=_SERIES[1], lw=2, label="cargas nodales"))
        if handles:
            fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.0),
                       ncol=len(handles), frameon=False)
    _frame(ax, pts, spec.get("zoom"))
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    n_el = sum(len(b.data) for b in m.cells if b.type in _OUTLINE or b.type in _LINES)
    ax.set_title(spec.get("title") or f"Malla: {n_el} elementos, {len(pts)} nodos")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def _instant_files(files: ResultFiles, at) -> list[tuple[str, float | None, Path]]:
    inst = files.instants()
    out = []
    for f in at:
        if f == "final":
            out.append(("final", None, files.final()))
            continue
        if not inst:
            raise ValueError("No hay VTK por instante: declare output.results.frequency "
                             "({at: [...]}, {every: N} o all_steps).")
        k = int(np.argmin([abs(fac - float(f)) for fac, _ in inst]))
        out.append((f, inst[k][0], inst[k][1]))
    return out


def plot_deformed(files: ResultFiles, spec: dict, out: Path) -> Path:
    plt = _pyplot()
    from matplotlib.collections import LineCollection, PolyCollection
    from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
    chosen = _instant_files(files, spec.get("at", ["final"]))
    meshes = [_read_vtu(p) for _, _, p in chosen]
    field = spec.get("field")
    # Escala común: el mayor desplazamiento dibujado mide el 5 % del tamaño.
    ref = np.asarray(meshes[0].points[:, :2], dtype=float)
    size = float(np.max(ref.max(axis=0) - ref.min(axis=0)))
    umax = max(float(np.max(np.linalg.norm(np.asarray(m.point_data["Displacements"])[:, :2], axis=1)))
               for m in meshes)
    sc = spec.get("scale", "auto")
    scale = (0.05 * size / umax if umax > 0 else 1.0) if sc == "auto" else float(sc)
    geo = [_polygons(m, np.asarray(m.point_data["Displacements"], dtype=float), scale) for m in meshes]
    values = None
    norm, cmap = None, None
    if field:
        values = [_cell_values(m, field, g[2] + g[4]) for m, g in zip(meshes, geo)]
        allv = np.concatenate(values)
        lo, hi = spec.get("limits") or [float(np.nanmin(allv)), float(np.nanmax(allv))]
        if lo < 0.0 < hi:
            cmap = LinearSegmentedColormap.from_list("div", _DIVERGING)
            norm = TwoSlopeNorm(vcenter=0.0, vmin=lo, vmax=hi)
        else:
            cmap = LinearSegmentedColormap.from_list("seq", _BLUE_RAMP)
            norm = Normalize(vmin=lo, vmax=hi if hi > lo else lo + 1.0)
    hist = None
    col = spec.get("label_column")
    if col:
        hist = files.history()
        if col not in hist:
            raise ValueError(f"label_column {col!r} no es una columna de la historia "
                             f"({sorted(hist)}).")
    zoom = spec.get("zoom")
    w, h = _panel_size(ref, zoom, n_panels=len(chosen))
    fig, axes = plt.subplots(len(chosen), 1, figsize=(w, h), squeeze=False)
    edges = spec.get("edges", "auto")
    if edges == "auto":
        edges = max(len(g[1]) for g in geo) <= _EDGES_AUTO_MAX
    mappable = None
    for ax, (req, fac, path), m, g, k in zip(axes[:, 0], chosen, meshes, geo, range(len(chosen))):
        pts, polys, poly_src, segs, seg_src = g
        if spec.get("show_undeformed", False):
            und = _polygons(m)
            ax.add_collection(PolyCollection(und[1], facecolors="none", edgecolors=_GRID,
                                             linewidths=0.4))
        if polys:
            pc = PolyCollection(polys, edgecolors=_EDGE if edges else "face",
                                linewidths=0.2 if edges else 0.0)
            if values is not None:
                pc.set_array(values[k][:len(polys)])
                pc.set_cmap(cmap)
                pc.set_norm(norm)
                mappable = pc
            else:
                pc.set_facecolor("#e3eefb")
            ax.add_collection(pc)
        if segs:
            lc = LineCollection(segs, linewidths=1.5, colors=_SERIES[0])
            if values is not None:
                lc.set_array(values[k][len(polys):])
                lc.set_cmap(cmap)
                lc.set_norm(norm)
                mappable = mappable or lc
            ax.add_collection(lc)
        _frame(ax, pts, zoom)
        if req == "final":
            label = "final"
        else:
            label = f"factor {fac:.4g}"
        if hist is not None:
            step = _step_of(path)
            if step is None:
                row = len(hist["step"]) - 1
            else:
                rows = np.where(hist["step"] == step)[0]
                row = int(rows[0]) if rows.size else None
            if row is not None and row >= 0:
                label = f"{col} = {hist[col][row]:.4g}"
        ax.set_title(label, loc="left")
    if mappable is not None:
        cb = fig.colorbar(mappable, ax=list(axes[:, 0]), shrink=0.9, pad=0.02)
        cb.set_label(spec.get("field_label") or field)
        cb.outline.set_visible(False)
    fig.suptitle(spec.get("title") or (f"Deformada × {scale:.3g}"
                                       + (f", {spec.get('field_label') or field}" if field else "")),
                 color=_INK, fontsize=10)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def _step_of(path: Path):
    tail = path.stem.rsplit("_", 1)[-1]
    return int(tail) if tail.isdigit() else None


def plot_curve(files: ResultFiles, spec: dict, out: Path) -> Path:
    plt = _pyplot()
    x = spec["x"]
    ys = [spec["y"]] if isinstance(spec["y"], str) else list(spec["y"])
    runs = [(files.base_name, files)]
    for other in spec.get("compare", []):
        path = Path(other)
        if not path.is_absolute():
            path = files.base_dir / path
        runs.append((path.stem, ResultFiles.of_yaml(path)))
    series = [(name, run, y) for name, run in runs for y in ys]
    labels = spec.get("labels")
    if labels is not None and len(labels) != len(series):
        raise ValueError(f"labels: {len(labels)} nombres para {len(series)} curvas.")
    if len(series) > len(_SERIES):
        raise ValueError(f"{len(series)} curvas en una figura; como mucho {len(_SERIES)} "
                         "(divida la figura).")
    fig, ax = plt.subplots(figsize=(6.0, 3.8))
    for k, (name, run, y) in enumerate(series):
        h = run.history()
        for c in (x, y):
            if c not in h:
                raise ValueError(f"{run.history_path.name}: no tiene la columna {c!r} "
                                 f"({sorted(h)}).")
        label = labels[k] if labels else (name if len(ys) == 1 else
                                          (y if len(runs) == 1 else f"{name}: {y}"))
        ax.plot(h[x], h[y], color=_SERIES[k], lw=1.6, label=label)
    ax.set_xlabel(spec.get("xlabel") or x)
    ax.set_ylabel(spec.get("ylabel") or (ys[0] if len(ys) == 1 else ""))
    ax.grid(True, color=_GRID, lw=0.6)
    ax.set_axisbelow(True)
    if len(series) > 1:
        ax.legend(frameon=False)
    if spec.get("title"):
        ax.set_title(spec["title"])
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


_DRAW = {"mesh": plot_mesh, "deformed": plot_deformed, "curve": plot_curve}


def make_plots(files: ResultFiles, plots: list) -> list[Path]:
    """Dibuja las figuras de ``plots`` (``output.plots``) con los resultados de
    ``files``; devuelve las rutas escritas."""
    written = []
    for i, spec in enumerate(plots or []):
        out = _output_path(files, spec, i)
        out.parent.mkdir(parents=True, exist_ok=True)
        written.append(_DRAW[spec["type"]](files, spec, out))
        _log.info(f"Figura escrita: {out}")
    return written
