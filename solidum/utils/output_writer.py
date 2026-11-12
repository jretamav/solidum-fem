"""Escritura de los resultados declarados en el bloque ``output`` del YAML (ADR 0036).

:class:`OutputWriter` se engancha al ``step_callback`` del análisis y escribe,
junto al YAML:

- el modelo sin deformar (``pre_process.export``): ``<nombre>_pre_proceso.vtu``;
- los resultados por paso en VTK según ``results.frequency``:
  ``"last_step"`` (por omisión, sólo el final), ``"all_steps"``,
  ``{every: N}`` (cada N pasos) o ``{at: [f1, f2, ...]}`` (el primer paso
  convergido con factor de carga, o tiempo, ≥ cada valor); con más de un
  paso, un ``<nombre>.pvd`` para ParaView; siempre ``<nombre>_final.vtu``;
- la **historia** (``history``): un CSV con una fila por paso convergido,
  ``step``, ``factor`` y las columnas declaradas (desplazamiento medio, mínimo
  o máximo de un nodo, grupo físico o coordenada; reacción sumada de un grupo
  de apoyos o de desplazamiento impuesto). Se escribe fila a fila, así que un
  análisis interrumpido conserva lo calculado;
- la salida en texto estilo FEAP (``text_export``).

La reacción de un paso es ``F_int(U) − factor·F`` en los grados de libertad
prescritos (sin restricciones multipunto de apoyo), la del ADR 0004 §3; sólo
en los análisis estáticos.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from solidum.logging import get_logger

_log = get_logger("output")

_FREQ_WORDS = ("last_step", "all_steps")
_REDUCE = {"mean": np.mean, "min": np.min, "max": np.max, "sum": np.sum}
_COLUMN_KINDS = ("displacement", "reaction")
_SELECTORS = ("node", "nodes", "group", "at")


def validate_output_block(output) -> list[str]:
    """Errores de forma de ``output.results.frequency`` y ``output.history``."""
    errors = []
    if not isinstance(output, dict):
        return errors
    res = output.get("results") or {}
    if isinstance(res, dict) and "frequency" in res:
        freq = res["frequency"]
        if isinstance(freq, str):
            if freq not in _FREQ_WORDS:
                errors.append(f"output.results.frequency: {freq!r} no admitido; se admiten "
                              f"{list(_FREQ_WORDS)}, {{every: N}} o {{at: [f1, f2, ...]}}.")
        elif isinstance(freq, dict) and len(freq) == 1 and "every" in freq:
            n = freq["every"]
            if isinstance(n, bool) or not isinstance(n, int) or n < 1:
                errors.append(f"output.results.frequency.every: se esperaba un entero ≥ 1, no {n!r}.")
        elif isinstance(freq, dict) and len(freq) == 1 and "at" in freq:
            at = freq["at"]
            if not (isinstance(at, list) and at and all(
                    isinstance(v, (int, float)) and not isinstance(v, bool) for v in at)):
                errors.append("output.results.frequency.at: se esperaba una lista de números.")
        else:
            errors.append("output.results.frequency: se admite 'last_step', 'all_steps', "
                          "{every: N} o {at: [f1, f2, ...]}.")
    hist = output.get("history")
    if hist is None:
        return errors
    if not isinstance(hist, dict) or not isinstance(hist.get("columns"), list) or not hist["columns"]:
        errors.append("output.history: se esperaba {columns: [...]} (y opcionalmente 'file').")
        return errors
    unknown = set(hist) - {"columns", "file"}
    if unknown:
        errors.append(f"output.history: claves no admitidas {sorted(unknown)}; se admiten 'columns' y 'file'.")
    names = set()
    for i, col in enumerate(hist["columns"]):
        ctx = f"output.history.columns[{i}]"
        if not isinstance(col, dict):
            errors.append(f"{ctx}: se esperaba un diccionario.")
            continue
        name = col.get("name")
        if not isinstance(name, str) or not name or name in ("step", "factor") or name in names:
            errors.append(f"{ctx}.name: se esperaba un nombre único (no 'step' ni 'factor').")
        names.add(name)
        kinds = [k for k in _COLUMN_KINDS if k in col]
        if len(kinds) != 1:
            errors.append(f"{ctx}: declare exactamente una magnitud: {list(_COLUMN_KINDS)}.")
            continue
        extra = set(col) - {"name", "scale", kinds[0]}
        if extra:
            errors.append(f"{ctx}: claves no admitidas {sorted(extra)}.")
        if "scale" in col and (isinstance(col["scale"], bool) or not isinstance(col["scale"], (int, float))):
            errors.append(f"{ctx}.scale: se esperaba un número.")
        sel = col[kinds[0]]
        if not isinstance(sel, dict):
            errors.append(f"{ctx}.{kinds[0]}: se esperaba un diccionario con 'dof' y uno de {list(_SELECTORS)}.")
            continue
        found = [s for s in _SELECTORS if s in sel]
        if len(found) != 1:
            errors.append(f"{ctx}.{kinds[0]}: declare exactamente uno de {list(_SELECTORS)}.")
        if not isinstance(sel.get("dof"), str):
            errors.append(f"{ctx}.{kinds[0]}.dof: falta el grado de libertad (p. ej. 'uy').")
        allowed = {"dof", "tol", *_SELECTORS} | ({"reduce"} if kinds[0] == "displacement" else set())
        if set(sel) - allowed:
            errors.append(f"{ctx}.{kinds[0]}: claves no admitidas {sorted(set(sel) - allowed)}.")
        if "reduce" in sel and sel["reduce"] not in ("mean", "min", "max"):
            errors.append(f"{ctx}.{kinds[0]}.reduce: se admite 'mean', 'min' o 'max'.")
    return errors


class OutputWriter:
    """Escribe los resultados declarados en ``output`` mientras corre el análisis.

    Parameters
    ----------
    domain, assembler
        Modelo y ensamblador del análisis.
    output_config : dict
        Bloque ``output`` del YAML.
    yaml_path : str | Path
        Archivo del modelo: los resultados van a su carpeta, con
        ``output.file_name`` (por omisión, el nombre del YAML) como base.
    F_ext : ndarray
        Carga de referencia del análisis (la que escala el factor de carga).
    static : bool
        ``True`` si el análisis es estático (las reacciones de la historia
        sólo se definen ahí).
    """

    def __init__(self, domain, assembler, output_config: dict, yaml_path, F_ext, *, static: bool = True):
        self.domain, self.assembler = domain, assembler
        cfg = output_config or {}
        yaml_path = Path(yaml_path)
        self.base_dir = yaml_path.parent
        self.base_name = cfg.get("file_name") or yaml_path.stem
        self.base = str(self.base_dir / self.base_name)
        self.F_ext = np.zeros(domain.total_dofs) if F_ext is None else np.asarray(F_ext, dtype=float)
        res = cfg.get("results") or {}
        self.freq = res.get("frequency", "last_step")
        self.n_vars = res.get("nodal_results")
        self.e_vars = res.get("element_results")
        self.pre = bool((cfg.get("pre_process") or {}).get("export", False))
        self.text_cfg = cfg.get("text_export") or {}
        self.export_text = bool(self.text_cfg.get("export", False))
        self.n_steps = 0
        self.last = None
        self.pvd: list = []
        self.text_data: list = []
        self._at_pending = sorted(float(v) for v in self.freq["at"]) if isinstance(self.freq, dict) and "at" in self.freq else []
        self.columns = self._resolve_columns(cfg.get("history"), static)
        self.history_path = None
        if self.columns:
            hist = cfg["history"]
            self.history_path = self.base_dir / (hist.get("file") or f"{self.base_name}_historia.csv")
            with open(self.history_path, "w", newline="", encoding="utf-8") as f:
                csv.writer(f).writerow(["step", "factor", *[c["name"] for c in self.columns]])

    # ------------------------------------------------------------------
    def _nodes_of(self, sel: dict, ctx: str):
        dom = self.domain
        if "node" in sel:
            ids = [sel["node"]]
        elif "nodes" in sel:
            ids = list(sel["nodes"])
        elif "group" in sel:
            groups = getattr(dom, "physical_groups", None) or {}
            if sel["group"] not in groups:
                raise ValueError(f"{ctx}: grupo físico {sel['group']!r} inexistente. Disponibles: {sorted(groups)}.")
            ids = list(groups[sel["group"]])
        else:
            p = np.asarray(sel["at"], dtype=float)
            tol = float(sel.get("tol", 1e-6))
            ids = [nid for nid, n in dom.nodes.items()
                   if len(n.coordinates) >= len(p)
                   and np.linalg.norm(np.asarray(n.coordinates[:len(p)], dtype=float) - p) <= tol]
            if not ids:
                raise ValueError(f"{ctx}: ningún nodo a distancia ≤ {tol} de {list(p)}.")
        missing = [i for i in ids if i not in dom.nodes]
        if missing:
            raise ValueError(f"{ctx}: nodos inexistentes {missing}.")
        return [dom.nodes[i] for i in ids]

    def _resolve_columns(self, hist, static):
        if not hist:
            return []
        cols = []
        for i, col in enumerate(hist["columns"]):
            kind = "displacement" if "displacement" in col else "reaction"
            ctx = f"output.history.columns[{i}] ({col['name']})"
            if kind == "reaction" and not static:
                raise ValueError(f"{ctx}: la reacción sólo se registra en análisis estáticos.")
            sel = col[kind]
            nodes = self._nodes_of(sel, ctx)
            dof = sel["dof"]
            lacking = [n.id for n in nodes if dof not in n.dofs]
            if lacking:
                raise ValueError(f"{ctx}: los nodos {lacking[:10]} no tienen el grado de libertad {dof!r}.")
            dofs = np.array([n.dofs[dof] for n in nodes], dtype=int)
            if kind == "reaction":
                free = [n.id for n in nodes if dof not in n.boundary_conditions]
                if free:
                    raise ValueError(f"{ctx}: la reacción necesita '{dof}' prescrito en todos los nodos "
                                     f"del selector; libres: {free[:10]}.")
            cols.append({"name": col["name"], "kind": kind, "dofs": dofs,
                         "reduce": _REDUCE[sel.get("reduce", "mean") if kind == "displacement" else "sum"],
                         "scale": float(col.get("scale", 1.0))})
        return cols

    # ------------------------------------------------------------------
    def start(self):
        if self.pre:
            from solidum.utils.vtk_exporter import VtkExporter
            ruta = self.base + "_pre_proceso.vtu"
            VtkExporter(self.domain).export(ruta, U=None, F_ext=self.F_ext, nodal_vars=self.n_vars,
                                            elem_vars=self.e_vars)
            _log.info(f"Modelo sin deformar exportado a {ruta}")

    def _want_vtk(self, step, factor) -> bool:
        f = self.freq
        if f == "all_steps":
            return True
        if isinstance(f, dict) and "every" in f:
            return step % int(f["every"]) == 0
        if self._at_pending and factor >= self._at_pending[0] - 1e-12 * max(1.0, abs(factor)):
            while self._at_pending and factor >= self._at_pending[0] - 1e-12 * max(1.0, abs(factor)):
                self._at_pending.pop(0)
            return True
        return False

    def step(self, step, U, factor):
        """Callback de paso convergido: ``(paso, U, factor de carga o tiempo)``."""
        self.n_steps += 1
        n = self.n_steps
        factor = float(factor)
        self.last = (n, np.array(U, copy=True), factor)
        if self._want_vtk(n, factor):
            from solidum.utils.vtk_exporter import VtkExporter
            name = f"{self.base_name}_{n:04d}.vtu"
            VtkExporter(self.domain).export(str(self.base_dir / name), U, self.F_ext * factor,
                                            self.n_vars, self.e_vars)
            self.pvd.append((factor, name))
        if self.columns:
            row = [n, factor]
            F_int = None
            for c in self.columns:
                if c["kind"] == "displacement":
                    val = c["reduce"](U[c["dofs"]])
                else:
                    if F_int is None:
                        _, F_int = self.assembler.assemble_non_linear_system(U)
                    val = np.sum(F_int[c["dofs"]] - factor * self.F_ext[c["dofs"]])
                row.append(f"{c['scale'] * float(val):.10g}")
            with open(self.history_path, "a", newline="", encoding="utf-8") as f:
                csv.writer(f).writerow(row)
        if self.export_text:
            self.text_data.append(self._text_step(n, U, factor))

    def finish(self, U_final=None):
        """Cierra la salida: resultado final, archivo PVD y texto."""
        from solidum.utils.vtk_exporter import VtkExporter
        if self.n_steps == 0 and U_final is not None:
            # Análisis sin pasos (lineal): un único estado con factor 1.
            self.step(1, U_final, 1.0)
        if self.last is None:
            return
        _, U, factor = self.last
        if len(self.pvd) > 1:
            with open(self.base + ".pvd", "w", encoding="utf-8") as f:
                f.write('<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" '
                        'byte_order="LittleEndian">\n<Collection>\n')
                for fac, name in self.pvd:
                    f.write(f'  <DataSet timestep="{fac:.6f}" group="" part="0" file="{name}"/>\n')
                f.write("</Collection>\n</VTKFile>\n")
        VtkExporter(self.domain).export(self.base + "_final.vtu", U, self.F_ext * factor, self.n_vars, self.e_vars)
        if self.export_text:
            self._write_text()
        _log.info(f"Resultados escritos en {self.base_dir} ({self.base_name}*)")

    # ------------------------------------------------------------------
    # Salida de texto estilo FEAP (trasladada de examples/ejecutar_yaml.py)
    def _sel(self, key, pool):
        v = self.text_cfg.get(key, [])
        return list(pool) if isinstance(v, str) and v.lower() == "all" else list(v)

    def _text_step(self, n, U, factor):
        from solidum.utils.vtk_exporter import element_average_stress, von_mises
        n_vars = self.text_cfg.get("nodal_results", [])
        e_vars = self.text_cfg.get("element_results", [])
        nodes = {}
        for nid in self._sel("nodes", self.domain.nodes):
            node = self.domain.nodes.get(nid)
            if node is None:
                continue
            entry = {}
            if not n_vars or "Displacements" in n_vars:
                entry["U"] = tuple(U[node.dofs[d]] if d in node.dofs else 0.0 for d in ("ux", "uy"))
            if not n_vars or "External_Forces" in n_vars:
                entry["F"] = tuple(self.F_ext[node.dofs[d]] * factor if d in node.dofs else 0.0
                                   for d in ("ux", "uy"))
            nodes[nid] = entry
        elems = {}
        for eid in self._sel("elements", self.domain.elements):
            el = self.domain.elements.get(eid)
            if el is None:
                continue
            entry = {}
            if not e_vars or "Von_Mises" in e_vars:
                entry["Von_Mises"] = von_mises(element_average_stress(el, U))
            if not e_vars or "Internal_State" in e_vars:
                svs = getattr(getattr(el, "state", None), "vars", None) or getattr(el, "state_vars", None)
                val = 0.0
                if svs and svs[0] is not None:
                    key = getattr(getattr(el, "material", None), "PRIMARY_STATE_VAR", None)
                    vals = [sv[key] if (key and key in sv and isinstance(sv[key], (int, float)))
                            else next((v for v in sv.values()
                                       if isinstance(v, (int, float)) and not isinstance(v, bool)), 0.0)
                            for sv in svs if sv is not None]
                    val = sum(vals) / len(vals) if vals else 0.0
                entry["Internal_State"] = val
            elems[eid] = entry
        return {"step": n, "factor": factor, "nodes": nodes, "elems": elems}

    def _write_text(self):
        n_vars = self.text_cfg.get("nodal_results", [])
        e_vars = self.text_cfg.get("element_results", [])
        ruta = self.base + "_text_export.txt"
        dom = self.domain
        with open(ruta, "w", encoding="utf-8") as f:
            f.write(" * * * * *  S O L I D U M   F E M   O U T P U T  * * * * *\n\n")
            f.write(f"     Problem Name            - - - - - - : {self.base_name}\n")
            f.write(f"     Number of Nodal Points  - - - - - - : {len(dom.nodes)}\n")
            f.write(f"     Number of Elements  - - - - - - - - : {len(dom.elements)}\n")
            f.write(f"     Total Degrees-of-Freedom  - - - - - : {dom.total_dofs}\n")
            has_disp = not n_vars or "Displacements" in n_vars
            has_force = not n_vars or "External_Forces" in n_vars
            has_vm = not e_vars or "Von_Mises" in e_vars
            has_is = not e_vars or "Internal_State" in e_vars
            for st in self.text_data:
                f.write("\n ---------------------------------------------------\n\n")
                f.write(f"     Step:{st['step']:8d}    Load Factor:{st['factor']:14.5E}\n")
                if has_disp:
                    f.write("\n     N o d a l   D i s p l a c e m e n t s\n\n      Node       1 Displ       2 Displ\n")
                    for nid, d in st["nodes"].items():
                        if "U" in d:
                            f.write(f"{nid:10d}{d['U'][0]:14.4E}{d['U'][1]:14.4E}\n")
                if has_force:
                    f.write("\n     N o d a l   F o r c e s\n\n      Node       1 Force       2 Force\n")
                    for nid, d in st["nodes"].items():
                        if "F" in d:
                            f.write(f"{nid:10d}{d['F'][0]:14.4E}{d['F'][1]:14.4E}\n")
                if st["elems"]:
                    f.write("\n     E l e m e n t   V a r i a b l e s\n\n")
                    f.write("      Elmt" + ("     Von_Mises" if has_vm else "") + ("     Int_State" if has_is else "") + "\n")
                    for eid, d in st["elems"].items():
                        line = f"{eid:10d}"
                        if has_vm and "Von_Mises" in d:
                            line += f"{d['Von_Mises']:14.4E}"
                        if has_is and "Internal_State" in d:
                            line += f"{d['Internal_State']:14.4E}"
                        f.write(line + "\n")
        _log.info(f"Salida de texto (estilo FEAP) escrita en {ruta}")
