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
  de apoyos o de desplazamiento impuesto; energías; esfuerzo, deformación,
  variable interna o fuerza interna de un elemento o de un conjunto de
  elementos). Se escribe fila a fila, así que un análisis interrumpido
  conserva lo calculado;
- la salida en texto estilo FEAP (``text_export``);
- las figuras de ``plots`` (ADR 0039), dibujadas al terminar con los archivos
  anteriores (:mod:`solidum.post`).

La reacción de un paso es ``F_int(U) − factor·F`` en los grados de libertad
prescritos (sin restricciones multipunto de apoyo), la del ADR 0004 §3; sólo
en los análisis estáticos.

Energías (análisis estáticos sin etapas; ADR 0038). En equilibrio, la fuerza
que el exterior (cargas y apoyos) ejerce sobre la estructura es
``F_int(U)``, así que:

- trabajo externo ``W_n = W_{n−1} + ½ (F_int,n + F_int,n−1)·(U_n − U_{n−1})``
  (regla del trapecio entre pasos convergidos; incluye el trabajo de los
  desplazamientos impuestos);
- energía elástica ``E_n = ½ U_n·F_int,n``: la que se recupera descargando
  hasta el origen, válida si todos los materiales descargan por la secante
  (elásticos y de daño: ``SECANT_UNLOADING``) y no hay esfuerzos iniciales;
- energía disipada ``D_n = W_n − E_n``.

Las magnitudes de elemento se leen del estado confirmado del paso (el mismo
que el análisis usó), sin volver a evaluar el material; la deformación, de
``B·u``.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from solidum.logging import get_logger

_log = get_logger("output")

_FREQ_WORDS = ("last_step", "all_steps")
_REDUCE = {"mean": np.mean, "min": np.min, "max": np.max, "sum": np.sum}
_COLUMN_KINDS = ("displacement", "reaction", "energy", "element")
_SELECTORS = ("node", "nodes", "group", "at")
#: Energías de la historia.
ENERGIES = ("external_work", "elastic", "dissipated")
_ELEMENT_SELECTORS = ("id", "ids", "group", "where")
_ELEMENT_KEYS = {"var", "component", "end", "reduce", *_ELEMENT_SELECTORS}
#: Fuerzas internas de los elementos estructurales 1D (ElementForces).
FORCE_COMPONENTS = ("N", "V", "M", "Vy", "Vz", "T", "My", "Mz")
_VOIGT = {1: ("xx",), 3: ("xx", "yy", "xy"), 6: ("xx", "yy", "zz", "xy", "yz", "xz")}
_DERIVED = {"stress": ("mises", "max_principal"), "strain": ("max_principal",)}


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
    if "plots" in output:
        from solidum.post.plots import validate_plots_block
        errors.extend(validate_plots_block(output["plots"]))
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
        if kinds[0] == "energy":
            if col["energy"] not in ENERGIES:
                errors.append(f"{ctx}.energy: {col['energy']!r} no admitida; una de {list(ENERGIES)}.")
            continue
        if kinds[0] == "element":
            errors.extend(_validate_element_column(col["element"], f"{ctx}.element"))
            continue
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


def _validate_element_column(sel, ctx: str) -> list[str]:
    errors = []
    if not isinstance(sel, dict):
        return [f"{ctx}: se esperaba un diccionario con 'var' y, opcionalmente, uno de "
                f"{list(_ELEMENT_SELECTORS)}."]
    if set(sel) - _ELEMENT_KEYS:
        errors.append(f"{ctx}: claves no admitidas {sorted(set(sel) - _ELEMENT_KEYS)}; se admiten "
                      f"{sorted(_ELEMENT_KEYS)}.")
    found = [k for k in _ELEMENT_SELECTORS if k in sel]
    if len(found) > 1:
        errors.append(f"{ctx}: declare como mucho uno de {list(_ELEMENT_SELECTORS)} (sin ninguno, "
                      "todos los elementos).")
    var = sel.get("var")
    if not isinstance(var, str) or not var:
        errors.append(f"{ctx}.var: falta la magnitud ('stress', 'strain', una variable interna "
                      f"como 'damage', o una fuerza interna {list(FORCE_COMPONENTS)}).")
        return errors
    comp = sel.get("component")
    if var in ("stress", "strain"):
        valid = sorted({c for cs in _VOIGT.values() for c in cs} | {"zz"} | set(_DERIVED[var]))
        if comp not in valid:
            errors.append(f"{ctx}.component: con var: {var} se requiere uno de {valid}.")
    elif comp is not None:
        errors.append(f"{ctx}.component: sólo con var: stress o strain.")
    if var in FORCE_COMPONENTS:
        if sel.get("end") not in ("i", "j"):
            errors.append(f"{ctx}.end: con una fuerza interna se requiere el extremo 'i' o 'j'.")
    elif "end" in sel:
        errors.append(f"{ctx}.end: sólo con una fuerza interna {list(FORCE_COMPONENTS)}.")
    if "reduce" in sel and sel["reduce"] not in ("mean", "min", "max"):
        errors.append(f"{ctx}.reduce: se admite 'mean', 'min' o 'max'.")
    if "where" in sel:
        from solidum.utils.expressions import POSITION, ExpressionError, compile_predicate
        try:
            compile_predicate(sel["where"], POSITION, {"near": lambda a, b, tol=None: True})
        except ExpressionError as exc:
            errors.append(f"{ctx}.where: {exc}")
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
        from solidum.core.load_vector import equivalent_loads_of
        self.equivalent_loads = equivalent_loads_of(F_ext)
        self.F_ext = np.zeros(domain.total_dofs) if F_ext is None else np.asarray(F_ext, dtype=float)
        self.static = bool(static)
        self._W = 0.0
        self._U_prev = None
        self._F_prev = None
        self._read_trial = False
        res = cfg.get("results") or {}
        self.freq = res.get("frequency", "last_step")
        self.n_vars = res.get("nodal_results")
        self.e_vars = res.get("element_results")
        self.pre = bool((cfg.get("pre_process") or {}).get("export", False))
        self.text_cfg = cfg.get("text_export") or {}
        self.export_text = bool(self.text_cfg.get("export", False))
        self.plots = cfg.get("plots") or []
        self.history_file = (cfg.get("history") or {}).get("file")
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
            kind = next(k for k in _COLUMN_KINDS if k in col)
            ctx = f"output.history.columns[{i}] ({col['name']})"
            if kind == "energy":
                cols.append(self._energy_column(col, ctx, static))
                continue
            if kind == "element":
                cols.append({"name": col["name"], "kind": "element",
                             "probe": _ElementProbe(self, col["element"], ctx, static),
                             "reduce": _REDUCE[col["element"].get("reduce", "mean")],
                             "scale": float(col.get("scale", 1.0))})
                continue
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

    def _energy_column(self, col, ctx, static):
        which = col["energy"]
        if not static:
            raise ValueError(f"{ctx}: las energías sólo se registran en análisis estáticos sin etapas.")
        if self._F_prev is None:
            K0, F0 = self.assembler.assemble_non_linear_system(np.zeros(self.domain.total_dofs))
            self._F_prev = np.array(F0, dtype=float)
            self._U_prev = np.zeros(self.domain.total_dofs)
            self._K0_scale = float(np.max(np.abs(K0.diagonal()))) if K0.shape[0] else 0.0
        if which != "external_work":
            not_secant = sorted({name for el in self.domain.elements.values()
                                 for name in _non_secant(el)})
            if not_secant:
                raise ValueError(
                    f"{ctx}: la energía {which!r} usa E = ½·U·F_int(U), que sólo es la energía "
                    f"elástica recuperable si todos los materiales descargan por la secante hasta "
                    f"el origen (elásticos y de daño); no lo hacen {not_secant[:6]}. El trabajo "
                    f"externo ('external_work') sí está disponible.")
            coords = np.array([list(n.coordinates) + [0.0] * (3 - len(n.coordinates))
                               for n in self.domain.nodes.values()])
            size = float(np.linalg.norm(coords.max(axis=0) - coords.min(axis=0))) if len(coords) else 0.0
            if np.max(np.abs(self._F_prev), initial=0.0) > 1e-10 * self._K0_scale * size:
                raise ValueError(
                    f"{ctx}: el modelo tiene fuerzas internas en U = 0 (esfuerzos iniciales, p. ej. "
                    f"cables pretensados): la energía {which!r} no se puede obtener de ½·U·F_int(U).")
        return {"name": col["name"], "kind": "energy", "energy": which,
                "scale": float(col.get("scale", 1.0))}

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
            if self._read_trial or any(c["kind"] in ("reaction", "energy") for c in self.columns):
                _, F_int = self.assembler.assemble_non_linear_system(U)
            energies = self._energies(U, F_int) if self._F_prev is not None else None
            for c in self.columns:
                if c["kind"] == "displacement":
                    val = c["reduce"](U[c["dofs"]])
                elif c["kind"] == "reaction":
                    val = np.sum(F_int[c["dofs"]] - factor * self.F_ext[c["dofs"]])
                elif c["kind"] == "energy":
                    val = energies[c["energy"]]
                else:
                    val = c["reduce"](c["probe"].values(U, factor, self._read_trial))
                row.append(f"{c['scale'] * float(val):.10g}")
            with open(self.history_path, "a", newline="", encoding="utf-8") as f:
                csv.writer(f).writerow(row)
        if self.export_text:
            self.text_data.append(self._text_step(n, U, factor))

    def _energies(self, U, F_int) -> dict:
        U = np.asarray(U, dtype=float)
        self._W += 0.5 * float((F_int + self._F_prev) @ (U - self._U_prev))
        self._U_prev, self._F_prev = U.copy(), np.array(F_int, dtype=float)
        E = 0.5 * float(U @ F_int)
        return {"external_work": self._W, "elastic": E, "dissipated": self._W - E}

    def finish(self, U_final=None):
        """Cierra la salida: resultado final, archivo PVD y texto."""
        from solidum.utils.vtk_exporter import VtkExporter
        if self.n_steps == 0 and U_final is not None:
            # Análisis sin pasos (lineal): un único estado con factor 1. El
            # solver lineal no confirma el estado: las magnitudes de elemento
            # se leen del de prueba, que deja el ensamblaje en U.
            self._read_trial = True
            self.step(1, U_final, 1.0)
        if self.last is None:
            return
        _, U, factor = self.last
        if self.pvd:
            # También con un solo instante: el .pvd es el índice factor → archivo
            # que usan las figuras (ADR 0039).
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
        if self.plots:
            from solidum.post.plots import ResultFiles, make_plots
            try:
                make_plots(ResultFiles(self.base_dir, self.base_name, self.history_file), self.plots)
            except (FileNotFoundError, ValueError, ImportError) as exc:
                # Los resultados ya están escritos: un fallo de una figura no
                # invalida el análisis; se corrige el YAML y se regeneran con
                # python -m solidum.post.
                _log.error(f"No se pudieron dibujar las figuras: {exc} Los resultados están "
                           f"escritos; corrija output.plots y regenérelas con "
                           f"'python -m solidum.post'.")

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


def _non_secant(el) -> list[str]:
    """Nombres de los materiales del elemento que no descargan por la secante
    (``SECANT_UNLOADING``); el propio elemento si no se le encuentra material."""
    refs = [getattr(el, kw, None) for kw in getattr(el, "REFERENCE_KWARGS", {"material": None})]
    refs = [r for r in refs if r is not None]
    if not refs:
        return [type(el).__name__]
    return [type(r).__name__ for r in refs if not getattr(r, "SECANT_UNLOADING", False)]


class _ElementProbe:
    """Valores de una magnitud de elemento en los elementos de un selector:
    por punto de integración (esfuerzo, deformación, variable interna) o por
    elemento (fuerza interna en un extremo)."""

    def __init__(self, writer, sel: dict, ctx: str, static: bool):
        self.writer = writer
        self.var = sel["var"]
        self.component = sel.get("component")
        self.end = 0 if sel.get("end", "i") == "i" else 1
        self.ctx = ctx
        elements = self._select(sel)
        explicit = any(k in sel for k in _ELEMENT_SELECTORS)
        ok = [el for el in elements if self._supports(el)]
        if explicit and len(ok) != len(elements):
            bad = [el.id for el in elements if el not in ok]
            raise ValueError(f"{ctx}: los elementos {bad[:10]} no tienen {self._what()}.")
        if not ok:
            raise ValueError(f"{ctx}: ningún elemento del selector tiene {self._what()}.")
        self.elements = ok
        if self.var in FORCE_COMPONENTS:
            loads = writer.equivalent_loads
            self.f_eq = {el.id: (None if loads is None else loads.get(el.id)) for el in ok}
            if not static and any(v is not None for v in self.f_eq.values()):
                raise ValueError(f"{ctx}: las fuerzas internas con carga distribuida sólo se "
                                 "registran en análisis estáticos.")
            return
        if self.var in ("stress", "strain"):
            dims = {el.STRAIN_DIM for el in ok}
            if len(dims) != 1:
                raise ValueError(f"{ctx}: elementos con distinto número de componentes {sorted(dims)}.")
            self.n_sigma = dims.pop()
            comps = _VOIGT[self.n_sigma]
            allowed = set(comps) | set(_DERIVED[self.var]) | (
                {"zz"} if self.var == "stress" and self.n_sigma == 3 else set())
            if self.component not in allowed:
                raise ValueError(f"{ctx}: component {self.component!r} no existe con "
                                 f"{self.n_sigma} componentes; se admiten {sorted(allowed)}.")

    # -- selección ---------------------------------------------------------
    def _select(self, sel):
        dom = self.writer.domain
        if "id" in sel:
            ids = [sel["id"]]
        elif "ids" in sel:
            ids = list(sel["ids"])
        elif "group" in sel:
            groups = getattr(dom, "physical_groups", None) or {}
            if sel["group"] not in groups:
                raise ValueError(f"{self.ctx}: grupo físico {sel['group']!r} inexistente. "
                                 f"Disponibles: {sorted(groups)}.")
            inside = set(groups[sel["group"]])
            ids = [eid for eid, el in dom.elements.items()
                   if all(nd.id in inside for nd in el.nodes)]
            if not ids:
                raise ValueError(f"{self.ctx}: ningún elemento tiene todos sus nodos en el grupo "
                                 f"{sel['group']!r}.")
        elif "where" in sel:
            from solidum.utils.expressions import POSITION, compile_predicate
            pts = np.array([_pad3(nd.coordinates) for nd in dom.nodes.values()])
            scale = float(np.linalg.norm(pts.max(axis=0) - pts.min(axis=0)))

            def near(a, b, tol=None):
                return abs(a - b) <= (1e-6 * scale if tol is None else tol)
            pred = compile_predicate(sel["where"], POSITION, {"near": near})
            ids = [eid for eid, el in dom.elements.items()
                   if pred(*np.mean([_pad3(nd.coordinates) for nd in el.nodes], axis=0))]
            if not ids:
                raise ValueError(f"{self.ctx}: ningún elemento tiene el centro donde {sel['where']!r}.")
        else:
            ids = list(dom.elements)
        missing = [i for i in ids if i not in dom.elements]
        if missing:
            raise ValueError(f"{self.ctx}: elementos inexistentes {missing[:10]}.")
        return [dom.elements[i] for i in ids]

    def _what(self) -> str:
        if self.var in FORCE_COMPONENTS:
            return f"la fuerza interna {self.var!r}"
        if self.var in ("stress", "strain"):
            return f"{'esfuerzos' if self.var == 'stress' else 'deformaciones'} por punto de integración"
        return f"la variable interna escalar {self.var!r}"

    def _supports(self, el) -> bool:
        if self.var in FORCE_COMPONENTS:
            try:
                ef = el.internal_forces(np.zeros(self.writer.domain.total_dofs))
            except Exception:  # noqa: BLE001 — un elemento sin el contrato no lo admite
                return False
            return ef is not None and self.var in ef.components
        # Sólo elementos mecánicos con esfuerzos por punto de integración
        # (sólidos y barras): no los térmicos (flujo) ni los marcos
        # (resultantes de sección).
        mechanical = (isinstance(getattr(el, "STRAIN_DIM", None), int)
                      and not {"rx", "ry", "rz"} & set(getattr(el, "DOF_NAMES", ())))
        state = getattr(el, "state", None)
        if not mechanical or state is None or not len(getattr(state, "stresses", None) or ()):
            return False
        if self.var == "stress":
            return True
        if self.var == "strain":
            return hasattr(el, "compute_gauss_state")
        schema = getattr(getattr(el, "material", None), "STATE_SCHEMA", None) or {}
        return tuple(schema.get(self.var, (None,))) == ()

    # -- valores -------------------------------------------------------------
    def values(self, U, factor, trial: bool) -> np.ndarray:
        if self.var in FORCE_COMPONENTS:
            out = []
            for el in self.elements:
                f_eq = self.f_eq[el.id]
                ef = el.internal_forces(U, equivalent_load=None if f_eq is None else factor * f_eq)
                out.append(float(ef.components[self.var][self.end]))
            return np.array(out)
        if self.var == "strain":
            eps = np.vstack([np.atleast_2d(el.compute_gauss_state(U)["strain"])
                             for el in self.elements])
            return _component(eps, self.component, strain=True)
        # Elementos de una familia del ensamblaje por lotes: lectura directa de
        # sus arreglos de estado; el resto, elemento a elemento.
        parts_sig, parts_out, loose = [], [], []
        for fam, rows in self._family_rows():
            if fam is None:
                loose = rows
                continue
            st = fam.state
            S = st.S_trial if trial else st.S_committed
            sig = (st.sig_trial if trial else st.sig_committed)[rows]
            if self.var == "stress":
                parts_sig.append(sig)
                need_szz = self.n_sigma == 3 and self.component in ("zz", "mises", "max_principal")
                parts_out.append(np.asarray(fam.material.batch_out_of_plane_stress(sig, S[rows]))
                                 if need_szz else np.zeros(len(rows)))
            else:
                k = fam.schema.names.index(self.var)
                parts_out.append(S[rows, fam.schema.offsets[k]])
        loose_sig, loose_out = [], []
        for el in loose:
            st = el.state
            stresses = st.stresses_trial if trial else st.stresses
            vars_ = st.vars_trial if trial else st.vars
            for g in range(len(stresses)):
                if self.var == "stress":
                    sig = np.atleast_1d(np.asarray(stresses[g], dtype=float))
                    loose_sig.append(sig)
                    szz = 0.0
                    if self.n_sigma == 3 and self.component in ("zz", "mises", "max_principal"):
                        szz = float(el.material.out_of_plane_stress(sig, vars_[g]))
                    loose_out.append(szz)
                else:
                    loose_out.append(float(vars_[g][self.var]))
        out = np.concatenate([*(np.asarray(p, dtype=float).reshape(-1) for p in parts_out),
                              np.asarray(loose_out, dtype=float)])
        if self.var != "stress":
            return out
        sig = np.vstack([*(np.atleast_2d(p) for p in parts_sig),
                         *(np.atleast_2d(r) for r in loose_sig)])
        return _component(sig, self.component, szz=out)

    def _family_rows(self):
        """``[(familia, filas)]`` de los elementos en familia y ``(None,
        [elementos sueltos])``; se recalcula si alguna familia dejó de estar
        adoptada."""
        from solidum.math.batch.postprocess import family_of
        plan = getattr(self, "_plan", None)
        if plan is not None and all(f is None or f.adopted for f, _ in plan):
            return plan
        groups, loose = {}, []
        for el in self.elements:
            fam = family_of(el)
            if fam is None:
                loose.append(el)
            else:
                groups.setdefault(fam, []).append(el)
        plan = []
        for fam, members in groups.items():
            index = {id(e): i for i, e in enumerate(fam.elements)}
            rows = np.concatenate([np.arange(fam.state.rows(index[id(e)]).start,
                                             fam.state.rows(index[id(e)]).stop)
                                   for e in members])
            plan.append((fam, rows))
        plan.append((None, loose))
        self._plan = plan
        return plan


def _pad3(c) -> list:
    c = [float(v) for v in c]
    return c + [0.0] * (3 - len(c))


def _component(v: np.ndarray, component: str, *, strain: bool = False,
               szz: np.ndarray | None = None) -> np.ndarray:
    """Componente o invariante de ``v`` ``(P, n)`` en Voigt del proyecto."""
    n = v.shape[1]
    if component in _VOIGT[n]:
        return v[:, _VOIGT[n].index(component)]
    if component == "zz":
        return szz
    # Tensor completo 3×3 por punto (en 2D, σzz fuera del plano o εzz = 0).
    T = np.zeros((v.shape[0], 3, 3))
    shear = 0.5 if strain else 1.0              # γ = 2ε en la deformación
    if n == 1:
        T[:, 0, 0] = v[:, 0]
    elif n == 3:
        T[:, 0, 0], T[:, 1, 1] = v[:, 0], v[:, 1]
        T[:, 2, 2] = 0.0 if szz is None else szz
        T[:, 0, 1] = T[:, 1, 0] = shear * v[:, 2]
    else:
        T[:, 0, 0], T[:, 1, 1], T[:, 2, 2] = v[:, 0], v[:, 1], v[:, 2]
        T[:, 0, 1] = T[:, 1, 0] = shear * v[:, 3]
        T[:, 1, 2] = T[:, 2, 1] = shear * v[:, 4]
        T[:, 0, 2] = T[:, 2, 0] = shear * v[:, 5]
    if component == "max_principal":
        if strain and n == 3:
            # En 2D la deformación fuera del plano no se conoce: principal en el plano.
            return np.linalg.eigvalsh(T[:, :2, :2])[:, -1]
        return np.linalg.eigvalsh(T)[:, -1]
    dev = T - np.trace(T, axis1=1, axis2=2)[:, None, None] * np.eye(3) / 3.0
    return np.sqrt(1.5 * np.einsum("pij,pij->p", dev, dev))
