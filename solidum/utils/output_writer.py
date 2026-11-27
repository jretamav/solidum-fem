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
  elementos; salto o tracción locales de una interfaz, ``var: jump |
  traction`` con ``component: n | s | s1 | s2`` (``n1``, ``n2`` en la
  adherencia 3D), y los elementos de una interfaz por su nombre,
  ``interface: <nombre>``, ADR 0043 §7). Se escribe fila a fila, así que un
  análisis interrumpido conserva lo calculado;
- la salida en texto estilo FEAP (``text_export``);
- las figuras de ``plots`` (ADR 0039), dibujadas al terminar con los archivos
  anteriores (:mod:`solidum.post`).

La reacción de un paso es ``F_int(U) − F`` en los grados de libertad
prescritos (sin restricciones multipunto de apoyo), la del ADR 0004 §3, con
``F`` la carga aplicada del paso: ``factor·F_ref`` o, en una etapa,
``F_c + λ·F_ref`` (ADR 0026); sólo en los análisis estáticos.

Con etapas (``stages``) la historia lleva además la columna ``stage``, el
número de la etapa (desde 1); ``factor`` es el λ de la etapa, que vuelve a
empezar en 0 en cada una. Las reacciones y las energías se registran si
ninguna etapa es en el tiempo; la reacción de un grado de libertad que en esa
etapa todavía está libre es una celda vacía («no aplica»).

Energías (análisis estáticos, también por etapas; ADR 0038). En equilibrio, la fuerza
que el exterior (cargas y apoyos) ejerce sobre la estructura es
``F_int(U)``, así que:

- trabajo externo ``W_n = W_{n−1} + ½ (F_int,n + F_int,n−1)·(U_n − U_{n−1})``
  (regla del trapecio entre pasos convergidos; incluye el trabajo de los
  desplazamientos impuestos);
- energía elástica (almacenada) ``E_n = Σ_e E_e``: ``½ u_e·f_e`` en cada
  elemento, la que se recupera descargando hasta el origen, válida si sus
  materiales descargan por la secante (elásticos y de daño:
  ``SECANT_UNLOADING``) y no hay esfuerzos iniciales; los elementos que
  conocen su energía almacenada la entregan (``Element.stored_energy``: las
  interfaces, ``Σ_g m_g·ψ_g`` de su ley, que cubre la fricción y la holgura
  inicial; ADR 0043 §5). Sin interfaces, ``E_n = ½ U_n·F_int,n``;
- energía disipada ``D_n = W_n − E_n``. Con la regla del trapecio, un paso
  que cruza un cambio brusco de rigidez (una interfaz que se cierra o empieza
  a deslizar dentro del paso) añade a ``W`` un error que aparece como
  disipación y baja con el tamaño del paso.

Las magnitudes de elemento se leen del estado confirmado del paso (el mismo
que el análisis usó), sin volver a evaluar el material; la deformación, de
``B·u``.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from solidum.core.interface_material import all_component_names
from solidum.core.point_quantity import JUMP, STRAIN, point_quantity
from solidum.logging import get_logger

_log = get_logger("output")

_FREQ_WORDS = ("last_step", "all_steps")
_REDUCE = {"mean": np.mean, "min": np.min, "max": np.max, "sum": np.sum}
_COLUMN_KINDS = ("displacement", "reaction", "energy", "element")
_SELECTORS = ("node", "nodes", "group", "at")
#: Energías de la historia.
ENERGIES = ("external_work", "elastic", "dissipated")
_ELEMENT_SELECTORS = ("id", "ids", "group", "where", "interface")
_ELEMENT_KEYS = {"var", "component", "end", "reduce", *_ELEMENT_SELECTORS}
#: Fuerzas internas de los elementos estructurales 1D (ElementForces).
FORCE_COMPONENTS = ("N", "V", "M", "Vy", "Vz", "T", "My", "Mz")
_VOIGT = {1: ("xx",), 3: ("xx", "yy", "xy"), 6: ("xx", "yy", "zz", "xy", "yz", "xz")}
_DERIVED = {"stress": ("mises", "max_principal"), "strain": ("max_principal",)}
#: Magnitudes por punto de los elementos de interfaz (ADR 0043 §7), en ejes locales.
INTERFACE_VARS = ("jump", "traction")
#: Componentes locales de una interfaz: de superficie (n, s / n, s1, s2) o de
#: la adherencia de una barra (n, s / n1, n2, s); de la fuente única de la
#: familia (``solidum.core.interface_material``).
INTERFACE_COMPONENTS = all_component_names()


def validate_output_block(output) -> list[str]:
    """Errores de forma de ``output.results.frequency`` y ``output.history``."""
    errors = []
    if not isinstance(output, dict):
        return errors
    res = output.get("results") or {}
    if isinstance(res, dict) and not isinstance(res.get("integration_points", False), bool):
        errors.append("output.results.integration_points: se esperaba true o false.")
    if isinstance(res, dict):
        for key in ("nodal_results", "element_results"):
            names = res.get(key)
            if names is not None and not (isinstance(names, list)
                                          and all(isinstance(v, str) for v in names)):
                errors.append(f"output.results.{key}: se esperaba una lista de nombres de campos.")
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
        if not isinstance(name, str) or not name or name in ("step", "stage", "factor") or name in names:
            errors.append(f"{ctx}.name: se esperaba un nombre único (no 'step', 'stage' ni 'factor').")
        names.add(name)
        kinds = [k for k in _COLUMN_KINDS if k in col]
        if len(kinds) != 1:
            errors.append(f"{ctx}: declare exactamente una magnitud: {list(_COLUMN_KINDS)}.")
            continue
        extra = set(col) - {"name", "scale", kinds[0]} - ({"elements"} if kinds[0] == "energy"
                                                          else set())
        if extra:
            errors.append(f"{ctx}: claves no admitidas {sorted(extra)}.")
        if "scale" in col and (isinstance(col["scale"], bool) or not isinstance(col["scale"], (int, float))):
            errors.append(f"{ctx}.scale: se esperaba un número.")
        if kinds[0] == "energy":
            if col["energy"] not in ENERGIES:
                errors.append(f"{ctx}.energy: {col['energy']!r} no admitida; una de {list(ENERGIES)}.")
            if "elements" in col:
                errors.extend(_validate_zone(col["elements"], f"{ctx}.elements"))
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


def _validate_zone(sel, ctx: str) -> list[str]:
    """Selector de elementos de una energía por zona: exactamente uno de
    ``id``, ``ids``, ``group`` o ``where``."""
    if not isinstance(sel, dict):
        return [f"{ctx}: se esperaba un diccionario con uno de {list(_ELEMENT_SELECTORS)}."]
    errors = []
    found = [k for k in _ELEMENT_SELECTORS if k in sel]
    if len(found) != 1 or set(sel) - set(_ELEMENT_SELECTORS):
        errors.append(f"{ctx}: declare exactamente uno de {list(_ELEMENT_SELECTORS)}.")
    if "where" in sel:
        from solidum.utils.expressions import POSITION, ExpressionError, compile_predicate
        try:
            compile_predicate(sel["where"], POSITION, {"near": lambda a, b, tol=None: True})
        except ExpressionError as exc:
            errors.append(f"{ctx}.where: {exc}")
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
        errors.append(f"{ctx}.var: falta la magnitud ('stress', 'strain', 'jump', 'traction', una "
                      f"variable interna como 'damage', o una fuerza interna "
                      f"{list(FORCE_COMPONENTS)}).")
        return errors
    comp = sel.get("component")
    if var in ("stress", "strain"):
        valid = sorted({c for cs in _VOIGT.values() for c in cs} | {"zz"} | set(_DERIVED[var]))
        if comp not in valid:
            errors.append(f"{ctx}.component: con var: {var} se requiere uno de {valid}.")
    elif var in INTERFACE_VARS:
        if comp not in INTERFACE_COMPONENTS:
            errors.append(f"{ctx}.component: con var: {var} se requiere una componente local "
                          f"{list(INTERFACE_COMPONENTS)}.")
    elif comp is not None:
        errors.append(f"{ctx}.component: sólo con var: stress, strain, jump o traction.")
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
        sólo se definen ahí); con etapas, si ninguna es en el tiempo.
    staged : bool
        ``True`` en un análisis por etapas: la historia lleva la columna
        ``stage`` (número de la etapa, desde 1) y la carga de cada etapa la da
        :meth:`begin_stage`.
    stage_prescribed : set, optional
        Pares ``(nodo, dof)`` que prescribe alguna etapa (ADR 0026 §7): su
        reacción se admite aunque estén libres en el nivel raíz.
    stage_distributed_loads : bool
        ``True`` si alguna etapa tiene carga distribuida: con una etapa en el
        tiempo, las fuerzas internas de marco de la historia serían sin su
        carga equivalente, y se rechazan.
    """

    def __init__(self, domain, assembler, output_config: dict, yaml_path, F_ext, *,
                 static: bool = True, staged: bool = False, stage_prescribed: set | None = None,
                 stage_distributed_loads: bool = False):
        self.domain, self.assembler = domain, assembler
        cfg = output_config or {}
        yaml_path = Path(yaml_path)
        self.base_dir = yaml_path.parent
        self.base_name = cfg.get("file_name") or yaml_path.stem
        self.base = str(self.base_dir / self.base_name)
        from solidum.core.load_vector import equivalent_loads_of
        self.equivalent_loads = equivalent_loads_of(F_ext)
        self.F_ext = np.zeros(domain.total_dofs) if F_ext is None else np.asarray(F_ext, dtype=float)
        # Carga aplicada de un paso: F_const + factor·F_ext (ADR 0026); sin
        # etapas, F_const = 0. Cargas equivalentes por elemento de cada parte
        # (ADR 0027), para las fuerzas internas de la historia.
        self.F_const = np.zeros(domain.total_dofs)
        self._eq_const = None
        self._eq_ref = self.equivalent_loads
        self.static = bool(static)
        self.staged = bool(staged)
        self.stage = None
        self._stage_prescribed = set(stage_prescribed or ())
        self._stage_distributed_loads = bool(stage_distributed_loads)
        self._W = 0.0
        self._U_prev = None
        self._F_prev = None
        # Elementos que entregan su energía almacenada; None mientras ninguna
        # columna pida 'elastic' o 'dissipated' (ADR 0043 §5).
        self._stored_elements = None
        self._read_trial = False
        res = cfg.get("results") or {}
        self.freq = res.get("frequency", "last_step")
        self.n_vars = res.get("nodal_results")
        self.e_vars = res.get("element_results")
        self._check_field_names()
        self.integration_points = bool(res.get("integration_points", False))
        self.pvd_points: list = []
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
                csv.writer(f).writerow(["step", *(["stage"] if self.staged else []), "factor",
                                        *[c["name"] for c in self.columns]])

    def begin_stage(self, number: int, F_constant=None, F_reference=None) -> None:
        """Empieza la etapa ``number`` (desde 1) de un análisis por etapas.

        En una etapa estática la carga aplicada de un paso es ``F_constant +
        λ·F_reference`` (ADR 0026): la de las etapas anteriores, constante, y
        la de ésta, que crece con λ; con ella se calculan las reacciones y las
        fuerzas internas con carga distribuida. En una etapa en el tiempo,
        ``None, None``. El trabajo externo y las energías siguen acumulándose
        de una etapa a otra: la regla del trapecio usa ``F_int``, que no
        depende de cómo se reparta la carga."""
        from solidum.core.load_vector import equivalent_loads_of
        n = self.domain.total_dofs
        self.stage = int(number)
        self.F_const = np.zeros(n) if F_constant is None else np.asarray(F_constant, dtype=float)
        self.F_ext = np.zeros(n) if F_reference is None else np.asarray(F_reference, dtype=float)
        self._eq_const = equivalent_loads_of(F_constant)
        self._eq_ref = equivalent_loads_of(F_reference)

    def _applied(self, factor: float) -> np.ndarray:
        """Carga aplicada del paso: ``F_const + factor·F_ext``."""
        return self.F_const + factor * self.F_ext

    def element_equivalent_load(self, element_id: int, factor: float):
        """Carga equivalente de la carga distribuida del elemento en el paso
        (ADR 0027): la constante más ``factor`` veces la de referencia;
        ``None`` si el elemento no tiene."""
        parts = []
        if self._eq_const is not None and (c := self._eq_const.get(element_id)) is not None:
            parts.append(np.asarray(c, dtype=float))
        if self._eq_ref is not None and (r := self._eq_ref.get(element_id)) is not None:
            parts.append(factor * np.asarray(r, dtype=float))
        return sum(parts) if parts else None

    # ------------------------------------------------------------------
    def _check_field_names(self) -> None:
        """Un nombre de ``nodal_results`` o ``element_results`` que no es un
        campo del VTK de este modelo es un error con la lista de los que hay
        (antes se ignoraba en silencio)."""
        from solidum.utils.vtk_exporter import (NODAL_FIELD_NAMES, STRESS_COMPONENT_NAMES,
                                                state_variable_names)
        if self.n_vars is not None:
            unknown = [v for v in self.n_vars if v not in NODAL_FIELD_NAMES]
            if unknown:
                raise ValueError(f"output.results.nodal_results: campos inexistentes {unknown}; "
                                 f"se admiten {list(NODAL_FIELD_NAMES)}.")
        if self.e_vars is not None:
            elements = list(self.domain.elements.values())
            states, angles = state_variable_names(elements)
            # Salto y tracción de las interfaces con celda (ADR 0043 §7).
            local = dict.fromkeys(
                f"{p}_{c}" for el in elements
                if point_quantity(el).kind == JUMP and getattr(el, "VTK_EXPORT", True) is not False
                for p in ("Jump", "Traction") for c in el.jump_components())
            known = ["Internal_State", *STRESS_COMPONENT_NAMES, "Von_Mises", "Flux",
                     *(f"{c}_{e}" for c in FORCE_COMPONENTS for e in ("i", "j")),
                     *local, *(s for s in states if s not in angles)]
            unknown = [v for v in self.e_vars if v not in known]
            if unknown:
                hint = (f" {sorted(angles & set(unknown))} son ángulos o códigos: no se promedian "
                        f"por celda; se leen en el archivo por puntos de integración "
                        f"(results.integration_points: true)." if angles & set(unknown) else "")
                raise ValueError(f"output.results.element_results: campos inexistentes en este "
                                 f"modelo {unknown}; se admiten {known}.{hint}")

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
        if "node" in sel or "nodes" in sel:
            # Un nodo duplicado por ``split`` (ADR 0043 §6): el id no dice de qué lado.
            for i in ids:
                copies = dom.split_copies(i)
                if copies is not None:
                    raise ValueError(
                        f"{ctx}: el nodo {i} se duplicó al separar una interfaz (nodos "
                        f"{list(copies)}, uno por lado); el id no dice de qué lado se lee. Use "
                        f"'group' (los grupos <nombre>_a y <nombre>_b de la interfaz) o 'at'.")
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
                raise ValueError(f"{ctx}: la reacción sólo se registra en análisis estáticos "
                                 f"(con etapas, si ninguna es en el tiempo).")
            sel = col[kind]
            nodes = self._nodes_of(sel, ctx)
            dof = sel["dof"]
            lacking = [n.id for n in nodes if dof not in n.dofs]
            if lacking:
                raise ValueError(f"{ctx}: los nodos {lacking[:10]} no tienen el grado de libertad {dof!r}.")
            dofs = np.array([n.dofs[dof] for n in nodes], dtype=int)
            if kind == "reaction":
                # Con etapas basta que alguna lo prescriba (ADR 0026 §7): donde
                # está libre, F_int − F es el residuo de equilibrio (≈ 0).
                free = [n.id for n in nodes if dof not in n.boundary_conditions
                        and (n.id, dof) not in self._stage_prescribed]
                if free:
                    where = " (en el nivel raíz o en alguna etapa)" if self.staged else ""
                    raise ValueError(f"{ctx}: la reacción necesita '{dof}' prescrito en todos los nodos "
                                     f"del selector{where}; libres: {free[:10]}.")
            cols.append({"name": col["name"], "kind": kind, "dofs": dofs, "nodes": nodes,
                         "dof": dof,
                         "reduce": _REDUCE[sel.get("reduce", "mean") if kind == "displacement" else "sum"],
                         "scale": float(col.get("scale", 1.0))})
        return cols

    @staticmethod
    def _reaction(c, F_int, F_applied):
        """Reacción sumada de la columna ``c`` en los nodos donde el grado de
        libertad está prescrito ahora; ``None`` si en ninguno (una etapa
        anterior a la que lo prescribe: no hay apoyo, y ``F_int − F`` sería
        sólo el residuo del equilibrio)."""
        on = np.array([c["dof"] in n.boundary_conditions for n in c["nodes"]], dtype=bool)
        if not on.any():
            return None
        dofs = c["dofs"][on]
        return np.sum(F_int[dofs] - F_applied[dofs])

    def _energy_column(self, col, ctx, static):
        which = col["energy"]
        if not static:
            raise ValueError(f"{ctx}: las energías sólo se registran en análisis estáticos (con "
                             f"etapas, si ninguna es en el tiempo: la energía cinética y la "
                             f"amortiguada no entran en la cuenta).")
        if "elements" in col:
            return {"name": col["name"], "kind": "zone_energy",
                    "zone": _ZoneEnergy(self, col["elements"], ctx, which),
                    "scale": float(col.get("scale", 1.0))}
        if self._F_prev is None:
            K0, F0 = self.assembler.assemble_non_linear_system(np.zeros(self.domain.total_dofs))
            self._F_prev = np.array(F0, dtype=float)
            self._U_prev = np.zeros(self.domain.total_dofs)
            self._K0_scale = float(np.max(np.abs(K0.diagonal()))) if K0.shape[0] else 0.0
        if which != "external_work":
            known, unknown = _stored_energy_split(self.domain.elements.values(),
                                                  self.domain.total_dofs)
            if unknown:
                raise ValueError(
                    f"{ctx}: la energía {which!r} necesita la energía almacenada de cada elemento. "
                    f"Los de interfaz la toman de su ley; los demás usan ½·u_e·f_e, que sólo es la "
                    f"energía elástica recuperable si sus materiales descargan por la secante "
                    f"hasta el origen (elásticos y de daño). No se conoce para {unknown[:6]}. El "
                    f"trabajo externo ('external_work') sí está disponible.")
            self._stored_elements = known
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
            VtkExporter(self.domain).export(str(self.base_dir / name), U, self._applied(factor),
                                            self.n_vars, self.e_vars)
            self.pvd.append((factor, name))
            if self.integration_points:
                name = f"{self.base_name}_puntos_{n:04d}.vtu"
                if VtkExporter(self.domain).export_integration_points(
                        str(self.base_dir / name), U, trial=self._read_trial):
                    self.pvd_points.append((factor, name))
        if self.columns:
            row = [n, *([self.stage] if self.staged else []), factor]
            F_int = None
            if self._read_trial:
                # Análisis sin estado confirmado (lineal): las magnitudes de
                # elemento se leen del estado de prueba, que deja este
                # ensamblaje en U.
                _, F_int = self.assembler.assemble_non_linear_system(U)
            elif any(c["kind"] in ("reaction", "energy") for c in self.columns):
                # La F_int del paso confirmado, sin reensamblar (ADR 0040).
                F_int = self.assembler.internal_forces_at(U)
            energies = self._energies(U, F_int) if self._F_prev is not None else None
            F_applied = self._applied(factor)
            for c in self.columns:
                if c["kind"] == "displacement":
                    val = c["reduce"](U[c["dofs"]])
                elif c["kind"] == "reaction":
                    val = self._reaction(c, F_int, F_applied)
                elif c["kind"] == "energy":
                    val = energies[c["energy"]]
                elif c["kind"] == "zone_energy":
                    val = c["zone"].value(U, self._read_trial)
                else:
                    val = c["reduce"](c["probe"].values(U, factor, self._read_trial))
                # Celda vacía: «no aplica» (las figuras la leen como NaN, un hueco).
                row.append("" if val is None else f"{c['scale'] * float(val):.10g}")
            with open(self.history_path, "a", newline="", encoding="utf-8") as f:
                csv.writer(f).writerow(row)
        if self.export_text:
            self.text_data.append(self._text_step(n, U, factor))

    def _energies(self, U, F_int) -> dict:
        U = np.asarray(U, dtype=float)
        self._W += 0.5 * float((F_int + self._F_prev) @ (U - self._U_prev))
        self._U_prev, self._F_prev = U.copy(), np.array(F_int, dtype=float)
        E = None
        if self._stored_elements is not None:
            # E = Σ_e E_e: ½·u_e·f_e salvo en los elementos que entregan su
            # energía almacenada (interfaces; ADR 0043 §5).
            E = 0.5 * float(U @ F_int)
            if self._stored_elements:
                E += _stored_energy_correction(self._stored_elements, U, self._read_trial)
        return {"external_work": self._W, "elastic": E,
                "dissipated": None if E is None else self._W - E}

    def finish_stopped(self, exc: BaseException) -> None:
        """Cierra la salida de un análisis detenido (``exc``) con el último
        paso convergido: VTK final, ``.pvd``, texto y figuras. La historia ya
        está escrita fila a fila. Si escribe, lo anota en ``exc`` (nota de
        PEP 678, que muestran la traza y ``python -m solidum``). No lanza: un
        fallo aquí se registra y no tapa el error que detuvo el análisis."""
        if self.last is None:
            _log.warning("El análisis se detuvo antes del primer paso convergido: no hay "
                         "resultados que escribir (salvo el pre-proceso).")
            return
        n, _, factor = self.last
        _log.warning(f"El análisis se detuvo ({type(exc).__name__}); se escriben los resultados "
                     f"del último paso convergido: paso {n}, factor {factor:.6g}.")
        try:
            self.finish()
        except Exception as err:  # noqa: BLE001 — no tapar el error del análisis
            _log.error(f"No se pudieron cerrar los resultados del análisis detenido: {err}")
            return
        note = (f"Resultados escritos hasta el último paso convergido (paso {n}, factor "
                f"{factor:.6g}): historia, VTK final y figuras en {self.base_dir}.")
        if hasattr(exc, "add_note"):
            exc.add_note(note)
        else:  # Python 3.10: PEP 678 llega en 3.11; la lista la lee python -m solidum
            exc.__notes__ = [*getattr(exc, "__notes__", ()), note]

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
            _write_pvd(self.base + ".pvd", self.pvd)
        if self.pvd_points:
            _write_pvd(self.base + "_puntos.pvd", self.pvd_points)
        VtkExporter(self.domain).export(self.base + "_final.vtu", U, self._applied(factor),
                                        self.n_vars, self.e_vars)
        if self.integration_points:
            VtkExporter(self.domain).export_integration_points(
                self.base + "_puntos_final.vtu", U, trial=self._read_trial)
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
                # solidum-post (= python -m solidum.post).
                _log.error(f"No se pudieron dibujar las figuras: {exc} Los resultados están "
                           f"escritos; corrija output.plots y regenérelas con "
                           f"'solidum-post modelo.yaml' (o 'python -m solidum.post modelo.yaml').")

    # ------------------------------------------------------------------
    # Salida de texto estilo FEAP (trasladada de examples/ejecutar_yaml.py)
    def _sel(self, key, pool):
        v = self.text_cfg.get(key, [])
        return list(pool) if isinstance(v, str) and v.lower() == "all" else list(v)

    def _text_step(self, n, U, factor):
        from solidum.utils.vtk_exporter import element_average_stress, von_mises
        n_vars = self.text_cfg.get("nodal_results", [])
        e_vars = self.text_cfg.get("element_results", [])
        F = self._applied(factor)
        nodes = {}
        for nid in self._sel("nodes", self.domain.nodes):
            node = self.domain.nodes.get(nid)
            if node is None:
                continue
            entry = {}
            if not n_vars or "Displacements" in n_vars:
                entry["U"] = tuple(U[node.dofs[d]] if d in node.dofs else 0.0 for d in ("ux", "uy"))
            if not n_vars or "External_Forces" in n_vars:
                entry["F"] = tuple(F[node.dofs[d]] if d in node.dofs else 0.0
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


def _write_pvd(path: str, entries: list) -> None:
    """Índice ``factor → archivo`` de ParaView (y de las figuras, ADR 0039)."""
    with open(path, "w", encoding="utf-8") as f:
        f.write('<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" '
                'byte_order="LittleEndian">\n<Collection>\n')
        for fac, name in entries:
            f.write(f'  <DataSet timestep="{fac:.6f}" group="" part="0" file="{name}"/>\n')
        f.write("</Collection>\n</VTKFile>\n")


def _non_secant(el) -> list[str]:
    """Nombres de los materiales del elemento que no descargan por la secante
    (``SECANT_UNLOADING``); el propio elemento si no se le encuentra material."""
    refs = [getattr(el, kw, None) for kw in getattr(el, "REFERENCE_KWARGS", {"material": None})]
    refs = [r for r in refs if r is not None]
    if not refs:
        return [type(el).__name__]
    return [type(r).__name__ for r in refs if not getattr(r, "SECANT_UNLOADING", False)]


def _stored_energy_split(elements, n_dofs: int) -> tuple[list, list[str]]:
    """Elementos que entregan su energía almacenada (``Element.stored_energy``,
    ADR 0043 §5) y nombres de los materiales cuya energía no se conoce: los de
    esos elementos que no la definen y, en los demás, los que no descargan por
    la secante (su término es ``½·u_e·f_e``, ADR 0038)."""
    U0 = np.zeros(n_dofs)
    known, unknown = [], set()
    for el in elements:
        try:
            energy = el.stored_energy(U0)
        except NotImplementedError:
            unknown.update(type(m).__name__ for m in el.materials())
            continue
        if energy is None:
            unknown.update(_non_secant(el))
            continue
        # Contrato de la pareja (ADR 0043 §5): quien entrega su energía da
        # también el f_e del mismo estado; y G_e (force_matrix) es fija, de
        # desplazamientos pequeños. Se comprueba aquí, al crear la columna, y
        # no en el primer paso.
        if not callable(getattr(el, "nodal_forces", None)):
            raise TypeError(f"{type(el).__name__}(id={el.id}) define stored_energy sin "
                            f"nodal_forces (ADR 0043 §5).")
        if getattr(el, "GEOMETRICALLY_NONLINEAR", False) and _is_interface(el):
            raise ValueError(
                f"{type(el).__name__}(id={el.id}): las energías de la historia usan la matriz "
                f"G_e de desplazamientos pequeños de las interfaces; con desplazamientos "
                f"grandes no están disponibles.")
        known.append(el)
    return known, sorted(unknown)


def _stored_energy_correction(elements, U: np.ndarray, trial: bool) -> float:
    """``Σ_e (E_e − ½·u_e·f_e)`` de los elementos que entregan su energía: lo
    que hay que sumar a ``½·U·F_int`` para tener ``Σ_e E_e``."""
    total = 0.0
    for el in elements:
        u_e = U[el.get_global_dof_indices()]
        total += el.stored_energy(U, trial) - 0.5 * float(u_e @ el.nodal_forces(trial))
    return total


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
            if not static and (any(v is not None for v in self.f_eq.values())
                               or writer._stage_distributed_loads):
                # Con etapas, las cargas de cada una no se conocen aquí: basta
                # que alguna tenga carga distribuida (la constante de una etapa
                # estática sigue en la etapa en el tiempo).
                raise ValueError(f"{ctx}: las fuerzas internas con carga distribuida sólo se "
                                 "registran en análisis estáticos (con etapas, si ninguna es en "
                                 "el tiempo).")
            return
        if self.var in ("stress", "strain"):
            dims = {point_quantity(el).dim for el in ok}
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
        return _select_elements(self.writer.domain, sel, self.ctx)

    def _what(self) -> str:
        if self.var in FORCE_COMPONENTS:
            return f"la fuerza interna {self.var!r}"
        if self.var in ("stress", "strain"):
            return f"{'esfuerzos' if self.var == 'stress' else 'deformaciones'} por punto de integración"
        if self.var in INTERFACE_VARS:
            word = "el salto" if self.var == "jump" else "la tracción"
            return (f"{word} de interfaz con la componente {self.component!r} (las interfaces de "
                    f"superficie tienen n, s en 2D y n, s1, s2 en 3D; la adherencia de una "
                    f"barra, n, s en 2D y n1, n2, s en 3D)")
        return f"la variable interna escalar {self.var!r}"

    def _supports(self, el) -> bool:
        if self.var in FORCE_COMPONENTS:
            try:
                ef = el.internal_forces(np.zeros(self.writer.domain.total_dofs))
            except Exception:  # noqa: BLE001 — un elemento sin el contrato no lo admite
                return False
            return ef is not None and self.var in ef.components
        if point_quantity(el).kind == JUMP:
            # Interfaces (ADR 0043 §7): salto y tracción locales por punto y
            # las variables escalares de su ley. JUMP implica LocalJumpElement
            # (lo exige el barrido de elementos).
            if self.var in INTERFACE_VARS:
                return self.component in el.jump_components()
            schema = getattr(getattr(el, "material", None), "STATE_SCHEMA", None) or {}
            return tuple(schema.get(self.var, (None,))) == ()
        if self.var in INTERFACE_VARS:
            return False
        # Sólo elementos de desplazamientos con esfuerzos por punto de integración
        # (sólidos y barras): no los térmicos (flujo), ni las interfaces
        # (salto y tracción en ejes locales), ni los marcos (resultantes de
        # sección). Clasificación única: solidum.core.point_quantity.
        mechanical = (point_quantity(el).kind == STRAIN
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
                # La del paso: constante de las etapas anteriores más factor
                # veces la que crece (ADR 0026 §5, 0027).
                ef = el.internal_forces(
                    U, equivalent_load=self.writer.element_equivalent_load(el.id, factor))
                out.append(float(ef.components[self.var][self.end]))
            return np.array(out)
        if self.var == "strain":
            eps = np.vstack([np.atleast_2d(el.compute_gauss_state(U)["strain"])
                             for el in self.elements])
            return _component(eps, self.component, strain=True)
        if self.var in INTERFACE_VARS:
            # El salto, de U (cinemática, como la deformación); la tracción,
            # del estado con que el análisis equilibró el paso.
            out = []
            for el in self.elements:
                k = el.jump_components().index(self.component)
                if self.var == "jump":
                    out.append(np.asarray(el.jumps(el.get_local_displacements(U)))[:, k])
                else:
                    st = el.state
                    out.append([float(np.atleast_1d(t)[k])
                                for t in (st.stresses_trial if trial else st.stresses)])
            return np.concatenate([np.asarray(v, dtype=float).reshape(-1) for v in out])
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


def _select_elements(dom, sel: dict, ctx: str) -> list:
    """Elementos de un selector ``id`` | ``ids`` | ``group`` (todos sus nodos
    en el grupo físico) | ``where`` (condición sobre el centro) |
    ``interface`` (los elementos de la interfaz de ese nombre, ADR 0043 §7);
    sin selector, todos."""
    if "id" in sel:
        ids = [sel["id"]]
    elif "ids" in sel:
        ids = list(sel["ids"])
    elif "interface" in sel:
        interfaces = getattr(dom, "interfaces", None) or {}
        if sel["interface"] not in interfaces:
            raise ValueError(f"{ctx}: interfaz {sel['interface']!r} inexistente. Disponibles (bloque "
                             f"'interfaces'): {sorted(interfaces)}.")
        ids = list(interfaces[sel["interface"]])
    elif "group" in sel:
        groups = getattr(dom, "physical_groups", None) or {}
        if sel["group"] not in groups:
            raise ValueError(f"{ctx}: grupo físico {sel['group']!r} inexistente. "
                             f"Disponibles: {sorted(groups)}.")
        inside = set(groups[sel["group"]])
        ids = [eid for eid, el in dom.elements.items()
               if all(nd.id in inside for nd in el.nodes)]
        if not ids:
            raise ValueError(f"{ctx}: ningún elemento tiene todos sus nodos en el grupo "
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
            raise ValueError(f"{ctx}: ningún elemento tiene el centro donde {sel['where']!r}.")
    else:
        ids = list(dom.elements)
    missing = [i for i in ids if i not in dom.elements]
    if missing:
        raise ValueError(f"{ctx}: elementos inexistentes {missing[:10]}.")
    return [dom.elements[i] for i in ids]


class _ZoneEnergy:
    """Energías de un conjunto de elementos (``energy`` con ``elements``).

    Con la fuerza interna de cada elemento ``f_e = Σ_g w_g·detJ_g·s·B_gᵀσ_g``
    (``s``: espesor en 2D) evaluada con los esfuerzos del estado con que el
    análisis equilibró el paso, las mismas fórmulas que las del modelo
    completo restringidas a la zona:

    - trabajo ``W_n = W_{n−1} + ½ Σ_e (f_e,n + f_e,n−1)·(u_e,n − u_e,n−1)``:
      el que recibe el material de la zona, ``∫∫σ:dε dV``; para todo el
      modelo es el trabajo externo;
    - energía elástica ``E_n = ½ Σ_e u_e,n·f_e,n`` (materiales que descargan
      por la secante);
    - energía disipada ``D_n = W_n − E_n``.

    La suma sobre una partición del modelo da las energías del modelo
    completo (hasta el redondeo). Sólo elementos cuya fuerza interna es
    ``∫Bᵀσ`` con la cinemática compilada (``BATCH_KINEMATICS``), los sólidos
    isoparamétricos, y los de interfaz, cuya fuerza es ``f_e = G_e·t_e``
    (``LocalJumpElement.force_matrix``) y cuya energía almacenada la da la
    ley (``stored_energy``, ADR 0043 §5) en vez de ``½·u_e·f_e``.
    """

    def __init__(self, writer, sel: dict, ctx: str, which: str):
        from solidum.math.batch.family import _has_own_state
        self.which = which
        elements = _select_elements(writer.domain, sel, ctx)
        bad = [el.id for el in elements if not _is_interface(el) and (
               getattr(type(el), "BATCH_KINEMATICS", None) is None
               or getattr(el, "points", None) is None or _has_own_state(type(el))
               or {"rx", "ry", "rz", "T"} & set(getattr(el, "DOF_NAMES", ())))]
        large = [el.id for el in elements
                 if _is_interface(el) and getattr(el, "GEOMETRICALLY_NONLINEAR", False)]
        if large:
            raise ValueError(f"{ctx}: las energías por zona usan la matriz G_e de "
                             f"desplazamientos pequeños de las interfaces; con desplazamientos "
                             f"grandes no están disponibles (elementos {large[:10]}).")
        if bad:
            raise ValueError(f"{ctx}: la energía por zona necesita sólidos isoparamétricos (fuerza "
                             f"interna ∫Bᵀσ) o elementos de interfaz; no lo son los elementos "
                             f"{bad[:10]}.")
        known = []
        if which != "external_work":
            known, unknown = _stored_energy_split(elements, writer.domain.total_dofs)
            if unknown:
                raise ValueError(
                    f"{ctx}: la energía {which!r} necesita la energía almacenada de cada elemento "
                    f"de la zona. Los de interfaz la toman de su ley; los sólidos usan ½·u·f_int, "
                    f"que sólo es la energía elástica recuperable si sus materiales descargan por "
                    f"la secante. No se conoce para {unknown[:6]}. El trabajo ('external_work') "
                    f"sí está disponible.")
        groups: dict = {}
        for el in elements:
            if _is_interface(el):
                key = (type(el), el.force_matrix().shape)
            else:
                pts, w = el.batch_quadrature()
                key = (type(el), pts.shape, len(el.get_global_dof_indices()))
            groups.setdefault(key, []).append(el)
        known_ids = {id(el) for el in known}
        self.groups = [_ZoneGroup(members, known_ids) for members in groups.values()]
        for g in self.groups:
            if np.any(g.stresses(trial=False) != 0.0):
                raise ValueError(f"{ctx}: la zona tiene esfuerzos iniciales: sus energías no se "
                                 f"pueden obtener de ½·u·f_int.")
        self.W = 0.0

    def value(self, U: np.ndarray, trial: bool) -> float:
        E = 0.0
        for g in self.groups:
            f = np.einsum("mij,mj->mi", g.G, g.stresses(trial))
            u = U[g.dofs]
            self.W += 0.5 * float(np.sum((f + g.f_prev) * (u - g.u_prev)))
            if g.stored.any():
                # Los que entregan su energía almacenada (ADR 0043 §5).
                E += 0.5 * float(np.sum(u[~g.stored] * f[~g.stored]))
                E += sum(el.stored_energy(U, trial)
                         for el, s in zip(g.members, g.stored) if s)
            else:
                E += 0.5 * float(np.sum(u * f))
            g.f_prev, g.u_prev = f, u
        return {"external_work": self.W, "elastic": E, "dissipated": self.W - E}[self.which]


def _is_interface(el) -> bool:
    """Elemento de salto local con ``f_e = G_e·t_e`` (ADR 0043): ``JUMP``
    implica ``LocalJumpElement`` (lo exige el barrido de elementos)."""
    return point_quantity(el).kind == JUMP


class _ZoneGroup:
    """Elementos de una zona con la misma clase y regla de integración:
    ``f_e = G_e·σ_e`` con ``G_e = [w_g·detJ_g·s·B_gᵀ]`` fijo (desplazamientos
    pequeños) y ``σ_e`` los esfuerzos de sus puntos en fila; en una interfaz,
    ``G_e = [m_g·(W·R·B)_gᵀ]`` y ``σ_e`` las tracciones. ``stored[m]``: el
    miembro ``m`` entrega su energía almacenada."""

    def __init__(self, members: list, stored_ids=frozenset()):
        from solidum.math.batch.postprocess import family_of
        self.members = members
        self.stored = np.array([id(el) in stored_ids for el in members], dtype=bool)
        el0 = members[0]
        cls = type(el0)
        self.dofs = np.array([el.get_global_dof_indices() for el in members], dtype=np.int64)
        n_dof = self.dofs.shape[1]
        if _is_interface(el0):
            G = np.array([el.force_matrix() for el in members])
        else:
            pts, w = el0.batch_quadrature()
            dim = pts.shape[1]
            n_sig = int(el0.STRAIN_DIM)
            scale_attr = getattr(cls, "BATCH_SCALE", None)
            B = np.zeros((n_sig, n_dof))
            G = np.zeros((len(members), n_dof, len(pts) * n_sig))
            for m, el in enumerate(members):
                coords = np.ascontiguousarray(el.get_coordinate_matrix(ndim=dim), dtype=np.float64)
                s = 1.0 if scale_attr is None else float(getattr(el, scale_attr))
                # La regla de cada elemento (dos de la misma clase pueden tener
                # reglas distintas con el mismo número de puntos).
                pts_e, w_e = el.batch_quadrature()
                for g, (pt, wg) in enumerate(zip(pts_e, w_e)):
                    detJ = cls.BATCH_KINEMATICS(np.ascontiguousarray(pt), coords, B)
                    G[m, :, g * n_sig:(g + 1) * n_sig] = (wg * detJ * s) * B.T
        self.G = G
        self.n_flat = G.shape[2]
        self.f_prev = np.zeros((len(members), n_dof))
        self.u_prev = np.zeros((len(members), n_dof))
        self._family_of = family_of
        self._plan = None

    def _gather_plan(self):
        """``[(familia, índices en la zona, filas de la familia)]`` y los
        elementos sueltos; se rehace si una familia dejó de estar adoptada."""
        plan = self._plan
        # Las familias se crean con el primer ensamblaje: un elemento suelto
        # al construir la zona puede pertenecer después a una.
        if (plan is not None and all(f.adopted for f, _, _ in plan[0])
                and all(self._family_of(self.members[m]) is None for m in plan[1])):
            return plan
        fams: dict = {}
        loose = []
        for m, el in enumerate(self.members):
            fam = self._family_of(el)
            if fam is None:
                loose.append(m)
            else:
                fams.setdefault(fam, []).append(m)
        fam_plan = []
        for fam, idx in fams.items():
            index = {id(e): i for i, e in enumerate(fam.elements)}
            rows = np.concatenate([np.arange(fam.state.rows(index[id(self.members[m])]).start,
                                             fam.state.rows(index[id(self.members[m])]).stop)
                                   for m in idx])
            fam_plan.append((fam, np.array(idx, dtype=np.int64), rows))
        self._plan = (fam_plan, loose)
        return self._plan

    def stresses(self, trial: bool) -> np.ndarray:
        fam_plan, loose = self._gather_plan()
        out = np.zeros((len(self.members), self.n_flat))
        for fam, idx, rows in fam_plan:
            sig = fam.state.sig_trial if trial else fam.state.sig_committed
            out[idx] = sig[rows].reshape(len(idx), self.n_flat)
        for m in loose:
            st = self.members[m].state
            s = st.stresses_trial if trial else st.stresses
            out[m] = np.concatenate([np.atleast_1d(np.asarray(x, dtype=float)) for x in s])
        return out


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
