"""Mallado con Gmsh desde un archivo de geometría ``.geo`` (ADR 0035).

El YAML declara la geometría y las opciones de la malla; Solidum llama a Gmsh,
escribe el ``.msh`` y lo lee con :class:`solidum.utils.gmsh_parser.GmshParser`.
Así un estudio de malla (h = 10, 5, 2.5 mm) o el paso de primer a segundo
orden es cambiar un dato del YAML, y el modelo es reproducible.

Parámetros del ``.geo``: el archivo los declara con su valor por omisión,
``DefineConstant[ h = 10 ];``, y ``parameters`` del YAML los sobrescribe
(``gmsh.parser.setNumber`` antes de ``gmsh.merge``; ``DefineConstant`` conserva
el valor ya definido).

Gmsh es una dependencia opcional (``pip install gmsh`` o el grupo ``mesh`` de
Solidum): sólo se importa al mallar.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from solidum.logging import get_logger

_log = get_logger("mesh.gmsh")

# Claves admitidas en el bloque ``mesh`` del YAML cuando se malla desde la geometría.
GEOMETRY_KEYS = ("geometry", "parameters", "order", "element_shape", "serendipity", "output")


@dataclass
class GeometryMeshOptions:
    """Opciones de la malla generada desde un ``.geo``.

    Attributes
    ----------
    geometry : Path
        Archivo ``.geo`` de Gmsh.
    parameters : dict[str, float]
        Valores de los parámetros del ``.geo`` (``DefineConstant``).
    order : int
        1 (elementos lineales) o 2 (cuadráticos).
    element_shape : str
        ``"tri"`` (triángulos y tetraedros, por omisión) o ``"quad"``
        (recombinación en cuadriláteros; en 3D, los hexaedros los define el
        ``.geo`` con mallas transfinitas y ``Recombine``).
    serendipity : bool
        Con ``order = 2``: ``True`` (por omisión) da cuadriláteros de 8 nodos y
        hexaedros de 20 (``Quad8``, ``Hex20``); ``False``, de 9 y 27.
    output : Path
        ``.msh`` que se escribe.
    """
    geometry: Path
    output: Path
    parameters: dict = field(default_factory=dict)
    order: int = 1
    element_shape: str = "tri"
    serendipity: bool = True


def validate_geometry_block(block: dict, ctx: str = "mesh") -> list[str]:
    """Errores de forma del bloque ``mesh`` con ``geometry`` (sin abrir Gmsh)."""
    errors = []
    unknown = set(block) - set(GEOMETRY_KEYS)
    if unknown:
        errors.append(f"{ctx}: claves no admitidas {sorted(unknown)}; se admiten {list(GEOMETRY_KEYS)}.")
    geo = block.get("geometry")
    if not isinstance(geo, str) or not geo.lower().endswith(".geo"):
        errors.append(f"{ctx}.geometry: se esperaba la ruta de un archivo .geo, no {geo!r}.")
    params = block.get("parameters", {}) or {}
    if not isinstance(params, dict):
        errors.append(f"{ctx}.parameters: se esperaba un diccionario {{nombre: valor}}.")
    else:
        for k, v in params.items():
            if isinstance(v, bool) or not isinstance(v, (int, float)):
                errors.append(f"{ctx}.parameters.{k}: se esperaba un número, no {v!r}.")
    if block.get("order", 1) not in (1, 2):
        errors.append(f"{ctx}.order: se admite 1 o 2, no {block.get('order')!r}.")
    if block.get("element_shape", "tri") not in ("tri", "quad"):
        errors.append(f"{ctx}.element_shape: se admite 'tri' o 'quad', no {block.get('element_shape')!r}.")
    if not isinstance(block.get("serendipity", True), bool):
        errors.append(f"{ctx}.serendipity: se esperaba true o false.")
    out = block.get("output")
    if out is not None and (not isinstance(out, str) or not out.lower().endswith(".msh")):
        errors.append(f"{ctx}.output: se esperaba la ruta de un archivo .msh, no {out!r}.")
    return errors


def options_from_block(block: dict, base_dir: Path) -> GeometryMeshOptions:
    geo = (base_dir / block["geometry"]).resolve()
    out = block.get("output") or geo.with_suffix(".msh").name
    return GeometryMeshOptions(
        geometry=geo, output=(base_dir / out).resolve(),
        parameters={k: float(v) for k, v in (block.get("parameters") or {}).items()},
        order=int(block.get("order", 1)), element_shape=block.get("element_shape", "tri"),
        serendipity=bool(block.get("serendipity", True)))


def mesh_from_geometry(opts: GeometryMeshOptions) -> Path:
    """Malla ``opts.geometry`` con Gmsh y escribe ``opts.output``; devuelve su ruta."""
    try:
        import gmsh
    except ImportError as exc:
        raise ImportError(
            "Mallar desde un .geo necesita Gmsh: pip install gmsh (o el grupo 'mesh' de "
            "Solidum). Sin Gmsh, genere el .msh aparte y declárelo con 'mesh: archivo.msh'."
        ) from exc
    if not opts.geometry.is_file():
        raise FileNotFoundError(f"mesh.geometry: no existe {opts.geometry}.")
    texto = opts.geometry.read_text(encoding="utf-8", errors="replace")
    ausentes = [k for k in opts.parameters if not re.search(rf"\b{re.escape(k)}\b", texto)]
    if ausentes:
        raise ValueError(f"mesh.parameters: {ausentes} no aparecen en {opts.geometry.name} "
                         "(¿errata? el .geo debe declararlos con DefineConstant).")
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
        for name, value in opts.parameters.items():
            gmsh.parser.setNumber(name, [value])
        gmsh.merge(str(opts.geometry))
        gmsh.model.geo.synchronize()
        if opts.element_shape == "quad":
            gmsh.option.setNumber("Mesh.RecombineAll", 1)
            gmsh.option.setNumber("Mesh.Algorithm", 8)          # Frontal-Delaunay para cuadriláteros
            gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)
        dim = max((d for d, _ in gmsh.model.getEntities()), default=0)
        if dim < 2:
            raise ValueError(f"mesh.geometry: {opts.geometry.name} no define superficies ni volúmenes.")
        gmsh.model.mesh.generate(dim)
        if opts.order == 2:
            gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 1 if opts.serendipity else 0)
            gmsh.model.mesh.setOrder(2)
        opts.output.parent.mkdir(parents=True, exist_ok=True)
        gmsh.write(str(opts.output))
        used = {k: gmsh.parser.getNumber(k) for k in opts.parameters}
        _log.info(f"Malla generada con Gmsh {gmsh.__version__} desde {opts.geometry.name} "
                  f"(orden {opts.order}, {opts.element_shape}, parámetros {used}) -> {opts.output}")
    finally:
        gmsh.finalize()
    return opts.output
