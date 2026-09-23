"""Utilidades compartidas por los ejemplos del manual de ejemplos (Example_manual).

Cada ejemplo del manual vive en su carpeta de ``examples/`` con un ``run.py``
que sigue el mismo contrato:

- ``calcular() -> dict``: resuelve el problema con Solidum y devuelve las
  cifras que el manual cita, **junto con su referencia analítica y el error
  relativo**. El manual no contiene ninguna cifra escrita a mano: las lee de
  aquí. Un test (``tests/test_examples_manual.py``) llama a esta función y
  comprueba que los errores declarados están dentro de su tolerancia, así que
  lo que el manual afirma se verifica en cada ejecución de la suite.
- ``figuras(res, carpeta)``: dibuja las figuras a partir de ``res``.
- ``main()``: ejecuta ambas y guarda ``resultados.json`` y las figuras.

Este módulo sólo aporta lo común: estilo de las figuras (coherente entre
ejemplos y legible en pantalla), guardado en PDF vectorial —lo que usa el
manual— y en PNG —lo que muestra GitHub en el README del ejemplo—, y la
resolución de un YAML estático conservando el dominio (``solidum.run_yaml``
sólo devuelve el resultado).
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import solidum  # noqa: E402
from solidum.logging import set_log_level  # noqa: E402
from solidum.math.assembly import Assembler  # noqa: E402
from solidum.utils.yaml_parser import YamlParser  # noqa: E402

# Los ejemplos resuelven decenas de modelos (barridos de malla o esbeltez): el
# registro INFO de cada solve sepultaría el resumen. Los avisos sí se ven.
set_log_level("WARNING")

# Paleta con buen contraste en pantalla (clara u oscura) y distinguible para
# daltonismo: numérico en azul, analítico en negro, secundarios en naranja y
# verde azulado.
COLOR_FE = "#1f5fa8"
COLOR_ANALITICO = "#222222"
COLOR_SECUNDARIO = "#d95f02"
COLOR_TERCIARIO = "#1b9e77"


def estilo_figuras() -> None:
    """Estilo común de matplotlib para todas las figuras del manual."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "figure.figsize": (6.4, 3.9),
        "figure.dpi": 110,
        "font.size": 10,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,        # fuentes incrustadas como TrueType
    })


def guardar_figura(fig, carpeta: str | Path, nombre: str) -> None:
    """Guarda ``nombre.pdf`` (vectorial, para el manual) y ``nombre.png``."""
    carpeta = Path(carpeta)
    fig.savefig(carpeta / f"{nombre}.pdf")
    fig.savefig(carpeta / f"{nombre}.png", dpi=150)
    import matplotlib.pyplot as plt
    plt.close(fig)


def _a_json(v):
    if isinstance(v, (np.floating, np.integer)):
        return v.item()
    if isinstance(v, np.ndarray):
        return v.tolist()
    if isinstance(v, dict):
        return {k: _a_json(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_a_json(x) for x in v]
    return v


def guardar_resultados(carpeta: str | Path, res: dict) -> None:
    """Escribe ``resultados.json`` (lo que el manual cita)."""
    path = Path(carpeta) / "resultados.json"
    path.write_text(json.dumps(_a_json(res), indent=2, ensure_ascii=False) + "\n",
                    encoding="utf-8")


def error_relativo(numerico: float, referencia: float) -> float:
    return abs(numerico - referencia) / abs(referencia)


_NODOS_QUAD4 = np.array([[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]])


def esfuerzos_nodales(domain, estados: dict) -> dict[int, np.ndarray]:
    """Esfuerzos nodales recuperados: extrapolación desde los puntos de Gauss
    de cada elemento a sus nodos y promedio entre los elementos que comparten
    el nodo.

    El esfuerzo del método de desplazamientos es más preciso en los puntos de
    Gauss que en los nodos, pero los máximos que interesan (borde de un
    agujero, superficie de un cilindro) están en la frontera. Para el Quad4
    con 2×2 puntos se extrapola con el campo bilineal que pasa por los cuatro
    valores de Gauss —en coordenadas naturales escaladas por √3, los puntos de
    Gauss son las esquinas de un cuadrado unitario—; en elementos con un solo
    punto, el valor constante.

    ``estados`` es la salida de ``solidum.math.batch.postprocess.gauss_states``.
    Devuelve ``{node_id: [σxx, σyy, σxy]}``.
    """
    suma: dict[int, np.ndarray] = {}
    cuenta: dict[int, int] = {}
    for eid, gs in estados.items():
        elem = domain.elements[eid]
        sig = np.asarray(gs["stress"])
        if sig.shape[0] == 1:
            nodales = np.repeat(sig, len(elem.nodes), axis=0)
        elif sig.shape[0] == 4 and len(elem.nodes) == 4:
            s = np.sign(np.asarray(gs["points_natural"]))          # (4 Gauss, 2)
            r3 = np.sqrt(3.0) * _NODOS_QUAD4                        # (4 nodos, 2)
            pesos = 0.25 * (1 + r3[:, None, 0] * s[None, :, 0]) * (1 + r3[:, None, 1] * s[None, :, 1])
            nodales = pesos @ sig
        else:
            raise NotImplementedError(f"recuperación nodal no implementada para {type(elem).__name__}")
        for node, val in zip(elem.nodes, nodales):
            suma[node.id] = suma.get(node.id, 0.0) + val
            cuenta[node.id] = cuenta.get(node.id, 0) + 1
    return {nid: suma[nid] / cuenta[nid] for nid in suma}


def resolver_yaml_estatico(ruta_yaml: str | Path):
    """Resuelve un YAML estático como ``solidum.run_yaml`` y devuelve
    ``(domain, assembler, resultado)``: el ejemplo necesita el dominio para
    extraer esfuerzos por punto de Gauss o fuerzas internas por elemento."""
    ruta_yaml = Path(ruta_yaml).resolve()
    cwd = os.getcwd()
    os.chdir(ruta_yaml.parent)   # mallas relativas al YAML
    try:
        parser = YamlParser(str(ruta_yaml))
        domain = parser.parse()
        domain.generate_equation_numbers()
        assembler = Assembler(domain)
        solver = parser.get_solver(assembler)
        F_ext = (parser.get_external_forces()
                 + parser.get_body_load(assembler)
                 + parser.get_thermal_loads())
        res = solidum.run(domain, assembler=assembler, solver=solver, F_applied=F_ext)
    finally:
        os.chdir(cwd)
    return domain, assembler, res
