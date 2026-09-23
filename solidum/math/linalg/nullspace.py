"""Modos de cuerpo rígido: casi-núcleo de ``K`` para el precondicionador AMG (ADR 0018).

Un multimalla algebraico por agregación suavizada (AMG) construye sus mallas
gruesas a partir de los vectores que el operador casi anula. En un problema
escalar (conducción) basta el vector constante, que es lo que AMG supone por
defecto. En elasticidad no: el núcleo de ``K`` sin apoyos son los **modos de
cuerpo rígido** —3 en el plano, 6 en el espacio—, y sin ellos AMG aproxima
mal los modos suaves de rotación. Medido sobre ``Hex8 40³`` (ADR 0018): con
los modos, CG converge en 11 iteraciones; sin ellos, en 112, más lento que
sin precondicionar.

Los modos se derivan **del nombre de los DOF**, sin información del tipo de
elemento, así que cubren cualquier combinación del catálogo:

- ``ux``, ``uy``, ``uz`` — una traslación unitaria por componente presente.
- Rotación infinitesimal alrededor del eje ``a``: ``u = e_a × (x − c)`` en
  los DOF de traslación y ``r_a = 1`` en el DOF de giro, si existe. Se
  incluye cuando el giro deja rastro en el modelo: cuando están presentes
  los dos desplazamientos perpendiculares al eje o el propio giro ``r_a``.
  En un modelo plano ``(ux, uy[, rz])`` sólo aparece la rotación alrededor
  de ``z``.
- Cualquier otro DOF (``T`` y futuros campos escalares): un modo constante
  por nombre. Es el núcleo del laplaciano, el de un problema de difusión sin
  condiciones de Dirichlet.

En un dominio mixto (mecánico + térmico) los bloques no se mezclan: los modos
mecánicos valen cero en ``T`` y el constante térmico vale cero en ``u``.

Son el núcleo **exacto** de ``K`` en la configuración de referencia de un
modelo sin apoyos (lo verifica ``tests/test_linalg_iterative.py`` sobre
sólidos 2D y 3D, armaduras, marcos 2D y 3D, térmico y dominio mixto); con
no linealidad o apoyos son un casi-núcleo, que es lo que AMG necesita.

Las coordenadas se centran en el centroide y cada modo de rotación se divide
por la dimensión característica del modelo, para que todas las columnas
tengan magnitud comparable con independencia de las unidades del usuario.

Además de AMG, los modos alimentan la detección de **mecanismos** antes de un
análisis estático (ADR 0019): un modo rígido que los apoyos no restringen
significa que la rigidez es singular. Por eso :func:`rigid_body_basis` guarda
junto a la matriz qué es cada columna (traslación, giro o campo escalar), el
centro y la escala, para poder describir el movimiento libre en palabras.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

_TRANSLATIONS = ("ux", "uy", "uz")
_ROTATIONS = ("rx", "ry", "rz")
_AXES = ("x", "y", "z")

# Rotación alrededor de cada eje: (DOF de giro, [(DOF de traslación, signo,
# índice de la coordenada que multiplica)]). Sale de u = e_a × r:
#   e_x × r = ( 0, −z,  y)
#   e_y × r = ( z,  0, −x)
#   e_z × r = (−y,  x,  0)
_ROTATION_FIELDS = {
    "x": ("rx", (("uy", -1.0, 2), ("uz", +1.0, 1))),
    "y": ("ry", (("ux", +1.0, 2), ("uz", -1.0, 0))),
    "z": ("rz", (("ux", -1.0, 1), ("uy", +1.0, 0))),
}


@dataclass(frozen=True)
class RigidBasis:
    """Modos de cuerpo rígido y lo necesario para interpretarlos.

    Attributes
    ----------
    B : np.ndarray, shape ``(total_dofs, m)``
        Una columna por modo.
    labels : tuple[tuple[str, str], ...]
        Qué es cada columna: ``("translation", "x")``, ``("rotation", "z")``
        o ``("scalar", "T")``.
    center : np.ndarray, shape ``(3,)``
        Centroide de los nodos: punto respecto al que giran los modos.
    length : float
        Dimensión característica por la que se dividen los modos de giro.
    """

    B: np.ndarray
    labels: tuple
    center: np.ndarray
    length: float


def _dof_arrays(nodes) -> dict:
    """``{nombre: (índices de ecuación, posición del nodo)}`` en una sola
    pasada por los nodos."""
    eq: dict[str, list[int]] = {}
    idx: dict[str, list[int]] = {}
    for i, node in enumerate(nodes):
        for name, k in node.dofs.items():
            eq.setdefault(name, []).append(k)
            idx.setdefault(name, []).append(i)
    return {name: (np.asarray(eq[name], dtype=np.int64), np.asarray(idx[name], dtype=np.int64))
            for name in eq}


def rigid_body_basis(domain) -> RigidBasis:
    """Modos de cuerpo rígido del dominio en el espacio **completo** de DOF,
    con sus etiquetas, centro y escala (ver el docstring del módulo).

    Parameters
    ----------
    domain
        Dominio con la numeración de ecuaciones ya generada. Se usan sólo
        ``node.coordinates`` y ``node.dofs``.
    """
    ndof = int(domain.total_dofs)
    nodes = [n for n in domain.nodes.values() if n.dofs]
    if ndof == 0 or not nodes:
        return RigidBasis(np.zeros((ndof, 0)), (), np.zeros(3), 1.0)

    coords = np.zeros((len(nodes), 3))
    for i, node in enumerate(nodes):
        c = np.asarray(node.coordinates, dtype=float).ravel()[:3]
        coords[i, : c.size] = c
    center = coords.mean(axis=0)
    rel = coords - center
    length = float(np.max(np.linalg.norm(rel, axis=1)))
    if length <= 0.0:
        length = 1.0

    arrays = _dof_arrays(nodes)
    present = set(arrays)

    columns: list[np.ndarray] = []
    labels: list[tuple[str, str]] = []

    for axis, name in zip(_AXES, _TRANSLATIONS):
        if name in present:
            col = np.zeros(ndof)
            col[arrays[name][0]] = 1.0
            columns.append(col)
            labels.append(("translation", axis))

    for axis in _AXES:
        rot_dof, parts = _ROTATION_FIELDS[axis]
        leaves_trace = rot_dof in present or all(p[0] in present for p in parts)
        if not leaves_trace:
            continue
        col = np.zeros(ndof)
        for dof_name, sign, coord_idx in parts:
            if dof_name in present:
                eq, idx = arrays[dof_name]
                col[eq] = sign * rel[idx, coord_idx] / length
        if rot_dof in present:
            col[arrays[rot_dof][0]] = 1.0 / length
        columns.append(col)
        labels.append(("rotation", axis))

    mechanical = set(_TRANSLATIONS) | set(_ROTATIONS)
    for name in sorted(present - mechanical):
        col = np.zeros(ndof)
        col[arrays[name][0]] = 1.0
        columns.append(col)
        labels.append(("scalar", name))

    B = np.column_stack(columns) if columns else np.zeros((ndof, 0))
    return RigidBasis(B, tuple(labels), center, length)


def rigid_body_modes(domain) -> np.ndarray:
    """Modos de cuerpo rígido del dominio en el espacio **completo** de DOF.

    Parameters
    ----------
    domain
        Dominio con la numeración de ecuaciones ya generada
        (``domain.total_dofs > 0``). Se usan sólo ``node.coordinates`` y
        ``node.dofs``.

    Returns
    -------
    np.ndarray, shape ``(total_dofs, m)``
        Una columna por modo; ``m = 0`` si el dominio no tiene DOF.
    """
    return rigid_body_basis(domain).B
