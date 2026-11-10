"""Movimientos de cuerpo rígido sobre los grados de libertad del modelo.

- :func:`dof_unit_vector`: 1 en el grado de libertad ``dof_name`` de cada nodo
  que lo tiene, 0 en los demás.
- :func:`rigid_body_vector`: traslación unitaria (``ux``, ``uy``, ``uz``) o giro
  unitario (``rx``, ``ry``, ``rz``) del modelo entero alrededor de un eje
  paralelo a ``x``, ``y`` o ``z`` que pasa por el punto ``about``.
- :func:`rotation_column`: la cinemática del giro, compartida con los modos de
  cuerpo rígido del precondicionador AMG y de la detección de mecanismos
  (:mod:`solidum.math.linalg.nullspace`, ADR 0018 y 0019).

Son los vectores de influencia ``r`` de la excitación por movimiento uniforme
del suelo (Chopra 2017, §9.4 y §13.1): el desplazamiento de cada grado de
libertad cuando el modelo se mueve como cuerpo rígido una unidad. Los apoyos se
mueven con el suelo, así que ``r`` vale lo mismo en los grados de libertad
restringidos. Los factores de participación modal son ``Γₙ = φₙᵀ·M·r``
(:func:`solidum.math.modal_response.participation_factors`).

Un giro ``θ = 1`` alrededor del eje ``e`` que pasa por ``c`` desplaza el punto
``x`` en ``e × (x − c)`` y gira los grados de libertad de rotación ``θ·e``. Los
elementos estándar reproducen exactamente los movimientos de cuerpo rígido, así
que ``K·r = 0`` en un modelo sin apoyos.
"""
from __future__ import annotations

import numpy as np

from solidum.core.dof_names import ROTATIONS, TRANSLATIONS

# Giro alrededor de cada eje: [(DOF de traslación, signo, índice de la
# coordenada que multiplica)]. Sale de u = e_a × r:
#   e_x × r = ( 0, −z,  y)
#   e_y × r = ( z,  0, −x)
#   e_z × r = (−y,  x,  0)
ROTATION_FIELDS = {
    "rx": (("uy", -1.0, 2), ("uz", +1.0, 1)),
    "ry": (("ux", +1.0, 2), ("uz", -1.0, 0)),
    "rz": (("ux", -1.0, 1), ("uy", +1.0, 0)),
}


def dof_index_arrays(nodes) -> dict:
    """``{nombre: (índices de ecuación, posición del nodo en nodes)}`` en una
    sola pasada por los nodos."""
    eq: dict[str, list[int]] = {}
    idx: dict[str, list[int]] = {}
    for i, node in enumerate(nodes):
        for name, k in node.dofs.items():
            eq.setdefault(name, []).append(k)
            idx.setdefault(name, []).append(i)
    return {name: (np.asarray(eq[name], dtype=np.int64), np.asarray(idx[name], dtype=np.int64))
            for name in eq}


def node_coordinates(nodes) -> np.ndarray:
    """Coordenadas de los nodos completadas a 3 componentes, shape ``(n, 3)``."""
    coords = np.zeros((len(nodes), 3))
    for i, node in enumerate(nodes):
        c = np.asarray(node.coordinates, dtype=float).ravel()[:3]
        coords[i, : c.size] = c
    return coords


def rotation_column(arrays: dict, coords: np.ndarray, center: np.ndarray,
                    motion: str, ndof: int) -> np.ndarray:
    """Giro unitario ``motion`` (``"rx"``, ``"ry"`` o ``"rz"``) alrededor de un
    eje que pasa por ``center`` (3 componentes): ``e × (x − c)`` en las
    traslaciones presentes y 1 en el propio giro, si existe."""
    col = np.zeros(ndof)
    for dof_name, sign, k in ROTATION_FIELDS[motion]:
        if dof_name in arrays:
            eq, idx = arrays[dof_name]
            col[eq] = sign * (coords[idx, k] - center[k])
    if motion in arrays:
        col[arrays[motion][0]] = 1.0
    return col


def dof_unit_vector(domain, dof_name: str, ndof: int | None = None) -> np.ndarray:
    """Vector con 1 en el grado de libertad ``dof_name`` de cada nodo que lo
    tiene (restringido o no) y 0 en los demás.

    Para una traslación (``ux``, ``uy``, ``uz``) es el movimiento de cuerpo
    rígido unitario en esa dirección.

    Raises
    ------
    ValueError
        Si ningún nodo tiene ese grado de libertad.
    """
    ndof = domain.total_dofs if ndof is None else int(ndof)
    r = np.zeros(ndof)
    for node in domain.nodes.values():
        idx = node.dofs.get(dof_name)
        if idx is not None and idx >= 0:
            r[idx] = 1.0
    if not r.any():
        raise ValueError(f"ningún nodo tiene el grado de libertad {dof_name!r}.")
    return r


def rigid_body_vector(domain, motion: str, *, about=None,
                      ndof: int | None = None) -> np.ndarray:
    """Movimiento de cuerpo rígido unitario del modelo.

    Parameters
    ----------
    domain
        Dominio con los grados de libertad ya numerados.
    motion : {"ux", "uy", "uz", "rx", "ry", "rz"}
        Traslación unitaria en la dirección indicada, o giro unitario (1 rad)
        alrededor de un eje paralelo a ``x``, ``y`` o ``z``. En un modelo
        plano sólo existe el giro ``rz``.
    about : array-like, optional
        Punto por el que pasa el eje del giro, con tantas coordenadas como el
        modelo. Obligatorio para un giro; no se admite con una traslación.
    ndof : int, optional
        Tamaño del vector; por omisión, ``domain.total_dofs``.

    Returns
    -------
    np.ndarray, shape (ndof,)
        ``r``: el desplazamiento de cada grado de libertad.
    """
    if motion in TRANSLATIONS:
        if about is not None:
            raise ValueError(
                f"rigid_body_vector: una traslación ({motion!r}) no tiene punto de "
                f"giro; 'about' sólo se da con un giro ({list(ROTATIONS)}).")
        return dof_unit_vector(domain, motion, ndof)
    if motion not in ROTATIONS:
        raise ValueError(
            f"rigid_body_vector: movimiento {motion!r} desconocido; se admiten "
            f"{list(TRANSLATIONS + ROTATIONS)}.")
    if about is None:
        raise ValueError(
            f"rigid_body_vector: un giro ({motion!r}) necesita el punto por el que "
            f"pasa su eje ('about').")

    nodes = [n for n in domain.nodes.values() if n.dofs]
    dim = max((len(n.coordinates) for n in nodes), default=0)
    center = np.asarray(about, dtype=float).ravel()
    if center.size != dim:
        raise ValueError(
            f"rigid_body_vector: 'about' debe tener {dim} coordenadas (las del "
            f"modelo); tiene {center.size}.")
    if dim == 2 and motion != "rz":
        raise ValueError(
            f"rigid_body_vector: en un modelo plano sólo existe el giro 'rz' "
            f"(alrededor del eje z); se pidió {motion!r}.")

    center3 = np.zeros(3)
    center3[:dim] = center
    ndof = domain.total_dofs if ndof is None else int(ndof)
    r = rotation_column(dof_index_arrays(nodes), node_coordinates(nodes), center3, motion, ndof)
    if not r.any():
        raise ValueError(
            f"rigid_body_vector: el giro {motion!r} no mueve ningún grado de "
            f"libertad del modelo.")
    return r
