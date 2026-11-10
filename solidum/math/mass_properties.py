"""Masa total y centro de masa del modelo, a partir de su matriz de masa.

Sirven para comprobar la entrada de masas (densidades, masas de nudo) antes de
un análisis dinámico, y para situar la resultante de las fuerzas de inercia
(excentricidades en planta, altura de la resultante).

Para cada traslación ``d`` presente en el modelo (``ux``, ``uy``, ``uz``):

- **Masa** ``m_d = r_dᵀ·M·r_d``, con ``r_d`` la traslación rígida unitaria
  (:func:`solidum.math.rigid_body.dof_unit_vector`). Incluye la masa de los
  grados de libertad restringidos: es la masa total del modelo en esa
  dirección.
- **Centro** de la masa ``d``: la línea de acción de la resultante de las
  fuerzas de inercia ante una aceleración rígida en ``d``. Sus coordenadas
  perpendiculares a ``d`` salen de los momentos de esas fuerzas, ``r_dᵀ·M·r_k``
  con ``r_k`` el giro rígido alrededor del eje ``k`` (perpendicular a ``d`` y
  a la coordenada buscada): ``(e_k × x)_d = ε_dkj·x_j`` da
  ``x̄_j = c_j + ε_dkj·(r_dᵀ·M·r_k(c))/m_d``. Es exacto con masa consistente o
  concentrada, también en barras inclinadas, porque los elementos reproducen
  exactamente los movimientos de cuerpo rígido. La coordenada a lo largo de
  ``d`` no afecta a la línea de acción y queda en ``NaN``.

Si la masa es la misma en todas las direcciones (densidades y masas de nudo
iguales en ``ux``, ``uy``, ``uz``), las estimaciones de cada coordenada
coinciden y ``MassProperties.center`` es el centro de masa. Si no (una masa
de nudo sólo en ``ux``), cada dirección tiene su centro y ``center`` es
``None``.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from solidum.core.dof_names import ROTATIONS, TRANSLATIONS
from solidum.math.rigid_body import dof_unit_vector, rigid_body_vector

# Coincidencia de las estimaciones de una coordenada del centro, relativa a la
# dimensión del modelo.
_CENTER_RTOL = 1.0e-9


def _levi_civita(i: int, j: int, k: int) -> float:
    return float((i - j) * (j - k) * (k - i)) / 2.0


@dataclass(frozen=True)
class MassProperties:
    """Masa total y centro de masa por dirección de traslación.

    Attributes
    ----------
    total_mass : dict
        ``{"ux": m_x, "uy": m_y[, "uz": m_z]}``: masa total en cada dirección.
    center_of_mass : dict
        ``{"ux": array, ...}``: punto de la línea de acción de la resultante de
        inercia de cada dirección, con tantas coordenadas como el modelo; la
        coordenada a lo largo de la dirección es ``NaN``. Todo ``NaN`` si la
        dirección no tiene masa.
    center : np.ndarray or None
        Centro de masa del modelo, si las estimaciones de cada coordenada
        coinciden entre direcciones (a ``1e-9`` del tamaño del modelo);
        ``None`` si la masa depende de la dirección.
    """

    total_mass: dict
    center_of_mass: dict
    center: np.ndarray | None


def _isotropic(domain, total: dict) -> bool:
    """``True`` si la masa no depende de la dirección. La de los elementos no
    depende (la densidad es escalar); sí una masa de nudo distinta en ``ux``,
    ``uy``, ``uz`` de un mismo nudo, o masas totales distintas."""
    masses = list(total.values())
    if masses and np.ptp(masses) > _CENTER_RTOL * max(masses):
        return False
    per_node: dict[int, dict] = {}
    for (node_id, dof), value in getattr(domain, "nodal_masses", {}).items():
        if dof in TRANSLATIONS:
            per_node.setdefault(node_id, {})[dof] = float(value)
    for node_id, values in per_node.items():
        present = [d for d in TRANSLATIONS if d in domain.nodes[node_id].dofs]
        if len({values.get(d, 0.0) for d in present}) > 1:
            return False
    return True


def _common_center(total: dict, centers: dict, dim: int, size: float):
    """Centro único: cada coordenada, la media de sus estimaciones si
    coinciden (en 3D cada coordenada tiene dos)."""
    points = [p for d, p in centers.items() if total[d] > 0.0]
    if not points:
        return None
    stacked = np.vstack(points)
    center = np.empty(dim)
    for j in range(dim):
        values = stacked[:, j][np.isfinite(stacked[:, j])]
        if values.size == 0 or np.ptp(values) > _CENTER_RTOL * size:
            return None
        center[j] = float(np.mean(values))
    return center


def mass_properties(assembler, lumping: str = "consistent") -> MassProperties:
    """Masa total y centro de masa del modelo con la matriz de masa
    ``assembler.assemble_mass_matrix(lumping)`` (incluye las masas de nudo).

    Ver el docstring del módulo para las definiciones.
    """
    domain = assembler.domain
    if domain.total_dofs == 0:
        domain.generate_equation_numbers()
    ndof = domain.total_dofs
    M = assembler.assemble_mass_matrix(lumping=lumping)
    nodes = [n for n in domain.nodes.values() if n.dofs]
    dim = max((len(n.coordinates) for n in nodes), default=0)
    coords = np.zeros((len(nodes), dim))
    for i, n in enumerate(nodes):
        coords[i, :len(n.coordinates)] = n.coordinates
    c = coords.mean(axis=0) if len(nodes) else np.zeros(dim)
    available = ("rz",) if dim == 2 else ROTATIONS

    rotations = {}
    total, centers = {}, {}
    for d_idx, d in enumerate(TRANSLATIONS):
        if not any(d in n.dofs for n in nodes):
            continue
        r_d = dof_unit_vector(domain, d, ndof)
        Mr = M @ r_d
        m_d = float(r_d @ Mr)
        total[d] = m_d
        center = np.full(dim, np.nan)
        if m_d > 0.0:
            for j in range(dim):
                if j == d_idx:
                    continue
                k = 3 - d_idx - j
                motion = ROTATIONS[k]
                if motion not in available:
                    continue
                if motion not in rotations:
                    rotations[motion] = rigid_body_vector(domain, motion, about=c, ndof=ndof)
                center[j] = c[j] + _levi_civita(d_idx, k, j) * float(Mr @ rotations[motion]) / m_d
        centers[d] = center
    size = float(np.max(np.ptp(coords, axis=0))) if len(nodes) else 0.0
    center = (_common_center(total, centers, dim, max(size, 1e-300))
              if _isotropic(domain, total) else None)
    return MassProperties(total_mass=total, center_of_mass=centers, center=center)
