"""Tipo `ElementForces` — esfuerzos internos de un elemento estructural 1D en ejes locales.

Vive en ``solidum.core`` porque es parte del contrato de :class:`Element`
(ADR 0002): cada subclase de elemento **estructural 1D** (truss, cable,
frame 2D/3D) devuelve un :class:`ElementForces` desde
``internal_forces(U)``. Mantenerlo aquí en lugar de en ``solidum.results``
evita que ``core/element.py`` importe la capa de resultados (capa más
alta), preservando la jerarquía conceptual ``core → ... → results``.

**Dominio de aplicación (ADR 0012)**: este contrato aplica **solo a
elementos estructurales 1D**. Los sólidos (2D y 3D) **no** lo exponen:
en un continuo el equivalente de "esfuerzo seccional" es el campo
tensorial ``σ(x)``, accesible vía
:meth:`Element.compute_gauss_state(U)` (devuelve ``{stress, strain,
points_natural, points_global}`` por punto de integración). Convertir
``σ(x)`` a un valor único por elemento (promedio, centroide, valor
representativo) introduce ambigüedad — cada consumidor escoge la
agregación que necesite a partir de ``compute_gauss_state``.

``solidum.results`` lo re-exporta como parte de la API pública, así
``from solidum.results import ElementForces`` sigue funcionando para
consumidores externos.

Convenciones de signo: ``Reglas.md §5``.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np


ElementKind = Literal["truss", "cable", "frame2d", "frame3d"]

# Claves válidas por familia. Un componente ausente del dict significa "no aplica"
# a ese tipo de elemento (p. ej. no hay M en un truss) — no se rellena con ceros.
_VALID_COMPONENTS: dict[ElementKind, frozenset[str]] = {
    "truss":   frozenset({"N"}),
    "cable":   frozenset({"N"}),
    "frame2d": frozenset({"N", "V", "M"}),
    "frame3d": frozenset({"N", "Vy", "Vz", "T", "My", "Mz"}),
}


@dataclass(frozen=True)
class ElementForces:
    """Esfuerzos internos en ejes locales, evaluados en los nodos i, j.

    Inmutable (shallow): ``frozen=True`` impide reasignar atributos; los
    arrays NumPy interiores siguen siendo mutables por contenido. Tratar
    como read-only por convención.

    Parameters
    ----------
    kind
        Familia del elemento. Determina qué claves pueden aparecer en ``components``.
    components
        Diccionario ``{nombre: array(2,)}`` con el valor en el nodo i (índice 0) y
        en el nodo j (índice 1). Las claves válidas dependen de ``kind``:

        - ``"truss"`` / ``"cable"``: ``{"N"}``.
        - ``"frame2d"``: ``{"N", "V", "M"}``. Convención de viga estructural
          (Reglas.md §5): ``V`` rota el diferencial en horario si es positivo,
          ``M`` positivo es sagging.
        - ``"frame3d"``: ``{"N", "Vy", "Vz", "T", "My", "Mz"}``. Convención
          stress-resultant / RHR pura (Reglas.md §5).
    """

    kind: ElementKind
    components: dict[str, np.ndarray]

    def __post_init__(self) -> None:
        valid = _VALID_COMPONENTS.get(self.kind)
        if valid is None:
            raise ValueError(f"ElementForces.kind desconocido: {self.kind!r}")

        extra = set(self.components) - valid
        if extra:
            raise ValueError(
                f"Componentes inválidos para kind={self.kind!r}: {sorted(extra)}. "
                f"Válidos: {sorted(valid)}"
            )

        for name, arr in self.components.items():
            if not isinstance(arr, np.ndarray) or arr.shape != (2,):
                raise ValueError(
                    f"components[{name!r}] debe ser np.ndarray shape (2,), "
                    f"recibido {type(arr).__name__} shape={getattr(arr, 'shape', None)}"
                )

    def at_node_i(self) -> dict[str, float]:
        """Valores en el nodo i como escalares."""
        return {k: float(v[0]) for k, v in self.components.items()}

    def at_node_j(self) -> dict[str, float]:
        """Valores en el nodo j como escalares."""
        return {k: float(v[1]) for k, v in self.components.items()}
