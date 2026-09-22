"""Ensamblaje por lotes (ADR 0014).

Agrupa los elementos en *familias* derivadas de sus contratos
declarativos, guarda el estado interno como arreglos y evalúa todos los
puntos de Gauss de una familia dentro de un único kernel compilado. El
camino por elemento (``compute_element_state`` / ``compute_state``) sigue
siendo el contrato obligatorio y la referencia física; el camino por
lotes lo reproduce a precisión de máquina y sólo se activa para los
componentes que lo declaran.

Piezas:

- :mod:`signatures` — firmas Numba de los kernels puntuales y de familia.
- :mod:`schema` — de ``STATE_SCHEMA`` a filas de arreglo.
- :mod:`state` — ``FamilyState`` y vistas compatibles con ``ElementState``.
- :mod:`kernels` — los bucles compilados sobre elementos y puntos de Gauss
  (ensamblaje serie y paralelo, post-proceso, reducción COO → CSR).
- :mod:`family` — construcción de familias, evaluación por trozos y
  post-proceso por familia (``gauss_state``).
- :mod:`postprocess` — de dominio a familias: ``gauss_states(domain, U)``.
"""
from solidum.math.batch.family import (
    Family,
    build_families,
    element_is_batchable,
    family_key,
    material_is_batchable,
)
from solidum.math.batch.postprocess import family_of, gauss_states, group_by_family
from solidum.math.batch.schema import StateSchema
from solidum.math.batch.signatures import FAMILY_SIG, GAUSS_SIG, KIN_SIG, MAT_SIG
from solidum.math.batch.state import BatchedElementState, FamilyState

__all__ = [
    "Family",
    "FamilyState",
    "BatchedElementState",
    "StateSchema",
    "build_families",
    "element_is_batchable",
    "material_is_batchable",
    "family_key",
    "family_of",
    "gauss_states",
    "group_by_family",
    "KIN_SIG",
    "MAT_SIG",
    "FAMILY_SIG",
    "GAUSS_SIG",
]
