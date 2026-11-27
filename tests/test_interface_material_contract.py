"""Contrato de la familia de materiales de interfaz sobre todo el registro
(ADR 0033 y ADR 0043 §4).

``STATE_SCHEMA`` lo leen la historia de la salida y la exportación por punto:
una ley que no lo declara, o que devuelve variables que no declara, dejaría sus
variables internas fuera de la salida sin aviso. El barrido recorre el
registro, así que la ley N+1 entra el día que se registra.
"""
import numpy as np
import pytest

import solidum  # noqa: F401 — puebla los registros
from solidum.registry import InterfaceMaterialRegistry

from _programa_principal import nombres_programa_principal

FABRICAS = {
    "BondSlipModelCode": dict(tau_max=10e6, s1=1e-3, s2=2e-3, s3=5e-3, alpha=0.4, tau_f=2e6,
                              k0=5e10, k_n=1e12),
    "BondSlipMultilinear": dict(points=[[1e-4, 5e6], [1e-3, 8e6]], k_n=1e12),
    "CohesiveDamageIsotropic": dict(sigma_t0=3e6, G_f=100.0, K_e=1e13, softening="linear"),
    "CohesiveDamageMixedMode": dict(sigma_t0=3e6, G_f=100.0, K_e=1e13, beta=1.0,
                                    softening="exponential"),
    "ContactPenalty": dict(K_n=1e12, mu=0.3, K_t=1e12),
}


def nombres():
    return nombres_programa_principal(InterfaceMaterialRegistry)


def test_toda_ley_registrada_tiene_fabrica():
    assert set(nombres()) == set(FABRICAS)


def _saltos(n):
    """Apertura, cierre con deslizamiento y descarga: recorre las ramas."""
    base = [np.r_[1e-5, np.zeros(n - 1)], np.r_[-1e-6, np.full(n - 1, 2e-4)],
            np.r_[4e-4, np.full(n - 1, -3e-4)], np.zeros(n)]
    return base


@pytest.mark.parametrize("nombre", sorted(FABRICAS))
def test_state_schema_declarado_y_coincide_con_el_estado(nombre):
    ley = InterfaceMaterialRegistry.create(nombre, **FABRICAS[nombre])
    schema = ley.STATE_SCHEMA
    assert isinstance(schema, dict) and schema, f"{nombre}: STATE_SCHEMA sin declarar"
    probadas = 0
    for n in (2, 3):
        if not ley.accepts(ley.GEOMETRY, n):
            continue
        probadas += 1
        estado = None
        for salto in _saltos(n):
            _, _, estado = ley.compute_traction(salto, estado)
            assert set(estado) == set(schema), (nombre, n, sorted(estado), sorted(schema))
            for clave, forma in schema.items():
                valor = np.asarray(estado[clave])
                if forma == ():
                    assert valor.ndim == 0, (nombre, clave, valor)
                else:
                    assert valor.ndim == 1, (nombre, clave, valor)
    assert probadas, f"{nombre}: no admite 2 ni 3 componentes"


@pytest.mark.parametrize("nombre", sorted(FABRICAS))
def test_una_ley_sin_rigidez_abierta_y_con_rigidez_cerrada_declara_is_unilateral(nombre):
    """Comprobación por comportamiento (revisión del arquitecto del
    2026-11-26): si la rigidez normal es nula con una apertura pequeña y no
    nula al cerrar, la ley es unilateral y debe declararlo. De esa
    declaración depende que los análisis de rigidez fija rechacen el
    contacto (``ensure_no_unilateral_interfaces``)."""
    ley = InterfaceMaterialRegistry.create(nombre, **FABRICAS[nombre])
    for n in (2, 3):
        if not ley.accepts(ley.GEOMETRY, n):
            continue
        _, T_abierta, _ = ley.compute_traction(np.r_[1e-9, np.zeros(n - 1)], None)
        _, T_cerrada, _ = ley.compute_traction(np.r_[-1e-9, np.zeros(n - 1)], None)
        if T_abierta[0, 0] == 0.0 and T_cerrada[0, 0] > 0.0:
            assert ley.IS_UNILATERAL, f"{nombre}: unilateral sin declarar IS_UNILATERAL"
        elif ley.IS_UNILATERAL:
            assert T_abierta[0, 0] == 0.0, f"{nombre}: declara IS_UNILATERAL con rigidez abierta"


@pytest.mark.parametrize("nombre", sorted(FABRICAS))
def test_la_variable_principal_esta_en_el_schema(nombre):
    ley = InterfaceMaterialRegistry.create(nombre, **FABRICAS[nombre])
    assert ley.PRIMARY_STATE_VAR in ley.STATE_SCHEMA
    assert ley.STATE_SCHEMA[ley.PRIMARY_STATE_VAR] == ()


# --- cada ley con cada elemento de interfaz ---------------------------------------

_CARAS = {"line": "INTERFAZ_2D", "line3": "INTERFAZ_2D6", "tri": "INTERFAZ_3D6",
          "quad": "INTERFAZ_3D8"}


def _clases_de_interfaz():
    from solidum.elements.interface import InterfaceElement
    from solidum.registry import ElementRegistry
    return sorted(n for n in nombres_programa_principal(ElementRegistry)
                  if issubclass(ElementRegistry.get(n), InterfaceElement))


@pytest.mark.parametrize("elemento", _clases_de_interfaz())
@pytest.mark.parametrize("nombre", sorted(FABRICAS))
def test_cada_ley_con_cada_elemento_se_acepta_o_se_rechaza_con_el_mensaje_de_la_familia(
        nombre, elemento):
    """Revisión del arquitecto del 2026-11-26: el barrido de elementos usa una
    sola ley. Aquí cada ley del registro con cada clase de interfaz de espesor
    nulo: una de superficie con las componentes del elemento se acepta y
    evalúa en reposo (fuerza nula, rigidez finita); las demás se rechazan con
    el mensaje de la familia."""
    import test_element_contract_sweep as barrido
    from solidum.registry import ElementRegistry
    cls = ElementRegistry.get(elemento)
    ley = InterfaceMaterialRegistry.create(nombre, **FABRICAS[nombre])
    nodos = barrido.nodos(getattr(barrido, _CARAS[cls.FACE_CELL]))
    kwargs = {"thickness": 2.0} if cls.JUMP_DIM == 2 else {}
    if ley.accepts("surface", cls.JUMP_DIM):
        el = cls(1, nodos, ley, **kwargs)
        K, F = el.compute_element_state(np.zeros(len(nodos) * cls.JUMP_DIM))
        assert np.all(np.isfinite(K)) and np.allclose(F, 0.0)
    else:
        with pytest.raises(ValueError, match="no es un material de interfaz de superficie"):
            cls(1, nodos, ley, **kwargs)
