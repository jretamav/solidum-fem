"""Aviso de resolvedor lento (decisión del usuario, 2026-11-24).

Sin Pardiso, el despachador usa SuperLU. Antes avisaba en toda matriz
simétrica definida positiva —también en una armadura de 3 nudos— y
recomendaba CHOLMOD con conda, que en Windows no se usa. Ahora sólo una
factorización que tarda más de ``SLOW_FACTORIZATION_SECONDS`` recomienda,
una vez por proceso, ``pip install pypardiso``.
"""
import logging
import warnings

import numpy as np
import pytest
import scipy.sparse as sp

import solidum.math.linalg.lu as lu_mod
from solidum.constants import SLOW_FACTORIZATION_SECONDS
from solidum.math.linalg import dispatcher
from solidum.math.linalg.base import StiffnessProperties
from solidum.math.linalg.lu import LUSolver


def _K(n=40):
    return sp.diags([-np.ones(n - 1), 4.0 * np.ones(n), -np.ones(n - 1)], [-1, 0, 1], format="csc")


@pytest.fixture
def avisos(monkeypatch):
    """Mensajes WARNING del registro de la capa algebraica, con el «una vez
    por proceso» reiniciado."""
    monkeypatch.setattr(lu_mod, "_slow_warning_emitted", False)
    registros = []

    class Recolector(logging.Handler):
        def emit(self, record):
            registros.append(record.getMessage())

    h = Recolector(level=logging.WARNING)
    logger = logging.getLogger("solidum.linalg")
    logger.addHandler(h)
    yield registros
    logger.removeHandler(h)


def test_el_umbral_es_un_segundo():
    assert SLOW_FACTORIZATION_SECONDS == 1.0


def test_factorizacion_lenta_sin_pardiso_recomienda_pip_install_pypardiso_una_vez(avisos, monkeypatch):
    monkeypatch.setattr(lu_mod, "SLOW_FACTORIZATION_SECONDS", -1.0)   # toda factorización «lenta»
    s = LUSolver(suggest_faster_backend=True)
    s.factorize(_K())
    s.factorize(_K())
    LUSolver(suggest_faster_backend=True).factorize(_K())
    assert len(avisos) == 1, avisos
    assert "pip install pypardiso" in avisos[0]
    assert "conda" not in avisos[0]
    assert "40 incógnitas" in avisos[0]


def test_factorizacion_rapida_no_avisa(avisos):
    """Con el umbral real, un modelo pequeño (como los del curso) no dice nada."""
    LUSolver(suggest_faster_backend=True).factorize(_K())
    assert avisos == []


def test_lu_pedido_o_de_respaldo_no_avisa_aunque_tarde(avisos, monkeypatch):
    monkeypatch.setattr(lu_mod, "SLOW_FACTORIZATION_SECONDS", -1.0)
    LUSolver().factorize(_K())
    assert avisos == []


@pytest.mark.parametrize("spd", [True, False])
def test_el_despachador_sin_pardiso_elige_lu_con_sugerencia_y_sin_runtimewarning(spd, monkeypatch):
    monkeypatch.setattr(dispatcher, "_HAS_PARDISO", False)
    monkeypatch.setattr(dispatcher, "_HAS_CHOLESKY", False)
    props = StiffnessProperties(is_symmetric=True, is_positive_definite=spd, size=100)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        s = dispatcher.select_solver(props)
    assert isinstance(s, LUSolver) and s.suggest_faster_backend
    assert not [x for x in w if issubclass(x.category, RuntimeWarning)], [str(x.message) for x in w]


def test_lu_pedido_en_el_yaml_no_lleva_sugerencia(monkeypatch):
    monkeypatch.setattr(dispatcher, "_HAS_PARDISO", False)
    props = StiffnessProperties(is_symmetric=True, is_positive_definite=True, size=100)
    s = dispatcher.select_solver(props, override="lu")
    assert isinstance(s, LUSolver) and not s.suggest_faster_backend
