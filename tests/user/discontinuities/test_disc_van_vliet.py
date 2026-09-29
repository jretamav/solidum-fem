"""Validación: ensayo de tracción de Van Vliet con grieta prescrita.

Reproduce el primer ejemplo numérico de Retama (2010), §8.1: probeta de
concreto en hueso de perro con carga excéntrica, grieta prescrita en una banda
de ``CST_Embedded2D`` a media altura y ley cohesiva lineal o exponencial. El
modelo está en ``_van_vliet.py``; la figura, en
``docs/user/discontinuities/validacion/van_vliet_fase1.png``.

Criterio de aceptación (validado por el usuario, 2026-09-29): el pico, la
deformación en el pico y la carga a 20, 50 y 100 µm frente a los valores
leídos de la fig. 8.3 de la tesis. Las tolerancias suman la incertidumbre de
lectura (±0.3 kN en el pico, ±0.5 kN en las ramas) y las diferencias de
modelado que la tesis no permite fijar (acero, malla, rótulas).

Diferencia abierta, fuera del criterio: más allá de ~100 µm la carga de
Solidum se queda por encima de la de la tesis (1.5 frente a 0.8 kN a 200 µm
con la ley exponencial). Las dos mitades giran sobre el borde de la banda
opuesto a la carga, como en la fig. 8.4; ese elemento queda comprimido y
nunca se agrieta, y la grieta cerca del pivote apenas abre. La cola depende
de la malla (1.19 kN a 200 µm con una banda más fina).
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))

from _van_vliet import TESIS, ensayo, resumen  # noqa: E402

TOLERANCIA = {"P_max": 0.8, 20: 1.5, 50: 1.0, 100: 0.7}      # [kN]
DELTA_PICO = (5.0, 9.0)                                       # [µm]


@pytest.fixture(scope="module", params=["exponential", "linear"])
def resultado(request):
    return request.param, resumen(ensayo(request.param))


def test_carga_pico(resultado):
    ley, r = resultado
    assert abs(r["P_max"] - TESIS[ley]["P_max"]) <= TOLERANCIA["P_max"], r


def test_deformacion_en_el_pico(resultado):
    ley, r = resultado
    assert DELTA_PICO[0] <= r["d_pico"] <= DELTA_PICO[1], r


@pytest.mark.parametrize("delta", [20, 50, 100])
def test_rama_de_ablandamiento(resultado, delta):
    ley, r = resultado
    assert abs(r[delta] - TESIS[ley][delta]) <= TOLERANCIA[delta], (delta, r)
