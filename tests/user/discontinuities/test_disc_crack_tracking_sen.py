"""Regresión: la grieta de la viga SEN con la malla libre de 5 mm no rodea un
nodo (spec ``CrackTracking``, acceptance.sen_malla_5mm_sin_rodeo).

Es el modelo del validador del 2026-11-18 (placas de reparto pegadas), 300
pasos (~35 s). Con el seguimiento global anterior (commit ``4f1b58a``) la
grieta se detenía en 17 elementos hacia el paso 245 y la carga subía sin
límite: 90.4 kN en el paso 300 (224 kN en el 800). Con la prolongación recta,
pico de 45.89 kN en el paso 98, 40 elementos agrietados y 7.45 kN en el paso
300 (medido el 2026-11-19).
"""
from __future__ import annotations

import logging
import shutil
from pathlib import Path

import pytest

import solidum

DATA = Path(__file__).resolve().parent / "data" / "viga_sen"


def test_la_grieta_pasa_y_la_carga_baja_despues_del_pico(tmp_path):
    pytest.importorskip("gmsh")
    for name in ("sen_h05.yaml", "viga_sen.geo"):
        shutil.copy(DATA / name, tmp_path / name)
    level = logging.getLogger("solidum").level
    logging.getLogger("solidum").setLevel(logging.ERROR)
    try:
        model = solidum.load_yaml(str(tmp_path / "sen_h05.yaml"))
        loads = []
        model.solve(step_callback=lambda step, U, lam: loads.append(float(lam)), write_output=False)
    finally:
        logging.getLogger("solidum").setLevel(level)
    embedded = [e for e in model.domain.elements.values() if type(e).__name__ == "CST_Embedded2D"]
    cracked = sum(e.discontinuity_state is not None for e in embedded)
    p_max = max(loads)
    assert len(loads) == 300
    assert p_max == pytest.approx(45.89e3, rel=0.01)
    assert cracked >= 38
    assert loads[-1] < 0.2 * p_max
