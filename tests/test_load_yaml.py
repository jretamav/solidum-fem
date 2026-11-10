"""``load_yaml``: el modelo de un YAML, listo para resolver y consultar.

Criterios: ``load_yaml(path).solve()`` da exactamente lo mismo que
``run_yaml(path)`` en los pipelines estático, modal y por etapas; el
``YamlModel`` conserva el dominio (grados de libertad numerados) y el
ensamblador con el que se resolvió, y ``domain.last_result`` es el resultado.
"""
import logging
import textwrap

import numpy as np
import pytest

import solidum
from solidum import StagedResult, YamlModel, load_yaml, run_yaml

_MODEL = """\
nodes:
  - {id: 1, coords: [0.0, 0.0]}
  - {id: 2, coords: [0.0, 1.5]}
  - {id: 3, coords: [0.0, 3.0]}
materials:
  - {id: 1, type: Elastic1D, E: 2.0e11, density: 7850.0}
elements:
  - {id: 1, type: Frame2DEuler, nodes: [1, 2], material: 1, A: 6.0e-3, I: 8.5e-5}
  - {id: 2, type: Frame2DEuler, nodes: [2, 3], material: 1, A: 6.0e-3, I: 8.5e-5}
boundary_conditions_by_node:
  - {node_id: 1, ux: 0.0, uy: 0.0, rz: 0.0}
"""

_CASES = {
    "estatico": _MODEL + """\
point_loads:
  - {node_id: 3, ux: 1.0e4, uy: -2.0e4}
solver: {type: LinearSolver}
""",
    "modal": _MODEL + """\
nodal_masses:
  - {node_id: 3, ux: 2.0e3, uy: 2.0e3}
solver: {type: ModalSolver, n_modes: 3}
""",
    "etapas": _MODEL + """\
stages:
  - name: vertical
    loads:
      point_loads:
        - {node_id: 3, uy: -2.0e4}
    solver: {type: NonlinearSolver, num_steps: 2}
  - name: lateral
    loads:
      point_loads:
        - {node_id: 3, ux: 1.0e4}
    solver: {type: NonlinearSolver, num_steps: 2}
""",
}


@pytest.fixture(autouse=True)
def _silencio():
    logging.disable(logging.WARNING)
    yield
    logging.disable(logging.NOTSET)


def _write(tmp_path, text):
    path = tmp_path / "modelo.yaml"
    path.write_text(textwrap.dedent(text), encoding="utf-8")
    return path


def _key_array(res):
    if isinstance(res, StagedResult):
        return np.concatenate([r.U for r in res])
    if isinstance(res, solidum.ModalResult):
        return res.frequencies_rad
    return res.U


@pytest.mark.parametrize("case", sorted(_CASES))
def test_igual_que_run_yaml(tmp_path, case):
    path = _write(tmp_path, _CASES[case])
    model = load_yaml(path)
    assert isinstance(model, YamlModel)
    assert model.domain.total_dofs == 9
    assert model.assembler.domain is model.domain
    res = model.solve()
    assert model.domain.last_result is (res.final if isinstance(res, StagedResult) else res)
    ref = run_yaml(path)
    assert type(res) is type(ref)
    np.testing.assert_array_equal(_key_array(res), _key_array(ref))


def test_el_modelo_se_consulta_despues_de_resolver(tmp_path):
    model = load_yaml(_write(tmp_path, _CASES["estatico"]))
    res = model.solve()
    top = model.domain.nodes[3]
    assert res.U[top.dofs["ux"]] > 0.0
    M = model.assembler.assemble_mass_matrix()
    assert M.shape == (9, 9)
