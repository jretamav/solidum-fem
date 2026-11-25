"""Residuo no finito en el corrector: el intento se abandona y el paso se repite más corto.

Un material cuyo problema local no converge devuelve un esfuerzo no finito
(``FixedSmearedCrack2D`` con aplastamiento, spec §8b y Q5). El corrector lo
detecta antes de resolver, abandona el intento con un motivo propio (no
«tangente singular») y el solver repite el paso más corto (bisección).
"""
import logging

import numpy as np
import pytest

import solidum  # noqa: F401
from solidum.elements.solid_2d import Quad4
from solidum.materials.elastic_2d import Elastic2D
from solidum.math.assembly import Assembler
from solidum.math.solvers import NonlinearSolver
from solidum.mesh import structured_mesh

E, NU = 1000.0, 0.0
U_FINAL = 1.0e-2
D_LIM = 4.0e-3              # incremento de ε_xx por encima del cual el punto «no converge»


class _ElasticoQueFallaConPasosGrandes(Elastic2D):
    """Elástico que devuelve un esfuerzo no finito si el incremento de ε_xx
    respecto del estado confirmado supera ``D_LIM`` (imita un problema local
    que sólo converge con pasos pequeños). Sin kernel por lotes: el camino
    por elemento llama a ``compute_state``."""

    STATE_SCHEMA = {"eps_xx": ()}
    BATCH_KERNEL = None

    def batch_kernel(self):
        return None

    def initial_state(self):
        return {"eps_xx": 0.0}

    def compute_state(self, strain, state_vars=None):
        sig, C, _ = super().compute_state(strain, None)
        e0 = 0.0 if not state_vars else float(state_vars.get("eps_xx", 0.0))
        if abs(strain[0] - e0) > D_LIM:
            sig = np.full(3, np.nan)
            C = np.full((3, 3), np.nan)
        return sig, C, {"eps_xx": float(strain[0])}


def _barra(num_steps):
    mat = _ElasticoQueFallaConPasosGrandes(E, NU)
    mesh = structured_mesh(Quad4, (1, 1), size=(1.0, 1.0), material=mat, thickness=1.0)
    dom = mesh.domain
    for nd in mesh.nodes_on("xmin"):
        nd.fix_dof("ux", 0.0)
    for nd in mesh.nodes_on("ymin"):
        nd.fix_dof("uy", 0.0)
    right = mesh.nodes_on("xmax")
    for nd in right:
        nd.fix_dof("ux", U_FINAL)
    dom.generate_equation_numbers(verbose=False)
    return dom, right


def test_residuo_no_finito_repite_el_paso_mas_corto_y_termina(caplog):
    """Con un paso de ε_xx = 10⁻² el punto falla; la bisección lo divide hasta
    pasos de ≤ 4·10⁻³ y el análisis termina en la solución elástica. El motivo
    registrado es el del residuo no finito, no una tangente singular."""
    dom, right = _barra(num_steps=1)
    asm = Assembler(dom)
    with caplog.at_level(logging.DEBUG, logger="solidum"):
        U = NonlinearSolver(asm, num_steps=1).solve(np.zeros(dom.total_dofs))
    for nd in right:
        assert U[nd.dofs["ux"]] == pytest.approx(U_FINAL, rel=1e-12)
    texto = caplog.text
    assert "el residuo no es finito" in texto and "bisección" in texto
    assert "tangente singular" not in texto


def test_residuo_no_finito_en_un_elemento_con_todos_sus_gdl_prescritos(caplog):
    """El punto que falla está en un elemento con todos sus grados de libertad
    prescritos: el NaN sólo llega a filas de reacciones, no a ‖TᵀR‖. El
    corrector mira el residuo completo y repite el paso (H2 de la
    verificación del bloque 2 del modelo 3; sin eso, el paso «convergía» con
    esfuerzos NaN y sin aviso)."""
    mat = Elastic2D(E, NU)
    mesh = structured_mesh(Quad4, (2, 1), size=(2.0, 1.0), material=mat, thickness=1.0)
    dom = mesh.domain
    falla = _ElasticoQueFallaConPasosGrandes(E, NU)
    for eid, el in list(dom.elements.items()):
        if np.mean([nd.coordinates[0] for nd in el.nodes]) < 1.0:
            dom.replace_element(Quad4(eid, el.nodes, falla, thickness=1.0))
    for nd in dom.nodes.values():
        x = nd.coordinates[0]
        nd.fix_dof("uy", 0.0)
        if abs(x) < 1e-12:
            nd.fix_dof("ux", 0.0)
        elif abs(x - 1.0) < 1e-12:
            nd.fix_dof("ux", U_FINAL)                 # el elemento izquierdo, todo prescrito
    dom.generate_equation_numbers(verbose=False)
    asm = Assembler(dom)
    with caplog.at_level(logging.DEBUG, logger="solidum"):
        U = NonlinearSolver(asm, num_steps=1).solve(np.zeros(dom.total_dofs))
    assert np.all(np.isfinite(U))
    assert "el residuo no es finito" in caplog.text and "bisección" in caplog.text


def test_el_corrector_marca_el_residuo_no_finito_con_su_motivo():
    """``CorrectorResult.failure_reason`` y la excepción tipada llevan el motivo."""
    from solidum.math.solvers.corrector import CorrectorResult
    r = CorrectorResult(converged=False, x=None, state=None, n_solves=0,
                        nonfinite_residual=True)
    assert "el residuo no es finito" in r.failure_reason()
    err = r.divergence_error(last_load_factor=0.5)
    assert "el residuo no es finito" in str(err)
