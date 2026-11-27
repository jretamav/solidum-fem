"""El contacto requiere un análisis no lineal (decisión del usuario del
2026-11-26, Diálogo Q4 de la spec ``Interface2D4``; ADR 0043).

Invariantes que protegen estas pruebas:
- Todo análisis de rigidez fija (estático lineal y su perturbación sobre una
  etapa, modal, armónico, espectral, Newmark y HHT lineales, diferencias
  centrales lineales) rechaza un modelo con una ley de interfaz unilateral,
  con un mensaje que nombra la ley y los elementos; los estáticos y los
  transitorios nombran además sus análisis no lineales.
- No rechaza lo que no es contacto: una ley de interfaz bilateral (cohesiva)
  ni un material unilateral de barra (cable).
- Los análisis no lineales siguen admitiendo el contacto.
"""
from __future__ import annotations

import inspect

import numpy as np
import pytest

import solidum
from solidum.materials.elastic_2d import Elastic2D
from solidum.materials.interface import CohesiveDamageIsotropic, ContactPenalty
from solidum.math.assembly import Assembler
from solidum.math.solvers import (
    CentralDifferenceSolver,
    LinearSolver,
    ModalSolver,
    NonlinearSolver,
)
from solidum.math.solvers.model_checks import (
    ensure_no_unilateral_interfaces,
    nonlinear_alternatives,
)
from solidum.mesh import structured_mesh
from solidum.mesh.interfaces import split_interfaces
from solidum.registry import SolverRegistry

from _programa_principal import nombres_programa_principal

MENSAJE = "el contacto requiere un análisis no lineal"


def modelo(law, density=7.8e3):
    """Placa 4 × 2 de Quad4 con una grieta interior x = 2, y ≤ 1, separada
    con la ley ``law``; apoyada abajo y con una carga a la derecha."""
    mat = Elastic2D(E=1000.0, nu=0.25, density=density)
    dom = structured_mesh("Quad4", lines=[np.linspace(0, 4, 5), np.linspace(0, 2, 5)],
                          material=mat).domain
    ids = sorted(i for i, n in dom.nodes.items()
                 if abs(n.coordinates[0] - 2.0) < 1e-9 and n.coordinates[1] <= 1.0 + 1e-9)
    dom.physical_groups, dom.physical_group_cells = {"grieta": ids}, {"grieta": []}
    made = split_interfaces(dom, "grieta", law, name="g")
    for n in dom.nodes.values():
        if abs(n.coordinates[1]) < 1e-9:
            n.boundary_conditions = {"ux": 0.0, "uy": 0.0}
    dom.generate_equation_numbers()
    asm = Assembler(dom)
    F = np.zeros(dom.total_dofs)
    for n in dom.nodes.values():
        if abs(n.coordinates[0] - 4.0) < 1e-9:
            F[n.dofs["ux"]] = -1.0
    return dom, asm, F, made


# Cómo resolver un solver de cada familia (PIPELINE_KIND) sin más datos: la
# comprobación ocurre antes de avanzar, así que basta un paso corto.
RESOLVER = {
    "static": lambda S, a, F: S(a).solve(F),
    "modal": lambda S, a, F: S(a, n_modes=2).solve(),
    "harmonic": lambda S, a, F: S(a, omega=np.array([1.0]), F_amplitude=F).solve(),
    "spectrum": lambda S, a, F: S(a, n_modes=2, direction={"dof_name": "ux"},
                                  spectrum=lambda T: 1.0).solve(),
    "transient": lambda S, a, F: S(a, 1.0e-3, 1.0e-3).solve(),
}


def _rigidez_fija():
    """Solvers del programa principal cuya instancia por omisión declara
    ``FIXED_STIFFNESS = True`` (el barrido de solvers exige que lo declaren)."""
    out = []
    for nombre in nombres_programa_principal(SolverRegistry):
        cls = SolverRegistry.get(nombre)
        if cls.PIPELINE_KIND not in RESOLVER:
            continue
        flag = inspect.getattr_static(cls, "FIXED_STIFFNESS", None)
        if flag is True or isinstance(flag, property):
            out.append(nombre)
    return out


def test_todo_solver_de_rigidez_fija_tiene_resolutor_en_la_prueba():
    """Un solver de rigidez fija de una familia nueva (p. ej. pandeo) no se
    salta en silencio: hay que añadir aquí cómo resolverlo."""
    for nombre in nombres_programa_principal(SolverRegistry):
        cls = SolverRegistry.get(nombre)
        flag = inspect.getattr_static(cls, "FIXED_STIFFNESS", None)
        if flag is True or isinstance(flag, property):
            assert cls.PIPELINE_KIND in RESOLVER, (
                f"{nombre}: añada a RESOLVER cómo resolver la familia {cls.PIPELINE_KIND!r}")


def test_la_lista_de_rigidez_fija_es_la_esperada():
    """Si un solver cambia de clase, la prueba de abajo lo cubre sola; esta
    lista sólo documenta los de hoy."""
    assert set(_rigidez_fija()) >= {"LinearSolver", "ModalSolver", "HarmonicSolver",
                                     "ResponseSpectrumSolver", "NewmarkSolver", "HHTSolver",
                                     "CentralDifferenceSolver"}


@pytest.mark.parametrize("nombre", _rigidez_fija())
def test_un_analisis_de_rigidez_fija_rechaza_el_contacto(nombre):
    """Recorre el registro: todo solver de rigidez fija rechaza el contacto,
    con un mensaje que nombra la ley, los elementos y, en estática y
    dinámica, los análisis no lineales (derivados del mismo atributo)."""
    cls = SolverRegistry.get(nombre)
    _, asm, F, made = modelo(ContactPenalty(K_n=1.0e5))
    with pytest.raises(ValueError, match=MENSAJE) as info:
        RESOLVER[cls.PIPELINE_KIND](cls, asm, F)
    texto = str(info.value)
    assert texto.startswith(nombre) and "ContactPenalty" in texto
    assert str(made.elements[0]) in texto
    alternativas = nonlinear_alternatives(cls.PIPELINE_KIND)
    if alternativas is None:
        assert "Análisis no lineales" not in texto
    else:
        assert alternativas in texto
        assert all(SolverRegistry.get(n.split()[0]).PIPELINE_KIND == cls.PIPELINE_KIND
                   for n in alternativas.replace(" o ", ", ").split(", "))


def test_la_perturbacion_lineal_sobre_una_etapa_rechaza_el_contacto():
    dom, asm, F, _ = modelo(ContactPenalty(K_n=1.0e5))
    base = solidum.run(dom, F_applied=F, solver=NonlinearSolver(asm, num_steps=1), assembler=asm)
    with pytest.raises(ValueError, match=MENSAJE):
        solidum.run(dom, F_applied=F, solver=LinearSolver(asm), assembler=asm, previous=base)


def test_los_analisis_no_lineales_admiten_el_contacto():
    dom, asm, F, made = modelo(ContactPenalty(K_n=1.0e5))
    res = solidum.run(dom, F_applied=F, solver=NonlinearSolver(asm, num_steps=2), assembler=asm)
    assert res.converged
    # La carga empuja la cara derecha hacia la grieta: se cierra.
    assert any(t[0] < 0.0 for e in made.elements for t in dom.elements[e].state.stresses)
    CentralDifferenceSolver(asm, 1.0e-4, 1.0e-4, nonlinear=True,
                            F_func=lambda t: F).solve()


def test_una_ley_de_interfaz_bilateral_se_admite():
    law = CohesiveDamageIsotropic(sigma_t0=3.0, G_f=0.1, K_e=1.0e4, softening="linear")
    dom, asm, F, _ = modelo(law)
    ensure_no_unilateral_interfaces(dom, "prueba", "static")
    U = LinearSolver(asm).solve(F)
    assert np.all(np.isfinite(U))
    ModalSolver(asm, n_modes=2).solve()


def test_un_cable_no_es_contacto():
    """``Material.IS_UNILATERAL`` (los cables) no entra en esta regla."""
    from solidum.core.domain import Domain
    from solidum.elements.cable import Cable2DCorot
    from solidum.materials.cable_1d import CableMaterial1D
    dom = Domain()
    dom.add_node(1, [0.0, 0.0])
    dom.add_node(2, [1.0, 0.0])
    dom.add_element(Cable2DCorot(1, [dom.nodes[1], dom.nodes[2]], CableMaterial1D(E=2.0e11), 1e-4))
    assert CableMaterial1D.IS_UNILATERAL
    ensure_no_unilateral_interfaces(dom, "prueba", "static")


def test_desde_el_yaml(tmp_path):
    ruta = tmp_path / "m.yaml"
    ruta.write_text("""\
materials:
  - {id: 1, type: Elastic2D, E: 1000.0, nu: 0.0}
interface_materials:
  - {id: 1, type: ContactPenalty, K_n: 1.0e5}
mesh_element: Quad4
mesh:
  structured:
    blocks:
      - {name: abajo, material: 1, x: {length: 2, n: 2}, y: {length: 1, n: 1}}
      - {name: arriba, material: 1, origin: [0, 1.1], x: {length: 2, n: 2}, y: {length: 1, n: 1}}
interfaces:
  - {name: junta, pair: [abajo.ymax, arriba.ymin], material: 1}
boundary_conditions_by_group:
  abajo.ymin: {ux: 0.0, uy: 0.0}
  arriba.ymax: {ux: 0.0, uy: -0.2}
solver: {type: LinearSolver}
""", encoding="utf-8")
    with pytest.raises(ValueError, match=MENSAJE):
        solidum.load_yaml(ruta).solve(write_output=True)
