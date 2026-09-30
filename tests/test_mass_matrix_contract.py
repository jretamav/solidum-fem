"""Contrato opcional de masa de un elemento (ADR 0009).

``Element.compute_mass_matrix`` es opcional: la base lanza
``NotImplementedError``. ``Assembler.assemble_mass_matrix`` no espera a ese
fallo a mitad de ensamblaje: antes de ensamblar recorre el dominio y reporta en
una sola ``ValueError`` qué tipos de elemento no implementan la masa y qué
materiales no declaran ``density``, para que un análisis modal o dinámico no
se interrumpa con un mensaje críptico.

Se prueba con una barra de juguete sin registrar que renuncia a su matriz de
masa, para no añadir clases a los barridos de contrato.
"""
import pytest

from solidum.core.domain import Domain
from solidum.core.element import Element
from solidum.elements.truss import Truss2D
from solidum.materials.elastic import Elastic1D
from solidum.math.assembly import Assembler

E, A, RHO = 200.0e9, 1.0e-3, 7850.0


class _BarraSinMasa(Truss2D):
    """``Truss2D`` que vuelve a la ``compute_mass_matrix`` de la base."""
    compute_mass_matrix = Element.compute_mass_matrix


def _dos_barras(segunda, material_segunda):
    """Barra ``Truss2D`` con densidad entre los nodos 1-2 y la barra
    ``segunda`` entre los nodos 2-3."""
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [1.0, 0.0])
    n3 = dom.add_node(3, [2.0, 0.0])
    dom.add_element(Truss2D(1, [n1, n2], Elastic1D(E=E, density=RHO), A=A))
    dom.add_element(segunda(2, [n2, n3], material_segunda, A=A))
    n1.fix_dof("ux", 0.0)
    for n in (n1, n2, n3):
        n.fix_dof("uy", 0.0)
    dom.generate_equation_numbers(verbose=False)
    return dom


def test_la_base_lanza_not_implemented():
    dom = _dos_barras(_BarraSinMasa, Elastic1D(E=E, density=RHO))
    with pytest.raises(NotImplementedError, match="_BarraSinMasa.*compute_mass_matrix"):
        dom.elements[2].compute_mass_matrix()


def test_el_ensamblador_nombra_al_elemento_sin_masa():
    dom = _dos_barras(_BarraSinMasa, Elastic1D(E=E, density=RHO))
    with pytest.raises(ValueError) as exc:
        Assembler(dom).assemble_mass_matrix()
    texto = str(exc.value)
    assert "Elementos sin `compute_mass_matrix` implementado" in texto
    assert "['_BarraSinMasa']" in texto          # sólo el tipo que falla
    assert "density" not in texto                # los materiales están bien


def test_el_mensaje_agrega_masa_y_densidad():
    """Los dos problemas se reportan juntos, en una sola pasada."""
    dom = _dos_barras(_BarraSinMasa, Elastic1D(E=E))
    with pytest.raises(ValueError) as exc:
        Assembler(dom).assemble_mass_matrix()
    texto = str(exc.value)
    assert "Materiales sin `density` declarada" in texto
    assert "['Elastic1D']" in texto
    assert "Elementos sin `compute_mass_matrix` implementado" in texto
    assert "['_BarraSinMasa']" in texto


def test_con_masa_en_todos_los_elementos_se_ensambla():
    dom = _dos_barras(Truss2D, Elastic1D(E=E, density=RHO))
    M = Assembler(dom).assemble_mass_matrix()
    assert M.shape == (dom.total_dofs, dom.total_dofs)
