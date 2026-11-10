"""Masa total y centro de masa del modelo (``solidum.math.mass_properties``) y
masas de nudo en el VTK.

Para cada traslación ``d``, la masa es ``m_d = r_dᵀ·M·r_d`` y el centro es la
línea de acción de la resultante de las fuerzas de inercia ante una
aceleración rígida en ``d``: sus coordenadas perpendiculares a ``d`` salen de
``r_dᵀ·M·r_giro``, exactas porque los elementos reproducen los giros rígidos.
Criterios:

1. Sólo masas de nudo: la masa total es su suma y el centro la media
   ponderada de sus posiciones, exactamente. Con masas distintas por
   dirección, el centro se da por dirección y no hay un punto único.
2. Marco plano con una barra inclinada y masa consistente y concentrada, y
   marco espacial: masa y centro iguales a los de las barras (ρ·A·L en el
   punto medio de cada una), a precisión de máquina.
3. Dos placas de densidad distinta (sólidos 2D): centro analítico.
4. Las masas de nudo se exportan al VTK (``Nodal_Mass``, ``Nodal_Inertia``).
"""
import math
import os
import tempfile

import numpy as np
import pytest

from solidum import Assembler, Domain
from solidum.elements.frame import Frame2DEuler
from solidum.elements.frame3d import Frame3D
from solidum.elements.solid_2d import Quad4
from solidum.materials.elastic import Elastic1D
from solidum.materials.elastic_2d import Elastic2D
from solidum.math.mass_properties import mass_properties

RHO, A, I = 7850.0, 6.0e-3, 8.5e-5
RTOL = 1e-12


def _bars_reference(segments, rho_a):
    """Masa y centro de un conjunto de barras rectas (masa ρA·L en el punto medio)."""
    masses = np.array([rho_a * np.linalg.norm(np.subtract(b, a)) for a, b in segments])
    mids = np.array([0.5 * (np.asarray(a) + np.asarray(b)) for a, b in segments])
    return masses.sum(), masses @ mids / masses.sum()


def test_masas_de_nudo():
    dom = Domain()
    mat = Elastic1D(E=2.0e11, density=0.0)
    pts = [(0.0, 0.0), (0.0, 3.0), (4.0, 3.0), (4.0, 6.5)]
    nodes = [dom.add_node(i + 1, list(p)) for i, p in enumerate(pts)]
    for k in range(3):
        dom.add_element(Frame2DEuler(k + 1, [nodes[k], nodes[k + 1]], material=mat, A=A, I=I))
    m = {2: 1.2e4, 3: 0.8e4, 4: 2.5e4}
    for nid, value in m.items():
        dom.add_nodal_mass(nid, ux=value, uy=value)
    dom.generate_equation_numbers()
    props = mass_properties(Assembler(dom))
    total = sum(m.values())
    xbar = sum(v * pts[n - 1][0] for n, v in m.items()) / total
    ybar = sum(v * pts[n - 1][1] for n, v in m.items()) / total
    assert props.total_mass == pytest.approx({"ux": total, "uy": total}, rel=RTOL)
    # Línea de acción: la coordenada a lo largo de la dirección no está definida.
    assert math.isnan(props.center_of_mass["ux"][0])
    assert props.center_of_mass["ux"][1] == pytest.approx(ybar, rel=RTOL)
    assert props.center_of_mass["uy"][0] == pytest.approx(xbar, rel=RTOL)
    np.testing.assert_allclose(props.center, [xbar, ybar], rtol=RTOL)

    # Masa sólo en ux en el nudo 4: cada dirección tiene su centro.
    dom.add_nodal_mass(4, ux=1.0e4)
    props = mass_properties(Assembler(dom))
    assert props.total_mass["ux"] == pytest.approx(total + 1.0e4, rel=RTOL)
    assert props.center is None


@pytest.mark.parametrize("lumping", ["consistent", "lumped"])
def test_marco_plano_con_barra_inclinada(lumping):
    dom = Domain()
    mat = Elastic1D(E=2.0e11, density=RHO)
    pts = [(0.0, 0.0), (0.0, 3.2), (2.5, 4.4), (5.0, 3.2), (5.0, 0.0)]
    nodes = [dom.add_node(i + 1, list(p)) for i, p in enumerate(pts)]
    pairs = [(0, 1), (1, 2), (2, 3), (3, 4), (0, 2)]
    for k, (a, b) in enumerate(pairs, start=1):
        dom.add_element(Frame2DEuler(k, [nodes[a], nodes[b]], material=mat, A=A, I=I))
    for dof in ("ux", "uy", "rz"):
        nodes[0].fix_dof(dof, 0.0)
    dom.generate_equation_numbers()
    total, center = _bars_reference([(pts[a], pts[b]) for a, b in pairs], RHO * A)
    props = mass_properties(Assembler(dom), lumping=lumping)
    assert props.total_mass["ux"] == pytest.approx(total, rel=RTOL)
    assert props.total_mass["uy"] == pytest.approx(total, rel=RTOL)
    np.testing.assert_allclose(props.center, center, rtol=RTOL)


def test_marco_espacial():
    dom = Domain()
    mat = Elastic1D(E=2.0e11, density=RHO)
    pts = [(0.0, 0.0, 0.0), (0.0, 0.0, 3.0), (4.0, 1.0, 3.5), (4.5, 5.0, 3.0), (1.0, 4.0, 0.0)]
    nodes = [dom.add_node(i + 1, list(p)) for i, p in enumerate(pts)]
    pairs = [(0, 1), (1, 2), (2, 3), (3, 4), (1, 3)]
    for k, (a, b) in enumerate(pairs, start=1):
        dom.add_element(Frame3D(k, [nodes[a], nodes[b]], mat, A=A, Iy=I, Iz=0.6 * I,
                                J=1.2 * I, nu=0.3))
    dom.generate_equation_numbers()
    total, center = _bars_reference([(pts[a], pts[b]) for a, b in pairs], RHO * A)
    props = mass_properties(Assembler(dom))
    for d in ("ux", "uy", "uz"):
        assert props.total_mass[d] == pytest.approx(total, rel=RTOL)
    np.testing.assert_allclose(props.center, center, rtol=RTOL)


def test_placas_de_densidad_distinta():
    dom = Domain()
    rho1, rho2, t = 2400.0, 7850.0, 0.2
    pts = [(0.0, 0.0), (2.0, 0.0), (5.0, 0.0), (0.0, 1.5), (2.0, 1.5), (5.0, 1.5)]
    n = [dom.add_node(i + 1, list(p)) for i, p in enumerate(pts)]
    dom.add_element(Quad4(1, [n[0], n[1], n[4], n[3]], Elastic2D(E=3e10, nu=0.2, density=rho1),
                          thickness=t))
    dom.add_element(Quad4(2, [n[1], n[2], n[5], n[4]], Elastic2D(E=2e11, nu=0.3, density=rho2),
                          thickness=t))
    dom.generate_equation_numbers()
    m1, m2 = rho1 * 2.0 * 1.5 * t, rho2 * 3.0 * 1.5 * t
    center = (m1 * np.array([1.0, 0.75]) + m2 * np.array([3.5, 0.75])) / (m1 + m2)
    for lumping in ("consistent", "lumped"):
        props = mass_properties(Assembler(dom), lumping=lumping)
        assert props.total_mass["ux"] == pytest.approx(m1 + m2, rel=RTOL)
        np.testing.assert_allclose(props.center, center, rtol=RTOL)


def test_masas_de_nudo_en_el_vtk():
    meshio = pytest.importorskip("meshio")
    from solidum.utils.vtk_exporter import VtkExporter

    dom = Domain()
    mat = Elastic1D(E=2.0e11, density=RHO)
    n1, n2 = dom.add_node(1, [0.0, 0.0]), dom.add_node(2, [0.0, 3.0])
    dom.add_element(Frame2DEuler(1, [n1, n2], material=mat, A=A, I=I))
    dom.add_nodal_mass(2, ux=1.5e3, uy=1.2e3, rz=40.0)
    dom.generate_equation_numbers()
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "masas.vtu")
        VtkExporter(dom).export(path, U=np.zeros(dom.total_dofs))
        mesh = meshio.read(path)
    np.testing.assert_array_equal(mesh.point_data["Nodal_Mass"][1], [1.5e3, 1.2e3, 0.0])
    np.testing.assert_array_equal(mesh.point_data["Nodal_Mass"][0], [0.0, 0.0, 0.0])
    np.testing.assert_array_equal(mesh.point_data["Nodal_Inertia"][1], [0.0, 0.0, 40.0])
