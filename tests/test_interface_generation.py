"""Creación de elementos de interfaz desde la malla (ADR 0043 §6): ``pair`` y
``split``, en Python y desde el bloque ``interfaces`` del YAML.

Invariantes que protegen estas pruebas:
- ``split``: nodos = originales + Σ(k − 1) (k = componentes de los sólidos
  alrededor de cada nodo al cortar el grupo); un elemento por lado separado;
  la punta no se duplica; una unión en T da 3 copias; ids deterministas; el
  lado A no depende de la orientación de las celdas del archivo de malla.
- Reasignar nodos no cambia la rigidez de ningún sólido (bit a bit), y el
  modelo separado con las copias atadas por MPC da el mismo U que el entero.
- Contabilidad: los grupos de nodos reciben las copias; cada celda de un grupo
  sigue siendo un lado de un solo elemento; los grupos ``<name>_a``/``_b``.
- Una referencia por id a un nodo duplicado (apoyo, carga, masa, MPC, historia)
  es error; un selector por grupo o coordenada toma todas las copias, y en una
  carga o masa eso es error (cada copia la recibiría entera).
- ``pair``: holgura nula o inicial, ``n`` de A hacia B; errores de pareja, de
  cara interior, de caras que no se miran, de holgura con una cohesiva y de
  caras cuadráticas; presión de contacto exacta (resortes en serie) desde el
  YAML con parámetros ≠ 1 (Reglas §6).
"""
from __future__ import annotations

import textwrap
from pathlib import Path

import numpy as np
import pytest

import solidum
from solidum.core.domain import Domain
from solidum.elements.interface import Interface2D4, Interface3D8, InterfaceElement
from solidum.elements.thermal.quad4_thermal import Quad4Thermal
from solidum.elements.truss import Truss2D
from solidum.materials.elastic import Elastic1D
from solidum.materials.elastic_2d import Elastic2D
from solidum.materials.elastic_3d import Elastic3D
from solidum.materials.interface import CohesiveDamageIsotropic, ContactPenalty
from solidum.materials.thermal_conduction import ThermalConduction
from solidum.math.assembly import Assembler
from solidum.math.solvers import NonlinearSolver
from solidum.mesh import structured_mesh
from solidum.mesh.interfaces import pair_interfaces, split_interfaces
from solidum.utils.yaml_parser import YamlParser, YamlValidationError

MAT = Elastic2D(E=1000.0, nu=0.25)
LAW = ContactPenalty(K_n=1.0e5)
THICK = 2.0


def placa(nx=4, ny=4, L=4.0, H=2.0, element="Quad4", dom=None, origin=(0.0, 0.0)):
    m = structured_mesh(element, lines=[origin[0] + np.linspace(0, L, nx + 1),
                                        origin[1] + np.linspace(0, H, ny + 1)],
                        material=MAT, thickness=THICK, domain=dom)
    return m.domain


def grupo(dom, name, pred, cells=None):
    """Grupo de nodos que cumplen ``pred(x, y)``; sin celdas, como los de una
    malla estructurada (o con las celdas dadas, como los de Gmsh)."""
    ids = sorted(i for i, n in dom.nodes.items() if pred(*map(float, n.coordinates[:2])))
    if not getattr(dom, "physical_groups", None):
        dom.physical_groups, dom.physical_group_cells = {}, {}
    dom.physical_groups[name] = ids
    dom.physical_group_cells[name] = list(cells or [])
    return ids


def at(dom, x, y):
    return [i for i, n in dom.nodes.items() if np.allclose(n.coordinates[:2], (x, y))]


def grieta_vertical(dom, y_max=1.0, x=2.0):
    return grupo(dom, "grieta", lambda X, Y: abs(X - x) < 1e-9 and Y <= y_max + 1e-9)


def interfaces(dom):
    return [e for e in dom.elements.values() if isinstance(e, InterfaceElement)]


def sides(dom):
    """frozenset(ids) → elementos sólidos que tienen ese lado."""
    out = {}
    for el in dom.elements.values():
        for local in getattr(el, "EDGE_NODES", None) or getattr(el, "FACE_NODES", None) or ():
            if isinstance(el, InterfaceElement):
                continue
            out.setdefault(frozenset(el.nodes[j].id for j in local), []).append(el)
    return out


# ---------------------------------------------------------------------------
# split
# ---------------------------------------------------------------------------

def test_split_duplica_la_boca_y_el_interior_pero_no_la_punta():
    dom = placa()
    grieta_vertical(dom)
    n0 = len(dom.nodes)
    r = split_interfaces(dom, "grieta", LAW, name="g")
    boca, medio, punta = at(dom, 2.0, 0.0), at(dom, 2.0, 0.5), at(dom, 2.0, 1.0)
    assert len(boca) == 2 and len(medio) == 2 and len(punta) == 1
    assert len(dom.nodes) == n0 + sum(len(c) - 1 for c in r.copies.values()) == n0 + 2
    assert len(r.elements) == 2 and all(type(dom.elements[e]) is Interface2D4 for e in r.elements)
    assert dom.split_nodes == r.copies


def test_cada_elemento_de_interfaz_une_las_copias_y_n_apunta_de_a_hacia_b():
    dom = placa()
    grieta_vertical(dom)
    split_interfaces(dom, "grieta", LAW, name="g")
    for e in interfaces(dom):
        XA = np.array([e.nodes[i].coordinates for i in e.A_NODES], float)
        XB = np.array([e.nodes[i].coordinates for i in e.B_NODES], float)
        np.testing.assert_array_equal(XA, XB)
        np.testing.assert_allclose(e.R[0][0], [1.0, 0.0])     # A a la izquierda
        assert e.thickness == THICK
        assert np.all(e.g0 == 0.0)
        # Cada copia de A es sólo de sólidos de la izquierda y cada copia de B,
        # sólo de la derecha (la punta, sin duplicar, es de los dos).
        copias = {i for ids in dom.split_nodes.values() for i in ids}
        for lado, nodos in ((-1.0, e.nodes[:2]), (1.0, e.nodes[2:])):
            for n in nodos:
                if n.id not in copias:
                    continue
                dueños = [el for el in dom.elements.values()
                          if not isinstance(el, InterfaceElement) and n in el.nodes]
                assert dueños
                for el in dueños:
                    xc = np.mean([float(m.coordinates[0]) for m in el.nodes])
                    assert lado * (xc - 2.0) > 0.0


def test_union_en_t_da_tres_copias():
    dom = placa(nx=4, ny=4)
    grupo(dom, "h", lambda x, y: abs(y - 1.0) < 1e-9 and 1.0 - 1e-9 <= x <= 3.0 + 1e-9)
    grupo(dom, "v", lambda x, y: abs(x - 2.0) < 1e-9 and y >= 1.0 - 1e-9)
    n0 = len(dom.nodes)
    r = split_interfaces(dom, ["h", "v"], LAW, name="t")
    assert len(at(dom, 2.0, 1.0)) == 3                      # unión de 3 ramas
    assert len(at(dom, 1.0, 1.0)) == len(at(dom, 3.0, 1.0)) == 1   # puntas
    assert len(at(dom, 2.0, 2.0)) == 2                      # boca en el borde
    assert len(dom.nodes) == n0 + 4
    assert len(r.elements) == 4


def test_ids_deterministas():
    def hecho():
        dom = placa()
        grieta_vertical(dom)
        r = split_interfaces(dom, "grieta", LAW, name="g")
        return r.copies, [[n.id for n in dom.elements[e].nodes] for e in r.elements]
    assert hecho() == hecho()


def test_la_rigidez_de_cada_solido_no_cambia_bit_a_bit():
    dom = placa()
    grieta_vertical(dom)
    antes = {i: e.compute_global_stiffness().copy() for i, e in dom.elements.items()}
    split_interfaces(dom, "grieta", LAW, name="g")
    for i, K in antes.items():
        assert np.array_equal(dom.elements[i].compute_global_stiffness(), K)


def _resolver(dom, carga):
    for n in dom.nodes.values():
        if abs(float(n.coordinates[0])) < 1e-12:
            n.fix_dof("ux", 0.0)
            n.fix_dof("uy", 0.0)
    dom.generate_equation_numbers()
    F = np.zeros(dom.total_dofs)
    for nid, (fx, fy) in carga.items():
        F[dom.nodes[nid].dofs["ux"]] += fx
        F[dom.nodes[nid].dofs["uy"]] += fy
    asm = Assembler(dom)
    return solidum.run(dom, assembler=asm, solver=NonlinearSolver(asm, num_steps=1),
                       F_applied=F).U


def test_el_modelo_separado_y_atado_por_mpc_da_el_mismo_u():
    carga = {}
    ref = placa()
    esquina = at(ref, 4.0, 2.0)[0]
    carga[esquina] = (30.0, -70.0)
    U_ref = _resolver(ref, carga)

    dom = placa()
    grieta_vertical(dom)
    r = split_interfaces(dom, "grieta", LAW, name="g")
    for p, ids in r.copies.items():
        for c in ids[1:]:
            for dof in ("ux", "uy"):
                dom.add_linear_constraint(slave=(c, dof), masters=[(p, dof)], coefficients=[1.0])
    U = _resolver(dom, carga)
    for nid, n in ref.nodes.items():
        for dof in ("ux", "uy"):
            assert U[dom.nodes[nid].dofs[dof]] == pytest.approx(U_ref[n.dofs[dof]], rel=1e-9,
                                                                abs=1e-12)
    # Sin las ataduras la grieta se abre: la prueba distingue.
    libre = placa()
    grieta_vertical(libre)
    split_interfaces(libre, "grieta", LAW, name="g")
    U_libre = _resolver(libre, carga)
    assert np.max(np.abs(U_libre[:len(U_ref)] - U_ref)) > 1e-6


def test_los_grupos_reciben_las_copias_y_las_celdas_pasan_al_dueno():
    dom = placa()
    grieta_vertical(dom)
    base = [tuple(at(dom, x, 0.0)[0] for x in (a, a + 1.0)) for a in (0.0, 1.0, 2.0, 3.0)]
    grupo(dom, "abajo", lambda x, y: abs(y) < 1e-9, cells=base)
    split_interfaces(dom, "grieta", LAW, name="g")
    boca = at(dom, 2.0, 0.0)
    assert set(boca) <= set(dom.physical_groups["abajo"])
    assert set(boca) <= set(dom.physical_groups["grieta"])
    idx = sides(dom)
    for cell in dom.physical_group_cells["abajo"]:
        assert len(idx[frozenset(cell)]) == 1
    # Los labios: un lado de un solo elemento cada uno, uno por grupo.
    for g in ("g_a", "g_b"):
        assert len(dom.physical_group_cells[g]) == 2
        for cell in dom.physical_group_cells[g]:
            assert len(idx[frozenset(cell)]) == 1
    for g in ("g_a", "g_b"):
        assert len(set(dom.physical_groups[g]) & set(boca)) == 1


def test_la_orientacion_de_las_celdas_no_cambia_la_interfaz():
    def con(celdas_al_reves):
        dom = placa()
        ids = [at(dom, 2.0, y)[0] for y in (0.0, 0.5, 1.0)]
        cells = list(zip(ids[:-1], ids[1:]))
        if celdas_al_reves:
            cells = [c[::-1] for c in cells]
        grupo(dom, "grieta", lambda x, y: abs(x - 2.0) < 1e-9 and y <= 1.0 + 1e-9, cells=cells)
        r = split_interfaces(dom, "grieta", LAW, name="g")
        return [[n.id for n in dom.elements[e].nodes] for e in r.elements]
    assert con(False) == con(True)


def test_las_celdas_de_gmsh_mandan_sobre_los_nodos():
    """Con celdas (Gmsh) sólo se separan los lados que son celdas del grupo,
    aunque otro lado tenga todos sus nodos en el grupo."""
    dom = placa()
    ids = [at(dom, 2.0, y)[0] for y in (0.0, 0.5, 1.0)]
    grupo(dom, "grieta", lambda x, y: abs(x - 2.0) < 1e-9 and y <= 1.0 + 1e-9,
          cells=[(ids[0], ids[1])])
    r = split_interfaces(dom, "grieta", LAW, name="g")
    assert len(r.elements) == 1


@pytest.mark.parametrize("pred,mensaje", [
    (lambda x, y: abs(y) < 1e-9, "contorno"),
    (lambda x, y: x <= 1.0 + 1e-9, "es una región"),
    (lambda x, y: abs(x - 2.0) < 1e-9 and 0.5 - 1e-9 <= y <= 1.0 + 1e-9, "no separa ningún nodo"),
    (lambda x, y: abs(x - 2.0) < 1e-9 and abs(y - 1.0) < 1e-9, "no contiene ningún lado"),
])
def test_split_rechaza_grupos_que_no_son_lados_interiores(pred, mensaje):
    dom = placa()
    grupo(dom, "x", pred)
    n0, e0 = len(dom.nodes), len(dom.elements)
    with pytest.raises(ValueError, match=mensaje):
        split_interfaces(dom, "x", LAW, name="g")
    assert (len(dom.nodes), len(dom.elements)) == (n0, e0)


def test_split_rechaza_una_barra_en_un_nodo_que_se_duplicaria():
    dom = placa()
    grieta_vertical(dom)
    b, c = at(dom, 2.0, 0.0)[0], at(dom, 3.0, 1.0)[0]
    dom.add_element(Truss2D(100, [dom.nodes[b], dom.nodes[c]], Elastic1D(E=200.0), A=0.1))
    with pytest.raises(ValueError, match="Truss2D 100"):
        split_interfaces(dom, "grieta", LAW, name="g")


def test_split_rechaza_elementos_termicos_en_los_nodos():
    dom = placa()
    grieta_vertical(dom)
    k = ThermalConduction(k=2.0)
    for i, el in list(dom.elements.items()):
        dom.add_element(Quad4Thermal(1000 + i, list(el.nodes), k))
    with pytest.raises(ValueError, match="térmicos"):
        split_interfaces(dom, "grieta", LAW, name="g")


def test_split_rechaza_caras_cuadraticas_3d():
    """Sin interfaz cuadrática 3D (spec Interface2D6, «Advertencias»)."""
    mat = Elastic3D(E=1000.0, nu=0.25)
    dom = structured_mesh("Hex20", lines=[np.linspace(0, 2, 3), np.linspace(0, 1, 2),
                                          np.linspace(0, 2, 3)], material=mat).domain
    ids = sorted(i for i, n in dom.nodes.items()
                 if abs(n.coordinates[0] - 1.0) < 1e-9 and n.coordinates[2] <= 1.0 + 1e-9)
    dom.physical_groups, dom.physical_group_cells = {"grieta": ids}, {"grieta": []}
    with pytest.raises(ValueError, match="no existen todavía"):
        split_interfaces(dom, "grieta", LAW, name="g")


@pytest.mark.parametrize("element", ["Quad8", "Quad9", "Tri6"])
def test_split_con_solidos_cuadraticos_crea_interface2d6(element):
    dom = placa(element=element)
    grieta_vertical(dom)
    n0 = len(dom.nodes)
    r = split_interfaces(dom, "grieta", LAW, name="g")
    # Grieta de y = 0 a 1 con lados de 0.5: 2 lados, 5 nodos (2 intermedios);
    # se duplican todos menos la punta.
    assert len(r.elements) == 2 and len(dom.nodes) == n0 + 4
    for eid in r.elements:
        e = dom.elements[eid]
        assert type(e).__name__ == "Interface2D6"
        # Convenio de la spec: A = (1, 2, 5), B frente a ella = (4, 3, 6).
        XA = np.array([e.nodes[i].coordinates for i in (0, 1, 4)], float)
        XB = np.array([e.nodes[i].coordinates for i in (3, 2, 5)], float)
        np.testing.assert_array_equal(XA, XB)
        # Esquinas primero y el intermedio en el centro; n = +x.
        np.testing.assert_allclose(XA[2], 0.5 * (XA[0] + XA[1]))
        np.testing.assert_allclose(e.R[:, 0], [[1.0, 0.0]] * 3)


def test_split_3d_con_hexaedros():
    mat = Elastic3D(E=1000.0, nu=0.25)
    m = structured_mesh("Hex8", lines=[np.linspace(0, 2, 3), np.linspace(0, 1, 2),
                                       np.linspace(0, 2, 3)], material=mat)
    dom = m.domain
    ids = sorted(i for i, n in dom.nodes.items()
                 if abs(n.coordinates[0] - 1.0) < 1e-9 and n.coordinates[2] <= 1.0 + 1e-9)
    dom.physical_groups, dom.physical_group_cells = {"grieta": ids}, {"grieta": []}
    n0 = len(dom.nodes)
    r = split_interfaces(dom, "grieta", LAW, name="g")
    assert len(r.elements) == 1 and type(dom.elements[r.elements[0]]) is Interface3D8
    assert len(dom.nodes) == n0 + 2                       # la arista de la boca (z = 0)
    e = dom.elements[r.elements[0]]
    np.testing.assert_allclose(e.R[0][0], [1.0, 0.0, 0.0], atol=1e-14)


# ---------------------------------------------------------------------------
# pair
# ---------------------------------------------------------------------------

def dos_bloques(gap=0.0, dx=0.0, element="Quad4"):
    dom = placa(nx=2, ny=1, L=2.0, H=1.0, element=element)
    placa(nx=2, ny=1, L=2.0, H=1.0, element=element, dom=dom, origin=(dx, 1.0 + gap))
    y_a, y_b = 1.0, 1.0 + gap
    # Los nodos de cada cara, por bloque (ids de la segunda malla, mayores).
    n_a = len(dom.nodes) // 2
    grupo(dom, "a", lambda x, y: abs(y - y_a) < 1e-9)
    grupo(dom, "b", lambda x, y: abs(y - y_b) < 1e-9)
    dom.physical_groups["a"] = [i for i in dom.physical_groups["a"] if i <= n_a]
    dom.physical_groups["b"] = [i for i in dom.physical_groups["b"] if i > n_a]
    return dom


@pytest.mark.parametrize("gap", [0.0, 0.1])
def test_pair_empareja_por_la_normal_de_a(gap):
    dom = dos_bloques(gap)
    r = pair_interfaces(dom, "a", "b", LAW, name="junta")
    assert len(r.elements) == 2
    for eid in r.elements:
        e = dom.elements[eid]
        np.testing.assert_allclose(e.R[0][0], [0.0, 1.0])
        np.testing.assert_allclose(e.g0, gap)
        XA = np.array([e.nodes[i].coordinates for i in e.A_NODES], float)
        XB = np.array([e.nodes[i].coordinates for i in e.B_NODES], float)
        np.testing.assert_allclose(XB - XA, [[0.0, gap]] * 2, atol=1e-14)
    assert set(dom.physical_groups["junta_a"]) == set(dom.physical_groups["a"])
    assert set(dom.physical_groups["junta_b"]) == set(dom.physical_groups["b"])


def test_pair_nodos_sin_pareja():
    dom = dos_bloques(0.1, dx=0.3)
    with pytest.raises(ValueError, match="sin pareja"):
        pair_interfaces(dom, "a", "b", LAW, name="junta")


def test_pair_rechaza_una_cara_interior():
    dom = placa()
    grupo(dom, "medio", lambda x, y: abs(y - 1.0) < 1e-9)
    grupo(dom, "arriba", lambda x, y: abs(y - 2.0) < 1e-9)
    with pytest.raises(ValueError, match="split"):
        pair_interfaces(dom, "medio", "arriba", LAW, name="junta")


def test_pair_rechaza_caras_que_no_se_miran():
    dom = dos_bloques(0.1)
    n_a = len(dom.nodes) // 2
    grupo(dom, "techo", lambda x, y: abs(y - 2.1) < 1e-9)
    dom.physical_groups["techo"] = [i for i in dom.physical_groups["techo"] if i > n_a]
    with pytest.raises(ValueError, match="no miran"):
        pair_interfaces(dom, "a", "techo", LAW, name="junta")


def test_pair_con_holgura_y_una_cohesiva_es_error():
    dom = dos_bloques(0.1)
    n_el = len(dom.elements)
    ley = CohesiveDamageIsotropic(sigma_t0=3.0, G_f=0.1, K_e=1e4, softening="linear")
    with pytest.raises(ValueError, match="no admite holgura"):
        pair_interfaces(dom, "a", "b", ley, name="junta")
    assert len(dom.elements) == n_el


def test_pair_rechaza_grupos_de_nombre_ocupado():
    dom = dos_bloques(0.1)
    grupo(dom, "junta_a", lambda x, y: y < 0.5)
    with pytest.raises(ValueError, match="ya existen"):
        pair_interfaces(dom, "a", "b", LAW, name="junta")


# ---------------------------------------------------------------------------
# Atomicidad, referencias previas, numeración y elección del elemento
# (revisión del arquitecto del 2026-11-26)
# ---------------------------------------------------------------------------

def _huella(dom):
    """Todo lo que ``split`` cambia si llega a mutar el dominio."""
    return (sorted((i, tuple(map(float, n.coordinates)), tuple(sorted(n.dofs.items())),
                    tuple(sorted(n.boundary_conditions.items()))) for i, n in dom.nodes.items()),
            {i: tuple(n.id for n in e.nodes) for i, e in dom.elements.items()},
            dom.topology_version, dict(dom.split_nodes), dict(dom.interfaces),
            {g: tuple(v) for g, v in dom.physical_groups.items()},
            {g: tuple(map(tuple, v)) for g, v in dom.physical_group_cells.items()},
            dom.total_dofs)


def _placa_cuarto_de_punto():
    """Quad8 con el intermedio de un lado de la grieta a un cuarto: la cara de
    la interfaz es degenerada (la tangente se anula en un extremo)."""
    def mapeo(X):
        x, y = X
        if abs(x - 2.0) < 1e-9 and abs(y - 0.75) < 1e-9:
            y = 0.625
        return np.array([x, y])
    m = structured_mesh("Quad8", lines=[np.linspace(0, 4, 5), np.linspace(0, 2, 5)],
                        mapping=mapeo, material=MAT, thickness=THICK)
    dom = m.domain
    grupo(dom, "grieta", lambda x, y: abs(x - 2.0) < 1e-9 and y <= 1.0 + 1e-9)
    return dom


@pytest.mark.parametrize("caso, mensaje", [
    ("integracion", "integration debe ser"),
    ("elemento_3d", "Interface3D8 une caras"),
    ("no_es_interfaz", "no es un elemento de interfaz"),
    ("apoyo_previo", "ya tienen apoyos"),
    ("mpc_previa", "restricciones multipunto"),
    ("masa_previa", "masas"),
    ("cara_degenerada", "degenerada"),
])
def test_split_fallido_no_cambia_el_dominio(caso, mensaje):
    """Cada causa de error deja el dominio como estaba: nodos, conectividad,
    ``topology_version``, ``split_nodes``, grupos y numeración."""
    if caso == "cara_degenerada":
        dom = _placa_cuarto_de_punto()
    else:
        dom = placa()
        grieta_vertical(dom)
    dom.generate_equation_numbers()
    p = at(dom, 2.0, 0.5)[0]
    kwargs = {}
    if caso == "integracion":
        kwargs["integration"] = "gaus"
    elif caso == "elemento_3d":
        kwargs["element"] = "Interface3D8"
    elif caso == "no_es_interfaz":
        kwargs["element"] = "Quad4"
    elif caso == "apoyo_previo":
        dom.nodes[p].fix_dof("ux", 0.0)
    elif caso == "mpc_previa":
        q = at(dom, 1.0, 0.5)[0]
        dom.add_linear_constraint((q, "ux"), [(p, "ux")], [1.0])
    elif caso == "masa_previa":
        dom.add_nodal_mass(p, ux=2.0)
    antes = _huella(dom)
    with pytest.raises(ValueError, match=mensaje):
        split_interfaces(dom, "grieta", LAW, name="g", **kwargs)
    assert _huella(dom) == antes


def test_split_despues_de_numerar_y_ensamblar_da_la_misma_rigidez():
    """Con el dominio numerado y un ensamblador ya construido, ``split``
    renumera y el mismo ensamblador rehace su topología: K es la del dominio
    separado antes de numerar, bit a bit, sin grados de libertad en −1."""
    tarde = placa()
    grieta_vertical(tarde)
    tarde.generate_equation_numbers()
    asm = Assembler(tarde)
    asm.assemble_system()
    n0 = asm.K_global.shape[0]
    split_interfaces(tarde, "grieta", LAW, name="g")
    assert all(v >= 0 for n in tarde.nodes.values() for v in n.dofs.values())
    asm.assemble_system()
    pronto = placa()
    grieta_vertical(pronto)
    split_interfaces(pronto, "grieta", LAW, name="g")
    pronto.generate_equation_numbers()
    ref = Assembler(pronto)
    ref.assemble_system()
    assert asm.K_global.shape[0] > n0
    assert (asm.K_global != ref.K_global).nnz == 0


def test_ensamblar_con_grados_de_libertad_sin_numerar_es_un_error_claro():
    """Un nodo añadido después de numerar queda con ecuación −1: el
    ensamblador lo dice en vez de leer en silencio el último grado de
    libertad."""
    dom = placa(nx=1, ny=1)
    dom.generate_equation_numbers()
    el = dom.elements[min(dom.elements)]
    old = el.nodes[0]
    dom.add_node(99, list(old.coordinates))
    dom.reassign_element_nodes(el.id, {old.id: 99})
    with pytest.raises(ValueError, match="sin numerar"):
        Assembler(dom).assemble_system()


def test_element_elige_la_clase_y_por_omision_la_del_programa_principal():
    dom = placa()
    grieta_vertical(dom)
    made = split_interfaces(dom, "grieta", LAW, name="g", element="Interface2D4")
    assert {type(dom.elements[i]) for i in made.elements} == {Interface2D4}
    dom2 = dos_bloques(0.0)
    made2 = pair_interfaces(dom2, "a", "b", LAW, name="j", element=Interface2D4)
    assert {type(dom2.elements[i]) for i in made2.elements} == {Interface2D4}


def test_barra_embebida_que_cruza_la_grieta_se_ata_a_los_nodos_de_su_lado(tmp_path):
    """El bloque ``interfaces`` va antes que ``constraints`` (ADR 0043 §6): cada
    punto de la barra se ata a los nodos del sólido que lo contiene, ya
    separado. Con el orden invertido, los puntos del lado B quedaban atados a
    los nodos del lado A (*medido* por el arquitecto)."""
    ruta = _yaml(tmp_path, """\
    materials:
      - {id: 1, type: Elastic2D, E: 1000.0, nu: 0.25}
      - {id: 2, type: Elastic1D, E: 2.0e5}
    interface_materials:
      - {id: 1, type: ContactPenalty, K_n: 1.0e5}
    mesh_element: Quad4
    mesh:
      structured:
        blocks:
          - {name: placa, material: 1, x: {length: 4, n: 4}, y: {length: 2, n: 4}}
        groups:
          grieta: "near(x, 2) and y < 1.0 + 1e-9"
    interfaces:
      - {name: g, split: grieta, material: 1}
    constraints:
      - {type: EmbeddedReinforcement, points: [[0.3, 0.55], [3.7, 0.55]], area: 0.01,
         material: 2}
    """)
    model = solidum.load_yaml(ruta)
    dom = model.domain
    assert dom.split_nodes
    (bar,) = model.parser.applied_constraints
    solidos = [set(n.id for n in e.nodes) for e in dom.elements.values()
               if not isinstance(e, InterfaceElement) and e.id not in bar.element_ids]
    atados = 0
    for c in dom.linear_constraints:
        if c["slave"][0] not in bar.node_ids:
            continue
        maestros = {m[0] for m in c["masters"]}
        assert any(maestros <= s for s in solidos), (c["slave"], maestros)
        atados += 1
    assert atados > 0


# ---------------------------------------------------------------------------
# Domain.reassign_element_nodes
# ---------------------------------------------------------------------------

def test_reasignar_exige_la_misma_posicion_y_cambia_la_huella():
    dom = placa(nx=1, ny=1)
    el = dom.elements[min(dom.elements)]
    old = el.nodes[0]
    dom.add_node(99, [float(c) + 0.5 for c in old.coordinates])
    with pytest.raises(ValueError, match="misma posición|no está en la posición"):
        dom.reassign_element_nodes(el.id, {old.id: 99})
    dom.add_node(98, list(old.coordinates))
    v = dom.topology_version
    dom.reassign_element_nodes(el.id, {old.id: 98})
    assert el.nodes[0] is dom.nodes[98] and dom.topology_version == v + 1
    assert set(dom.nodes[98].dofs) == {"ux", "uy"}
    with pytest.raises(ValueError, match="no tiene el nodo"):
        dom.reassign_element_nodes(el.id, {old.id: 98})


def test_reasignar_no_deja_nodos_repetidos():
    """El mapa debe ser inyectivo y no llevar un nodo a otro que el elemento
    ya tiene (revisión del arquitecto del 2026-11-26)."""
    dom = placa()
    grieta_vertical(dom)
    made = split_interfaces(dom, "grieta", LAW, name="g")
    el = dom.elements[made.elements[0]]           # caras A y B con nodos coincidentes
    a = el.nodes[el.A_NODES[0]]
    b = el.nodes[el.B_NODES[0]]
    assert a.id != b.id and list(a.coordinates) == list(b.coordinates)
    dom.add_node(9998, list(a.coordinates))
    v = dom.topology_version
    with pytest.raises(ValueError, match="mismo nuevo"):
        dom.reassign_element_nodes(el.id, {a.id: 9998, b.id: 9998})
    with pytest.raises(ValueError, match="ya tiene los nodos"):
        dom.reassign_element_nodes(el.id, {a.id: b.id})
    assert dom.topology_version == v and el.nodes[el.A_NODES[0]] is a


def test_un_solido_no_acepta_una_ley_de_interfaz():
    """Hueco anterior al bloque: la validación de compatibilidad salía antes
    de tiempo con un material sin ``STRAIN_DIM``."""
    from solidum.elements.solid_2d import Quad4
    dom = placa(nx=1, ny=1)
    with pytest.raises(ValueError, match="es una ley de interfaz"):
        Quad4(99, list(dom.elements[min(dom.elements)].nodes), LAW)


# ---------------------------------------------------------------------------
# Bloque YAML
# ---------------------------------------------------------------------------

E_B, K_N, T_B = 1000.0, 1.0e5, 3.0
BLOQUES = f"""\
    materials:
      - {{id: 1, type: Elastic2D, E: {E_B}, nu: 0.0}}
    interface_materials:
      - {{id: 1, type: ContactPenalty, K_n: {K_N}}}
      - {{id: 2, type: CohesiveDamageIsotropic, sigma_t0: 3.0, G_f: 0.1, K_e: 1.0e4,
         softening: linear}}
    mesh_element: Quad4
    mesh_thickness: {T_B}
    mesh:
      structured:
        blocks:
          - {{name: abajo, material: 1, x: {{length: 2, n: 4}}, y: {{length: 1, n: 2}}}}
          - {{name: arriba, material: 1, origin: [0, 1.1], x: {{length: 2, n: 4}},
             y: {{length: 1, n: 2}}}}
        groups:
          Origen: "near(x, 0) and near(y, 0)"
    """.rstrip(" ")


def _yaml(tmp_path, texto):
    ruta = tmp_path / "m.yaml"
    ruta.write_text(textwrap.dedent(texto), encoding="utf-8")
    return ruta


def test_presion_de_contacto_exacta_desde_el_yaml(tmp_path):
    """Holgura 0.1 y bajada 0.2: la presión es la de tres resortes en serie,
    σ = 0.1/(2/E + 1/K_n) (*calculado*), con espesor 3 en la fuerza."""
    ruta = _yaml(tmp_path, BLOQUES + """\
    interfaces:
      - {name: junta, pair: [abajo.ymax, arriba.ymin], material: 1}
    boundary_conditions_by_group:
      abajo.ymin: {uy: 0.0}
      Origen: {ux: 0.0}
      arriba.ymax: {uy: -0.2, ux: 0.0}
    solver: {type: NonlinearSolver, num_steps: 4}
    """)
    model = solidum.load_yaml(ruta)
    res = model.solve()
    sigma = 0.1 / (2.0 / E_B + 1.0 / K_N)
    (made,) = model.parser.applied_interfaces
    for eid in made.elements:
        for t in model.domain.elements[eid].state.stresses:
            np.testing.assert_allclose(t, [-sigma, 0.0], rtol=1e-9, atol=1e-9 * sigma)
    abajo = model.domain.physical_groups["abajo.ymin"]
    R = sum(res.R[model.domain.nodes[i].dofs["uy"]] for i in abajo)
    assert R == pytest.approx(sigma * 2.0 * T_B, rel=1e-9)


def test_split_desde_el_yaml_y_apoyo_por_grupo_en_las_dos_copias(tmp_path):
    ruta = _yaml(tmp_path, """\
    materials:
      - {id: 1, type: Elastic2D, E: 1000.0, nu: 0.25}
    interface_materials:
      - {id: 1, type: ContactPenalty, K_n: 1.0e5}
    mesh_element: Quad4
    mesh:
      structured:
        x: {length: 4, n: 4}
        y: {length: 2, n: 4}
        groups:
          grieta: "near(x, 2) and y <= 1"
          boca: "near(x, 2) and near(y, 0)"
    interfaces:
      - {name: g, split: grieta, material: 1, integration: gauss}
    boundary_conditions_by_group:
      boca: {uy: 0.0}
    """)
    parser = YamlParser(str(ruta))
    dom = parser.parse()
    boca = dom.physical_groups["boca"]
    assert len(boca) == 2
    assert all(dom.nodes[i].boundary_conditions == {"uy": 0.0} for i in boca)
    (made,) = parser.applied_interfaces
    assert all(dom.elements[e].integration == "gauss" for e in made.elements)


@pytest.mark.parametrize("bloque,mensaje", [
    ("boundary_conditions: [{node_id: @, ux: 0.0}]", "se duplicó"),
    ("point_loads: [{node_id: @, uy: 1.0}]", "se duplicó"),
    ("nodal_masses: [{node_id: @, ux: 1.0}]", "se duplicó"),
    ("linear_constraints: [{slave: {node: @, dof: ux}, masters: [{node: 1, dof: ux}], "
     "coefficients: [1.0]}]", "se duplicó"),
    ("point_loads_by_group: {boca: {uy: 1.0}}", "recibiría la carga entera"),
    ("point_loads_by_coord: [{coord: y, val: 0.0, uy: 1.0}]", "recibiría la carga entera"),
    ("nodal_masses_by_group: {boca: {ux: 1.0}}", "recibiría la masa entera"),
])
def test_referencias_ambiguas_a_un_nodo_duplicado(tmp_path, bloque, mensaje):
    base = """\
    materials:
      - {id: 1, type: Elastic2D, E: 1000.0, nu: 0.25}
    interface_materials:
      - {id: 1, type: ContactPenalty, K_n: 1.0e5}
    mesh_element: Quad4
    mesh:
      structured:
        x: {length: 4, n: 4}
        y: {length: 2, n: 4}
        groups:
          grieta: "near(x, 2) and y <= 1"
          boca: "near(x, 2) and near(y, 0)"
    interfaces:
      - {name: g, split: grieta, material: 1}
    """.rstrip(" ")
    # El id del nodo de la boca, leído de un modelo igual sin la entrada.
    dom = YamlParser(str(_yaml(tmp_path, base))).parse()
    boca = dom.physical_groups["boca"][0]
    ruta = _yaml(tmp_path, base + "    " + bloque.replace("@", str(boca)) + "\n")
    with pytest.raises(YamlValidationError, match=mensaje):
        YamlParser(str(ruta)).parse()


def test_la_historia_rechaza_el_id_de_un_nodo_duplicado(tmp_path):
    ruta = _yaml(tmp_path, """\
    materials:
      - {id: 1, type: Elastic2D, E: 1000.0, nu: 0.25}
    interface_materials:
      - {id: 1, type: ContactPenalty, K_n: 1.0e5}
    mesh_element: Quad4
    mesh:
      structured:
        x: {length: 4, n: 4}
        y: {length: 2, n: 4}
        groups:
          grieta: "near(x, 2) and y <= 1"
          boca: "near(x, 2) and near(y, 0)"
    interfaces:
      - {name: g, split: grieta, material: 1}
    boundary_conditions_by_group:
      xmin: {ux: 0.0, uy: 0.0}
      xmax: {ux: 0.01}
    solver: {type: NonlinearSolver, num_steps: 1}
    output:
      history:
        columns:
          - {name: u, displacement: {node: @, dof: ux}}
    """)
    texto = ruta.read_text(encoding="utf-8")
    sin = YamlParser(str(_yaml(tmp_path, texto.split("output:")[0]))).parse()
    boca = sin.physical_groups["boca"][0]
    ruta = _yaml(tmp_path, texto.replace("@", str(boca)))
    with pytest.raises(ValueError, match="se duplicó"):
        solidum.run_yaml(ruta, write_output=True)


@pytest.mark.parametrize("entrada,mensaje", [
    ("{name: j, material: 1}", "exactamente uno"),
    ("{name: j, pair: [abajo.ymax], material: 1}", "dos nombres de grupo"),
    ("{name: j, split: grieta}", "falta 'material'"),
    ("{name: j, split: grieta, material: 7}", "no está declarado"),
    ("{name: j, split: grieta, material: 1, integration: lobatto}", "integration"),
    ("{pair: [abajo.ymax, arriba.ymin], material: 1}", "falta 'name'"),
    ("{name: j, split: grieta, material: 1, espesor: 2}", "claves desconocidas"),
])
def test_forma_del_bloque_interfaces(tmp_path, entrada, mensaje):
    ruta = _yaml(tmp_path, BLOQUES + f"    interfaces:\n      - {entrada}\n")
    with pytest.raises(YamlValidationError, match=mensaje):
        YamlParser(str(ruta)).parse()


def test_holgura_con_una_cohesiva_desde_el_yaml(tmp_path):
    ruta = _yaml(tmp_path, BLOQUES + """\
    interfaces:
      - {name: junta, pair: [abajo.ymax, arriba.ymin], material: 2}
    """)
    with pytest.raises(YamlValidationError, match=r"interfaces\[0\] \(junta\).*no admite holgura"):
        YamlParser(str(ruta)).parse()


def test_un_elemento_de_interfaz_no_es_receta_de_un_grupo_de_la_malla(tmp_path):
    ruta = _yaml(tmp_path, """\
    materials:
      - {id: 1, type: Elastic2D, E: 1000.0, nu: 0.25}
    interface_materials:
      - {id: 1, type: ContactPenalty, K_n: 1.0e5}
    mesh:
      structured:
        element: Interface2D4
        material: 1
        x: {length: 4, n: 4}
        y: {length: 2, n: 4}
    """)
    with pytest.raises(YamlValidationError, match="es un elemento de interfaz"):
        YamlParser(str(ruta)).parse()


def test_interfaces_sin_malla_es_error(tmp_path):
    ruta = _yaml(tmp_path, """\
    nodes: [{id: 1, coords: [0, 0]}, {id: 2, coords: [1, 0]}]
    materials:
      - {id: 1, type: Elastic1D, E: 1.0}
    interface_materials:
      - {id: 1, type: ContactPenalty, K_n: 1.0e5}
    elements:
      - {id: 1, type: Truss2D, nodes: [1, 2], material: 1, A: 1.0}
    interfaces:
      - {name: j, split: grieta, material: 1}
    """)
    with pytest.raises(YamlValidationError, match="necesitan una malla"):
        YamlParser(str(ruta)).parse()


# ---------------------------------------------------------------------------
# Mallas de Gmsh (celdas del grupo dadas por la malla)
# ---------------------------------------------------------------------------

GMSH_DATA = (Path(__file__).resolve().parent / "data" / "gmsh")


@pytest.fixture
def gmsh():
    return pytest.importorskip("gmsh")


def test_split_de_una_curva_incrustada_en_gmsh(tmp_path, gmsh):
    ruta = _yaml(tmp_path, f"""\
    materials:
      - {{id: 1, type: Elastic2D, E: 1000.0, nu: 0.25}}
    interface_materials:
      - {{id: 1, type: ContactPenalty, K_n: 1.0e5}}
    mesh:
      geometry: {(GMSH_DATA / "grieta.geo").as_posix()}
      output: m.msh
    interfaces:
      - {{name: g, split: Grieta, material: 1}}
    """)
    parser = YamlParser(str(ruta))
    dom = parser.parse()
    (made,) = parser.applied_interfaces
    assert len(made.elements) == 2                      # h = 0.25 en 0.5 de grieta
    assert len(at(dom, 1.0, 0.0)) == 2 and len(at(dom, 1.0, 0.5)) == 1
    idx = sides(dom)
    for g in ("Grieta", "Abajo", "g_a", "g_b"):
        for cell in dom.physical_group_cells[g]:
            assert len(idx[frozenset(cell)]) == 1, g     # labios y contorno: un solo dueño
    assert len(dom.physical_group_cells["Grieta"]) == 4  # los dos labios


def test_pair_de_dos_superficies_de_gmsh_con_nodos_coincidentes(tmp_path, gmsh):
    """Compresión uniforme (ν = 0) a través de la junta: la presión es la de
    tres resortes en serie (*calculado*), con nodos coincidentes (g_0 = 0)."""
    ruta = _yaml(tmp_path, f"""\
    materials:
      - {{id: 1, type: Elastic2D, E: {E_B}, nu: 0.0}}
    interface_materials:
      - {{id: 1, type: ContactPenalty, K_n: {K_N}}}
    mesh_thickness: {T_B}
    mesh:
      geometry: {(GMSH_DATA / "bloques.geo").as_posix()}
      output: m.msh
    interfaces:
      - {{name: junta, pair: [TechoAbajo, PisoArriba], material: 1}}
    boundary_conditions_by_group:
      Base: {{uy: 0.0}}
      Origen: {{ux: 0.0}}
      Tapa: {{uy: -0.01}}
    solver: {{type: NonlinearSolver, num_steps: 2}}
    """)
    model = solidum.load_yaml(ruta)
    (made,) = model.parser.applied_interfaces
    assert len(made.elements) == 4
    model.solve()
    sigma = 0.01 / (2.0 / E_B + 1.0 / K_N)
    for eid in made.elements:
        e = model.domain.elements[eid]
        assert np.all(e.g0 == 0.0)
        for t in e.state.stresses:
            np.testing.assert_allclose(t, [-sigma, 0.0], rtol=1e-8, atol=1e-8 * sigma)


def test_espesores_distintos_a_los_dos_lados_son_error():
    """Spec Interface2D4 §6: el espesor de la interfaz 2D es el de los sólidos
    de los dos lados; si difieren, no está definido."""
    dom = dos_bloques(0.1)
    n_a = len(dom.nodes) // 2
    for el in list(dom.elements.values()):
        if min(n.id for n in el.nodes) > n_a:
            dom.replace_element(type(el)(el.id, list(el.nodes), MAT, thickness=THICK / 2))
    with pytest.raises(ValueError, match="espesores distintos"):
        pair_interfaces(dom, "a", "b", LAW, name="junta")
    split = placa()
    grieta_vertical(split)
    for el in list(split.elements.values()):
        if np.mean([float(n.coordinates[0]) for n in el.nodes]) > 2.0:
            split.replace_element(type(el)(el.id, list(el.nodes), MAT, thickness=THICK / 2))
    n0 = len(split.nodes)
    with pytest.raises(ValueError, match="espesores distintos"):
        split_interfaces(split, "grieta", LAW, name="g")
    assert len(split.nodes) == n0
