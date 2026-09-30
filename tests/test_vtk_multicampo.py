"""Exportador VTK en dominios multicampo y con restricciones multipunto.

Auditoría de las MPC del 2026-09-30 (hallazgos TE-5, TE-6 y la parte VTK de
TE-9):

- **Suavizado nodal de σ en un dominio multicampo** ``(ux, uy, T)`` con
  ``Quad4`` y ``Quad4Thermal`` superpuestos sobre los mismos nodos (ADR 0016).
  Hasta el 2026-09-30 el suavizado filtraba las celdas sólo por tipo
  (``line``): las celdas térmicas, que comparten el tipo ``quad``, entraban
  con σ = 0 y ``Sigma_*_nodal`` salía dividido por dos. En un modelo sólo
  térmico se exportaban además campos ``Sigma_*_nodal`` nulos. Ahora sólo
  promedian las celdas de sólidos mecánicos, y el caso sólo mecánico no
  cambia bit a bit.
- **MPC** (``u_s = u_m`` entre nodos coincidentes, oráculo del nodo fundido):
  ``Displacements``, ``Temperature`` y los campos por celda del modelo atado
  son los del modelo fundido, también en los nodos esclavos. El suavizado
  nodal en la interfaz es la limitación documentada (TE-6): cada copia
  promedia sólo las celdas de su lado, y el valor fundido es la media de las
  dos cuando cada lado aporta el mismo número de celdas.
"""
import numpy as np
import pytest

import solidum
from solidum.core.domain import Domain
from solidum.elements.solid_2d.quad4 import Quad4
from solidum.elements.solid_2d.tri3 import Tri3
from solidum.elements.thermal.quad4_thermal import Quad4Thermal
from solidum.elements.truss import Truss2D
from solidum.materials.elastic import Elastic1D
from solidum.materials.elastic_2d import Elastic2D
from solidum.materials.thermal_conduction import ThermalConduction
from solidum.utils.vtk_exporter import STRESS_COMPONENT_NAMES, VtkExporter

meshio = pytest.importorskip("meshio")

NX, NY, LX, LY, IC = 4, 2, 2.0, 1.0, 2     # IC: columna de la interfaz atada
ESP = 0.2
DUP = 100000                               # desplazamiento de ids de las copias

NODAL_SIGMA = tuple(f"{n}_nodal" for n in STRESS_COMPONENT_NAMES) + ("Von_Mises_nodal",)


def _placa(*, mecanico=True, termico=True, mpc=False):
    """Placa de NX × NY celdas con ``Quad4`` y/o ``Quad4Thermal`` sobre los
    mismos nodos. Con ``mpc`` la columna ``IC`` de nodos se duplica para las
    celdas de su derecha y cada copia se ata a su original DOF a DOF
    (``u_s = u_m``, ``T_s = T_m``): cinemáticamente es el modelo fundido.

    Borde izquierdo con ``ux = 0`` y ``T = 100``, ``uy = 0`` en la esquina;
    borde derecho con carga vertical y aporte de calor.
    """
    dom = Domain()
    mat_m = Elastic2D(E=30e9, nu=0.2)
    mat_t = ThermalConduction(k=45.0, c=460.0, density=7850.0, dim=2)
    nid, dup = {}, {}
    for j in range(NY + 1):
        for i in range(NX + 1):
            n = j * (NX + 1) + i + 1
            dom.add_node(n, [LX * i / NX, LY * j / NY])
            nid[i, j] = n
    if mpc:
        for j in range(NY + 1):
            dup[IC, j] = DUP + nid[IC, j]
            dom.add_node(dup[IC, j], [LX * IC / NX, LY * j / NY])
    eid = 1
    for j in range(NY):
        for i in range(NX):
            esq = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)]
            ids = [dup[c] if (mpc and c[0] == IC and i >= IC) else nid[c] for c in esq]
            nodos = [dom.nodes[k] for k in ids]
            if mecanico:
                dom.add_element(Quad4(eid, nodos, mat_m, thickness=ESP))
                eid += 1
            if termico:
                dom.add_element(Quad4Thermal(eid, nodos, mat_t, thickness=ESP))
                eid += 1
    dofs = (("ux", "uy") if mecanico else ()) + (("T",) if termico else ())
    for j in range(NY + 1):
        if mecanico:
            dom.nodes[nid[0, j]].fix_dof("ux", 0.0)
        if termico:
            dom.nodes[nid[0, j]].fix_dof("T", 100.0)
    if mecanico:
        dom.nodes[nid[0, 0]].fix_dof("uy", 0.0)
    if mpc:
        for key, n in dup.items():
            for d in dofs:
                dom.add_linear_constraint((n, d), [(nid[key], d)], [1.0])
    dom.generate_equation_numbers(verbose=False)
    F = np.zeros(dom.total_dofs)
    for j in range(NY + 1):
        if mecanico:
            F[dom.nodes[nid[NX, j]].dofs["uy"]] = -1e5
        if termico:
            F[dom.nodes[nid[NX, j]].dofs["T"]] = 40.0
    return dom, nid, dup, F


def _exportar(dom, U, path):
    VtkExporter(dom).export(str(path), U=U)
    return meshio.read(str(path))


def _por_nodo(dom, mesh, campo):
    """``{id de nodo: valor}`` de un campo nodal (orden de ``domain.nodes``)."""
    return dict(zip(dom.nodes, mesh.point_data[campo]))


def _mapear_U(origen, U, destino):
    """Copia a ``destino`` los DOF de ``U`` (numeración de ``origen``) que
    ambos dominios comparten por (nodo, nombre de DOF)."""
    V = np.zeros(destino.total_dofs)
    for n, node in destino.nodes.items():
        for d, k in node.dofs.items():
            V[k] = U[origen.nodes[n].dofs[d]]
    return V


# ---------------------------------------------------------------------------
# TE-5: suavizado nodal de σ sólo sobre sólidos mecánicos
# ---------------------------------------------------------------------------

def test_suavizado_multicampo_igual_al_solo_mecanico(tmp_path):
    """Con térmicos superpuestos, ``Sigma_*_nodal`` es el del modelo sólo
    mecánico con el mismo campo de desplazamientos (antes salía a la mitad)."""
    mix, _, _, F = _placa()
    r = solidum.run(mix, F_applied=F)
    mec, _, _, _ = _placa(termico=False)
    U_mec = _mapear_U(mix, r.U, mec)

    m_mix = _exportar(mix, r.U, tmp_path / "mix.vtu")
    m_mec = _exportar(mec, U_mec, tmp_path / "mec.vtu")

    assert np.max(np.abs(m_mec.point_data["Sigma_XX_nodal"])) > 1e6
    for campo in NODAL_SIGMA:
        np.testing.assert_allclose(m_mix.point_data[campo], m_mec.point_data[campo],
                                   rtol=1e-12, atol=1e-6, err_msg=campo)


def test_solo_termico_no_exporta_sigma_nodal(tmp_path):
    """Un modelo sin sólidos mecánicos no exporta ``Sigma_*_nodal`` (antes,
    campos nulos); ``Temperature`` y ``Flux`` sí."""
    dom, _, _, F = _placa(mecanico=False)
    r = solidum.run(dom, F_applied=F)
    m = _exportar(dom, r.U, tmp_path / "termico.vtu")
    for campo in NODAL_SIGMA:
        assert campo not in m.point_data
    assert "Temperature" in m.point_data
    assert "Flux" in m.cell_data


def test_solo_mecanico_sin_cambios_bit_a_bit(tmp_path):
    """Sin térmicos, el suavizado es exactamente la media de orden 0 de las
    celdas no lineales de siempre, en el mismo orden de acumulación: sólidos
    de dos tipos de celda (``quad`` y ``triangle``) y una barra (``line``,
    fuera del promedio), con un campo de σ no uniforme."""
    dom = Domain()
    mat = Elastic2D(E=30e9, nu=0.2)
    xy = {1: (0, 0), 2: (1, 0), 3: (2, 0), 4: (0, 1), 5: (1, 1), 6: (2, 1), 7: (3, 1)}
    for n, c in xy.items():
        dom.add_node(n, list(map(float, c)))
    N = dom.nodes
    dom.add_element(Quad4(1, [N[1], N[2], N[5], N[4]], mat, thickness=ESP))
    dom.add_element(Tri3(2, [N[2], N[3], N[6]], mat, thickness=ESP))
    dom.add_element(Tri3(3, [N[2], N[6], N[5]], mat, thickness=ESP))
    dom.add_element(Truss2D(4, [N[6], N[7]], Elastic1D(E=200e9), A=1e-3))
    dom.generate_equation_numbers(verbose=False)
    rng = np.random.default_rng(7)
    U = 1e-4 * rng.standard_normal(dom.total_dofs)

    m = _exportar(dom, U, tmp_path / "mec.vtu")

    # Referencia: el algoritmo previo a la corrección, sobre los σ por celda.
    n_nodes = len(dom.nodes)
    acc = np.zeros((n_nodes, 6))
    cnt = np.zeros(n_nodes, dtype=np.int64)
    for k, block in enumerate(m.cells):
        if block.type == "line":
            continue
        s = np.column_stack([m.cell_data[c][k] for c in STRESS_COMPONENT_NAMES])
        conn = np.asarray(block.data, dtype=np.int64)
        np.add.at(acc, conn.ravel(), np.repeat(s, conn.shape[1], axis=0))
        cnt += np.bincount(conn.ravel(), minlength=n_nodes)
    mask = cnt > 0
    acc[mask] /= cnt[mask][:, None]

    for k, name in enumerate(STRESS_COMPONENT_NAMES):
        np.testing.assert_array_equal(m.point_data[f"{name}_nodal"], acc[:, k])
    assert np.all(m.point_data["Sigma_XX_nodal"][list(dom.nodes).index(7)] == 0.0)


# ---------------------------------------------------------------------------
# TE-9 (parte VTK) y TE-6: modelo atado por MPC frente al fundido
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def fundido_y_atado(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("vtk_mpc")
    out = {}
    for mpc in (False, True):
        dom, nid, dup, F = _placa(mpc=mpc)
        r = solidum.run(dom, F_applied=F)
        out[mpc] = (dom, nid, dup, _exportar(dom, r.U, tmp / f"mpc_{mpc}.vtu"))
    return out


def _original(n):
    return n - DUP if n > DUP else n


@pytest.mark.parametrize("campo", ["Displacements", "Temperature"])
def test_mpc_campos_nodales_iguales_al_fundido(fundido_y_atado, campo):
    """``Displacements`` y ``Temperature`` del modelo atado coinciden con los
    del fundido en todos los nodos, esclavos incluidos."""
    d0, _, _, m0 = fundido_y_atado[False]
    d1, _, dup, m1 = fundido_y_atado[True]
    v0, v1 = _por_nodo(d0, m0, campo), _por_nodo(d1, m1, campo)
    escala = max(np.max(np.abs(v)) for v in v0.values())
    assert escala > 0.0
    for n in dup.values():
        assert n in v1
    for n, v in v1.items():
        np.testing.assert_allclose(v, v0[_original(n)], rtol=0.0, atol=1e-10 * escala,
                                   err_msg=f"{campo} nodo {n}")


@pytest.mark.parametrize("campo", ["Sigma_XX", "Sigma_YY", "Tau_XY", "Von_Mises", "Flux"])
def test_mpc_campos_por_celda_iguales_al_fundido(fundido_y_atado, campo):
    """Los campos por celda (σ de los ``Quad4``, flujo de los térmicos) no
    dependen de que el nodo esté fundido o atado."""
    _, _, _, m0 = fundido_y_atado[False]
    _, _, _, m1 = fundido_y_atado[True]
    a0 = np.concatenate([np.asarray(b) for b in m0.cell_data[campo]])
    a1 = np.concatenate([np.asarray(b) for b in m1.cell_data[campo]])
    escala = np.max(np.abs(a0))
    assert escala > 0.0
    np.testing.assert_allclose(a1, a0, rtol=0.0, atol=1e-10 * escala)


def test_mpc_suavizado_nodal_limitacion_en_la_interfaz(fundido_y_atado):
    """TE-6, limitación documentada: fuera de la interfaz ``Sigma_XX_nodal``
    coincide con el fundido; en la interfaz cada copia promedia sólo las
    celdas de su lado, las dos difieren, y el valor fundido es su media (cada
    lado aporta el mismo número de celdas a cada nodo de la columna). La
    diferencia entre copias se comprueba en las fibras extremas (la fila
    central está cerca de la fibra neutra, con σ_xx ≈ 0)."""
    d0, nid, _, m0 = fundido_y_atado[False]
    d1, _, dup, m1 = fundido_y_atado[True]
    s0 = _por_nodo(d0, m0, "Sigma_XX_nodal")
    s1 = _por_nodo(d1, m1, "Sigma_XX_nodal")
    escala = max(abs(v) for v in s0.values())
    tol = 1e-10 * escala
    for (i, j), n in nid.items():
        if i != IC:
            assert abs(s1[n] - s0[n]) <= tol, f"nodo {n}"
    for j in range(NY + 1):
        m, s = nid[IC, j], dup[IC, j]
        if j in (0, NY):
            assert abs(s1[m] - s1[s]) > 1e-2 * escala, f"interfaz j={j}"
        assert abs(0.5 * (s1[m] + s1[s]) - s0[m]) <= tol, f"interfaz j={j}"


def test_mpc_transitorio_termico_temperatura_en_esclavos(tmp_path):
    """La temperatura final de un transitorio térmico con MPC se exporta en
    los nodos esclavos igual que en su maestro."""
    dom, _, dup, _ = _placa(mecanico=False, mpc=True)
    r = solidum.run_thermal_transient(dom, dt=5.0, n_steps=4, T_initial=20.0)
    m = _exportar(dom, r.T_final, tmp_path / "transitorio.vtu")
    T = _por_nodo(dom, m, "Temperature")
    assert max(abs(T[n] - 20.0) for n in T) > 1.0
    for n in dup.values():
        assert abs(T[n] - T[_original(n)]) <= 1e-10 * max(abs(v) for v in T.values())
