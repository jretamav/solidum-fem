"""Variables internas en la salida VTK (ADR 0036, ampliación del 2026-11-25).

Barra de tres bloques en tracción uniaxial (ν = 0, esfuerzo plano, espesor
10): los dos bloques exteriores de grieta distribuida fija
(``FixedSmearedCrack2D``) y la columna central del modelo de daño de tracción
y compresión (``TensionCompressionDamage2D``), debilitada al 95 %, que es la
que se agrieta. Se comprueba:

1. El VTK por celda lleva **todas** las variables internas escalares de los
   materiales con su nombre (antes sólo la principal, ``Internal_State``),
   ``NaN`` en las celdas cuyo material no la tiene, y no promedia los ángulos
   (``STATE_ANGLES``: la θ de la grieta fija).
2. El archivo por puntos de integración (``integration_points: true``) lleva,
   en cada punto, las variables del estado confirmado —también θ—, el
   esfuerzo (σ_xx·A = P en la columna, uniaxial) y el desplazamiento
   interpolado.
3. ``element_results`` con un nombre que no es un campo del modelo es un
   error con la lista de los que hay (antes se ignoraba en silencio); un
   ángulo, con la indicación del archivo por puntos.
"""
import csv
import textwrap

import numpy as np
import pytest

import solidum

meshio = pytest.importorskip("meshio")

E, FT, GF, T, H = 30000.0, 3.0, 0.1, 10.0, 20.0

MODELO = f"""\
materials:
  - {{id: 1, type: FixedSmearedCrack2D, E: {E}, nu: 0.0, f_t: {FT}, G_f: {GF}}}
  - {{id: 2, type: TensionCompressionDamage2D, E: {E}, nu: 0.0, f_t: {0.95 * FT}, G_f: {GF},
      f_c: 30.0, G_c: 10.0}}
mesh_element: Quad4
mesh_thickness: {T}
mesh:
  structured:
    blocks:
      - {{name: a, material: 1, x: {{length: 50, n: 5}}, y: {{length: {H}, n: 2}}}}
      - {{name: w, material: 2, origin: [50, 0], x: {{length: 10, n: 1}}, y: {{length: {H}, n: 2}}}}
      - {{name: b, material: 1, origin: [60, 0], x: {{length: 40, n: 4}}, y: {{length: {H}, n: 2}}}}
    groups:
      Origen: "near(x, 0) and near(y, 0)"
boundary_conditions_by_group:
  a.xmin: {{ux: 0.0}}
  Origen: {{uy: 0.0}}
  b.xmax: {{ux: 0.4}}
solver: {{type: NonlinearSolver, num_steps: 100}}
output:
  file_name: barra
  log_level: ERROR
  results:
    frequency: {{at: [0.5]}}
    integration_points: true
  history:
    columns:
      - {{name: P, reaction: {{group: b.xmax, dof: ux}}}}
"""


@pytest.fixture(scope="module")
def resuelto(tmp_path_factory):
    d = tmp_path_factory.mktemp("variables")
    ruta = d / "barra.yaml"
    ruta.write_text(MODELO, encoding="utf-8")
    modelo = solidum.load_yaml(ruta)
    modelo.solve(write_output=True)
    with open(d / "barra_historia.csv", newline="", encoding="utf-8") as f:
        filas = list(csv.reader(f))
    P = float(filas[-1][filas[0].index("P")])
    return d, modelo, P


def _cell(m, name):
    return np.concatenate([np.asarray(a, dtype=float) for a in m.cell_data[name]])


def test_el_vtk_por_celda_lleva_todas_las_variables_internas(resuelto):
    d, modelo, _ = resuelto
    m = meshio.read(d / "barra_final.vtu")
    for name in ("damage_t", "damage_c", "kappa_t", "cracked", "w_max_n", "damage_n", "slip"):
        assert name in m.cell_data, name
    centro = np.array([np.mean([n.coordinates[0] for n in e.nodes])
                       for e in modelo.domain.elements.values()])
    columna = (centro > 50) & (centro < 60)
    dt, kt, cr = _cell(m, "damage_t"), _cell(m, "kappa_t"), _cell(m, "cracked")
    # Cada variable sólo en las celdas de su material; NaN en las demás.
    assert np.all(np.isfinite(kt[columna])) and np.all(np.isnan(kt[~columna]))
    assert np.all(np.isfinite(cr[~columna])) and np.all(np.isnan(cr[columna]))
    # Un nombre que comparten dos materiales (d⁺ del modelo 2, daño de la
    # segunda grieta del 3) lleva en cada celda el valor de su material.
    assert np.all(np.isfinite(dt))
    assert np.all(dt[columna] > 0.99)                    # la columna se separó
    assert np.all(cr[~columna] == 0.0)                   # el resto no se agrietó
    # La principal sigue siendo Internal_State (d⁺ en la columna).
    np.testing.assert_allclose(_cell(m, "Internal_State")[columna], dt[columna])


def test_los_angulos_no_se_promedian_por_celda(resuelto):
    d, _, _ = resuelto
    m = meshio.read(d / "barra_final.vtu")
    assert "theta" not in m.cell_data


def test_el_archivo_por_puntos_lleva_el_estado_de_cada_punto(resuelto):
    d, modelo, P = resuelto
    p = meshio.read(d / "barra_puntos_final.vtu")
    pd = p.point_data
    els = list(modelo.domain.elements.values())
    assert p.points.shape[0] == 4 * len(els)             # Quad4 2×2
    assert "theta" in pd and "kappa_t" in pd
    ids, gp = pd["Element_Id"].astype(int), pd["Point"].astype(int)
    for el in els:
        for g, sv in enumerate(el.state.vars):
            k = int(np.flatnonzero((ids == el.id) & (gp == g + 1))[0])
            for name, value in sv.items():
                np.testing.assert_allclose(pd[name][k], value, rtol=0, atol=0, err_msg=name)
            others = [n for n in ("theta", "kappa_t") if n not in sv]
            assert all(np.isnan(pd[n][k]) for n in others)
    # Uniaxial: σ_xx·A = P en todos los puntos de la columna (espesor ≠ 1).
    en_columna = (p.points[:, 0] > 50) & (p.points[:, 0] < 60)
    np.testing.assert_allclose(pd["Sigma_XX"][en_columna] * H * T, P, rtol=1e-6)


def test_el_desplazamiento_por_puntos_interpola_el_nodal(resuelto):
    d, _, _ = resuelto
    p = meshio.read(d / "barra_puntos_final.vtu")
    m = meshio.read(d / "barra_final.vtu")
    abajo = np.isclose(m.points[:, 1], 0.0)
    orden = np.argsort(m.points[abajo, 0])
    x_n = m.points[abajo, 0][orden]
    u_n = np.asarray(m.point_data["Displacements"])[abajo, 0][orden]
    # ux sólo depende de x dentro de cada elemento (tracción uniaxial, ν = 0).
    np.testing.assert_allclose(p.point_data["Displacements"][:, 0],
                               np.interp(p.points[:, 0], x_n, u_n), rtol=1e-9, atol=1e-14)


def test_el_archivo_por_puntos_tiene_su_indice_de_instantes(resuelto):
    d, _, _ = resuelto
    assert (d / "barra_puntos.pvd").is_file()
    texto = (d / "barra_puntos.pvd").read_text(encoding="utf-8")
    assert "barra_puntos_00" in texto


def test_sin_pedirlo_no_se_escribe_el_archivo_por_puntos(tmp_path):
    ruta = tmp_path / "barra.yaml"
    ruta.write_text(MODELO.replace("    integration_points: true\n", "")
                    .replace("num_steps: 100", "num_steps: 5"), encoding="utf-8")
    solidum.run_yaml(ruta, write_output=True)
    assert not list(tmp_path.glob("*_puntos*"))


@pytest.mark.parametrize("campos,mensaje", [
    ('["Von_Mises", "damage_c"]', None),
    ('["Von_Mises", "dano"]', "campos inexistentes en este modelo"),
    ('["theta"]', "son ángulos"),
])
def test_element_results_valida_los_nombres(tmp_path, campos, mensaje):
    ruta = tmp_path / "barra.yaml"
    ruta.write_text(MODELO.replace("    integration_points: true\n",
                                   f"    element_results: {campos}\n")
                    .replace("num_steps: 100", "num_steps: 5"), encoding="utf-8")
    if mensaje is None:
        solidum.run_yaml(ruta, write_output=True)
        m = meshio.read(tmp_path / "barra_final.vtu")
        assert set(m.cell_data) == {"Von_Mises", "damage_c"}
    else:
        with pytest.raises(ValueError, match=mensaje):
            solidum.run_yaml(ruta, write_output=True)


def test_nodal_results_valida_los_nombres(tmp_path):
    ruta = tmp_path / "barra.yaml"
    ruta.write_text(MODELO.replace("    integration_points: true\n",
                                   "    nodal_results: [Desplazamientos]\n"), encoding="utf-8")
    with pytest.raises(ValueError, match="Displacements"):
        solidum.run_yaml(ruta, write_output=True)


def test_element_results_admite_variables_que_aun_no_existen_al_empezar(tmp_path):
    """La plasticidad crea su estado en el primer paso: ``alpha`` es un campo
    válido aunque al crear el escritor el estado todavía esté vacío."""
    ruta = tmp_path / "placa.yaml"
    ruta.write_text(textwrap.dedent("""\
        materials:
          - {id: 1, type: VonMises2D, E: 200000.0, nu: 0.3, sigma_y: 250.0}
        mesh_element: Quad4
        mesh:
          structured: {x: {length: 1, n: 2}, y: {length: 1, n: 2}}
        boundary_conditions_by_group:
          xmin: {ux: 0.0, uy: 0.0}
          xmax: {ux: 0.003}
        solver: {type: NonlinearSolver, num_steps: 2}
        output:
          file_name: placa
          results:
            element_results: [alpha]
        """), encoding="utf-8")
    solidum.run_yaml(ruta, write_output=True)
    m = meshio.read(tmp_path / "placa_final.vtu")
    assert set(m.cell_data) == {"alpha"} and _cell(m, "alpha").max() > 0.0
