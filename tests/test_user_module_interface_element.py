"""Un elemento de interfaz de un módulo de usuario (ADR 0020 y ADR 0043 §1).

Promesa del ADR 0043: ``InterfaceElement`` es la base pública de los elementos
de espesor nulo y un módulo de usuario la hereda por la ruta pública
``solidum.elements.interface``, declarando sólo su cara y sus nodos. Con un
módulo de juguete (escrito en un directorio temporal y cargado en un
subproceso, para no contaminar el registro de la sesión de pytest):

1. sus fuentes sólo usan el contrato público (el detector de
   ``tests/test_user_module_public_imports.py`` no encuentra nada);
2. ``import solidum`` no lo carga; ``load_user_module`` lo registra y la
   clasificación por punto lo trata como salto (``JUMP``);
3. con los mismos nodos da la misma rigidez que ``Interface2D4``, bit a bit;
4. pasa la batería del barrido de elementos con su propio arnés;
5. el generador lo crea con ``element:`` (por omisión elige el del programa
   principal) y la salida lo trata como interfaz: celda VTK con ``Jump_*``,
   columna ``var: jump`` y energía almacenada de la historia.
"""
import os
import subprocess
import sys
import textwrap
from pathlib import Path

from test_user_module_public_imports import importaciones_privadas

ROOT = Path(__file__).resolve().parents[1]
TESTS = Path(__file__).resolve().parent
NOMBRE = "juguete_interfaz"

FUENTE = textwrap.dedent('''\
    """Elemento de interfaz de juguete: hereda la base pública."""
    from solidum.elements.interface import InterfaceElement
    from solidum.registry import ElementRegistry


    @ElementRegistry.register
    class InterfazJuguete2D(InterfaceElement):
        DOF_NAMES = ["ux", "uy"]
        JUMP_DIM = 2
        FACE_CELL = "line"
        N_FACE_NODES = 2
        A_NODES = (0, 1)
        B_NODES = (3, 2)
        N_INTEGRATION_POINTS = 2
        VTK_CELL_TYPE = "quad"
    ''')
ARCHIVOS = {"__init__.py": '"""Módulo de usuario de juguete."""\n', "interfaz.py": FUENTE}

GUION = textwrap.dedent('''\
    import sys
    from pathlib import Path

    import numpy as np

    import solidum
    import solidum.user

    solidum.user.__path__.append(sys.argv[1])
    carpeta = Path(sys.argv[2])
    assert not any(m.startswith("solidum.user.juguete_interfaz") for m in sys.modules), \\
        "import solidum cargó el módulo de usuario"
    solidum.load_user_module("juguete_interfaz")
    from solidum.core.point_quantity import JUMP, point_quantity
    from solidum.elements.interface import Interface2D4
    from solidum.materials.elastic_2d import Elastic2D
    from solidum.materials.interface import ContactPenalty
    from solidum.math.assembly import Assembler
    from solidum.mesh import split_interfaces, structured_mesh
    from solidum.registry import ElementRegistry

    cls = ElementRegistry.get("InterfazJuguete2D")
    assert point_quantity(cls).kind == JUMP
    assert solidum.user.user_module_of(cls) == "juguete_interfaz"
    assert solidum.user.user_module_of(Interface2D4) is None

    def placa():
        dom = structured_mesh("Quad4", lines=[np.linspace(0, 4, 5), np.linspace(0, 2, 5)],
                              material=Elastic2D(E=1000.0, nu=0.25), thickness=2.0).domain
        ids = sorted(i for i, n in dom.nodes.items()
                     if abs(n.coordinates[0] - 2.0) < 1e-9 and n.coordinates[1] <= 1.0 + 1e-9)
        dom.physical_groups, dom.physical_group_cells = {"grieta": ids}, {"grieta": []}
        return dom

    ley = ContactPenalty(K_n=1.0e5, mu=0.3, K_t=1.0e5)
    rigideces = {}
    for elemento in ("InterfazJuguete2D", None):
        dom = placa()
        made = split_interfaces(dom, "grieta", ley, name="g", element=elemento)
        tipos = {type(dom.elements[i]).__name__ for i in made.elements}
        # Por omisión, el del programa principal (aunque el de juguete encaja).
        assert tipos == ({"InterfazJuguete2D"} if elemento else {"Interface2D4"}), tipos
        dom.generate_equation_numbers()
        asm = Assembler(dom)
        asm.assemble_system()
        rigideces[elemento] = asm.K_global
    assert (rigideces["InterfazJuguete2D"] != rigideces[None]).nnz == 0, "K distinta"

    modelo = carpeta / "m.yaml"
    modelo.write_text("""\\
    user_modules: [juguete_interfaz]
    materials:
      - {id: 1, type: Elastic2D, E: 1000.0, nu: 0.0}
    interface_materials:
      - {id: 1, type: ContactPenalty, K_n: 1.0e5}
    mesh_element: Quad4
    mesh:
      structured:
        blocks:
          - {name: abajo, material: 1, x: {length: 2, n: 2}, y: {length: 1, n: 1}}
          - {name: arriba, material: 1, origin: [0, 1.1], x: {length: 2, n: 2},
             y: {length: 1, n: 1}}
        groups:
          Origen: "near(x, 0) and near(y, 0)"
    interfaces:
      - {name: junta, pair: [abajo.ymax, arriba.ymin], material: 1, element: InterfazJuguete2D}
    boundary_conditions_by_group:
      abajo.ymin: {uy: 0.0}
      Origen: {ux: 0.0}
      arriba.ymax: {uy: -0.2, ux: 0.0}
    solver: {type: NonlinearSolver, num_steps: 4}
    output:
      file_name: junta
      history:
        columns:
          - {name: gn, element: {interface: junta, var: jump, component: n}}
          - {name: E, energy: elastic}
    """, encoding="utf-8")
    m = solidum.load_yaml(modelo)
    tipos = {type(m.domain.elements[i]).__name__ for i in m.domain.interfaces["junta"]}
    assert tipos == {"InterfazJuguete2D"}, tipos
    m.solve(write_output=True)
    import meshio
    vtu = meshio.read(carpeta / "junta_final.vtu")
    assert "Jump_n" in vtu.cell_data and "Traction_n" in vtu.cell_data
    filas = (carpeta / "junta_historia.csv").read_text(encoding="utf-8").splitlines()
    cab = filas[0].split(",")
    ultima = dict(zip(cab, map(float, filas[-1].split(","))))
    E_B, K_N = 1000.0, 1.0e5
    sigma = 0.1 / (2.0 / E_B + 1.0 / K_N)
    assert abs(ultima["gn"] - (-sigma / K_N)) < 1e-9 * 0.1, ultima
    assert ultima["E"] > 0.0
    print("OK")
    ''')

BATERIA = textwrap.dedent('''\
    import sys

    sys.path.insert(0, {tests!r})
    import solidum
    import solidum.user

    solidum.user.__path__.append({modulos!r})
    solidum.load_user_module("juguete_interfaz")

    import test_element_contract_sweep as barrido
    from solidum.registry import ElementRegistry

    ArnesJuguete = barrido.arnes(
        lambda: barrido.nombres_de_paquete(ElementRegistry, "solidum.user.juguete_interfaz"),
        {{"InterfazJuguete2D": (lambda c: c(1, barrido.nodos(barrido.INTERFAZ_2D),
                                            barrido.MAT_INTERFAZ, thickness=barrido.ESPESOR),
                                 None)}})


    class TestJugueteRegistro(ArnesJuguete, barrido.TestRegistroCompleto):
        pass


    class TestJugueteContratoDeclarativo(ArnesJuguete, barrido.TestContratoDeclarativo):
        pass


    class TestJugueteRigidezElemental(ArnesJuguete, barrido.TestRigidezElemental):
        pass


    class TestJugueteMasa(ArnesJuguete, barrido.TestMasa):
        pass


    class TestJugueteValidacionDefensiva(ArnesJuguete, barrido.TestValidacionDefensiva):
        pass
    ''')


_ENTORNO = {**os.environ, "PYTHONIOENCODING": "utf-8"}


def _modulo(tmp_path) -> Path:
    raiz = tmp_path / "modulos"
    paquete = raiz / NOMBRE
    paquete.mkdir(parents=True)
    for archivo, texto in ARCHIVOS.items():
        (paquete / archivo).write_text(texto, encoding="utf-8")
    return raiz


def test_el_modulo_de_juguete_solo_usa_el_contrato_publico():
    assert importaciones_privadas(FUENTE, f"solidum.user.{NOMBRE}") == []


def test_hereda_la_base_y_el_generador_y_la_salida_lo_tratan_como_interfaz(tmp_path):
    raiz = _modulo(tmp_path)
    guion = tmp_path / "prueba.py"
    guion.write_text(GUION, encoding="utf-8")
    r = subprocess.run([sys.executable, str(guion), str(raiz), str(tmp_path)], cwd=ROOT,
                       capture_output=True, text=True, encoding="utf-8", errors="replace",
                       env=_ENTORNO, timeout=300)
    assert r.returncode == 0, r.stdout + r.stderr
    assert r.stdout.strip().endswith("OK")


def test_pasa_la_bateria_del_barrido_de_elementos(tmp_path):
    raiz = _modulo(tmp_path)
    archivo = tmp_path / "test_bateria_juguete.py"
    archivo.write_text(BATERIA.format(tests=str(TESTS), modulos=str(raiz)), encoding="utf-8")
    r = subprocess.run([sys.executable, "-m", "pytest", str(archivo), "-q", "-n", "0",
                        "-p", "no:cacheprovider", "--rootdir", str(tmp_path)],
                       cwd=tmp_path, capture_output=True, text=True, encoding="utf-8",
                       errors="replace", env=_ENTORNO, timeout=600)
    assert r.returncode == 0, r.stdout[-4000:] + r.stderr[-2000:]
    assert " passed" in r.stdout and "skipped" not in r.stdout.split("passed")[0][-30:]
