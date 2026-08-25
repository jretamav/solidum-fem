"""Tests de la vía YAML para el análisis térmico (Etapa 8).

Cablea la familia térmica al parser: bloque ``thermal_materials`` paralelo a
``materials``, cargas ``thermal_loads`` (fuente volumétrica y flujo de
frontera) y despacho al pipeline ``thermal_transient``.

**Por qué existe este archivo.** Los tests de los componentes térmicos
construyen el modelo desde código Python. Que un componente funcione desde
Python no implica que sea *usable*: la vía normal del programa es el YAML, y
el catálogo lo anuncia como seleccionable. Estos tests blindan esa promesa —
si alguien añade un material térmico y olvida el bloque del parser, el
componente quedaría inaccesible sin que ningún test lo notara.

Se verifica además que **el YAML mecánico no cambia de comportamiento**: la
resolución del material contra dos familias y la suma del término térmico al
vector de cargas no deben alterar un modelo que no declare nada térmico.
"""
from __future__ import annotations

import numpy as np
import pytest

import solidum
from solidum.math.assembly import Assembler
from solidum.results import SolveResult, ThermalTransientResult
from solidum.utils.yaml_parser import YamlParser, YamlValidationError

K_ACERO = 45.0
C_ACERO = 460.0
RHO_ACERO = 7850.0
ALPHA_ACERO = K_ACERO / (RHO_ACERO * C_ACERO)


@pytest.fixture
def escribir(tmp_path):
    """Escribe un YAML en el tmp_path del test y devuelve su ruta."""
    def _escribir(texto: str, nombre: str = "modelo.yaml") -> str:
        ruta = tmp_path / nombre
        ruta.write_text(texto, encoding="utf-8")
        return str(ruta)
    return _escribir


# ----------------------------------------------------------------------
# Modelos de referencia
# ----------------------------------------------------------------------

def _pared_plana(nx: int = 2, L: float = 2.0, H: float = 0.5,
                 thickness: float = 0.3, con_capacidad: bool = False,
                 solver_block: str = "solver:\n  type: LinearSolver") -> str:
    """Fila de `nx` Quad4Thermal con Dirichlet 100 → 20 en los extremos."""
    nodos, elems, bcs = [], [], []
    nid = {}
    c = 1
    for j in range(2):
        for i in range(nx + 1):
            nodos.append(f"  - {{id: {c}, coords: [{L*i/nx}, {H*j}]}}")
            nid[(i, j)] = c
            c += 1
    for i in range(nx):
        ns = [nid[(i, 0)], nid[(i + 1, 0)], nid[(i + 1, 1)], nid[(i, 1)]]
        elems.append(
            f"  - {{id: {i+1}, type: Quad4Thermal, nodes: {ns}, "
            f"material: 1, thickness: {thickness}}}"
        )
    for j in range(2):
        bcs.append(f"  - {{node_id: {nid[(0, j)]}, T: 100.0}}")
        bcs.append(f"  - {{node_id: {nid[(nx, j)]}, T: 20.0}}")

    props = f"k: {K_ACERO}"
    if con_capacidad:
        props += f", c: {C_ACERO}, density: {RHO_ACERO}"

    return (
        "nodes:\n" + "\n".join(nodos) + "\n"
        "thermal_materials:\n"
        f"  - {{id: 1, type: ThermalConduction, {props}}}\n"
        "elements:\n" + "\n".join(elems) + "\n"
        "boundary_conditions:\n" + "\n".join(bcs) + "\n"
        + solver_block + "\n"
    )


def _dofs_fila(dom, nx: int):
    """DOFs de temperatura de la fila inferior, de izquierda a derecha."""
    # Los nodos se numeraron 1..nx+1 en la fila j=0.
    return [dom.nodes[i + 1].dofs["T"] for i in range(nx + 1)]


# ----------------------------------------------------------------------
# Estacionario
# ----------------------------------------------------------------------

class TestEstacionarioPorYaml:
    """El régimen permanente lo resuelve el `LinearSolver` sin modificación."""

    def test_pared_plana_reproduce_el_perfil_lineal_analitico(self, escribir):
        nx = 6
        ruta = escribir(_pared_plana(nx=nx))
        r = solidum.run_yaml(ruta)

        assert isinstance(r, SolveResult)
        # Reparsear para recuperar la numeración de DOFs y leer el perfil.
        parser = YamlParser(ruta)
        dom = parser.parse()
        dom.generate_equation_numbers()

        perfil = r.U[_dofs_fila(dom, nx)]
        exacto = 100.0 + (20.0 - 100.0) * np.linspace(0.0, 1.0, nx + 1)
        assert np.allclose(perfil, exacto, rtol=0, atol=1e-9)

    def test_el_material_termico_se_construye_desde_su_bloque(self, escribir):
        from solidum.materials.thermal_conduction import ThermalConduction

        ruta = escribir(_pared_plana())
        parser = YamlParser(ruta)
        parser.parse()
        mat = parser.thermal_materials[1]
        assert isinstance(mat, ThermalConduction)
        assert mat.FLUX_DIM == 2
        assert np.allclose(mat.conductivity, K_ACERO * np.eye(2))

    def test_el_estacionario_no_exige_capacidad_calorifica(self, escribir):
        """`ρ` y `c` sobran aquí: el mensaje del material ya lo decía."""
        ruta = escribir(_pared_plana(con_capacidad=False))
        r = solidum.run_yaml(ruta)
        assert r.converged

    def test_dirichlet_por_nombre_de_dof_no_necesito_caso_especial(self, escribir):
        """`dof: T` funciona por el mecanismo genérico ya existente.

        Las BCs del parser se aplican por *nombre* de DOF, no por una lista
        cerrada de nombres mecánicos. Fue lo único de la cadena térmica que no
        hubo que tocar, y merece quedar blindado.
        """
        nx = 4
        ruta = escribir(_pared_plana(nx=nx))
        parser = YamlParser(ruta)
        dom = parser.parse()
        assert dom.nodes[1].boundary_conditions.get("T") == 100.0
        assert dom.nodes[nx + 1].boundary_conditions.get("T") == 20.0


# ----------------------------------------------------------------------
# Transitorio
# ----------------------------------------------------------------------

class TestTransitorioPorYaml:

    @staticmethod
    def _bloque_solver(dt: float, n_steps: int, extra: str = "") -> str:
        return (
            "solver:\n"
            "  type: ThetaMethodSolver\n"
            f"  dt: {dt}\n"
            f"  n_steps: {n_steps}\n"
            "  T_initial: 20.0\n" + extra
        )

    def test_despacha_al_pipeline_termico_y_devuelve_su_resultado(self, escribir):
        ruta = escribir(_pared_plana(
            con_capacidad=True,
            solver_block=self._bloque_solver(2000.0, 50),
        ))
        r = solidum.run_yaml(ruta)
        assert isinstance(r, ThermalTransientResult)
        assert r.n_steps == 50
        assert r.theta == 1.0        # default L-estable
        assert r.order == 1

    def test_converge_al_mismo_estacionario_que_el_linear_solver(self, escribir):
        """El transitorio integrado hasta `t ≫ L²/α` debe igualar al estático.

        Es el cross-check de las dos vías YAML entre sí: si el parser cableara
        mal las cargas o las BCs en una de las dos, divergirían.
        """
        nx, L = 6, 2.0
        t_dif = L * L / ALPHA_ACERO
        n = 300
        estatico = solidum.run_yaml(escribir(_pared_plana(nx=nx, L=L),
                                              "estatico.yaml"))
        transitorio = solidum.run_yaml(escribir(
            _pared_plana(nx=nx, L=L, con_capacidad=True,
                         solver_block=self._bloque_solver(5 * t_dif / n, n)),
            "transitorio.yaml",
        ))
        assert np.allclose(transitorio.T_final, estatico.U, rtol=0, atol=1e-9)

    def test_theta_configurable_desde_yaml(self, escribir):
        ruta = escribir(_pared_plana(
            con_capacidad=True,
            solver_block=self._bloque_solver(2000.0, 10, "  theta: 0.5\n"),
        ))
        r = solidum.run_yaml(ruta)
        assert r.theta == 0.5
        assert r.order == 2      # Crank-Nicolson gana un orden

    def test_output_every_configurable_desde_yaml(self, escribir):
        ruta = escribir(_pared_plana(
            con_capacidad=True,
            solver_block=self._bloque_solver(1000.0, 20, "  output_every: 5\n"),
        ))
        r = solidum.run_yaml(ruta)
        assert r.n_steps == 20               # pasos integrados
        assert len(r.t_history) < 21         # columnas almacenadas

    def test_el_transitorio_si_exige_capacidad_calorifica(self, escribir):
        """Sin `ρc` debe fallar con el mensaje térmico, no con el mecánico."""
        ruta = escribir(_pared_plana(
            con_capacidad=False,
            solver_block=self._bloque_solver(1000.0, 5),
        ))
        with pytest.raises(ValueError, match="LinearSolver"):
            solidum.run_yaml(ruta)


# ----------------------------------------------------------------------
# Cargas térmicas
# ----------------------------------------------------------------------

class TestCargasTermicas:

    @staticmethod
    def _un_elemento(thermal_loads: str, thickness: float = 2.0,
                      solver_block: str = "solver:\n  type: LinearSolver") -> str:
        return (
            "nodes:\n"
            "  - {id: 1, coords: [0.0, 0.0]}\n"
            "  - {id: 2, coords: [1.0, 0.0]}\n"
            "  - {id: 3, coords: [1.0, 1.0]}\n"
            "  - {id: 4, coords: [0.0, 1.0]}\n"
            "thermal_materials:\n"
            f"  - {{id: 1, type: ThermalConduction, k: {K_ACERO}}}\n"
            "elements:\n"
            f"  - {{id: 1, type: Quad4Thermal, nodes: [1,2,3,4], "
            f"material: 1, thickness: {thickness}}}\n"
            "boundary_conditions:\n"
            "  - {node_id: 1, T: 0.0}\n"
            "  - {node_id: 2, T: 0.0}\n"
            + thermal_loads + solver_block + "\n"
        )

    def test_fuente_volumetrica_integra_Q_por_volumen(self, escribir):
        Q, thickness = 1000.0, 2.0
        ruta = escribir(self._un_elemento(
            "thermal_loads:\n  body_source:\n    - {Q: %g}\n" % Q,
            thickness=thickness,
        ))
        parser = YamlParser(ruta)
        dom = parser.parse()
        dom.generate_equation_numbers()
        F = parser.get_thermal_loads()
        # Volumen = 1 × 1 × thickness.
        assert F.sum() == pytest.approx(Q * 1.0 * 1.0 * thickness, rel=1e-12)

    def test_flujo_saliente_extrae_energia(self, escribir):
        """`q̄ > 0` ⇒ vector NEGATIVO (Reglas.md §5).

        El signo sale de la forma débil, donde el término de frontera aparece
        como `−∫ w·q̄ dΓ`. Es la convención que más fácilmente se invierte por
        descuido, y la que un test debe fijar.
        """
        q_bar, thickness = 50.0, 2.0
        ruta = escribir(self._un_elemento(
            "thermal_loads:\n  boundary_flux:\n"
            "    - {element: 1, edge: 2, q: %g}\n" % q_bar,
            thickness=thickness,
        ))
        parser = YamlParser(ruta)
        dom = parser.parse()
        dom.generate_equation_numbers()
        F = parser.get_thermal_loads()
        # Borde 2 = nodos 3-4, longitud 1.
        assert F.sum() == pytest.approx(-q_bar * 1.0 * thickness, rel=1e-12)
        assert np.all(F <= 0.0)

    def test_flujo_entrante_se_declara_con_signo_negativo(self, escribir):
        ruta = escribir(self._un_elemento(
            "thermal_loads:\n  boundary_flux:\n"
            "    - {element: 1, edge: 2, q: -50.0}\n",
        ))
        parser = YamlParser(ruta)
        dom = parser.parse()
        dom.generate_equation_numbers()
        F = parser.get_thermal_loads()
        assert F.sum() > 0.0

    def test_fuente_y_flujo_combinados_contra_solucion_analitica(self, escribir):
        """Pared 1D con generación interna y enfriamiento en la cara opuesta.

        Balance: `d/dy(k·dT/dy) + Q = 0` con `T(0) = 0` y `−k·T'(L) = q̄`.
        Integrando: `T(y) = −(Q/2k)·y² + C₁·y` con `C₁ = (Q·L − q̄)/k`.
        Para `Q = 1000`, `k = 45`, `L = 1`, `q̄ = 50` ⇒ `T(L) = 10` exacto.
        """
        ruta = escribir(self._un_elemento(
            "thermal_loads:\n"
            "  body_source:\n    - {Q: 1000.0}\n"
            "  boundary_flux:\n    - {element: 1, edge: 2, q: 50.0}\n",
        ))
        r = solidum.run_yaml(ruta)
        parser = YamlParser(ruta)
        dom = parser.parse()
        dom.generate_equation_numbers()
        T_arriba = r.U[dom.nodes[3].dofs["T"]]
        assert T_arriba == pytest.approx(10.0, rel=1e-9)

    def test_body_source_acotado_a_elementos_concretos(self, escribir):
        """Sin `elements:` la fuente va a todo el dominio; con él, sólo a esos."""
        nx = 4
        base = _pared_plana(nx=nx)
        con_todo = base.replace(
            "solver:", "thermal_loads:\n  body_source:\n    - {Q: 500.0}\nsolver:")
        con_uno = base.replace(
            "solver:",
            "thermal_loads:\n  body_source:\n"
            "    - {Q: 500.0, elements: [1]}\nsolver:")

        def _suma(texto, nombre):
            ruta = escribir(texto, nombre)
            parser = YamlParser(ruta)
            dom = parser.parse()
            dom.generate_equation_numbers()
            return parser.get_thermal_loads().sum()

        total = _suma(con_todo, "todo.yaml")
        parcial = _suma(con_uno, "uno.yaml")
        assert parcial == pytest.approx(total / nx, rel=1e-12)

    def test_forma_escalar_equivale_a_lista_de_uno(self, escribir):
        """`body_source: {Q: ...}` y `body_source: [{Q: ...}]` son lo mismo."""
        def _suma(bloque, nombre):
            ruta = escribir(self._un_elemento(bloque), nombre)
            parser = YamlParser(ruta)
            dom = parser.parse()
            dom.generate_equation_numbers()
            return parser.get_thermal_loads().sum()

        escalar = _suma("thermal_loads:\n  body_source: {Q: 700.0}\n", "a.yaml")
        lista = _suma("thermal_loads:\n  body_source:\n    - {Q: 700.0}\n", "b.yaml")
        assert escalar == pytest.approx(lista, rel=1e-12)

    def test_las_cargas_termicas_llegan_al_transitorio(self, escribir):
        """El solver transitorio recibe `F_func` derivado del bloque.

        `F_func` es un callable y el YAML no puede expresarlo; el parser
        envuelve el vector constante. Sin este cableado el transitorio con
        fuente daría un resultado silenciosamente distinto al estacionario.
        """
        cargas = ("thermal_loads:\n"
                  "  body_source:\n    - {Q: 1000.0}\n"
                  "  boundary_flux:\n    - {element: 1, edge: 2, q: 50.0}\n")

        estatico = solidum.run_yaml(escribir(
            self._un_elemento(cargas), "est.yaml"))

        # El mismo modelo con ρc declarada e integrado hasta t ≫ L²/α, que con
        # L = 1 m y el acero vale ≈ 8.0e4 s. Por debajo de esa escala el campo
        # aún evoluciona y la comparación no significaría nada.
        t_difusion = 1.0 / ALPHA_ACERO
        n_steps = 200
        dt = 40.0 * t_difusion / n_steps
        yaml_transitorio = (
            "nodes:\n"
            "  - {id: 1, coords: [0.0, 0.0]}\n"
            "  - {id: 2, coords: [1.0, 0.0]}\n"
            "  - {id: 3, coords: [1.0, 1.0]}\n"
            "  - {id: 4, coords: [0.0, 1.0]}\n"
            "thermal_materials:\n"
            f"  - {{id: 1, type: ThermalConduction, k: {K_ACERO}, "
            f"c: {C_ACERO}, density: {RHO_ACERO}}}\n"
            "elements:\n"
            "  - {id: 1, type: Quad4Thermal, nodes: [1,2,3,4], "
            "material: 1, thickness: 2.0}\n"
            "boundary_conditions:\n"
            "  - {node_id: 1, T: 0.0}\n"
            "  - {node_id: 2, T: 0.0}\n"
            + cargas +
            "solver:\n"
            "  type: ThetaMethodSolver\n"
            f"  dt: {dt}\n"
            f"  n_steps: {n_steps}\n"
            "  T_initial: 0.0\n"
        )
        transitorio = solidum.run_yaml(escribir(yaml_transitorio, "tra.yaml"))

        assert np.allclose(transitorio.T_final, estatico.U, rtol=0, atol=1e-6)

    def test_flujo_sin_edge_ni_face_es_rechazado(self, escribir):
        ruta = escribir(self._un_elemento(
            "thermal_loads:\n  boundary_flux:\n    - {element: 1, q: 50.0}\n",
        ))
        parser = YamlParser(ruta)
        dom = parser.parse()
        dom.generate_equation_numbers()
        with pytest.raises(ValueError, match="edge.*face|face.*edge"):
            parser.get_thermal_loads()


# ----------------------------------------------------------------------
# Validación y mensajes
# ----------------------------------------------------------------------

class TestValidacion:

    def test_tipo_termico_desconocido_lista_los_disponibles(self, escribir):
        ruta = escribir(_pared_plana().replace(
            "type: ThermalConduction", "type: NoExiste"))
        with pytest.raises(YamlValidationError, match="ThermalConduction"):
            solidum.run_yaml(ruta)

    def test_el_error_distingue_la_familia_termica_de_la_mecanica(self, escribir):
        """No debe ofrecer materiales mecánicos como alternativa a uno térmico."""
        ruta = escribir(_pared_plana().replace(
            "type: ThermalConduction", "type: NoExiste"))
        with pytest.raises(YamlValidationError) as exc:
            solidum.run_yaml(ruta)
        assert "material térmico desconocido" in str(exc.value)
        assert "Elastic2D" not in str(exc.value)

    def test_referencia_a_material_inexistente_nombra_los_dos_bloques(self, escribir):
        ruta = escribir(_pared_plana().replace("material: 1", "material: 99"))
        with pytest.raises(YamlValidationError) as exc:
            solidum.run_yaml(ruta)
        assert "thermal_materials" in str(exc.value)

    def test_modelo_termico_sin_bloque_materials_es_valido(self, escribir):
        """`thermal_materials` sola basta: no hay que declarar `materials: []`."""
        ruta = escribir(_pared_plana())
        assert "materials:\n" not in ruta      # sanity del fixture
        r = solidum.run_yaml(ruta)
        assert r.converged

    def test_modelo_sin_ningun_material_sigue_rechazandose(self, escribir):
        ruta = escribir(_pared_plana().replace("thermal_materials:", "otros:"))
        with pytest.raises(YamlValidationError, match="thermal_materials"):
            solidum.run_yaml(ruta)

    def test_parametro_no_aceptado_por_el_elemento_es_rechazado(self, escribir):
        """La validación de kwargs por firma vale igual para los térmicos."""
        ruta = escribir(_pared_plana().replace(
            "thickness: 0.3", "thickness: 0.3, inventado: 1.0"))
        with pytest.raises(YamlValidationError, match="inventado"):
            solidum.run_yaml(ruta)


# ----------------------------------------------------------------------
# Cuadratura
# ----------------------------------------------------------------------

class TestCuadratura:
    """Los térmicos reciben la **clave**, los mecánicos la regla materializada.

    El elemento térmico resuelve la cuadratura en su propio constructor para
    poder guardar `quadrature_key` y reportarla en diagnósticos. Pasarle el
    objeto ya materializado, como se hace con los mecánicos, rompería esa
    invariante — de ahí la bifurcación en el parser, que estos tests fijan.
    """

    def test_clave_de_cuadratura_se_propaga_al_elemento(self, escribir):
        ruta = escribir(_pared_plana().replace(
            "thickness: 0.3", "thickness: 0.3, quadrature: '3x3'"))
        parser = YamlParser(ruta)
        dom = parser.parse()
        elem = dom.elements[1]
        assert elem.quadrature_key == "3x3"
        assert elem.N_INTEGRATION_POINTS == 9

    def test_cuadratura_por_defecto_si_no_se_declara(self, escribir):
        ruta = escribir(_pared_plana())
        parser = YamlParser(ruta)
        dom = parser.parse()
        assert dom.elements[1].quadrature_key == "2x2"

    def test_el_mecanico_sigue_recibiendo_la_regla_materializada(self, escribir):
        """Blindaje de no-regresión de la bifurcación introducida."""
        ruta = escribir(
            "nodes:\n"
            "  - {id: 1, coords: [0.0, 0.0]}\n"
            "  - {id: 2, coords: [1.0, 0.0]}\n"
            "  - {id: 3, coords: [1.0, 1.0]}\n"
            "  - {id: 4, coords: [0.0, 1.0]}\n"
            "materials:\n"
            "  - {id: 1, type: Elastic2D, E: 210.0e9, nu: 0.3}\n"
            "elements:\n"
            "  - {id: 1, type: Quad4, nodes: [1,2,3,4], material: 1, "
            "quadrature: '3x3'}\n"
            "boundary_conditions:\n"
            "  - {node_id: 1, ux: 0.0, uy: 0.0}\n"
            "  - {node_id: 2, uy: 0.0}\n"
            "solver:\n  type: LinearSolver\n"
        )
        parser = YamlParser(ruta)
        dom = parser.parse()
        assert dom.elements[1].N_INTEGRATION_POINTS == 9


# ----------------------------------------------------------------------
# No regresión del camino mecánico
# ----------------------------------------------------------------------

class TestNoRegresionMecanica:
    """El cableado térmico no debe alterar un modelo puramente mecánico."""

    @staticmethod
    def _voladizo() -> str:
        # Las cargas se nombran por el DOF sobre el que actúan (`ux`), no por
        # una componente de fuerza (`fx`): el parser las aplica por nombre de
        # DOF, el mismo mecanismo genérico que permite `T:` en el térmico.
        return (
            "nodes:\n"
            "  - {id: 1, coords: [0.0, 0.0]}\n"
            "  - {id: 2, coords: [1.0, 0.0]}\n"
            "  - {id: 3, coords: [1.0, 1.0]}\n"
            "  - {id: 4, coords: [0.0, 1.0]}\n"
            "materials:\n"
            "  - {id: 1, type: Elastic2D, E: 210.0e9, nu: 0.3}\n"
            "elements:\n"
            "  - {id: 1, type: Quad4, nodes: [1,2,3,4], material: 1}\n"
            "boundary_conditions:\n"
            "  - {node_id: 1, ux: 0.0, uy: 0.0}\n"
            "  - {node_id: 4, ux: 0.0, uy: 0.0}\n"
            "point_loads:\n"
            "  - {node_id: 2, ux: 1000.0}\n"
            "solver:\n  type: LinearSolver\n"
        )

    def test_modelo_mecanico_sin_bloques_termicos_no_cambia(self, escribir):
        r = solidum.run_yaml(escribir(self._voladizo()))
        assert r.converged
        assert np.any(r.U != 0.0)

    def test_get_thermal_loads_es_cero_sin_bloque(self, escribir):
        """El término que se suma al vector de cargas debe ser inocuo."""
        ruta = escribir(self._voladizo())
        parser = YamlParser(ruta)
        dom = parser.parse()
        dom.generate_equation_numbers()
        F = parser.get_thermal_loads()
        assert F.shape == (dom.total_dofs,)
        assert not np.any(F)
