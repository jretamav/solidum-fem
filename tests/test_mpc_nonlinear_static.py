"""Los solvers estáticos no lineales convergen con restricciones multipunto.

Con una restricción afín ``u_s = g + Σ α_i·u_mi`` (ADR 0004) el equilibrio
reducido es ``Tᵀ·R = 0``: la fuerza del DOF esclavo se transfiere a sus
maestros. Hasta el 2026-09-30 ``NonlinearSolver`` y ``ArcLengthSolver`` (y
sus variantes, que heredan ``ArcProblem``) medían la convergencia con
``‖R[libres]‖`` —el residuo en los DOF no esclavos—, que en la solución exacta
no se anula en los maestros de una MPC: el Newton llegaba a la solución
(corrección ~1e-25) pero su medida del residuo nunca bajaba y el solver
acababa en ``OscillatingNewtonError``. Con sólo condiciones de Dirichlet
``T`` es la identidad sobre los libres y ambas medidas coinciden, por eso no
se veía. Lo destapó la batería de modelos válidos de la revisión de la red de
seguridad del ADR 0019. ``NewtonNewmarkSolver`` ya medía ``‖TᵀR‖``.

Dos familias de comprobaciones:

- Con material elástico lineal, cada solver reproduce la solución de
  ``LinearSolver`` en λ = 1 (apoyo oblicuo y placas atadas).
- **Oráculo del nodo fundido** (auditoría de las MPC del 2026-09-30): atar dos
  nodos coincidentes con ``u_s = u_m`` es cinemáticamente idéntico a fundirlos
  en uno, así que el modelo atado y el de malla conforme dan el mismo campo de
  desplazamientos hasta el redondeo, también con plasticidad, con cadenas de
  restricciones, con el maestro apoyado y con no linealidad geométrica. Con un
  término independiente ``g`` la parte atada se traslada rígidamente ``λ·g``:
  ``g``, como un desplazamiento prescrito, se aplica proporcionalmente al
  factor de carga.
"""
import numpy as np
import pytest

from solidum.core.domain import Domain
from solidum.elements.solid_2d.quad4 import Quad4
from solidum.elements.truss import Truss2D, Truss2DCorot
from solidum.materials.elastic import Elastic1D
from solidum.materials.elastic_2d import Elastic2D
from solidum.materials.von_mises_2d import VonMises2D
from solidum.math.assembly import Assembler
from solidum.math.convergence import ConvergenceCriterion
from solidum.math.solvers import (
    ArcLengthSolver,
    IndirectDisplacementSolver,
    LinearSolver,
    NonlinearSolver,
)

E, NU, T_ESP = 30.0e9, 0.2, 0.2
NX, NY, LX, LY = 8, 4, 2.0, 1.0


def _placa(dom, x0, nid0, eid0, mat):
    """Malla Quad4 de NX × NY sobre [x0, x0 + LX] × [0, LY]; devuelve el mapa
    (i, j) → nodo."""
    nodos = {}
    for j in range(NY + 1):
        for i in range(NX + 1):
            nid = nid0 + j * (NX + 1) + i
            nodos[i, j] = dom.add_node(nid, [x0 + LX * i / NX, LY * j / NY])
    eid = eid0
    for j in range(NY):
        for i in range(NX):
            conn = [nodos[i, j], nodos[i + 1, j], nodos[i + 1, j + 1], nodos[i, j + 1]]
            dom.add_element(Quad4(eid, conn, mat, thickness=T_ESP))
            eid += 1
    return nodos


def _apoyo_oblicuo():
    """Placa apoyada en una esquina y con un rodillo oblicuo a 30° en la otra,
    impuesto con una MPC (u_y = tan30°·u_x). Carga vertical en el borde
    superior."""
    dom = Domain()
    mat = Elastic2D(E=E, nu=NU)
    n = _placa(dom, 0.0, 1, 1, mat)
    n[0, 0].fix_dof("ux", 0.0)
    n[0, 0].fix_dof("uy", 0.0)
    c = n[NX, 0]
    dom.add_linear_constraint((c.id, "uy"), [(c.id, "ux")], [np.tan(np.radians(30.0))])
    dom.generate_equation_numbers(verbose=False)
    F = np.zeros(dom.total_dofs)
    for i in range(NX + 1):
        F[n[i, NY].dofs["uy"]] = -1.0e5
    return dom, F, n[NX // 2, NY]


def _placas_atadas():
    """Dos placas no conformes, separadas 1 mm, unidas por MPC nudo a nudo
    (u_esclavo = u_maestro) en su borde común; la primera empotrada en x = 0 y
    la segunda cargada en su borde libre."""
    dom = Domain()
    mat = Elastic2D(E=E, nu=NU)
    a = _placa(dom, 0.0, 1, 1, mat)
    b = _placa(dom, LX + 1.0e-3, 1001, 1001, mat)
    for j in range(NY + 1):
        a[0, j].fix_dof("ux", 0.0)
        a[0, j].fix_dof("uy", 0.0)
        for dof in ("ux", "uy"):
            dom.add_linear_constraint((b[0, j].id, dof), [(a[NX, j].id, dof)], [1.0])
    dom.generate_equation_numbers(verbose=False)
    F = np.zeros(dom.total_dofs)
    for j in range(NY + 1):
        F[b[NX, j].dofs["ux"]] = 1.0e5
        F[b[NX, j].dofs["uy"]] = -1.0e5
    return dom, F, b[NX, NY]


MODELOS = {"apoyo_oblicuo": _apoyo_oblicuo, "placas_atadas": _placas_atadas}


def _solvers(control_node):
    return {
        "NonlinearSolver": lambda asm: NonlinearSolver(asm, num_steps=2),
        "ArcLengthSolver": lambda asm: ArcLengthSolver(asm, initial_dlambda=0.3, max_lambda=1.0),
        "IndirectDisplacementSolver": lambda asm: IndirectDisplacementSolver(
            asm, initial_dlambda=0.3, max_lambda=1.0, control=[(control_node, "uy", -1.0)]),
    }


@pytest.mark.parametrize("modelo", list(MODELOS))
@pytest.mark.parametrize("solver", ["NonlinearSolver", "ArcLengthSolver", "IndirectDisplacementSolver"])
def test_con_mpc_reproduce_la_solucion_lineal(modelo, solver):
    dom, F, ctrl = MODELOS[modelo]()
    U_lin = LinearSolver(Assembler(dom)).solve(F)
    dom2, F2, ctrl2 = MODELOS[modelo]()
    U = _solvers(ctrl2)[solver](Assembler(dom2)).solve(F2)
    np.testing.assert_allclose(U, U_lin, rtol=0.0, atol=1e-6 * np.abs(U_lin).max())


@pytest.mark.parametrize("modelo", list(MODELOS))
def test_la_mpc_se_cumple_en_la_solucion(modelo):
    """La solución del Newton cumple la restricción afín exactamente."""
    dom, F, _ = MODELOS[modelo]()
    asm = Assembler(dom)
    U = NonlinearSolver(asm, num_steps=2).solve(F)
    cs = asm.constraint_set
    T, g = cs.build(dom.total_dofs)
    np.testing.assert_allclose(T @ U[cs.free_dofs(dom.total_dofs)] + g, U,
                               rtol=0.0, atol=1e-12 * np.abs(U).max())


# --- Oráculo del nodo fundido --------------------------------------------------

NXO, NYO, LXO, LYO = 4, 2, 1.0, 0.5
P_PUNTA = 2.0e6 / 3.0  # por nodo del borde libre: plastifica con VonMises2D
DESFASE = (0.3, 0.1)   # traslación de la segunda placa que cierra g = −desfase

MATERIALES = {
    "elastico": lambda: Elastic2D(E=200.0e9, nu=0.3),
    "von_mises": lambda: VonMises2D(E=200.0e9, nu=0.3, sigma_y=250.0e6, H=2.0e9,
                                    hypothesis="plane_stress"),
}


def _estricto():
    return ConvergenceCriterion(rtol_force=1e-10, rtol_disp=1e-10)


def _placas(atadas, material, *, cadena=False, desfase=(0.0, 0.0), apoyo_en_interfaz=False):
    """Voladizo [0, 2·LXO] × [0, LYO] de dos placas Quad4, empotrado en x = 0 y
    con carga vertical en el borde libre.

    ``atadas=False``: malla conforme, las placas comparten los nodos de x = LXO.
    ``atadas=True``: la segunda placa tiene nodos propios en la interfaz,
    esclavos (``ux``, ``uy``) de los de la primera; con ``desfase`` la segunda
    placa está trasladada y la restricción lleva ``g = −desfase``.
    ``cadena`` cuelga del nodo inferior de la interfaz una barra ``Truss2D``
    que baja a un apoyo. En el modelo atado la barra tiene su propio nodo
    superior ``t``, esclavo del nodo de interfaz de la segunda placa (con
    ``g = +desfase``), que a su vez es esclavo del de la primera: la cadena
    ``t ← b ← a`` se cierra en ``t ← a`` con ``g = 0``.
    ``apoyo_en_interfaz`` apoya ``uy`` en el nodo inferior de la interfaz (el
    maestro, en el modelo atado).

    Devuelve ``(dom, F, a, b)``, con ``a`` y ``b`` los mapas (i, j) → nodo."""
    dom = Domain()
    mat = MATERIALES[material]()
    nid = [0]
    eid = [0]

    def nodo(x, y):
        nid[0] += 1
        return dom.add_node(nid[0], [x, y])

    def placa(x0, y0=0.0, izquierda=None):
        n = {}
        for j in range(NYO + 1):
            for i in range(NXO + 1):
                if i == 0 and izquierda is not None:
                    n[i, j] = izquierda[j]
                else:
                    n[i, j] = nodo(x0 + LXO * i / NXO, y0 + LYO * j / NYO)
        for j in range(NYO):
            for i in range(NXO):
                eid[0] += 1
                dom.add_element(Quad4(eid[0], [n[i, j], n[i + 1, j], n[i + 1, j + 1], n[i, j + 1]],
                                      mat, thickness=T_ESP))
        return n

    a = placa(0.0)
    if not atadas:
        b = placa(LXO, izquierda=[a[NXO, j] for j in range(NYO + 1)])
    else:
        b = placa(LXO + desfase[0], y0=desfase[1])
        for j in range(NYO + 1):
            for d, s in (("ux", desfase[0]), ("uy", desfase[1])):
                dom.add_linear_constraint((b[0, j].id, d), [(a[NXO, j].id, d)], [1.0], g=-s)
    if cadena:
        pie = nodo(LXO, -1.0)
        cabeza = nodo(LXO, 0.0) if atadas else a[NXO, 0]
        dom.add_element(Truss2D(10_000, [pie, cabeza], Elastic1D(E=200.0e9), A=1.0e-3))
        pie.fix_dof("ux", 0.0)
        pie.fix_dof("uy", 0.0)
        if atadas:
            for d, s in (("ux", desfase[0]), ("uy", desfase[1])):
                dom.add_linear_constraint((cabeza.id, d), [(b[0, 0].id, d)], [1.0], g=+s)
    for j in range(NYO + 1):
        a[0, j].fix_dof("ux", 0.0)
        a[0, j].fix_dof("uy", 0.0)
    if apoyo_en_interfaz:
        a[NXO, 0].fix_dof("uy", 0.0)
    dom.generate_equation_numbers(verbose=False)
    F = np.zeros(dom.total_dofs)
    for j in range(NYO + 1):
        F[b[NXO, j].dofs["uy"]] = -P_PUNTA
    return dom, F, a, b


def _solver_oraculo(nombre, asm, b):
    if nombre == "NonlinearSolver":
        return NonlinearSolver(asm, _estricto(), num_steps=10)
    if nombre == "ArcLengthSolver":
        return ArcLengthSolver(asm, _estricto(), initial_dlambda=0.1, max_lambda=1.0, max_steps=300)
    return IndirectDisplacementSolver(asm, _estricto(), initial_dlambda=0.1, max_lambda=1.0,
                                      max_steps=400, control=[(b[NXO, NYO], "uy", -1.0)])


def _desplazamientos(a, b, U):
    """(ux, uy) de los nodos de cada placa, en el orden de sus mapas."""
    return np.array([[U[n.dofs["ux"]], U[n.dofs["uy"]]]
                     for placa in (a, b) for _, n in sorted(placa.items())])


def _resolver(nombre, atadas, material, **kw):
    dom, F, a, b = _placas(atadas, material, **kw)
    asm = Assembler(dom)
    historia = []
    solver = _solver_oraculo(nombre, asm, b)
    U = solver.solve(F, step_callback=lambda k, U, lam: historia.append((lam, U.copy())))
    return dom, asm, a, b, U, historia


CASOS_ORACULO = [
    ("elastico", "NonlinearSolver"),
    ("elastico", "ArcLengthSolver"),
    ("elastico", "IndirectDisplacementSolver"),
    ("von_mises", "NonlinearSolver"),
    ("von_mises", "IndirectDisplacementSolver"),
]


@pytest.mark.parametrize("material,solver", CASOS_ORACULO)
def test_atadas_igual_que_fundidas(material, solver):
    """El modelo atado reproduce el de malla conforme hasta el redondeo.

    Con plasticidad se omite ``ArcLengthSolver``: su longitud de arco se mide
    con ``‖ΔU‖`` sobre todos los DOF, esclavos incluidos, así que el modelo
    atado recorre otra secuencia de pasos, y con un material dependiente de la
    trayectoria llega a otro equilibrio igual de válido (spec del solver)."""
    *_, a0, b0, U0, _ = _resolver(solver, False, material)
    *_, a1, b1, U1, _ = _resolver(solver, True, material)
    u0, u1 = _desplazamientos(a0, b0, U0), _desplazamientos(a1, b1, U1)
    np.testing.assert_allclose(u1, u0, rtol=0.0, atol=1e-10 * np.abs(u0).max())


def test_la_plasticidad_esta_activa_en_el_oraculo():
    """Control del oráculo: con VonMises2D la respuesta se aparta de la elástica,
    así que la comparación anterior ejercita de verdad la plasticidad."""
    *_, a, b, U_el, _ = _resolver("NonlinearSolver", False, "elastico")
    *_, _, _, U_vm, _ = _resolver("NonlinearSolver", False, "von_mises")
    u_el, u_vm = _desplazamientos(a, b, U_el), _desplazamientos(a, b, U_vm)
    assert np.abs(u_vm - u_el).max() > 0.5 * np.abs(u_el).max()


@pytest.mark.parametrize("material,solver", [
    ("elastico", "NonlinearSolver"),
    ("elastico", "ArcLengthSolver"),
    ("elastico", "IndirectDisplacementSolver"),
    ("von_mises", "NonlinearSolver"),
])
def test_termino_independiente_escala_con_lambda(material, solver):
    """Con ``g = −desfase`` la segunda placa se traslada rígidamente, sin
    esfuerzos: sus desplazamientos son los del modelo fundido menos el desfase,
    y en cada paso convergido ``U = T·U_libre + λ·g``."""
    *_, a0, b0, U0, _ = _resolver(solver, False, material)
    dom, asm, a1, b1, U1, historia = _resolver(solver, True, material, desfase=DESFASE)
    esperado = _desplazamientos(a0, b0, U0)
    esperado[len(a0):] -= np.array(DESFASE)
    escala = np.abs(_desplazamientos(a0, b0, U0)).max()
    np.testing.assert_allclose(_desplazamientos(a1, b1, U1), esperado, rtol=0.0, atol=1e-9 * escala)
    cs = asm.constraint_set
    T, g = cs.build(dom.total_dofs)
    libres = cs.free_dofs(dom.total_dofs)
    assert historia and historia[-1][0] == pytest.approx(1.0)
    for lam, U in historia:
        np.testing.assert_allclose(T @ U[libres] + lam * g, U, rtol=0.0, atol=1e-12 * max(escala, 1.0))


@pytest.mark.parametrize("material", list(MATERIALES))
@pytest.mark.parametrize("desfase", [(0.0, 0.0), DESFASE], ids=["g_nulo", "g_compuesto"])
def test_cadena_de_restricciones(material, desfase):
    """Barra colgada por una cadena t ← b ← a (cierre transitivo del ADR 0004):
    igual que colgarla del nodo fundido. Con desfase, los términos
    independientes de la cadena se componen (+desfase − desfase = 0) y la
    segunda placa sólo se traslada."""
    *_, a0, b0, U0, _ = _resolver("NonlinearSolver", False, material, cadena=True)
    *_, a1, b1, U1, _ = _resolver("NonlinearSolver", True, material, cadena=True, desfase=desfase)
    esperado = _desplazamientos(a0, b0, U0)
    escala = np.abs(esperado).max()
    esperado[len(a0):] -= np.array(desfase)
    np.testing.assert_allclose(_desplazamientos(a1, b1, U1), esperado, rtol=0.0, atol=1e-9 * escala)


def test_maestro_apoyado():
    """Maestro con apoyo: el cierre transitivo hace Dirichlet al esclavo, y los
    desplazamientos siguen siendo los del nodo fundido apoyado."""
    *_, a0, b0, U0, _ = _resolver("NonlinearSolver", False, "von_mises", apoyo_en_interfaz=True)
    *_, a1, b1, U1, _ = _resolver("NonlinearSolver", True, "von_mises", apoyo_en_interfaz=True)
    u0 = _desplazamientos(a0, b0, U0)
    np.testing.assert_allclose(_desplazamientos(a1, b1, U1), u0, rtol=0.0, atol=1e-10 * np.abs(u0).max())


# --- No linealidad geométrica ----------------------------------------------------

H_ARCO, SEMILUZ, EA, P_REF = 0.1, 1.0, 1.0e7, 2000.0


def _arco(atado, carga_en_esclavo=False):
    """Arco rebajado de dos barras Truss2DCorot con la clave fija en ``ux``
    (simetría). Atado: cada barra con su propio nodo de clave, el de la
    derecha esclavo del de la izquierda."""
    dom = Domain()
    mat = Elastic1D(E=EA)
    izq = dom.add_node(1, [-SEMILUZ, 0.0])
    der = dom.add_node(2, [SEMILUZ, 0.0])
    clave = dom.add_node(3, [0.0, H_ARCO])
    clave2 = dom.add_node(4, [0.0, H_ARCO]) if atado else clave
    dom.add_element(Truss2DCorot(1, [izq, clave], mat, A=1.0))
    dom.add_element(Truss2DCorot(2, [clave2, der], mat, A=1.0))
    for n in (izq, der):
        n.fix_dof("ux", 0.0)
        n.fix_dof("uy", 0.0)
    clave.fix_dof("ux", 0.0)
    if atado:
        for d in ("ux", "uy"):
            dom.add_linear_constraint((clave2.id, d), [(clave.id, d)], [1.0])
    dom.generate_equation_numbers(verbose=False)
    F = np.zeros(dom.total_dofs)
    F[(clave2 if carga_en_esclavo else clave).dofs["uy"]] = -P_REF
    return dom, F, clave


@pytest.mark.parametrize("carga_en_esclavo", [False, True])
def test_arco_corrotacional_atado_sigue_la_trayectoria_del_fundido(carga_en_esclavo):
    """Con control en la clave, el arco atado recorre paso a paso la misma
    trayectoria de equilibrio que el fundido, con la carga en el maestro o en
    el esclavo."""
    trazas = []
    for atado in (False, True):
        dom, F, clave = _arco(atado, carga_en_esclavo and atado)
        solver = IndirectDisplacementSolver(Assembler(dom), _estricto(), initial_dlambda=0.1,
                                            max_lambda=3.0, max_steps=60,
                                            control=[(clave, "uy", -1.0)])
        traza = []
        solver.solve(F, step_callback=lambda k, U, lam, t=traza, c=clave: t.append((lam, U[c.dofs["uy"]])))
        trazas.append(np.array(traza))
    fundido, atado = trazas
    assert len(fundido) == len(atado) > 10
    np.testing.assert_allclose(atado, fundido, rtol=0.0, atol=1e-12)
    # El trazado entra en la rama no lineal: el factor de carga pasa por un máximo.
    assert fundido[:, 0].max() > fundido[-1, 0] + 1e-6
