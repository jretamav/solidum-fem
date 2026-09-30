"""Los solvers dinámicos con restricciones multipunto (auditoría de las MPC, 2026-09-30).

Las restricciones afines ``u_s = g + Σ α_i·u_mi`` (ADR 0004) entran en los
ocho solvers dinámicos por el mismo operador de eliminación directa que los
apoyos: ``M_red = TᵀMT``, ``C_red = TᵀCT``, ``K_red = TᵀKT``, ``F_red = TᵀF``
y la expansión ``u = T·u_libre + g`` (``u̇``, ``ü``, modos y amplitudes
armónicas con ``T``, sin ``g``). Tres familias de comprobaciones:

- **Oráculo del nodo fundido.** Atar dos nodos coincidentes con ``u_s = u_m``
  en todos sus DOFs es cinemáticamente idéntico a fundirlos en uno: el modelo
  atado y el de malla conforme dan la misma respuesta hasta el redondeo. Se
  usa una ménsula de ``Frame2DEuler`` inclinada 30° (acopla ``ux`` y ``uy``)
  partida en su nodo central, con el esclavo a cada lado del corte, para los
  cinco transitorios (desplazamientos, velocidades, aceleraciones e historia
  de fuerzas internas), el modal (frecuencias y modos expandidos), el
  armónico y el espectral (combinación, factores de participación y masas
  efectivas).
- **Diferencias centradas con una MPC de varios maestros** (nodo colgante
  ``u_h = ½·u_a + ½·u_b``). Con masa concentrada ``TᵀMT = Σ_s m_s·T_sᵀT_s``
  deja de ser diagonal y el esquema explícito no puede invertirla: el solver
  debe rechazarla nombrando la causa, con el mismo umbral relativo a escala
  de metros y a escala de 10 µm. Hasta el 2026-09-30 el umbral era absoluto
  por debajo de 1 kg y, a microescala, el acoplamiento se descartaba en
  silencio (la masa del esclavo contaba ¼ en cada maestro en lugar de ½); a
  escala de metros el mensaje culpaba a un Frame3D oblicuo. El mismo modelo
  con Newmark concentrada conserva la cantidad de movimiento
  ``Σ m_i·ü_i = ΣF`` (oráculo exacto: la traslación rígida está en el espacio
  de ``T``).
- **Condiciones iniciales.** Los transitorios reconstruyen los esclavos como
  ``T·u0[libres] + g`` y ``T·u̇0[libres]``; un valor incoherente en un esclavo
  se descarta y el solver lo avisa en el log. Las filas Dirichlet no avisan
  (su proyección es deliberada y tiene test en ``test_dynamic_audit_fixes``).
"""
import contextlib
import logging

import numpy as np
import pytest

from solidum.core.domain import Domain
from solidum.elements.frame.euler import Frame2DEuler
from solidum.elements.solid_2d.quad4 import Quad4
from solidum.materials.elastic import Elastic1D
from solidum.materials.elastic_2d import Elastic2D
from solidum.math.assembly import Assembler
from solidum.math.solvers import (
    CentralDifferenceSolver,
    HarmonicSolver,
    HHTSolver,
    ModalSolver,
    NewmarkSolver,
    NewtonHHTSolver,
    NewtonNewmarkSolver,
    ResponseSpectrumSolver,
)

E, A, I, RHO = 2.0e11, 0.01, 1.0e-4, 7850.0
NEL, L = 6, 3.0
K_MID = NEL // 2
DUP_ID = 100
LADOS = ("derecha", "izquierda")


# ----------------------------------------------------------------------------
# Modelos
# ----------------------------------------------------------------------------

def _mensula(partida: bool, lado: str = "derecha"):
    """Ménsula de ``NEL`` Frame2DEuler inclinada 30°, empotrada en el nodo 1.

    Si ``partida``, el nodo central se duplica (nodo ``DUP_ID``, coincidente),
    los elementos a su derecha usan el duplicado y una MPC por DOF ata ambos
    nodos; ``lado`` indica cuál es el esclavo: ``"derecha"`` el duplicado,
    ``"izquierda"`` el nodo original.
    """
    dom = Domain()
    mat = Elastic1D(E=E, density=RHO)
    c, s = np.cos(np.radians(30.0)), np.sin(np.radians(30.0))
    nodos = {i: dom.add_node(i + 1, [L * i / NEL * c, L * i / NEL * s])
             for i in range(NEL + 1)}
    dup = None
    if partida:
        dup = dom.add_node(DUP_ID, list(np.asarray(nodos[K_MID].coordinates, dtype=float)))
    for e in range(NEL):
        ni, nj = nodos[e], nodos[e + 1]
        if partida and e == K_MID:
            ni = dup
        dom.add_element(Frame2DEuler(e + 1, [ni, nj], mat, A=A, I=I))
    for d in ("ux", "uy", "rz"):
        nodos[0].fix_dof(d, 0.0)
    if partida:
        for d in ("ux", "uy", "rz"):
            if lado == "derecha":
                dom.add_linear_constraint((DUP_ID, d), [(nodos[K_MID].id, d)], [1.0])
            else:
                dom.add_linear_constraint((nodos[K_MID].id, d), [(DUP_ID, d)], [1.0])
    dom.generate_equation_numbers(verbose=False)
    return dom, nodos


def _indices_fundido(dom_f, dom_s):
    """Para cada DOF del modelo fundido, su índice en el modelo atado."""
    idx = np.zeros(dom_f.total_dofs, dtype=int)
    for nid, node in dom_f.nodes.items():
        for d, i in node.dofs.items():
            idx[i] = dom_s.nodes[nid].dofs[d]
    return idx


def _indices_duplicado(dom_f, dom_s):
    """DOFs del duplicado (modelo atado) y del nodo central (modelo fundido)."""
    ds = np.array([dom_s.nodes[DUP_ID].dofs[d] for d in ("ux", "uy", "rz")])
    df = np.array([dom_f.nodes[K_MID + 1].dofs[d] for d in ("ux", "uy", "rz")])
    return ds, df


def _carga_extremo(dom, nodos):
    ndof = dom.total_dofs
    tip = nodos[NEL]

    def F(t):
        f = np.zeros(ndof)
        f[tip.dofs["uy"]] = -1.0e4 * np.sin(40.0 * t)
        f[tip.dofs["ux"]] = 3.0e3 * np.sin(40.0 * t)
        return f

    return F


def _assert_rel(a, b, rtol, msg=""):
    a, b = np.asarray(a), np.asarray(b)
    den = max(float(np.abs(b).max()), 1e-300)
    err = float(np.abs(a - b).max()) / den
    assert err <= rtol, f"{msg}: error relativo {err:.2e} > {rtol:.0e}"


# ----------------------------------------------------------------------------
# Transitorios: oráculo del nodo fundido
# ----------------------------------------------------------------------------

T_END, DT = 0.01, 1.0e-4
T_END_CD, DT_CD = 2.0e-3, 2.0e-6
RAYLEIGH = {"alpha": 2.0, "beta": 1.0e-5}

TRANSITORIOS = {
    "NewmarkSolver": lambda asm, **kw: NewmarkSolver(asm, T_END, DT, **kw),
    "HHTSolver": lambda asm, **kw: HHTSolver(asm, T_END, DT, alpha=-0.1, **kw),
    "NewtonNewmarkSolver": lambda asm, **kw: NewtonNewmarkSolver(asm, T_END, DT, **kw),
    "NewtonHHTSolver": lambda asm, **kw: NewtonHHTSolver(asm, T_END, DT, alpha=-0.1, **kw),
    "CentralDifferenceSolver": lambda asm, **kw: CentralDifferenceSolver(asm, T_END_CD, DT_CD, **kw),
}

# (lumping, rayleigh): el explícito sólo admite masa concentrada.
VARIANTES = {
    "consistente": ("consistent", None),
    "concentrada_rayleigh": ("lumped", RAYLEIGH),
}


def _transitorio(nombre, partida, lado, lumping, rayleigh, **extra):
    dom, nodos = _mensula(partida, lado)
    kw = dict(F_func=_carga_extremo(dom, nodos), rayleigh=rayleigh,
              record_internal_forces=True, **extra)
    if nombre != "CentralDifferenceSolver":
        kw["lumping"] = lumping
    return dom, TRANSITORIOS[nombre](Assembler(dom), **kw).solve()


@pytest.mark.parametrize("lado", LADOS)
@pytest.mark.parametrize("solver,variante", [
    (s, v) for s in TRANSITORIOS for v in VARIANTES
    if not (s == "CentralDifferenceSolver" and VARIANTES[v][0] != "lumped")
])
def test_transitorio_reproduce_el_nodo_fundido(solver, variante, lado):
    lumping, rayleigh = VARIANTES[variante]
    dom_f, r_f = _transitorio(solver, False, lado, lumping, rayleigh)
    dom_s, r_s = _transitorio(solver, True, lado, lumping, rayleigh)
    idx = _indices_fundido(dom_f, dom_s)
    ds, df = _indices_duplicado(dom_f, dom_s)

    _assert_rel(r_s.u_history[idx], r_f.u_history, 1e-10, "u")
    _assert_rel(r_s.udot_history[idx], r_f.udot_history, 1e-9, "u̇")
    # ü: el lado izquierdo da ~2e-9 por el redondeo de la factorización.
    _assert_rel(r_s.uddot_history[idx], r_f.uddot_history, 1e-7, "ü")
    # El nodo duplicado se mueve con el central (u, u̇, ü expandidos con T).
    _assert_rel(r_s.u_history[ds], r_f.u_history[df], 1e-10, "u duplicado")
    _assert_rel(r_s.udot_history[ds], r_f.udot_history[df], 1e-9, "u̇ duplicado")
    _assert_rel(r_s.uddot_history[ds], r_f.uddot_history[df], 1e-7, "ü duplicado")

    # Historia de fuerzas internas N, V, M elemento a elemento.
    hist_f, hist_s = r_f.element_forces_history, r_s.element_forces_history
    assert set(hist_f) == set(hist_s)
    for eid, lista_f in hist_f.items():
        assert len(hist_s[eid]) == len(lista_f)
        for c in ("N", "V", "M"):
            serie_f = np.array([ef.components[c] for ef in lista_f])
            serie_s = np.array([ef.components[c] for ef in hist_s[eid]])
            escala = max(float(np.abs(serie_f).max()), 1e-9)
            err = float(np.abs(serie_s - serie_f).max()) / escala
            assert err <= 1e-8, f"elemento {eid}, {c}: error relativo {err:.2e}"


def test_carga_en_el_nodo_esclavo():
    """Una carga aplicada en el DOF esclavo se transfiere a su maestro (TᵀF)."""
    def resolver(partida):
        dom, nodos = _mensula(partida, "derecha")
        nodo = dom.nodes[DUP_ID] if partida else nodos[K_MID]
        ndof = dom.total_dofs

        def F(t):
            f = np.zeros(ndof)
            f[nodo.dofs["uy"]] = -2.0e4 * np.sin(60.0 * t)
            return f

        return dom, NewmarkSolver(Assembler(dom), T_END, DT, F_func=F).solve()

    dom_f, r_f = resolver(False)
    dom_s, r_s = resolver(True)
    _assert_rel(r_s.u_history[_indices_fundido(dom_f, dom_s)], r_f.u_history, 1e-10, "u")


def test_mpc_con_termino_independiente_lineal_y_newton_coinciden():
    """Con ``g ≠ 0`` Newmark lineal y Newton-Newmark dan la misma historia."""
    def modelo():
        dom, nodos = _mensula(False)
        tip = nodos[NEL]
        # Rodillo oblicuo afín en el extremo: u_y = 0.3·u_x + 1 mm.
        dom.add_linear_constraint((tip.id, "uy"), [(tip.id, "ux")], [0.3], g=1.0e-3)
        dom.generate_equation_numbers(verbose=False)
        return dom, nodos

    dom_a, nodos_a = modelo()
    r_a = NewmarkSolver(Assembler(dom_a), T_END, DT, F_func=_carga_extremo(dom_a, nodos_a)).solve()
    dom_b, nodos_b = modelo()
    r_b = NewtonNewmarkSolver(Assembler(dom_b), T_END, DT, F_func=_carga_extremo(dom_b, nodos_b)).solve()
    _assert_rel(r_b.u_history, r_a.u_history, 1e-9, "u")
    tip = nodos_a[NEL]
    np.testing.assert_allclose(
        r_a.u_history[tip.dofs["uy"]], 0.3 * r_a.u_history[tip.dofs["ux"]] + 1.0e-3,
        rtol=0.0, atol=1e-15)


# ----------------------------------------------------------------------------
# Modal, armónico y espectral: oráculo del nodo fundido
# ----------------------------------------------------------------------------

N_MODOS = 6


@pytest.mark.parametrize("lado", LADOS)
@pytest.mark.parametrize("lumping", ["consistent", "lumped"])
def test_modal_reproduce_el_nodo_fundido(lumping, lado):
    dom_f, _ = _mensula(False)
    r_f = ModalSolver(Assembler(dom_f), n_modes=N_MODOS, lumping=lumping).solve()
    dom_s, _ = _mensula(True, lado)
    asm_s = Assembler(dom_s)
    r_s = ModalSolver(asm_s, n_modes=N_MODOS, lumping=lumping).solve()

    _assert_rel(r_s.frequencies_rad, r_f.frequencies_rad, 1e-10, "ω")
    idx = _indices_fundido(dom_f, dom_s)
    ds, df = _indices_duplicado(dom_f, dom_s)
    for j in range(N_MODOS):
        a, b = r_s.modes[idx, j], r_f.modes[:, j]
        signo = np.sign(a @ b)
        _assert_rel(signo * a, b, 1e-8, f"modo {j + 1}")
        _assert_rel(signo * r_s.modes[ds, j], r_f.modes[df, j], 1e-8, f"modo {j + 1}, duplicado")
    # Los modos expandidos con T son M-ortonormales con la masa global del
    # modelo atado (ΦᵀMΦ = φ_redᵀ·TᵀMT·φ_red = I).
    M = asm_s.assemble_mass_matrix(lumping=lumping)
    np.testing.assert_allclose(r_s.modes.T @ (M @ r_s.modes), np.eye(N_MODOS), atol=1e-9)


@pytest.mark.parametrize("lado", LADOS)
@pytest.mark.parametrize("rayleigh", [None, RAYLEIGH], ids=["sin_rayleigh", "rayleigh"])
@pytest.mark.parametrize("lumping", ["consistent", "lumped"])
def test_armonico_reproduce_el_nodo_fundido(lumping, rayleigh, lado):
    def resolver(partida):
        dom, nodos = _mensula(partida, lado)
        F = np.zeros(dom.total_dofs, dtype=complex)
        F[nodos[NEL].dofs["uy"]] = -1.0e4
        F[nodos[NEL].dofs["ux"]] = 3.0e3j
        r = HarmonicSolver(Assembler(dom), omega_min=5.0, omega_max=2000.0, n_omega=25,
                           F_amplitude=F, rayleigh=rayleigh, lumping=lumping).solve()
        return dom, r

    dom_f, r_f = resolver(False)
    dom_s, r_s = resolver(True)
    ds, df = _indices_duplicado(dom_f, dom_s)
    _assert_rel(r_s.u_complex[_indices_fundido(dom_f, dom_s)], r_f.u_complex, 1e-10, "û")
    _assert_rel(r_s.u_complex[ds], r_f.u_complex[df], 1e-10, "û duplicado")


@pytest.mark.parametrize("lado", LADOS)
@pytest.mark.parametrize("direccion", ["ux", "uy"])
@pytest.mark.parametrize("combinacion", ["SRSS", "CQC"])
def test_espectral_reproduce_el_nodo_fundido(combinacion, direccion, lado):
    def resolver(partida):
        dom, _ = _mensula(partida, lado)
        r = ResponseSpectrumSolver(
            Assembler(dom), n_modes=N_MODOS, direction={"dof_name": direccion},
            spectrum={"type": "constant_acceleration", "Sa": 3.0},
            combination=combinacion).solve()
        return dom, r

    dom_f, r_f = resolver(False)
    dom_s, r_s = resolver(True)
    idx = _indices_fundido(dom_f, dom_s)
    ds, df = _indices_duplicado(dom_f, dom_s)
    _assert_rel(r_s.u_combined[idx], r_f.u_combined, 1e-9, "u combinada")
    _assert_rel(r_s.u_combined[ds], r_f.u_combined[df], 1e-9, "u combinada, duplicado")
    # Γ_n·φ_n no depende del signo del modo.
    _assert_rel(r_s.u_per_mode[idx], r_f.u_per_mode, 1e-9, "u por modo")
    # Γ_n = φ_nᵀ·M·r con φ_n = T·φ_red: la masa del duplicado entra una vez.
    _assert_rel(np.abs(r_s.participation_factors), np.abs(r_f.participation_factors), 1e-9, "|Γ|")
    _assert_rel(r_s.effective_masses, r_f.effective_masses, 1e-9, "masa efectiva")
    assert r_s.total_mass == pytest.approx(r_f.total_mass, rel=1e-12)
    np.testing.assert_allclose(r_s.cumulative_effective_mass_ratio(),
                               r_f.cumulative_effective_mass_ratio(), rtol=1e-9)


# ----------------------------------------------------------------------------
# Diferencias centradas y MPC con varios maestros
# ----------------------------------------------------------------------------

def _nodo_colgante(escala: float, colgante: bool = True, apoyos: bool = False):
    """Refinamiento no conforme con Quad4: el nodo 7, en el punto medio del
    lado 2-3 del elemento grueso, cuelga de sus extremos con
    ``u_7 = ½·u_2 + ½·u_3`` (en ``ux`` y ``uy``). ``escala`` multiplica todas
    las longitudes (1.0 → metros, 1e-5 → 10 µm; masas nodales ~1e-13 kg)."""
    dom = Domain()
    mat = Elastic2D(E=1.7e11, nu=0.28, density=2330.0)
    P = {1: (0, 0), 2: (1, 0), 3: (1, 1), 4: (0, 1), 5: (2, 0), 6: (2, .5),
         7: (1, .5), 8: (2, 1)}
    n = {k: dom.add_node(k, [escala * v[0], escala * v[1]]) for k, v in P.items()}
    esp = 0.1 * escala
    dom.add_element(Quad4(1, [n[1], n[2], n[3], n[4]], mat, thickness=esp))
    dom.add_element(Quad4(2, [n[2], n[5], n[6], n[7]], mat, thickness=esp))
    dom.add_element(Quad4(3, [n[7], n[6], n[8], n[3]], mat, thickness=esp))
    if apoyos:
        for k in (1, 4):
            n[k].fix_dof("ux", 0.0)
            n[k].fix_dof("uy", 0.0)
    if colgante:
        for d in ("ux", "uy"):
            dom.add_linear_constraint((7, d), [(2, d), (3, d)], [0.5, 0.5])
    dom.generate_equation_numbers(verbose=False)
    return dom, n


def _pasos_colgante(escala):
    # Δt ≈ 0.05·L/c con c = √(E/ρ) ≈ 8500 m/s.
    dt = 0.05 * escala / 8500.0
    return 20 * dt, dt


@pytest.mark.parametrize("escala", [1.0, 1.0e-5], ids=["1m", "10um"])
def test_central_difference_rechaza_mpc_multimaestro_con_masa(escala):
    """El nodo colgante hace TᵀMT no diagonal: rechazo con la causa correcta,
    también a microescala (umbral relativo, no absoluto)."""
    dom, n = _nodo_colgante(escala)
    f = np.zeros(dom.total_dofs)
    f[n[8].dofs["ux"]] = 1.0e-3 * escala ** 2
    t_end, dt = _pasos_colgante(escala)
    solver = CentralDifferenceSolver(Assembler(dom), t_end, dt, F_func=lambda t: f)
    with pytest.raises(ValueError, match="varios maestros") as exc:
        solver.solve()
    texto = str(exc.value)
    assert "no es estrictamente diagonal" in texto
    assert "nodo 7 (ux)" in texto and "nodo 7 (uy)" in texto
    assert "Frame3D" not in texto


@pytest.mark.parametrize("escala", [1.0, 1.0e-5], ids=["1m", "10um"])
def test_central_difference_sin_mpc_acepta_el_mismo_modelo(escala):
    """Control: sin el nodo colgante el explícito conserva Σm·ü = ΣF."""
    dom, n = _nodo_colgante(escala, colgante=False)
    F0 = 1.0e-3 * escala ** 2
    f = np.zeros(dom.total_dofs)
    f[n[8].dofs["ux"]] = F0
    t_end, dt = _pasos_colgante(escala)
    asm = Assembler(dom)
    r = CentralDifferenceSolver(asm, t_end, dt, F_func=lambda t: f).solve()
    m = asm.assemble_mass_matrix(lumping="lumped").diagonal()
    ux = np.array([nd.dofs["ux"] for nd in dom.nodes.values()])
    np.testing.assert_allclose((m[ux, None] * r.uddot_history[ux]).sum(axis=0) / F0, 1.0,
                               rtol=1e-10)


@pytest.mark.parametrize("escala", [1.0, 1.0e-5], ids=["1m", "10um"])
def test_newmark_concentrada_con_mpc_multimaestro_conserva_cantidad_de_movimiento(escala):
    """Oráculo exacto sin modelo fundido: sólido libre, fuerza constante; con la
    masa reducida completa TᵀMT, Σ_i m_i·ü_x,i = ΣF_x en todo instante (la
    traslación rígida está en el espacio de T porque ½ + ½ = 1)."""
    dom, n = _nodo_colgante(escala)
    F0 = 1.0e-3 * escala ** 2
    f = np.zeros(dom.total_dofs)
    f[n[8].dofs["ux"]] = F0
    t_end, dt = _pasos_colgante(escala)
    asm = Assembler(dom)
    r = NewmarkSolver(asm, t_end, dt, F_func=lambda t: f, lumping="lumped").solve()
    m = asm.assemble_mass_matrix(lumping="lumped").diagonal()
    ux = np.array([nd.dofs["ux"] for nd in dom.nodes.values()])
    np.testing.assert_allclose((m[ux, None] * r.uddot_history[ux]).sum(axis=0) / F0, 1.0,
                               rtol=1e-8)
    # La MPC se cumple en la historia.
    for d in ("ux", "uy"):
        np.testing.assert_allclose(
            r.u_history[n[7].dofs[d]],
            0.5 * (r.u_history[n[2].dofs[d]] + r.u_history[n[3].dofs[d]]),
            rtol=0.0, atol=1e-12 * np.abs(r.u_history).max())


def test_central_difference_acepta_rodillo_oblicuo_de_un_maestro():
    """Un apoyo oblicuo ``u_y = tanθ·u_x`` (un solo maestro, coeficiente ≠ 1)
    mantiene TᵀMT diagonal: el explícito lo admite y cumple la restricción."""
    dom, n = _nodo_colgante(1.0, colgante=False, apoyos=True)
    tan30 = np.tan(np.radians(30.0))
    dom.add_linear_constraint((5, "uy"), [(5, "ux")], [tan30])
    dom.generate_equation_numbers(verbose=False)
    f = np.zeros(dom.total_dofs)
    f[n[8].dofs["uy"]] = -1.0e4
    t_end, dt = _pasos_colgante(1.0)
    r = CentralDifferenceSolver(Assembler(dom), 10 * t_end, dt, F_func=lambda t: f).solve()
    np.testing.assert_allclose(r.u_history[n[5].dofs["uy"]], tan30 * r.u_history[n[5].dofs["ux"]],
                               rtol=0.0, atol=1e-15)
    assert np.abs(r.u_history).max() > 0.0


def test_central_difference_admite_esclavo_multimaestro_sin_masa():
    """Un esclavo con varios maestros pero sin masa no acopla TᵀMT: no se
    rechaza de antemano por tener varias entradas en su fila de T."""
    import scipy.sparse as sp
    M = sp.diags([2.0, 3.0, 0.0]).tocsr()
    T = sp.csr_matrix(np.array([[1.0, 0.0], [0.0, 1.0], [0.5, 0.5]]))
    M_red = (T.T @ M @ T).tocsr()
    inv = CentralDifferenceSolver._invert_diagonal_mass(M_red, T=T, M=M)
    np.testing.assert_allclose(inv, [0.5, 1.0 / 3.0])
    # Con masa en el esclavo sí se rechaza, y el mensaje lo nombra.
    M2 = sp.diags([2.0, 3.0, 1.0]).tocsr()
    with pytest.raises(ValueError, match="varios maestros"):
        CentralDifferenceSolver._invert_diagonal_mass((T.T @ M2 @ T).tocsr(), T=T, M=M2)


def test_central_difference_masa_no_diagonal_sin_mpc_conserva_el_mensaje_de_frame3d():
    """Sin filas multimaestro en T, la causa diagnosticada sigue siendo la del
    Frame3D oblicuo (masa global bloque-diagonal)."""
    import scipy.sparse as sp
    M = sp.csr_matrix(np.array([[1.0e-6, 1.0e-9], [1.0e-9, 1.0e-6]]))
    T = sp.identity(2, format="csr")
    with pytest.raises(ValueError, match="Frame3D"):
        CentralDifferenceSolver._invert_diagonal_mass(M, T=T, M=M)


# ----------------------------------------------------------------------------
# Condiciones iniciales incoherentes con una MPC
# ----------------------------------------------------------------------------

class _CapturaLog(logging.Handler):
    """Captura los avisos del logger ``solidum`` (``propagate = False``, ADR
    0005: ``caplog`` no los ve)."""

    def __init__(self):
        super().__init__(level=logging.WARNING)
        self.mensajes: list[str] = []

    def emit(self, record):
        self.mensajes.append(record.getMessage())


@contextlib.contextmanager
def _capturar_avisos():
    handler = _CapturaLog()
    logger = logging.getLogger("solidum")
    logger.addHandler(handler)
    try:
        yield handler
    finally:
        logger.removeHandler(handler)


def _velocidad_inicial(dom, nodos, en_duplicado: bool):
    """Perfil v_y(x) = x/L en todos los nodos; en el duplicado, el del nodo
    central si ``en_duplicado`` (coherente) o cero (incoherente)."""
    v0 = np.zeros(dom.total_dofs)
    for i, nodo in nodos.items():
        v0[nodo.dofs["uy"]] = i / NEL
    if DUP_ID in dom.nodes and en_duplicado:
        v0[dom.nodes[DUP_ID].dofs["uy"]] = K_MID / NEL
    return v0


def _corto(nombre):
    """Transitorio corto para las pruebas de condiciones iniciales."""
    if nombre == "CentralDifferenceSolver":
        return lambda asm, **kw: CentralDifferenceSolver(asm, 20 * DT_CD, DT_CD, **kw)
    cls = {"NewmarkSolver": NewmarkSolver, "HHTSolver": HHTSolver,
           "NewtonNewmarkSolver": NewtonNewmarkSolver, "NewtonHHTSolver": NewtonHHTSolver}[nombre]
    return lambda asm, **kw: cls(asm, 5 * DT, DT, lumping="lumped", **kw)


@pytest.mark.parametrize("solver", list(TRANSITORIOS))
def test_velocidad_inicial_incoherente_con_la_mpc_avisa(solver):
    """``u̇0`` en el nodo central pero no en su duplicado (esclavo a la
    izquierda): el esclavo se reconstruye desde el maestro (0) y el solver
    avisa nombrando el DOF."""
    dom, nodos = _mensula(True, "izquierda")
    v0 = _velocidad_inicial(dom, nodos, en_duplicado=False)
    with _capturar_avisos() as log:
        r = _corto(solver)(Assembler(dom), u0_dot=v0).solve()
    avisos = [m for m in log.mensajes if "restricción multipunto" in m]
    assert len(avisos) == 1, log.mensajes
    assert "u0_dot" in avisos[0]
    assert f"nodo {nodos[K_MID].id} (uy)" in avisos[0]
    # Conducta documentada: el esclavo toma T·u̇0[libres] (el valor del maestro).
    k_uy = nodos[K_MID].dofs["uy"]
    assert r.udot_history[k_uy, 0] == pytest.approx(0.0, abs=1e-15)


@pytest.mark.parametrize("solver", list(TRANSITORIOS))
def test_condiciones_iniciales_coherentes_no_avisan(solver):
    dom, nodos = _mensula(True, "izquierda")
    v0 = _velocidad_inicial(dom, nodos, en_duplicado=True)
    u0 = 1.0e-4 * v0
    with _capturar_avisos() as log:
        _corto(solver)(Assembler(dom), u0=u0, u0_dot=v0).solve()
    assert not [m for m in log.mensajes if "restricción multipunto" in m], log.mensajes


def test_desplazamiento_inicial_incoherente_avisa():
    dom, nodos = _mensula(True, "derecha")
    u0 = np.zeros(dom.total_dofs)
    u0[dom.nodes[DUP_ID].dofs["ux"]] = 1.0e-4   # sólo en el esclavo
    with _capturar_avisos() as log:
        NewmarkSolver(Assembler(dom), 5 * DT, DT, u0=u0).solve()
    avisos = [m for m in log.mensajes if "restricción multipunto" in m]
    assert len(avisos) == 1 and "condición inicial u0 " in avisos[0]
    assert f"nodo {DUP_ID} (ux)" in avisos[0]


def test_condicion_inicial_en_un_apoyo_dirichlet_no_avisa():
    """Las filas Dirichlet (sin maestros) quedan fuera del aviso: su proyección
    sobre el valor prescrito es deliberada."""
    dom, nodos = _mensula(True, "derecha")
    v0 = _velocidad_inicial(dom, nodos, en_duplicado=True)
    v0[nodos[0].dofs["ux"]] = 5.0   # apoyo empotrado
    with _capturar_avisos() as log:
        NewmarkSolver(Assembler(dom), 5 * DT, DT, u0_dot=v0).solve()
    assert not [m for m in log.mensajes if "restricción multipunto" in m], log.mensajes
