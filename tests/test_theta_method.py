"""Tests de ``ThetaMethodSolver`` — integración temporal θ (Etapa 8).

Cubre los ``acceptance`` de ``docs/specs/ThetaMethodSolver.md``.

El eje del blindaje no es que el solver "dé un número": es que **distinga
realmente los dos esquemas** que ofrece. Un θ-method mal implementado
—signos cruzados en el término ``(1−θ)``, factor ``Δt`` mal colocado—
todavía converge al estacionario correcto, porque en el límite ``t → ∞``
el término temporal desaparece. Lo que no sobrevive a un error así es la
**tasa de convergencia temporal**: por eso el test de orden (2 para
Crank-Nicolson, 1 para el resto) es el más discriminante de la suite.

Materiales con coeficientes físicos reales (acero) en vez de unitarios,
para que un factor perdido no se camufle multiplicando por 1.
"""
from __future__ import annotations

import contextlib
import logging

import numpy as np
import pytest
from scipy.linalg import expm

from solidum.core.domain import Domain
from solidum.core.node import Node
from solidum.elements.thermal.hex8_thermal import Hex8Thermal
from solidum.elements.thermal.quad4_thermal import Quad4Thermal
from solidum.materials.thermal_conduction import ThermalConduction
from solidum.math.assembly import Assembler
from solidum.math.solvers.linear import LinearSolver
from solidum.math.solvers.theta_method import ThetaMethodSolver
from solidum.results import ThermalTransientResult

# Acero estructural: coeficientes no unitarios a propósito.
K_ACERO = 45.0       # W/(m·K)
C_ACERO = 460.0      # J/(kg·K)
RHO_ACERO = 7850.0   # kg/m³
ALPHA_ACERO = K_ACERO / (RHO_ACERO * C_ACERO)   # ≈ 1.246e-5 m²/s



class _CapturaLog(logging.Handler):
    """Captura los mensajes del logger ``solidum``.

    ``caplog`` de pytest no sirve aquí: el logger del proyecto declara
    ``propagate = False`` (ADR 0005) para no contaminar la aplicación que
    embeba Solidum, así que sus registros nunca llegan al logger raíz que
    ``caplog`` intercepta. Se engancha un handler directamente.
    """

    def __init__(self):
        super().__init__(level=logging.WARNING)
        self.mensajes: list[str] = []

    def emit(self, record):
        self.mensajes.append(record.getMessage())

    @property
    def texto(self) -> str:
        return " ".join(self.mensajes)


@contextlib.contextmanager
def capturar_avisos():
    handler = _CapturaLog()
    logger = logging.getLogger("solidum")
    logger.addHandler(handler)
    try:
        yield handler
    finally:
        logger.removeHandler(handler)


# ----------------------------------------------------------------------
# Utilidades de malla
# ----------------------------------------------------------------------

def _material(dim: int = 2) -> ThermalConduction:
    return ThermalConduction(
        k=K_ACERO, c=C_ACERO, density=RHO_ACERO, dim=dim,
    )


def barra_2d(nx: int, L: float, *, H: float = 0.1, thickness: float = 1.0,
             T_izq: float | None = None, T_der: float | None = None):
    """Una fila de ``Quad4Thermal`` a lo largo de ``x``.

    Modelo unidimensional en la práctica: las caras superior e inferior
    quedan sin BC, es decir adiabáticas (Neumann homogéneo), así que el
    campo sólo varía con ``x``.
    """
    mat = _material(2)
    dom = Domain()
    nid: dict[tuple[int, int], int] = {}
    c = 1
    for j in range(2):
        for i in range(nx + 1):
            dom.add_node(c, [L * i / nx, H * j])
            nid[(i, j)] = c
            c += 1
    for i in range(nx):
        ns = [dom.nodes[nid[(i, 0)]], dom.nodes[nid[(i + 1, 0)]],
              dom.nodes[nid[(i + 1, 1)]], dom.nodes[nid[(i, 1)]]]
        dom.add_element(Quad4Thermal(i + 1, ns, mat, thickness=thickness))
    if T_izq is not None:
        for j in range(2):
            dom.nodes[nid[(0, j)]].fix_dof("T", T_izq)
    if T_der is not None:
        for j in range(2):
            dom.nodes[nid[(nx, j)]].fix_dof("T", T_der)
    dom.generate_equation_numbers()
    return dom, nid


def dofs_fila(dom: Domain, nid: dict, nx: int, fila: int = 0) -> list[int]:
    return [dom.nodes[nid[(i, fila)]].dofs["T"] for i in range(nx + 1)]


# ----------------------------------------------------------------------
# Verificación — contra solución analítica o el solver estacionario
# ----------------------------------------------------------------------

class TestOrdenDeConvergenciaTemporal:
    """El test que realmente distingue Crank-Nicolson de Euler implícito.

    Se compara contra la solución **exacta del sistema semidiscreto**
    ``T(t) = exp(−C⁻¹K t)·T₀``, no contra la solución analítica del
    continuo. La razón es que el error total mezcla dos fuentes —
    discretización espacial (fija) y temporal (la que se refina) — y el
    error espacial actúa como suelo que enmascara la tasa temporal en
    cuanto el paso se hace pequeño. Comparando contra la solución exacta
    en el tiempo del **mismo** sistema de ODEs que el solver integra, sólo
    queda el error temporal y la pendiente log-log es la del esquema.
    """

    @staticmethod
    def _pendiente(theta: float, n_refinamientos: int = 4) -> float:
        nx = 6
        dom, _ = barra_2d(nx, L=0.2, T_izq=0.0, T_der=0.0)
        asm = Assembler(dom)
        asm.assemble_system()

        ndof = asm.ndof
        cs = asm.constraint_set
        T_op, _g = cs.build(ndof)
        libres = cs.free_dofs(ndof)
        K_ff = (T_op.T @ asm.K_global @ T_op).toarray()
        C_ff = (T_op.T @ asm.assemble_mass_matrix(lumping="lumped")
                @ T_op).toarray()

        T0, t_end = 100.0, 20.0
        exacta = expm(-np.linalg.solve(C_ff, K_ff) * t_end) @ np.full(
            K_ff.shape[0], T0,
        )

        errs, dts = [], []
        for m in range(n_refinamientos):
            n = 20 * 2 ** m
            dt = t_end / n
            r = ThetaMethodSolver(
                asm, dt=dt, n_steps=n, T_initial=T0, theta=theta,
                lumping="lumped",
            ).solve()
            errs.append(np.max(np.abs(r.T_final[libres] - exacta)))
            dts.append(dt)
        return float(np.polyfit(np.log(dts), np.log(errs), 1)[0])

    def test_euler_implicito_es_primer_orden(self):
        p = self._pendiente(1.0)
        assert abs(p - 1.0) < 0.15, f"pendiente {p:.4f}, esperada ≈ 1"

    def test_crank_nicolson_es_segundo_orden(self):
        p = self._pendiente(0.5)
        assert abs(p - 2.0) < 0.15, f"pendiente {p:.4f}, esperada ≈ 2"

    def test_galerkin_es_primer_orden(self):
        """θ = 2/3 es incondicionalmente estable pero sólo de 1er orden.

        Confirma que la ganancia de orden es exclusiva de ``θ = 1/2``: el
        término de error de truncación de primer orden sólo se cancela en
        el punto medio exacto, no en cualquier θ ≥ 1/2.
        """
        p = self._pendiente(2.0 / 3.0)
        assert abs(p - 1.0) < 0.15, f"pendiente {p:.4f}, esperada ≈ 1"

    def test_crank_nicolson_es_mas_preciso_que_euler_al_mismo_paso(self):
        """Contrapartida del default robusto, medida en vez de afirmada."""
        nx = 6
        dom, _ = barra_2d(nx, L=0.2, T_izq=0.0, T_der=0.0)
        asm = Assembler(dom)
        asm.assemble_system()
        cs = asm.constraint_set
        T_op, _ = cs.build(asm.ndof)
        libres = cs.free_dofs(asm.ndof)
        K_ff = (T_op.T @ asm.K_global @ T_op).toarray()
        C_ff = (T_op.T @ asm.assemble_mass_matrix(lumping="lumped")
                @ T_op).toarray()
        T0, t_end, n = 100.0, 20.0, 40
        exacta = expm(-np.linalg.solve(C_ff, K_ff) * t_end) @ np.full(
            K_ff.shape[0], T0,
        )
        err = {}
        for th in (1.0, 0.5):
            r = ThetaMethodSolver(asm, dt=t_end / n, n_steps=n,
                                   T_initial=T0, theta=th).solve()
            err[th] = np.max(np.abs(r.T_final[libres] - exacta))
        assert err[0.5] < err[1.0] / 50.0, (
            f"CN {err[0.5]:.3e} no es sustancialmente mejor que "
            f"Euler {err[1.0]:.3e} al mismo Δt"
        )


class TestConvergenciaAlEstacionario:
    """Integrar hasta ``t ≫ L²/α`` debe reproducir el ``LinearSolver``.

    El tiempo de difusión global ``L²/α`` es la escala que gobierna cuándo
    el transitorio ha terminado. Integrar por debajo de ella no es un error
    del esquema sino un análisis inacabado — es justo lo que el diagnóstico
    de ``Δt`` característico busca prevenir, visto desde el otro extremo.
    """

    def test_perfil_lineal_coincide_con_el_solver_estacionario(self):
        nx, L = 8, 2.0
        dom, nid = barra_2d(nx, L=L, T_izq=100.0, T_der=20.0)
        asm = Assembler(dom)
        T_estac = LinearSolver(asm).solve(np.zeros(dom.total_dofs))

        t_difusion = L * L / ALPHA_ACERO           # ≈ 3.21e5 s
        n_steps = 400
        dt = 5.0 * t_difusion / n_steps            # t_end = 5·L²/α
        r = ThetaMethodSolver(
            asm, dt=dt, n_steps=n_steps, T_initial=20.0, theta=1.0,
        ).solve()

        assert np.allclose(r.T_final, T_estac, rtol=0, atol=1e-9)

    def test_el_perfil_convergido_es_lineal_en_x(self):
        nx, L = 8, 2.0
        dom, nid = barra_2d(nx, L=L, T_izq=100.0, T_der=20.0)
        asm = Assembler(dom)
        t_difusion = L * L / ALPHA_ACERO
        r = ThetaMethodSolver(
            asm, dt=5.0 * t_difusion / 400, n_steps=400,
            T_initial=20.0, theta=1.0,
        ).solve()
        perfil = r.T_final[dofs_fila(dom, nid, nx)]
        exacto = 100.0 + (20.0 - 100.0) * np.linspace(0.0, 1.0, nx + 1)
        assert np.allclose(perfil, exacto, rtol=0, atol=1e-9)

    def test_transitorio_inacabado_no_ha_convergido_todavia(self):
        """Contraste explícito: por debajo de ``L²/α`` el campo aún evoluciona.

        No es un fallo del esquema. Se blinda para que un futuro cambio no
        haga "converger" antes de tiempo — lo que sí indicaría un error.
        """
        nx, L = 8, 2.0
        dom, _ = barra_2d(nx, L=L, T_izq=100.0, T_der=20.0)
        asm = Assembler(dom)
        T_estac = LinearSolver(asm).solve(np.zeros(dom.total_dofs))
        t_difusion = L * L / ALPHA_ACERO
        r = ThetaMethodSolver(
            asm, dt=0.6 * t_difusion / 400, n_steps=400,
            T_initial=20.0, theta=1.0,
        ).solve()
        assert np.max(np.abs(r.T_final - T_estac)) > 1e-3


class TestConduccionSemiInfinita:
    """Carslaw-Jaeger §2.4 — ``T(x,t) = T_s + (T_0−T_s)·erf(x/(2√(αt)))``.

    Es el único test de la suite que compara contra una solución analítica
    del **continuo**, no del sistema semidiscreto: valida la cadena completa
    elemento + material + integrador contra física conocida.

    Hipótesis de medio semi-infinito: sólo vale mientras el frente térmico
    no alcance el extremo opuesto. Se comprueba explícitamente que la
    profundidad de penetración ``2√(αt)`` es pequeña frente a ``L``.
    """

    def test_perfil_coincide_con_la_funcion_error(self):
        erf = pytest.importorskip("scipy.special").erf

        L, nx = 1.0, 60
        T0, Ts, t_obs = 0.0, 100.0, 100.0
        dom, nid = barra_2d(nx, L=L, T_izq=Ts, T_der=T0)
        asm = Assembler(dom)

        penetracion = 2.0 * np.sqrt(ALPHA_ACERO * t_obs)
        assert penetracion < 0.15 * L, (
            "el frente alcanza el extremo: la hipótesis semi-infinita no vale"
        )

        n = 400
        r = ThetaMethodSolver(
            asm, dt=t_obs / n, n_steps=n, T_initial=T0, theta=1.0,
        ).solve()

        xs = np.linspace(0.0, L, nx + 1)
        num = r.T_final[dofs_fila(dom, nid, nx)]
        ana = Ts + (T0 - Ts) * erf(xs / (2.0 * np.sqrt(ALPHA_ACERO * t_obs)))

        # Error normalizado al salto de temperatura impuesto: en la cola
        # ambos campos son ~0 y un error relativo puntual carece de sentido.
        err = np.max(np.abs(num - ana)) / abs(Ts - T0)
        assert err < 2.0e-2, f"error normalizado {err:.4e}"


# ----------------------------------------------------------------------
# Estabilidad, L-estabilidad y principio del máximo
# ----------------------------------------------------------------------

class TestEstabilidad:

    def test_theta_1_no_diverge_con_paso_enorme(self):
        """Incondicionalmente estable: ``Δt`` tres órdenes sobre ``h²/α``."""
        nx = 8
        dom, _ = barra_2d(nx, L=2.0, T_izq=100.0, T_der=20.0)
        asm = Assembler(dom)
        dt_car = ThetaMethodSolver(
            asm, dt=1.0, n_steps=1, T_initial=20.0,
        ).characteristic_dt()
        r = ThetaMethodSolver(
            asm, dt=1000.0 * dt_car, n_steps=30, T_initial=20.0, theta=1.0,
        ).solve()
        assert r.converged
        assert np.all(np.isfinite(r.T_history))
        mn, mx = r.extremes()
        assert 20.0 - 1e-9 <= mn and mx <= 100.0 + 1e-9

    def test_theta_0_diverge_con_el_mismo_paso(self):
        """Contraste explícito: el explícito sí tiene límite de estabilidad."""
        nx = 8
        dom, _ = barra_2d(nx, L=2.0, T_izq=100.0, T_der=20.0)
        asm = Assembler(dom)
        dt_car = ThetaMethodSolver(
            asm, dt=1.0, n_steps=1, T_initial=20.0,
        ).characteristic_dt()
        r = ThetaMethodSolver(
            asm, dt=1000.0 * dt_car, n_steps=30, T_initial=20.0, theta=0.0,
        ).solve()
        assert not r.converged, "θ=0 debería divergir con este Δt"

    def test_theta_1_converge_monotonamente_al_estacionario(self):
        nx = 8
        dom, nid = barra_2d(nx, L=2.0, T_izq=100.0, T_der=20.0)
        asm = Assembler(dom)
        r = ThetaMethodSolver(
            asm, dt=2000.0, n_steps=200, T_initial=20.0, theta=1.0,
        ).solve()
        dof = dom.nodes[nid[(2, 0)]].dofs["T"]
        historia = r.temperature_at(dof)
        assert np.all(np.diff(historia) >= -1e-12), (
            "el nodo debe calentarse monótonamente hacia el estacionario"
        )


class TestLEstabilidadAnteEscalon:
    """Documenta —no corrige— la diferencia entre A-estable y L-estable.

    Es el corazón físico de la spec: un esquema estable puede producir un
    resultado **cualitativamente imposible**. Crank-Nicolson amplifica los
    modos altos con factor → −1 en vez de → 0, así que un escalón de
    temperatura le produce oscilación amortiguada con valores fuera del
    rango de los datos.

    Los tests fijan el comportamiento observado para que un cambio futuro
    en la formulación se note: si Euler implícito empezara a oscilar, sería
    un error; si Crank-Nicolson dejara de hacerlo con este Δt, también
    (indicaría que el término ``(1−θ)`` se perdió).
    """

    @staticmethod
    def _correr(theta: float, lumping: str = "lumped"):
        nx = 10
        dom, nid = barra_2d(nx, L=0.5, T_izq=100.0)   # extremo dcho adiabático
        asm = Assembler(dom)
        dt_car = ThetaMethodSolver(
            asm, dt=1.0, n_steps=1, T_initial=0.0,
        ).characteristic_dt()
        r = ThetaMethodSolver(
            asm, dt=5.0 * dt_car, n_steps=25, T_initial=0.0,
            theta=theta, lumping=lumping,
        ).solve()
        return dom, nid, r

    def test_euler_implicito_respeta_el_principio_del_maximo(self):
        """Ningún nodo, en ningún paso, fuera de ``[T_0, T_pared]``."""
        _dom, _nid, r = self._correr(1.0)
        mn, mx = r.extremes()
        assert mn >= 0.0 - 1e-12, f"T mínima {mn} < T_0 = 0"
        assert mx <= 100.0 + 1e-12, f"T máxima {mx} > T_pared = 100"

    def test_euler_implicito_es_monotono_en_un_nodo_interior(self):
        dom, nid, r = self._correr(1.0)
        h = r.temperature_at(dom.nodes[nid[(3, 0)]].dofs["T"])
        assert np.all(np.diff(h) >= -1e-12), "Euler implícito no debe oscilar"

    def test_crank_nicolson_se_sale_del_rango_de_los_datos(self):
        """Sobrepasa ``T_pared``: viola el principio del máximo.

        No es un fallo del solver — es la propiedad conocida de un esquema
        A-estable no L-estable, y la razón documentada de que el default sea
        θ = 1. El test la cuantifica para que el usuario que elija θ = 0.5
        sepa exactamente a qué se expone.
        """
        _dom, _nid, r = self._correr(0.5)
        _mn, mx = r.extremes()
        assert mx > 100.0 + 1.0, (
            f"T máxima {mx:.4f} no supera T_pared: se esperaba sobrepaso"
        )

    def test_crank_nicolson_oscila_en_un_nodo_interior(self):
        dom, nid, r = self._correr(0.5)
        h = r.temperature_at(dom.nodes[nid[(3, 0)]].dofs["T"])
        assert np.min(np.diff(h)) < -1e-6, (
            "se esperaba al menos un decremento (oscilación) con θ = 0.5"
        )

    def test_capacidad_consistente_no_salva_a_crank_nicolson(self):
        """Los dos ejes del problema son independientes.

        ``lumped`` ataca la oscilación espacial; ``θ = 1`` la temporal. Con
        Crank-Nicolson la oscilación temporal persiste aunque la capacidad
        sea lumped o consistente: son fenómenos distintos y el default de
        cada eje se justifica por separado.
        """
        for lumping in ("lumped", "consistent"):
            _d, _n, r = self._correr(0.5, lumping=lumping)
            _mn, mx = r.extremes()
            assert mx > 100.0 + 1.0, f"sin sobrepaso con lumping={lumping}"


# ----------------------------------------------------------------------
# Condición inicial
# ----------------------------------------------------------------------

class TestCondicionInicial:

    def test_escalar_y_vector_equivalente_dan_el_mismo_resultado(self):
        nx = 8
        dom, _ = barra_2d(nx, L=2.0, T_izq=100.0, T_der=20.0)
        asm = Assembler(dom)
        rA = ThetaMethodSolver(
            asm, dt=100.0, n_steps=10, T_initial=20.0,
        ).solve()
        rB = ThetaMethodSolver(
            asm, dt=100.0, n_steps=10,
            T_initial=np.full(dom.total_dofs, 20.0),
        ).solve()
        assert np.allclose(rA.T_history, rB.T_history, rtol=0, atol=1e-14)

    def test_campo_inicial_no_uniforme_arranca_del_estado_dado(self):
        """Caso de uso: continuar desde un estacionario previo."""
        nx = 8
        dom, nid = barra_2d(nx, L=2.0, T_izq=100.0, T_der=20.0)
        asm = Assembler(dom)
        T_estac = LinearSolver(asm).solve(np.zeros(dom.total_dofs))
        r = ThetaMethodSolver(
            asm, dt=1000.0, n_steps=5, T_initial=T_estac, theta=1.0,
        ).solve()
        # Arrancar del estacionario ⇒ el campo no se mueve.
        assert np.allclose(r.T_final, T_estac, rtol=0, atol=1e-9)

    def test_vector_de_tamano_incorrecto_es_rechazado(self):
        nx = 8
        dom, _ = barra_2d(nx, L=2.0, T_izq=100.0, T_der=20.0)
        asm = Assembler(dom)
        s = ThetaMethodSolver(
            asm, dt=100.0, n_steps=2, T_initial=np.zeros(3),
        )
        with pytest.raises(ValueError, match="T_initial es un vector"):
            s.solve()

    def test_avisa_si_contradice_dirichlet_y_continua(self):
        """Un choque térmico es modelización legítima, no un error."""
        nx = 8
        dom, _ = barra_2d(nx, L=2.0, T_izq=100.0, T_der=20.0)
        asm = Assembler(dom)
        with capturar_avisos() as log:
            r = ThetaMethodSolver(
                asm, dt=100.0, n_steps=3, T_initial=999.0,
            ).solve()
        assert r.converged
        assert "T_initial contradice" in log.texto

    def test_el_aviso_menciona_la_oscilacion_si_no_es_L_estable(self):
        """Con θ ≠ 1 el aviso añade el riesgo concreto que corre el usuario."""
        nx = 8
        dom, _ = barra_2d(nx, L=2.0, T_izq=100.0, T_der=20.0)
        asm = Assembler(dom)
        with capturar_avisos() as log:
            ThetaMethodSolver(
                asm, dt=100.0, n_steps=3, T_initial=999.0, theta=0.5,
            ).solve()
        assert "L-estable" in log.texto

    def test_con_theta_1_el_aviso_no_menciona_oscilacion(self):
        """El aviso es específico: θ = 1 no puede oscilar, no debe insinuarlo."""
        nx = 8
        dom, _ = barra_2d(nx, L=2.0, T_izq=100.0, T_der=20.0)
        asm = Assembler(dom)
        with capturar_avisos() as log:
            ThetaMethodSolver(
                asm, dt=100.0, n_steps=3, T_initial=999.0, theta=1.0,
            ).solve()
        assert "T_initial contradice" in log.texto
        assert "L-estable" not in log.texto

    def test_sin_contradiccion_no_hay_aviso(self):
        nx = 8
        dom, _ = barra_2d(nx, L=2.0, T_izq=100.0, T_der=100.0)
        asm = Assembler(dom)
        with capturar_avisos() as log:
            ThetaMethodSolver(
                asm, dt=100.0, n_steps=3, T_initial=100.0,
            ).solve()
        assert "T_initial contradice" not in log.texto


# ----------------------------------------------------------------------
# Dirichlet variable en el tiempo
# ----------------------------------------------------------------------

class TestDirichletVariable:

    @staticmethod
    def _modelo_sinusoidal(omega: float, amplitud: float = 50.0):
        nx = 8
        dom, nid = barra_2d(nx, L=0.3, T_izq=0.0, T_der=0.0)
        asm = Assembler(dom)
        ndof = dom.total_dofs
        izq = [dom.nodes[nid[(0, j)]].dofs["T"] for j in range(2)]

        def g_de_t(t: float) -> np.ndarray:
            g = np.zeros(ndof)
            for d in izq:
                g[d] = amplitud * np.sin(omega * t)
            return g

        return dom, nid, asm, izq, g_de_t

    def test_la_frontera_sigue_exactamente_la_excitacion(self):
        omega = 2.0 * np.pi / 500.0
        dom, nid, asm, izq, g = self._modelo_sinusoidal(omega)
        s = ThetaMethodSolver(
            asm, dt=5.0, n_steps=300, T_initial=0.0, theta=1.0,
            dirichlet_func=g,
        )
        r = s.solve()
        esperado = 50.0 * np.sin(omega * r.t_history)
        assert np.allclose(r.temperature_at(izq[0]), esperado,
                            rtol=0, atol=1e-12)

    def test_la_factorizacion_se_reutiliza_en_todos_los_pasos(self):
        """Sólo cambia el acoplamiento ``K_fp·ḡ``; ``A_ff`` es constante."""
        omega = 2.0 * np.pi / 500.0
        _d, _n, asm, _izq, g = self._modelo_sinusoidal(omega)
        s = ThetaMethodSolver(
            asm, dt=5.0, n_steps=300, T_initial=0.0, dirichlet_func=g,
        )
        s.solve()
        assert s._n_factorizations == 1

    def test_el_interior_oscila_con_amplitud_atenuada(self):
        """Física esperada: la onda térmica se amortigua con la profundidad."""
        omega = 2.0 * np.pi / 500.0
        dom, nid, asm, izq, g = self._modelo_sinusoidal(omega)
        r = ThetaMethodSolver(
            asm, dt=5.0, n_steps=300, T_initial=0.0, dirichlet_func=g,
        ).solve()
        # Tras el arranque, medir amplitud pico a pico en régimen.
        interior = r.temperature_at(dom.nodes[nid[(4, 0)]].dofs["T"])[100:]
        amp_int = 0.5 * (interior.max() - interior.min())
        assert 0.0 < amp_int < 50.0, (
            f"amplitud interior {amp_int:.4f} debe estar atenuada respecto "
            "a los 50 de la frontera"
        )

    def test_dirichlet_constante_via_funcion_iguala_al_declarado(self):
        """``dirichlet_func`` devolviendo el ``g`` del modelo es un no-op.

        Blinda el término de capacidad ``C_fp·(ḡ_{n+1} − ḡ_n)``: si el valor
        prescrito no varía, ese término se anula y el resultado debe ser
        idéntico al del camino sin ``dirichlet_func``.
        """
        nx = 8
        dom, _ = barra_2d(nx, L=2.0, T_izq=100.0, T_der=20.0)
        asm = Assembler(dom)
        asm.assemble_system()
        _T_op, g0 = asm.constraint_set.build(asm.ndof)

        r_sin = ThetaMethodSolver(
            asm, dt=2000.0, n_steps=20, T_initial=20.0,
        ).solve()
        r_con = ThetaMethodSolver(
            asm, dt=2000.0, n_steps=20, T_initial=20.0,
            dirichlet_func=lambda t: g0,
        ).solve()
        assert np.allclose(r_sin.T_history, r_con.T_history,
                            rtol=0, atol=1e-12)

    def test_vector_de_tamano_incorrecto_es_rechazado(self):
        nx = 8
        dom, _ = barra_2d(nx, L=2.0, T_izq=100.0, T_der=20.0)
        asm = Assembler(dom)
        s = ThetaMethodSolver(
            asm, dt=100.0, n_steps=2, T_initial=20.0,
            dirichlet_func=lambda t: np.zeros(3),
        )
        with pytest.raises(ValueError, match="dirichlet_func"):
            s.solve()


# ----------------------------------------------------------------------
# Cargas
# ----------------------------------------------------------------------

class TestCargas:

    def test_fuente_constante_lleva_al_estacionario_del_linear_solver(self):
        """Con ``F`` constante el transitorio converge a ``K·T = F``."""
        nx, L = 6, 0.5
        dom, nid = barra_2d(nx, L=L, T_izq=0.0, T_der=0.0)
        asm = Assembler(dom)
        asm.assemble_system()
        ndof = asm.ndof

        Q = 5.0e4   # W/m³
        F = np.zeros(ndof)
        for el in dom.elements.values():
            fe = el.compute_body_source(Q)
            for a, n in enumerate(el.nodes):
                F[n.dofs["T"]] += fe[a]

        T_estac = LinearSolver(asm).solve(F)
        t_dif = L * L / ALPHA_ACERO
        r = ThetaMethodSolver(
            asm, dt=5.0 * t_dif / 300, n_steps=300, T_initial=0.0,
            theta=1.0, F_func=lambda t: F,
        ).solve()
        assert np.allclose(r.T_final, T_estac, rtol=0, atol=1e-8)

    def test_sin_F_func_el_campo_lo_mueve_solo_dirichlet(self):
        nx = 8
        dom, _ = barra_2d(nx, L=2.0, T_izq=100.0, T_der=20.0)
        asm = Assembler(dom)
        r1 = ThetaMethodSolver(
            asm, dt=1000.0, n_steps=10, T_initial=20.0,
        ).solve()
        ndof = dom.total_dofs
        r2 = ThetaMethodSolver(
            asm, dt=1000.0, n_steps=10, T_initial=20.0,
            F_func=lambda t: np.zeros(ndof),
        ).solve()
        assert np.allclose(r1.T_history, r2.T_history, rtol=0, atol=1e-14)

    def test_F_func_de_tamano_incorrecto_es_rechazado(self):
        nx = 8
        dom, _ = barra_2d(nx, L=2.0, T_izq=100.0, T_der=20.0)
        asm = Assembler(dom)
        s = ThetaMethodSolver(
            asm, dt=100.0, n_steps=2, T_initial=20.0,
            F_func=lambda t: np.zeros(3),
        )
        with pytest.raises(ValueError, match="F_func"):
            s.solve()


# ----------------------------------------------------------------------
# Diagnósticos declarativos
# ----------------------------------------------------------------------

class TestReporteDelOrdenEfectivo:
    """Contrapartida honesta de elegir robustez como default.

    El coste en precisión de ``θ = 1`` debe ser **consultable**, no una
    penalización silenciosa: por eso ``order`` es propiedad del solver
    (antes de correr) y campo del resultado (después).
    """

    @pytest.mark.parametrize("theta, orden", [
        (1.0, 1), (0.5, 2), (2.0 / 3.0, 1), (0.0, 1), (0.75, 1),
    ])
    def test_orden_reportado(self, theta, orden):
        nx = 4
        dom, _ = barra_2d(nx, L=1.0, T_izq=0.0, T_der=0.0)
        asm = Assembler(dom)
        s = ThetaMethodSolver(
            asm, dt=1.0, n_steps=1, T_initial=0.0, theta=theta,
        )
        assert s.order == orden

    def test_el_resultado_transporta_el_orden_y_los_parametros(self):
        nx = 4
        dom, _ = barra_2d(nx, L=1.0, T_izq=0.0, T_der=100.0)
        asm = Assembler(dom)
        r = ThetaMethodSolver(
            asm, dt=7.0, n_steps=3, T_initial=0.0, theta=0.5,
        ).solve()
        assert isinstance(r, ThermalTransientResult)
        assert r.order == 2
        assert r.theta == 0.5
        assert r.dt == 7.0
        assert r.n_steps == 3

    @pytest.mark.parametrize("theta, incond, l_estable", [
        (1.0, True, True), (0.5, True, False),
        (2.0 / 3.0, True, False), (0.0, False, False), (0.25, False, False),
    ])
    def test_banderas_de_estabilidad(self, theta, incond, l_estable):
        nx = 4
        dom, _ = barra_2d(nx, L=1.0, T_izq=0.0, T_der=0.0)
        asm = Assembler(dom)
        s = ThetaMethodSolver(
            asm, dt=1.0, n_steps=1, T_initial=0.0, theta=theta,
        )
        assert s.is_unconditionally_stable is incond
        assert s.is_l_stable is l_estable


class TestPasoCaracteristico:
    """``Δt ~ h²/α`` — información, no restricción."""

    def test_coincide_con_h_cuadrado_sobre_alfa(self):
        nx, L, H = 4, 0.4, 0.1
        dom, _ = barra_2d(nx, L=L, H=H, T_izq=0.0, T_der=100.0)
        asm = Assembler(dom)
        s = ThetaMethodSolver(asm, dt=1.0, n_steps=1, T_initial=0.0)
        h = np.hypot(L / nx, H)          # diagonal del bounding box
        assert s.characteristic_dt() == pytest.approx(
            h * h / ALPHA_ACERO, rel=1e-12,
        )

    def test_escala_inversamente_con_la_difusividad(self):
        """Un material que difunde el doble de rápido, la mitad de Δt."""
        nx = 4
        dom, _ = barra_2d(nx, L=0.4, T_izq=0.0, T_der=100.0)
        asm = Assembler(dom)
        dt1 = ThetaMethodSolver(
            asm, dt=1.0, n_steps=1, T_initial=0.0,
        ).characteristic_dt()

        # Mismo modelo con conductividad doble ⇒ difusividad doble.
        mat2 = ThermalConduction(k=2 * K_ACERO, c=C_ACERO,
                                  density=RHO_ACERO, dim=2)
        for el in dom.elements.values():
            el.material = mat2
        dt2 = ThetaMethodSolver(
            asm, dt=1.0, n_steps=1, T_initial=0.0,
        ).characteristic_dt()
        assert dt2 == pytest.approx(dt1 / 2.0, rel=1e-12)

    def test_no_estimable_sin_capacidad_declarada(self):
        """Material sin ``c``/``density`` ⇒ ``None``, no una excepción.

        Es un diagnóstico opcional: su ausencia no debe impedir un análisis
        estacionario ni romper la construcción del solver.
        """
        mat = ThermalConduction(k=K_ACERO, dim=2)   # sin c ni density
        dom = Domain()
        for i, xy in enumerate([[0, 0], [1, 0], [1, 1], [0, 1]], start=1):
            dom.add_node(i, [float(xy[0]), float(xy[1])])
        dom.add_element(Quad4Thermal(1, list(dom.nodes.values()), mat))
        dom.nodes[1].fix_dof("T", 0.0)
        dom.generate_equation_numbers()
        asm = Assembler(dom)
        s = ThetaMethodSolver(asm, dt=1.0, n_steps=1, T_initial=0.0)
        assert s.characteristic_dt() is None

    def test_avisa_si_el_paso_se_salta_el_transitorio(self):
        nx = 8
        dom, _ = barra_2d(nx, L=2.0, T_izq=100.0, T_der=20.0)
        asm = Assembler(dom)
        dt_car = ThetaMethodSolver(
            asm, dt=1.0, n_steps=1, T_initial=20.0,
        ).characteristic_dt()
        with capturar_avisos() as log:
            ThetaMethodSolver(
                asm, dt=100.0 * dt_car, n_steps=3, T_initial=20.0,
            ).solve()
        assert "frente térmico" in log.texto

    def test_no_avisa_si_el_paso_es_razonable(self):
        nx = 8
        dom, _ = barra_2d(nx, L=2.0, T_izq=100.0, T_der=20.0)
        asm = Assembler(dom)
        dt_car = ThetaMethodSolver(
            asm, dt=1.0, n_steps=1, T_initial=20.0,
        ).characteristic_dt()
        with capturar_avisos() as log:
            ThetaMethodSolver(
                asm, dt=0.5 * dt_car, n_steps=3, T_initial=20.0,
            ).solve()
        assert "frente térmico" not in log.texto


# ----------------------------------------------------------------------
# Almacenamiento y resultado
# ----------------------------------------------------------------------

class TestSubmuestreoDeSalida:

    def test_output_every_reduce_columnas_sin_alterar_la_integracion(self):
        nx = 8
        dom, _ = barra_2d(nx, L=2.0, T_izq=100.0, T_der=20.0)
        asm = Assembler(dom)
        completo = ThetaMethodSolver(
            asm, dt=100.0, n_steps=10, T_initial=20.0, output_every=1,
        ).solve()
        reducido = ThetaMethodSolver(
            asm, dt=100.0, n_steps=10, T_initial=20.0, output_every=3,
        ).solve()
        assert completo.T_history.shape[1] == 11
        assert reducido.T_history.shape[1] < 11
        # El campo final es idéntico: el submuestreo afecta al almacenamiento,
        # no al cálculo.
        assert np.allclose(reducido.T_final, completo.T_final,
                            rtol=0, atol=1e-13)

    def test_el_instante_final_siempre_se_almacena(self):
        """Aunque ``n_steps`` no sea múltiplo de ``output_every``."""
        nx = 8
        dom, _ = barra_2d(nx, L=2.0, T_izq=100.0, T_der=20.0)
        asm = Assembler(dom)
        r = ThetaMethodSolver(
            asm, dt=100.0, n_steps=10, T_initial=20.0, output_every=3,
        ).solve()
        assert r.t_history[0] == 0.0
        assert r.t_history[-1] == pytest.approx(1000.0)

    def test_n_steps_cuenta_pasos_integrados_no_columnas(self):
        nx = 8
        dom, _ = barra_2d(nx, L=2.0, T_izq=100.0, T_der=20.0)
        asm = Assembler(dom)
        r = ThetaMethodSolver(
            asm, dt=100.0, n_steps=10, T_initial=20.0, output_every=4,
        ).solve()
        assert r.n_steps == 10
        assert len(r.t_history) < 11


class TestApiDelResultado:

    @staticmethod
    def _resultado():
        nx = 8
        dom, nid = barra_2d(nx, L=2.0, T_izq=100.0, T_der=20.0)
        asm = Assembler(dom)
        r = ThetaMethodSolver(
            asm, dt=1000.0, n_steps=10, T_initial=20.0,
        ).solve()
        return dom, nid, r

    def test_T_final_es_la_ultima_columna(self):
        _d, _n, r = self._resultado()
        assert np.array_equal(r.T_final, r.T_history[:, -1])

    def test_temperature_at_devuelve_la_historia_de_un_dof(self):
        dom, nid, r = self._resultado()
        dof = dom.nodes[nid[(3, 0)]].dofs["T"]
        h = r.temperature_at(dof)
        assert h.shape == r.t_history.shape
        assert np.array_equal(h, r.T_history[dof, :])

    def test_extremes_devuelve_min_y_max_globales(self):
        _d, _n, r = self._resultado()
        mn, mx = r.extremes()
        assert mn == pytest.approx(float(np.min(r.T_history)))
        assert mx == pytest.approx(float(np.max(r.T_history)))

    def test_el_resultado_es_inmutable(self):
        _d, _n, r = self._resultado()
        with pytest.raises(Exception):
            r.theta = 0.5

    def test_los_dofs_prescritos_llevan_el_valor_impuesto(self):
        dom, nid, r = self._resultado()
        dof_izq = dom.nodes[nid[(0, 0)]].dofs["T"]
        assert np.allclose(r.temperature_at(dof_izq), 100.0)


# ----------------------------------------------------------------------
# Validación de entradas
# ----------------------------------------------------------------------

class TestRechazoDeEntradasInvalidas:

    @staticmethod
    def _asm():
        dom, _ = barra_2d(4, L=1.0, T_izq=0.0, T_der=100.0)
        return Assembler(dom)

    @pytest.mark.parametrize("theta", [-0.1, 1.1, 2.0, -1.0])
    def test_theta_fuera_de_rango(self, theta):
        with pytest.raises(ValueError, match="theta"):
            ThetaMethodSolver(
                self._asm(), dt=1.0, n_steps=1, T_initial=0.0, theta=theta,
            )

    @pytest.mark.parametrize("dt", [0.0, -1.0, -1e-9])
    def test_dt_no_positivo(self, dt):
        with pytest.raises(ValueError, match="dt"):
            ThetaMethodSolver(self._asm(), dt=dt, n_steps=1, T_initial=0.0)

    @pytest.mark.parametrize("n", [0, -1])
    def test_n_steps_invalido(self, n):
        with pytest.raises(ValueError, match="n_steps"):
            ThetaMethodSolver(self._asm(), dt=1.0, n_steps=n, T_initial=0.0)

    @pytest.mark.parametrize("k", [0, -2])
    def test_output_every_invalido(self, k):
        with pytest.raises(ValueError, match="output_every"):
            ThetaMethodSolver(
                self._asm(), dt=1.0, n_steps=1, T_initial=0.0, output_every=k,
            )

    def test_sin_dirichlet_el_mensaje_nombra_el_modo_de_temperatura_uniforme(self):
        """Diagnóstico físico, no un fallo algebraico genérico."""
        mat = _material(2)
        dom = Domain()
        for i, xy in enumerate([[0, 0], [1, 0], [1, 1], [0, 1]], start=1):
            dom.add_node(i, [float(xy[0]), float(xy[1])])
        dom.add_element(Quad4Thermal(1, list(dom.nodes.values()), mat))
        dom.generate_equation_numbers()
        s = ThetaMethodSolver(
            Assembler(dom), dt=1.0, n_steps=2, T_initial=0.0,
        )
        with pytest.raises(ValueError, match="temperatura uniforme"):
            s.solve()

    def test_material_sin_capacidad_falla_ruidosamente(self):
        """Sin ``ρc`` el transitorio debe abortar, no producir basura."""
        mat = ThermalConduction(k=K_ACERO, dim=2)   # sin c ni density
        dom = Domain()
        for i, xy in enumerate([[0, 0], [1, 0], [1, 1], [0, 1]], start=1):
            dom.add_node(i, [float(xy[0]), float(xy[1])])
        dom.add_element(Quad4Thermal(1, list(dom.nodes.values()), mat))
        dom.nodes[1].fix_dof("T", 0.0)
        dom.generate_equation_numbers()
        s = ThetaMethodSolver(
            Assembler(dom), dt=1.0, n_steps=2, T_initial=0.0,
        )
        with pytest.raises(ValueError):
            s.solve()

    def test_material_sin_capacidad_orienta_hacia_el_solver_estacionario(self):
        """El transitorio necesita ``ρc``; el estacionario no.

        El mensaje debe orientar hacia el ``LinearSolver`` en vez de exigir
        datos que el usuario quizá no necesite para su análisis. El
        ``Assembler`` delega en ``ThermalMaterial.volumetric_capacity``, que
        conoce la física de su familia; sin esa delegación se emitiría el
        mensaje mecánico genérico, que aconseja ``density = 0.0`` — consejo
        físicamente incorrecto aquí, porque daría capacidad calorífica nula.
        """
        mat = ThermalConduction(k=K_ACERO, dim=2)   # sin c ni density
        dom = Domain()
        for i, xy in enumerate([[0, 0], [1, 0], [1, 1], [0, 1]], start=1):
            dom.add_node(i, [float(xy[0]), float(xy[1])])
        dom.add_element(Quad4Thermal(1, list(dom.nodes.values()), mat))
        dom.nodes[1].fix_dof("T", 0.0)
        dom.generate_equation_numbers()
        s = ThetaMethodSolver(
            Assembler(dom), dt=1.0, n_steps=2, T_initial=0.0,
        )
        with pytest.raises(ValueError, match="LinearSolver"):
            s.solve()

    def test_un_material_mecanico_sin_density_conserva_su_mensaje(self):
        """La delegación no altera el dominio mecánico.

        El ``Assembler`` sólo cede el mensaje a materiales que exponen
        ``volumetric_capacity``. Un material mecánico sin ``density`` sigue
        recibiendo el diagnóstico de ADR 0008 —donde ``density = 0.0`` **sí**
        es legítimo (penalty, restricción)—, intacto. Es la garantía de que
        el arreglo del mensaje térmico no tuvo efectos colaterales.
        """
        from solidum.elements.solid_2d import Quad4
        from solidum.materials.elastic_2d import Elastic2D

        mat = Elastic2D(E=210e9, nu=0.3)      # sin density (ADR 0008)
        assert not hasattr(mat, "volumetric_capacity")

        dom = Domain()
        for i, xy in enumerate([[0, 0], [1, 0], [1, 1], [0, 1]], start=1):
            dom.add_node(i, [float(xy[0]), float(xy[1])])
        dom.add_element(Quad4(1, list(dom.nodes.values()), mat))
        dom.generate_equation_numbers()

        with pytest.raises(ValueError, match="ADR 0008"):
            Assembler(dom).assemble_mass_matrix()


# ----------------------------------------------------------------------
# 3D y cross-check
# ----------------------------------------------------------------------

class TestTransitorio3D:
    """El integrador es agnóstico a la dimensión del elemento.

    Mismo problema resuelto con una fila de ``Quad4Thermal`` y con una capa
    de ``Hex8Thermal`` con las caras ``z`` adiabáticas: el transitorio debe
    coincidir paso a paso. Es el análogo temporal del cross-check 2D↔3D
    estacionario, y detecta un error en la tercera dimensión sin necesitar
    un mallador de geometría curva.
    """

    @staticmethod
    def _malla_3d(nx: int, L: float, H: float = 0.1, W: float = 0.1,
                   T_izq: float = 100.0, T_der: float = 20.0):
        mat = _material(3)
        dom = Domain()
        nid: dict[tuple[int, int, int], int] = {}
        c = 1
        for kz in range(2):
            for j in range(2):
                for i in range(nx + 1):
                    dom.add_node(c, [L * i / nx, H * j, W * kz])
                    nid[(i, j, kz)] = c
                    c += 1
        for i in range(nx):
            ns = [
                dom.nodes[nid[(i, 0, 0)]], dom.nodes[nid[(i + 1, 0, 0)]],
                dom.nodes[nid[(i + 1, 1, 0)]], dom.nodes[nid[(i, 1, 0)]],
                dom.nodes[nid[(i, 0, 1)]], dom.nodes[nid[(i + 1, 0, 1)]],
                dom.nodes[nid[(i + 1, 1, 1)]], dom.nodes[nid[(i, 1, 1)]],
            ]
            dom.add_element(Hex8Thermal(i + 1, ns, mat))
        for kz in range(2):
            for j in range(2):
                dom.nodes[nid[(0, j, kz)]].fix_dof("T", T_izq)
                dom.nodes[nid[(nx, j, kz)]].fix_dof("T", T_der)
        dom.generate_equation_numbers()
        return dom, nid

    def test_cross_check_2d_vs_3d_en_todo_el_transitorio(self):
        nx, L, W = 6, 2.0, 0.1
        dt, n = 2000.0, 40

        dom2, nid2 = barra_2d(nx, L=L, H=0.1, thickness=W,
                               T_izq=100.0, T_der=20.0)
        r2 = ThetaMethodSolver(
            Assembler(dom2), dt=dt, n_steps=n, T_initial=20.0, theta=1.0,
        ).solve()

        dom3, nid3 = _M3 = self._malla_3d(nx, L=L, H=0.1, W=W)
        r3 = ThetaMethodSolver(
            Assembler(dom3), dt=dt, n_steps=n, T_initial=20.0, theta=1.0,
        ).solve()

        for i in range(nx + 1):
            h2 = r2.temperature_at(dom2.nodes[nid2[(i, 0)]].dofs["T"])
            h3 = r3.temperature_at(dom3.nodes[nid3[(i, 0, 0)]].dofs["T"])
            assert np.allclose(h2, h3, rtol=0, atol=1e-9), (
                f"discrepancia 2D/3D en el nodo x={i}"
            )

    def test_3d_converge_al_estacionario_lineal(self):
        nx, L = 6, 2.0
        dom, nid = self._malla_3d(nx, L=L)
        asm = Assembler(dom)
        t_dif = L * L / ALPHA_ACERO
        r = ThetaMethodSolver(
            asm, dt=5.0 * t_dif / 300, n_steps=300, T_initial=20.0,
            theta=1.0,
        ).solve()
        perfil = np.array([
            r.T_final[dom.nodes[nid[(i, 0, 0)]].dofs["T"]]
            for i in range(nx + 1)
        ])
        exacto = 100.0 + (20.0 - 100.0) * np.linspace(0.0, 1.0, nx + 1)
        assert np.allclose(perfil, exacto, rtol=0, atol=1e-9)


# ----------------------------------------------------------------------
# Integración con la infraestructura del proyecto
# ----------------------------------------------------------------------

class TestIntegracionConElProyecto:

    def test_esta_registrado_en_el_solver_registry(self):
        from solidum.registry import SolverRegistry
        assert SolverRegistry.get("ThetaMethodSolver") is ThetaMethodSolver

    def test_declara_su_pipeline(self):
        assert ThetaMethodSolver.PIPELINE_KIND == "thermal_transient"

    def test_el_pipeline_esta_reconocido_por_run_yaml(self):
        from solidum.entry import _KNOWN_PIPELINE_KINDS
        assert "thermal_transient" in _KNOWN_PIPELINE_KINDS

    def test_run_thermal_transient_construye_el_solver_por_defecto(self):
        from solidum.entry import run_thermal_transient
        nx = 8
        dom, _ = barra_2d(nx, L=2.0, T_izq=100.0, T_der=20.0)
        r = run_thermal_transient(
            dom, dt=1000.0, n_steps=10, T_initial=20.0,
        )
        assert isinstance(r, ThermalTransientResult)
        assert dom.last_result is r
        assert r.theta == 1.0            # default L-estable

    def test_run_thermal_transient_exige_los_parametros_minimos(self):
        from solidum.entry import run_thermal_transient
        dom, _ = barra_2d(4, L=1.0, T_izq=0.0, T_der=100.0)
        with pytest.raises(ValueError, match="run_thermal_transient"):
            run_thermal_transient(dom, dt=1.0)

    def test_es_exportado_en_la_raiz_del_paquete(self):
        import solidum
        assert solidum.ThetaMethodSolver is ThetaMethodSolver
        assert solidum.ThermalTransientResult is ThermalTransientResult

    def test_no_acepta_rayleigh(self):
        """El amortiguamiento de Rayleigh no existe en 1er orden.

        Blinda una decisión de la spec: si alguien lo añadiera "por
        simetría" con ``NewmarkSolver`` estaría introduciendo un modelo de
        disipación que no corresponde a esta ecuación.
        """
        dom, _ = barra_2d(4, L=1.0, T_izq=0.0, T_der=100.0)
        with pytest.raises(TypeError):
            ThetaMethodSolver(
                Assembler(dom), dt=1.0, n_steps=1, T_initial=0.0,
                rayleigh={"alpha": 0.1, "beta": 0.0},
            )
