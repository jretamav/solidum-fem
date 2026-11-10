# solidum_fem/solidum/math/solvers/modal.py
"""``ModalSolver`` — problema de valores característicos ``K·φ = ω²·M·φ``
(ADR 0009 fase 1).
"""
import numpy as np

from solidum.constants import MODAL_ARPACK_TOLERANCE, MODAL_EIGENVALUE_ZERO_RTOL
from solidum.core.dof_names import TRANSLATIONS
from solidum.math.linalg import EigenSolver
from solidum.math.linalg.eigen import SemidefiniteMassError
from solidum.math.modal_response import participating_mass, participation_factors
from solidum.math.rigid_body import dof_unit_vector
from solidum.math.solvers._shared import LinearizedStateMixin, _dof_keys, _dof_labels, _log
from solidum.math.solvers.diagnostics import MasslessDofsError
from solidum.registry import SolverRegistry
from solidum.results import ModalResult


@SolverRegistry.register
class ModalSolver(LinearizedStateMixin):
    """Análisis modal: frecuencias propias y modos de vibración (ADR 0009 fase 1).

    PIPELINE_KIND = "modal" — ``run_yaml`` despacha a ``run_modal``.

    Resuelve el problema generalizado simétrico ``K·φ = ω²·M·φ`` sobre los
    DOFs libres tras imposición de Dirichlet por eliminación directa. La
    capa algebraica delega en :class:`EigenSolver` (ARPACK Lanczos con
    shift-invert centrado en ``sigma``).

    Grados de libertad sin masa (ADR 0024): masas sólo en algunos nudos o
    elementos con ``density = 0`` dejan ``M`` semidefinida. Esos grados de
    libertad se condensan estáticamente, de forma exacta, y sus componentes
    de los modos se recuperan. Hay tantas frecuencias finitas como grados
    de libertad con masa; si el problema no tiene solución, se lanza
    :class:`~solidum.math.solvers.diagnostics.MasslessDofsError`, que los
    nombra. El resultado los lista en ``condensed_dofs``.

    Sobre el estado de una etapa (ADR 0031): con ``continue_from(previous)``
    (o ``run_modal(..., previous=...)``) usa la rigidez tangente del estado
    que dejó la etapa, no la elástica: los periodos con P-Δ tras la gravedad,
    o los de la estructura agrietada tras un empuje. Un autovalor negativo
    (estado inestable en ese modo) se reporta con su signo en
    ``eigenvalues``, con frecuencia y periodo ``NaN`` y un aviso.

    Participación modal: el resultado trae, para cada traslación del modelo
    (``ux``, ``uy``, ``uz``), los factores de participación ``Γₙ = φₙᵀ·M·r``
    y las masas efectivas ``Γₙ²``, con ``r`` la traslación rígida unitaria
    (Chopra 2017, §13.1). Para otro vector ``r`` (un giro, una dirección
    oblicua), :func:`solidum.math.modal_response.participation_factors` con
    :func:`solidum.math.rigid_body.rigid_body_vector`.

    Parameters
    ----------
    assembler : Assembler
        Vinculación al modelo. Se reusan sus cachés de topología COO y
        ConstraintSet.
    n_modes : int
        Número de modos a calcular. Restricción ARPACK: ``1 ≤ n_modes < N_libre``.
    sigma : float, default 0.0
        Shift (rad²/s²) para shift-invert. ``0`` → frecuencias más bajas
        (uso por defecto en ingeniería estructural). Cualquier otro valor
        extrae los autovalores más cercanos a él.
    which : str, default "LM"
        Estrategia ARPACK. Con shift-invert el default es el adecuado;
        override solo para diagnóstico.
    tolerance : float, default 1.0e-9
        Tolerancia ARPACK sobre el residuo del autovalor.
    lumping : str, default "consistent"
        Discretización de la masa. Valores: ``"consistent"`` (default) o
        ``"lumped"`` (ADR 0009 fases 1 y 2, ambas cerradas). El esquema
        lumped delega en ``solidum.math.mass_lumping.lump_hrz`` para sólidos
        isoparamétricos y en la fórmula nodal directa de cada elemento
        para vigas/marcos.
    linear_algebra : str, default "auto"
        Reservado. ARPACK usa su factorización LU interna sin override;
        el parámetro vive aquí para coherencia con los demás solvers y
        para no romper la firma cuando se cablee a la capa algebraica
        (ADR 0003).
    """

    PIPELINE_KIND = "modal"

    def __init__(
        self,
        assembler,
        n_modes: int,
        *,
        sigma: float = 0.0,
        which: str = "LM",
        tolerance: float = MODAL_ARPACK_TOLERANCE,
        lumping: str = "consistent",
        linear_algebra: str = "auto",
    ):
        self.assembler = assembler
        self.n_modes = int(n_modes)
        self.sigma = float(sigma)
        self.which = str(which)
        self.tolerance = float(tolerance)
        self.lumping = str(lumping)
        self.linear_algebra = str(linear_algebra)

    def solve(self) -> ModalResult:
        _log.info("--- INICIANDO SOLVER MODAL ---")

        K_global = self._state_stiffness()
        self.K_used = K_global
        M_global = self.assembler.assemble_mass_matrix(lumping=self.lumping)

        K_red, M_red, T = self.assembler.reduce_pair(K_global, M_global)

        eig = EigenSolver(
            sigma=self.sigma, which=self.which, tol=self.tolerance
        )
        try:
            lambdas, phi_red = eig.solve(K_red, M_red, self.n_modes)
        except SemidefiniteMassError as exc:
            free = self.assembler.constraint_set.free_dofs(self.assembler.ndof)
            raise MasslessDofsError(
                cause=exc.cause,
                massless=exc.massless,
                massless_dofs=_dof_labels(self.assembler.domain, free[exc.massless]),
                n_mass=exc.n_mass,
                n_modes=self.n_modes,
            ) from exc
        if eig.last_massless.size:
            _log.info(
                f"  {eig.last_massless.size} grado(s) de libertad libre(s) sin masa: "
                f"condensación estática exacta (ADR 0024); "
                f"{K_red.shape[0] - eig.last_massless.size} con masa."
            )

        # Autovalores casi nulos por ruido numérico (modos de cuerpo rígido o
        # casi): cero, con T = ∞. Uno negativo de verdad (rigidez tangente
        # indefinida: el estado de la etapa es inestable en ese modo) no se
        # recorta (ADR 0031): queda con su signo en ``eigenvalues`` y con
        # frecuencia y periodo NaN.
        lam_max = float(np.max(np.abs(lambdas))) if lambdas.size else 0.0
        lambdas = np.where(
            np.abs(lambdas) < MODAL_EIGENVALUE_ZERO_RTOL * lam_max, 0.0, lambdas,
        )
        unstable = np.flatnonzero(lambdas < 0.0)
        if unstable.size:
            _log.warning(
                "ModalSolver: %d autovalor(es) negativo(s) (modos %s, λ = %s): la "
                "rigidez es indefinida y el estado es inestable en esos modos (p. ej. "
                "una carga axial por encima de la de pandeo). Se reportan con su signo "
                "en `eigenvalues`; su frecuencia y su periodo son NaN.",
                unstable.size, ", ".join(str(i + 1) for i in unstable[:6]),
                ", ".join(f"{lambdas[i]:.3e}" for i in unstable[:6]))

        with np.errstate(invalid="ignore"):
            omegas = np.where(lambdas >= 0.0, np.sqrt(np.abs(lambdas)), np.nan)
        frequencies_hz = omegas / (2.0 * np.pi)
        with np.errstate(divide="ignore", invalid="ignore"):
            periods = np.where(
                frequencies_hz > 0.0,
                1.0 / frequencies_hz,
                np.where(np.isnan(frequencies_hz), np.nan, np.inf),
            )

        modes = np.asarray(T @ phi_red)

        n_found = int(omegas.size)
        if n_found:
            _log.info(
                f"  -> {n_found} modo(s) calculado(s). "
                f"ω₁={omegas[0]:.4e} rad/s ({frequencies_hz[0]:.4e} Hz)"
            )

        # Participación en cada traslación del modelo: r = traslación rígida
        # unitaria (los apoyos se mueven con el suelo), Γₙ = φₙᵀ·M·r.
        # Masa de referencia de la razón acumulada: la que puede participar,
        # bᵀ·M_LL⁻¹·b (decisión del usuario del 2026-11-09, ADR 0009 §8).
        domain = self.assembler.domain
        participation, effective, total, participating = {}, {}, {}, {}
        for name in TRANSLATIONS:
            if not any(name in node.dofs for node in domain.nodes.values()):
                continue
            r = dof_unit_vector(domain, name, self.assembler.ndof)
            Mr = M_global @ r
            gamma = participation_factors(modes, M_global, r)
            participation[name] = gamma
            effective[name] = gamma * gamma
            total[name] = float(r @ Mr)
            participating[name] = participating_mass(M_red, T.T @ Mr)
        if effective:
            _log.info("  -> masa efectiva acumulada: " + ", ".join(
                f"{name} {100.0 * float(np.sum(m)) / participating[name]:.1f} %"
                if participating[name] > 0.0 else f"{name} —"
                for name, m in effective.items()))
        free = self.assembler.constraint_set.free_dofs(self.assembler.ndof)
        condensed = tuple(_dof_keys(domain, free[eig.last_massless]))

        # ``converged`` refleja lo que reporta ARPACK: ``False`` si se detuvo
        # con menos modos de los pedidos (``n_modes`` es entonces el número
        # realmente devuelto, no el solicitado).
        return ModalResult(
            frequencies_rad=omegas,
            frequencies_hz=frequencies_hz,
            periods=periods,
            modes=modes,
            n_modes=n_found,
            converged=bool(eig.last_converged),
            participation_factors=participation,
            effective_masses=effective,
            condensed_dofs=condensed,
            total_mass=total,
            participating_mass=participating,
            eigenvalues=np.asarray(lambdas, dtype=float),
        )
