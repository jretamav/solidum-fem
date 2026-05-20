# fenix_fem/solidum/math/solvers/central_difference.py
"""``CentralDifferenceSolver`` — integración explícita por diferencias centradas
para ``M·ü + C·u̇ + K·u = F(t)`` (ADR 0009 fase 5).

Esquema **explícito de segundo orden** sin sistema lineal a resolver en cada
paso — la única "inversión" es ``M⁻¹``, trivial cuando M es diagonal
(ADR 0009 fase 2, mass lumping). Apropiado para:

- Wave propagation en sólidos (ondas elásticas, impacto).
- Dinámica no lineal de respuesta rápida (rotura por impacto, explosión).
- Problemas donde el paso temporal Δt natural cae por debajo del paso de
  estabilidad incondicional de Newmark.

**Condición de estabilidad** (Courant-Friedrichs-Lewy, sistema no
amortiguado):

.. math:: \\Delta t \\le \\frac{2}{\\omega_\\max}

donde ``ω_max`` es la frecuencia natural máxima del sistema discreto. Para
una barra axial elástica con ``L_el = L/N`` y ``c_p = √(E/ρ)``:

.. math:: \\Delta t_\\text{crit} \\approx \\frac{L_\\text{el}}{c_p}

(condición CFL). Si se excede, la solución diverge — el solver detecta la
explosión exponencial y lanza ``RuntimeError`` con mensaje informativo.

**Esquema predictor-corrector (leapfrog Belytschko-Liu-Moran)**:

.. code-block:: text

    Inicialización:
        a_0 = M⁻¹·(F(0) − F_int(u_0) − C·v_0 − F_dir)
        v_{1/2} = v_0 + (Δt/2)·a_0

    Paso n → n+1:
        u_{n+1}   = u_n + Δt·v_{n+1/2}
        F_int_{n+1} = K·u_{n+1}    (lineal)
                    = ensamble(u_{n+1})  (no lineal)
        a_{n+1}   = M⁻¹·(F(t_{n+1}) − F_int_{n+1} − C·v_{n+1/2} − F_dir)
        v_{n+1}   = v_{n+1/2} + (Δt/2)·a_{n+1}
        v_{n+3/2} = v_{n+1} + (Δt/2)·a_{n+1}    [para el paso siguiente]

Sin Newton interno — para problemas no lineales basta con reevaluar
``F_int(u_{n+1})`` con el estado actual; no se requiere convergencia
iterativa porque la ecuación de aceleración es explícita.

Referencias:

- Belytschko, T., Liu, W. K., Moran, B. (2014). *Nonlinear Finite Elements
  for Continua and Structures*, §6.2.
- Hughes, T. J. R. (2000). *The Finite Element Method*, §9.1.
- Crisfield, M. A. (1997). *Non-Linear Finite Element Analysis of Solids
  and Structures*, vol. 2, §24.6.
"""
from __future__ import annotations

from typing import Callable

import numpy as np
import scipy.sparse as sp

from solidum.constants import LUMPED_MASS_OFF_DIAGONAL_RTOL
from solidum.math.damping import resolve_rayleigh_config
from solidum.math.solvers._shared import _log
from solidum.registry import SolverRegistry
from solidum.results import TransientResult


@SolverRegistry.register
class CentralDifferenceSolver:
    """Integración explícita por diferencias centradas (ADR 0009 fase 5).

    Parameters
    ----------
    assembler : Assembler
        Vinculación al modelo.
    t_end : float
        Tiempo final (s).
    dt : float
        Paso temporal constante (s). Debe satisfacer la condición CFL
        ``Δt < 2/ω_max``; el solver detecta divergencia a posteriori.
    nonlinear : bool, default False
        Si ``True``, recalcula ``F_int`` cada paso vía
        ``Assembler.assemble_non_linear_system`` (descarta la tangente).
        Si ``False``, usa ``K·u`` con ``K`` constante ensamblada una vez.
    rayleigh : dict or None, default None
        Amortiguamiento Rayleigh ``C = α·M + β·K``. Mismo contrato que
        :class:`NewmarkSolver`. Con amortiguamiento, el ``Δt_crit`` cae
        ligeramente — responsabilidad del usuario.
    u0, u0_dot : np.ndarray or None
        Condiciones iniciales globales (shape ``(ndof,)``). Cero por
        defecto.
    F_func : callable or None
        ``F_func(t) -> np.ndarray`` evaluado en cada paso. ``None`` ⇒
        vibración libre.
    lumping : str, default "lumped"
        Masa lumped por defecto — el esquema explícito **pierde su razón
        de ser** sin M diagonal (cada paso requeriría una factorización).
        Si se pasa ``"consistent"`` el solver lanza ``ValueError``: usar
        Newmark implícito si la formulación requiere masa consistente.
    divergence_threshold : float, default 1.0e8
        Factor multiplicativo sobre la escala inicial de ``u`` por encima
        del cual se considera divergencia. Determina el momento en que se
        aborta con ``RuntimeError`` (típicamente porque ``Δt > Δt_crit``).

    Notes
    -----
    El solver no calcula ``Δt_crit`` automáticamente. Una estimación
    canónica por elemento es ``Δt_crit_el = L_el / c_p`` con
    ``c_p = √(E/ρ)`` (barras y sólidos 2D); el ``Δt_crit`` del sistema es
    el mínimo entre elementos. Para Frame3D con lumped bloque-diagonal,
    cada paso requiere invertir bloques 6×6 nodales (no implementado en
    esta versión inicial): el solver acepta solamente lumped estrictamente
    diagonal — Frame3D oblicuo se rechaza con ``ValueError``.
    """

    PIPELINE_KIND = "transient"

    def __init__(
        self,
        assembler,
        t_end: float,
        dt: float,
        *,
        nonlinear: bool = False,
        rayleigh: dict | None = None,
        u0: np.ndarray | None = None,
        u0_dot: np.ndarray | None = None,
        F_func: Callable[[float], np.ndarray] | None = None,
        lumping: str = "lumped",
        divergence_threshold: float = 1.0e8,
        linear_algebra: str = "auto",  # noqa: ARG002 — aceptado por compatibilidad de firma
    ):
        if dt <= 0.0:
            raise ValueError(f"CentralDifferenceSolver: dt={dt} debe ser positivo.")
        if t_end <= 0.0:
            raise ValueError(
                f"CentralDifferenceSolver: t_end={t_end} debe ser positivo."
            )
        if lumping != "lumped":
            raise ValueError(
                f"CentralDifferenceSolver: lumping='{lumping}' no soportado. "
                f"El esquema explícito requiere M diagonal — usa "
                f"lumping='lumped' (ADR 0009 fase 2). Si la formulación "
                f"requiere masa consistente, usa NewmarkSolver implícito."
            )
        self.assembler = assembler
        self.t_end = float(t_end)
        self.dt = float(dt)
        self.nonlinear = bool(nonlinear)
        self.rayleigh_cfg = rayleigh
        self.u0 = u0
        self.u0_dot = u0_dot
        self.F_func = F_func
        self.lumping = lumping
        self.divergence_threshold = float(divergence_threshold)

    @staticmethod
    def _invert_diagonal_mass(M_red: sp.spmatrix) -> np.ndarray:
        """Devuelve el vector ``1/diag(M_red)`` si M_red es diagonal pura.

        El esquema explícito sólo es eficiente con M estrictamente diagonal
        (ADR 0009 fase 2 garantiza esto para todos los elementos excepto
        Frame3D oblicuo). Si M_red tiene off-diagonals no nulos, el coste
        de invertir un bloque 6×6 nodal cada paso anula la ventaja del
        explícito frente a Newmark — se rechaza con ``ValueError``.
        """
        M_dense = M_red.toarray() if sp.issparse(M_red) else np.asarray(M_red)
        diag = np.diag(M_dense)
        off = M_dense - np.diag(diag)
        off_max = float(np.max(np.abs(off))) if off.size else 0.0
        diag_max = float(np.max(np.abs(diag))) if diag.size else 1.0
        if off_max > LUMPED_MASS_OFF_DIAGONAL_RTOL * max(diag_max, 1.0):
            raise ValueError(
                "CentralDifferenceSolver: M_lumped no es estrictamente "
                "diagonal tras reducir Dirichlet — probablemente Frame3D "
                "con eje oblicuo a globales (ADR 0009 fase 2 docs). "
                "Alinea el eje del frame con un eje global o usa "
                f"NewmarkSolver. Off-diagonal max = {off_max:.3e} vs "
                f"diag max = {diag_max:.3e}."
            )
        if np.min(diag) <= 0.0:
            raise ValueError(
                "CentralDifferenceSolver: M_lumped tiene entradas no "
                f"positivas en la diagonal (min={np.min(diag):.3e}). "
                "Material sin density o densidad negativa."
            )
        return 1.0 / diag

    def solve(self) -> TransientResult:
        _log.info("--- INICIANDO SOLVER CENTRAL DIFFERENCES ---")

        # Ensamblar K (lineal) y M (lumped por construcción de este solver).
        self.assembler.assemble_system()
        K = self.assembler.K_global
        M = self.assembler.assemble_mass_matrix(lumping=self.lumping)

        # Amortiguamiento Rayleigh.
        alpha_r, beta_r = resolve_rayleigh_config(
            self.rayleigh_cfg, source=type(self).__name__,
        )

        # Reducción por Dirichlet.
        cs = self.assembler.constraint_set
        T, g = cs.build(self.assembler.ndof)
        free_dofs = cs.free_dofs(self.assembler.ndof)

        K_red = (T.T @ K @ T).tocsr()
        M_red = (T.T @ M @ T).tocsr()
        C_red = alpha_r * M_red + beta_r * K_red
        F_dir = T.T @ (K @ g)

        # M⁻¹ diagonal — vector denso, multiplicación trivial.
        M_inv_diag = self._invert_diagonal_mass(M_red)

        ndof = self.assembler.ndof
        n_free = K_red.shape[0]

        # Condiciones iniciales.
        u0_global = (np.zeros(ndof) if self.u0 is None
                      else np.asarray(self.u0, dtype=float).reshape(ndof))
        u0_dot_global = (np.zeros(ndof) if self.u0_dot is None
                          else np.asarray(self.u0_dot, dtype=float).reshape(ndof))
        u_free = u0_global[free_dofs].copy()
        v_free = u0_dot_global[free_dofs].copy()

        # F_int(u_0) y a_0.
        F0_global = (np.zeros(ndof) if self.F_func is None
                      else np.asarray(self.F_func(0.0), dtype=float).reshape(ndof))
        F0_red = T.T @ F0_global

        F_int_red = self._compute_internal_forces_reduced(u_free, T, g, K_red, F_dir)

        a_free = M_inv_diag * (F0_red - F_int_red - C_red @ v_free - F_dir)

        # Inicialización leapfrog: v_{1/2} = v_0 + (Δt/2)·a_0.
        dt = self.dt
        half_dt = 0.5 * dt
        v_half = v_free + half_dt * a_free

        # Historiales.
        n_steps = int(np.ceil(self.t_end / dt))
        t_history = np.linspace(0.0, n_steps * dt, n_steps + 1)
        u_history = np.zeros((ndof, n_steps + 1))
        udot_history = np.zeros((ndof, n_steps + 1))
        uddot_history = np.zeros((ndof, n_steps + 1))
        u_history[:, 0] = T @ u_free + g
        udot_history[:, 0] = T @ v_free
        uddot_history[:, 0] = T @ a_free

        # Escala inicial para detección de divergencia.
        initial_scale = max(
            float(np.linalg.norm(u_free)),
            float(np.linalg.norm(F0_red)) / max(float(np.max(np.diag(M_red.toarray()))), 1.0),
            1.0,
        )
        threshold = self.divergence_threshold * initial_scale

        # Bucle temporal.
        for step in range(n_steps):
            t_next = t_history[step + 1]

            # u_{n+1} = u_n + Δt·v_{n+1/2}.
            u_free = u_free + dt * v_half

            # F_int(u_{n+1}).
            F_int_red = self._compute_internal_forces_reduced(u_free, T, g, K_red, F_dir)

            # F(t_{n+1}).
            F_global = (np.zeros(ndof) if self.F_func is None
                         else np.asarray(self.F_func(t_next), dtype=float).reshape(ndof))
            F_red = T.T @ F_global

            # a_{n+1} = M⁻¹·(F − F_int − C·v_{n+1/2} − F_dir).
            a_free = M_inv_diag * (F_red - F_int_red - C_red @ v_half - F_dir)

            # v_{n+1} (output) = v_{n+1/2} + (Δt/2)·a_{n+1}.
            v_free = v_half + half_dt * a_free
            # v_{n+3/2} (próximo paso) = v_{n+1} + (Δt/2)·a_{n+1}.
            v_half = v_free + half_dt * a_free

            # Detección de divergencia.
            u_norm = float(np.linalg.norm(u_free))
            if not np.isfinite(u_norm) or u_norm > threshold:
                raise RuntimeError(
                    f"CentralDifferenceSolver: divergencia detectada en paso "
                    f"{step + 1}/{n_steps} (t={t_next:.3e}). "
                    f"|u| = {u_norm:.3e} > {threshold:.3e}. Probable "
                    f"violación de la condición CFL Δt < 2/ω_max. "
                    f"Reduce dt (actual {dt:.3e}) o usa NewmarkSolver "
                    f"implícito (incondicionalmente estable)."
                )

            # Volcado a historial global.
            u_history[:, step + 1] = T @ u_free + g
            udot_history[:, step + 1] = T @ v_free
            uddot_history[:, step + 1] = T @ a_free

        _log.info(
            f"  -> {n_steps} pasos completados. "
            f"Rayleigh: α={alpha_r:.4e}, β={beta_r:.4e}. "
            f"Modo: {'no lineal' if self.nonlinear else 'lineal'}."
        )

        return TransientResult(
            t_history=t_history,
            u_history=u_history,
            udot_history=udot_history,
            uddot_history=uddot_history,
            n_steps=n_steps,
            alpha_rayleigh=alpha_r,
            beta_rayleigh=beta_r,
        )

    def _compute_internal_forces_reduced(
        self,
        u_free: np.ndarray,
        T,
        g: np.ndarray,
        K_red: sp.spmatrix,
        F_dir: np.ndarray,
    ) -> np.ndarray:
        """``F_int`` proyectado a DOFs libres, **excluyendo la contribución del apoyo**.

        Ambas ramas devuelven ``F_int`` sin el término que mueve a los DOFs
        prescritos (``F_dir = T.T·K·g``). El cómputo de la aceleración resta
        ``-F_dir`` una vez fuera de este método, de modo que la rama lineal
        y la no lineal pasan por el mismo balance.

        - Lineal: ``F_int(u_global) = K·u_global = K·T·u_free + K·g``;
          proyectado a libres queda ``K_red·u_free + F_dir``. Aquí devolvemos
          sólo ``K_red·u_free`` (la contribución del apoyo viene de ``-F_dir``
          en el llamador).
        - No lineal: ``Assembler.assemble_non_linear_system(u_global)`` recalcula
          ``F_int`` sobre el estado actual e incluye **implícitamente** la
          contribución del apoyo (porque ``u_global = T·u_free + g``). Restamos
          ``F_dir`` aquí dentro para alinearnos con la rama lineal — si no se
          restara, el llamador contaría el efecto del apoyo dos veces y la
          dinámica con apoyos prescritos no nulos saldría corrompida.
        """
        if not self.nonlinear:
            return K_red @ u_free
        # Reconstruir u_global con apoyos prescritos y recalcular F_int.
        u_global = T @ u_free + g
        _K_t, F_int_global = self.assembler.assemble_non_linear_system(u_global)
        return T.T @ F_int_global - F_dir
