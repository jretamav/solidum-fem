"""``ContactPenalty`` — contacto unilateral por penalización con fricción de
Coulomb opcional (material de interfaz de superficie, ADR 0033 y 0043).

Spec ``docs/specs/ContactPenalty.md``.

Componentes en los ejes locales de la interfaz, normal primero: ``(g_n, g_s)``
en 2D y ``(g_n, g_s1, g_s2)`` en 3D. ``g_n = g_0 + Δu_n`` es la **holgura
normal actual** que calcula el elemento (``g_n > 0`` separación); las
tangenciales, el desplazamiento relativo desde la referencia. ``t`` es la
tracción sobre la cara A (normal exterior ``n``): en compresión
``t = −p_N·n``, y ``t_s`` tiene el sentido del deslizamiento de B respecto
de A.

- **Normal**: ``t_n = K_n·min(g_n, 0)``, presión ``p_N = −t_n ≥ 0``.
- **Tangencial**: sin fricción (``mu = 0``): ``t_s = 0``. Con fricción:
  ``t_s = K_t·(g_s − g_s^p)``, ``‖t_s‖ ≤ μ·p_N``, retorno radial (solución
  cerrada de Euler implícito; el límite no depende del deslizamiento).
- **Abierta** (``g_n > 0``): ``t = 0`` y ``g_s^p = g_s`` (el punto de
  adherencia se reinicia).
- **Holgura nula** (``g_n = 0``): rama cerrada, para que la primera
  iteración tenga rigidez (spec §11).
"""
from __future__ import annotations

import math

import numpy as np

from solidum.core.interface_material import InterfaceMaterial
from solidum.registry import InterfaceMaterialRegistry

OPEN, STICK, SLIP = 0, 1, 2


@InterfaceMaterialRegistry.register
class ContactPenalty(InterfaceMaterial):
    """Contacto unilateral por penalización, con fricción de Coulomb opcional.

    Parameters
    ----------
    K_n : float
        Rigidez de penalización normal [Pa/m]. Obligatoria, sin valor por
        omisión (decisión del usuario del 2026-11-26). Guía: ``K_n ≈ 10·E/h``
        con ``E`` y ``h`` del sólido vecino; la interpenetración es
        ``p_N/K_n``.
    mu : float, optional
        Coeficiente de fricción de Coulomb (≥ 0). Sin él, contacto sin
        fricción (decisión del usuario del 2026-11-26).
    K_t : float, optional
        Rigidez de adherencia tangencial [Pa/m]. Obligatoria si ``mu > 0`` y
        rechazada si ``mu = 0``. Guía: empezar con ``K_t = K_n``; el
        deslizamiento previo al umbral es ``μ·p_N/K_t``.
    """

    GEOMETRY = "surface"
    JUMP_DIM = None
    PRIMARY_STATE_VAR = "contact_pressure"
    # Por instancia: True sólo sin fricción (spec §10). La de la clase es la
    # cota conservadora.
    IS_SYMMETRIC = False
    IS_UNILATERAL = True
    ACCEPTS_INITIAL_GAP = True
    SECANT_UNLOADING = False
    STATE_SCHEMA = {
        "slip": (None,),             # g_s^p, 1 (2D) o 2 (3D) componentes
        "slip_accumulated": (),      # ∫‖ġ_s^p‖ con la superficie cerrada
        "contact_pressure": (),      # p_N
        "status": (),                # 0 abierta, 1 cerrada sin deslizar, 2 deslizando
    }
    # El régimen es un código: el VTK por celda no lo promedia.
    STATE_CODES = ("status",)

    def __init__(self, K_n: float, mu: float = 0.0, K_t: float | None = None):
        name = type(self).__name__
        if not (isinstance(K_n, (int, float)) and math.isfinite(K_n) and K_n > 0.0):
            raise ValueError(f"{name}: K_n debe ser un número finito > 0 (recibido {K_n!r}).")
        if not (isinstance(mu, (int, float)) and math.isfinite(mu) and mu >= 0.0):
            raise ValueError(f"{name}: mu debe ser un número finito ≥ 0 (recibido {mu!r}).")
        if mu > 0.0:
            if K_t is None:
                raise ValueError(
                    f"{name}: con fricción (mu = {mu}) hace falta K_t, la rigidez de "
                    f"adherencia tangencial [Pa/m].")
            if not (isinstance(K_t, (int, float)) and math.isfinite(K_t) and K_t > 0.0):
                raise ValueError(f"{name}: K_t debe ser un número finito > 0 (recibido {K_t!r}).")
        elif K_t is not None:
            raise ValueError(
                f"{name}: K_t sólo interviene con fricción; sin mu (o con mu = 0) el "
                f"contacto no tiene rigidez tangencial.")
        self.K_n = float(K_n)
        self.mu = float(mu)
        self.K_t = None if K_t is None else float(K_t)
        self.IS_SYMMETRIC = self.mu == 0.0

    # -- ley -------------------------------------------------------------------

    def compute_traction(self, jump, state_vars=None):
        g = np.asarray(jump, dtype=np.float64).reshape(-1)
        n = g.size
        if n not in (2, 3):
            raise ValueError(
                f"{type(self).__name__}: el salto debe tener 2 (2D) o 3 (3D) "
                f"componentes; tiene {n}.")
        m = n - 1
        st = state_vars or {}
        slip_old = np.asarray(st.get("slip", np.zeros(m)), dtype=np.float64).reshape(-1)
        if slip_old.size != m:
            raise ValueError(
                f"{type(self).__name__}: el estado tiene {slip_old.size} componentes de "
                f"deslizamiento y el salto {m} tangenciales.")
        acc_old = float(st.get("slip_accumulated", 0.0))

        g_n, g_s = float(g[0]), g[1:]
        traction = np.zeros(n)
        tangent = np.zeros((n, n))

        if g_n > 0.0:
            # Abierta: sin tracción ni rigidez; todo el desplazamiento
            # tangencial es deslizamiento (el punto de adherencia se reinicia).
            new_state = {"slip": g_s.copy(), "slip_accumulated": acc_old,
                         "contact_pressure": 0.0, "status": float(OPEN)}
            return traction, tangent, new_state

        # Cerrada (g_n ≤ 0; en g_n = 0, rama cerrada: rigidez de arranque).
        p = max(-self.K_n * g_n, 0.0)
        traction[0] = -p
        tangent[0, 0] = self.K_n

        if self.mu == 0.0:
            # Sin fricción: desliza libre, t_s = 0.
            slip_new = g_s.copy()
            status = STICK
        else:
            bound = self.mu * p
            t_trial = self.K_t * (g_s - slip_old)
            norm_trial = float(np.linalg.norm(t_trial))
            if norm_trial <= bound:
                # Adherencia (incluye t_trial = 0 con p = 0: rigidez de arranque).
                traction[1:] = t_trial
                tangent[1:, 1:] = self.K_t * np.eye(m)
                slip_new = slip_old.copy()
                status = STICK
            else:
                # Deslizamiento: retorno radial al círculo ‖t_s‖ = μ·p_N.
                direction = t_trial / norm_trial
                traction[1:] = bound * direction
                slip_new = slip_old + ((norm_trial - bound) / self.K_t) * direction
                # ∂t_s/∂g_n = −μ·K_n·m (el límite crece con la compresión).
                tangent[1:, 0] = -self.mu * self.K_n * direction
                # ∂t_s/∂g_s = (μ·p_N·K_t/‖t^tr‖)(I − m⊗m); nula en 2D.
                tangent[1:, 1:] = (bound * self.K_t / norm_trial) * (
                    np.eye(m) - np.outer(direction, direction))
                status = SLIP

        acc_new = acc_old + float(np.linalg.norm(slip_new - slip_old))
        new_state = {"slip": slip_new, "slip_accumulated": acc_new,
                     "contact_pressure": p, "status": float(status)}
        return traction, tangent, new_state

    def stored_energy(self, jump, state_vars=None) -> float:
        """``ψ = ½·K_n·min(g_n, 0)² + ½·K_t·‖g_s − g_s^p‖²`` (spec §5), con
        ``g_s^p`` el del estado del paso (el que devolvió
        :meth:`compute_traction` para este salto)."""
        g = np.asarray(jump, dtype=np.float64).reshape(-1)
        closure = min(float(g[0]), 0.0)
        psi = 0.5 * self.K_n * closure * closure
        if self.mu > 0.0 and g[0] <= 0.0:
            st = state_vars or {}
            slip = np.asarray(st.get("slip", np.zeros(g.size - 1)), dtype=np.float64)
            elastic = g[1:] - slip
            psi += 0.5 * self.K_t * float(elastic @ elastic)
        return psi
