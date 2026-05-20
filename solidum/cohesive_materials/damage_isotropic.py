# solidum_fem/solidum/cohesive_materials/damage_isotropic.py
"""``CohesiveDamageIsotropic`` — material cohesivo traction-jump con daño escalar Modo-I.

Traducción a código de Retama (2010), Cap. 3: daño escalar isótropo activado
por tracción (Rankine), softening lineal o exponencial gobernado por la
energía de fractura ``G_f``. Ver ``docs/specs/CohesiveDamageIsotropic.md``.
"""
from __future__ import annotations

import numpy as np

from solidum.constants import DAMAGE_MAX
from solidum.core.cohesive_material import CohesiveMaterial
from solidum.registry import CohesiveMaterialRegistry


@CohesiveMaterialRegistry.register
class CohesiveDamageIsotropic(CohesiveMaterial):
    """
    Material cohesivo de daño isótropo escalar, Modo-I, traction-jump.

    Opera sobre el salto ``[[u]] ∈ ℝ²`` en el frame local ``(n, s)`` de ``Γ_d``
    y devuelve la tracción ``t ∈ ℝ²`` en los mismos ejes. La componente normal
    rige activación, evolución del daño y energía disipada; la tangencial
    desliza libremente sin disipación (Modo-I puro). El modo mixto I–II queda
    diferido (fase G del ADR 0010).

    Parameters
    ----------
    sigma_t0 : float
        Resistencia a tracción ``σ_{t0}`` [Pa]. Umbral por encima del cual la
        tracción normal inicia el daño.
    G_f : float
        Energía de fractura ``G_F`` [N/m]. Área bajo la curva ``t–[[u_n]]``.
    K_e : float
        Rigidez elástica del salto en dirección normal ``K_e`` [Pa/m]. Penalty:
        debe ser suficientemente grande para que ``κ_0 = σ_{t0}/K_e ≪`` escala
        del problema, pero no tanto que degrade el condicionamiento. Guía
        práctica: ``K_e ≈ 10·E_bulk/ℓ_c``. Sin default automático (ADR 0010
        separa contratos bulk / cohesivo).
    softening : str
        Forma de la curva: ``'linear'`` o ``'exponential'``.

    Notes
    -----
    Forma cerrada de ``ω(κ)`` para softening lineal (κ ∈ [κ_0, w_c])::

        ω(κ) = 1 − σ_{t0}·(w_c − κ) / [K_e·κ·(w_c − κ_0)]
        w_c = 2·G_F / σ_{t0}                                    (apertura crítica)
        dω/dκ = σ_{t0}·w_c / [K_e·κ²·(w_c − κ_0)]

    Forma cerrada para softening exponencial (κ > κ_0)::

        T_soft(κ) = σ_{t0}·exp[−σ_{t0}·(κ − κ_0) / H],   H = G_F − σ_{t0}·κ_0/2
        ω(κ) = 1 − T_soft(κ) / (K_e·κ)
        dω/dκ = (1 − ω)·(1/κ + σ_{t0}/H)

    Tangente algorítmica en frame local (Modo-I) ``T_tan = α·(n⊗n)``, simétrica.
    Por construcción ``T_tan[1,1] = 0`` (rigidez tangencial nula): la grieta
    desliza libre en ``s``. Modo mixto introducirá rigidez en esa componente.

    Para deducción física, condiciones de Kuhn-Tucker, validación energética
    y benchmarks, ver ``docs/specs/CohesiveDamageIsotropic.md``.
    """
    JUMP_DIM = 2
    PRIMARY_STATE_VAR = 'damage'
    IS_SYMMETRIC = True

    SOFTENING_LINEAR = 'linear'
    SOFTENING_EXPONENTIAL = 'exponential'
    _VALID_SOFTENING = (SOFTENING_LINEAR, SOFTENING_EXPONENTIAL)

    def __init__(self, sigma_t0: float, G_f: float, K_e: float, softening: str):
        if sigma_t0 <= 0.0:
            raise ValueError(
                f"CohesiveDamageIsotropic: sigma_t0 debe ser > 0 (recibido {sigma_t0})."
            )
        if G_f <= 0.0:
            raise ValueError(
                f"CohesiveDamageIsotropic: G_f debe ser > 0 (recibido {G_f})."
            )
        if K_e <= 0.0:
            raise ValueError(
                f"CohesiveDamageIsotropic: K_e debe ser > 0 (recibido {K_e})."
            )
        if softening not in self._VALID_SOFTENING:
            raise ValueError(
                f"CohesiveDamageIsotropic: softening debe ser uno de "
                f"{self._VALID_SOFTENING}, recibido {softening!r}."
            )

        self.sigma_t0 = sigma_t0
        self.G_f = G_f
        self.K_e = K_e
        self.softening = softening

        # Magnitudes derivadas
        self.kappa_0 = sigma_t0 / K_e           # umbral del salto equivalente

        if softening == self.SOFTENING_LINEAR:
            self.w_c = 2.0 * G_f / sigma_t0     # apertura crítica
            self.H = None
        else:  # exponential
            self.H = G_f - 0.5 * sigma_t0 * self.kappa_0
            if self.H <= 0.0:
                raise ValueError(
                    f"CohesiveDamageIsotropic: softening exponencial requiere "
                    f"K_e > σ_t0² / (2·G_F)  (H = G_F − σ_t0·κ_0/2 > 0). "
                    f"Con los parámetros recibidos, H = {self.H:.3e} ≤ 0."
                )
            self.w_c = None

    def compute_traction(self, jump: np.ndarray, state_vars=None):
        # Estado committed (κ histórico). Si no hay estado previo, el material
        # parte virgen en κ_0; toda apertura por debajo del umbral es elástica.
        kappa_old = self.kappa_0 if state_vars is None else state_vars.get('kappa', self.kappa_0)

        # Frame local: n = (1, 0). [[u_n]] es la componente 0 del salto.
        u_n = float(jump[0])
        u_n_pos = u_n if u_n > 0.0 else 0.0  # ⟨[[u_n]]⟩ (McAuley)

        # Régimen carga / descarga (Kuhn-Tucker)
        if u_n_pos > kappa_old:
            kappa_new = u_n_pos
            loading = True
        else:
            kappa_new = kappa_old
            loading = False

        # ω físico (puede llegar a 1.0 exacto en lineal con κ ≥ w_c, o
        # asintóticamente en exponencial); cap_tangent señala al llamador que
        # debe imponer la rigidez residual numérica en la tangente para evitar
        # singularidad del Newton, pero la tracción se calcula con ω real.
        omega, domega_dkappa, cap_tangent = self._damage(kappa_new)

        # Tracción: t_n = (1 − ω)·K_e·[[u_n]]; t_s = 0 en Modo-I.
        # ω = 1 ⇒ t_n = 0 exactamente (grieta totalmente abierta, físico).
        traction = np.array([(1.0 - omega) * self.K_e * u_n, 0.0])

        # Tangente algorítmica en frame local (rank-1 sobre n⊗n). Cap residual
        # ``(1 − DAMAGE_MAX)·K_e`` sólo en la *tangente*, para mantener el
        # sistema lineal de Newton no singular cuando ω → 1. No contamina la
        # tracción ni la energía disipada (ver caveat numérico §12 de la spec).
        # En carga activa fuera de la zona capada se usa la consistente:
        # ``K_e·[(1 − ω) − u_n·dω/dκ]`` con dκ/d[[u_n]] = 1 (κ_new = u_n > 0).
        consistent = loading and (omega > 0.0) and (not cap_tangent)
        if consistent:
            stiffness_nn = self.K_e * ((1.0 - omega) - u_n * domega_dkappa)
        else:
            stiffness_factor = max(1.0 - omega, 1.0 - DAMAGE_MAX)
            stiffness_nn = stiffness_factor * self.K_e

        tangent = np.zeros((2, 2))
        tangent[0, 0] = stiffness_nn

        new_state = {'kappa': kappa_new, 'damage': omega}
        return traction, tangent, new_state

    def _damage(self, kappa: float):
        """Devuelve ``(ω(κ), dω/dκ, cap_tangent)`` para el ``κ`` corriente.

        ``ω`` es el valor *físico* (sin truncar por ``DAMAGE_MAX``). El flag
        ``cap_tangent`` indica que ``ω ≥ DAMAGE_MAX`` (zona numéricamente
        singular en rigidez) y el llamador debe usar la rigidez residual.
        ``dω/dκ`` es respecto al historial; el llamador la compone con
        ``∂κ/∂[[u_n]] = 1`` (carga activa) al construir la tangente.
        """
        if kappa <= self.kappa_0:
            return 0.0, 0.0, False
        if self.softening == self.SOFTENING_LINEAR:
            return self._damage_linear(kappa)
        return self._damage_exponential(kappa)

    def _damage_linear(self, kappa: float):
        # Saturación geométrica: por encima de la apertura crítica la grieta
        # está totalmente abierta (ω = 1.0 exacto, sin residuo).
        if kappa >= self.w_c:
            return 1.0, 0.0, True

        s = self.sigma_t0
        wc = self.w_c
        k0 = self.kappa_0
        Ke = self.K_e
        omega = 1.0 - s * (wc - kappa) / (Ke * kappa * (wc - k0))
        domega = s * wc / (Ke * kappa * kappa * (wc - k0))
        return omega, domega, omega >= DAMAGE_MAX

    def _damage_exponential(self, kappa: float):
        s = self.sigma_t0
        k0 = self.kappa_0
        Ke = self.K_e
        H = self.H
        T_soft = s * np.exp(-s * (kappa - k0) / H)
        omega = 1.0 - T_soft / (Ke * kappa)
        domega = (1.0 - omega) * (1.0 / kappa + s / H)
        return omega, domega, omega >= DAMAGE_MAX
