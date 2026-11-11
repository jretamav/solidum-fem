"""``CohesiveDamageIsotropic`` — ley cohesiva tracción-salto con daño escalar
de modo I (material de interfaz de superficie, ADR 0033).

Daño escalar isótropo activado por tracción (Rankine), ablandamiento lineal o
exponencial gobernado por la energía de fractura ``G_f``, con rigidez de
penalización ``K_e`` (Alfaiate, Wells y Sluys 2002; Retama 2010, cap. 3).
Spec ``docs/specs/CohesiveDamageIsotropic.md``.
"""
from __future__ import annotations

import numpy as np

from solidum.core.interface_material import InterfaceMaterial
from solidum.materials.interface._softening import SofteningCurve
from solidum.registry import InterfaceMaterialRegistry


@InterfaceMaterialRegistry.register
class CohesiveDamageIsotropic(InterfaceMaterial):
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

    Es el modelo de daño isótropo con penalización de Alfaiate, Wells y Sluys
    (2002, ecs. 8-17), que Retama (2010) sigue en las ecs. 3.2-3.16; la única
    diferencia es ``H`` en lugar de ``G_F`` en el exponente, para que la
    energía disipada sea exactamente ``G_F``.

    Tangente en frame local (Modo-I) ``T_tan = α·(n⊗n)``, simétrica:

    - carga activa: ``α = dT_soft/dκ`` (tangente consistente, negativa en
      todo el ablandamiento y nula con la grieta totalmente abierta);
    - rama elástica, descarga y recarga: secante ``α = (1 − ω)·K_e``;
    - cierre (``[[u_n]] < 0``): ``t_n = K_e·[[u_n]]``, sin daño; la
      penalización impide la interpenetración (Alfaiate et al. 2002, p. 667).

    Por construcción ``T_tan[1,1] = 0`` (rigidez tangencial nula, Retama 2010,
    p. 67): la grieta desliza libre en ``s``. Modo mixto introducirá rigidez
    en esa componente.

    Para deducción física, condiciones de Kuhn-Tucker, validación energética
    y benchmarks, ver ``docs/specs/CohesiveDamageIsotropic.md``.
    """
    GEOMETRY = "surface"
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
        self.sigma_t0 = sigma_t0
        self.G_f = G_f
        self.K_e = K_e
        self.softening = softening

        # Envolvente compartida con CohesiveDamageMixedMode (_softening.py);
        # valida la forma y K_e > σ_t0²/(2·G_F): si no se cumple, la rama
        # elástica consumiría toda la energía y ω saltaría a 1 en cuanto
        # κ > κ_0, disipando más que G_F sin aviso.
        self._curve = SofteningCurve(sigma_t0, G_f, K_e, softening,
                                     shapes=self._VALID_SOFTENING,
                                     owner="CohesiveDamageIsotropic")
        self.kappa_0 = self._curve.kappa_0       # umbral del salto equivalente
        # Magnitudes de la spec de esta ley: apertura crítica (lineal) y
        # energía del ablandamiento (exponencial).
        self.w_c = 2.0 * G_f / sigma_t0 if softening == self.SOFTENING_LINEAR else None
        self.H = self._curve.H if softening == self.SOFTENING_EXPONENTIAL else None

    def compute_traction(self, jump: np.ndarray, state_vars=None):
        # Estado committed (κ histórico). Si no hay estado previo, el material
        # parte virgen en κ_0; toda apertura por debajo del umbral es elástica.
        kappa_old = self.kappa_0 if state_vars is None else state_vars.get('kappa', self.kappa_0)

        # Frame local: n = (1, 0). [[u_n]] es la componente 0 del salto.
        u_n = float(jump[0])

        if u_n < 0.0:
            # Cierre de la grieta: se recupera la relación elástica inicial,
            # sin daño, y la penalización K_e impide la interpenetración de
            # las caras (Alfaiate, Wells y Sluys 2002, p. 667: "if crack
            # closure occurs, the initial elastic constitutive relation is
            # recovered"). El daño no evoluciona en compresión.
            kappa_new = kappa_old
            traction_n = self.K_e * u_n
            stiffness_nn = self.K_e
        else:
            # Kuhn-Tucker sobre ⟨[[u_n]]⟩ (Retama 2010, ec. 3.13-3.15).
            loading = u_n > kappa_old
            kappa_new = u_n if loading else kappa_old
            T_soft, dT_soft = self._envelope(kappa_new)
            if loading and kappa_new > self.kappa_0:
                # Carga sobre la envolvente de ablandamiento: t_n = T_soft(κ)
                # con κ = [[u_n]], y la tangente consistente es su pendiente
                # (ec. 3.11 con ∂κ/∂[[u_n]] = 1). Es negativa en toda la rama
                # y se usa hasta la apertura total, donde vale 0: el sistema
                # local del elemento no se vuelve singular porque su K_jj
                # incluye el término del volumen.
                traction_n = T_soft
                stiffness_nn = dT_soft
            else:
                # Rama elástica (κ = κ_0) o descarga/recarga por debajo de κ:
                # secante al origen S(κ) = T_soft(κ)/κ = (1 − ω)·K_e (ec. 3.12).
                secant = T_soft / kappa_new
                traction_n = secant * u_n
                stiffness_nn = secant

        # ω físico, sólo informativo (estado y exportación): 1 − S(κ)/K_e.
        omega = 1.0 - self._envelope(kappa_new)[0] / (self.K_e * kappa_new)

        traction = np.array([traction_n, 0.0])     # t_s = 0: Modo I puro
        tangent = np.zeros((2, 2))
        tangent[0, 0] = stiffness_nn

        new_state = {'kappa': kappa_new, 'damage': omega}
        return traction, tangent, new_state

    def _envelope(self, kappa: float):
        """Envolvente ``(T_soft(κ), dT_soft/dκ)`` de la tracción normal.

        Para ``κ ≤ κ_0`` es la rama elástica ``K_e·κ`` (pendiente ``K_e``);
        para ``κ > κ_0`` la de ablandamiento, que empieza en ``σ_t0`` y
        encierra con la rama elástica un área ``G_F``. Se evalúa sin pasar
        por ``1 − ω``: con ``K_e`` de penalización, ``ω`` está a ``~κ_0/κ``
        de 1 y la resta perdería cifras.
        """
        return self._curve(kappa)
