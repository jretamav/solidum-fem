# solidum_fem/solidum/materials/thermal_conduction.py
"""Conducción de calor de Fourier (Etapa 8, spec ``docs/specs/ThermalConduction.md``).

Ley constitutiva ``q = -k·∇T`` con ``k`` tensor de conductividad. Lineal y
sin variables internas: el análogo térmico de ``Elastic2D``/``Elastic3D``.
"""
import numpy as np

from solidum.core.thermal_material import ThermalMaterial
from solidum.registry import ThermalMaterialRegistry


@ThermalMaterialRegistry.register
class ThermalConduction(ThermalMaterial):
    """Material de conducción de calor con conductividad tensorial.

    El contrato guarda siempre ``k`` como tensor; el caso isótropo es el
    particular ``k = k·I`` y el constructor lo expande. Evita romper el
    contrato cuando entren materiales ortótropos (madera, composites, roca
    estratificada), habituales en el dominio del usuario.

    Parameters
    ----------
    k : float or array-like
        Conductividad térmica [W/(m·K)]. Escalar ⇒ isótropo, se expande a
        ``k·I``. Matriz ``(d, d)`` ⇒ ortótropo/anisótropo; debe ser
        simétrica (relaciones recíprocas de Onsager) y definida positiva
        (segunda ley: no puede fluir calor hacia lo más caliente).
    c : float, optional
        Calor específico [J/(kg·K)]. Opcional al construir; obligatorio
        sólo en análisis transitorio.
    density : float, optional
        Densidad [kg/m³]. Opcional al construir; obligatoria sólo en
        análisis transitorio (criterio del ADR 0008).
    dim : int, default 2
        Dimensión del problema (2 ó 3). Determina ``FLUX_DIM`` y el tamaño
        al que se expande un ``k`` escalar. Se ignora si ``k`` llega ya
        como matriz, cuyo tamaño manda.

    Notes
    -----
    El signo negativo de ``q = -k·∇T`` es física, no convención: el calor
    fluye de mayor a menor temperatura. De ahí se sigue que ``k`` deba ser
    definida positiva — si no lo fuera, existiría una dirección en la que
    el calor fluiría espontáneamente hacia la región más caliente.
    """

    PRIMARY_STATE_VAR = None
    IS_SYMMETRIC = True

    # Tolerancia relativa para la comprobación de simetría de k. Escalada
    # con la magnitud del propio tensor para que el criterio sea
    # adimensional y no dependa de las unidades del usuario (ADR 0006).
    _SYMMETRY_RTOL = 1.0e-12

    def __init__(
        self,
        k,
        c: float | None = None,
        density: float | None = None,
        dim: int = 2,
    ):
        k_arr = np.atleast_1d(np.asarray(k, dtype=float))

        if k_arr.ndim == 1 and k_arr.size == 1:
            # Escalar ⇒ isótropo. Se expande a k·I con la dimensión pedida.
            if dim not in (2, 3):
                raise ValueError(
                    f"ThermalConduction: dim={dim} inválido; esperado 2 ó 3."
                )
            k_scalar = float(k_arr[0])
            if k_scalar <= 0.0:
                raise ValueError(
                    f"ThermalConduction: k={k_scalar} debe ser estrictamente "
                    "positivo (el calor fluye de mayor a menor temperatura)."
                )
            self.conductivity = k_scalar * np.eye(dim)
            self.is_isotropic = True
        elif k_arr.ndim == 2:
            if k_arr.shape[0] != k_arr.shape[1]:
                raise ValueError(
                    f"ThermalConduction: k con forma {k_arr.shape} debe ser "
                    "cuadrada (2×2 en 2D, 3×3 en 3D)."
                )
            if k_arr.shape[0] not in (2, 3):
                raise ValueError(
                    f"ThermalConduction: k de tamaño {k_arr.shape[0]} no "
                    "soportado; esperado 2×2 ó 3×3."
                )
            scale = np.max(np.abs(k_arr))
            if not np.allclose(
                k_arr, k_arr.T, rtol=0.0, atol=self._SYMMETRY_RTOL * scale
            ):
                raise ValueError(
                    "ThermalConduction: el tensor de conductividad k debe ser "
                    "simétrico (relaciones recíprocas de Onsager). "
                    f"Recibido:\n{k_arr}"
                )
            eigvals = np.linalg.eigvalsh(k_arr)
            if np.min(eigvals) <= 0.0:
                raise ValueError(
                    "ThermalConduction: el tensor de conductividad k debe ser "
                    "definido positivo (segunda ley de la termodinámica). "
                    f"Autovalores: {eigvals}."
                )
            self.conductivity = k_arr.copy()
            self.is_isotropic = bool(
                np.allclose(k_arr, k_arr[0, 0] * np.eye(k_arr.shape[0]))
            )
        else:
            raise ValueError(
                f"ThermalConduction: k con forma {k_arr.shape} no interpretable; "
                "esperado escalar (isótropo) o matriz cuadrada 2×2 / 3×3."
            )

        if c is not None and c <= 0.0:
            raise ValueError(
                f"ThermalConduction: c={c} debe ser estrictamente positivo."
            )
        if density is not None and density <= 0.0:
            raise ValueError(
                f"ThermalConduction: density={density} debe ser estrictamente "
                "positiva."
            )

        self.specific_heat = None if c is None else float(c)
        self.density = None if density is None else float(density)

    @property
    def FLUX_DIM(self) -> int:  # noqa: N802 — nombre de contrato, no de método
        """Dimensión del gradiente y del flujo, deducida del tensor."""
        return self.conductivity.shape[0]

    def compute_flux(self, grad_T):
        """Devuelve ``(q, k)`` con ``q = -k·∇T`` (ley de Fourier)."""
        grad = np.asarray(grad_T, dtype=float)
        if grad.shape != (self.FLUX_DIM,):
            raise ValueError(
                f"ThermalConduction: ∇T de forma {grad.shape} incompatible con "
                f"FLUX_DIM={self.FLUX_DIM}. El material se construyó para un "
                f"problema {self.FLUX_DIM}D."
            )
        return -self.conductivity @ grad, self.conductivity

    @property
    def thermal_diffusivity(self) -> float:
        """Difusividad térmica ``α = k/(ρc)`` [m²/s] del caso isótropo.

        Magnitud derivada, no parámetro: gobierna la velocidad de
        propagación del frente térmico y sirve para estimar el paso de
        tiempo característico ``Δt ~ h²/α``.

        Raises
        ------
        ValueError
            Si el material es anisótropo (la difusividad sería tensorial y
            no un escalar), o si falta ``c`` o ``density``.
        """
        if not self.is_isotropic:
            raise ValueError(
                "ThermalConduction: la difusividad escalar no está definida "
                "para conductividad anisótropa; sería un tensor. Usa los "
                "autovalores de `conductivity` para estimar los extremos."
            )
        rho_c = self.volumetric_capacity(consumer="la difusividad térmica")
        return float(self.conductivity[0, 0]) / rho_c
