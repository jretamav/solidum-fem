"""Piezas comunes de la constitutiva de los marcos (ADR 0021).

Los marcos Euler-Bernoulli (``Frame2DEuler``, ``Frame2DEulerCorot``,
``Frame3D``) aceptan **exactamente uno** de dos caminos, decidido al
construir:

- ``material`` + propiedades geométricas: la forma cerrada elástica de
  siempre. Sólo admite materiales con ``IS_LINEAR_ELASTIC = True``; con un
  material inelástico escalaría toda la rigidez local con la ``E_t`` de la
  deformación axial (flector nulo con fluencia axial, de signo cambiado con
  ablandamiento: auditoría 2026-09-29). La misma guarda protege a
  ``Frame2DTimoshenko`` hasta la etapa 2.
- ``section`` (+ ``n_ip``): elemento en desplazamientos que integra la ley de
  sección en ``n_ip`` puntos de Gauss-Legendre a lo largo del elemento
  (spec ``FiberSection`` §5). El estado de la sección en el punto ``g`` vive
  en ``state.vars[g]`` y sus fuerzas internas ``s_g`` en ``state.stresses[g]``.

:class:`FrameSectionMixin` reúne lo que comparten: la validación de los
argumentos, la guarda de material elástico lineal, el estado por punto de
Gauss, la integración ``K = Σ w_g Bᵀ k_s B``, ``F = Σ w_g Bᵀ s`` y la
evaluación sin huella en el trial que usan las fuerzas internas. Cada marco
aporta su cinemática (las matrices ``B_g``) y su transformación a globales.
"""
from __future__ import annotations

from typing import ClassVar, Tuple

import numpy as np

from solidum.core.element_state import ElementState
from solidum.core.section import Section

# Puntos de Gauss-Legendre a lo largo del elemento (spec §5): 3 por defecto;
# 2 es el mínimo (integra exactamente la rigidez elástica, cuyo integrando es
# cuadrático); con 1 el modo de flexión antisimétrica no tiene energía.
DEFAULT_N_IP = 3
MIN_N_IP = 2


def gauss_legendre_unit(n: int) -> Tuple[np.ndarray, np.ndarray]:
    """Puntos ``ξ ∈ (0, 1)`` y pesos (suma 1) de Gauss-Legendre de ``n``
    puntos sobre el intervalo unidad: ``∫₀ᴸ f dx ≈ L·Σ w_g f(ξ_g L)``."""
    xi, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (xi + 1.0), 0.5 * w


def hermite_second_derivatives(xi: float, L: float) -> np.ndarray:
    """``[N₁'', N₂'', N₃'', N₄'']`` de las funciones cúbicas de Hermite en
    ``x = ξ·L`` (spec §5), para ``v = N₁v₁ + N₂θ₁ + N₃v₂ + N₄θ₂``::

        N₁'' = (12ξ − 6)/L²,  N₂'' = (6ξ − 4)/L,
        N₃'' = (6 − 12ξ)/L²,  N₄'' = (6ξ − 2)/L.
    """
    return np.array([
        (12.0 * xi - 6.0) / L**2,
        (6.0 * xi - 4.0) / L,
        (6.0 - 12.0 * xi) / L**2,
        (6.0 * xi - 2.0) / L,
    ])


class FrameSectionMixin:
    """Constitutiva de un marco: forma cerrada con material elástico lineal o
    integración de una sección (ADR 0021). Va delante de ``Element`` en la
    herencia.

    Contrato de la subclase
    -----------------------
    SECTION_DIM : ClassVar[int]
        Número de fuerzas internas de la sección que consume: 2 (``[N, M_z]``)
        en los marcos planos, 3 (``[N, M_y, M_z]``) en ``Frame3D``.
    SECTION_RESPONSE : ClassVar[str]
        Método de la sección que evalúa: ``"response_2d"`` o ``"response_3d"``.

    El constructor de la subclase llama a :meth:`_setup_constitution` antes
    de ``Element.__init__`` (que construye el estado con
    ``N_INTEGRATION_POINTS``).
    """

    SECTION_DIM: ClassVar[int] = 2
    SECTION_RESPONSE: ClassVar[str] = "response_2d"

    # Camino con material mientras la subclase no declare una sección.
    section = None

    # ------------------------------------------------------------------
    # Construcción
    # ------------------------------------------------------------------

    def _setup_constitution(self, element_id: int, material, section, n_ip,
                            geometry: dict, material_only: dict | None = None) -> None:
        """Valida que se declare exactamente uno de los dos caminos.

        ``geometry`` son las propiedades geométricas obligatorias del camino
        con material (``{"A": A, "I": I}``…); ``material_only``, parámetros
        opcionales que sólo tienen sentido con material (``{"nu": nu}``).
        Con ``section`` fija ``self.section``, ``self.n_ip``, la cuadratura
        y el atributo de instancia ``N_INTEGRATION_POINTS = n_ip``.
        """
        where = f"{type(self).__name__}(id={element_id})"
        if section is not None:
            if material is not None:
                raise ValueError(
                    f"{where}: declare 'material' con las propiedades geométricas "
                    f"{list(geometry)} o 'section', no ambos."
                )
            given = [k for k, v in {**geometry, **(material_only or {})}.items()
                     if v is not None]
            if given:
                raise ValueError(
                    f"{where}: con 'section' no se declaran {given}: la sección "
                    f"integra la rigidez y la masa de sus fibras."
                )
            if not isinstance(section, Section):
                raise ValueError(
                    f"{where}: 'section' debe ser una sección de marco (p. ej. "
                    f"FiberSection); se recibió {type(section).__name__}."
                )
            if (isinstance(n_ip, (bool, np.bool_))
                    or not isinstance(n_ip, (int, np.integer)) or n_ip < MIN_N_IP):
                raise ValueError(
                    f"{where}: n_ip={n_ip!r} debe ser un entero ≥ {MIN_N_IP}. Con 1 "
                    f"punto de Gauss el modo de flexión antisimétrica no tiene "
                    f"energía; con 2 la rigidez elástica se integra exactamente."
                )
            self.section = section
            self.n_ip = int(n_ip)
            self.N_INTEGRATION_POINTS = self.n_ip
            self._gauss_xi, self._gauss_w = gauss_legendre_unit(self.n_ip)
            return

        if material is None:
            raise ValueError(
                f"{where}: falta la constitutiva. Declare 'material' con las "
                f"propiedades geométricas {list(geometry)}, o 'section'."
            )
        missing = [k for k, v in geometry.items() if v is None]
        if missing:
            raise ValueError(
                f"{where}: con 'material' son obligatorias las propiedades "
                f"geométricas {missing}."
            )
        if n_ip != DEFAULT_N_IP:
            raise ValueError(
                f"{where}: n_ip={n_ip!r} sólo aplica con 'section'; con 'material' "
                f"el marco usa su forma cerrada, sin puntos de Gauss."
            )
        self.section = None

    def _validate_material_compatibility(self) -> None:
        """Validación de la base (dimensión, unilateralidad) y, después, la
        guarda de material elástico lineal del camino con material."""
        super()._validate_material_compatibility()
        material = getattr(self, "material", None)
        if material is None or getattr(material, "IS_LINEAR_ELASTIC", False):
            return
        name = type(self).__name__
        if "section" in getattr(type(self), "REFERENCE_KWARGS", {}):
            hint = (
                "Para un material inelástico declare una sección por fibras "
                "(FiberSection, bloque 'sections' del YAML) y pásela con "
                "'section': la plasticidad o el daño avanzan fibra a fibra."
            )
        else:
            from solidum.registry import ElementRegistry
            with_section = [n for n in ElementRegistry.names()
                            if "section" in ElementRegistry.get(n).REFERENCE_KWARGS]
            hint = (
                f"{name} todavía no admite sección por fibras; para plasticidad "
                f"distribuida use {with_section} con 'section' (FiberSection)."
            )
        raise ValueError(
            f"{name}(id={self.id}): el material {type(material).__name__} no es "
            f"elástico lineal (IS_LINEAR_ELASTIC = False). Con 'material' y "
            f"propiedades geométricas el marco usa su forma cerrada elástica, que "
            f"con un material inelástico escalaría toda la rigidez con la E_t de "
            f"la deformación axial y daría flectores incorrectos (ADR 0021). {hint}"
        )

    def _init_state(self) -> None:
        """Con sección: ``N_INTEGRATION_POINTS`` puntos, cada uno con la lista
        de estados iniciales de la sección y ``s_g = 0``."""
        section = getattr(self, "section", None)
        if section is None:
            super()._init_state()
            return
        n = self.N_INTEGRATION_POINTS
        self.state = ElementState(n, init_stress=np.zeros(self.SECTION_DIM))
        for g in range(n):
            self.state.vars[g] = section.initial_states()
            self.state.vars_trial[g] = section.initial_states()

    # ------------------------------------------------------------------
    # Contratos genéricos
    # ------------------------------------------------------------------

    def materials(self) -> list:
        """Materiales de las fibras con sección; ``[material]`` sin ella."""
        if self.section is not None:
            return self.section.materials()
        return super().materials()

    # ------------------------------------------------------------------
    # Integración de la sección
    # ------------------------------------------------------------------

    def _integrate_section(self, B: np.ndarray, weights: np.ndarray,
                           v: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """``K = Σ w_g B_gᵀ k_s,g B_g`` y ``F = Σ w_g B_gᵀ s_g`` con
        ``e_g = B_g·v``; guarda en el trial los estados y ``s_g`` de cada
        punto. ``B``: ``(n_ip, SECTION_DIM, n)``; ``weights`` ya incluye la
        longitud (``w_g·L``)."""
        response = getattr(self.section, self.SECTION_RESPONSE)
        state = self.state
        n = B.shape[2]
        F = np.zeros(n)
        K = np.zeros((n, n))
        for g in range(B.shape[0]):
            Bg = B[g]
            s, k_s, new_states = response(Bg @ v, state.vars[g])
            state.vars_trial[g] = new_states
            state.stresses_trial[g] = s
            F += weights[g] * (Bg.T @ s)
            K += weights[g] * (Bg.T @ k_s @ Bg)
        # Simetría exacta: k_s es simétrica, pero el producto triple puede
        # redondear distinto las entradas (i, j) y (j, i).
        return 0.5 * (K + K.T), F

    def _section_forces_at(self, B: np.ndarray, v: np.ndarray):
        """Deformaciones generalizadas ``e_g`` y fuerzas internas ``s_g`` de
        cada punto de Gauss para ``v``, a partir del estado committed y sin
        tocar el trial. Arreglos ``(n_ip, SECTION_DIM)``."""
        response = getattr(self.section, self.SECTION_RESPONSE)
        e = np.array([B[g] @ v for g in range(B.shape[0])])
        s = np.array([response(e[g], self.state.vars[g])[0] for g in range(B.shape[0])])
        return e, s

    def _evaluate_without_trace(self, u_e: np.ndarray):
        """``compute_element_state(u_e)`` sin dejar huella en el trial: la
        evaluación para post-proceso (fuerzas internas) no es un iterado del
        paso en curso, como ``compute_global_stiffness``."""
        snapshot = self.state.snapshot_trial()
        try:
            return self.compute_element_state(u_e)
        finally:
            self.state.restore_trial(snapshot)

    def _axial_area(self) -> float:
        """Área bruta: ``A`` con material, ``Σ A_i`` con sección. La usa la
        carga de una fuerza de cuerpo uniforme (``q = A·b``)."""
        if self.section is None:
            return self.A
        return self.section.total_area()
