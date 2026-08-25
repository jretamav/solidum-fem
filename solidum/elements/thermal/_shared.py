# solidum_fem/solidum/elements/thermal/_shared.py
"""Base común de los elementos térmicos (Etapa 8).

Centraliza lo que comparten ``Quad4Thermal`` y ``Hex8Thermal``: validación
del material, bucle de Gauss para la matriz de conductividad, matriz de
capacidad (consistente y lumped), fuente volumétrica y post-proceso.

Las subclases sólo declaran su geometría —funciones de forma, gradiente,
cuadratura por defecto, topología de la frontera— y heredan el resto.

Se introduce con el primer elemento térmico en vez de esperar al segundo
porque la ecuación discretizada ``C·Ṫ + K·T = F`` es idéntica en 2D y 3D:
la dimensión sólo cambia el tamaño de ``B``, no la formulación. No es una
abstracción especulativa sino la constatación de que el segundo caso real
—el ``Hex8Thermal`` de esta misma etapa— difiere únicamente en geometría.
"""
from __future__ import annotations

from abc import abstractmethod
from typing import ClassVar, List

import numpy as np

from solidum.core.element import Element, validate_lumping_kwarg
from solidum.core.node import Node
from solidum.core.thermal_material import ThermalMaterial
from solidum.math.mass_lumping import lump_hrz


class _ThermalSolid(Element):
    """Elemento de conducción de calor: un DOF escalar ``T`` por nodo.

    Contrato de subclases
    ---------------------
    ``FLUX_DIM`` : dimensión del gradiente (2 ó 3).
    ``N_NODES`` : número de nodos.
    ``DEFAULT_QUADRATURE`` : clave en ``QuadratureRegistry``.
    ``_shape_functions(p)`` : ``N`` en el punto natural ``p``.
    ``_gradient(p, coords)`` : ``(dN_dx, detJ)`` — la matriz ``B`` térmica.
    ``_boundary_load(index, q_bar)`` : vector nodal del flujo de frontera.

    Notas
    -----
    ``STRAIN_DIM`` se declara como ``None`` para desactivar la validación
    dimensional de la base mecánica, que compara contra el ``STRAIN_DIM``
    del material. Aquí la compatibilidad se comprueba contra ``FLUX_DIM``
    en :meth:`_validate_material_compatibility`, que se sobreescribe.
    """

    DOF_NAMES: ClassVar[List[str]] = ["T"]
    STRAIN_DIM = None            # no aplica: el campo es escalar
    FLUX_DIM: ClassVar[int]
    N_NODES: ClassVar[int]
    DEFAULT_QUADRATURE: ClassVar[str]

    # Espesor: sólo lo usan los elementos planos. Los 3D lo dejan en 1.0 y
    # el volumen sale de la geometría.
    thickness: float = 1.0

    def __init__(
        self,
        element_id: int,
        nodes: List[Node],
        material: ThermalMaterial,
        quadrature: str | None = None,
    ):
        from solidum.registry import QuadratureRegistry

        if len(nodes) != self.N_NODES:
            raise ValueError(
                f"{type(self).__name__}(id={element_id}): requiere "
                f"{self.N_NODES} nodos, recibió {len(nodes)}."
            )

        key = quadrature if quadrature is not None else self.DEFAULT_QUADRATURE
        self.points, self.weights = QuadratureRegistry.get(key)
        self.quadrature_key = key
        self.N_INTEGRATION_POINTS = len(self.points)

        super().__init__(element_id, nodes, material)

    # ------------------------------------------------------------------
    # Contrato geométrico que cada subclase concreta
    # ------------------------------------------------------------------

    @abstractmethod
    def _shape_functions(self, point) -> np.ndarray:
        """Funciones de forma evaluadas en el punto natural."""

    @abstractmethod
    def _gradient(self, point, coords):
        """Devuelve ``(dN_dx, detJ)`` — la matriz ``B`` térmica y el jacobiano."""

    # ------------------------------------------------------------------
    # Validación e inicialización
    # ------------------------------------------------------------------

    def _validate_material_compatibility(self) -> None:
        """Comprueba familia y dimensión del material térmico.

        Sustituye la validación de la base mecánica, que compara
        ``STRAIN_DIM``. Un material mecánico pasado por error se detecta
        aquí con mensaje explícito, en vez de fallar más tarde con un
        ``AttributeError`` sobre ``compute_flux``.
        """
        mat = self.material
        if mat is None:
            raise ValueError(
                f"{type(self).__name__}(id={self.id}): requiere un material "
                "térmico; no admite construirse sin material."
            )
        if not isinstance(mat, ThermalMaterial):
            raise TypeError(
                f"{type(self).__name__}(id={self.id}): el material "
                f"{type(mat).__name__} no pertenece a la familia térmica "
                "(ThermalMaterial). Los elementos térmicos consumen materiales "
                "que relacionan flujo con gradiente de temperatura "
                "(compute_flux), no esfuerzo con deformación. Use "
                "ThermalConduction o similar."
            )
        mat_dim = getattr(mat, "FLUX_DIM", None)
        if mat_dim is not None and mat_dim != self.FLUX_DIM:
            raise ValueError(
                f"{type(self).__name__}(id={self.id}): incompatibilidad "
                f"dimensional. El elemento requiere FLUX_DIM={self.FLUX_DIM} "
                f"pero el material {type(mat).__name__} declara "
                f"FLUX_DIM={mat_dim}. Construya el material con la dimensión "
                "del problema (dim=2/3, o un tensor k del tamaño adecuado)."
            )

    def _init_state(self) -> None:
        """La conducción de Fourier no tiene variables internas.

        Se omite el ``ElementState`` de la base mecánica: no hay historia
        que promover trial → committed. Cuando entre un material térmico
        con memoria (``k(T)``, cambio de fase), este método se reimplementa
        con el caso real delante.
        """
        self.state = None

    def commit_state(self) -> None:
        """Sin estado que confirmar (modelo lineal sin memoria)."""
        return

    # ------------------------------------------------------------------
    # Matrices elementales
    # ------------------------------------------------------------------

    def compute_conductivity_matrix(self) -> np.ndarray:
        """Matriz de conductividad ``K_e = ∫ Bᵀ k B t dΩ``.

        Simétrica y **semidefinida positiva**: tiene el modo nulo del campo
        de temperatura uniforme, que no produce gradiente y por tanto
        tampoco flujo. Es el análogo térmico de los modos de sólido rígido.
        """
        n = self.N_NODES
        K_e = np.zeros((n, n))
        coords = self.get_coordinate_matrix(ndim=self.FLUX_DIM)
        k_tensor = self.material.conductivity

        for p, w in zip(self.points, self.weights):
            B, detJ = self._gradient(p, coords)
            K_e += (B.T @ k_tensor @ B) * (detJ * w * self.thickness)

        return K_e

    def compute_capacity_matrix(self, lumping: str = "lumped") -> np.ndarray:
        """Matriz de capacidad ``C_e = ∫ ρc Nᵀ N t dΩ``.

        Parameters
        ----------
        lumping : {"lumped", "consistent"}
            **El default es ``"lumped"``, al revés que en dinámica
            estructural.** La razón es física: la capacidad consistente
            produce oscilaciones espurias ante un frente térmico abrupto,
            con temperaturas que pueden salirse del rango de los datos —
            violando el principio del máximo de la ecuación de difusión.
            La lumped amortigua esas oscilaciones y preserva la monotonía.
            Ver ``docs/specs/Quad4Thermal.md`` §5.
        """
        validate_lumping_kwarg(lumping, type(self).__name__)

        rho_c = self.material.volumetric_capacity(
            consumer="la matriz de capacidad calorífica"
        )

        n = self.N_NODES
        C_e = np.zeros((n, n))
        coords = self.get_coordinate_matrix(ndim=self.FLUX_DIM)
        volume = 0.0

        for p, w in zip(self.points, self.weights):
            N = self._shape_functions(p)
            _, detJ = self._gradient(p, coords)
            dV = detJ * w * self.thickness
            C_e += np.outer(N, N) * (rho_c * dV)
            volume += dV

        if lumping == "lumped":
            # HRZ canónico. El campo es escalar —un DOF por nodo—, así que
            # hay una sola "dirección" y la capacidad total del elemento es
            # ρc·V_e. El esquema escala la diagonal para conservarla, que es
            # exactamente lo que se necesita: la energía almacenable no debe
            # depender de cómo se distribuya entre nodos.
            C_e = lump_hrz(
                C_e,
                total_mass=rho_c * volume,
                n_translational_dirs=1,
            )

        return C_e

    # La base mecánica pide `compute_element_state`; en el térmico el
    # sistema es lineal y sin historia, así que la matriz es la de
    # conductividad y no hay fuerzas internas dependientes del estado.
    def compute_element_state(self, u_e: np.ndarray):
        K_e = self.compute_conductivity_matrix()
        return K_e, K_e @ u_e

    def compute_global_stiffness(self) -> np.ndarray:
        """Alias térmico: la 'rigidez' del problema es la conductividad."""
        return self.compute_conductivity_matrix()

    def compute_mass_matrix(self, lumping: str = "lumped") -> np.ndarray:
        """Alias térmico de la capacidad, para el pipeline transitorio."""
        return self.compute_capacity_matrix(lumping=lumping)

    # ------------------------------------------------------------------
    # Cargas
    # ------------------------------------------------------------------

    def compute_body_source(self, Q: float) -> np.ndarray:
        """Vector nodal de la fuente volumétrica ``f = ∫ Q Nᵀ t dΩ``.

        ``Q`` en [W/m³]. Para ``Q`` uniforme la suma de las componentes es
        exactamente ``Q·V_e`` (``Q·A_e·t`` en 2D), invariante ante
        distorsión de la geometría.

        La fuente es carga del elemento y no propiedad del material: la
        generación por hidratación o efecto Joule varía en espacio y
        tiempo, y como propiedad exigiría duplicar materiales.
        """
        n = self.N_NODES
        f = np.zeros(n)
        coords = self.get_coordinate_matrix(ndim=self.FLUX_DIM)

        for p, w in zip(self.points, self.weights):
            N = self._shape_functions(p)
            _, detJ = self._gradient(p, coords)
            f += N * (float(Q) * detJ * w * self.thickness)

        return f

    # ------------------------------------------------------------------
    # Post-proceso
    # ------------------------------------------------------------------

    def compute_gauss_state(self, T_global: np.ndarray) -> dict:
        """Gradiente y flujo en cada punto de Gauss.

        Paralelo al contrato de los sólidos (ADR 0012). No expone
        ``internal_forces``: ese método corresponde a elementos
        estructurales 1D por el cierre por dominio del propio ADR.
        """
        T_e = self.get_local_displacements(T_global)
        coords = self.get_coordinate_matrix(ndim=self.FLUX_DIM)

        n_g = len(self.points)
        grads = np.zeros((n_g, self.FLUX_DIM))
        fluxes = np.zeros((n_g, self.FLUX_DIM))
        pts_global = np.zeros((n_g, self.FLUX_DIM))

        for idx, p in enumerate(self.points):
            B, _ = self._gradient(p, coords)
            grad_T = B @ T_e
            q, _ = self.material.compute_flux(grad_T)
            grads[idx] = grad_T
            fluxes[idx] = q
            pts_global[idx] = self._shape_functions(p) @ coords

        return {
            "points_natural": np.asarray(self.points, dtype=float),
            "points_global": pts_global,
            "grad_T": grads,
            "flux": fluxes,
        }
