# fenix_fem/solidum/core/element.py
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import ClassVar, List, Tuple

import numpy as np

from solidum.core.element_forces import ElementForces
from solidum.core.element_state import ElementState
from solidum.core.material import Material
from solidum.core.node import Node


# Valores válidos del kwarg ``lumping`` de :meth:`Element.compute_mass_matrix`
# (ADR 0009 fase 2). Centralizado aquí para que la validación no se duplique
# en cada subclase concreta (auditoría H-1.6).
SUPPORTED_LUMPING: frozenset[str] = frozenset({"consistent", "lumped"})


def validate_lumping_kwarg(lumping: str, elem_class_name: str) -> None:
    """Valida el kwarg ``lumping`` contra el set canónico del proyecto.

    Lanza ``NotImplementedError`` con mensaje uniforme si el valor no
    está soportado. El primer argumento es el string recibido; el
    segundo es el nombre de la clase concreta, para que el mensaje le
    sea útil al usuario.
    """
    if lumping not in SUPPORTED_LUMPING:
        raise NotImplementedError(
            f"{elem_class_name}.compute_mass_matrix: lumping={lumping!r} "
            f"no soportado. Valores admitidos: {sorted(SUPPORTED_LUMPING)}."
        )


class Element(ABC):
    """Clase base abstracta para todos los elementos finitos de Solidum FEM.

    Contrato declarativo de subclases
    ---------------------------------
    DOF_NAMES : ClassVar[List[str]]
        Lista canónica de DOFs por nodo (ej. Quad4: ['ux', 'uy']).
    STRAIN_DIM : ClassVar[int]
        Dimensión de strain que el elemento requiere del material:
            1  → 1D (Truss, Frame; epsilon escalar)
            3  → Voigt 2D plano
            6  → Voigt 3D
        Validado contra `material.STRAIN_DIM` al construir.
    N_INTEGRATION_POINTS : ClassVar[int], default=1
        Puntos de Gauss del elemento. Determina el tamaño del ElementState.
    PRESERVES_SYMMETRY : ClassVar[bool], default=True
        Indica si la matriz tangente del elemento es simétrica cuando lo es
        la del material. ``True`` para todos los elementos derivados de un
        funcional de energía con cargas conservativas (formulación
        variacional estándar). ``False`` para follower loads (presión que
        sigue a la superficie deformada) y formulaciones corrotacionales con
        rotaciones finitas no simetrizadas. La capa algebraica (ADR 0003)
        lo agrega con ``material.IS_SYMMETRIC`` para elegir backend.
    ACCEPTS_UNILATERAL : ClassVar[bool], default=False
        Indica si el elemento es robusto frente a materiales unilaterales
        (``material.IS_UNILATERAL = True``), cuyo módulo tangente puede
        colapsar a cero en algún régimen. Solo elementos preparados para
        manejar rigidez nula localmente (típicamente corotacionales con
        K_G no degenerado) deben sobreescribir a ``True``. Validado en
        construcción para evitar que un truss ordinario acepte un material
        de cable.

    Métodos abstractos a implementar
    --------------------------------
    compute_element_state(u_e) → (K_e, F_int_e)

    Métodos opcionales (implementan análisis específicos)
    -----------------------------------------------------
    compute_mass_matrix(lumping)
        Matriz de masa elemental en ejes globales. Necesaria para análisis
        modal y dinámico (ADR 0009). Si no se sobreescribe, la base lanza
        ``NotImplementedError`` con mensaje claro al ser invocada. El
        contrato declarativo acepta el parámetro ``lumping`` desde el día
        uno (fase 1: solo ``"consistent"``; fases futuras: ``"lumped"``,
        ``"hrz"``, ...) para que la extensión sea aditiva.
    compute_body_load(b), compute_edge_traction(edge, t_vec)
        Cargas distribuidas consistentes. Opcionales según el tipo.

    Métodos heredados (no sobreescribir salvo necesidad)
    ----------------------------------------------------
    commit_state(), state_vars, stresses, get_global_dof_indices,
    get_local_displacements, compute_global_stiffness.
    """

    DOF_NAMES: ClassVar[List[str]]
    STRAIN_DIM: ClassVar[int]
    N_INTEGRATION_POINTS: ClassVar[int] = 1
    PRESERVES_SYMMETRY: ClassVar[bool] = True
    ACCEPTS_UNILATERAL: ClassVar[bool] = False

    def __init__(self, element_id: int, nodes: List[Node],
                 material: Material | None = None):
        self.id = element_id
        self.nodes = nodes
        self.material = material

        self._validate_material_compatibility()
        self._register_dofs()
        self._init_state()

    # ------------------------------------------------------------------
    # Inicialización compartida
    # ------------------------------------------------------------------

    def _validate_material_compatibility(self) -> None:
        """Atrapa al construir el mismatch material/elemento (1D vs 2D, etc).

        Sin esto el error aparece como `matmul` críptico tras varias iteraciones
        del solver. Con esto, falla inmediatamente con un mensaje claro.
        """
        if self.material is None:
            return
        mat_dim = getattr(self.material, "STRAIN_DIM", None)
        if mat_dim is None:
            return
        if mat_dim != self.STRAIN_DIM:
            raise ValueError(
                f"{type(self).__name__}(id={self.id}): incompatibilidad dimensional. "
                f"El elemento requiere STRAIN_DIM={self.STRAIN_DIM} pero el material "
                f"{type(self.material).__name__} declara STRAIN_DIM={mat_dim}. "
                f"Use un material 1D con elementos Truss/Frame, 2D con Quad/Tri, etc."
            )
        if getattr(self.material, "IS_UNILATERAL", False) and not self.ACCEPTS_UNILATERAL:
            raise ValueError(
                f"{type(self).__name__}(id={self.id}): material unilateral "
                f"{type(self.material).__name__} no admitido. Su módulo tangente "
                f"puede colapsar a cero y degenerar la matriz global en este "
                f"elemento. Use un elemento que lo acepte (p. ej. Cable2DCorot o "
                f"Cable3DCorot) o sustituya el material por uno bilateral."
            )

    def _register_dofs(self) -> None:
        for node in self.nodes:
            for dof_name in self.DOF_NAMES:
                node.add_dof(dof_name)

    def _init_state(self) -> None:
        init_stress = 0.0 if self.STRAIN_DIM == 1 else np.zeros(self.STRAIN_DIM)
        self.state = ElementState(self.N_INTEGRATION_POINTS, init_stress=init_stress)

    # ------------------------------------------------------------------
    # Acceso al estado interno (heredado por todas las subclases)
    # ------------------------------------------------------------------

    @property
    def state_vars(self):
        return self.state.vars

    @property
    def stresses(self):
        return self.state.stresses

    def commit_state(self) -> None:
        """Confirma las variables internas trial como estado convergido.

        Se llama una vez por paso de carga, después de que el solver converja.
        Para materiales puramente elásticos es esencialmente un no-op pero el
        método debe existir.
        """
        self.state.commit()

    # ------------------------------------------------------------------
    # Contrato abstracto
    # ------------------------------------------------------------------

    @abstractmethod
    def compute_element_state(self, u_e: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Calcula la matriz tangente local K_e y el vector de fuerzas internas F_int_e.

        Actualiza variables internas en estado *trial* (sin commitear).

        Returns
        -------
        K_e      : ndarray (ndof_e × ndof_e)
        F_int_e  : ndarray (ndof_e,)
        """

    # ------------------------------------------------------------------
    # Hook de preparación de paso (ADR 0010 §5)
    # ------------------------------------------------------------------

    def prepare_step(self, U_committed: np.ndarray) -> None:
        """Prepara el elemento al inicio de un paso de carga, antes del Newton.

        Hook invocado por los solvers no lineales una vez por paso, con el
        campo de desplazamientos **convergido** del paso anterior. La
        implementación base es no-op; lo sobreescriben los elementos que
        necesitan tomar decisiones discretas estables entre pasos (p. ej.
        activación de discontinuidad embebida en ``CST_Embedded2D``, ADR
        0010 §5) y que se romperían si las hicieran dentro del Newton (donde
        los predictores lineales oscilan sobre el umbral y producen
        chattering).

        Parameters
        ----------
        U_committed : np.ndarray
            Campo global de desplazamientos convergido del paso anterior. El
            elemento extrae los DOFs locales con ``get_local_displacements``
            si los necesita.
        """
        # No-op por defecto. Subclases pueden sobreescribir.

    # ------------------------------------------------------------------
    # Contrato opcional para análisis modal y dinámico (ADR 0009)
    # ------------------------------------------------------------------

    def compute_mass_matrix(self, lumping: str = "consistent") -> np.ndarray:
        """Matriz de masa elemental en ejes globales, ordenada como ``compute_element_state``.

        El contrato base lanza ``NotImplementedError`` para señalar que la
        subclase no la implementa. El ``Assembler.assemble_mass_matrix``
        intercepta el fallo y reporta de forma agregada qué tipos de
        elemento carecen de masa, evitando que un análisis modal se cuelgue
        en runtime con un mensaje críptico.

        Parameters
        ----------
        lumping : {"consistent"}, default "consistent"
            Estrategia de discretización de la inercia (ADR 0009 §1). En la
            fase 1 solo se admite ``"consistent"`` (masa de Galerkin
            ``∫ρ·NᵀN dΩ``); las subclases deben lanzar
            ``NotImplementedError`` para otros valores hasta que la rama
            correspondiente se implemente (e.g. ``"lumped"``, ``"hrz"``).
            El parámetro vive en la firma desde el día uno para que añadir
            esas ramas sea aditivo y no rompa llamadores existentes.

        Returns
        -------
        np.ndarray, shape (ndof_e, ndof_e)
            Matriz de masa simétrica, positiva semi-definida.
        """
        raise NotImplementedError(
            f"{type(self).__name__} no implementa `compute_mass_matrix`. "
            f"El análisis solicitado (modal o dinámico) requiere la matriz "
            f"de masa elemental. Implementa el método en la subclase o "
            f"retira este elemento del modelo."
        )

    # ------------------------------------------------------------------
    # API pública de resultados (ADR 0002)
    # ------------------------------------------------------------------

    def internal_forces(self, U_global: np.ndarray) -> ElementForces | None:
        """Esfuerzos internos en ejes locales — **solo elementos estructurales 1D**.

        Contrato canónico del ADR 0002 con **dominio acotado por ADR 0012**:
        aplica a elementos estructurales 1D (trusses, cables, frames 2D/3D)
        que sí tienen "esfuerzos seccionales" discretos por extremo (N, V,
        M, T). Devuelve un :class:`ElementForces` inmutable. Convenciones
        de signo en Reglas.md §5:

        - Frames 2D, trusses y cables: convención de viga estructural
          (``V`` con signo invertido respecto al interno; ``M`` positivo
          = sagging con ``+y`` arriba).
        - Frames 3D: stress-resultant / RHR pura sin inversión de signos.

        Los elementos corotacionales deben usar su ``state`` ya comprometido
        tras la solución (no recomputar desde cero).

        Returns
        -------
        ElementForces | None
            ``None`` para **sólidos 2D y 3D** (Tri3, Quad4, Tri6, Quad8,
            Quad9, CST_Embedded2D, Hex8, Tet4). En un sólido continuo el
            equivalente de "esfuerzo seccional" es el campo tensorial
            ``σ(x)``, accesible vía :meth:`compute_gauss_state(U)`. La
            conversión a un único valor representativo por elemento
            (promedio, centroide, ponderado por volumen) la decide cada
            consumidor según su semántica — no la inventa este contrato.

            Las subclases barra/viga/cable sobrescriben este método y
            devuelven el ``ElementForces`` tipado.

        Notes
        -----
        Cierre del ADR 0002 con dominio explícito (ADR 0012, 2026-05-19):
        no es deuda pendiente, es decisión arquitectural. Sólidos exponen
        :meth:`compute_gauss_state` como API canónica de salida; estructurales
        exponen :meth:`internal_forces`. Las dos APIs coexisten porque
        responden a primitivas físicamente distintas — ``σ(x)`` campo
        tensorial vs ``(N, V, M, T)`` resultantes seccionales.

        Adicionalmente, varios elementos exponen un método legado
        :meth:`compute_internal_forces` que devuelve un ``dict`` ad-hoc
        por familia (``{N, V, M}`` para frames, ``{stress, strain}``
        promediado para sólidos). Es la API histórica para post-proceso
        libre (VTK exporter, scripts) y **no es el contrato del ADR 0002**.
        Coexisten ambas.
        """
        return None

    # ------------------------------------------------------------------
    # Utilidades heredadas
    # ------------------------------------------------------------------

    def get_global_dof_indices(self) -> List[int]:
        """Mapea DOF_NAMES × nodos a índices globales de ecuación."""
        indices = []
        for node in self.nodes:
            for dof_name in self.DOF_NAMES:
                indices.append(node.dofs[dof_name])
        return indices

    def get_local_displacements(self, U_global: np.ndarray) -> np.ndarray:
        """Extrae el vector u_e del campo de desplazamientos global."""
        return U_global[self.get_global_dof_indices()]

    def compute_global_stiffness(self) -> np.ndarray:
        """Matriz de rigidez K_e evaluada en el estado no deformado (u_e = 0).

        Punto de partida del análisis lineal y delegación garantizada en la
        misma `compute_element_state`, asegurando consistencia lineal/no-lineal.
        """
        ndof_e = len(self.DOF_NAMES) * len(self.nodes)
        K_e, _ = self.compute_element_state(np.zeros(ndof_e))
        return K_e

    def get_coordinate_matrix(self, ndim: int = 2) -> np.ndarray:
        """Coordenadas de los nodos en una matriz NumPy (n_nodos × ndim)."""
        coords = np.zeros((len(self.nodes), ndim))
        for i, node in enumerate(self.nodes):
            for j in range(min(ndim, len(node.coordinates))):
                coords[i, j] = node.coordinates[j]
        return coords
