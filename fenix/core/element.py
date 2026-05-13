# fenix_fem/fenix/core/element.py
from abc import ABC, abstractmethod
from typing import ClassVar, List, Optional, Tuple

import numpy as np

from fenix.core.element_state import ElementState
from fenix.core.material import Material
from fenix.core.node import Node
from fenix.results import ElementForces


class Element(ABC):
    """Clase base abstracta para todos los elementos finitos de Fenix FEM.

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
                 material: Optional[Material] = None):
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
    # API pública de resultados (ADR 0002)
    # ------------------------------------------------------------------

    def internal_forces(self, U_global: np.ndarray) -> Optional[ElementForces]:
        """Esfuerzos internos en ejes locales del elemento, en los nodos i, j.

        Convenciones de signo fijadas en Reglas.md §5:
        - Frames 2D, trusses y cables: convención de viga estructural.
        - Frames 3D: stress-resultant / RHR pura.

        Los elementos corotacionales deben usar su ``state`` ya comprometido
        tras la solución (no recomputar desde cero).

        Returns
        -------
        ElementForces | None
            ``None`` por defecto (elementos sólidos y cualquier elemento para
            el que N/V/M no aplican). Subclases barra/viga/cable sobreescriben.
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
