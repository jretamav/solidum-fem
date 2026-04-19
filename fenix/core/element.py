# fenix_fem/fenix/core/element.py
from abc import ABC, abstractmethod
from typing import ClassVar, List, Tuple
import numpy as np
from fenix.core.node import Node

class Element(ABC):
    """
    Clase base abstracta para todos los elementos finitos de Fenix FEM.

    Contrato obligatorio para subclases
    ------------------------------------
    DOF_NAMES : ClassVar[List[str]]
        Lista canónica de DOF por nodo (sin repetir). Ejemplo Quad4: ['ux', 'uy'].
        El número de DOFs del elemento es len(DOF_NAMES) × len(nodes).
    compute_element_state(u_e)  — rigidez tangente + fuerzas internas dado u_e.
    commit_state()              — confirma variables internas tras convergencia del paso.

    Implementación heredada (no sobreescribir salvo necesidad)
    -----------------------------------------------------------
    get_global_dof_indices()    — mapea DOF_NAMES × nodos a índices globales.
    compute_global_stiffness()  — llama compute_element_state(zeros); usada en análisis lineal.
    """

    # Cada subclase concreta DEBE declarar este atributo de clase.
    # Ejemplo: DOF_NAMES = ['ux', 'uy']
    DOF_NAMES: ClassVar[List[str]]

    def __init__(self, element_id: int, nodes: List[Node]):
        self.id = element_id
        self.nodes = nodes

    # ------------------------------------------------------------------
    # Contrato abstracto — toda subclase DEBE implementar estos métodos
    # ------------------------------------------------------------------

    @abstractmethod
    def compute_element_state(self, u_e: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calcula la matriz de rigidez tangente local K_e y el vector de fuerzas
        internas F_int_e dado el vector de desplazamientos locales u_e.

        Actualiza las variables internas en estado *trial* (sin commitear).

        Returns
        -------
        K_e      : ndarray (ndof_e × ndof_e) — matriz de rigidez tangente local
        F_int_e  : ndarray (ndof_e,)          — vector de fuerzas internas locales
        """

    @abstractmethod
    def commit_state(self) -> None:
        """
        Confirma las variables internas del estado trial como estado convergido.
        Debe llamarse una vez por paso de carga, después de que el solver converja.
        Para materiales puramente elásticos puede ser un no-op, pero debe existir.
        """

    # ------------------------------------------------------------------
    # Implementación base — disponible para todas las subclases
    # ------------------------------------------------------------------

    def get_global_dof_indices(self) -> List[int]:
        """
        Mapea DOF_NAMES × nodos a sus índices globales de ecuación.

        Itera sobre cada nodo y para cada nombre en DOF_NAMES extrae el índice
        del diccionario node.dofs. No asume DOFs uniformes: cada nodo puede
        contener DOFs distintos siempre que incluya los de DOF_NAMES.

        Returns
        -------
        indices : List[int] — longitud = len(DOF_NAMES) × len(nodes)
        """
        indices = []
        for node in self.nodes:
            for dof_name in self.DOF_NAMES:
                indices.append(node.dofs[dof_name])
        return indices

    def compute_global_stiffness(self) -> np.ndarray:
        """
        Devuelve la matriz de rigidez global K_e evaluada en el estado no deformado
        (u_e = 0). Usada por el LinearSolver y como punto de partida en análisis lineal.

        No sobreescribir: la implementación base delega en compute_element_state(),
        garantizando consistencia entre el análisis lineal y el no-lineal.
        """
        ndof_e = len(self.DOF_NAMES) * len(self.nodes)
        K_e, _ = self.compute_element_state(np.zeros(ndof_e))
        return K_e

    def get_coordinate_matrix(self, ndim: int = 2) -> np.ndarray:
        """Extrae las coordenadas de los nodos en una matriz NumPy (n_nodos × ndim)."""
        coords = np.zeros((len(self.nodes), ndim))
        for i, node in enumerate(self.nodes):
            for j in range(min(ndim, len(node.coordinates))):
                coords[i, j] = node.coordinates[j]
        return coords
