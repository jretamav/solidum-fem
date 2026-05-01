# fenix_fem/fenix/core/domain.py
from typing import Dict, Optional, TYPE_CHECKING
from fenix.core.node import Node
from fenix.logging import get_logger

if TYPE_CHECKING:
    from fenix.results import SolveResult

_log = get_logger("core.domain")


class Domain:
    def __init__(self):
        self.nodes: Dict[int, Node] = {}
        self.elements = {}
        self.total_dofs: int = 0
        # Resultado agregado de la última solución (ADR 0002). Lo asigna el
        # entrypoint fenix.run tras solver.solve vía build_solve_result.
        self.last_result: Optional["SolveResult"] = None
        # Restricciones afines lineales (MPC) declaradas por el usuario
        # (ADR 0004 fase 2). Cada entrada: dict con 'slave' = (node_id,
        # dof_name), 'masters' = [(node_id, dof_name), ...], 'coefficients',
        # y opcional 'g'. El Assembler las traduce a DOFs globales y las
        # acumula en su ConstraintSet.
        self.linear_constraints: list[dict] = []

    def add_linear_constraint(
        self,
        slave: tuple[int, str],
        masters: list[tuple[int, str]],
        coefficients: list[float],
        g: float = 0.0,
    ) -> None:
        """Declara una restricción afín ``u_s = g + Σ α_i · u_mi``.

        Parameters
        ----------
        slave
            Par ``(node_id, dof_name)`` del DOF esclavo.
        masters
            Lista de pares ``(node_id, dof_name)`` de los DOFs maestro.
        coefficients
            Coeficientes ``α_i`` en el mismo orden que ``masters``.
        g
            Término independiente (cero por defecto).

        Casos de uso típicos: apoyo en plano oblicuo, periodicidad, unión
        rígida, simetría no alineada con ejes globales.
        """
        if len(masters) != len(coefficients):
            raise ValueError(
                "masters y coefficients deben tener la misma longitud "
                f"({len(masters)} vs {len(coefficients)})."
            )
        self.linear_constraints.append({
            "slave": slave,
            "masters": list(masters),
            "coefficients": list(coefficients),
            "g": float(g),
        })

    def add_node(self, node_id: int, coordinates: list[float]) -> Node:
        if node_id in self.nodes:
            raise ValueError(f"El nodo {node_id} ya existe.")
        new_node = Node(node_id, coordinates)
        self.nodes[node_id] = new_node
        return new_node

    def get_node(self, node_id: int) -> Node:
        return self.nodes.get(node_id)

    def add_element(self, element) -> None:
        if element.id in self.elements:
            raise ValueError(f"El elemento {element.id} ya existe.")
        self.elements[element.id] = element

    def generate_equation_numbers(self, verbose: bool = False):
        eq_number = 0
        for node_id, node in self.nodes.items():
            for dof_name in node.dofs.keys():
                node.dofs[dof_name] = eq_number
                eq_number += 1
        self.total_dofs = eq_number
        if verbose:
            _log.info(f"Numeración completada. Grados de libertad totales: {self.total_dofs}")
