# solidum_fem/solidum/core/domain.py
from __future__ import annotations

import math
from numbers import Real
from typing import Dict, TYPE_CHECKING

from solidum.core.dof_names import MECHANICAL, ROTATIONS, TRANSLATIONS
from solidum.core.node import Node
from solidum.logging import get_logger

if TYPE_CHECKING:
    from solidum.results import ModalResult, SolveResult

_log = get_logger("core.domain")


class Domain:
    def __init__(self):
        self.nodes: Dict[int, Node] = {}
        self.elements = {}
        self.total_dofs: int = 0
        # Resultado agregado de la última solución (ADR 0002). Para análisis
        # estáticos lo asigna ``solidum.run`` vía ``build_solve_result``
        # (``SolveResult``); para análisis modal lo asigna ``solidum.run_modal``
        # directamente (``ModalResult``). Ambos son frozen dataclasses.
        self.last_result: SolveResult | ModalResult | None = None
        # Restricciones afines lineales (MPC) declaradas por el usuario
        # (ADR 0004 fase 2). Cada entrada: dict con 'slave' = (node_id,
        # dof_name), 'masters' = [(node_id, dof_name), ...], 'coefficients',
        # y opcional 'g'. El Assembler las traduce a DOFs globales y las
        # acumula en su ConstraintSet.
        self.linear_constraints: list[dict] = []
        # Nodos que crea Solidum, no el usuario: los de una barra embebida
        # (ADR 0032). Los selectores por coordenada del YAML no los toman:
        # una carga o un apoyo por coordenada es del sólido.
        self.generated_nodes: set[int] = set()
        # Masas concentradas en nudos (ADR 0024): (node_id, dof) → masa [kg]
        # en un grado de libertad de traslación, inercia rotacional [kg·m²]
        # en uno de giro. El Assembler las suma a la diagonal de M.
        self.nodal_masses: dict[tuple[int, str], float] = {}
        # ``True`` en cuanto un análisis compromete estado interno en los
        # elementos (``Assembler.commit_all_states``): el modelo ya no está en
        # su estado inicial (ADR 0026).
        self.committed_history: bool = False

    def add_nodal_mass(self, node_id: int, **masses: float) -> None:
        """Masa concentrada en un nudo (ADR 0024).

        Ejemplo::

            domain.add_nodal_mass(5, ux=12.0e3, uy=12.0e3)   # masa [kg]
            domain.add_nodal_mass(5, rz=4.0e3)               # inercia [kg·m²]

        Cada argumento nombra un grado de libertad mecánico del nodo:
        ``ux``, ``uy``, ``uz`` reciben una masa; ``rx``, ``ry``, ``rz``, una
        inercia rotacional (masa·longitud²). La masa entra en la diagonal de
        la matriz de masa con los dos esquemas (``consistent`` y ``lumped``),
        así que la ven el análisis modal, el espectral y los dinámicos.

        Es **sólo inercia**: ``gravity`` no genera su peso. Si su peso actúa,
        declárelo como carga.

        Varias llamadas sobre el mismo grado de libertad se suman. Que el
        nodo tenga ese grado de libertad se comprueba al ensamblar la matriz
        de masa, porque los grados de libertad los crean los elementos.

        Raises
        ------
        ValueError
            Si el nodo no existe, si no se da ningún grado de libertad, si
            alguno no es mecánico (``T``) o si algún valor no es un número
            finito no negativo. Todos los errores se informan juntos.
        """
        if node_id not in self.nodes:
            raise ValueError(f"add_nodal_mass: no hay ningún nodo {node_id} en el modelo.")
        if not masses:
            raise ValueError(
                "add_nodal_mass: indica al menos un grado de libertad, p. ej. "
                f"add_nodal_mass({node_id}, ux=1000.0)."
            )
        errors = []
        for dof, value in masses.items():
            if dof not in MECHANICAL:
                errors.append(
                    f"'{dof}' no es un grado de libertad mecánico; admitidos: "
                    f"{list(TRANSLATIONS)} (masa) y {list(ROTATIONS)} (inercia "
                    f"rotacional). Un grado de libertad escalar como 'T' no tiene "
                    f"masa: su análogo es una capacidad calorífica."
                )
            if (not isinstance(value, Real) or isinstance(value, bool)
                    or not math.isfinite(float(value)) or float(value) < 0.0):
                errors.append(f"{dof}={value!r}: debe ser un número finito no negativo.")
        if errors:
            raise ValueError(
                f"add_nodal_mass (nodo {node_id}):\n" + "\n".join(f"  - {e}" for e in errors)
            )
        for dof, value in masses.items():
            key = (node_id, dof)
            self.nodal_masses[key] = self.nodal_masses.get(key, 0.0) + float(value)

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
