# solidum_fem/solidum/core/domain.py
from __future__ import annotations

import math
import weakref
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
        self._topology_version = 0
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
        # Nodos duplicados al separar caras interiores (``split``, ADR 0043
        # §6): id original → ids de todas sus copias, el original primero.
        # Una referencia por id a cualquiera de ellos es ambigua (¿qué lado?)
        # y el lector del YAML la rechaza nombrando las copias.
        self.split_nodes: dict[int, tuple[int, ...]] = {}
        # Interfaces creadas desde la malla (``pair``/``split``, ADR 0043 §6):
        # nombre → ids de sus elementos. Las selecciona la salida por nombre
        # (``interface: <nombre>`` en la historia y en las energías por zona).
        self.interfaces: dict[str, tuple[int, ...]] = {}
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
        """Añade un elemento y lo enlaza con este dominio (``element.domain``,
        pieza P10, ADR 0042). Es, con :meth:`replace_element`, el único camino
        para cambiar ``elements``.

        Raises
        ------
        ValueError
            Si ya hay un elemento con ese id, o si el elemento pertenece a otro
            dominio vivo.
        """
        if element.id in self.elements:
            raise ValueError(f"El elemento {element.id} ya existe.")
        self._check_not_in_other_domain(element)
        self.elements[element.id] = element
        element._domain_ref = weakref.ref(self)
        self._topology_version += 1

    def replace_element(self, element) -> None:
        """Sustituye el elemento que tiene el mismo id (p. ej. para cambiar su
        material o su tipo). El nuevo queda enlazado con este dominio y el
        anterior, desenlazado (``domain`` = ``None``).

        Raises
        ------
        ValueError
            Si no hay ningún elemento con ese id, o si el nuevo pertenece a
            otro dominio vivo.
        """
        old = self.elements.get(element.id)
        if old is None:
            raise ValueError(f"replace_element: no hay ningún elemento {element.id} que sustituir.")
        self._check_not_in_other_domain(element)
        self.elements[element.id] = element
        if old is not element:
            old._domain_ref = None
        element._domain_ref = weakref.ref(self)
        self._topology_version += 1

    def reassign_element_nodes(self, element_id: int, mapping: dict[int, int]) -> None:
        """Sustituye nodos de un elemento por otros nodos del dominio en la
        **misma posición** (ADR 0043 §6): ``mapping`` lleva id viejo → id
        nuevo. Es, con :meth:`add_element` y :meth:`replace_element`, el único
        camino para cambiar la conectividad de un elemento del dominio.

        Lo usa la separación de caras interiores (``split``): las copias de un
        nodo están donde el original, así que la geometría que el elemento
        calculó al construirse sigue valiendo y su rigidez no cambia. Registra
        en los nodos nuevos los grados de libertad del elemento y aumenta
        :attr:`topology_version`.

        Raises
        ------
        ValueError
            Si no hay ningún elemento ``element_id``, si el elemento no tiene
            alguno de los nodos viejos, si algún nodo nuevo no existe o si no
            está en la misma posición que el viejo. No cambia nada si falla.
        """
        el = self.elements.get(element_id)
        if el is None:
            raise ValueError(f"reassign_element_nodes: no hay ningún elemento {element_id}.")
        ids = [n.id for n in el.nodes]
        for old_id, new_id in mapping.items():
            if old_id not in ids:
                raise ValueError(f"reassign_element_nodes: el elemento {element_id} no tiene el "
                                 f"nodo {old_id} (tiene {ids}).")
            new = self.nodes.get(new_id)
            if new is None:
                raise ValueError(f"reassign_element_nodes: no hay ningún nodo {new_id}.")
            old = el.nodes[ids.index(old_id)]
            if tuple(map(float, new.coordinates)) != tuple(map(float, old.coordinates)):
                raise ValueError(
                    f"reassign_element_nodes: el nodo {new_id} {list(new.coordinates)} no está "
                    f"en la posición del nodo {old_id} {list(old.coordinates)}; el elemento "
                    f"{element_id} guarda la geometría que calculó al construirse.")
        new_ids = list(mapping.values())
        if len(set(new_ids)) != len(new_ids):
            raise ValueError(f"reassign_element_nodes: el mapa {mapping} lleva dos nodos viejos "
                             f"al mismo nuevo; el elemento {element_id} quedaría con un nodo "
                             f"repetido.")
        repeated = sorted(set(new_ids) & (set(ids) - set(mapping)))
        if repeated:
            raise ValueError(f"reassign_element_nodes: el elemento {element_id} ya tiene los "
                             f"nodos {repeated}; quedarían repetidos.")
        nodes = [self.nodes[mapping[n.id]] if n.id in mapping else n for n in el.nodes]
        for node in nodes:
            for dof_name in el.DOF_NAMES:
                node.add_dof(dof_name)
        el.nodes = nodes
        self._topology_version += 1

    def record_split_nodes(self, copies: dict[int, tuple[int, ...]]) -> None:
        """Anota los nodos que duplicó ``split`` (ADR 0043 §6): id original →
        ids de todas sus copias, con el original primero."""
        self.split_nodes.update(copies)
        origin = self.__dict__.setdefault("_split_origin", {})
        for p, ids in copies.items():
            for c in ids:
                origin[c] = p

    def split_copies(self, node_id: int) -> tuple[int, ...] | None:
        """Ids de todas las copias si ``node_id`` es uno de los nodos que
        duplicó ``split`` (el original o una copia); si no, ``None``."""
        origin = self.__dict__.get("_split_origin", {})
        if len(origin) != sum(len(ids) for ids in self.split_nodes.values()):
            # split_nodes es público: si se escribió sin record_split_nodes,
            # el mapa inverso se rehace desde él (la fuente es split_nodes).
            origin = {c: p for p, ids in self.split_nodes.items() for c in ids}
            self.__dict__["_split_origin"] = origin
        p = origin.get(node_id)
        return None if p is None else self.split_nodes[p]

    @property
    def topology_version(self) -> int:
        """Contador que aumenta con cada ``add_element``, ``replace_element`` y
        ``reassign_element_nodes`` (ADR 0042 y 0043): huella barata para los
        datos que un componente deriva de los elementos del dominio (p. ej. el
        contorno libre) y guarda."""
        return self._topology_version

    def _check_not_in_other_domain(self, element) -> None:
        other = element.domain
        if other is not None and other is not self:
            raise ValueError(
                f"El elemento {element.id} ya pertenece a otro dominio; un elemento sólo puede "
                f"estar en un dominio (cree uno nuevo, o use copy.deepcopy del dominio entero).")

    def __setstate__(self, state):
        """Al deserializar o copiar con ``copy.deepcopy``, los elementos vuelven
        a enlazarse con este dominio (no se serializa la referencia, P10).

        Sólo los que no tienen dominio vivo: una copia superficial
        (``copy.copy``) comparte los elementos con el original, que los
        conserva."""
        self.__dict__.update(state)
        ref = weakref.ref(self)
        for element in self.elements.values():
            if element.domain is None:
                element._domain_ref = ref

    def generate_equation_numbers(self, verbose: bool = False):
        eq_number = 0
        for node_id, node in self.nodes.items():
            for dof_name in node.dofs.keys():
                node.dofs[dof_name] = eq_number
                eq_number += 1
        self.total_dofs = eq_number
        if verbose:
            _log.info(f"Numeración completada. Grados de libertad totales: {self.total_dofs}")

    def renumber_if_numbered(self) -> None:
        """Rehace la numeración si el dominio ya estaba numerado. La llaman los
        generadores que añaden nodos o reasignan elementos sobre un dominio ya
        construido (refuerzo embebido, ADR 0032; ``split``, ADR 0043 §6): sin
        ella, los grados de libertad nuevos quedarían con ecuación −1."""
        if self.total_dofs:
            self.generate_equation_numbers()
