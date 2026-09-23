# solidum_fem/solidum/math/assembly.py
from typing import Callable

import numpy as np
import scipy.sparse as sp

from solidum.bc.constraints import ConstraintSet
from solidum.constants import BATCH_ASSEMBLY_DEFAULT, BATCH_PARALLEL_DEFAULT
from solidum.core.domain import Domain
from solidum.core.element import Element
from solidum.logging import get_logger
from solidum.math.batch.family import Family, build_families
from solidum.math.batch.kernels import coo_to_csr_reduce, coo_to_csr_reduce_parallel

_log = get_logger("assembly")


class Assembler:
    """Ensambla las matrices y vectores globales del dominio.

    Dos caminos coexisten (ADR 0014):

    - **por lotes** (``batch=True``, default de ``BATCH_ASSEMBLY_DEFAULT``):
      los elementos cuya clase declara ``BATCH_KINEMATICS`` y cuyo material
      declara ``STATE_SCHEMA`` + kernel puntual se agrupan en *familias*
      (misma clase, misma instancia de material, misma cuadratura, mismo
      número de nodos) y se evalúan dentro de un único kernel compilado;
      su estado interno vive en arreglos (``FamilyState``) y ``elem.state``
      pasa a ser una vista sobre ellos.
    - **por elemento**: el bucle de siempre sobre ``compute_element_state``.
      Lo siguen los elementos sin kernel (1D, embebido, materiales sin
      esquema) y todo el modelo cuando ``batch=False``.

    Ambos caminos ejecutan las mismas funciones compiladas por punto de
    Gauss y producen ``K`` y ``F_int`` idénticos a precisión de máquina; lo
    garantiza ``tests/test_batch_assembly.py`` sobre los registros
    completos. Los solvers no distinguen el camino: consumen
    ``assemble_non_linear_system``, ``commit_all_states`` y ``reduce`` como
    siempre.

    Parameters
    ----------
    domain
        Dominio con nodos y elementos.
    batch
        ``True`` / ``False`` para forzar el camino; ``None`` lee
        ``solidum.constants.BATCH_ASSEMBLY_DEFAULT``.
    batch_memory_budget
        Presupuesto (bytes) para los temporales de un trozo de familia;
        ``None`` ⇒ ``BATCH_MEMORY_BUDGET_BYTES``.
    parallel
        Variante del kernel de familia: ``True`` reparte los elementos de
        cada trozo entre los hilos de Numba (``prange``), ``False`` usa la
        variante serie; ``None`` lee ``BATCH_PARALLEL_DEFAULT``. El
        resultado es bit a bit el mismo. El número de hilos lo fija Numba
        (``numba.set_num_threads`` / ``NUMBA_NUM_THREADS``).
    """

    def __init__(self, domain: Domain, batch: bool | None = None,
                 batch_memory_budget: int | None = None,
                 parallel: bool | None = None):
        self.domain = domain
        self.ndof = 0
        self.K_global = None
        self.F_global = None

        self.batch = BATCH_ASSEMBLY_DEFAULT if batch is None else bool(batch)
        self.batch_memory_budget = batch_memory_budget
        self.parallel = BATCH_PARALLEL_DEFAULT if parallel is None else bool(parallel)

        # Variables para caché de topología COO
        self._topology_built = False
        self._coo_rows = None
        self._coo_cols = None
        self._total_entries = 0
        self._elem_dof_indices = []
        # Posición del bloque de cada elemento (orden del dominio) en `data`.
        self._elem_ptr = np.zeros(0, dtype=np.int64)
        # Familias por lotes y posiciones (orden del dominio) de los
        # elementos que siguen el camino por elemento.
        self._families: list[Family] = []
        self._loose: list[int] = []

        # Caché del ConstraintSet (ADR 0004 fase 1).
        self._constraint_set: ConstraintSet | None = None

        # Casi-núcleo del sistema reducido para el precondicionador AMG
        # (ADR 0018), con la huella (topología, restricciones) con que se
        # calculó.
        self._near_nullspace: np.ndarray | None = None
        self._near_nullspace_key: tuple | None = None

        # Caché de la matriz de masa global (ADR 0009). M es lineal y no
        # cambia entre llamadas: se ensambla una vez por análisis y se reusa.
        self._M_global: sp.csr_matrix | None = None
        self._M_lumping: str | None = None

        # Huellas para detectar cambios en el modelo entre llamadas: sin
        # ellas un apoyo añadido o un elemento nuevo tras el primer
        # ensamblaje se ignoraba en silencio (auditoría 2026-09-22).
        self._topology_key: tuple | None = None
        self._constraint_fp: int | None = None

        # Cargas nodales equivalentes por elemento (ejes globales) del último
        # ``assemble_body_load`` / ``assemble_self_weight``. Las consume
        # ``build_solve_result`` para que ``internal_forces`` devuelva las
        # fuerzas internas de extremo (``F_int − f_eq``) y no ``K·u``.
        self.equivalent_loads: dict[int, np.ndarray] = {}

    def _get_element_global_indices(self, element) -> list:
        """Extrae los índices globales de DOFs del elemento en el orden correcto.

        Delega en element.get_global_dof_indices(), que itera sobre DOF_NAMES × nodos
        sin asumir uniformidad de DOFs entre nodos.
        """
        return element.get_global_dof_indices()

    def _build_topology(self):
        """Precalcula y cachea la topología de ensamblaje (filas, columnas y
        mapeo) una sola vez, y construye las familias por lotes.

        El vector COO ``data`` se ordena **por familia**: cada familia ocupa
        un bloque contiguo (sus elementos en orden del dominio) y los
        elementos del camino por elemento van detrás. La conversión a CSR
        suma duplicados con independencia del orden, así que el patrón de
        ``K`` no cambia; lo que se gana es volcar cada trozo de familia con
        una asignación de rodaja en vez de un índice disperso.
        """
        if self.domain.total_dofs == 0:
            self.domain.generate_equation_numbers()

        self.ndof = self.domain.total_dofs
        elements = list(self.domain.elements.values())
        self._elem_dof_indices = [self._get_element_global_indices(e) for e in elements]
        self._total_entries = sum(len(idx) ** 2 for idx in self._elem_dof_indices)

        self._release_families()
        if self.batch:
            self._families, self._loose = build_families(
                elements, self._elem_dof_indices, self.batch_memory_budget,
                parallel=self.parallel)
        else:
            self._families, self._loose = [], list(range(len(elements)))

        self._coo_rows = np.zeros(self._total_entries, dtype=np.int32)
        self._coo_cols = np.zeros(self._total_entries, dtype=np.int32)
        self._elem_ptr = np.zeros(len(elements), dtype=np.int64)

        ptr = 0

        def _place(pos: int) -> None:
            nonlocal ptr
            global_indices = self._elem_dof_indices[pos]
            n_idx = len(global_indices)
            self._coo_rows[ptr:ptr + n_idx ** 2] = np.repeat(global_indices, n_idx)
            self._coo_cols[ptr:ptr + n_idx ** 2] = np.tile(global_indices, n_idx)
            self._elem_ptr[pos] = ptr
            ptr += n_idx ** 2

        for fam in self._families:
            fam.ptr_start = ptr
            for pos in fam.positions:
                _place(int(pos))
            fam.adopt_states()
        for pos in self._loose:
            _place(pos)

        self._build_csr_map()

        self._topology_built = True
        self._topology_key = self._current_topology_key()
        # La masa cacheada corresponde a la topología anterior.
        self._M_global = None
        self._M_lumping = None

    def _build_csr_map(self) -> None:
        """Estructura CSR de ``K`` y mapa inverso COO → CSR, calculados una vez.

        ``scipy.sparse.coo_matrix(...).tocsr()`` ordena y suma duplicados en
        **cada** ensamblaje (unos 115 ms para 8 000 Hex8, más que el propio
        kernel elástico por lotes). Con el patrón fijo, basta calcular una
        vez qué entradas COO alimentan cada entrada CSR (``_csr_src_ptr``,
        ``_csr_src_idx``: la lista de fuentes de cada posición CSR, en orden
        creciente de índice COO) y reducir con un kernel compilado, serie o
        paralelo según ``parallel``. La suma por posición sigue el mismo
        orden que ``np.bincount``, así que el resultado es idéntico bit a
        bit; la diferencia es que cada posición CSR se calcula de forma
        independiente y se reparte entre hilos. Las filas y columnas COO se
        descartan después: el mapa inverso (4 B por entrada más 8 B por
        entrada CSR) ocupa menos que ellas.
        """
        n = self._total_entries
        ndof = self.ndof
        rows = self._coo_rows
        cols = self._coo_cols
        pattern = sp.csr_matrix(
            (np.ones(n, dtype=np.float64), (rows, cols)), shape=(ndof, ndof)
        )
        pattern.sum_duplicates()
        pattern.sort_indices()
        self._csr_indptr = pattern.indptr.astype(np.int32, copy=False)
        self._csr_indices = pattern.indices.astype(np.int32, copy=False)
        self._nnz = int(self._csr_indices.shape[0])
        csr_keys = (np.repeat(np.arange(ndof, dtype=np.int64), np.diff(self._csr_indptr))
                    * ndof + self._csr_indices)
        coo_keys = rows.astype(np.int64) * ndof + cols
        csr_map = np.searchsorted(csr_keys, coo_keys)
        # Mapa inverso: fuentes COO de cada posición CSR, estables en el
        # orden COO (mismo orden de acumulación que ``np.bincount``).
        order = np.argsort(csr_map, kind="stable")
        counts = np.bincount(csr_map, minlength=self._nnz)
        self._csr_src_ptr = np.zeros(self._nnz + 1, dtype=np.int64)
        np.cumsum(counts, out=self._csr_src_ptr[1:])
        self._csr_src_idx = np.ascontiguousarray(order, dtype=np.int32)
        self._coo_rows = None
        self._coo_cols = None
        # Vector COO reutilizado en cada ensamblaje: todas sus entradas se
        # sobreescriben (bloques completos por familia y por elemento), así
        # que no hace falta ponerlo a cero.
        self._data = np.empty(n, dtype=np.float64)

    def _to_csr(self, data: np.ndarray) -> sp.csr_matrix:
        """Matriz CSR a partir del vector de entradas COO (patrón cacheado)."""
        values = np.empty(self._nnz, dtype=np.float64)
        reduce = coo_to_csr_reduce_parallel if self.parallel else coo_to_csr_reduce
        reduce(np.ascontiguousarray(data, dtype=np.float64),
               self._csr_src_ptr, self._csr_src_idx, values)
        return sp.csr_matrix(
            (values, self._csr_indices, self._csr_indptr),
            shape=(self.ndof, self.ndof),
        )

    def _release_families(self) -> None:
        """Devuelve el estado de las familias a los ``ElementState`` clásicos
        y olvida las familias (antes de reconstruir la topología)."""
        for fam in self._families:
            fam.release()
        self._families = []
        self._loose = []

    @property
    def families(self) -> list[Family]:
        """Familias por lotes construidas con la topología vigente (vacía si
        ``batch=False`` o si aún no se ha ensamblado)."""
        return list(self._families)

    def _current_topology_key(self) -> tuple:
        return (len(self.domain.elements), self.domain.total_dofs, self.batch, self.parallel)

    def _ensure_topology(self) -> None:
        """Construye la topología COO si falta o si el modelo cambió de
        tamaño (elementos añadidos, renumeración) desde la última vez."""
        if not self._topology_built or self._topology_key != self._current_topology_key():
            self._build_topology()

    def invalidate(self) -> None:
        """Descarta todas las cachés (topología, familias, restricciones, masa).

        Las huellas de topología y de restricciones ya detectan por sí solas
        elementos nuevos y apoyos/MPC añadidos o cambiados; este método es
        el punto explícito para cambios que no dejan huella (densidad de un
        material, un material o un espesor sustituidos en un elemento
        existente, ...). El estado interno de las familias vuelve a los
        elementos, así que no se pierde historia."""
        self._release_families()
        self._topology_built = False
        self._topology_key = None
        self._constraint_set = None
        self._constraint_fp = None
        self._M_global = None
        self._M_lumping = None

    def _elements(self) -> list:
        return list(self.domain.elements.values())

    def assemble_system(self):
        """Matriz de rigidez global en la configuración de referencia
        (``u = 0``): ``compute_global_stiffness`` en el camino por elemento y
        el kernel de familia con ``U = 0`` en el camino por lotes. En ambos
        casos el estado *trial* **no se modifica**: es una evaluación
        auxiliar y el trial previo se restaura al salir."""
        self._ensure_topology()

        self.F_global = np.zeros(self.ndof)
        data = self._data

        if self._families:
            U0 = np.zeros(self.ndof)
            F_discard = np.zeros(self.ndof)
            for fam in self._families:
                fam.evaluate(U0, data, F_discard, keep_trial=True)

        elements = self._elements()
        for pos in self._loose:
            K_e = elements[pos].compute_global_stiffness()
            n_idx = len(self._elem_dof_indices[pos])
            p0 = self._elem_ptr[pos]
            data[p0:p0 + n_idx ** 2] = K_e.ravel()

        # CSR: eficiente para solve y para modificar entradas diagonales existentes
        self.K_global = self._to_csr(data)

    def assemble_non_linear_system(self, U_current: np.ndarray):
        """Construye la matriz tangente global y el vector de fuerzas internas a máxima velocidad."""
        self._ensure_topology()

        F_int_global = np.zeros(self.ndof)
        data = self._data
        U_current = np.asarray(U_current, dtype=np.float64)

        for fam in self._families:
            fam.evaluate(U_current, data, F_int_global)

        elements = self._elements()
        for pos in self._loose:
            element = elements[pos]
            global_indices = self._elem_dof_indices[pos]
            u_local = U_current[global_indices]

            K_e, F_int_e = element.compute_element_state(u_local)

            F_int_global[global_indices] += F_int_e

            n_idx = len(global_indices)
            p0 = self._elem_ptr[pos]
            data[p0:p0 + n_idx ** 2] = K_e.ravel()

        return self._to_csr(data), F_int_global

    # ------------------------------------------------------------------
    # Imposición de Dirichlet por eliminación directa (ADR 0004 fase 1).
    # ------------------------------------------------------------------

    def _constraint_fingerprint(self) -> int:
        """Huella de las restricciones declaradas en el dominio (Dirichlet y
        MPC). O(n_restricciones) por llamada; despreciable frente a un
        ensamblaje y evita servir un ``ConstraintSet`` obsoleto."""
        dirichlet = tuple(
            (nid, dof, float(val))
            for nid, node in self.domain.nodes.items()
            for dof, val in node.boundary_conditions.items()
        )
        linear = tuple(
            (tuple(spec['slave']), tuple(tuple(m) for m in spec['masters']),
             tuple(float(c) for c in spec['coefficients']), float(spec.get('g', 0.0)))
            for spec in getattr(self.domain, 'linear_constraints', [])
        )
        return hash((dirichlet, linear, self.domain.total_dofs))

    @property
    def constraint_set(self) -> ConstraintSet:
        """ConstraintSet derivado de ``Node.boundary_conditions`` y de las
        restricciones lineales del dominio. Se reconstruye automáticamente si
        cambian (huella), así que añadir un apoyo tras el primer ``solve``
        surte efecto en el siguiente."""
        fp = self._constraint_fingerprint()
        if self._constraint_set is None or fp != self._constraint_fp:
            self._constraint_set = self._build_constraint_set()
            self._constraint_fp = fp
        return self._constraint_set

    def _build_constraint_set(self) -> ConstraintSet:
        self._ensure_topology()
        cs = ConstraintSet()
        for node in self.domain.nodes.values():
            for dof_name, value in node.boundary_conditions.items():
                cs.add_dirichlet(node.dofs[dof_name], value)
        for spec in getattr(self.domain, "linear_constraints", []):
            slave_node_id, slave_dof_name = spec["slave"]
            slave_node = self.domain.get_node(slave_node_id)
            if slave_node is None or slave_dof_name not in slave_node.dofs:
                raise ValueError(
                    f"Restricción lineal: nodo {slave_node_id} o DOF "
                    f"'{slave_dof_name}' inexistente."
                )
            slave = slave_node.dofs[slave_dof_name]
            masters: list[int] = []
            for nid, dn in spec["masters"]:
                node = self.domain.get_node(nid)
                if node is None or dn not in node.dofs:
                    raise ValueError(
                        f"Restricción lineal: nodo maestro {nid} o DOF "
                        f"'{dn}' inexistente."
                    )
                masters.append(node.dofs[dn])
            cs.add_linear(slave, masters, spec["coefficients"], spec.get("g", 0.0))
        return cs

    def reduce(
        self,
        K: sp.spmatrix,
        F: np.ndarray,
        U_current: np.ndarray | None = None,
        load_factor: float = 1.0,
    ):
        """Reduce ``K, F`` a la subred de DOFs libres (eliminación directa).

        Modelo afín ``u = T·u_libre + g`` con ``T`` sparse de shape
        ``(ndof, n_libre)`` que selecciona DOFs libres y, en filas esclavas,
        lleva los coeficientes ``α_si`` de las restricciones lineales. ``g``
        es no nulo solo en filas esclavas con término independiente.

        Para sistemas incrementales (Newton, corrector arc-length) se
        interpreta el resultado como ``δu = T·δu_libre + g_inc`` con
        ``g_inc = T·U_current[free] + load_factor·g_indep − U_current``,
        que en DOFs libres se anula y en DOFs esclavos representa la
        corrección necesaria para satisfacer la restricción.

        Returns
        -------
        K_red, F_red, T, g
            ``K_red = TᵀKT``, ``F_red = Tᵀ(F − K·g)``, el operador ``T`` y
            el vector ``g``. La reconstrucción se hace con :meth:`expand`.
        """
        cs = self.constraint_set
        T, g_indep = cs.build(self.ndof)

        if U_current is None:
            g = load_factor * g_indep
        else:
            free_dofs = cs.free_dofs(self.ndof)
            u_free = U_current[free_dofs]
            g = T @ u_free + load_factor * g_indep - U_current

        K_red = T.T @ K @ T
        F_red = T.T @ (F - K @ g)
        return K_red, F_red, T, g

    def expand(
        self,
        u_red: np.ndarray,
        T: sp.spmatrix,
        g: np.ndarray,
    ) -> np.ndarray:
        """Reconstruye ``u`` completo: ``u = T·u_red + g``."""
        return T @ u_red + g

    # ------------------------------------------------------------------
    # Matriz de masa global (ADR 0009)
    # ------------------------------------------------------------------

    def assemble_mass_matrix(self, lumping: str = "consistent") -> sp.csr_matrix:
        """Ensambla la matriz de masa global ``M`` (ADR 0009).

        Reutiliza la topología cacheada por ``_build_topology`` (posición de
        cada elemento en el vector de entradas y mapa COO → CSR) — ``M`` y
        ``K`` comparten el mismo patrón de sparsity porque se ensamblan sobre
        los mismos pares de DOFs elemento a elemento. El coste de ensamblaje
        de ``M`` es comparable al de ``K`` por elemento.

        La masa lineal es constante en el tiempo y se cachea por análisis:
        llamadas sucesivas con el mismo ``lumping`` devuelven el resultado
        previo sin recomputar.

        Pre-validación agregada (espíritu de ADR 0008 y del YAML parser):

        - Materiales que no declaran ``density`` se acumulan y reportan
          juntos con ``ValueError``.
        - Elementos cuya subclase no override ``compute_mass_matrix`` (la
          base lanza ``NotImplementedError``) se reportan en la misma
          excepción.

        Parameters
        ----------
        lumping : {"consistent"}, default "consistent"
            Estrategia de discretización de la inercia. Se propaga literal
            a cada elemento; cada subclase valida los valores que admite.

        Returns
        -------
        scipy.sparse.csr_matrix, shape (ndof, ndof)
            Matriz de masa global, simétrica, positiva (semi)definida.

        Raises
        ------
        ValueError
            Si algún material no declara ``density`` o algún elemento no
            implementa ``compute_mass_matrix``. El mensaje agrupa ambos
            tipos de problema.
        """
        self._ensure_topology()

        if self._M_global is not None and self._M_lumping == lumping:
            return self._M_global

        missing_density: list[str] = []
        missing_method: list[str] = []
        for element in self.domain.elements.values():
            mat = element.material
            if mat is None or getattr(mat, "density", None) is None:
                # Un material que sabe explicar su propio requisito lo hace
                # él mismo. El mensaje genérico de abajo está escrito para el
                # dominio mecánico —donde `density = 0.0` es legítimo en un
                # material sin masa por diseño (penalty, restricción)— y ese
                # consejo es FÍSICAMENTE INCORRECTO en otras familias: en un
                # material térmico produciría `ρc = 0`, capacidad calorífica
                # nula, que hace singular la matriz de capacidad y deja el
                # transitorio irresoluble. Delegar preserva el diagnóstico
                # accionable de cada familia sin que el ensamblador tenga que
                # conocerlas (Reglas.md §1: el coste del componente N+1).
                explicar = getattr(mat, "volumetric_capacity", None)
                if callable(explicar):
                    explicar(consumer="el ensamblaje de la matriz de capacidad")
                missing_density.append(
                    type(mat).__name__ if mat is not None else "<sin material>"
                )
            if type(element).compute_mass_matrix is Element.compute_mass_matrix:
                missing_method.append(type(element).__name__)

        if missing_density or missing_method:
            parts = ["assemble_mass_matrix: el cálculo no puede completarse."]
            if missing_density:
                parts.append(
                    "- Materiales sin `density` declarada (ADR 0008): "
                    f"{sorted(set(missing_density))}. Añade el atributo "
                    "`density` (kg/m³ en las unidades del problema) en cada "
                    "uno; usa 0.0 explícitamente si el material es sin masa "
                    "por diseño (penalty, restricción)."
                )
            if missing_method:
                parts.append(
                    "- Elementos sin `compute_mass_matrix` implementado "
                    f"(ADR 0009): {sorted(set(missing_method))}. La subclase "
                    "debe sobreescribir el método para participar en análisis "
                    "modal o dinámico."
                )
            raise ValueError("\n".join(parts))

        data = self._data
        for i, element in enumerate(self.domain.elements.values()):
            M_e = element.compute_mass_matrix(lumping=lumping)
            n_idx = len(self._elem_dof_indices[i])
            p0 = self._elem_ptr[i]
            data[p0:p0 + n_idx ** 2] = M_e.ravel()

        self._M_global = self._to_csr(data)
        self._M_lumping = lumping
        return self._M_global

    def near_nullspace(self) -> np.ndarray:
        """Modos de cuerpo rígido restringidos a los DOF libres (ADR 0018).

        Casi-núcleo del sistema reducido ``K_red = TᵀKT`` para el
        precondicionador AMG del backend iterativo. Se derivan de
        coordenadas y nombres de DOF (:func:`rigid_body_modes`) y se
        restringen a los DOF libres en el orden de las columnas de ``T``:
        con restricciones afines el valor en un esclavo lo reconstruye
        ``T`` a partir de sus maestros, así que basta la restricción.

        Se calcula sólo si alguien lo pide —los solvers lo pasan como
        proveedor perezoso en ``StiffnessProperties.near_nullspace`` y sólo
        AMG lo invoca— y se cachea mientras no cambien la topología ni las
        restricciones.
        """
        from solidum.math.linalg.nullspace import rigid_body_modes

        self._ensure_topology()
        cs = self.constraint_set
        key = (self._topology_key, self._constraint_fp)
        if self._near_nullspace is None or self._near_nullspace_key != key:
            B = rigid_body_modes(self.domain)
            free = cs.free_dofs(self.ndof)
            self._near_nullspace = np.ascontiguousarray(B[free])
            self._near_nullspace_key = key
        return self._near_nullspace

    def reduce_pair(
        self,
        K: sp.spmatrix,
        M: sp.spmatrix,
    ) -> tuple[sp.spmatrix, sp.spmatrix, sp.spmatrix]:
        """Reduce simultáneamente ``K`` y ``M`` por eliminación directa.

        Variante de :meth:`reduce` para el problema generalizado modal
        ``K·φ = ω²·M·φ`` (ADR 0009 §5). Los términos ``g_indep`` de
        restricciones lineales no homogéneas no aplican: los modos son
        solución del problema homogéneo asociado, así que el operador ``T``
        de selección/coeficientes basta para mapear DOFs libres a globales.

        Returns
        -------
        K_red, M_red, T
            Matrices reducidas ``K_red = TᵀKT`` y ``M_red = TᵀMT``, junto al
            operador ``T``. La expansión de modos al espacio completo se
            obtiene como ``φ = T·φ_red`` (en DOFs prescritos el modo es 0).
        """
        cs = self.constraint_set
        T, _g_indep = cs.build(self.ndof)
        K_red = T.T @ K @ T
        M_red = T.T @ M @ T
        return K_red, M_red, T

    def commit_all_states(self):
        """Confirma las variables internas de todos los elementos tras la
        convergencia del paso: copia trial → committed por familia (una
        asignación de arreglo) y ``commit_state`` en los elementos del
        camino por elemento."""
        if not self._topology_built or self._topology_key != self._current_topology_key():
            for elem in self.domain.elements.values():
                elem.commit_state()
            return
        for fam in self._families:
            fam.commit()
        elements = self._elements()
        for pos in self._loose:
            elements[pos].commit_state()

    def prepare_all_steps(self, U_committed: np.ndarray) -> None:
        """Invoca ``prepare_step(U_committed)`` en todos los elementos del dominio.

        Hook de preparación de paso (ADR 0010 §5). Lo llaman los solvers no
        lineales una vez por paso, con el campo convergido del paso anterior,
        antes del primer ensamblaje del Newton. La implementación base de
        ``Element.prepare_step`` es no-op; los elementos con discontinuidad
        embebida lo sobreescriben para chequear activación.
        """
        for elem in self.domain.elements.values():
            elem.prepare_step(U_committed)

    def apply_point_load(self, node_id: int, dof_name: str, value: float):
        node = self.domain.get_node(node_id)
        if node and dof_name in node.dofs:
            idx = node.dofs[dof_name]
            self.F_global[idx] += value

    @staticmethod
    def _promote_to_3d(v, name: str) -> np.ndarray:
        """Acepta vector 2D ó 3D y promueve a R³ con padding cero si es 2D."""
        arr = np.asarray(v, dtype=float).ravel()
        if arr.size == 2:
            return np.array([arr[0], arr[1], 0.0])
        if arr.size == 3:
            return arr
        raise ValueError(
            f"{name}: dimensión inesperada ({arr.size}). Esperado 2 ó 3 componentes."
        )

    def _iterate_body_force(self, get_b_for_element: Callable) -> np.ndarray:
        """Recorrido interno común a `assemble_body_load` y `assemble_self_weight`.

        Itera sobre todos los elementos que declaran `compute_body_load` y
        delega en `get_b_for_element(element) -> np.ndarray` la obtención del
        vector de fuerza de cuerpo apropiado para cada uno (en sus
        componentes 2D ó 3D según `'uz' in element.DOF_NAMES`).
        """
        self._ensure_topology()

        F_body = np.zeros(self.ndof)
        equivalent_loads: dict[int, np.ndarray] = {}
        for i, element in enumerate(self.domain.elements.values()):
            method = getattr(element, "compute_body_load", None)
            if method is None:
                continue
            b_elem = get_b_for_element(element)
            f_e = np.asarray(method(b_elem), dtype=float)
            F_body[self._elem_dof_indices[i]] += f_e
            equivalent_loads[element.id] = f_e
        self.equivalent_loads = equivalent_loads
        return F_body

    def assemble_body_load(self, b) -> np.ndarray:
        """Acumula la integral consistente de fuerza de cuerpo uniforme.

        Aplica el mismo vector `b` (fuerza por unidad de volumen, en ejes
        globales) a todos los elementos. Útil cuando el usuario precalcula
        `ρ·g` para una estructura monomaterial, o cuando la fuerza de cuerpo
        no proviene de gravedad (campo electromagnético, centrífuga, etc.).

        Para peso propio en estructuras multimaterial, preferir
        :meth:`assemble_self_weight` que usa `material.density` por elemento
        (ADR 0008).

        Promoción 2D↔3D: se acepta `b` con 2 ó 3 componentes; cada elemento
        recibe el slice apropiado según su dimensionalidad
        (`'uz' in DOF_NAMES` ⇒ 3D, en otro caso 2D). La omisión silenciosa
        de elementos sin `compute_body_load` es deliberada — el método es
        opcional, no parte del contrato base de `Element`.

        Parameters
        ----------
        b : np.ndarray or list, shape (2,) or (3,)
            Fuerza de cuerpo por unidad de volumen en ejes globales.

        Returns
        -------
        np.ndarray, shape (ndof,)
            Vector global de fuerzas nodales consistentes con `∫NᵀbA dx`.
        """
        b3 = self._promote_to_3d(b, "assemble_body_load")

        def _b_for_element(element):
            return b3 if 'uz' in element.DOF_NAMES else b3[:2]

        return self._iterate_body_force(_b_for_element)

    def assemble_self_weight(self, g) -> np.ndarray:
        """Acumula peso propio usando `material.density` de cada elemento (ADR 0008).

        Cada elemento ve `b_element = material.density · g`. Comportamiento
        según el valor de `density` del material asignado:

        - `density is None` (no declarado al construir el material): el
          cálculo no puede completarse físicamente — falta información de
          masa. La función falla con `ValueError` listando los materiales
          afectados por nombre. Sin posibilidad de fallo silencioso.
        - `density == 0.0` declarado explícitamente: caso legítimo de
          material sin masa (penalty, restricción, fixture). Aporte nulo
          al vector global; emite `WARNING` informativo.
        - `density > 0`: aporte normal `ρ·g·integral consistente`.

        Promoción 2D↔3D: idéntica a `assemble_body_load`.

        Parameters
        ----------
        g : np.ndarray or list, shape (2,) or (3,)
            Vector aceleración gravitatoria en ejes globales (m/s² o las
            unidades correspondientes del problema). Típicamente
            `[0, -9.81]` en 2D ó `[0, 0, -9.81]` en 3D.

        Returns
        -------
        np.ndarray, shape (ndof,)
            Vector global de fuerzas nodales consistentes con peso propio.

        Raises
        ------
        ValueError
            Si algún material del modelo no declara `density` y se invoca
            este cálculo.
        """
        g3 = self._promote_to_3d(g, "assemble_self_weight")
        warned_materials: set[int] = set()

        missing: list[str] = []

        def _b_for_element(element):
            density = element.material.density
            if density is None:
                missing.append(type(element.material).__name__)
                return 0.0 * (g3 if 'uz' in element.DOF_NAMES else g3[:2])  # placeholder; falla más abajo
            if density == 0.0 and id(element.material) not in warned_materials:
                warned_materials.add(id(element.material))
                _log.warning(
                    f"assemble_self_weight: material '{type(element.material).__name__}' "
                    f"declara density=0.0 explícitamente; su aporte al peso propio "
                    f"es nulo. Caso legítimo si el material representa restricción/penalty "
                    f"sin masa física; revísalo si esperabas peso propio aquí."
                )
            g_elem = g3 if 'uz' in element.DOF_NAMES else g3[:2]
            return density * g_elem

        F = self._iterate_body_force(_b_for_element)
        if missing:
            unique = sorted(set(missing))
            raise ValueError(
                "assemble_self_weight: los siguientes materiales no declaran "
                "density y son requeridos por el cálculo de peso propio: "
                f"{unique}. Añade el atributo `density` (kg/m³ en las unidades "
                "del problema) en cada uno de ellos. Si el material no tiene "
                "masa física (penalty/restricción), declarar `density=0.0` "
                "explícitamente para silenciar este error."
            )
        return F
