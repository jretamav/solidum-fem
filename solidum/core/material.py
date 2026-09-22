# solidum_fem/solidum/core/material.py
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Callable, ClassVar

import numpy as np

from solidum.constants import ADMISSIBILITY_TOL_ABS, ADMISSIBILITY_TOL_REL


class BatchMaterialHooks:
    """Contrato opcional del camino por lotes (ADR 0014), común a las
    familias de materiales (mecánicos y térmicos).

    STATE_SCHEMA : ClassVar[dict | None]
        Variables internas como ``{nombre: forma}`` en el orden en que se
        almacenan por filas: ``{"eps_p": (4,), "alpha": ()}``. ``{}`` para
        un material sin historia. **Obligatorio en todo material del
        catálogo** (lo exige el barrido de contratos): con él, el estado
        interno se guarda como arreglos y el commit deja de copiar
        diccionarios, tenga o no el material kernel por lotes. ``None``
        (el default de la base) significa "no declarado" y excluye al
        material del camino por lotes sin más consecuencia.

    BATCH_KERNEL : ClassVar[callable | None]
        Función ``@njit`` con la firma ``MAT_SIG`` de
        ``solidum.math.batch.signatures``. Opcional. Un material cuya
        constitutiva dependa de la instancia (hipótesis plane stress /
        plane strain) sobreescribe :meth:`batch_kernel` en vez del
        atributo.

    Los tres métodos siguientes alimentan al kernel con constantes de la
    familia: :meth:`batch_params` (escalares en el orden que el kernel
    documenta), :meth:`batch_matrix` (la matriz ``C`` del material) e
    :meth:`initial_state` (valores iniciales de las variables internas,
    p. ej. ``κ = κ₀`` en daño). :meth:`batch_report` recibe el número de
    puntos que el kernel marcó con ``flag = 1`` tras cada evaluación.
    """

    STATE_SCHEMA: ClassVar[dict | None] = None
    BATCH_KERNEL: ClassVar[Callable | None] = None

    def batch_kernel(self) -> Callable | None:
        return type(self).BATCH_KERNEL

    def batch_params(self) -> np.ndarray:
        return np.zeros(0)

    def batch_matrix(self) -> np.ndarray:
        raise NotImplementedError(
            f"{type(self).__name__} declara un kernel por lotes pero no "
            f"implementa batch_matrix()."
        )

    def initial_state(self) -> dict | None:
        """Estado inicial por punto de Gauss: ceros para todas las variables
        del esquema (``None`` si no hay)."""
        schema = self.STATE_SCHEMA
        if not schema:
            return None
        out = {}
        for name, shape in schema.items():
            shape = tuple(shape)
            out[name] = 0.0 if shape == () else np.zeros(shape)
        return out

    def batch_report(self, n_flagged: int) -> None:
        """Hook tras cada evaluación por lotes con puntos marcados."""
        return None

    def batch_out_of_plane_stress(self, sigma: np.ndarray, S: np.ndarray) -> np.ndarray:
        """``σ_zz`` de ``P`` puntos a la vez: ``(P,)`` a partir de los
        esfuerzos ``sigma`` ``(P, n_sigma)`` y de las filas de estado
        committed ``S`` ``(P, n_state)``. Es la versión por arreglos de
        :meth:`Material.out_of_plane_stress` para el post-proceso por
        lotes (exportador VTK).

        Por defecto: ceros si el material no sobreescribe
        ``out_of_plane_stress`` (plane stress, 1D, 3D, térmicos); si lo
        sobreescribe y no aporta versión vectorizada, se evalúa punto a
        punto desempaquetando cada fila (correcto, más lento). Los
        materiales 2D con plane strain del catálogo sobreescriben este
        método con la fórmula vectorizada.
        """
        P = int(sigma.shape[0])
        fn = getattr(type(self), "out_of_plane_stress", None)
        base = getattr(Material, "out_of_plane_stress", None)
        if fn is None or fn is base:
            return np.zeros(P)
        from solidum.math.batch.schema import StateSchema
        schema = StateSchema(self.STATE_SCHEMA or {})
        return np.array([float(self.out_of_plane_stress(sigma[i], schema.unpack(S[i])))
                         for i in range(P)])


class Material(BatchMaterialHooks, ABC):
    """Clase base abstracta para todos los materiales de Solidum FEM.

    Contrato de subclases
    ---------------------
    STRAIN_DIM : ClassVar[int]
        Dimensión del vector de deformación que el material espera/produce:
            1  → escalar 1D (Truss, Frame)
            3  → Voigt 2D plano [εxx, εyy, γxy]
            6  → Voigt 3D       [εxx, εyy, εzz, γxy, γyz, γxz]
        Los elementos validan al construirse que `material.STRAIN_DIM`
        coincide con la dimensión que ellos requieren, atajando errores
        que de otro modo aparecerían como `matmul` críptico en runtime.

    PRIMARY_STATE_VAR : ClassVar[str | None]
        Nombre de la clave principal dentro del dict de state_vars
        (ej. 'alpha', 'damage'). Se usa en VtkExporter para exportar
        la variable de estado sin hardcodear nombres. `None` para
        materiales puramente elásticos sin variables internas.

    IS_SYMMETRIC : ClassVar[bool]
        Indica si la contribución del material a la matriz tangente es
        simétrica. ``True`` para elasticidad, daño isótropo, plasticidad
        asociada (J2, etc.). ``False`` para plasticidad no asociada,
        daño anisótropo con back-stress, viscoplasticidad con
        endurecimiento cinemático no asociado. La capa algebraica
        (ADR 0003) lo agrega para elegir entre Cholesky/LDLᵀ/LU.

    IS_UNILATERAL : ClassVar[bool]
        ``True`` para materiales con respuesta unilateral cuyo módulo
        tangente puede colapsar a cero en algún régimen (p. ej. cables
        sin rigidez en compresión, gaps de solo compresión). Solo se
        aceptan en elementos que declaren ``ACCEPTS_UNILATERAL = True``
        — habitualmente formulaciones corotacionales preparadas para
        manejar rigidez nula sin degenerar la matriz global.

    STATE_SCHEMA : ClassVar[dict]
        Variables internas como ``{nombre: forma}`` (ADR 0014). Obligatorio;
        ``{}`` si el material no tiene historia. Ver
        :class:`BatchMaterialHooks` para el resto del contrato por lotes
        (``BATCH_KERNEL``, ``batch_params``, ``batch_matrix``,
        ``initial_state``), todo opcional.
    """

    STRAIN_DIM: ClassVar[int]
    PRIMARY_STATE_VAR: ClassVar[str | None] = None
    IS_SYMMETRIC: ClassVar[bool] = True
    IS_UNILATERAL: ClassVar[bool] = False

    # Densidad física del medio continuo (kg/m³, en las unidades del usuario).
    # Opcional al construir: análisis que no invocan peso propio o matriz de
    # masa no necesitan declararla. Default `None` indica "no declarada".
    # Cuando un consumidor que requiere masa la lee y encuentra `None`, falla
    # con `ValueError` identificando el material — sin posibilidad de fallo
    # silencioso de masa cero. Ver ADR 0008.
    density: float | None = None

    @abstractmethod
    def compute_state(self, strain, state_vars=None):
        """Calcula esfuerzo, módulo tangente/secante y nuevas variables internas.

        Semántica trial/commit
        ----------------------
        ``compute_state`` se llama múltiples veces por paso de carga: una vez
        por cada iteración del solver no lineal (Newton-Raphson, arc-length,
        Newton dinámico). El ``new_state_vars`` devuelto es el estado **trial**
        del paso — refleja la actualización tentativa correspondiente al
        ``strain`` recibido, pero **no se considera definitivo hasta que el
        solver invoque ``commit_all_states()`` al converger el paso**.

        Quien gestiona la promoción trial → committed es el ``ElementState``
        (ver ``solidum.core.element``): cada iteración Newton produce un nuevo
        ``state_vars`` trial; el commit congela los trial actuales como
        committed para que en el siguiente paso se construya el predictor
        elástico a partir del estado convergido (``ε^p_n``, ``α_n``, ``κ_n``,
        ``d_n``, ...).

        Si el paso no converge, el solver llama a la lógica de rollback que
        descarta los trial — los committed sobreviven. Esto es esencial para
        bisección adaptativa y para line search: una llamada a
        ``compute_state`` con un ``strain`` distinto **no contamina** el
        estado committed, sólo el trial corriente.

        Las llamadas a ``compute_state`` deben ser puras respecto al estado
        que el material guarda internamente — toda la historia entra por
        ``state_vars`` y sale por ``new_state_vars``. La instancia ``self``
        sólo contiene parámetros constantes del modelo (E, ν, σ_y, ...).

        Returns
        -------
        (stress, tangent_modulus, new_state_vars)
            ``new_state_vars`` es **trial**; el solver decide cuándo
            promoverlo a committed.
        """
        pass

    # ------------------------------------------------------------------
    # Política de tolerancias del criterio de admisibilidad (ADR 0006)
    # ------------------------------------------------------------------

    def admissibility_scale(self, state_vars=None) -> float:
        """Escala característica del criterio de admisibilidad, **en unidades
        de esfuerzo**.

        Se usa para normalizar la tolerancia de la comparación ``f ≤ tol``
        (fluencia, umbral de daño, etc.) según el estado corriente del
        material, garantizando invariancia bajo cambio de unidades y
        adaptatividad ante endurecimiento o evolución de variables internas.

        Returns
        -------
        float
            Esfuerzo característico positivo. El default ``1.0`` solo lo
            usan materiales sin check de admisibilidad (elásticos puros,
            cables); para esos materiales la magnitud es irrelevante porque
            ``admissibility_scale`` no se invoca.
        """
        return 1.0

    def admissibility_tol(self, state_vars=None) -> float:
        """Tolerancia numérica del criterio de admisibilidad, en unidades de esfuerzo.

        Política única del proyecto (ADR 0006), combinada absoluta + relativa::

            tol = ADMISSIBILITY_TOL_ABS + ADMISSIBILITY_TOL_REL · escala

        donde ``escala = admissibility_scale(state_vars)``. Es el **único punto**
        donde vive la fórmula; modificarla aquí propaga a todos los consumidores,
        incluidos los que no pueden invocar ``is_admissible`` por estar dentro
        de un kernel ``@njit`` y precomputar la tolerancia fuera (ver VonMises2D).
        """
        return ADMISSIBILITY_TOL_ABS + ADMISSIBILITY_TOL_REL * self.admissibility_scale(state_vars)

    def is_admissible(self, f: float, state_vars=None) -> bool:
        """Decide si el estado de prueba está dentro de la región admisible.

        Aplica la política estándar combinada absoluta + relativa::

            f ≤ ADMISSIBILITY_TOL_ABS + ADMISSIBILITY_TOL_REL · escala

        donde ``escala = admissibility_scale(state_vars)``. Tanto ``f`` como
        la escala deben venir en unidades de esfuerzo. El término absoluto
        actúa como piso para estados donde ``escala → 0``.
        """
        return f <= self.admissibility_tol(state_vars)

    # ------------------------------------------------------------------
    # Componente fuera del plano (post-proceso de materiales 2D)
    # ------------------------------------------------------------------

    def out_of_plane_stress(self, sigma, state_vars=None) -> float:
        """``σ_zz`` asociada a un estado Voigt 2D ``σ = [σ_xx, σ_yy, σ_xy]``.

        Sólo tiene contenido en materiales 2D bajo *plane strain*, donde
        ``σ_zz ≠ 0`` pero no forma parte del vector que devuelve
        ``compute_state``. Los consumidores de post-proceso (invariantes 3D
        como Von Mises en el exportador VTK) la necesitan para no calcular
        un invariante de plane stress sobre un estado de plane strain.

        El default ``0.0`` cubre plane stress, los materiales 1D y los 3D
        (cuyo ``σ_zz`` ya viaja en el vector Voigt 6D). Las subclases 2D
        con hipótesis plane strain lo sobreescriben con su propia
        constitutiva; ver ``solidum.materials._plane_strain``.
        """
        return 0.0
