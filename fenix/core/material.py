# fenix_fem/fenix/core/material.py
from abc import ABC, abstractmethod
from typing import ClassVar, Optional

from fenix.constants import ADMISSIBILITY_TOL_ABS, ADMISSIBILITY_TOL_REL


class Material(ABC):
    """Clase base abstracta para todos los materiales de Fenix FEM.

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

    PRIMARY_STATE_VAR : ClassVar[Optional[str]]
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
    """

    STRAIN_DIM: ClassVar[int]
    PRIMARY_STATE_VAR: ClassVar[Optional[str]] = None
    IS_SYMMETRIC: ClassVar[bool] = True
    IS_UNILATERAL: ClassVar[bool] = False

    # Densidad física del medio continuo (kg/m³, en las unidades del usuario).
    # OBLIGATORIA en cada subclase concreta: el constructor de todo material
    # debe declarar `density` como argumento requerido. Valor `0.0` se admite
    # para casos legítimos (materiales de penalty, restricción pura, fixtures
    # de test que ignoran masa) pero solo cuando el usuario lo declare
    # explícitamente — nunca por omisión silenciosa. Ver ADR 0008.
    density: float

    @abstractmethod
    def compute_state(self, strain, state_vars=None):
        """Calcula esfuerzo, módulo tangente/secante y nuevas variables internas.

        Returns
        -------
        (stress, tangent_modulus, new_state_vars)
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
