# solidum_fem/solidum/materials/orthotropic_2d.py
"""Elasticidad lineal ortótropa 2D bajo hipótesis plane stress.

Primer material anisótropo del catálogo. Ver docs/specs/Orthotropic2D.md y
ADR 0013 (dónde vive la orientación material).
"""
import numpy as np

from solidum.core.material import Material
from solidum.registry import MaterialRegistry


@MaterialRegistry.register
class Orthotropic2D(Material):
    """Elástico lineal ortótropo 2D, plane stress, con orientación de fibra.

    Cuatro constantes independientes ``(E1, E2, G12, nu12)`` más un ángulo
    ``theta`` que sitúa los ejes principales del material respecto a los ejes
    globales. La quinta constante ``nu21`` queda determinada por reciprocidad.

    Convención de subíndices de Poisson (Jones/Tsai)
    -----------------------------------------------
    ``nu12 = -eps_22 / eps_11`` bajo carga uniaxial en la dirección **1**:
    **primer índice = dirección de carga, segundo = dirección de la
    contracción medida**. Existe en la literatura la convención opuesta;
    confundirlas **invierte los papeles de E1 y E2 sin fallo ruidoso**. Fijada
    en la spec y en ADR 0013 §5.

    Ejes del material
    -----------------
    En bambú, la dirección **1** es la longitudinal (haces vasculares, la
    rígida) y la **2** la transversal en el plano de análisis. En ejes del
    material no hay acoplamiento entre esfuerzos normales y distorsión
    angular; al rotar ``theta`` aparece acoplamiento tracción-cortante
    (``C[0,2]``, ``C[1,2]`` no nulos): una tracción pura genera distorsión.

    La orientación vive **en el material** por decisión del ADR 0013: es un
    atajo consciente y desechable, válido mientras la fibra sea uniforme
    (probetas, tiras a ángulo constante). Cuando aparezca un caso con
    orientación variable en la malla, migra al elemento vía un parámetro
    ``orientation`` aditivo en ``compute_state``.

    Parameters
    ----------
    E1, E2 : float
        Módulos de Young en las direcciones 1 y 2. Estrictamente positivos.
    G12 : float
        Módulo de cortante en el plano 1-2. **Independiente** de E y nu — no
        se deriva de ellos como en isótropo; debe medirse.
    nu12 : float
        Coeficiente de Poisson mayor. Admisible ``|nu12| < sqrt(E1/E2)``, que
        **permite valores mayores que 0.5** en materiales muy ortótropos
        (ver `_validate`).
    theta : float, default 0.0
        Ángulo en **grados** de los ejes del material respecto a los globales,
        positivo antihorario (Reglas.md §5).
    density : float, optional
        Densidad del medio (ADR 0008). Sólo la exigen peso propio y masa.
    """

    STRAIN_DIM = 3
    PRIMARY_STATE_VAR = None
    IS_SYMMETRIC = True

    def __init__(self, E1: float, E2: float, G12: float, nu12: float,
                 theta: float = 0.0, density: float | None = None):
        self._validate(E1, E2, G12, nu12, density)

        self.E1, self.E2, self.G12 = float(E1), float(E2), float(G12)
        self.nu12 = float(nu12)
        self.theta = float(theta)
        self.density = density

        # Reciprocidad: la simetría de la matriz de flexibilidad (consecuencia
        # de la existencia de una energía de deformación) fija nu21. Por eso
        # hay cuatro constantes independientes y no cinco.
        self.nu21 = self.nu12 * self.E2 / self.E1

        self.C_mat = self._build_C_material()
        self.C = self._rotate(self.C_mat, self.theta)

    # ------------------------------------------------------------------
    # Validación
    # ------------------------------------------------------------------

    @staticmethod
    def _validate(E1, E2, G12, nu12, density) -> None:
        """Admisibilidad ortótropa — NO es la del isótropo.

        `Elastic2D` valida ``-1 < nu < 0.5``, criterio correcto para isótropo e
        **incorrecto aquí**. La condición ortótropa es que C sea definida
        positiva, que se reduce a ``1 - nu12*nu21 > 0``, es decir
        ``|nu12| < sqrt(E1/E2)``. Con E1/E2 = 20 el límite es 4.47: un nu12 de
        0.6 es físicamente legítimo en bambú y un validador heredado del
        isótropo lo rechazaría.
        """
        for nombre, valor in (("E1", E1), ("E2", E2), ("G12", G12)):
            if valor <= 0.0:
                raise ValueError(
                    f"Orthotropic2D: {nombre}={valor} debe ser estrictamente "
                    f"positivo (requisito de definición positiva de C)."
                )

        limite = np.sqrt(E1 / E2)
        if abs(nu12) >= limite:
            raise ValueError(
                f"Orthotropic2D: nu12={nu12} viola la admisibilidad ortótropa "
                f"|nu12| < sqrt(E1/E2) = {limite:.4f}. El criterio NO es el "
                f"isótropo (-1, 0.5): en un material muy ortótropo nu12 > 0.5 "
                f"puede ser legítimo. Revise E1, E2 y nu12, y compruebe que "
                f"nu12 sigue la convención Jones/Tsai (carga en 1, contracción "
                f"medida en 2)."
            )

        if density is not None and density < 0.0:
            raise ValueError(
                f"Orthotropic2D: density={density} no puede ser negativa."
            )

    # ------------------------------------------------------------------
    # Construcción de la constitutiva
    # ------------------------------------------------------------------

    def _build_C_material(self) -> np.ndarray:
        """Rigidez en ejes del material — sin acoplamiento (bloque cortante nulo).

        Inversa de la flexibilidad. El factor ``1 - nu12*nu21`` juega el papel
        estructural que ``1 - nu**2`` tiene en el isótropo.
        """
        denom = 1.0 - self.nu12 * self.nu21
        return np.array([
            [self.E1 / denom,             self.nu12 * self.E2 / denom, 0.0],
            [self.nu12 * self.E2 / denom, self.E2 / denom,             0.0],
            [0.0,                         0.0,                    self.G12],
        ])

    @staticmethod
    def _rotate(C_mat: np.ndarray, theta_deg: float) -> np.ndarray:
        """Lleva C de ejes del material a globales: ``C_glob = T^T C_mat T``.

        `T` es la transformación de **deformaciones** en Voigt engineering
        (gamma = 2*eps). Con engineering, las matrices de transformación de
        esfuerzo y de deformación **no coinciden** —difieren en factores 2 en
        el bloque cortante— y la forma ``T^T C T`` es válida precisamente para
        la T de deformaciones escrita aquí. Es el error clásico al programar
        esto; lo blinda el test de invariancia bajo rotación.
        """
        t = np.radians(theta_deg)
        c, s = np.cos(t), np.sin(t)
        T = np.array([
            [c * c,      s * s,      c * s],
            [s * s,      c * c,     -c * s],
            [-2.0 * c * s, 2.0 * c * s, c * c - s * s],
        ])
        return T.T @ C_mat @ T

    # ------------------------------------------------------------------
    # Contrato Material
    # ------------------------------------------------------------------

    def compute_state(self, strain: np.ndarray, state_vars=None):
        """Material sin historia: la tangente es constante en todo régimen.

        `C` se precalculó en el constructor (igual que hace `Elastic2D` con su
        hipótesis), así que aquí no se recomputa ni la rotación ni la inversa.
        """
        return self.C @ strain, self.C, state_vars
