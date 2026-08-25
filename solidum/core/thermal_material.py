# solidum_fem/solidum/core/thermal_material.py
"""Clase base abstracta para materiales térmicos (Etapa 8).

Familia paralela a ``Material``: relaciona el vector de flujo de calor
``q`` con el gradiente de temperatura ``∇T`` —dos vectores del espacio
físico— en vez de ``σ`` con ``ε`` en notación Voigt. Ni el gradiente ni el
flujo son tensores simétricos, así que la maquinaria de Voigt no aplica.

Es la misma razón que motivó separar ``CohesiveMaterial`` (ADR 0010):
compartir registro y clase base obligaría al parser YAML y a los elementos
a discriminar por tipo en cada uso.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import ClassVar


class ThermalMaterial(ABC):
    """Clase base de los materiales térmicos de Solidum FEM (Etapa 8).

    Contrato de subclases
    ---------------------
    FLUX_DIM : ClassVar[int]
        Dimensión del gradiente de temperatura y del flujo:
            2 → problema plano, ``∇T = (∂T/∂x, ∂T/∂y)``
            3 → problema tridimensional
        Los elementos térmicos validan al construirse que coincide con la
        dimensión que su formulación requiere. A diferencia de
        ``STRAIN_DIM``, no es una dimensión de Voigt sino la dimensión
        física del espacio.

    PRIMARY_STATE_VAR : ClassVar[str | None]
        Nombre de la variable interna exportada al post-proceso. ``None``
        en modelos lineales sin memoria, que es el caso de la conducción
        de Fourier.

    IS_SYMMETRIC : ClassVar[bool]
        ``True`` si la contribución del material a la matriz del elemento
        es simétrica. Lo es siempre que el tensor de conductividad lo sea
        —requisito de las relaciones recíprocas de Onsager—, de modo que
        el despachador algebraico (ADR 0003) puede elegir Cholesky.

    Notas sobre unidades
    --------------------
    Los términos térmicos tienen unidades propias: conductividad ``k``
    [W/(m·K)], calor específico ``c`` [J/(kg·K)], flujo ``q`` [W/m²],
    fuente volumétrica ``Q`` [W/m³]. La ``density`` comparte significado y
    unidades con la de los materiales mecánicos [kg/m³].

    Igual que en el resto del catálogo, el sistema **no convierte
    unidades**: la consistencia es responsabilidad del usuario
    (``Reglas.md §5``).
    """

    FLUX_DIM: ClassVar[int]
    PRIMARY_STATE_VAR: ClassVar[str | None] = None
    IS_SYMMETRIC: ClassVar[bool] = True

    # Propiedades de la capacidad calorífica volumétrica ρ·c. Ambas son
    # opcionales al construir, siguiendo el criterio del ADR 0008 para
    # ``density``: un análisis estacionario sólo usa la conductividad, y
    # obligar a declarar propiedades que no entran en ninguna ecuación
    # forzaría a inventar valores. Cuando un consumidor transitorio las
    # encuentra ``None``, falla con ``ValueError`` accionable — nunca con
    # capacidad cero silenciosa.
    density: float | None = None
    specific_heat: float | None = None

    @abstractmethod
    def compute_flux(self, grad_T):
        """Calcula el flujo de calor y el tensor de conductividad.

        Parameters
        ----------
        grad_T : np.ndarray
            Gradiente de temperatura ``∇T ∈ ℝ^{FLUX_DIM}`` en ejes
            globales, tal como lo entrega la matriz ``B`` del elemento.

        Returns
        -------
        (flux, conductivity)
            ``flux``: ``np.ndarray`` de tamaño ``FLUX_DIM`` — el vector
            ``q = -k·∇T``. El signo negativo es la segunda ley de la
            termodinámica, no una convención elegible.
            ``conductivity``: ``np.ndarray`` ``(FLUX_DIM, FLUX_DIM)`` — el
            tensor ``k``, análogo al módulo tangente del contrato mecánico.

        Notes
        -----
        Sin parámetro ``state_vars``: la conducción de Fourier es lineal y
        sin memoria. Cuando entre un material térmico con historia
        (conductividad dependiente de la temperatura, cambio de fase), la
        firma se ampliará entonces, con el caso real delante.
        """
        pass

    def volumetric_capacity(self, consumer: str = "el análisis transitorio") -> float:
        """Devuelve ``ρ·c`` [J/(m³·K)], validando que ambas estén declaradas.

        Parameters
        ----------
        consumer : str
            Descripción de quién solicita la capacidad, para que el mensaje
            de error diga en qué contexto falta el dato.

        Raises
        ------
        ValueError
            Si falta ``density``, ``specific_heat`` o ambas. El mensaje
            nombra el material, los parámetros ausentes y **la salida** —
            no sólo lo que falta, sino qué hacer al respecto.
        """
        faltan = []
        if self.density is None:
            faltan.append("'density'")
        if self.specific_heat is None:
            faltan.append("'c'")

        if faltan:
            lista = " y ".join(faltan)
            if len(faltan) > 1:
                sujeto, verbo, imperativo = "los parámetros", "declararon", "Decláralos"
            else:
                sujeto, verbo, imperativo = "el parámetro", "declaró", "Decláralo"
            raise ValueError(
                f"{type(self).__name__}: {consumer} requiere {sujeto} {lista}, "
                f"que no se {verbo} al construir el material. "
                f"{imperativo} en el material, o usa un solver estacionario "
                f"(LinearSolver), que no necesita la capacidad calorífica."
            )

        return self.density * self.specific_heat
