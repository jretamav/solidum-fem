"""Registries de Solidum FEM (Material, Element, Solver, Quadrature).

Soporta tres formas de registro:

    @MaterialRegistry.register            # nombre = cls.__name__
    class Foo(Material): ...

    @MaterialRegistry.register("Alias")   # nombre explícito
    class Foo(Material): ...

    MaterialRegistry.register("Foo", Foo) # forma legacy explícita

Las clases se descubren automáticamente al importar `solidum` gracias a
`solidum.autodiscover.initialize()`, eliminando la necesidad de mantener
listas de imports en `registry_initialization.py`.
"""
from __future__ import annotations

from typing import Any, Callable, Dict, Type


class _BaseRegistry:
    """Base genérica para registries con decorador.

    Cada subclase debe declarar su propio dict `_items` para no compartir
    almacenamiento.  El método `register` admite tres formas (ver módulo).
    """

    _items: Dict[str, Type] = {}
    _kind: str = "ítem"

    @classmethod
    def register(cls, name_or_class: str | Type | None = None,
                 klass: Type | None = None) -> Type | Callable[[Type], Type]:
        # Forma legacy: register("Foo", FooClass)
        if klass is not None:
            cls._items[name_or_class] = klass
            return klass

        # Forma decorador-sin-paréntesis: @register
        if isinstance(name_or_class, type):
            cls._items[name_or_class.__name__] = name_or_class
            return name_or_class

        # Forma decorador-con-paréntesis: @register("Alias") o @register()
        def decorator(target: Type) -> Type:
            registered_name = name_or_class if isinstance(name_or_class, str) else target.__name__
            cls._items[registered_name] = target
            return target
        return decorator

    @classmethod
    def create(cls, name: str, **kwargs) -> Any:
        if name not in cls._items:
            raise ValueError(
                f"{cls._kind} '{name}' no registrado. "
                f"Disponibles: {sorted(cls._items.keys())}"
            )
        return cls._items[name](**kwargs)

    @classmethod
    def names(cls) -> list:
        return sorted(cls._items.keys())

    @classmethod
    def get(cls, name: str) -> Type:
        """Devuelve la clase registrada bajo ``name`` sin instanciarla.

        Útil para introspección (firmas de constructor, atributos de clase
        declarativos como ``PIPELINE_KIND``) cuando ``create(name, **kwargs)``
        no aplica porque se quiere consultar la clase sin construirla.

        Raises
        ------
        KeyError
            Si ``name`` no está registrado.
        """
        if name not in cls._items:
            raise KeyError(
                f"{cls._kind} '{name}' no registrado. "
                f"Disponibles: {sorted(cls._items.keys())}"
            )
        return cls._items[name]


class MaterialRegistry(_BaseRegistry):
    _items: Dict[str, Type] = {}
    _kind = "Material"


class CohesiveMaterialRegistry(_BaseRegistry):
    """Registry paralelo a ``MaterialRegistry`` para materiales cohesivos
    traction-jump (ADR 0010). Separado intencionalmente: los cohesivos
    operan sobre ``[[u]]`` y devuelven ``t`` sobre ``Γ_d``, no sobre Voigt
    de ``ε`` en el bulk. Mezclarlos forzaría al parser YAML a discriminar
    por tipo en cada uso."""
    _items: Dict[str, Type] = {}
    _kind = "MaterialCohesivo"


class ElementRegistry(_BaseRegistry):
    _items: Dict[str, Type] = {}
    _kind = "Elemento"


class SolverRegistry(_BaseRegistry):
    _items: Dict[str, Type] = {}
    _kind = "Solucionador"


class QuadratureRegistry:
    """Registry para reglas de integración (puntos, pesos)."""
    _rules: Dict[str, tuple] = {}

    @classmethod
    def register(cls, name: str, points: list, weights: list) -> None:
        cls._rules[name] = (points, weights)

    @classmethod
    def get(cls, name: str) -> tuple:
        if name not in cls._rules:
            raise ValueError(
                f"Cuadratura '{name}' no registrada. "
                f"Disponibles: {sorted(cls._rules.keys())}"
            )
        return cls._rules[name]
