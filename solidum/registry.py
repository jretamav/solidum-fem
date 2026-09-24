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

Un registro que declara ``YAML_SECTION`` es una **familia de material**: el
lector YAML construye sus objetos a partir de esa sección sin conocerla por
nombre (ADR 0020, P1). Así un módulo de usuario declara su propia familia
sin tocar el programa principal.
"""
from __future__ import annotations

from typing import Any, Callable, ClassVar, Dict, Type


def _store(items: Dict[str, Type], name: str, klass: Type, kind: str) -> None:
    """Inscribe ``klass`` bajo ``name`` avisando si pisa a otra clase distinta.

    Con el descubrimiento recursivo el ganador de una colisión dependía del
    orden de importación, en silencio. Re-registrar el mismo objeto (p. ej.
    un módulo recargado) no avisa."""
    previous = items.get(name)
    if previous is not None and previous is not klass:
        from solidum.logging import get_logger
        get_logger("registry").warning(
            f"{kind} '{name}' ya estaba registrado ({previous.__module__}."
            f"{previous.__qualname__}); lo sustituye {klass.__module__}."
            f"{klass.__qualname__}."
        )
    items[name] = klass


class Registry:
    """Base genérica para registries con decorador.

    Cada subclase recibe su propio dict ``_items`` (``__init_subclass__``): si
    lo heredara, compartiría en silencio el almacenamiento de la base. El
    método ``register`` admite tres formas (ver módulo).

    Metadatos de familia (ADR 0020, P1), opcionales:

    ``YAML_SECTION``
        Sección de primer nivel del YAML con los objetos de la familia (lista
        de ``{id, type, …}``). Declararla convierte al registro en una familia
        de material que el lector YAML recorre sin conocerla por nombre.
    ``YAML_LABEL``
        Nombre de un objeto de la familia en los mensajes de error
        (``"material térmico"``).

    Y, para cualquier registro, ``SPEC_KIND``: el ``kind`` de las specs
    (``docs/specs/*.md``) cuyos componentes viven en él. ``tools/spec.py``
    localiza el registro por este valor, sin lista fija.
    """

    _items: Dict[str, Type] = {}
    _kind: str = "ítem"
    YAML_SECTION: ClassVar[str | None] = None
    YAML_LABEL: ClassVar[str] = "objeto"
    SPEC_KIND: ClassVar[str | None] = None
    _families: ClassVar[list] = []
    _spec_kinds: ClassVar[dict] = {}

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        if "_items" not in cls.__dict__:
            cls._items = {}
        kind = cls.__dict__.get("SPEC_KIND")
        if kind is not None:
            other = Registry._spec_kinds.get(kind)
            if other is not None and other.__qualname__ != cls.__qualname__:
                raise ValueError(
                    f"Tipo de spec '{kind}' ya declarado por "
                    f"{other.__module__}.{other.__qualname__}."
                )
            Registry._spec_kinds[kind] = cls
        section = cls.__dict__.get("YAML_SECTION")
        if section is not None:
            for other in Registry._families:
                if other.YAML_SECTION == section and other.__qualname__ != cls.__qualname__:
                    raise ValueError(
                        f"Sección YAML '{section}' ya declarada por "
                        f"{other.__module__}.{other.__qualname__}."
                    )
            # Recargar el módulo que define la familia la sustituye.
            Registry._families[:] = [f for f in Registry._families
                                     if f.__qualname__ != cls.__qualname__]
            Registry._families.append(cls)

    @staticmethod
    def families() -> list:
        """Familias de material (registros con ``YAML_SECTION``), en orden de
        definición."""
        return list(Registry._families)

    @staticmethod
    def for_spec_kind(kind: str):
        """Registro cuyo ``SPEC_KIND`` es ``kind``, o ``None``."""
        return Registry._spec_kinds.get(kind)

    @staticmethod
    def spec_kinds() -> list:
        """Todos los ``kind`` de spec declarados por algún registro."""
        return sorted(Registry._spec_kinds)

    @classmethod
    def register(cls, name_or_class: str | Type | None = None,
                 klass: Type | None = None) -> Type | Callable[[Type], Type]:
        # Forma legacy: register("Foo", FooClass)
        if klass is not None:
            _store(cls._items, name_or_class, klass, cls._kind)
            return klass

        # Forma decorador-sin-paréntesis: @register
        if isinstance(name_or_class, type):
            _store(cls._items, name_or_class.__name__, name_or_class, cls._kind)
            return name_or_class

        # Forma decorador-con-paréntesis: @register("Alias") o @register()
        def decorator(target: Type) -> Type:
            registered_name = name_or_class if isinstance(name_or_class, str) else target.__name__
            _store(cls._items, registered_name, target, cls._kind)
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


class MaterialRegistry(Registry):
    _items: Dict[str, Type] = {}
    _kind = "Material"
    SPEC_KIND = "material"
    YAML_SECTION = "materials"
    YAML_LABEL = "material"


class ThermalMaterialRegistry(Registry):
    """Registry paralelo a ``MaterialRegistry`` para materiales térmicos
    (Etapa 8). Separado intencionalmente: un material térmico relaciona el flujo de calor ``q`` con el gradiente
    ``∇T`` —dos vectores del espacio físico, sin notación Voigt ni tensor
    simétrico de por medio— en vez de ``σ`` con ``ε``. Compartir registro
    obligaría al parser YAML y a los elementos a discriminar por tipo en
    cada uso."""
    _items: Dict[str, Type] = {}
    _kind = "MaterialTermico"
    SPEC_KIND = "thermal_material"
    YAML_SECTION = "thermal_materials"
    YAML_LABEL = "material térmico"


class ElementRegistry(Registry):
    _items: Dict[str, Type] = {}
    _kind = "Elemento"
    SPEC_KIND = "element"


class SolverRegistry(Registry):
    _items: Dict[str, Type] = {}
    _kind = "Solucionador"
    SPEC_KIND = "solver"


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
