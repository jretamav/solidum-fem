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

Un registro que declara ``YAML_SECTION`` es una **familia**: el lector YAML
construye sus objetos a partir de esa sección sin conocerla por nombre (ADR
0020, P1). Así un módulo de usuario declara su propia familia sin tocar el
programa principal. Las familias del programa principal son las de
materiales (mecánicos y térmicos) y la de secciones de marco (ADR 0021).

Un objeto de familia se construye con ``cls(**datos)``; si su clase declara
el método de clase ``from_yaml(data, resolve)``, el lector lo usa en su
lugar y le pasa ``resolve(registro, id, contexto)`` para las referencias a
objetos de otras familias ya construidas (ADR 0021 §2).
"""
from __future__ import annotations

from typing import Any, Callable, ClassVar, Dict, Type


def _store(items: Dict[str, Type], name: str, klass: Type, label: str) -> None:
    """Inscribe ``klass`` bajo ``name`` avisando si pisa a otra clase distinta.

    Con el descubrimiento recursivo el ganador de una colisión dependía del
    orden de importación, en silencio. Re-registrar el mismo objeto (p. ej.
    un módulo recargado) no avisa."""
    previous = items.get(name)
    if previous is not None and previous is not klass:
        from solidum.logging import get_logger
        get_logger("registry").warning(
            f"{label} '{name}' ya estaba registrado ({previous.__module__}."
            f"{previous.__qualname__}); lo sustituye {klass.__module__}."
            f"{klass.__qualname__}."
        )
    items[name] = klass


class Registry:
    """Base genérica para registries con decorador.

    Cada subclase recibe su propio dict ``_items`` (``__init_subclass__``): si
    lo heredara, compartiría en silencio el almacenamiento de la base. Por eso
    una subclase no lo declara. El método ``register`` admite tres formas (ver
    módulo).

    Una subclase declara sólo atributos públicos (ADR 0020, P1); son todo lo
    que un módulo de usuario necesita para definir su propio registro:

    ``YAML_LABEL``
        Nombre de un objeto del registro en los mensajes de error, tanto del
        lector YAML como del propio registro (``"material térmico"``). En
        minúscula: los mensajes lo insertan a media frase.
    ``SPEC_KIND``
        El ``kind`` de las specs (``docs/specs/*.md``) cuyos componentes viven
        en él. ``tools/spec.py`` localiza el registro por este valor, sin
        lista fija.
    ``YAML_SECTION``
        Sólo en una familia: sección de primer nivel del YAML con sus objetos
        (lista de ``{id, type, …}``). Declararla convierte al registro en una
        familia que el lector YAML recorre sin conocerla por nombre.
    """

    _items: Dict[str, Type] = {}
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
        """Familias (registros con ``YAML_SECTION``), en orden de definición:
        el orden en que el lector YAML las construye."""
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
    def _label(cls) -> str:
        """``YAML_LABEL`` con mayúscula inicial, para abrir un mensaje."""
        return cls.YAML_LABEL[:1].upper() + cls.YAML_LABEL[1:]

    @classmethod
    def register(cls, name_or_class: str | Type | None = None,
                 klass: Type | None = None) -> Type | Callable[[Type], Type]:
        # Forma legacy: register("Foo", FooClass)
        if klass is not None:
            _store(cls._items, name_or_class, klass, cls._label())
            return klass

        # Forma decorador-sin-paréntesis: @register
        if isinstance(name_or_class, type):
            _store(cls._items, name_or_class.__name__, name_or_class, cls._label())
            return name_or_class

        # Forma decorador-con-paréntesis: @register("Alias") o @register()
        def decorator(target: Type) -> Type:
            registered_name = name_or_class if isinstance(name_or_class, str) else target.__name__
            _store(cls._items, registered_name, target, cls._label())
            return target
        return decorator

    @classmethod
    def create(cls, name: str, **kwargs) -> Any:
        if name not in cls._items:
            raise ValueError(
                f"{cls._label()} '{name}' no registrado. "
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
                f"{cls._label()} '{name}' no registrado. "
                f"Disponibles: {sorted(cls._items.keys())}"
            )
        return cls._items[name]


class MaterialRegistry(Registry):
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
    SPEC_KIND = "thermal_material"
    YAML_SECTION = "thermal_materials"
    YAML_LABEL = "material térmico"


class SectionRegistry(Registry):
    """Familia de secciones de marco (ADR 0021).

    Una sección se declara una vez en la sección YAML ``sections`` y la
    comparten muchos marcos, como un material. Sus objetos referencian
    materiales dentro de sus parches: por eso se define **después** de
    ``MaterialRegistry`` (las familias se construyen en orden de definición)
    y su clase construye el objeto con el gancho genérico
    ``from_yaml(data, resolve)`` del lector YAML."""
    SPEC_KIND = "section"
    YAML_SECTION = "sections"
    YAML_LABEL = "sección"


class ElementRegistry(Registry):
    SPEC_KIND = "element"
    YAML_LABEL = "elemento"


class SolverRegistry(Registry):
    SPEC_KIND = "solver"
    YAML_LABEL = "solver"


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
