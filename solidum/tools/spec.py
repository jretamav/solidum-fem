"""Parser y validador de specs de componentes (docs/specs/*.md).

Una spec es la única entrada del usuario al desarrollo: declara el contrato
físico de un elemento/material/solver. Este módulo extrae el bloque ```yaml
de la spec, valida su schema y lo confronta con la clase registrada.

Las specs del programa principal (``docs/specs/``) declaran además el bloque
``admission``: el criterio de admisión de Reglas.md §1 aplicado al componente,
con sus referencias en las listas cerradas de textos y de códigos generales.
Las de un módulo de usuario (``docs/user/<módulo>/specs/``, ADR 0020) no están
obligadas a declararlo.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import yaml

VALID_STATUSES = {"draft", "implemented", "validated"}

# Criterio de admisión al programa principal (Reglas.md §1): una formulación
# es estándar si aparece en un texto de referencia del área y está
# implementada en un código general establecido. Listas cerradas:
# identificador canónico (el que se escribe en la spec) → nombre con que
# Reglas.md §1 lo enumera. tests/test_specs.py comprueba que los nombres
# coinciden con Reglas.md; ampliar una lista es decisión del usuario y se hace
# primero allí.
ADMISSION_TEXTBOOKS: Dict[str, str] = {
    "zienkiewicz_taylor": "Zienkiewicz-Taylor",
    "hughes": "Hughes",
    "bathe": "Bathe",
    "belytschko_liu_moran": "Belytschko-Liu-Moran",
    "de_souza_neto": "de Souza Neto-Perić-Owen",
    "de_borst_crisfield": "de Borst-Crisfield-Remmers-Verhoosel",
    # Dinámica estructural (Chopra, Dynamics of Structures): añadido por el
    # usuario el 2026-09-30 para cubrir la combinación modal y los
    # integradores temporales, que los textos de MEF sólo tratan en parte.
    "chopra": "Chopra",
}
ADMISSION_GENERAL_CODES: Dict[str, str] = {
    "abaqus": "Abaqus",
    "ansys": "ANSYS",
    "diana": "DIANA",
    "code_aster": "Code_Aster",
    "opensees": "OpenSees",
}
# Clave del bloque ``admission`` → lista cerrada que la gobierna.
ADMISSION_LISTS: Dict[str, Dict[str, str]] = {
    "textbooks": ADMISSION_TEXTBOOKS,
    "general_codes": ADMISSION_GENERAL_CODES,
}


def valid_kinds() -> set:
    """``kind`` admitidos: el ``SPEC_KIND`` de cada registro (ADR 0020). Un
    módulo de usuario que declare su propio registro añade el suyo; para
    validar sus specs hay que cargarlo antes."""
    import solidum  # noqa: F401 — dispara autodiscover
    from solidum.registry import Registry
    return set(Registry.spec_kinds())

# Campo físico que gobierna un elemento. Determina qué declara su interfaz:
# un elemento mecánico tiene `strain_dim` (dimensión Voigt de la deformación
# que entrega al material); uno térmico no tiene deformación alguna, sino un
# gradiente escalar. Declararlo explícitamente evita que la validación
# presuponga que todo elemento es mecánico.
VALID_FIELDS = {"displacement", "temperature"}
DEFAULT_FIELD = "displacement"

_YAML_BLOCK_RE = re.compile(r"```yaml\s*\n(.*?)\n```", re.DOTALL)
_IMPLEMENTATION_SECTION_RE = re.compile(
    r"^##\s+Implementaci[oó]n\s*$(.*?)(?=^##\s|\Z)",
    re.MULTILINE | re.DOTALL | re.IGNORECASE,
)


@dataclass
class Spec:
    path: Path
    contract: Dict[str, Any]
    implementation_filled: bool

    @property
    def name(self) -> str:
        return self.contract["name"]

    @property
    def kind(self) -> str:
        return self.contract["kind"]

    @property
    def status(self) -> str:
        return self.contract["status"]


class SpecError(ValueError):
    """Error de schema o coherencia en una spec."""


def parse_spec(path: Path) -> Spec:
    text = Path(path).read_text(encoding="utf-8")
    match = _YAML_BLOCK_RE.search(text)
    if not match:
        raise SpecError(f"{path}: no se encontró bloque ```yaml con el contrato")
    contract = yaml.safe_load(match.group(1))
    if not isinstance(contract, dict):
        raise SpecError(f"{path}: el bloque YAML no es un mapping")

    impl_match = _IMPLEMENTATION_SECTION_RE.search(text)
    impl_body = impl_match.group(1).strip() if impl_match else ""
    # Consideramos "rellena" si hay al menos un archivo o clase declarados
    # (la plantilla deja guiones "—" como placeholder).
    # Acepta "Archivo: path", "**Archivo**: path", "- Archivo: path", etc.
    # Lo que indica 'rellena' es que tras "Archivo"/"Clase" y dos puntos haya
    # algo distinto del placeholder "—" o "-" solo.
    filled = bool(
        re.search(
            r"(Archivo|Clase)\*{0,2}\s*:\s*(?![—\-]\s*$)\S+",
            impl_body,
            flags=re.MULTILINE,
        )
    )

    return Spec(path=Path(path), contract=contract, implementation_filled=filled)


def is_user_module_spec(path: Path) -> bool:
    """¿Es ``path`` la spec de un módulo de usuario? Viven bajo
    ``docs/user/<módulo>/specs/`` (ADR 0020); el resto son del programa
    principal."""
    parts = Path(path).resolve().parts
    return any(a == "docs" and b == "user" for a, b in zip(parts, parts[1:]))


def validate_schema(spec: Spec, *, main_program: Optional[bool] = None) -> None:
    """Valida el schema del contrato.

    ``main_program`` indica si la spec es del programa principal, que exige
    el bloque ``admission``; por omisión se deduce de la ruta
    (:func:`is_user_module_spec`). En una spec de módulo de usuario el bloque
    es opcional —una pieza estándar puede vivir en un módulo mientras sólo él
    la consuma (Reglas.md §4)—, pero si se declara se valida igual.
    """
    c = spec.contract
    for field in ("name", "kind", "status", "interface"):
        if field not in c:
            raise SpecError(f"{spec.path}: falta campo obligatorio '{field}'")

    kinds = valid_kinds()
    if c["kind"] not in kinds:
        raise SpecError(
            f"{spec.path}: kind='{c['kind']}' inválido (esperado {sorted(kinds)})"
        )
    if c["status"] not in VALID_STATUSES:
        raise SpecError(
            f"{spec.path}: status='{c['status']}' inválido "
            f"(esperado {sorted(VALID_STATUSES)})"
        )

    if main_program is None:
        main_program = not is_user_module_spec(spec.path)
    _validate_admission(spec, required=main_program)

    if c["kind"] == "element":
        _validate_element_interface(spec)


def admission_entries(spec: Spec) -> Dict[str, List[Tuple[str, Optional[str]]]]:
    """Bloque ``admission`` normalizado: ``{lista: [(identificador, nota)]}``.

    Cada entrada de la spec es el identificador solo (``- hughes``) o un
    mapping de una clave con su nota breve (``- hughes: "§9.3"``); la nota
    sale ``None`` si falta. Lanza :class:`SpecError` si el formato o algún
    identificador no son válidos; ``{}`` si la spec no declara el bloque.
    """
    block = spec.contract.get("admission")
    if block is None:
        return {}
    where = f"{spec.path}: admission"
    if not isinstance(block, dict):
        raise SpecError(f"{where} debe ser un mapping con {sorted(ADMISSION_LISTS)}")
    unknown = sorted(set(block) - set(ADMISSION_LISTS))
    if unknown:
        raise SpecError(
            f"{where}: claves desconocidas {unknown} (esperadas {sorted(ADMISSION_LISTS)})"
        )

    result: Dict[str, List[Tuple[str, Optional[str]]]] = {}
    for key, closed in ADMISSION_LISTS.items():
        entries = block.get(key)
        valid = ", ".join(f"{k} ({v})" for k, v in closed.items())
        if not isinstance(entries, list) or not entries:
            raise SpecError(
                f"{where}.{key} debe ser una lista con al menos una entrada de la "
                f"lista cerrada de Reglas.md §1: {valid}"
            )
        normalized: List[Tuple[str, Optional[str]]] = []
        for entry in entries:
            if isinstance(entry, str):
                ref, note = entry, None
            elif isinstance(entry, dict) and len(entry) == 1:
                (ref, note), = entry.items()
                if note is not None and not isinstance(note, str):
                    raise SpecError(
                        f"{where}.{key}: la nota de '{ref}' debe ser texto, "
                        f"no {type(note).__name__}"
                    )
            else:
                raise SpecError(
                    f"{where}.{key}: entrada {entry!r} inválida; se espera "
                    f"'identificador' o 'identificador: \"nota\"'"
                )
            if ref not in closed:
                raise SpecError(
                    f"{where}.{key}: '{ref}' no está en la lista cerrada de "
                    f"Reglas.md §1. Válidos: {valid}"
                )
            if any(ref == seen for seen, _ in normalized):
                raise SpecError(f"{where}.{key}: '{ref}' repetido")
            normalized.append((ref, note or None))
        result[key] = normalized
    return result


def _validate_admission(spec: Spec, *, required: bool) -> None:
    """Criterio de admisión de Reglas.md §1: al menos un texto de referencia
    y un código general, de las listas cerradas. Obligatorio en el programa
    principal: si no se pueden rellenar las dos listas, el componente no es
    estándar y va a un módulo de usuario (ADR 0020)."""
    if required and spec.contract.get("admission") is None:
        raise SpecError(
            f"{spec.path}: falta el bloque 'admission' (criterio de admisión al "
            "programa principal, Reglas.md §1): declara 'textbooks' y "
            "'general_codes' con referencias de las listas cerradas. Si no se "
            "pueden rellenar las dos, el componente va a un módulo de usuario."
        )
    admission_entries(spec)


def _validate_element_interface(spec: Spec) -> None:
    """Valida la interfaz de un elemento según el campo físico que resuelve.

    Lo común a todo elemento —nodos, DOFs, puntos de integración— se exige
    siempre. Lo que depende de la física se pide por campo: un elemento de
    desplazamientos declara `strain_dim` (dimensión Voigt de la deformación
    que entrega al material), mientras que uno térmico no tiene deformación
    y declara `flux_dim` (dimensión del gradiente de temperatura, 2 ó 3).
    """
    iface = spec.contract.get("interface") or {}

    field = iface.get("field", DEFAULT_FIELD)
    if field not in VALID_FIELDS:
        raise SpecError(
            f"{spec.path}: interface.field='{field}' inválido "
            f"(esperado {sorted(VALID_FIELDS)})"
        )

    # Campos comunes a cualquier elemento, sea cual sea su física.
    required: dict[str, type] = {
        "dof_names": list,
        "n_nodes": int,
        "n_integration_points": int,
    }
    # Campo dimensional propio de cada física.
    dim_key = {"displacement": "strain_dim", "temperature": "flux_dim"}[field]
    required[dim_key] = int

    for key, typ in required.items():
        if key not in iface:
            raise SpecError(f"{spec.path}: interface.{key} ausente")
        if not isinstance(iface[key], typ):
            raise SpecError(
                f"{spec.path}: interface.{key} debe ser {typ.__name__}, "
                f"got {type(iface[key]).__name__}"
            )
    if not all(isinstance(d, str) for d in iface["dof_names"]):
        raise SpecError(f"{spec.path}: interface.dof_names debe ser lista de str")
    for key in ("n_nodes", "n_integration_points", dim_key):
        if iface[key] < 1:
            raise SpecError(f"{spec.path}: interface.{key} debe ser ≥ 1")


def cross_check_with_registry(spec: Spec) -> List[str]:
    """Verifica coherencia spec ↔ clase registrada.

    Devuelve lista de inconsistencias (vacía si todo ok). Solo aplica cuando
    status ∈ {implemented, validated}: en draft la clase puede aún no existir.
    """
    if spec.status == "draft":
        return []

    errors_pre: List[str] = []
    if not spec.implementation_filled:
        errors_pre.append(
            f"{spec.path}: status={spec.status} requiere sección "
            "'## Implementación' rellena (Archivo: ... / Clase: ...)"
        )

    import solidum  # noqa: F401 — dispara autodiscover
    from solidum.registry import Registry

    registry = Registry.for_spec_kind(spec.kind)

    errors: List[str] = list(errors_pre)
    if spec.name not in registry._items:
        errors.append(
            f"{spec.path}: status={spec.status} pero '{spec.name}' no está "
            f"registrado en {registry.__name__}"
        )
        return errors

    if spec.kind == "element":
        cls = registry._items[spec.name]
        iface = spec.contract["interface"]
        # El atributo dimensional depende del campo físico: los elementos de
        # desplazamientos declaran STRAIN_DIM; los térmicos, FLUX_DIM.
        field = iface.get("field", DEFAULT_FIELD)
        dim_attr, dim_key = {
            "displacement": ("STRAIN_DIM", "strain_dim"),
            "temperature": ("FLUX_DIM", "flux_dim"),
        }[field]
        for attr, key in (
            ("DOF_NAMES", "dof_names"),
            ("N_NODES", "n_nodes"),
            (dim_attr, dim_key),
            ("N_INTEGRATION_POINTS", "n_integration_points"),
        ):
            if not hasattr(cls, attr):
                # N_NODES no siempre se declara; skip si ausente
                if attr == "N_NODES":
                    continue
                errors.append(f"{spec.path}: {cls.__name__} no declara {attr}")
                continue
            code_val = getattr(cls, attr)
            spec_val = iface[key]
            if list(code_val) != list(spec_val) if attr == "DOF_NAMES" else code_val != spec_val:
                errors.append(
                    f"{spec.path}: {attr}={code_val!r} (código) vs "
                    f"{key}={spec_val!r} (spec)"
                )
    return errors


def collect_specs(specs_dir: Path) -> List[Path]:
    """Lista specs reales (excluye plantillas `_template_*.md`)."""
    return sorted(
        p for p in Path(specs_dir).glob("*.md") if not p.name.startswith("_template")
    )
