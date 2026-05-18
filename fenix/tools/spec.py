"""Parser y validador de specs de componentes (docs/specs/*.md).

Una spec es la única entrada del usuario al desarrollo: declara el contrato
físico de un elemento/material/solver. Este módulo extrae el bloque ```yaml
de la spec, valida su schema y lo confronta con la clase registrada.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List

import yaml

VALID_KINDS = {"element", "material", "cohesive_material", "solver"}
VALID_STATUSES = {"draft", "implemented", "validated"}

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


def validate_schema(spec: Spec) -> None:
    c = spec.contract
    for field in ("name", "kind", "status", "interface"):
        if field not in c:
            raise SpecError(f"{spec.path}: falta campo obligatorio '{field}'")

    if c["kind"] not in VALID_KINDS:
        raise SpecError(
            f"{spec.path}: kind='{c['kind']}' inválido (esperado {sorted(VALID_KINDS)})"
        )
    if c["status"] not in VALID_STATUSES:
        raise SpecError(
            f"{spec.path}: status='{c['status']}' inválido "
            f"(esperado {sorted(VALID_STATUSES)})"
        )

    if c["kind"] == "element":
        _validate_element_interface(spec)


def _validate_element_interface(spec: Spec) -> None:
    iface = spec.contract.get("interface") or {}
    required = {
        "dof_names": list,
        "n_nodes": int,
        "strain_dim": int,
        "n_integration_points": int,
    }
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
    for key in ("n_nodes", "strain_dim", "n_integration_points"):
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

    import fenix  # noqa: F401 — dispara autodiscover
    from fenix.registry import (
        CohesiveMaterialRegistry,
        ElementRegistry,
        MaterialRegistry,
        SolverRegistry,
    )

    registry = {
        "element": ElementRegistry,
        "material": MaterialRegistry,
        "cohesive_material": CohesiveMaterialRegistry,
        "solver": SolverRegistry,
    }[spec.kind]

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
        for attr, key in (
            ("DOF_NAMES", "dof_names"),
            ("N_NODES", "n_nodes"),
            ("STRAIN_DIM", "strain_dim"),
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
    elif spec.kind == "cohesive_material":
        cls = registry._items[spec.name]
        iface = spec.contract["interface"]
        for attr, key in (
            ("JUMP_DIM", "jump_dim"),
            ("PRIMARY_STATE_VAR", "primary_state_var"),
            ("IS_SYMMETRIC", "is_symmetric"),
        ):
            if key not in iface:
                continue  # campos opcionales en la spec; lo declarado debe coincidir
            if not hasattr(cls, attr):
                errors.append(f"{spec.path}: {cls.__name__} no declara {attr}")
                continue
            code_val = getattr(cls, attr)
            spec_val = iface[key]
            if code_val != spec_val:
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
