"""YAML loading and validation."""

from __future__ import annotations

from pathlib import Path
import math
import re
from typing import Any

import yaml

from .exceptions import ConfigError
from .local_nist import parse_composition, read_elements, validate_density
from .models import LayerSpec, MaterialSpec
from .physics import validate_positive_finite


SAFE_ID = re.compile(r"^[A-Za-z0-9_-]+$")
DISALLOWED_INPUT_KEYS = {
    "coefficient_file",
    "attenuation_file",
    "absorption_file",
    "spectrum_file",
    "incident_spectrum",
    "input_spectrum",
    "online_result_file",
    "results_file",
}
WORKFLOWS = {
    "coefficients",
    "attenuation_vs_thickness",
    "energy_absorption_vs_thickness",
    "multilayer_attenuation_profile",
}


def load_yaml(path: str | Path) -> dict[str, Any]:
    config_path = Path(path)
    try:
        data = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise ConfigError(f"Configuration file is missing: {config_path}") from exc
    except yaml.YAMLError as exc:
        raise ConfigError(f"Invalid YAML in {config_path}: {exc}") from exc
    if not isinstance(data, dict):
        raise ConfigError("YAML configuration must be a mapping")
    return data


def resolve_output_dir(value: object) -> Path:
    if value is None or str(value).strip() == "":
        raise ConfigError("output_dir is required")
    path = Path(str(value))
    return path if path.is_absolute() else Path.cwd() / path


def reject_disallowed_inputs(node: object, path: str = "$") -> None:
    if isinstance(node, dict):
        for key, value in node.items():
            if key in DISALLOWED_INPUT_KEYS:
                raise ConfigError(f"{path}.{key} is not allowed in YAML workflows")
            if isinstance(value, str) and re.search(r"(^|[\\/])results[\\/].*\.txt$", value):
                raise ConfigError(f"{path}.{key} must not point to historical results TXT files")
            reject_disallowed_inputs(value, f"{path}.{key}")
    elif isinstance(node, list):
        for index, value in enumerate(node):
            reject_disallowed_inputs(value, f"{path}[{index}]")


def validate_common(data: dict[str, Any]) -> tuple[str, Path]:
    reject_disallowed_inputs(data)
    if data.get("schema_version") != 1:
        raise ConfigError("schema_version must be 1")
    workflow = data.get("workflow")
    if workflow not in WORKFLOWS:
        raise ConfigError(f"workflow must be one of {sorted(WORKFLOWS)}")
    output_dir = resolve_output_dir(data.get("output_dir"))
    return str(workflow), output_dir


def _safe_id(value: object, field: str = "id") -> str:
    text = str(value or "").strip()
    if not SAFE_ID.fullmatch(text):
        raise ConfigError(f"{field} may contain only letters, digits, underscore, and hyphen")
    return text


def _finite_positive_mapping(composition: dict[str, Any]) -> dict[str, float]:
    try:
        parsed = parse_composition(composition)
    except Exception as exc:
        raise ConfigError(str(exc)) from exc
    return parsed


def material_from_dict(data: dict[str, Any], *, require_density: bool, allow_null_density: bool) -> MaterialSpec:
    if not isinstance(data, dict):
        raise ConfigError("material entries must be mappings")
    material_id = _safe_id(data.get("id"))
    basis = data.get("composition_basis")
    if basis not in {"atomic_ratio", "mass_fraction"}:
        raise ConfigError("composition_basis must be atomic_ratio or mass_fraction")
    composition = data.get("composition")
    if not isinstance(composition, dict):
        raise ConfigError("composition must be a mapping")
    parsed_composition = _finite_positive_mapping(composition)
    density_value = data.get("density_g_cm3")
    if density_value is None and not allow_null_density:
        raise ConfigError("density_g_cm3 may be null only in the coefficients workflow")
    try:
        density = validate_density(density_value, required=require_density)
    except Exception as exc:
        raise ConfigError(str(exc)) from exc
    density_source = str(data.get("density_source") or "not specified by user")
    # Force element-range validation before any output directory is created.
    elements = read_elements()
    unknown = [symbol for symbol in parsed_composition if symbol not in elements]
    if unknown:
        raise ConfigError(f"Unknown elements: {unknown}")
    unsupported = [symbol for symbol in parsed_composition if elements[symbol].atomic_number > 92]
    if unsupported:
        raise ConfigError(f"NIST v1.4 coefficient tables are limited to Z=1-92: {unsupported}")
    return MaterialSpec(
        id=material_id,
        composition_basis=str(basis),
        composition=parsed_composition,
        density_g_cm3=density,
        density_source=density_source,
    )


def layer_from_dict(data: dict[str, Any]) -> LayerSpec:
    material = material_from_dict(data, require_density=True, allow_null_density=False)
    try:
        thickness = validate_positive_finite(data.get("thickness_cm"), "thickness_cm")
    except Exception as exc:
        raise ConfigError(str(exc)) from exc
    return LayerSpec(
        id=material.id,
        composition_basis=material.composition_basis,
        composition=material.composition,
        density_g_cm3=material.density_g_cm3,
        density_source=material.density_source,
        thickness_cm=thickness,
    )


def energies_from_config(value: object, *, allow_null: bool) -> list[float] | None:
    if value is None:
        if allow_null:
            return None
        raise ConfigError("energies_kev is required")
    if not isinstance(value, list) or not value:
        raise ConfigError("energies_kev must be a non-empty list")
    energies: list[float] = []
    for item in value:
        try:
            number = float(item)
        except (TypeError, ValueError) as exc:
            raise ConfigError("energies_kev must contain numbers") from exc
        if not math.isfinite(number) or number <= 0:
            raise ConfigError("energies_kev must contain only positive finite values")
        energies.append(number)
    if len(set(energies)) != len(energies):
        raise ConfigError("energies_kev must not contain duplicates")
    return energies


def energy_range_from_config(value: object) -> tuple[float, float] | None:
    if value is None:
        return None
    if not isinstance(value, dict):
        raise ConfigError("energy_range_kev must be a mapping with start and stop")
    allowed = {"start", "stop"}
    unknown = set(value) - allowed
    if unknown:
        raise ConfigError(f"energy_range_kev contains unsupported keys: {sorted(unknown)}")
    if "start" not in value or "stop" not in value:
        raise ConfigError("energy_range_kev requires start and stop")
    try:
        start = validate_positive_finite(value.get("start"), "energy_range_kev.start")
        stop = validate_positive_finite(value.get("stop"), "energy_range_kev.stop")
    except Exception as exc:
        raise ConfigError(str(exc)) from exc
    if start >= stop:
        raise ConfigError("energy_range_kev.start must be smaller than energy_range_kev.stop")
    return start, stop


def thickness_from_config(data: object) -> tuple[float, float]:
    if not isinstance(data, dict):
        raise ConfigError("thickness must be a mapping")
    try:
        maximum = validate_positive_finite(data.get("maximum_cm"), "maximum_cm")
        step = validate_positive_finite(data.get("step_cm"), "step_cm")
    except Exception as exc:
        raise ConfigError(str(exc)) from exc
    return maximum, step


def validate_coefficients_config(data: dict[str, Any]) -> dict[str, Any]:
    workflow, output_dir = validate_common(data)
    if workflow != "coefficients":
        raise ConfigError("Expected workflow: coefficients")
    energies = energies_from_config(data.get("energies_kev"), allow_null=True)
    energy_range = energy_range_from_config(data.get("energy_range_kev"))
    if energies is not None and energy_range is not None:
        raise ConfigError("coefficients workflow accepts either energies_kev or energy_range_kev, not both")
    materials_raw = data.get("materials")
    if not isinstance(materials_raw, list) or not materials_raw:
        raise ConfigError("materials must be a non-empty list")
    materials = [
        material_from_dict(item, require_density=False, allow_null_density=True)
        for item in materials_raw
    ]
    ids = [material.id for material in materials]
    if len(set(ids)) != len(ids):
        raise ConfigError("material ids must be unique")
    online_comparison = bool(data.get("online_comparison", False))
    if online_comparison and energy_range is not None:
        raise ConfigError("online_comparison does not support energy_range_kev; use energies_kev or null")
    if online_comparison:
        unsupported_online = [
            material.id for material in materials
            if material.composition_basis != "atomic_ratio"
        ]
        if unsupported_online:
            raise ConfigError(
                "online_comparison supports atomic_ratio materials only; "
                f"unsupported materials: {unsupported_online}"
            )
    precision = str(data.get("precision", "low"))
    valid_precisions = {"direct", "super", "high", "medium", "low", "fast"}
    if precision not in valid_precisions:
        raise ConfigError(f"precision must be one of {sorted(valid_precisions)}")
    return {
        "workflow": workflow,
        "output_dir": output_dir,
        "materials": materials,
        "energies_kev": energies,
        "energy_range_kev": energy_range,
        "online_comparison": online_comparison,
        "precision": precision,
    }


def validate_single_material_thickness_config(data: dict[str, Any], expected_workflow: str) -> dict[str, Any]:
    workflow, output_dir = validate_common(data)
    if workflow != expected_workflow:
        raise ConfigError(f"Expected workflow: {expected_workflow}")
    material = material_from_dict(data.get("material"), require_density=True, allow_null_density=False)
    maximum, step = thickness_from_config(data.get("thickness"))
    return {
        "workflow": workflow,
        "output_dir": output_dir,
        "material": material,
        "energies_kev": energies_from_config(data.get("energies_kev"), allow_null=False),
        "maximum_cm": maximum,
        "step_cm": step,
    }


def validate_multilayer_config(data: dict[str, Any]) -> dict[str, Any]:
    workflow, output_dir = validate_common(data)
    if workflow != "multilayer_attenuation_profile":
        raise ConfigError("Expected workflow: multilayer_attenuation_profile")
    depth_sampling = data.get("depth_sampling")
    if not isinstance(depth_sampling, dict):
        raise ConfigError("depth_sampling must be a mapping")
    try:
        step = validate_positive_finite(depth_sampling.get("step_cm"), "depth_sampling.step_cm")
    except Exception as exc:
        raise ConfigError(str(exc)) from exc
    layers_raw = data.get("layers")
    if not isinstance(layers_raw, list) or not layers_raw:
        raise ConfigError("layers must be a non-empty list")
    layers = [layer_from_dict(item) for item in layers_raw]
    ids = [layer.id for layer in layers]
    if len(set(ids)) != len(ids):
        raise ConfigError("layer ids must be unique")
    return {
        "workflow": workflow,
        "output_dir": output_dir,
        "layers": layers,
        "energies_kev": energies_from_config(data.get("energies_kev"), allow_null=False),
        "step_cm": step,
    }
