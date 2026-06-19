"""Version-neutral access to local and online attenuation coefficient files."""

from __future__ import annotations

from pathlib import Path
import math

import numpy as np

from .attenuation_io_local import (
    interpolate_linear_attenuation as interpolate_local,
    read_local_coefficient_table,
)
from .attenuation_io_online import (
    interpolate_linear_attenuation as interpolate_online,
    read_online_coefficient_table,
)


LOCAL_SCHEMA = "nist_attenuation_local"
ONLINE_SCHEMA = "nist_xcom_attenuation_online"


def coefficient_schema(path: str | Path) -> str:
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.startswith("#"):
            break
        if line.startswith("# schema:"):
            return line.split(":", 1)[1].strip()
    raise ValueError(f"Coefficient schema metadata is missing: {path}")


def read_coefficient_table(path: str | Path):
    schema = coefficient_schema(path)
    if schema == LOCAL_SCHEMA:
        return read_local_coefficient_table(path)
    if schema == ONLINE_SCHEMA:
        return read_online_coefficient_table(path)
    raise ValueError(f"Unsupported coefficient schema {schema!r}: {path}")


def interpolate_linear_attenuation(table, target_energy_mev, *, density_g_cm3=None):
    schema = table.metadata.get("schema")
    if schema == LOCAL_SCHEMA:
        return interpolate_local(table, target_energy_mev, density_g_cm3=density_g_cm3)
    if schema == ONLINE_SCHEMA:
        return interpolate_online(table, target_energy_mev, density_g_cm3=density_g_cm3)
    raise ValueError(f"Unsupported coefficient schema {schema!r}")


def embedded_density_g_cm3(table) -> float | None:
    value = table.metadata.get("density_g_cm3", "none")
    if value.lower() == "none":
        return None
    density = float(value)
    if not math.isfinite(density) or density <= 0:
        raise ValueError(f"Invalid embedded density in {table.path}")
    return density


def resolve_density_g_cm3(table, density_g_cm3: float | None) -> tuple[float, str]:
    embedded = embedded_density_g_cm3(table)
    if embedded is not None:
        if density_g_cm3 is not None:
            requested = float(density_g_cm3)
            if not math.isclose(requested, embedded, rel_tol=1e-12, abs_tol=0.0):
                raise ValueError(
                    f"Provided density {requested:g} g/cm3 conflicts with embedded density "
                    f"{embedded:g} g/cm3 in {table.path}"
                )
        return embedded, "embedded_in_coefficient_file"
    if density_g_cm3 is None:
        raise ValueError(
            f"density_g_cm3 is required because {table.path} contains mass coefficients only"
        )
    density = float(density_g_cm3)
    if not math.isfinite(density) or density <= 0:
        raise ValueError("density_g_cm3 must be finite and positive")
    return density, "function_argument"


def _segments(energy: np.ndarray) -> list[tuple[int, int]]:
    segments: list[tuple[int, int]] = []
    start = 0
    for index in range(1, len(energy)):
        if energy[index] == energy[index - 1]:
            segments.append((start, index - 1))
            start = index
    segments.append((start, len(energy) - 1))
    return segments


def interpolate_table_rows(table, target_energy_mev, row_values) -> np.ndarray:
    target = np.atleast_1d(np.asarray(target_energy_mev, dtype=float))
    values = np.asarray(row_values, dtype=float)
    if values.shape != table.energy_mev.shape:
        raise ValueError("row_values must match the coefficient table length")
    if np.any(~np.isfinite(target)) or np.any(target <= 0):
        raise ValueError("Target energies must be finite and positive")
    if np.any(target < table.energy_mev[0]) or np.any(target > table.energy_mev[-1]):
        raise ValueError("Target energy is outside the coefficient table; extrapolation is forbidden")
    if np.any(~np.isfinite(values)) or np.any(values <= 0):
        raise ValueError("Coefficient rows must be finite and positive")

    result = np.full(target.shape, np.nan, dtype=float)
    for index, energy in enumerate(target):
        exact = np.flatnonzero(np.isclose(table.energy_mev, energy, rtol=0.0, atol=1e-12))
        if exact.size:
            above = exact[table.edge_side[exact] == "above"]
            result[index] = values[above[-1] if above.size else exact[-1]]
            continue
        for start, end in _segments(table.energy_mev):
            x = table.energy_mev[start : end + 1]
            if x[0] < energy < x[-1]:
                y = values[start : end + 1]
                result[index] = math.exp(np.interp(math.log(energy), np.log(x), np.log(y)))
                break
        if not math.isfinite(result[index]):
            raise ValueError(f"No valid interpolation segment for {energy:g} MeV")
    return result


def interpolate_linear_coefficient(
    table,
    target_energy_mev,
    *,
    expected_kind: str,
    density_g_cm3: float | None = None,
) -> tuple[np.ndarray, float, str]:
    if table.coefficient_kind != expected_kind:
        raise ValueError(
            f"Expected coefficient_kind={expected_kind!r}, received {table.coefficient_kind!r}"
        )
    density, density_source = resolve_density_g_cm3(table, density_g_cm3)
    mass_values = interpolate_table_rows(table, target_energy_mev, table.mass_coefficient)
    linear_values = mass_values * density
    if table.linear_coefficient is not None:
        stored = interpolate_table_rows(table, target_energy_mev, table.linear_coefficient)
        if not np.allclose(linear_values, stored, rtol=5e-12, atol=0.0):
            raise ValueError(f"Stored linear coefficients are inconsistent in {table.path}")
        linear_values = stored
    return linear_values, density, density_source


def provenance_metadata(table) -> dict[str, str]:
    keys = (
        "source_database",
        "source_version",
        "source_version_date",
        "source_url",
        "request_endpoint",
        "source_snapshot_accessed",
        "source_accessed_at_utc",
        "source_version_verified",
        "source_integrity_check",
        "response_sha256",
        "local_data_correction",
        "warning",
        "composition_basis",
        "composition_json",
        "mass_fractions_json",
        "citation",
        "nist_database_disclaimer_url",
    )
    return {key: table.metadata[key] for key in keys if key in table.metadata}


def write_dataframe_with_metadata(path, frame, metadata) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="\n") as handle:
        for key, value in metadata.items():
            handle.write(f"# {key}: {value}\n")
        frame.to_csv(handle, sep="\t", index=False, float_format="%.12e")
    return output
