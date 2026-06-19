"""Beer-Lambert transmission operations for online XCOM files."""

from __future__ import annotations

from pathlib import Path
import math

import numpy as np
import pandas as pd

from .attenuation_io_online import (
    CoefficientTable,
    interpolate_linear_attenuation,
    read_online_coefficient_table,
)


def validate_spectrum(energy_kev: object, intensity: object) -> tuple[np.ndarray, np.ndarray]:
    energy = np.asarray(energy_kev, dtype=float)
    values = np.asarray(intensity, dtype=float)
    if energy.ndim != 1 or values.ndim != 1 or energy.size != values.size or energy.size == 0:
        raise ValueError("energy_kev and intensity must be non-empty one-dimensional arrays of equal length")
    if np.any(~np.isfinite(energy)) or np.any(~np.isfinite(values)):
        raise ValueError("Spectrum contains non-finite values")
    if np.any(energy <= 0) or np.any(values < 0):
        raise ValueError("Spectrum energies must be positive and intensities must be non-negative")
    if np.any(np.diff(energy) <= 0):
        raise ValueError("Spectrum energies must be strictly increasing")
    return energy, values


def compute_transmission(
    coefficient_file: str | Path,
    energy_kev: object,
    intensity: object,
    *,
    thickness_cm: float,
    density_g_cm3: float | None = None,
) -> tuple[pd.DataFrame, CoefficientTable]:
    energy, incident = validate_spectrum(energy_kev, intensity)
    thickness = float(thickness_cm)
    if not math.isfinite(thickness) or thickness < 0:
        raise ValueError("thickness_cm must be finite and non-negative")
    table = read_online_coefficient_table(coefficient_file)
    linear_mu = interpolate_linear_attenuation(
        table,
        energy / 1000.0,
        density_g_cm3=density_g_cm3,
    )
    transmission = np.exp(-linear_mu * thickness)
    result = pd.DataFrame(
        {
            "Energy_keV": energy,
            "I0": incident,
            "Linear_mu_cm_inverse": linear_mu,
            "Transmission_fraction": transmission,
            "Attenuation_fraction": 1.0 - transmission,
            "I_after": incident * transmission,
        }
    )
    return result, table


def provenance_metadata(table: CoefficientTable) -> dict[str, str]:
    keys = (
        "source_database",
        "source_version",
        "source_version_date",
        "source_url",
        "request_endpoint",
        "source_accessed_at_utc",
        "source_version_verified",
        "response_sha256",
        "citation",
        "nist_database_disclaimer_url",
    )
    return {key: table.metadata[key] for key in keys if key in table.metadata}


def write_dataframe_with_metadata(
    path: str | Path,
    frame: pd.DataFrame,
    metadata: dict[str, object],
) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="\n") as handle:
        for key, value in metadata.items():
            handle.write(f"# {key}: {value}\n")
        frame.to_csv(handle, sep="\t", index=False, float_format="%.12e")
    return output
