"""TXT output helpers with metadata headers."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from .interpolation import build_compound_dense_grid, interpolate_element_to_grid
from .models import CoefficientResult
from .provenance import runtime_metadata


def write_dataframe(path: Path, frame: pd.DataFrame, metadata: dict[str, object]) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    merged = {**metadata, **runtime_metadata()}
    with output.open("w", encoding="utf-8", newline="\n") as handle:
        for key, value in merged.items():
            handle.write(f"# {key}: {value}\n")
        frame.to_csv(handle, sep="\t", index=False, float_format="%.12e")
    return output


def coefficient_columns(result: CoefficientResult) -> list[str]:
    if result.coefficient_kind == "attenuation":
        columns = ["Energy_keV", "Edge_side", "Mass_mu_over_rho_cm2_g"]
        if result.linear_coefficient is not None:
            columns.append("Linear_mu_cm_inverse")
        return columns
    columns = ["Energy_keV", "Edge_side", "Mass_mu_en_over_rho_approx_cm2_g"]
    if result.linear_coefficient is not None:
        columns.append("Linear_mu_en_approx_cm_inverse")
    return columns


def coefficient_frame(result: CoefficientResult) -> pd.DataFrame:
    columns = coefficient_columns(result)
    data: dict[str, object] = {
        "Energy_keV": [row.energy_kev for row in result.rows],
        "Edge_side": [row.edge_side for row in result.rows],
        columns[2]: result.mass_coefficient,
    }
    if result.linear_coefficient is not None:
        data[columns[3]] = result.linear_coefficient
    return pd.DataFrame(data)


def write_coefficient_result(path: Path, result: CoefficientResult) -> Path:
    columns = coefficient_columns(result)
    metadata = {
        **result.metadata,
        "workflow_columns": " ".join(columns),
        "warnings": " | ".join(result.warnings) if result.warnings else "none",
        "limitations": " | ".join(result.limitations) if result.limitations else "none",
    }
    return write_dataframe(path, coefficient_frame(result), metadata)


def write_dense_coefficient_result(
    path: Path,
    result: CoefficientResult,
    precision: str = "low",
) -> Path:
    """Write coefficient TXT with dense log-log interpolation for smooth Origin plots.

    Strategy: interpolate each element independently to a unified dense grid
    (preserving all absorption edges), then combine via mass-fraction additivity.
    """
    tables = result.element_tables
    if not tables:
        from .local_nist import load_material_tables
        _, tables = load_material_tables(result.material)

    weights = result.mass_fractions

    # Build unified dense grid from ALL element tables (preserves every edge)
    grid_kev, edge_sides = build_compound_dense_grid(tables, precision)

    # Interpolate each element to the dense grid, then combine
    dense_mass = np.zeros(len(grid_kev), dtype=float)
    for symbol, weight in weights.items():
        table = tables[symbol]
        source = (
            table.mass_attenuation
            if result.coefficient_kind == "attenuation"
            else table.mass_energy_absorption
        )
        dense_mass += weight * interpolate_element_to_grid(table, source, grid_kev, edge_sides)

    # Build data
    data: dict[str, object] = {
        "Energy_keV": grid_kev,
        "Edge_side": edge_sides,
    }
    if result.coefficient_kind == "attenuation":
        data["Mass_mu_over_rho_cm2_g"] = dense_mass
    else:
        data["Mass_mu_en_over_rho_approx_cm2_g"] = dense_mass

    # Linear coefficients if density available
    density = result.material.density_g_cm3
    if density is not None:
        if result.coefficient_kind == "attenuation":
            data["Linear_mu_cm_inverse"] = dense_mass * density
        else:
            data["Linear_mu_en_approx_cm_inverse"] = dense_mass * density

    frame = pd.DataFrame(data)

    columns = list(data.keys())
    metadata = {
        **result.metadata,
        "interpolation": "per-element dense log-log, then mass-fraction additivity",
        "grid_strategy": "all NIST points + all element edges + gap fill",
        "precision": precision,
        "grid_points": str(len(grid_kev)),
        "workflow_columns": " ".join(columns),
        "warnings": " | ".join(result.warnings) if result.warnings else "none",
        "limitations": " | ".join(result.limitations) if result.limitations else "none",
    }
    return write_dataframe(path, frame, metadata)


def write_text_table(
    path: Path,
    columns: Sequence[str],
    rows: Iterable[Sequence[object]],
    metadata: dict[str, object],
) -> Path:
    frame = pd.DataFrame(list(rows), columns=list(columns))
    return write_dataframe(path, frame, metadata)

