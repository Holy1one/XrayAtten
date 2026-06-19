"""Estimate photon energy absorption versus material thickness.

This workflow uses total attenuation ``mu`` and mass energy-absorption
``mu_en`` coefficients together. For a monoenergetic narrow primary beam, the
first-collision absorbed-energy fraction is estimated as::

    (mu_en / mu) * (1 - exp(-mu * thickness))

It is not a full photon/electron transport calculation and does not include
scattered-photon buildup or geometry-dependent escape and reabsorption.
"""

from __future__ import annotations

from pathlib import Path
import json
import math
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from atten_coeff_local import compute_coefficients
from attenuation_vs_thickness import build_thickness_grid
from support.coefficient_files import (
    interpolate_linear_coefficient,
    provenance_metadata,
    read_coefficient_table,
    write_dataframe_with_metadata,
)


BASE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = BASE_DIR / "results"

DEFAULT_COMPOSITION = {
    "Si": 50,
    "B": 20,
    "Li": 20,
    "Mg": 10,
    "Y": 2,
    "Cs": 10,
    "Ce": 1,
    "O": 155,
    "F": 9,
}
DEFAULT_DENSITY_G_CM3 = 2.5
DEFAULT_NAME = "Si50B20Li20Mg10Y2Cs10Ce1O155F9"


def _metadata_json(table, key: str) -> object | None:
    value = table.metadata.get(key)
    return None if value is None else json.loads(value)


def _validate_coefficient_pair(attenuation_table, absorption_table) -> None:
    if attenuation_table.coefficient_kind != "attenuation":
        raise ValueError("The attenuation file must contain total attenuation coefficients")
    if absorption_table.coefficient_kind != "energy_absorption_approx":
        raise ValueError(
            "The absorption file must contain local mass energy-absorption coefficients; "
            "NIST XCOM online attenuation output cannot be used as mu_en/rho"
        )
    if absorption_table.metadata.get("schema") != "nist_attenuation_local":
        raise ValueError("Energy-absorption coefficients are currently supported only by local")

    for key in ("composition_json", "mass_fractions_json"):
        left = _metadata_json(attenuation_table, key)
        right = _metadata_json(absorption_table, key)
        if left is not None and right is not None and left != right:
            raise ValueError(f"Attenuation and absorption files describe different materials ({key})")


def main(
    attenuation_coefficient_file: str | Path,
    energy_absorption_coefficient_file: str | Path,
    *,
    energies_kev: list[float] | np.ndarray,
    maximum_thickness_cm: float = 1.0,
    step_cm: float = 1e-4,
    density_g_cm3: float | None = None,
    output_dir: str | Path | None = None,
) -> dict[float, pd.DataFrame]:
    maximum = float(maximum_thickness_cm)
    step = float(step_cm)
    if not math.isfinite(maximum) or maximum <= 0:
        raise ValueError("maximum_thickness_cm must be finite and positive")
    if not math.isfinite(step) or step <= 0:
        raise ValueError("step_cm must be finite and positive")
    energy = np.asarray(energies_kev, dtype=float)
    if energy.ndim != 1 or energy.size == 0 or np.any(~np.isfinite(energy)) or np.any(energy <= 0):
        raise ValueError("energies_kev must be a non-empty array of positive finite values")
    if np.any(np.diff(energy) <= 0):
        raise ValueError("energies_kev must be strictly increasing without duplicates")

    attenuation_table = read_coefficient_table(attenuation_coefficient_file)
    absorption_table = read_coefficient_table(energy_absorption_coefficient_file)
    _validate_coefficient_pair(attenuation_table, absorption_table)

    linear_mu, attenuation_density, attenuation_density_source = interpolate_linear_coefficient(
        attenuation_table,
        energy / 1000.0,
        expected_kind="attenuation",
        density_g_cm3=density_g_cm3,
    )
    linear_mu_en, absorption_density, absorption_density_source = interpolate_linear_coefficient(
        absorption_table,
        energy / 1000.0,
        expected_kind="energy_absorption_approx",
        density_g_cm3=density_g_cm3,
    )
    if not math.isclose(attenuation_density, absorption_density, rel_tol=1e-12, abs_tol=0.0):
        raise ValueError("Attenuation and energy-absorption files use different densities")
    if np.any(linear_mu_en > linear_mu * (1.0 + 1e-10)):
        raise ValueError("mu_en must not exceed total attenuation mu")

    warnings.warn(
        "The plotted quantity is a first-collision, narrow-beam absorbed-energy estimate. "
        "It is not a full transport result with photon buildup or geometry-dependent reabsorption.",
        RuntimeWarning,
        stacklevel=2,
    )
    thickness = build_thickness_grid(maximum, step)
    directory = (
        Path(output_dir)
        if output_dir is not None
        else RESULTS_DIR / "energy_absorption" / "local_approx"
    )
    directory.mkdir(parents=True, exist_ok=True)

    figure, axis = plt.subplots(figsize=(9, 6))
    frames: dict[float, pd.DataFrame] = {}
    for energy_value, mu, mu_en in zip(energy, linear_mu, linear_mu_en):
        primary_transmission = np.exp(-mu * thickness)
        primary_removal = 1.0 - primary_transmission
        absorption_to_removal_ratio = mu_en / mu
        absorbed_energy = absorption_to_removal_ratio * primary_removal
        removed_but_not_absorbed = primary_removal - absorbed_energy
        frame = pd.DataFrame(
            {
                "Thickness_cm": thickness,
                "Linear_mu_cm_inverse": np.full_like(thickness, mu),
                "Linear_mu_en_cm_inverse": np.full_like(thickness, mu_en),
                "Primary_transmission_fraction": primary_transmission,
                "Primary_removal_fraction": primary_removal,
                "First_collision_absorbed_energy_fraction": absorbed_energy,
                "First_collision_absorbed_energy_percent": absorbed_energy * 100.0,
                "Removed_primary_energy_not_locally_absorbed_fraction": removed_but_not_absorbed,
            }
        )
        metadata = {
            "schema": "xray_first_collision_energy_absorption_vs_thickness_v1",
            "attenuation_coefficient_file": str(Path(attenuation_coefficient_file).resolve()),
            "energy_absorption_coefficient_file": str(Path(energy_absorption_coefficient_file).resolve()),
            "energy_keV": f"{energy_value:.12g}",
            "effective_density_g_cm3": f"{attenuation_density:.16g}",
            "attenuation_density_source": attenuation_density_source,
            "absorption_density_source": absorption_density_source,
            "equation": "absorbed_energy_fraction = (mu_en / mu) * (1 - exp(-mu * thickness_cm))",
            "model_derivation": "integral from 0 to thickness of mu_en * exp(-mu * x) dx",
            "scope": "monoenergetic narrow primary beam; first-collision energy absorption",
            "limitations": "no scattered-photon buildup; no geometry-dependent secondary-photon escape or reabsorption; no coupled photon-electron transport",
            "mu_en_material_model": "elemental mass-fraction additivity approximation for custom compounds",
            "nist_mu_en_definition_url": "https://physics.nist.gov/PhysRefData/XrayMassCoef/chap3.html",
            "nist_xcom_scope_url": "https://physics.nist.gov/PhysRefData/Xcom/Text/intro.html",
            **{f"absorption_{key}": value for key, value in provenance_metadata(absorption_table).items()},
        }
        write_dataframe_with_metadata(
            directory / f"{energy_value:g}keV_energy_absorption_local_approx.txt",
            frame,
            metadata,
        )
        axis.plot(thickness, absorbed_energy * 100.0, label=f"{energy_value:g} keV")
        frames[float(energy_value)] = frame

    axis.set_xlabel("Thickness (cm)")
    axis.set_ylabel("First-collision absorbed energy estimate (%)")
    axis.set_xlim(0, maximum)
    axis.set_ylim(0, 105)
    axis.grid(True, linestyle="--", alpha=0.6)
    axis.legend()
    figure.tight_layout()
    figure.savefig(directory / "energy_absorption_vs_thickness_local_approx.png", dpi=300)
    figure.savefig(directory / "energy_absorption_vs_thickness_local_approx.pdf")
    plt.close(figure)
    return frames


def prepare_default_coefficient_files() -> tuple[Path, Path]:
    coefficient_dir = RESULTS_DIR / "coefficients" / "local"
    attenuation = compute_coefficients(
        DEFAULT_COMPOSITION,
        composition_basis="atomic_ratio",
        coefficient_kind="attenuation",
        density_g_cm3=DEFAULT_DENSITY_G_CM3,
        output_dir=coefficient_dir,
        output_name=f"{DEFAULT_NAME}_attenuation_local",
    )
    absorption = compute_coefficients(
        DEFAULT_COMPOSITION,
        composition_basis="atomic_ratio",
        coefficient_kind="energy_absorption_approx",
        density_g_cm3=DEFAULT_DENSITY_G_CM3,
        output_dir=coefficient_dir,
        output_name=f"{DEFAULT_NAME}_energy_absorption_approx_local",
    )
    return Path(attenuation["txt_path"]), Path(absorption["txt_path"])


if __name__ == "__main__":
    attenuation_file, absorption_file = prepare_default_coefficient_files()
    main(
        attenuation_file,
        absorption_file,
        energies_kev=list(range(10, 71, 10)),
        maximum_thickness_cm=1.0,
        step_cm=0.001,
    )
