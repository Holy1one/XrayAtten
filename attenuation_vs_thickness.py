"""Calculate attenuation versus thickness at selected photon energies."""

from __future__ import annotations

from pathlib import Path
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from support.coefficient_files import (
	interpolate_linear_coefficient,
	provenance_metadata,
	read_coefficient_table,
	write_dataframe_with_metadata,
)


BASE_DIR = Path(__file__).resolve().parent


def build_thickness_grid(maximum_thickness_cm: float, step_cm: float) -> np.ndarray:
	grid = np.arange(0.0, maximum_thickness_cm, step_cm, dtype=float)
	if grid.size == 0 or not math.isclose(grid[-1], maximum_thickness_cm, rel_tol=0.0, abs_tol=1e-14):
		grid = np.append(grid, maximum_thickness_cm)
	return grid


def main(
	coefficient_file: str | Path,
	*,
	energies_kev: list[float] | np.ndarray,
	maximum_thickness_cm: float = 1.0,
	step_cm: float = 1e-4,
	density_g_cm3: float | None = None,
	output_dir: str | Path | None = None,
):
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

	table = read_coefficient_table(coefficient_file)
	linear_mu, effective_density, density_source = interpolate_linear_coefficient(
		table,
		energy / 1000.0,
		expected_kind="attenuation",
		density_g_cm3=density_g_cm3,
	)
	thickness = build_thickness_grid(maximum, step)
	version = "online" if table.metadata.get("schema") == "nist_xcom_attenuation_online" else "local"
	directory = (
		Path(output_dir)
		if output_dir is not None
		else BASE_DIR / "results" / "thickness" / version
	)
	directory.mkdir(parents=True, exist_ok=True)

	figure, axis = plt.subplots(figsize=(9, 6))
	frames: dict[float, pd.DataFrame] = {}
	for energy_value, mu in zip(energy, linear_mu):
		transmission = np.exp(-mu * thickness)
		efficiency = 1.0 - transmission
		frame = pd.DataFrame(
			{
				"Thickness_cm": thickness,
				"Linear_mu_cm_inverse": np.full_like(thickness, mu),
				"Transmission_fraction": transmission,
				"Attenuation_fraction": efficiency,
				"Attenuation_percent": efficiency * 100.0,
			}
		)
		metadata = {
			"schema": "xray_attenuation_vs_thickness",
			"coefficient_schema": table.metadata.get("schema", "unknown"),
			"coefficient_file": str(Path(coefficient_file).resolve()),
			"energy_keV": f"{energy_value:.12g}",
			"effective_density_g_cm3": f"{effective_density:.16g}",
			"density_source": density_source,
			"equation": "attenuation_fraction = 1 - exp(-linear_mu_cm_inverse * thickness_cm)",
			**provenance_metadata(table),
		}
		write_dataframe_with_metadata(directory / f"{energy_value:g}keV_attenuation_{version}.txt", frame, metadata)
		axis.plot(thickness, efficiency * 100.0, label=f"{energy_value:g} keV")
		frames[float(energy_value)] = frame

	axis.set_xlabel("Thickness (cm)")
	axis.set_ylabel("Primary-beam attenuation (%)")
	axis.set_xlim(0, maximum)
	axis.set_ylim(0, 105)
	axis.grid(True, linestyle="--", alpha=0.6)
	axis.legend()
	figure.tight_layout()
	figure.savefig(directory / f"attenuation_vs_thickness_{version}.png", dpi=300)
	figure.savefig(directory / f"attenuation_vs_thickness_{version}.pdf")
	plt.close(figure)
	return frames


if __name__ == "__main__":
	main(
		BASE_DIR / "results" / "coefficients" / "online" / "Si50B20Li20Mg10Y2Cs10Ce1O155F9_attenuation_online.txt",
		energies_kev=list(range(10, 71, 10)),
		maximum_thickness_cm=1.0,
		step_cm=0.001,
	)
