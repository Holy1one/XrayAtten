"""Sequential X-ray spectrum transmission through multiple material layers."""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BASE_DIR = Path(__file__).resolve().parent

from support.coefficient_files import (
	interpolate_linear_attenuation,
	provenance_metadata,
	read_coefficient_table,
	write_dataframe_with_metadata,
)


def simulate_layers(
	layers: list[dict],
	energy_kev,
	incident_intensity,
	*,
	output_dir: str | Path = BASE_DIR / "results" / "multilayer" / "online",
) -> list[pd.DataFrame]:
	if not layers:
		raise ValueError("At least one layer is required")
	directory = Path(output_dir)
	directory.mkdir(parents=True, exist_ok=True)

	current_intensity = incident_intensity
	results: list[pd.DataFrame] = []
	figure, axis = plt.subplots(figsize=(9, 6))
	axis.plot(energy_kev, incident_intensity, label="Incident spectrum", linewidth=1.5)

	for index, layer in enumerate(layers, start=1):
		name = str(layer["name"])
		coefficient_file = Path(layer["coefficient_file"])
		thickness_cm = float(layer["thickness_cm"])
		density = layer.get("density_g_cm3")
		table = read_coefficient_table(coefficient_file)
		linear_mu = interpolate_linear_attenuation(
			table,
			np.asarray(energy_kev, dtype=float) / 1000.0,
			density_g_cm3=density,
		)
		transmission = np.exp(-linear_mu * thickness_cm)
		frame = pd.DataFrame(
			{
				"Energy_keV": np.asarray(energy_kev, dtype=float),
				"I0": np.asarray(current_intensity, dtype=float),
				"Linear_mu_cm_inverse": linear_mu,
				"Transmission_fraction": transmission,
				"Attenuation_fraction": 1.0 - transmission,
				"I_after": np.asarray(current_intensity, dtype=float) * transmission,
			}
		)
		metadata = {
			"schema": "xray_multilayer_spectrum_stage",
			"coefficient_schema": table.metadata.get("schema", "unknown"),
			"layer_index": index,
			"layer_name": name,
			"coefficient_file": str(coefficient_file.resolve()),
			"thickness_cm": f"{thickness_cm:.16g}",
			"density_override_g_cm3": "none" if density is None else f"{float(density):.16g}",
			"equation": "I_after = I0 * exp(-linear_mu_cm_inverse * thickness_cm)",
			**provenance_metadata(table),
		}
		write_dataframe_with_metadata(directory / f"layer_{index}_{name}.txt", frame, metadata)
		axis.plot(frame["Energy_keV"], frame["I_after"], label=f"After layer {index}: {name}")
		current_intensity = frame["I_after"].to_numpy(dtype=float)
		results.append(frame)

	axis.set_xlabel("Energy (keV)")
	axis.set_ylabel("Intensity")
	axis.set_title("Sequential X-ray transmission through material layers")
	axis.grid(True, linestyle="--", alpha=0.5)
	axis.legend()
	figure.tight_layout()
	figure.savefig(directory / "multilayer_spectrum.png", dpi=300)
	plt.close(figure)
	return results


def main() -> None:
	spectrum = pd.read_csv(BASE_DIR / "data" / "spectra" / "70.csv")
	energy = spectrum.iloc[:, 0].to_numpy(dtype=float)
	intensity = spectrum.iloc[:, 1].to_numpy(dtype=float)
	valid = (energy >= 1.0) & (energy <= 70.0)
	layers = [
		{
			"name": "Red",
			"coefficient_file": BASE_DIR / "results" / "coefficients" / "online" / "C420H400Br42P20Mn8Sb2_attenuation_online.txt",
			"thickness_cm": 0.1,
		},
		{
			"name": "Green",
			"coefficient_file": BASE_DIR / "results" / "coefficients" / "online" / "Si50Ca10F30Na10Al20Mn11O1311_attenuation_online.txt",
			"thickness_cm": 0.1,
		},
		{
			"name": "Blue",
			"coefficient_file": BASE_DIR / "results" / "coefficients" / "online" / "Si50B20Li20Mg10Y2Cs10Ce1O155F9_attenuation_online.txt",
			"thickness_cm": 0.1,
		},
	]
	simulate_layers(
		layers,
		energy[valid],
		intensity[valid],
	)


if __name__ == "__main__":
	main()
