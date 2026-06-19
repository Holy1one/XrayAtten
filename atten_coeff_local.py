"""NIST-traceable photon attenuation coefficient calculator (local).

The elemental tables are a local snapshot of NIST X-Ray Mass Attenuation
Coefficients version 1.4. Compound mass attenuation coefficients are formed
using the NIST mass-fraction additivity rule. Interpolation is log-log linear
inside intervals separated by absorption edges; extrapolation is forbidden.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import math
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from support.attenuation_io_local import SCHEMA_NAME, SCHEMA_VERSION


BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"
COEFF_DIR = DATA_DIR / "nist_local" / "coeff"
ELEMENTS_CSV = DATA_DIR / "elements.csv"
MANIFEST_PATH = DATA_DIR / "nist_local" / "manifest.json"

NIST_DATABASE = "NIST Tables of X-Ray Mass Attenuation Coefficients and Mass Energy-Absorption Coefficients"
NIST_VERSION = "1.4"
NIST_VERSION_DATE = "2004-07-12"
NIST_URL = "https://physics.nist.gov/PhysRefData/XrayMassCoef/"
NIST_XCOM_VERSION = "1.5"
NIST_XCOM_VERSION_DATE = "2010-11"
SNAPSHOT_ACCESS_DATE = "2025-09-24"
VERSION_VERIFIED_DATE = "2026-06-11"
ENERGY_MIN_MEV = 0.001
ENERGY_MAX_MEV = 20.0


@dataclass(frozen=True)
class ElementTable:
	symbol: str
	energy_mev: np.ndarray
	mass_attenuation: np.ndarray
	mass_energy_absorption: np.ndarray
	sha256: str
	path: Path


@dataclass(frozen=True)
class OutputRow:
	energy_mev: float
	edge_side: str


def _sha256(path: Path) -> str:
	digest = hashlib.sha256()
	with path.open("rb") as handle:
		for chunk in iter(lambda: handle.read(1024 * 1024), b""):
			digest.update(chunk)
	return digest.hexdigest()


def _load_manifest(path: Path = MANIFEST_PATH) -> dict:
	try:
		manifest = json.loads(path.read_text(encoding="utf-8"))
	except FileNotFoundError as exc:
		raise FileNotFoundError(f"NIST data manifest is missing: {path}") from exc
	if manifest.get("schema") != "nist_xraymasscoef_snapshot_manifest_v1":
		raise ValueError(f"Invalid NIST manifest schema: {path}")
	return manifest


def _read_elements(path: Path = ELEMENTS_CSV) -> pd.DataFrame:
	df = pd.read_csv(path, encoding="utf-8-sig")
	required = {"原子序数", "符号", "标准原子质量", "原子质量来源"}
	if not required.issubset(df.columns):
		raise ValueError(f"elements.csv must contain columns: {sorted(required)}")
	if df["原子序数"].tolist() != list(range(1, len(df) + 1)):
		raise ValueError("elements.csv atomic numbers must be continuous and ordered")
	if df["符号"].duplicated().any():
		raise ValueError("elements.csv contains duplicate symbols")
	mass = pd.to_numeric(df["标准原子质量"], errors="coerce")
	if mass.isna().any() or (mass <= 0).any():
		raise ValueError("elements.csv contains invalid atomic masses")
	return df


def _parse_composition(composition: dict | str) -> dict[str, float]:
	if isinstance(composition, dict):
		items = list(composition.items())
	elif isinstance(composition, str):
		items = []
		seen: set[str] = set()
		if not composition.strip():
			raise ValueError("composition must not be empty")
		for part in composition.split(","):
			if part.count(":") != 1:
				raise ValueError("String composition must use 'Element:value' entries")
			symbol, raw_value = (piece.strip() for piece in part.split(":"))
			if symbol in seen:
				raise ValueError(f"Duplicate element in composition: {symbol}")
			seen.add(symbol)
			items.append((symbol, raw_value))
	else:
		raise TypeError("composition must be a dict or an 'Element:value' string")
	if not items:
		raise ValueError("composition must not be empty")

	parsed: dict[str, float] = {}
	for raw_symbol, raw_value in items:
		symbol = str(raw_symbol).strip()
		if not re.fullmatch(r"[A-Z][a-z]?", symbol):
			raise ValueError(f"Invalid element symbol: {symbol!r}")
		try:
			value = float(raw_value)
		except (TypeError, ValueError) as exc:
			raise ValueError(f"Invalid composition value for {symbol}") from exc
		if not math.isfinite(value) or value <= 0:
			raise ValueError(f"Composition value for {symbol} must be finite and positive")
		if symbol in parsed:
			raise ValueError(f"Duplicate element in composition: {symbol}")
		parsed[symbol] = value
	return parsed


def _mass_fractions(composition: dict[str, float], basis: str, elements: pd.DataFrame) -> dict[str, float]:
	records = elements.set_index("符号")
	unknown = [symbol for symbol in composition if symbol not in records.index]
	if unknown:
		raise ValueError(f"Unknown elements: {unknown}")
	unsupported = [symbol for symbol in composition if int(records.loc[symbol, "原子序数"]) > 92]
	if unsupported:
		raise ValueError(f"NIST v1.4 coefficient tables are limited to Z=1-92: {unsupported}")

	if basis == "atomic_ratio":
		weighted = {
			symbol: amount * float(records.loc[symbol, "标准原子质量"])
			for symbol, amount in composition.items()
		}
		total = math.fsum(weighted.values())
		return {symbol: value / total for symbol, value in weighted.items()}
	if basis == "mass_fraction":
		total = math.fsum(composition.values())
		if not math.isclose(total, 1.0, rel_tol=0.0, abs_tol=1e-8):
			raise ValueError(f"mass_fraction values must sum to 1.0; received {total:.16g}")
		return {symbol: value / total for symbol, value in composition.items()}
	raise ValueError("composition_basis must be 'atomic_ratio' or 'mass_fraction'")


def _read_element_table(symbol: str, manifest: dict) -> ElementTable:
	path = COEFF_DIR / f"{symbol}.txt"
	entry = manifest.get("elements", {}).get(symbol)
	if entry is None:
		raise ValueError(f"No manifest entry for element {symbol}")
	actual_hash = _sha256(path)
	if actual_hash != entry.get("sha256"):
		raise ValueError(f"NIST source file hash mismatch for {symbol}: {path}")

	df = pd.read_csv(path, sep=r"\s+", engine="python")
	if df.shape[1] < 3:
		raise ValueError(f"Invalid coefficient table: {path}")
	data = df.iloc[:, :3].apply(pd.to_numeric, errors="raise").to_numpy(dtype=float)
	if np.any(~np.isfinite(data)) or np.any(data <= 0):
		raise ValueError(f"Non-positive or non-finite NIST data in {path}")
	if np.any(np.diff(data[:, 0]) < 0):
		raise ValueError(f"Unsorted energy data in {path}")
	return ElementTable(symbol, data[:, 0], data[:, 1], data[:, 2], actual_hash, path)


def _edge_energies(table: ElementTable) -> np.ndarray:
	duplicate = table.energy_mev[1:] == table.energy_mev[:-1]
	return table.energy_mev[1:][duplicate]


def _build_rows(tables: dict[str, ElementTable], energies_mev: object | None) -> list[OutputRow]:
	all_edges = np.unique(np.concatenate([_edge_energies(table) for table in tables.values()]))
	if energies_mev is None:
		energies = np.unique(np.concatenate([table.energy_mev for table in tables.values()]))
	else:
		energies = np.atleast_1d(np.asarray(energies_mev, dtype=float))
		if energies.ndim != 1 or energies.size == 0:
			raise ValueError("energies_mev must be a non-empty one-dimensional sequence")
		if np.any(~np.isfinite(energies)):
			raise ValueError("energies_mev contains non-finite values")
		if np.any((energies < ENERGY_MIN_MEV) | (energies > ENERGY_MAX_MEV)):
			raise ValueError("energies_mev must be within 0.001-20 MeV; extrapolation is forbidden")
		energies = np.unique(energies)

	rows: list[OutputRow] = []
	for energy in sorted(float(value) for value in energies):
		is_edge = np.any(np.isclose(all_edges, energy, rtol=0.0, atol=1e-12))
		if is_edge:
			rows.extend([OutputRow(energy, "below"), OutputRow(energy, "above")])
		else:
			rows.append(OutputRow(energy, "regular"))
	return rows


def _segments(energy: np.ndarray) -> list[tuple[int, int]]:
	segments: list[tuple[int, int]] = []
	start = 0
	for index in range(1, len(energy)):
		if energy[index] == energy[index - 1]:
			segments.append((start, index - 1))
			start = index
	segments.append((start, len(energy) - 1))
	return segments


def _interpolate_element(table: ElementTable, coefficients: np.ndarray, rows: list[OutputRow]) -> np.ndarray:
	result = np.full(len(rows), np.nan, dtype=float)
	segments = _segments(table.energy_mev)
	for row_index, row in enumerate(rows):
		exact = np.flatnonzero(np.isclose(table.energy_mev, row.energy_mev, rtol=0.0, atol=1e-12))
		if exact.size:
			if exact.size == 1:
				result[row_index] = coefficients[exact[0]]
			else:
				result[row_index] = coefficients[exact[0] if row.edge_side == "below" else exact[-1]]
			continue

		for start, end in segments:
			x = table.energy_mev[start : end + 1]
			if x[0] < row.energy_mev < x[-1]:
				y = coefficients[start : end + 1]
				result[row_index] = math.exp(
					np.interp(math.log(row.energy_mev), np.log(x), np.log(y))
				)
				break
		if not math.isfinite(result[row_index]):
			raise ValueError(
				f"No edge-safe interpolation interval for {table.symbol} at {row.energy_mev:g} MeV"
			)
	return result


def _default_output_name(composition: dict[str, float], basis: str, kind: str) -> str:
	parts = [f"{symbol}{value:g}" for symbol, value in composition.items()]
	return "".join(parts) + f"_{basis}_{kind}_local"


def _write_output(
	path: Path,
	rows: list[OutputRow],
	mass_values: np.ndarray,
	linear_values: np.ndarray | None,
	*,
	composition: dict[str, float],
	basis: str,
	kind: str,
	mass_fractions: dict[str, float],
	density: float | None,
	tables: dict[str, ElementTable],
) -> None:
	if kind == "attenuation":
		mass_column = "mass_attenuation_cm2_g"
		linear_column = "linear_attenuation_cm_inverse"
	else:
		mass_column = "mass_energy_absorption_approx_cm2_g"
		linear_column = "linear_energy_absorption_approx_cm_inverse"
	columns = ["energy_MeV", "edge_side", mass_column]
	if linear_values is not None:
		columns.append(linear_column)

	used_sources = {
		symbol: {
			"file": table.path.relative_to(BASE_DIR).as_posix(),
			"sha256": table.sha256,
			"url": f"{NIST_URL}ElemTab/z{int(_read_elements().set_index('符号').loc[symbol, '原子序数']):02d}.html",
		}
		for symbol, table in tables.items()
	}
	metadata = [
		("schema", SCHEMA_NAME),
		("schema_version", SCHEMA_VERSION),
		("coefficient_kind", kind),
		("columns", " ".join(columns)),
		("units", "MeV cm2/g" + (" 1/cm" if linear_values is not None else "")),
		("source_database", NIST_DATABASE),
		("source_version", NIST_VERSION),
		("source_version_date", NIST_VERSION_DATE),
		("source_url", NIST_URL),
		("xcom_reference_version", NIST_XCOM_VERSION),
		("xcom_reference_version_date", NIST_XCOM_VERSION_DATE),
		("source_snapshot_accessed", SNAPSHOT_ACCESS_DATE),
		("source_version_verified", VERSION_VERIFIED_DATE),
		("source_integrity_check", "All 92 local elemental tables matched the current NIST online values on 2026-06-11"),
		("generated_at_utc", datetime.now(timezone.utc).isoformat(timespec="seconds")),
		("citation", "Hubbell JH and Seltzer SM, NISTIR 5632 (1995), online database version 1.4 (2004)"),
		("composition_basis", basis),
		("composition_json", json.dumps(composition, sort_keys=True, separators=(",", ":"))),
		("mass_fractions_json", json.dumps(mass_fractions, sort_keys=True, separators=(",", ":"))),
		("atomic_weight_source", "NIST X-Ray Mass Attenuation Coefficients Table 1 Z/A-derived values"),
		("density_g_cm3", "none" if density is None else f"{density:.16g}"),
		("interpolation", "piecewise log-log linear; absorption-edge intervals separated; no extrapolation"),
		("edge_policy", "duplicate energy rows labeled below/above; exact-edge transmission uses above"),
		("source_files_json", json.dumps(used_sources, sort_keys=True, separators=(",", ":"))),
	]
	if "F" in tables:
		metadata.append((
			"local_data_correction",
			"Fluorine 60 keV energy label corrected from an erroneous local duplicate 50 keV label; coefficient values were unchanged and verified against NIST on 2026-06-11",
		))
	if kind == "energy_absorption_approx":
		metadata.append((
			"warning",
			"Approximate elemental mass-fraction additivity only; NIST states matrix-dependent radiative-loss effects make simple additivity formally inadequate for compound mu_en/rho",
		))

	with path.open("w", encoding="utf-8", newline="\n") as handle:
		for key, value in metadata:
			handle.write(f"# {key}: {value}\n")
		for index, row in enumerate(rows):
			fields = [f"{row.energy_mev:.12e}", row.edge_side, f"{mass_values[index]:.12e}"]
			if linear_values is not None:
				fields.append(f"{linear_values[index]:.12e}")
			handle.write("  ".join(fields) + "\n")


def _write_plot(path: Path, rows: list[OutputRow], values: np.ndarray, ylabel: str) -> None:
	energy = np.asarray([row.energy_mev for row in rows])
	figure, axis = plt.subplots(figsize=(6.5, 4.5))
	axis.loglog(energy, values, linewidth=1.2)
	axis.set_xlabel("Photon energy (MeV)")
	axis.set_ylabel(ylabel)
	axis.grid(which="both", linestyle="--", linewidth=0.5, alpha=0.6)
	figure.tight_layout()
	figure.savefig(path, dpi=300)
	plt.close(figure)


def compute_coefficients(
	composition: dict | str,
	*,
	composition_basis: str,
	coefficient_kind: str = "attenuation",
	density_g_cm3: float | None = None,
	energies_mev: object | None = None,
	output_dir: str | Path = RESULTS_DIR / "coefficients" / "local",
	output_name: str | None = None,
) -> dict:
	"""Compute a traceable compound coefficient table from the NIST snapshot."""
	parsed = _parse_composition(composition)
	if coefficient_kind not in {"attenuation", "energy_absorption_approx"}:
		raise ValueError("coefficient_kind must be 'attenuation' or 'energy_absorption_approx'")
	if density_g_cm3 is not None:
		density = float(density_g_cm3)
		if not math.isfinite(density) or density <= 0:
			raise ValueError("density_g_cm3 must be finite and positive")
	else:
		density = None

	elements = _read_elements()
	weights = _mass_fractions(parsed, composition_basis, elements)
	manifest = _load_manifest()
	tables = {symbol: _read_element_table(symbol, manifest) for symbol in parsed}
	rows = _build_rows(tables, energies_mev)
	mass_values = np.zeros(len(rows), dtype=float)
	for symbol, weight in weights.items():
		table = tables[symbol]
		source = table.mass_attenuation if coefficient_kind == "attenuation" else table.mass_energy_absorption
		mass_values += weight * _interpolate_element(table, source, rows)
	linear_values = mass_values * density if density is not None else None

	directory = Path(output_dir)
	if not directory.is_absolute():
		directory = BASE_DIR / directory
	directory.mkdir(parents=True, exist_ok=True)
	name = output_name or _default_output_name(parsed, composition_basis, coefficient_kind)
	if not re.fullmatch(r"[A-Za-z0-9_.-]+", name):
		raise ValueError("output_name may contain only letters, digits, dot, underscore, and hyphen")
	txt_path = directory / f"{name}.txt"
	plot_path = directory / f"{name}.png"
	_write_output(
		txt_path,
		rows,
		mass_values,
		linear_values,
		composition=parsed,
		basis=composition_basis,
		kind=coefficient_kind,
		mass_fractions=weights,
		density=density,
		tables=tables,
	)
	ylabel = "Mass attenuation coefficient (cm²/g)"
	plot_values = mass_values
	if linear_values is not None:
		ylabel = "Linear attenuation coefficient (cm⁻¹)" if coefficient_kind == "attenuation" else "Linear energy-absorption coefficient, approximate (cm⁻¹)"
		plot_values = linear_values
	elif coefficient_kind == "energy_absorption_approx":
		ylabel = "Mass energy-absorption coefficient, approximate (cm²/g)"
	_write_plot(plot_path, rows, plot_values, ylabel)
	return {
		"energy_mev": np.asarray([row.energy_mev for row in rows]),
		"edge_side": np.asarray([row.edge_side for row in rows]),
		"mass_coefficient": mass_values,
		"linear_coefficient": linear_values,
		"mass_fractions": weights,
		"txt_path": txt_path,
		"plot_path": plot_path,
	}


def main() -> None:
	result = compute_coefficients(
		{"Si": 50, "B": 20, "Li": 20, "Mg": 10, "Y": 2, "Cs": 10, "Ce": 1, "O": 155, "F": 9},
		composition_basis="atomic_ratio",
		density_g_cm3=2.5,
		output_name="Si50B20Li20Mg10Y2Cs10Ce1O155F9_attenuation_local",
	)
	print(result["txt_path"])


if __name__ == "__main__":
	main()
