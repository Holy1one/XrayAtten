"""Strict reader and interpolation helpers for online XCOM text files."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json
import math

import numpy as np


SCHEMA_NAME = "nist_xcom_attenuation_online"
SCHEMA_VERSION = "2"


@dataclass(frozen=True)
class CoefficientTable:
	path: Path
	metadata: dict[str, str]
	columns: tuple[str, ...]
	energy_mev: np.ndarray
	edge_side: np.ndarray
	mass_coefficient: np.ndarray
	linear_coefficient: np.ndarray | None
	coefficient_kind: str


def _parse_metadata(lines: list[str]) -> tuple[dict[str, str], int]:
	metadata: dict[str, str] = {}
	data_start = 0
	for index, line in enumerate(lines):
		if not line.startswith("#"):
			data_start = index
			break
		text = line[1:].strip()
		if ":" in text:
			key, value = text.split(":", 1)
			metadata[key.strip()] = value.strip()
	else:
		data_start = len(lines)
	return metadata, data_start


def read_online_coefficient_table(path: str | Path) -> CoefficientTable:
	"""Read an online XCOM coefficient table."""
	file_path = Path(path)
	lines = file_path.read_text(encoding="utf-8").splitlines()
	metadata, data_start = _parse_metadata(lines)

	if metadata.get("schema") != SCHEMA_NAME:
		raise ValueError(
			f"{file_path} is not a {SCHEMA_NAME} file; legacy coefficient files are not supported"
		)
	if metadata.get("schema_version") != SCHEMA_VERSION:
		raise ValueError(f"Unsupported schema_version in {file_path}")

	columns = tuple(metadata.get("columns", "").split())
	if len(columns) not in (2, 3):
		raise ValueError(f"Invalid columns metadata in {file_path}")
	if columns[:2] != ("energy_MeV", "mass_attenuation_cm2_g"):
		raise ValueError(f"Invalid coefficient columns in {file_path}")

	kind = metadata.get("coefficient_kind", "")
	if kind != "attenuation":
		raise ValueError(f"Unknown coefficient_kind in {file_path}: {kind!r}")
	if len(columns) == 3 and columns[2] != "linear_attenuation_cm_inverse":
		raise ValueError(f"Coefficient columns do not match coefficient_kind in {file_path}")

	energy: list[float] = []
	mass: list[float] = []
	linear: list[float] = []
	for line_number, line in enumerate(lines[data_start:], start=data_start + 1):
		if not line.strip() or line.lstrip().startswith("#"):
			continue
		parts = line.split()
		if len(parts) != len(columns):
			raise ValueError(f"Invalid data row at {file_path}:{line_number}")
		try:
			energy_value = float(parts[0])
			mass_value = float(parts[1])
			linear_value = float(parts[2]) if len(columns) == 3 else None
		except ValueError as exc:
			raise ValueError(f"Non-numeric coefficient at {file_path}:{line_number}") from exc
		if not math.isfinite(energy_value) or not math.isfinite(mass_value):
			raise ValueError(f"Non-finite value at {file_path}:{line_number}")
		if energy_value <= 0 or mass_value <= 0:
			raise ValueError(f"Non-positive value at {file_path}:{line_number}")
		if linear_value is not None and (not math.isfinite(linear_value) or linear_value <= 0):
			raise ValueError(f"Invalid linear coefficient at {file_path}:{line_number}")
		energy.append(energy_value)
		mass.append(mass_value)
		if linear_value is not None:
			linear.append(linear_value)

	if not energy:
		raise ValueError(f"No coefficient rows found in {file_path}")
	energy_array = np.asarray(energy, dtype=float)
	if np.any(np.diff(energy_array) < 0):
		raise ValueError(f"Energy values are not sorted in {file_path}")
	edge_side = _derive_edge_sides(energy_array, file_path)

	table = CoefficientTable(
		path=file_path,
		metadata=metadata,
		columns=columns,
		energy_mev=energy_array,
		edge_side=edge_side,
		mass_coefficient=np.asarray(mass, dtype=float),
		linear_coefficient=np.asarray(linear, dtype=float) if linear else None,
		coefficient_kind=kind,
	)
	_validate_linear_relation(table)
	return table


def _derive_edge_sides(energy: np.ndarray, path: Path) -> np.ndarray:
	sides = np.full(energy.shape, "regular", dtype="<U7")
	index = 0
	while index < energy.size:
		end = index + 1
		while end < energy.size and energy[end] == energy[index]:
			end += 1
		count = end - index
		if count == 2:
			sides[index] = "below"
			sides[index + 1] = "above"
		elif count > 2:
			raise ValueError(f"Energy {energy[index]:g} MeV occurs more than twice in {path}")
		index = end
	return sides


def _metadata_density(table: CoefficientTable) -> float | None:
	value = table.metadata.get("density_g_cm3", "none")
	if value.lower() == "none":
		return None
	try:
		density = float(value)
	except ValueError as exc:
		raise ValueError(f"Invalid density metadata in {table.path}") from exc
	if not math.isfinite(density) or density <= 0:
		raise ValueError(f"Invalid density metadata in {table.path}")
	return density


def _validate_linear_relation(table: CoefficientTable) -> None:
	if table.linear_coefficient is None:
		return
	density = _metadata_density(table)
	if density is None:
		raise ValueError(f"Linear coefficient column requires density metadata in {table.path}")
	expected = table.mass_coefficient * density
	if not np.allclose(table.linear_coefficient, expected, rtol=5e-12, atol=0.0):
		raise ValueError(f"Linear coefficient is inconsistent with mass coefficient and density in {table.path}")


def _coefficient_rows_cm_inverse(
	table: CoefficientTable, density_g_cm3: float | None
) -> np.ndarray:
	if table.linear_coefficient is not None:
		if density_g_cm3 is not None:
			requested = float(density_g_cm3)
			stored = _metadata_density(table)
			if stored is None or not math.isclose(requested, stored, rel_tol=1e-12, abs_tol=0.0):
				raise ValueError("Provided density conflicts with the density embedded in the online file")
		return table.linear_coefficient
	if density_g_cm3 is None:
		raise ValueError("A positive density_g_cm3 is required for a mass-only attenuation table")
	density = float(density_g_cm3)
	if not math.isfinite(density) or density <= 0:
		raise ValueError("density_g_cm3 must be finite and positive")
	return table.mass_coefficient * density


def _segments(energy: np.ndarray) -> list[tuple[int, int]]:
	segments: list[tuple[int, int]] = []
	start = 0
	for index in range(1, len(energy)):
		if energy[index] == energy[index - 1]:
			segments.append((start, index - 1))
			start = index
	segments.append((start, len(energy) - 1))
	return segments


def interpolate_linear_attenuation(
	table: CoefficientTable,
	target_energy_mev: np.ndarray | list[float] | float,
	*,
	density_g_cm3: float | None = None,
) -> np.ndarray:
	"""Return linear attenuation using edge-aware log-log interpolation.

	At an exact absorption-edge energy, the above-edge value is selected.
	"""
	target = np.atleast_1d(np.asarray(target_energy_mev, dtype=float))
	if np.any(~np.isfinite(target)) or np.any(target <= 0):
		raise ValueError("Target energies must be finite and positive")
	if np.any(target < table.energy_mev[0]) or np.any(target > table.energy_mev[-1]):
		raise ValueError("Target energy is outside the coefficient table; extrapolation is forbidden")

	values = _coefficient_rows_cm_inverse(table, density_g_cm3)
	result = np.full(target.shape, np.nan, dtype=float)

	for index, value in enumerate(target):
		exact = np.flatnonzero(np.isclose(table.energy_mev, value, rtol=0.0, atol=1e-12))
		if exact.size:
			above = exact[table.edge_side[exact] == "above"]
			result[index] = values[above[-1] if above.size else exact[-1]]
			continue
		for start, end in _segments(table.energy_mev):
			x = table.energy_mev[start : end + 1]
			if x[0] < value < x[-1]:
				y = values[start : end + 1]
				result[index] = math.exp(
					np.interp(math.log(value), np.log(x), np.log(y))
				)
				break
		if not math.isfinite(result[index]):
			raise ValueError(f"No valid interpolation segment for {value:g} MeV")
	return result


def metadata_json(table: CoefficientTable, key: str) -> object:
	"""Decode a JSON-valued metadata field."""
	try:
		return json.loads(table.metadata[key])
	except (KeyError, json.JSONDecodeError) as exc:
		raise ValueError(f"Invalid JSON metadata field {key!r} in {table.path}") from exc
