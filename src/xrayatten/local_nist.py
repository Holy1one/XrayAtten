"""Local NIST v1.4 snapshot access and coefficient calculations."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Mapping

import numpy as np

from .exceptions import DataIntegrityError
from .interpolation import build_output_rows, interpolate_element
from .models import CoefficientResult, CoefficientTable, ElementRecord, MaterialSpec, OutputRow
from .physics import validate_positive_finite


PACKAGE_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = PACKAGE_DIR.parent.parent
DATA_DIR = PACKAGE_DIR / "data"
NIST_DIR = DATA_DIR / "nist_local"
COEFF_DIR = NIST_DIR / "coeff"
ELEMENTS_CSV = DATA_DIR / "elements.csv"
MANIFEST_PATH = NIST_DIR / "manifest.json"

DATA_SOURCE_LOCAL = "nist_local_v1_4"
NIST_DATABASE = "NIST Tables of X-Ray Mass Attenuation Coefficients and Mass Energy-Absorption Coefficients"
NIST_VERSION = "1.4"
NIST_VERSION_DATE = "2004-07-12"
NIST_URL = "https://physics.nist.gov/PhysRefData/XrayMassCoef/"
ENERGY_MIN_MEV = 0.001
ENERGY_MAX_MEV = 20.0
APPROX_WARNING = "elemental mass-fraction additivity approximation for custom-compound mu_en/rho"
INTERPOLATION_RULE = "piecewise log-log interpolation; absorption-edge intervals separated; no extrapolation"
EDGE_POLICY = "duplicate energy rows labeled below/above; exact absorption-edge hits use above-edge values"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_manifest(path: Path = MANIFEST_PATH) -> dict:
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise DataIntegrityError(f"NIST data manifest is missing: {path}") from exc
    if manifest.get("schema") != "nist_xraymasscoef_snapshot_manifest_v1":
        raise DataIntegrityError(f"Invalid NIST manifest schema: {path}")
    return manifest


def read_elements(path: Path = ELEMENTS_CSV) -> dict[str, ElementRecord]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    required = {"原子序数", "符号", "标准原子质量", "原子质量来源"}
    if not rows or not required.issubset(rows[0].keys()):
        raise DataIntegrityError(f"elements.csv must contain columns: {sorted(required)}")

    records: dict[str, ElementRecord] = {}
    atomic_numbers: list[int] = []
    for row in rows:
        try:
            atomic_number = int(row["原子序数"])
            atomic_weight = float(row["标准原子质量"])
        except (TypeError, ValueError) as exc:
            raise DataIntegrityError("elements.csv contains invalid numeric values") from exc
        symbol = row["符号"].strip()
        if symbol in records:
            raise DataIntegrityError(f"elements.csv contains duplicate symbol: {symbol}")
        if atomic_weight <= 0 or not math.isfinite(atomic_weight):
            raise DataIntegrityError(f"Invalid atomic weight for {symbol}")
        atomic_numbers.append(atomic_number)
        records[symbol] = ElementRecord(
            atomic_number=atomic_number,
            symbol=symbol,
            atomic_weight=atomic_weight,
            atomic_weight_source=row["原子质量来源"].strip(),
        )
    if atomic_numbers != list(range(1, len(atomic_numbers) + 1)):
        raise DataIntegrityError("elements.csv atomic numbers must be continuous and ordered")
    return records


def parse_composition(composition: Mapping[str, object] | str) -> dict[str, float]:
    if isinstance(composition, str):
        if not composition.strip():
            raise ValueError("composition must not be empty")
        items: list[tuple[str, object]] = []
        for part in composition.split(","):
            if part.count(":") != 1:
                raise ValueError("String composition must use 'Element:value' entries")
            symbol, value = (piece.strip() for piece in part.split(":", 1))
            items.append((symbol, value))
    elif isinstance(composition, Mapping):
        items = list(composition.items())
    else:
        raise TypeError("composition must be a mapping or an 'Element:value' string")
    if not items:
        raise ValueError("composition must not be empty")

    parsed: dict[str, float] = {}
    for raw_symbol, raw_value in items:
        symbol = str(raw_symbol).strip()
        if not re.fullmatch(r"[A-Z][a-z]?", symbol):
            raise ValueError(f"Invalid element symbol: {symbol!r}")
        if symbol in parsed:
            raise ValueError(f"Duplicate element in composition: {symbol}")
        try:
            value = float(raw_value)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Invalid composition value for {symbol}") from exc
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f"Composition value for {symbol} must be finite and positive")
        parsed[symbol] = value
    return parsed


def mass_fractions(
    composition: Mapping[str, object] | str,
    basis: str,
    elements: dict[str, ElementRecord] | None = None,
) -> tuple[dict[str, float], dict[str, float]]:
    parsed = parse_composition(composition)
    element_records = elements or read_elements()
    unknown = [symbol for symbol in parsed if symbol not in element_records]
    if unknown:
        raise ValueError(f"Unknown elements: {unknown}")
    unsupported = [symbol for symbol in parsed if element_records[symbol].atomic_number > 92]
    if unsupported:
        raise ValueError(f"NIST v1.4 coefficient tables are limited to Z=1-92: {unsupported}")

    if basis == "atomic_ratio":
        weighted = {
            symbol: amount * element_records[symbol].atomic_weight
            for symbol, amount in parsed.items()
        }
        total = math.fsum(weighted.values())
        return parsed, {symbol: value / total for symbol, value in weighted.items()}
    if basis == "mass_fraction":
        total = math.fsum(parsed.values())
        if not math.isclose(total, 1.0, rel_tol=0.0, abs_tol=1e-8):
            raise ValueError(f"mass_fraction values must sum to 1.0; received {total:.16g}")
        return parsed, {symbol: value / total for symbol, value in parsed.items()}
    raise ValueError("composition_basis must be 'atomic_ratio' or 'mass_fraction'")


def validate_density(density_g_cm3: object | None, *, required: bool) -> float | None:
    if density_g_cm3 is None:
        if required:
            raise ValueError("density_g_cm3 is required and must be positive")
        return None
    return validate_positive_finite(density_g_cm3, "density_g_cm3")


def read_element_table(symbol: str, manifest: dict | None = None) -> CoefficientTable:
    source_manifest = manifest or load_manifest()
    entry = source_manifest.get("elements", {}).get(symbol)
    if entry is None:
        raise DataIntegrityError(f"No manifest entry for element {symbol}")
    path = COEFF_DIR / f"{symbol}.txt"
    actual_hash = sha256_file(path)
    if actual_hash != entry.get("sha256"):
        raise DataIntegrityError(f"NIST source file hash mismatch for {symbol}: {path}")
    try:
        data = np.loadtxt(path, skiprows=1, dtype=float)
    except Exception as exc:
        raise DataIntegrityError(f"Invalid coefficient table: {path}") from exc
    if data.ndim != 2 or data.shape[1] < 3:
        raise DataIntegrityError(f"Invalid coefficient table shape: {path}")
    data = data[:, :3]
    if np.any(~np.isfinite(data)) or np.any(data <= 0):
        raise DataIntegrityError(f"Non-positive or non-finite NIST data in {path}")
    if np.any(np.diff(data[:, 0]) < 0):
        raise DataIntegrityError(f"Unsorted energy data in {path}")
    return CoefficientTable(
        symbol=symbol,
        energy_mev=data[:, 0],
        mass_attenuation=data[:, 1],
        mass_energy_absorption=data[:, 2],
        sha256=actual_hash,
        path=path,
    )


def load_material_tables(material: MaterialSpec) -> tuple[dict[str, float], dict[str, CoefficientTable]]:
    elements = read_elements()
    _, weights = mass_fractions(material.composition, material.composition_basis, elements)
    manifest = load_manifest()
    tables = {symbol: read_element_table(symbol, manifest) for symbol in weights}
    return weights, tables


def source_files_metadata(tables: dict[str, CoefficientTable]) -> dict[str, dict[str, str]]:
    manifest = load_manifest()
    elements = read_elements()
    metadata: dict[str, dict[str, str]] = {}
    for symbol, table in tables.items():
        record = elements[symbol]
        manifest_entry = manifest["elements"][symbol]
        metadata[symbol] = {
            "file": table.path.relative_to(PROJECT_ROOT).as_posix(),
            "sha256": table.sha256,
            "url": manifest_entry.get("url", f"{NIST_URL}ElemTab/z{record.atomic_number:02d}.html"),
        }
    return metadata


def base_metadata(
    material: MaterialSpec,
    mass_fraction_values: dict[str, float],
    tables: dict[str, CoefficientTable],
    coefficient_kind: str,
) -> dict[str, str]:
    manifest = load_manifest()
    metadata = {
        "data_source": DATA_SOURCE_LOCAL,
        "coefficient_kind": coefficient_kind,
        "source_database": manifest.get("database", NIST_DATABASE),
        "source_version": manifest.get("version", NIST_VERSION),
        "source_version_date": manifest.get("version_date", NIST_VERSION_DATE),
        "source_manifest": MANIFEST_PATH.relative_to(PROJECT_ROOT).as_posix(),
        "source_manifest_sha256": sha256_file(MANIFEST_PATH),
        "source_snapshot_accessed": str(manifest.get("snapshot_accessed", "")),
        "source_version_verified": str(manifest.get("version_verified", "")),
        "elements_source_file": ELEMENTS_CSV.relative_to(PROJECT_ROOT).as_posix(),
        "elements_source_sha256": sha256_file(ELEMENTS_CSV),
        "composition_basis": material.composition_basis,
        "composition_json": json.dumps(material.composition, sort_keys=True, separators=(",", ":")),
        "mass_fractions_json": json.dumps(mass_fraction_values, sort_keys=True, separators=(",", ":")),
        "density_g_cm3": "none" if material.density_g_cm3 is None else f"{material.density_g_cm3:.16g}",
        "density_source": material.density_source or "not specified by user",
        "interpolation": INTERPOLATION_RULE,
        "edge_policy": EDGE_POLICY,
        "source_files_json": json.dumps(source_files_metadata(tables), sort_keys=True, separators=(",", ":")),
        "scope": "monoenergetic narrow primary-beam attenuation using local NIST v1.4 snapshot",
    }
    if "F" in tables:
        metadata["local_data_correction"] = (
            "Fluorine 60 keV energy label corrected from an erroneous local duplicate 50 keV label; "
            "coefficient values were unchanged and verified against NIST on 2026-06-11"
        )
    if coefficient_kind == "energy_absorption_approx":
        metadata["approximation"] = APPROX_WARNING
        metadata["warning"] = APPROX_WARNING
    return metadata


def compute_coefficients(
    material: MaterialSpec,
    *,
    coefficient_kind: str = "attenuation",
    energies_kev: object | None = None,
) -> CoefficientResult:
    if coefficient_kind not in {"attenuation", "energy_absorption_approx"}:
        raise ValueError("coefficient_kind must be 'attenuation' or 'energy_absorption_approx'")
    density = validate_density(material.density_g_cm3, required=False)
    normalized_material = MaterialSpec(
        id=material.id,
        composition_basis=material.composition_basis,
        composition=parse_composition(material.composition),
        density_g_cm3=density,
        density_source=material.density_source or "not specified by user",
    )
    weights, tables = load_material_tables(normalized_material)
    rows = build_output_rows(tables, energies_kev)
    mass_values = np.zeros(len(rows), dtype=float)
    for symbol, weight in weights.items():
        table = tables[symbol]
        source = table.mass_attenuation if coefficient_kind == "attenuation" else table.mass_energy_absorption
        mass_values += weight * interpolate_element(table, source, rows)
    linear_values = mass_values * density if density is not None else None
    warnings: list[str] = []
    limitations: list[str] = []
    if coefficient_kind == "energy_absorption_approx":
        warnings.append(APPROX_WARNING)
        limitations.extend(
            [
                "No scattered-photon buildup is included.",
                "No fluorescence escape/reabsorption or bremsstrahlung transport is included.",
                "No coupled photon-electron Monte Carlo transport is included.",
            ]
        )
    metadata = base_metadata(normalized_material, weights, tables, coefficient_kind)
    return CoefficientResult(
        material=normalized_material,
        coefficient_kind=coefficient_kind,
        rows=rows,
        mass_coefficient=mass_values,
        linear_coefficient=linear_values,
        mass_fractions=weights,
        metadata=metadata,
        warnings=warnings,
        limitations=limitations,
        element_tables=tables,
    )


def coefficient_at_energy_kev(
    material: MaterialSpec,
    coefficient_kind: str,
    energy_kev: float,
) -> tuple[float, float | None, dict[str, str], dict[str, float]]:
    result = compute_coefficients(material, coefficient_kind=coefficient_kind, energies_kev=[energy_kev])
    above = [index for index, row in enumerate(result.rows) if row.edge_side == "above"]
    index = above[-1] if above else len(result.rows) - 1
    linear = None if result.linear_coefficient is None else float(result.linear_coefficient[index])
    return float(result.mass_coefficient[index]), linear, result.metadata, result.mass_fractions


def validate_data() -> dict[str, object]:
    manifest = load_manifest()
    elements = read_elements()
    table_entries = manifest.get("elements", {})
    if len(table_entries) != 92:
        raise DataIntegrityError(f"Expected 92 NIST element tables, found {len(table_entries)}")
    expected_symbols = [record.symbol for record in elements.values() if record.atomic_number <= 92]
    missing = [symbol for symbol in expected_symbols if symbol not in table_entries]
    if missing:
        raise DataIntegrityError(f"Manifest is missing NIST elements: {missing}")
    for symbol in expected_symbols:
        table = read_element_table(symbol, manifest)
        if table.energy_mev[0] < ENERGY_MIN_MEV or table.energy_mev[-1] > ENERGY_MAX_MEV:
            raise DataIntegrityError(f"Unexpected energy range for {symbol}")
    corrections = manifest.get("local_corrections", [])
    if not any(item.get("symbol") == "F" and item.get("corrected_value") == 0.06 for item in corrections):
        raise DataIntegrityError("Manifest must retain the fluorine 60 keV correction note")
    return {
        "schema": manifest["schema"],
        "database": manifest.get("database"),
        "version": manifest.get("version"),
        "elements": len(expected_symbols),
        "manifest_sha256": sha256_file(MANIFEST_PATH),
    }
