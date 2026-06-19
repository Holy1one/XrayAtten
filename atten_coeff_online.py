"""Fetch traceable compound attenuation coefficients from NIST XCOM v1.5.

The script submits compound formulae to the public NIST XCOM web form and
stores total attenuation coefficients including coherent scattering. It does
not use the local ``coeff/`` snapshot and does not silently fall back to it.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
import csv
import hashlib
import json
import math
import re
import time
from typing import Mapping, Sequence

from bs4 import BeautifulSoup
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import requests

from support.attenuation_io_online import SCHEMA_NAME, SCHEMA_VERSION


BASE_DIR = Path(__file__).resolve().parent
ELEMENTS_CSV = BASE_DIR / "data" / "elements.csv"
OUTPUT_DIR = BASE_DIR / "results" / "coefficients" / "online"

NIST_DATABASE = "NIST XCOM Photon Cross Sections Database (SRD 8)"
NIST_VERSION = "1.5"
NIST_VERSION_DATE = "2010-11"
NIST_LANDING_URL = "https://physics.nist.gov/PhysRefData/Xcom/html/xcom1.html"
NIST_REQUEST_ENDPOINT = "https://physics.nist.gov/cgi-bin/Xcom/xcom3_2"
NIST_VERSION_URL = "https://physics.nist.gov/PhysRefData/Xcom/Text/version.shtml"
NIST_HELP_URL = "https://physics.nist.gov/PhysRefData/Xcom/Text/chap4.html"
NIST_DISCLAIMER_URL = "https://www.nist.gov/pml/database-disclaimer"
VERSION_VERIFIED_DATE = "2026-06-12"

ENERGY_MIN_MEV = 0.001
ENERGY_MAX_MEV = 20.0
CUSTOM_ENERGIES_MEV = tuple(value / 1000.0 for value in range(1, 81))
MAX_ADDITIONAL_ENERGIES = 100
DEFAULT_TIMEOUT_SECONDS = 60.0
DEFAULT_RETRIES = 3
DEFAULT_REQUEST_INTERVAL_SECONDS = 1.0


@dataclass(frozen=True)
class MaterialConfig:
    key: str
    composition: Mapping[str, int | float | str]
    include_linear: bool
    density_g_cm3: float | None
    density_source: str
    output_name: str
    xcom_title: str


@dataclass(frozen=True)
class XcomFetch:
    content: bytes
    status_code: int
    accessed_at_utc: str
    response_sha256: str


@dataclass(frozen=True)
class XcomTable:
    energy_mev: np.ndarray
    edge_side: np.ndarray
    edge_label: np.ndarray
    mass_attenuation_cm2_g: np.ndarray
    constituents_by_atomic_number: dict[str, float]


PROJECT_DENSITY_SOURCE = (
    "Project material configuration supplied by the user; not provided or inferred by NIST XCOM"
)

MATERIALS: tuple[MaterialConfig, ...] = (
    MaterialConfig(
        key="green_glass",
        composition={"Si": 50, "Ca": 10, "F": 30, "Na": 10, "Al": 20, "Mn": 1.1, "O": 131.1},
        include_linear=True,
        density_g_cm3=2.2,
        density_source=PROJECT_DENSITY_SOURCE,
        output_name="Si50Ca10F30Na10Al20Mn11O1311_attenuation_online",
        xcom_title="GreenGlass",
    ),
    MaterialConfig(
        key="organic_green",
        composition={"C": 420, "H": 400, "Br": 42, "P": 20, "Mn": 8, "Sb": 2},
        include_linear=True,
        density_g_cm3=1.3,
        density_source=PROJECT_DENSITY_SOURCE,
        output_name="C420H400Br42P20Mn8Sb2_attenuation_online",
        xcom_title="OrganicGreen",
    ),
    MaterialConfig(
        key="blue_glass",
        composition={"Si": 50, "B": 20, "Li": 20, "Mg": 10, "Y": 2, "Cs": 10, "Ce": 1, "O": 155, "F": 9},
        include_linear=True,
        density_g_cm3=2.5,
        density_source=PROJECT_DENSITY_SOURCE,
        output_name="Si50B20Li20Mg10Y2Cs10Ce1O155F9_attenuation_online",
        xcom_title="BlueGlass",
    ),
)


class XcomError(RuntimeError):
    """Raised when NIST XCOM cannot return a valid coefficient table."""


def _element_records(path: Path = ELEMENTS_CSV) -> tuple[set[str], dict[int, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or "原子序数" not in rows[0] or "符号" not in rows[0]:
        raise ValueError(f"Invalid element table: {path}")
    symbols = {row["符号"].strip() for row in rows}
    by_atomic_number = {int(row["原子序数"]): row["符号"].strip() for row in rows}
    return symbols, by_atomic_number


def _parse_composition(composition: Mapping[str, object] | str) -> dict[str, Fraction]:
    if isinstance(composition, str):
        if not composition.strip():
            raise ValueError("composition must not be empty")
        items: list[tuple[str, object]] = []
        seen: set[str] = set()
        for part in composition.split(","):
            if part.count(":") != 1:
                raise ValueError("String composition must use 'Element:value' entries")
            symbol, value = (piece.strip() for piece in part.split(":", 1))
            if symbol in seen:
                raise ValueError(f"Duplicate element in composition: {symbol}")
            seen.add(symbol)
            items.append((symbol, value))
    elif isinstance(composition, Mapping):
        items = list(composition.items())
    else:
        raise TypeError("composition must be a mapping or an 'Element:value' string")
    if not items:
        raise ValueError("composition must not be empty")

    known_symbols, _ = _element_records()
    parsed: dict[str, Fraction] = {}
    for raw_symbol, raw_value in items:
        symbol = str(raw_symbol).strip()
        if not re.fullmatch(r"[A-Z][a-z]?", symbol) or symbol not in known_symbols:
            raise ValueError(f"Unknown or invalid element symbol: {symbol!r}")
        if symbol in parsed:
            raise ValueError(f"Duplicate element in composition: {symbol}")
        try:
            value = Fraction(str(raw_value))
        except (ValueError, ZeroDivisionError) as exc:
            raise ValueError(f"Invalid composition value for {symbol}") from exc
        if value <= 0:
            raise ValueError(f"Composition value for {symbol} must be positive")
        parsed[symbol] = value
    return parsed


def _lcm(left: int, right: int) -> int:
    return abs(left * right) // math.gcd(left, right)


def composition_to_xcom_formula(composition: Mapping[str, object] | str) -> tuple[str, dict[str, Fraction]]:
    """Convert an atomic-ratio composition to an integer XCOM formula."""
    parsed = _parse_composition(composition)
    denominator = 1
    for amount in parsed.values():
        denominator = _lcm(denominator, amount.denominator)
    integer_amounts = {
        symbol: amount.numerator * (denominator // amount.denominator)
        for symbol, amount in parsed.items()
    }
    formula = "".join(f"{symbol}{amount}" for symbol, amount in integer_amounts.items())
    return formula, parsed


def _fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return format(value.numerator / value.denominator, ".16g")


def _validate_energy_grid(energies_mev: Sequence[float]) -> tuple[float, ...]:
    values = np.asarray(tuple(energies_mev), dtype=float)
    if values.ndim != 1 or values.size == 0:
        raise ValueError("custom_energies_mev must be a non-empty one-dimensional sequence")
    if values.size > MAX_ADDITIONAL_ENERGIES:
        raise ValueError(f"NIST XCOM accepts at most {MAX_ADDITIONAL_ENERGIES} additional energies")
    if np.any(~np.isfinite(values)):
        raise ValueError("custom_energies_mev contains non-finite values")
    if np.any((values < ENERGY_MIN_MEV) | (values > ENERGY_MAX_MEV)):
        raise ValueError("custom energies must be within 0.001-20 MeV")
    unique = np.unique(values)
    return tuple(float(value) for value in unique)


def _validate_linear_configuration(
    include_linear: bool,
    density_g_cm3: float | None,
    density_source: str | None,
) -> tuple[float | None, str]:
    if not isinstance(include_linear, bool):
        raise TypeError("include_linear must be a boolean")
    if not include_linear:
        if density_g_cm3 is not None:
            raise ValueError("density_g_cm3 must be omitted when include_linear is False")
        return None, "none"
    if density_g_cm3 is None:
        raise ValueError("A material density is required when include_linear is True")
    density = float(density_g_cm3)
    if not math.isfinite(density) or density <= 0:
        raise ValueError("density_g_cm3 must be finite and positive")
    source = str(density_source or "").strip()
    if not source:
        raise ValueError("density_source is required when include_linear is True")
    return density, source


def _energy_payload(energies_mev: Sequence[float]) -> str:
    return "\n".join(format(value, ".4g") for value in energies_mev)


def _request_xcom(
    formula: str,
    title: str,
    custom_energies_mev: Sequence[float],
    *,
    include_standard_grid: bool,
    session: requests.Session | None = None,
    timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
    retries: int = DEFAULT_RETRIES,
) -> XcomFetch:
    if not re.fullmatch(r"[A-Za-z0-9 ]+", title):
        raise ValueError("XCOM title may contain only ASCII letters, digits, and spaces")
    timeout = float(timeout_seconds)
    if not math.isfinite(timeout) or timeout <= 0:
        raise ValueError("timeout_seconds must be finite and positive")
    if retries < 1:
        raise ValueError("retries must be at least 1")

    payload = {
        "Formula": formula,
        "Name": title,
        "Graph0": "on",
        "Energies": _energy_payload(custom_energies_mev),
        "NumAdd": "1",
        "WindowXmin": format(ENERGY_MIN_MEV, ".4g"),
        "WindowXmax": format(ENERGY_MAX_MEV, ".4g"),
        "ResizeFlag": "on",
    }
    if include_standard_grid:
        payload["Output"] = "on"

    client = session or requests.Session()
    close_client = session is None
    client.headers.setdefault(
        "User-Agent",
        "NIST-XCOM-attenuation-online/1.0 (scientific data retrieval; no authentication)",
    )
    last_error: Exception | None = None
    try:
        for attempt in range(retries):
            try:
                response = client.post(NIST_REQUEST_ENDPOINT, data=payload, timeout=timeout)
                if response.status_code == 200:
                    accessed = datetime.now(timezone.utc).isoformat(timespec="seconds")
                    return XcomFetch(
                        content=response.content,
                        status_code=response.status_code,
                        accessed_at_utc=accessed,
                        response_sha256=hashlib.sha256(response.content).hexdigest(),
                    )
                last_error = XcomError(f"NIST XCOM returned HTTP {response.status_code}")
                if response.status_code < 500 and response.status_code != 429:
                    break
            except requests.RequestException as exc:
                last_error = exc
            if attempt + 1 < retries:
                time.sleep(2**attempt)
    finally:
        if close_client:
            client.close()
    raise XcomError(f"NIST XCOM request failed after {retries} attempts: {last_error}")


def _parse_float(text: str, field: str) -> float:
    try:
        value = float(text.strip())
    except ValueError as exc:
        raise XcomError(f"Invalid numeric {field} in NIST XCOM response: {text!r}") from exc
    if not math.isfinite(value):
        raise XcomError(f"Non-finite {field} in NIST XCOM response")
    return value


def parse_xcom_html(
    html: bytes | str,
    *,
    energy_min_mev: float = ENERGY_MIN_MEV,
    energy_max_mev: float = ENERGY_MAX_MEV,
) -> XcomTable:
    """Parse XCOM HTML and retain only total attenuation with coherent scattering."""
    soup = BeautifulSoup(html, "html.parser")
    title = soup.title.get_text(" ", strip=True) if soup.title else ""
    if "XCOM: Error" in title:
        heading = soup.find("h2")
        detail = heading.get_text(" ", strip=True) if heading else title
        raise XcomError(f"NIST XCOM rejected the request: {detail}")

    constituents: dict[str, float] = {}
    for pre in soup.find_all("pre"):
        for atomic_number, fraction in re.findall(
            r"Z\s*=\s*(\d+)\s*:\s*([0-9.+\-Ee]+)", pre.get_text(" ", strip=True)
        ):
            constituents[str(int(atomic_number))] = _parse_float(fraction, "weight fraction")

    raw_rows: list[tuple[float, str, float]] = []
    for row in soup.select("table tr"):
        cells = row.find_all("td")
        if len(cells) != 9:
            continue
        try:
            energy = float(cells[1].get_text(" ", strip=True))
            mass_attenuation = float(cells[7].get_text(" ", strip=True))
        except ValueError:
            continue
        if not energy_min_mev <= energy <= energy_max_mev:
            continue
        edge_label = " ".join(cells[0].get_text(" ", strip=True).split())
        raw_rows.append((energy, edge_label, mass_attenuation))

    if not raw_rows:
        raise XcomError("NIST XCOM response did not contain a coefficient table")
    if not constituents:
        raise XcomError("NIST XCOM response did not contain constituent weight fractions")

    energy_values: list[float] = []
    edge_sides: list[str] = []
    edge_labels: list[str] = []
    mass_values: list[float] = []
    index = 0
    while index < len(raw_rows):
        end = index + 1
        while end < len(raw_rows) and raw_rows[end][0] == raw_rows[index][0]:
            end += 1
        group = raw_rows[index:end]
        if len(group) == 1:
            selected = [(group[0], "regular")]
        elif len(group) == 2 and any(item[1] for item in group):
            selected = [(group[0], "below"), (group[1], "above")]
        elif all(math.isclose(item[2], group[0][2], rel_tol=0.0, abs_tol=0.0) for item in group):
            selected = [(group[0], "regular")]
        else:
            raise XcomError(
                f"Ambiguous duplicate XCOM rows at {group[0][0]:.12g} MeV"
            )
        for (energy, label, mass_attenuation), side in selected:
            if energy <= 0 or mass_attenuation <= 0:
                raise XcomError("NIST XCOM returned a non-positive coefficient")
            energy_values.append(energy)
            edge_sides.append(side)
            edge_labels.append(label)
            mass_values.append(mass_attenuation)
        index = end

    energy_array = np.asarray(energy_values, dtype=float)
    if np.any(np.diff(energy_array) < 0):
        raise XcomError("NIST XCOM returned unsorted energies")
    return XcomTable(
        energy_mev=energy_array,
        edge_side=np.asarray(edge_sides, dtype=str),
        edge_label=np.asarray(edge_labels, dtype=str),
        mass_attenuation_cm2_g=np.asarray(mass_values, dtype=float),
        constituents_by_atomic_number=constituents,
    )


def _validate_requested_energies(table: XcomTable, requested: Sequence[float]) -> None:
    for energy in requested:
        if not np.any(np.isclose(table.energy_mev, energy, rtol=0.0, atol=5e-8)):
            raise XcomError(f"NIST XCOM response omitted requested energy {energy:.12g} MeV")


def _constituents_metadata(table: XcomTable) -> dict[str, dict[str, object]]:
    _, symbols_by_atomic_number = _element_records()
    return {
        atomic_number: {
            "symbol": symbols_by_atomic_number.get(int(atomic_number), "unknown"),
            "fraction_by_weight": fraction,
        }
        for atomic_number, fraction in table.constituents_by_atomic_number.items()
    }


def _write_output(
    path: Path,
    table: XcomTable,
    linear_values: np.ndarray | None,
    *,
    composition: dict[str, Fraction],
    xcom_formula: str,
    density_g_cm3: float | None,
    density_source: str,
    fetch: XcomFetch,
    custom_energies_mev: Sequence[float],
    include_standard_grid: bool,
) -> None:
    columns = ["energy_MeV", "mass_attenuation_cm2_g"]
    if linear_values is not None:
        columns.append("linear_attenuation_cm_inverse")
    composition_json = {
        symbol: _fraction_text(amount) for symbol, amount in composition.items()
    }
    metadata = [
        ("schema", SCHEMA_NAME),
        ("schema_version", SCHEMA_VERSION),
        ("coefficient_kind", "attenuation"),
        ("columns", " ".join(columns)),
        ("units", "MeV cm2/g" + (" 1/cm" if linear_values is not None else "")),
        ("source_database", NIST_DATABASE),
        ("source_version", NIST_VERSION),
        ("source_version_date", NIST_VERSION_DATE),
        ("source_version_verified", VERSION_VERIFIED_DATE),
        ("source_url", NIST_LANDING_URL),
        ("request_endpoint", NIST_REQUEST_ENDPOINT),
        ("source_accessed_at_utc", fetch.accessed_at_utc),
        ("generated_at_utc", datetime.now(timezone.utc).isoformat(timespec="seconds")),
        ("http_status", str(fetch.status_code)),
        ("response_sha256", fetch.response_sha256),
        ("citation", "Berger MJ et al., XCOM: Photon Cross Sections Database, NIST SRD 8, version 1.5 (2010)"),
        ("requested_composition_json", json.dumps(composition_json, separators=(",", ":"))),
        ("composition_basis", "atomic_ratio"),
        ("xcom_formula", xcom_formula),
        ("xcom_constituents_json", json.dumps(_constituents_metadata(table), separators=(",", ":"))),
        ("density_g_cm3", "none" if density_g_cm3 is None else format(density_g_cm3, ".16g")),
        ("density_source", density_source),
        ("linear_relation", "not_requested" if linear_values is None else "linear_attenuation = mass_attenuation * density"),
        ("energy_range_MeV", f"{ENERGY_MIN_MEV:g} {ENERGY_MAX_MEV:g}"),
        ("include_xcom_standard_grid", str(include_standard_grid).lower()),
        ("custom_energies_MeV_json", json.dumps(list(custom_energies_mev), separators=(",", ":"))),
        ("xcom_energy_precision", "XCOM uses at most 4 significant figures for additional energies"),
        ("edge_policy", "Repeated energy rows are below-edge then above-edge; the reader derives edge_side in memory"),
        ("public_access_statement", "Retrieved through the public NIST HTTPS form without authentication; this is not a license statement"),
        ("nist_version_history_url", NIST_VERSION_URL),
        ("nist_help_url", NIST_HELP_URL),
        ("nist_database_disclaimer_url", NIST_DISCLAIMER_URL),
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for key, value in metadata:
            handle.write(f"# {key}: {value}\n")
        for index, energy in enumerate(table.energy_mev):
            fields = [
                f"{energy:.12e}",
                f"{table.mass_attenuation_cm2_g[index]:.12e}",
            ]
            if linear_values is not None:
                fields.append(f"{linear_values[index]:.12e}")
            handle.write("  ".join(fields) + "\n")


def _write_plot(path: Path, table: XcomTable, values: np.ndarray, ylabel: str) -> None:
    figure, axis = plt.subplots(figsize=(6.5, 4.5))
    axis.loglog(table.energy_mev, values, linewidth=1.2)
    axis.set_xlabel("Photon energy (MeV)")
    axis.set_ylabel(ylabel)
    axis.grid(which="both", linestyle="--", linewidth=0.5, alpha=0.6)
    figure.tight_layout()
    figure.savefig(path, dpi=300)
    plt.close(figure)


def compute_coefficients(
    composition: Mapping[str, object] | str,
    *,
    include_linear: bool | None = None,
    density_g_cm3: float | None = None,
    density_source: str | None = None,
    custom_energies_mev: Sequence[float] = CUSTOM_ENERGIES_MEV,
    include_standard_grid: bool = True,
    output_dir: str | Path = OUTPUT_DIR,
    output_name: str | None = None,
    xcom_title: str = "Material",
    session: requests.Session | None = None,
    timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
    retries: int = DEFAULT_RETRIES,
    composition_basis: str = "atomic_ratio",
    coefficient_kind: str = "attenuation",
    energies_mev: Sequence[float] | None = None,
) -> dict[str, object]:
    """Fetch and write an XCOM compound attenuation table.

    ``composition_basis``, ``coefficient_kind``, and ``energies_mev`` retain a
    narrow compatibility path for callers familiar with the previous API.
    """
    if composition_basis != "atomic_ratio":
        raise ValueError("Online XCOM compound input supports atomic_ratio composition only")
    if coefficient_kind != "attenuation":
        raise ValueError("Online XCOM v1.5 output does not provide mass energy-absorption coefficients")
    if energies_mev is not None:
        custom_energies_mev = energies_mev
    energies = _validate_energy_grid(custom_energies_mev)
    if include_linear is None:
        include_linear = density_g_cm3 is not None
    density, density_source_value = _validate_linear_configuration(
        include_linear, density_g_cm3, density_source
    )
    formula, parsed = composition_to_xcom_formula(composition)

    fetch = _request_xcom(
        formula,
        xcom_title,
        energies,
        include_standard_grid=include_standard_grid,
        session=session,
        timeout_seconds=timeout_seconds,
        retries=retries,
    )
    table = parse_xcom_html(fetch.content)
    _validate_requested_energies(table, energies)
    linear_values = (
        table.mass_attenuation_cm2_g * density if density is not None else None
    )

    directory = Path(output_dir)
    if not directory.is_absolute():
        directory = BASE_DIR / directory
    directory.mkdir(parents=True, exist_ok=True)
    name = output_name or f"{formula}_attenuation_online"
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", name):
        raise ValueError("output_name may contain only letters, digits, dot, underscore, and hyphen")
    txt_path = directory / f"{name}.txt"
    plot_path = directory / f"{name}.png"
    _write_output(
        txt_path,
        table,
        linear_values,
        composition=parsed,
        xcom_formula=formula,
        density_g_cm3=density,
        density_source=density_source_value,
        fetch=fetch,
        custom_energies_mev=energies,
        include_standard_grid=include_standard_grid,
    )
    plot_values = linear_values if linear_values is not None else table.mass_attenuation_cm2_g
    ylabel = (
        "Linear attenuation coefficient (1/cm)"
        if linear_values is not None
        else "Mass attenuation coefficient (cm^2/g)"
    )
    _write_plot(plot_path, table, plot_values, ylabel)
    return {
        "energy_mev": table.energy_mev,
        "edge_side": table.edge_side,
        "edge_label": table.edge_label,
        "mass_coefficient": table.mass_attenuation_cm2_g,
        "linear_coefficient": linear_values,
        "xcom_formula": formula,
        "xcom_constituents": table.constituents_by_atomic_number,
        "txt_path": txt_path,
        "plot_path": plot_path,
        "response_sha256": fetch.response_sha256,
    }


def compute_registered_material(
    material: MaterialConfig,
    *,
    output_dir: str | Path = OUTPUT_DIR,
    session: requests.Session | None = None,
) -> dict[str, object]:
    return compute_coefficients(
        material.composition,
        include_linear=material.include_linear,
        density_g_cm3=material.density_g_cm3,
        density_source=material.density_source,
        custom_energies_mev=CUSTOM_ENERGIES_MEV,
        include_standard_grid=True,
        output_dir=output_dir,
        output_name=material.output_name,
        xcom_title=material.xcom_title,
        session=session,
    )


def main() -> None:
    results: list[dict[str, object]] = []
    with requests.Session() as session:
        for index, material in enumerate(MATERIALS):
            result = compute_registered_material(material, session=session)
            results.append(result)
            print(f"[{material.key}] wrote {result['txt_path']} and {result['plot_path']}")
            if index + 1 < len(MATERIALS):
                time.sleep(DEFAULT_REQUEST_INTERVAL_SECONDS)
    print(f"Completed {len(results)} NIST XCOM material requests.")


if __name__ == "__main__":
    main()
