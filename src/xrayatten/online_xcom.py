"""Optional online NIST XCOM v1.5 attenuation comparison."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import re
import time
from typing import Mapping, Sequence

import numpy as np
import pandas as pd

from .exceptions import OnlineXcomError
from .models import MaterialSpec
from .output import write_dataframe
from .plotting import plt


DATA_SOURCE_ONLINE = "nist_xcom_v1_5"
NIST_DATABASE = "NIST XCOM Photon Cross Sections Database (SRD 8)"
NIST_VERSION = "1.5"
NIST_VERSION_DATE = "2010-11"
NIST_LANDING_URL = "https://physics.nist.gov/PhysRefData/Xcom/html/xcom1.html"
NIST_REQUEST_ENDPOINT = "https://physics.nist.gov/cgi-bin/Xcom/xcom3_2"
ENERGY_MIN_MEV = 0.001
ENERGY_MAX_MEV = 20.0
DEFAULT_TIMEOUT_SECONDS = 60.0
DEFAULT_RETRIES = 3


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


def _fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return format(value.numerator / value.denominator, ".16g")


def _lcm(left: int, right: int) -> int:
    return abs(left * right) // math.gcd(left, right)


def _parse_atomic_ratio(composition: Mapping[str, object]) -> dict[str, Fraction]:
    parsed: dict[str, Fraction] = {}
    for raw_symbol, raw_value in composition.items():
        symbol = str(raw_symbol).strip()
        if not re.fullmatch(r"[A-Z][a-z]?", symbol):
            raise ValueError(f"Invalid element symbol: {symbol!r}")
        if symbol in parsed:
            raise ValueError(f"Duplicate element in composition: {symbol}")
        try:
            value = Fraction(str(raw_value))
        except (ValueError, ZeroDivisionError) as exc:
            raise ValueError(f"Invalid composition value for {symbol}") from exc
        if value <= 0:
            raise ValueError(f"Composition value for {symbol} must be positive")
        parsed[symbol] = value
    if not parsed:
        raise ValueError("composition must not be empty")
    return parsed


def composition_to_xcom_formula(composition: Mapping[str, object]) -> tuple[str, dict[str, Fraction]]:
    parsed = _parse_atomic_ratio(composition)
    denominator = 1
    for amount in parsed.values():
        denominator = _lcm(denominator, amount.denominator)
    integers = {
        symbol: amount.numerator * (denominator // amount.denominator)
        for symbol, amount in parsed.items()
    }
    return "".join(f"{symbol}{amount}" for symbol, amount in integers.items()), parsed


def _validate_energy_grid(energies_kev: Sequence[float] | None) -> tuple[float, ...]:
    if energies_kev is None:
        values = np.asarray([value / 1000.0 for value in range(1, 81)], dtype=float)
    else:
        values = np.asarray(energies_kev, dtype=float) / 1000.0
    if values.ndim != 1 or values.size == 0:
        raise ValueError("online comparison energies must be a non-empty one-dimensional sequence")
    if np.any(~np.isfinite(values)) or np.any(values < ENERGY_MIN_MEV) or np.any(values > ENERGY_MAX_MEV):
        raise ValueError("online comparison energies must be within 1-20000 keV")
    if values.size > 100:
        raise ValueError("NIST XCOM accepts at most 100 additional energies")
    return tuple(float(value) for value in np.unique(values))


def _payload_energies(energies_mev: Sequence[float]) -> str:
    return "\n".join(format(value, ".4g") for value in energies_mev)


def _request_xcom(
    formula: str,
    title: str,
    energies_mev: Sequence[float],
    *,
    session=None,
    timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
    retries: int = DEFAULT_RETRIES,
) -> XcomFetch:
    import requests

    payload = {
        "Formula": formula,
        "Name": title,
        "Graph0": "on",
        "Energies": _payload_energies(energies_mev),
        "NumAdd": "1",
        "WindowXmin": format(ENERGY_MIN_MEV, ".4g"),
        "WindowXmax": format(ENERGY_MAX_MEV, ".4g"),
        "ResizeFlag": "on",
    }
    client = session or requests.Session()
    close_client = session is None
    client.headers.setdefault("User-Agent", "xray-attenuation-online-comparison/0.1")
    last_error: Exception | None = None
    try:
        for attempt in range(retries):
            try:
                response = client.post(NIST_REQUEST_ENDPOINT, data=payload, timeout=timeout_seconds)
                if response.status_code == 200:
                    accessed = datetime.now(timezone.utc).isoformat(timespec="seconds")
                    return XcomFetch(
                        content=response.content,
                        status_code=response.status_code,
                        accessed_at_utc=accessed,
                        response_sha256=hashlib.sha256(response.content).hexdigest(),
                    )
                last_error = OnlineXcomError(f"NIST XCOM returned HTTP {response.status_code}")
            except requests.RequestException as exc:
                last_error = exc
            if attempt + 1 < retries:
                time.sleep(2**attempt)
    finally:
        if close_client:
            client.close()
    raise OnlineXcomError(f"NIST XCOM request failed after {retries} attempts: {last_error}")


def _parse_float(text: str, field: str) -> float:
    try:
        value = float(text.strip())
    except ValueError as exc:
        raise OnlineXcomError(f"Invalid numeric {field} in NIST XCOM response: {text!r}") from exc
    if not math.isfinite(value):
        raise OnlineXcomError(f"Non-finite {field} in NIST XCOM response")
    return value


def parse_xcom_html(html: bytes | str) -> XcomTable:
    from bs4 import BeautifulSoup

    soup = BeautifulSoup(html, "html.parser")
    title = soup.title.get_text(" ", strip=True) if soup.title else ""
    if "XCOM: Error" in title:
        heading = soup.find("h2")
        detail = heading.get_text(" ", strip=True) if heading else title
        raise OnlineXcomError(f"NIST XCOM rejected the request: {detail}")

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
        if ENERGY_MIN_MEV <= energy <= ENERGY_MAX_MEV:
            edge_label = " ".join(cells[0].get_text(" ", strip=True).split())
            raw_rows.append((energy, edge_label, mass_attenuation))

    if not raw_rows:
        raise OnlineXcomError("NIST XCOM response did not contain a coefficient table")
    if not constituents:
        raise OnlineXcomError("NIST XCOM response did not contain constituent weight fractions")

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
        else:
            raise OnlineXcomError(f"Ambiguous duplicate XCOM rows at {group[0][0]:.12g} MeV")
        for (energy, label, mass_attenuation), side in selected:
            if energy <= 0 or mass_attenuation <= 0:
                raise OnlineXcomError("NIST XCOM returned a non-positive coefficient")
            energy_values.append(energy)
            edge_sides.append(side)
            edge_labels.append(label)
            mass_values.append(mass_attenuation)
        index = end
    energy_array = np.asarray(energy_values, dtype=float)
    if np.any(np.diff(energy_array) < 0):
        raise OnlineXcomError("NIST XCOM returned unsorted energies")
    return XcomTable(
        energy_mev=energy_array,
        edge_side=np.asarray(edge_sides, dtype=str),
        edge_label=np.asarray(edge_labels, dtype=str),
        mass_attenuation_cm2_g=np.asarray(mass_values, dtype=float),
        constituents_by_atomic_number=constituents,
    )


def _validate_requested_energies(table: XcomTable, energies_mev: Sequence[float]) -> None:
    for energy in energies_mev:
        if not np.any(np.isclose(table.energy_mev, energy, rtol=0.0, atol=5e-8)):
            raise OnlineXcomError(f"NIST XCOM response omitted requested energy {energy:.12g} MeV")


def _write_plot(path: Path, table: XcomTable, linear_values: np.ndarray | None) -> Path:
    values = linear_values if linear_values is not None else table.mass_attenuation_cm2_g
    ylabel = "Linear attenuation coefficient (1/cm)" if linear_values is not None else "Mass attenuation coefficient (cm2/g)"
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.loglog(table.energy_mev * 1000.0, values, linewidth=1.2)
    ax.set_xlabel("Photon energy (keV)")
    ax.set_ylabel(ylabel)
    ax.grid(which="both", linestyle="--", linewidth=0.5, alpha=0.6)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return path


def compute_online_attenuation(
    material: MaterialSpec,
    *,
    energies_kev: Sequence[float] | None,
    output_dir: Path,
    session=None,
) -> list[Path]:
    if material.composition_basis != "atomic_ratio":
        raise ValueError("Online XCOM comparison supports atomic_ratio only")
    formula, parsed = composition_to_xcom_formula(material.composition)
    energies_mev = _validate_energy_grid(energies_kev)
    fetch = _request_xcom(formula, material.id, energies_mev, session=session)
    table = parse_xcom_html(fetch.content)
    _validate_requested_energies(table, energies_mev)
    linear = None
    if material.density_g_cm3 is not None:
        linear = table.mass_attenuation_cm2_g * material.density_g_cm3
    frame = pd.DataFrame(
        {
            "Energy_keV": table.energy_mev * 1000.0,
            "Edge_side": table.edge_side,
            "Mass_mu_over_rho_cm2_g": table.mass_attenuation_cm2_g,
        }
    )
    if linear is not None:
        frame["Linear_mu_cm_inverse"] = linear
    metadata = {
        "data_source": DATA_SOURCE_ONLINE,
        "coefficient_kind": "attenuation",
        "source_database": NIST_DATABASE,
        "source_version": NIST_VERSION,
        "source_version_date": NIST_VERSION_DATE,
        "source_url": NIST_LANDING_URL,
        "request_endpoint": NIST_REQUEST_ENDPOINT,
        "source_accessed_at_utc": fetch.accessed_at_utc,
        "response_sha256": fetch.response_sha256,
        "requested_composition_json": json.dumps({k: _fraction_text(v) for k, v in parsed.items()}, separators=(",", ":")),
        "composition_basis": "atomic_ratio",
        "xcom_formula": formula,
        "density_g_cm3": "none" if material.density_g_cm3 is None else f"{material.density_g_cm3:.16g}",
        "density_source": material.density_source,
        "comparison_scope": "online NIST XCOM v1.5 attenuation comparison only; not an official workflow input",
    }
    txt = output_dir / f"{material.id}_attenuation_online.txt"
    png = output_dir / f"{material.id}_attenuation_online.png"
    return [write_dataframe(txt, frame, metadata), _write_plot(png, table, linear)]


def run_online_comparison(
    materials: Sequence[MaterialSpec],
    *,
    energies_kev: Sequence[float] | None,
    output_dir: Path,
    session=None,
) -> list[Path]:
    paths: list[Path] = []
    for material in materials:
        paths.extend(
            compute_online_attenuation(
                material,
                energies_kev=energies_kev,
                output_dir=output_dir,
                session=session,
            )
        )
    return paths
