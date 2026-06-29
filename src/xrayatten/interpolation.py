"""Edge-aware piecewise log-log interpolation."""

from __future__ import annotations

import math
from typing import Sequence

import numpy as np

from .exceptions import InterpolationError
from .models import CoefficientTable, OutputRow


def edge_energies(table: CoefficientTable) -> np.ndarray:
    duplicate = table.energy_mev[1:] == table.energy_mev[:-1]
    return table.energy_mev[1:][duplicate]


def segments(energy: np.ndarray) -> list[tuple[int, int]]:
    pieces: list[tuple[int, int]] = []
    start = 0
    for index in range(1, len(energy)):
        if energy[index] == energy[index - 1]:
            pieces.append((start, index - 1))
            start = index
    pieces.append((start, len(energy) - 1))
    return pieces


def build_output_rows(tables: dict[str, CoefficientTable], energies_kev: object | None) -> list[OutputRow]:
    edges = [edge_energies(table) for table in tables.values()]
    all_edges = np.unique(np.concatenate(edges)) if edges else np.asarray([], dtype=float)
    if energies_kev is None:
        energies = np.unique(np.concatenate([table.energy_mev for table in tables.values()]))
    else:
        values = np.atleast_1d(np.asarray(energies_kev, dtype=float))
        if values.ndim != 1 or values.size == 0:
            raise ValueError("energies_kev must be a non-empty one-dimensional sequence")
        if np.any(~np.isfinite(values)) or np.any(values <= 0):
            raise ValueError("energies_kev must contain only positive finite values")
        energies = np.unique(values / 1000.0)

    rows: list[OutputRow] = []
    for energy in sorted(float(value) for value in energies):
        is_edge = np.any(np.isclose(all_edges, energy, rtol=0.0, atol=1e-12))
        if is_edge:
            rows.append(OutputRow(energy, "below"))
            rows.append(OutputRow(energy, "above"))
        else:
            rows.append(OutputRow(energy, "regular"))
    return rows


def _interpolate_to_kev_grid(
    table: CoefficientTable,
    coefficients: np.ndarray,
    grid_kev: np.ndarray,
    edge_sides: Sequence[str],
) -> np.ndarray:
    nist_kev = table.energy_mev * 1000.0
    grid = np.asarray(grid_kev, dtype=float)
    if grid.ndim != 1 or len(grid) != len(edge_sides):
        raise InterpolationError("grid_kev and edge_sides must be one-dimensional arrays of the same length")
    if np.any(~np.isfinite(grid)) or np.any(grid <= 0):
        raise InterpolationError("interpolation grid energies must be positive finite values")

    idx = np.searchsorted(nist_kev, grid)
    idx_clamped = np.clip(idx, 0, len(nist_kev) - 1)
    is_exact = np.abs(grid - nist_kev[idx_clamped]) < 1e-6

    result = np.empty(len(grid), dtype=float)
    result.fill(np.nan)

    for row_index in np.where(is_exact)[0]:
        matches = np.where(np.abs(nist_kev - grid[row_index]) < 1e-6)[0]
        if len(matches) >= 2:
            result[row_index] = coefficients[matches[-1] if edge_sides[row_index] == "above" else matches[0]]
        elif len(matches) == 1:
            result[row_index] = coefficients[matches[0]]

    non_exact = ~is_exact
    if np.any(non_exact):
        query = grid[non_exact]
        segment = np.searchsorted(nist_kev, query, side="right") - 1
        segment = np.clip(segment, 0, len(nist_kev) - 2)

        e_lo = nist_kev[segment]
        e_hi = nist_kev[segment + 1]
        c_lo = coefficients[segment]
        c_hi = coefficients[segment + 1]
        in_segment = (query > e_lo) & (query < e_hi)

        with np.errstate(divide="ignore", invalid="ignore"):
            log_ratio = np.log(query / e_lo) / np.log(e_hi / e_lo)
            values = np.exp(np.log(c_lo) + log_ratio * np.log(c_hi / c_lo))
        values[~in_segment] = np.nan
        result[non_exact] = values

    if np.any(~np.isfinite(result)):
        bad_index = int(np.flatnonzero(~np.isfinite(result))[0])
        raise InterpolationError(
            f"No edge-safe interpolation interval for {table.symbol} at {grid[bad_index]:g} keV"
        )
    return result


def interpolate_element(table: CoefficientTable, coefficients: np.ndarray, rows: list[OutputRow]) -> np.ndarray:
    grid_kev = np.asarray([row.energy_kev for row in rows], dtype=float)
    edge_sides = [row.edge_side for row in rows]
    return _interpolate_to_kev_grid(table, coefficients, grid_kev, edge_sides)


# ── Precision-based dense grid for Origin-ready output ─────────────────

# Base step specs: (lower_keV, upper_keV, step_keV) — "super" precision
_BASE_GRID_SPECS: list[tuple[float, float, float]] = [
    (1.0, 10.0, 0.00001),
    (10.0, 100.0, 0.0001),
    (100.0, 1000.0, 0.001),
    (1000.0, 10000.0, 0.01),
    (10000.0, 20000.0, 0.1),
]

def _collect_all_nist_points(element_tables: dict[str, CoefficientTable]) -> np.ndarray:
    """Collect all unique NIST energy points (MeV) from all element tables."""
    all_e = np.concatenate([table.energy_mev for table in element_tables.values()])
    return np.unique(all_e)


PRECISION_MULTIPLIERS: dict[str, int] = {
    "direct": 0,  # NIST grid points only, no fill
    "super": 1,
    "high": 10,
    "medium": 100,
    "low": 1000,
    "fast": 10000,
}

# Cache: (precision, frozenset_of_symbols) -> (grid, edge_sides)
_GRID_CACHE: dict[tuple, tuple[np.ndarray, list[str]]] = {}


def get_precision_multiplier(level: str) -> int:
    if level not in PRECISION_MULTIPLIERS:
        raise ValueError(f"precision must be one of {list(PRECISION_MULTIPLIERS.keys())}")
    return PRECISION_MULTIPLIERS[level]


def build_compound_dense_grid(
    element_tables: dict[str, CoefficientTable],
    precision: str = "low",
) -> tuple[np.ndarray, list[str]]:
    """Build a dense energy grid for compound coefficient output.

    Grid construction (in order):
    1. ALL NIST data points from all elements (exact, no approximation)
    2. ALL absorption edges from all elements (exact, below/above pairs)
    3. Fixed-step fill points at the specified precision between breakpoints

    These three sets are merged and deduplicated. Every NIST point and
    every absorption edge appears as an exact value in the final grid.

    Returns (energies_kev, edge_sides).
    """
    cache_key = (precision, frozenset(element_tables.keys()))
    if cache_key in _GRID_CACHE:
        return _GRID_CACHE[cache_key]

    # 1. Collect all NIST data points (keV) from all elements
    nist_kev_set: set[float] = set()
    for table in element_tables.values():
        for e in table.energy_mev:
            nist_kev_set.add(round(float(e) * 1000.0, 10))

    # 2. Collect all absorption edges (keV) — these need below/above pairs
    edge_kev_set: set[float] = set()
    for table in element_tables.values():
        for e in edge_energies(table):
            edge_kev_set.add(round(float(e) * 1000.0, 10))

    # 3. Build fixed-step fill points between consecutive breakpoints
    #    Breakpoints = NIST points (which include edges as duplicates)
    breakpoints_kev = sorted(nist_kev_set)
    mult = get_precision_multiplier(precision)

    fill_points: set[float] = set()
    if mult > 0:
        for i in range(len(breakpoints_kev) - 1):
            e_lo = breakpoints_kev[i]
            e_hi = breakpoints_kev[i + 1]
            gap = e_hi - e_lo
            if gap <= 0:
                continue
            base_step = 0.01
            for lo, hi, step in _BASE_GRID_SPECS:
                if lo <= e_lo < hi:
                    base_step = step
                    break
            step = base_step * mult
            if gap > step * 1.5:
                first = math.ceil(e_lo / step) * step
                if first < e_hi:
                    n = int((e_hi - first) / step)
                    for k in range(n):
                        val = round(first + k * step, 10)
                        if e_lo < val < e_hi:
                            fill_points.add(val)

    # 4. Merge all points, with edge energies duplicated (below + above)
    non_edge_points = (nist_kev_set | fill_points) - edge_kev_set
    # Edge energies need two entries each
    edge_list: list[float] = []
    for e in sorted(edge_kev_set):
        edge_list.append(e)
        edge_list.append(e)
    all_points = sorted(non_edge_points) + edge_list
    grid = np.array(sorted(all_points))

    # 5. Build edge_sides: mark below/above at every edge energy
    edge_sides: list[str] = ["regular"] * len(grid)
    for e_edge in sorted(edge_kev_set):
        # Find all grid indices at this exact energy
        idxs = np.where(np.abs(grid - e_edge) < 1e-9)[0]
        if len(idxs) >= 2:
            edge_sides[idxs[0]] = "below"
            edge_sides[idxs[1]] = "above"
        elif len(idxs) == 1:
            edge_sides[idxs[0]] = "below"

    _GRID_CACHE[cache_key] = (grid, edge_sides)
    return grid, edge_sides


def build_compound_range_grid(
    element_tables: dict[str, CoefficientTable],
    precision: str,
    start_kev: float,
    stop_kev: float,
) -> tuple[np.ndarray, list[str]]:
    """Return the precision grid rows inside an inclusive energy range."""
    if not math.isfinite(start_kev) or not math.isfinite(stop_kev) or start_kev <= 0 or stop_kev <= 0:
        raise InterpolationError("energy_range_kev values must be positive finite values")
    if start_kev >= stop_kev:
        raise InterpolationError("energy_range_kev.start must be smaller than energy_range_kev.stop")

    grid, edge_sides = build_compound_dense_grid(element_tables, precision)
    min_grid = float(grid[0])
    max_grid = float(grid[-1])
    tolerance = 1e-9
    if start_kev < min_grid - tolerance or stop_kev > max_grid + tolerance:
        raise InterpolationError(
            f"energy_range_kev must stay within the local NIST range {min_grid:g}-{max_grid:g} keV"
        )

    mask = (grid >= start_kev - tolerance) & (grid <= stop_kev + tolerance)
    if not np.any(mask):
        raise InterpolationError(
            "energy_range_kev contains no points from the selected precision grid; "
            "use a wider range or a higher precision"
        )
    return grid[mask], [side for side, keep in zip(edge_sides, mask) if keep]


def interpolate_element_to_grid(
    table: CoefficientTable,
    coefficients: np.ndarray,
    grid_kev: np.ndarray,
    edge_sides: list[str],
) -> np.ndarray:
    """Interpolate one element to a keV grid using the precision-grid logic."""
    return _interpolate_to_kev_grid(table, coefficients, grid_kev, edge_sides)


# ── Legacy single-element dense interpolation (kept for backward compat) ──

def dense_loglog_interpolation(
    energy_mev: np.ndarray,
    coefficients: np.ndarray,
    *,
    points_per_segment: int = 200,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Legacy single-element dense interpolation."""
    pieces = segments(energy_mev)

    breakpoints: set[float] = set()
    for start, end in pieces:
        breakpoints.add(float(energy_mev[start]))
        breakpoints.add(float(energy_mev[end]))
    for i in range(1, len(energy_mev)):
        if energy_mev[i] == energy_mev[i - 1]:
            breakpoints.add(float(energy_mev[i]))

    sorted_bp = sorted(breakpoints)

    dense_energies: list[float] = []
    dense_values: list[float] = []
    dense_labels: list[str] = []

    for i in range(len(sorted_bp) - 1):
        e_lo, e_hi = sorted_bp[i], sorted_bp[i + 1]

        is_edge = False
        for start, end in pieces:
            x = energy_mev[start : end + 1]
            if len(x) >= 2 and x[0] == x[1] and abs(e_lo - x[0]) < 1e-12:
                is_edge = True
                break

        if is_edge and (not dense_energies or abs(dense_energies[-1] - e_lo * 1000) > 1e-6):
            for idx in range(len(energy_mev)):
                if abs(energy_mev[idx] - e_lo) < 1e-12:
                    dense_energies.append(e_lo * 1000.0)
                    dense_values.append(float(coefficients[idx]))
                    dense_labels.append("below")
                    break
            for idx in range(len(energy_mev) - 1, -1, -1):
                if abs(energy_mev[idx] - e_lo) < 1e-12:
                    dense_energies.append(e_lo * 1000.0)
                    dense_values.append(float(coefficients[idx]))
                    dense_labels.append("above")
                    break
        elif not dense_energies or abs(dense_energies[-1] - e_lo * 1000) > 1e-6:
            dense_energies.append(e_lo * 1000.0)
            val = _interp_at(e_lo, energy_mev, coefficients, pieces)
            dense_values.append(val)
            dense_labels.append("regular")

        if e_hi > e_lo * (1.0 + 1e-12):
            fill = np.logspace(math.log10(e_lo), math.log10(e_hi), points_per_segment + 2)[1:-1]
            for e in fill:
                e_f = float(e)
                if dense_energies and abs(dense_energies[-1] - e_f * 1000) / max(abs(dense_energies[-1]), 1e-30) < 1e-8:
                    continue
                val = _interp_at(e_f, energy_mev, coefficients, pieces)
                dense_energies.append(e_f * 1000.0)
                dense_values.append(val)
                dense_labels.append("regular")

    e_last = sorted_bp[-1]
    if not dense_energies or abs(dense_energies[-1] - e_last * 1000) > 1e-6:
        dense_energies.append(e_last * 1000.0)
        val = _interp_at(e_last, energy_mev, coefficients, pieces)
        dense_values.append(val)
        dense_labels.append("regular")

    order = np.argsort(dense_energies)
    return (
        np.asarray(dense_energies)[order],
        np.asarray(dense_values)[order],
        [dense_labels[i] for i in order],
    )


def _interp_at(energy_mev: float, table_energy: np.ndarray, coefficients: np.ndarray, pieces: list[tuple[int, int]]) -> float:
    """Interpolate coefficient at a single energy using log-log piecewise."""
    exact = np.flatnonzero(np.isclose(table_energy, energy_mev, rtol=0.0, atol=1e-12))
    if exact.size:
        return float(coefficients[exact[0]])
    for start, end in pieces:
        x = table_energy[start : end + 1]
        if x[0] < energy_mev < x[-1]:
            y = coefficients[start : end + 1]
            return float(math.exp(np.interp(math.log(energy_mev), np.log(x), np.log(y))))
    return float("nan")
