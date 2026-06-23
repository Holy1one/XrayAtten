"""Matplotlib plotting helpers for workflow PNG outputs.

All plots read directly from TXT files — no additional processing.
"""

from __future__ import annotations

import os
from pathlib import Path
import tempfile

import matplotlib

if "MPLCONFIGDIR" not in os.environ:
    _MPLCONFIGDIR = Path(tempfile.gettempdir()) / "xrayatten-matplotlib"
    _MPLCONFIGDIR.mkdir(parents=True, exist_ok=True)
    os.environ["MPLCONFIGDIR"] = str(_MPLCONFIGDIR)
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_coefficient_from_txt(txt_path: Path, png_path: Path, max_points: int = 5000) -> Path:
    """Plot coefficient curve directly from a TXT file.

    For large TXT files (dense interpolation), downsample for plotting
    while the TXT retains full resolution for Origin.
    """
    frame = pd.read_csv(txt_path, comment="#", sep="\t")
    energy = frame["Energy_keV"].values

    # Pick the value column (3rd column)
    value_col = frame.columns[2]
    values = frame[value_col].values

    if "mu_en" in value_col:
        ylabel = "Mass energy-absorption coefficient, approximate (cm2/g)"
    else:
        ylabel = "Mass attenuation coefficient (cm2/g)"

    # Check if linear coefficients are available (4th column)
    if len(frame.columns) > 3:
        values = frame[frame.columns[3]].values
        if "mu_en" in frame.columns[3]:
            ylabel = "Linear energy-absorption coefficient, approximate (1/cm)"
        else:
            ylabel = "Linear attenuation coefficient (1/cm)"

    # Downsample for plotting if too many points
    if len(energy) > max_points:
        step = max(1, len(energy) // max_points)
        mask = np.zeros(len(energy), dtype=bool)
        mask[::step] = True
        # Always include edge rows
        if "Edge_side" in frame.columns:
            mask[frame["Edge_side"].values != "regular"] = True
        energy = energy[mask]
        values = values[mask]

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.loglog(energy, values, linewidth=0.8)
    ax.set_xlabel("Photon energy (keV)")
    ax.set_ylabel(ylabel)
    ax.grid(which="both", linestyle="--", linewidth=0.5, alpha=0.6)
    fig.tight_layout()
    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=300)
    plt.close(fig)
    return png_path


def plot_thickness(path: Path, frames: dict[float, pd.DataFrame], y_column: str, ylabel: str) -> Path:
    fig, ax = plt.subplots(figsize=(8, 5.5))
    max_thickness = 0.0
    for energy_kev, frame in frames.items():
        ax.plot(frame["Thickness_cm"], frame[y_column], label=f"{energy_kev:g} keV")
        max_thickness = max(max_thickness, float(frame["Thickness_cm"].max()))
    ax.set_xlabel("Thickness (cm)")
    ax.set_ylabel(ylabel)
    ax.set_xlim(0, max_thickness)
    ax.set_ylim(bottom=0)
    ax.grid(True, linestyle="--", alpha=0.6)
    ax.legend()
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return path


def plot_multilayer(path: Path, frames: dict[float, pd.DataFrame], layers: list) -> Path:
    total_thickness = sum(layer.thickness_cm for layer in layers)
    boundaries = [0.0]
    cumulative = 0.0
    for layer in layers:
        cumulative += layer.thickness_cm
        boundaries.append(total_thickness - cumulative)
    fig, ax = plt.subplots(figsize=(7.5, 6.0))
    for energy_kev, frame in frames.items():
        ax.plot(
            frame["Attenuation_fraction"] * 100.0,
            frame["Stack_coordinate_from_bottom_cm"],
            label=f"{energy_kev:g} keV",
        )
    for y in boundaries:
        ax.axhline(y, color="0.45", linestyle="--", linewidth=0.8, alpha=0.7)
    upper = total_thickness
    x_text = 1.0
    for layer in layers:
        lower = upper - layer.thickness_cm
        ax.text(x_text, (upper + lower) / 2.0, layer.id, va="center", ha="left", fontsize=8)
        upper = lower
    ax.set_xlabel("Primary-beam attenuation (%)")
    ax.set_ylabel("Stack coordinate from bottom (cm)")
    ax.set_ylim(0, total_thickness)
    ax.set_title(
        "Primary-beam attenuation profiles through multilayer stack\n"
        f"Total thickness = {total_thickness:.3f} cm"
    )
    ax.grid(True, linestyle=":", alpha=0.4)
    ax.legend(fontsize=8)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return path
