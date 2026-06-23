"""YAML-dispatched scientific workflows."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from .config import (
    load_yaml,
    validate_coefficients_config,
    validate_common,
    validate_multilayer_config,
    validate_single_material_thickness_config,
)
from .exceptions import ConfigError
from .local_nist import (
    APPROX_WARNING,
    DATA_SOURCE_LOCAL,
    EDGE_POLICY,
    INTERPOLATION_RULE,
    compute_coefficients,
    coefficient_at_energy_kev,
)
from .models import LayerSpec, MaterialSpec, WorkflowResult
from .output import write_coefficient_result, write_dense_coefficient_result, write_dataframe
from .physics import (
    attenuation_fraction,
    first_collision_absorbed_energy_fraction,
    thickness_grid,
    transmission_fraction,
)
from .plotting import plot_coefficient_from_txt, plot_multilayer, plot_thickness


def _energy_name(energy_kev: float) -> str:
    return f"{energy_kev:g}keV"


def run_config(path: str | Path) -> WorkflowResult:
    data = load_yaml(path)
    workflow, _ = validate_common(data)
    if workflow == "coefficients":
        return run_coefficients(validate_coefficients_config(data))
    if workflow == "attenuation_vs_thickness":
        return run_attenuation_vs_thickness(
            validate_single_material_thickness_config(data, "attenuation_vs_thickness")
        )
    if workflow == "energy_absorption_vs_thickness":
        return run_energy_absorption_vs_thickness(
            validate_single_material_thickness_config(data, "energy_absorption_vs_thickness")
        )
    if workflow == "multilayer_attenuation_profile":
        return run_multilayer_attenuation_profile(validate_multilayer_config(data))
    raise ConfigError(f"Unsupported workflow: {workflow}")


def run_coefficients(config: dict) -> WorkflowResult:
    output_dir: Path = config["output_dir"]
    energies_kev = config["energies_kev"]
    precision: str = config.get("precision", "low")
    output_paths: list[Path] = []
    for material in config["materials"]:
        attenuation = compute_coefficients(material, coefficient_kind="attenuation", energies_kev=energies_kev)
        attenuation_dir = output_dir / "attenuation"
        attenuation_txt = attenuation_dir / f"{material.id}_attenuation_local.txt"
        attenuation_png = attenuation_dir / f"{material.id}_attenuation_local.png"
        output_paths.append(write_dense_coefficient_result(attenuation_txt, attenuation, precision=precision))
        output_paths.append(plot_coefficient_from_txt(attenuation_txt, attenuation_png))

        absorption = compute_coefficients(
            material,
            coefficient_kind="energy_absorption_approx",
            energies_kev=energies_kev,
        )
        absorption_dir = output_dir / "energy_absorption"
        absorption_txt = absorption_dir / f"{material.id}_energy_absorption_approx_local.txt"
        absorption_png = absorption_dir / f"{material.id}_energy_absorption_approx_local.png"
        output_paths.append(write_dense_coefficient_result(absorption_txt, absorption, precision=precision))
        output_paths.append(plot_coefficient_from_txt(absorption_txt, absorption_png))

    if config["online_comparison"]:
        from .online_xcom import run_online_comparison

        try:
            online_paths = run_online_comparison(
                config["materials"],
                energies_kev=energies_kev,
                output_dir=output_dir / "online_comparison" / "attenuation",
            )
        except Exception as exc:
            raise RuntimeError(
                "Local coefficient outputs were written successfully, but online "
                "NIST XCOM comparison failed. Re-run with online_comparison: false "
                "to keep only local outputs, or retry when the network/service is available."
            ) from exc
        output_paths.extend(online_paths)
    return WorkflowResult("coefficients", output_dir, output_paths)


def run_attenuation_vs_thickness(config: dict) -> WorkflowResult:
    output_dir: Path = config["output_dir"]
    material: MaterialSpec = config["material"]
    thickness = thickness_grid(config["maximum_cm"], config["step_cm"])
    frames: dict[float, pd.DataFrame] = {}
    output_paths: list[Path] = []
    for energy_kev in config["energies_kev"]:
        mass_mu, linear_mu, metadata, _ = coefficient_at_energy_kev(material, "attenuation", energy_kev)
        if linear_mu is None:
            raise ValueError("density_g_cm3 is required for attenuation_vs_thickness")
        transmission = transmission_fraction(linear_mu, thickness)
        attenuation = 1.0 - transmission
        frame = pd.DataFrame(
            {
                "Thickness_cm": thickness,
                "Mass_mu_over_rho_cm2_g": np.full_like(thickness, mass_mu),
                "Linear_mu_cm_inverse": np.full_like(thickness, linear_mu),
                "Transmission_fraction": transmission,
                "Attenuation_fraction": attenuation,
            }
        )
        file_metadata = {
            **metadata,
            "workflow": "attenuation_vs_thickness",
            "energy_keV": f"{energy_kev:.12g}",
            "equation": "T = exp(-mu*x); R = 1 - exp(-mu*x)",
            "scope": "monoenergetic narrow primary-beam attenuation",
        }
        path = output_dir / f"{_energy_name(energy_kev)}_attenuation.txt"
        output_paths.append(write_dataframe(path, frame, file_metadata))
        frames[float(energy_kev)] = frame
    output_paths.append(
        plot_thickness(
            output_dir / "attenuation_vs_thickness.png",
            frames,
            "Attenuation_fraction",
            "Primary-beam attenuation (fraction)",
        )
    )
    return WorkflowResult("attenuation_vs_thickness", output_dir, output_paths)


def run_energy_absorption_vs_thickness(config: dict) -> WorkflowResult:
    output_dir: Path = config["output_dir"]
    material: MaterialSpec = config["material"]
    thickness = thickness_grid(config["maximum_cm"], config["step_cm"])
    frames: dict[float, pd.DataFrame] = {}
    output_paths: list[Path] = []
    for energy_kev in config["energies_kev"]:
        mass_mu, linear_mu, attenuation_metadata, _ = coefficient_at_energy_kev(material, "attenuation", energy_kev)
        mass_mu_en, linear_mu_en, absorption_metadata, _ = coefficient_at_energy_kev(
            material, "energy_absorption_approx", energy_kev
        )
        if linear_mu is None or linear_mu_en is None:
            raise ValueError("density_g_cm3 is required for energy_absorption_vs_thickness")
        if linear_mu_en > linear_mu * (1.0 + 1e-10):
            raise ValueError("mu_en must not exceed total attenuation mu")
        primary_transmission = transmission_fraction(linear_mu, thickness)
        primary_removal = 1.0 - primary_transmission
        absorbed = first_collision_absorbed_energy_fraction(linear_mu, linear_mu_en, thickness)
        removed_not_absorbed = primary_removal - absorbed
        frame = pd.DataFrame(
            {
                "Thickness_cm": thickness,
                "Transmission_fraction": primary_transmission,
                "Attenuation_fraction": primary_removal,
                "First_collision_absorbed_energy_fraction": absorbed,
                "Removed_not_absorbed_fraction": removed_not_absorbed,
            }
        )
        file_metadata = {
            **attenuation_metadata,
            "workflow": "energy_absorption_vs_thickness",
            "energy_keV": f"{energy_kev:.12g}",
            "coefficient_kind": "energy_absorption_approx",
            "mu_en_metadata": json.dumps(absorption_metadata, sort_keys=True, separators=(",", ":")),
            "quantity_name": "First-collision absorbed energy estimate",
            "equation": "A_E = (mu_en / mu) * (1 - exp(-mu*x))",
            "approximation": APPROX_WARNING,
            "limitations": (
                "No scattered-photon buildup; no fluorescence escape/reabsorption; "
                "no coupled photon-electron transport; not detector efficiency."
            ),
        }
        path = output_dir / f"{_energy_name(energy_kev)}_energy_absorption.txt"
        output_paths.append(write_dataframe(path, frame, file_metadata))
        frames[float(energy_kev)] = frame
    output_paths.append(
        plot_thickness(
            output_dir / "energy_absorption_vs_thickness.png",
            frames,
            "First_collision_absorbed_energy_fraction",
            "First-collision absorbed energy estimate (fraction)",
        )
    )
    return WorkflowResult("energy_absorption_vs_thickness", output_dir, output_paths)


def _layer_sample_rows(layers: list[LayerSpec], step_cm: float) -> list[tuple[int, LayerSpec, float, float, str]]:
    rows: list[tuple[int, LayerSpec, float, float, str]] = []
    cumulative = 0.0
    for layer_index, layer in enumerate(layers, start=1):
        layer_start = cumulative
        intervals = max(1, int(math.ceil(layer.thickness_cm / step_cm)))
        positions = np.linspace(0.0, layer.thickness_cm, intervals + 1)
        for point_index, depth_in_layer in enumerate(positions):
            global_depth = layer_start + float(depth_in_layer)
            if point_index == 0 and layer_index == 1:
                label = "stack_top"
            elif point_index == 0:
                label = "layer_entry"
            elif point_index == len(positions) - 1 and layer_index == len(layers):
                label = "stack_bottom"
            elif point_index == len(positions) - 1:
                label = "layer_exit"
            else:
                label = "interior"
            rows.append((layer_index, layer, global_depth, float(depth_in_layer), label))
        cumulative += layer.thickness_cm
    return rows


def _multilayer_metadata(
    energy_kev: float,
    layers: list[LayerSpec],
    layer_mass_fractions: list[dict[str, float]],
    step_cm: float,
    base_metadata: dict[str, str],
) -> dict[str, object]:
    total_thickness = sum(layer.thickness_cm for layer in layers)
    return {
        **base_metadata,
        "workflow": "multilayer_attenuation_profile",
        "energy_keV": f"{energy_kev:.12g}",
        "data_source": DATA_SOURCE_LOCAL,
        "interpolation": INTERPOLATION_RULE,
        "edge_policy": EDGE_POLICY,
        "layer_order": "top_to_bottom",
        "layer_compositions_json": json.dumps([layer.composition for layer in layers], sort_keys=True, separators=(",", ":")),
        "layer_mass_fractions_json": json.dumps(layer_mass_fractions, sort_keys=True, separators=(",", ":")),
        "layer_densities_g_cm3_json": json.dumps([layer.density_g_cm3 for layer in layers], separators=(",", ":")),
        "layer_thicknesses_cm_json": json.dumps([layer.thickness_cm for layer in layers], separators=(",", ":")),
        "total_thickness_cm": f"{total_thickness:.16g}",
        "depth_sampling_step_cm": f"{step_cm:.16g}",
        "equation": "tau(E,s)=sum(mu_i(E)*d_i); T=exp(-tau); R=1-T",
        "incident_intensity_normalization": "I0 = 1",
        "scope": "monoenergetic narrow primary-beam attenuation through a fixed multilayer stack",
        "limitations": (
            "No scattered-photon buildup, no fluorescence escape/reabsorption, "
            "and no coupled photon-electron transport are included."
        ),
    }


def run_multilayer_attenuation_profile(config: dict) -> WorkflowResult:
    output_dir: Path = config["output_dir"]
    layers: list[LayerSpec] = config["layers"]
    step_cm: float = config["step_cm"]
    energies_kev: list[float] = config["energies_kev"]
    total_thickness = sum(layer.thickness_cm for layer in layers)
    spatial_rows = _layer_sample_rows(layers, step_cm)
    frames: dict[float, pd.DataFrame] = {}
    output_paths: list[Path] = []

    for energy_kev in energies_kev:
        mass_mu_values: list[float] = []
        linear_mu_values: list[float] = []
        layer_mass_fractions: list[dict[str, float]] = []
        base_metadata: dict[str, str] | None = None
        for layer in layers:
            mass_mu, linear_mu, metadata, mass_fraction_values = coefficient_at_energy_kev(
                layer, "attenuation", energy_kev
            )
            if linear_mu is None:
                raise ValueError("All layers require density_g_cm3")
            mass_mu_values.append(mass_mu)
            linear_mu_values.append(linear_mu)
            layer_mass_fractions.append(mass_fraction_values)
            if base_metadata is None:
                base_metadata = metadata

        completed_tau_by_layer: list[float] = [0.0]
        for linear_mu, layer in zip(linear_mu_values, layers):
            completed_tau_by_layer.append(completed_tau_by_layer[-1] + linear_mu * layer.thickness_cm)

        data_rows: list[dict[str, object]] = []
        for layer_index, layer, global_depth, depth_in_layer, label in spatial_rows:
            idx = layer_index - 1
            tau_in_layer = linear_mu_values[idx] * depth_in_layer
            cumulative_tau = completed_tau_by_layer[idx] + tau_in_layer
            transmission = math.exp(-cumulative_tau)
            attenuation = 1.0 - transmission
            data_rows.append(
                {
                    "Depth_from_incident_surface_cm": global_depth,
                    "Stack_coordinate_from_bottom_cm": total_thickness - global_depth,
                    "Layer_id": layer.id,
                    "Interface_position": label,
                    "Cumulative_optical_depth": cumulative_tau,
                    "Transmission_fraction": transmission,
                    "Attenuation_fraction": attenuation,
                }
            )
        frame = pd.DataFrame(data_rows)
        metadata = _multilayer_metadata(
            energy_kev,
            layers,
            layer_mass_fractions,
            step_cm,
            base_metadata or {},
        )
        path = output_dir / f"{_energy_name(energy_kev)}_multilayer_attenuation_profile.txt"
        output_paths.append(write_dataframe(path, frame, metadata))
        frames[float(energy_kev)] = frame
    output_paths.append(plot_multilayer(output_dir / "multilayer_attenuation_profile.png", frames, layers))
    return WorkflowResult("multilayer_attenuation_profile", output_dir, output_paths)
