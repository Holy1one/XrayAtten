from __future__ import annotations

import pandas as pd
import pytest

from conftest import write_yaml
from xrayatten.local_nist import compute_coefficients
from xrayatten.models import MaterialSpec
from xrayatten.workflows import run_config


def test_coefficients_workflow_outputs_local_files(workspace_tmp, example_blue_material):
    config = write_yaml(
        workspace_tmp / "coeff.yaml",
        f"""
schema_version: 1
workflow: coefficients
output_dir: {workspace_tmp / "out"}
energies_kev: null
online_comparison: false
materials:
  - id: blue
    composition_basis: atomic_ratio
    composition:
      Si: 1
      O: 2
    density_g_cm3: 2.2
  - id: mass_only
    composition_basis: atomic_ratio
    composition:
      Si: 1
      O: 2
    density_g_cm3: null
""",
    )
    result = run_config(config)
    assert result.workflow == "coefficients"
    assert (workspace_tmp / "out" / "attenuation" / "blue_attenuation_local.txt").exists()
    assert (workspace_tmp / "out" / "energy_absorption" / "blue_energy_absorption_approx_local.png").exists()
    assert not (workspace_tmp / "out" / "online_comparison").exists()
    mass_only = pd.read_csv(workspace_tmp / "out" / "attenuation" / "mass_only_attenuation_local.txt", comment="#", sep="\t")
    assert "Linear_mu_cm_inverse" not in mass_only.columns


def test_coefficients_online_failure_keeps_local_outputs(workspace_tmp, monkeypatch):
    import xrayatten.online_xcom as online_xcom

    def fail_online(*_args, **_kwargs):
        raise RuntimeError("network unavailable")

    monkeypatch.setattr(online_xcom, "run_online_comparison", fail_online)
    config = write_yaml(
        workspace_tmp / "coeff_online.yaml",
        f"""
schema_version: 1
workflow: coefficients
output_dir: {workspace_tmp / "out_online"}
energies_kev: null
online_comparison: true
materials:
  - id: silica
    composition_basis: atomic_ratio
    composition:
      Si: 1
      O: 2
    density_g_cm3: 2.2
""",
    )
    with pytest.raises(RuntimeError, match="Local coefficient outputs were written successfully"):
        run_config(config)
    assert (workspace_tmp / "out_online" / "attenuation" / "silica_attenuation_local.txt").exists()
    assert (
        workspace_tmp / "out_online" / "energy_absorption" / "silica_energy_absorption_approx_local.txt"
    ).exists()


def test_coefficients_workflow_accepts_selected_energy_points(workspace_tmp):
    config = write_yaml(
        workspace_tmp / "coeff_points.yaml",
        f"""
schema_version: 1
workflow: coefficients
output_dir: {workspace_tmp / "out_points"}
energies_kev: [10, 20]
online_comparison: false
materials:
  - id: silica
    composition_basis: atomic_ratio
    composition:
      Si: 1
      O: 2
    density_g_cm3: 2.2
""",
    )
    run_config(config)

    frame = pd.read_csv(workspace_tmp / "out_points" / "attenuation" / "silica_attenuation_local.txt", comment="#", sep="\t")
    assert frame["Energy_keV"].tolist() == [10.0, 20.0]
    assert frame["Edge_side"].tolist() == ["regular", "regular"]


def test_coefficients_workflow_accepts_precision_energy_range(workspace_tmp):
    config = write_yaml(
        workspace_tmp / "coeff_range.yaml",
        f"""
schema_version: 1
workflow: coefficients
output_dir: {workspace_tmp / "out_range"}
energies_kev: null
energy_range_kev:
  start: 10
  stop: 12
precision: fast
online_comparison: false
materials:
  - id: silica
    composition_basis: atomic_ratio
    composition:
      Si: 1
      O: 2
    density_g_cm3: 2.2
""",
    )
    run_config(config)

    frame = pd.read_csv(workspace_tmp / "out_range" / "attenuation" / "silica_attenuation_local.txt", comment="#", sep="\t")
    assert frame["Energy_keV"].tolist() == [10.0, 11.0, 12.0]

    material = MaterialSpec(
        id="silica",
        composition_basis="atomic_ratio",
        composition={"Si": 1, "O": 2},
        density_g_cm3=2.2,
    )
    expected = compute_coefficients(material, coefficient_kind="attenuation", energies_kev=[11.0])
    actual = frame.loc[frame["Energy_keV"] == 11.0, "Mass_mu_over_rho_cm2_g"].iloc[0]
    assert actual == pytest.approx(expected.mass_coefficient[0])


def test_coefficients_range_rejects_empty_precision_grid(workspace_tmp):
    config = write_yaml(
        workspace_tmp / "coeff_empty_range.yaml",
        f"""
schema_version: 1
workflow: coefficients
output_dir: {workspace_tmp / "out_empty_range"}
energies_kev: null
energy_range_kev:
  start: 10.01
  stop: 10.02
precision: fast
online_comparison: false
materials:
  - id: silica
    composition_basis: atomic_ratio
    composition:
      Si: 1
      O: 2
    density_g_cm3: 2.2
""",
    )
    with pytest.raises(Exception, match="no points"):
        run_config(config)
    assert not (workspace_tmp / "out_empty_range").exists()
