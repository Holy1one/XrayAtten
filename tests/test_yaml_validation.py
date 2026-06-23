from __future__ import annotations

import pytest

from conftest import write_yaml
from xrayatten.config import resolve_output_dir
from xrayatten.workflows import run_config


def test_yaml_rejects_coefficient_file_key_before_output(workspace_tmp):
    config = write_yaml(
        workspace_tmp / "bad.yaml",
        f"""
schema_version: 1
workflow: attenuation_vs_thickness
output_dir: {workspace_tmp / "out"}
coefficient_file: results/old.txt
material:
  id: silica
  composition_basis: atomic_ratio
  composition:
    Si: 1
    O: 2
  density_g_cm3: 2.2
energies_kev: [10]
thickness:
  maximum_cm: 0.01
  step_cm: 0.001
""",
    )
    with pytest.raises(Exception, match="coefficient_file"):
        run_config(config)
    assert not (workspace_tmp / "out").exists()


def test_coefficients_rejects_explicit_energies(workspace_tmp):
    config = write_yaml(
        workspace_tmp / "bad_coeff.yaml",
        f"""
schema_version: 1
workflow: coefficients
output_dir: {workspace_tmp / "out"}
energies_kev: [10, 20]
materials:
  - id: silica
    composition_basis: atomic_ratio
    composition:
      Si: 1
      O: 2
    density_g_cm3: 2.2
""",
    )
    with pytest.raises(Exception, match="energies_kev"):
        run_config(config)
    assert not (workspace_tmp / "out").exists()


def test_relative_output_dir_resolves_from_current_working_directory(workspace_tmp, monkeypatch):
    monkeypatch.chdir(workspace_tmp)
    assert resolve_output_dir("out") == workspace_tmp / "out"


def test_online_comparison_rejects_mass_fraction_before_output(workspace_tmp):
    config = write_yaml(
        workspace_tmp / "bad_online.yaml",
        f"""
schema_version: 1
workflow: coefficients
output_dir: {workspace_tmp / "out_online_bad"}
energies_kev: null
online_comparison: true
materials:
  - id: water
    composition_basis: mass_fraction
    composition:
      H: 0.1118983441
      O: 0.8881016559
    density_g_cm3: 1.0
""",
    )
    with pytest.raises(Exception, match="atomic_ratio"):
        run_config(config)
    assert not (workspace_tmp / "out_online_bad").exists()
