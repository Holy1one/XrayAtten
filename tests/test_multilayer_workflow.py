from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from conftest import write_yaml
from xrayatten.workflows import run_config


def _multilayer_yaml(
    path,
    output_dir,
    *,
    duplicate_id=False,
    step=0.05,
    bad_thickness=False,
    include_energies=True,
):
    second_id = "top" if duplicate_id else "bottom"
    thickness = -0.1 if bad_thickness else 0.1
    energies = "energies_kev: [15, 30]\n" if include_energies else ""
    return write_yaml(
        path,
        f"""
schema_version: 1
workflow: multilayer_attenuation_profile
output_dir: {output_dir}
{energies}
depth_sampling:
  step_cm: {step}
layers:
  - id: top
    composition_basis: atomic_ratio
    composition:
      Si: 1
      O: 2
    density_g_cm3: 2.2
    thickness_cm: 0.1
  - id: {second_id}
    composition_basis: atomic_ratio
    composition:
      C: 1
    density_g_cm3: 1.3
    thickness_cm: {thickness}
""",
    )


def test_multilayer_profile_outputs_configured_energy_files(workspace_tmp):
    config = _multilayer_yaml(workspace_tmp / "multi.yaml", workspace_tmp / "out")
    run_config(config)
    txt_files = sorted((workspace_tmp / "out").glob("*keV_multilayer_attenuation_profile.txt"))
    assert len(txt_files) == 2
    assert {path.name for path in txt_files} == {
        "15keV_multilayer_attenuation_profile.txt",
        "30keV_multilayer_attenuation_profile.txt",
    }
    assert (workspace_tmp / "out" / "multilayer_attenuation_profile.png").exists()
    frame = pd.read_csv(workspace_tmp / "out" / "15keV_multilayer_attenuation_profile.txt", comment="#", sep="\t")
    total = 0.2
    assert frame["Depth_from_incident_surface_cm"].iloc[0] == pytest.approx(0.0)
    assert frame["Stack_coordinate_from_bottom_cm"].iloc[0] == pytest.approx(total)
    assert frame["Transmission_fraction"].iloc[0] == pytest.approx(1.0)
    assert frame["Attenuation_fraction"].iloc[0] == pytest.approx(0.0)
    assert frame["Depth_from_incident_surface_cm"].iloc[-1] == pytest.approx(total)
    assert frame["Stack_coordinate_from_bottom_cm"].iloc[-1] == pytest.approx(0.0)
    assert np.allclose(
        frame["Stack_coordinate_from_bottom_cm"],
        total - frame["Depth_from_incident_surface_cm"],
    )
    assert np.all(np.diff(frame["Cumulative_optical_depth"]) >= -1e-15)
    assert np.all(np.diff(frame["Transmission_fraction"]) <= 1e-15)
    assert np.all(np.diff(frame["Attenuation_fraction"]) >= -1e-15)
    assert np.allclose(frame["Transmission_fraction"] + frame["Attenuation_fraction"], 1.0)
    assert {"stack_top", "layer_exit", "layer_entry", "stack_bottom"}.issubset(set(frame["Interface_position"]))


def test_multilayer_validation_errors(workspace_tmp):
    with pytest.raises(Exception, match="unique"):
        run_config(_multilayer_yaml(workspace_tmp / "dup.yaml", workspace_tmp / "dup_out", duplicate_id=True))
    with pytest.raises(Exception, match="thickness_cm"):
        run_config(_multilayer_yaml(workspace_tmp / "bad_thick.yaml", workspace_tmp / "bad_thick_out", bad_thickness=True))
    with pytest.raises(Exception, match="depth_sampling.step_cm"):
        run_config(_multilayer_yaml(workspace_tmp / "bad_step.yaml", workspace_tmp / "bad_step_out", step=0))
    with pytest.raises(Exception, match="energies_kev"):
        run_config(
            _multilayer_yaml(
                workspace_tmp / "missing_energy.yaml",
                workspace_tmp / "missing_energy_out",
                include_energies=False,
            )
        )
