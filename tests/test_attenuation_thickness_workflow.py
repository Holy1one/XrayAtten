from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from conftest import write_yaml
from xrayatten.workflows import run_config


def test_attenuation_vs_thickness_workflow(workspace_tmp):
    config = write_yaml(
        workspace_tmp / "atten.yaml",
        f"""
schema_version: 1
workflow: attenuation_vs_thickness
output_dir: {workspace_tmp / "out"}
material:
  id: silica
  composition_basis: atomic_ratio
  composition:
    Si: 1
    O: 2
  density_g_cm3: 2.2
energies_kev: [10, 20]
thickness:
  maximum_cm: 0.01
  step_cm: 0.004
""",
    )
    run_config(config)
    frame = pd.read_csv(workspace_tmp / "out" / "10keV_attenuation.txt", comment="#", sep="\t")
    assert frame["Thickness_cm"].iloc[-1] == pytest.approx(0.01)
    assert np.allclose(frame["Transmission_fraction"] + frame["Attenuation_fraction"], 1.0)
    assert np.all(np.diff(frame["Attenuation_fraction"]) >= -1e-15)
    assert (workspace_tmp / "out" / "attenuation_vs_thickness.png").exists()


def test_attenuation_requires_density_before_output(workspace_tmp):
    config = write_yaml(
        workspace_tmp / "bad.yaml",
        f"""
schema_version: 1
workflow: attenuation_vs_thickness
output_dir: {workspace_tmp / "bad_out"}
material:
  id: silica
  composition_basis: atomic_ratio
  composition:
    Si: 1
    O: 2
  density_g_cm3: null
energies_kev: [10]
thickness:
  maximum_cm: 0.01
  step_cm: 0.004
""",
    )
    with pytest.raises(Exception, match="density_g_cm3"):
        run_config(config)
    assert not (workspace_tmp / "bad_out").exists()
