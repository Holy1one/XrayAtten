from __future__ import annotations

import numpy as np
import pandas as pd

from conftest import write_yaml
from xrayatten.workflows import run_config


def test_energy_absorption_vs_thickness_workflow(workspace_tmp):
    config = write_yaml(
        workspace_tmp / "abs.yaml",
        f"""
schema_version: 1
workflow: energy_absorption_vs_thickness
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
    path = workspace_tmp / "out" / "10keV_energy_absorption.txt"
    text = path.read_text(encoding="utf-8")
    assert "First-collision absorbed energy estimate" in text
    assert "elemental mass-fraction additivity approximation" in text
    frame = pd.read_csv(path, comment="#", sep="\t")
    assert np.all(frame["First_collision_absorbed_energy_fraction"] <= frame["Attenuation_fraction"] + 1e-15)
    balance = (
        frame["Transmission_fraction"]
        + frame["First_collision_absorbed_energy_fraction"]
        + frame["Removed_not_absorbed_fraction"]
    )
    assert np.allclose(balance, 1.0)
    assert (workspace_tmp / "out" / "energy_absorption_vs_thickness.png").exists()
