from __future__ import annotations

import numpy as np

from xrayatten.local_nist import load_manifest, read_element_table, read_elements, validate_data


def test_local_manifest_and_tables_validate():
    result = validate_data()
    assert result["schema"] == "nist_xraymasscoef_snapshot_manifest_v1"
    assert result["elements"] == 92
    assert result["version"] == "1.4"


def test_elements_table_and_fluorine_correction():
    elements = read_elements()
    assert len(elements) == 118
    assert elements["H"].atomic_number == 1
    assert elements["U"].atomic_number == 92

    manifest = load_manifest()
    assert any(item["symbol"] == "F" and item["corrected_value"] == 0.06 for item in manifest["local_corrections"])
    fluorine = read_element_table("F", manifest)
    assert np.any(np.isclose(fluorine.energy_mev, 0.06))
    assert np.all(np.diff(fluorine.energy_mev) >= 0)

