from __future__ import annotations

import numpy as np
import pytest

from xrayatten.local_nist import compute_coefficients
from xrayatten.interpolation import interpolate_element, interpolate_element_to_grid
from xrayatten.models import MaterialSpec


def test_exact_absorption_edge_outputs_below_and_above():
    material = MaterialSpec(id="si", composition_basis="atomic_ratio", composition={"Si": 1})
    result = compute_coefficients(material, coefficient_kind="attenuation", energies_kev=[1.8389])
    assert [row.edge_side for row in result.rows] == ["below", "above"]
    assert result.mass_coefficient[0] == pytest.approx(309.2)
    assert result.mass_coefficient[1] == pytest.approx(3192.0)


def test_extrapolation_is_forbidden():
    material = MaterialSpec(id="si", composition_basis="atomic_ratio", composition={"Si": 1})
    with pytest.raises(ValueError):
        compute_coefficients(material, coefficient_kind="attenuation", energies_kev=[0.5])


def test_point_and_grid_interpolation_use_same_logic():
    material = MaterialSpec(id="si", composition_basis="atomic_ratio", composition={"Si": 1})
    result = compute_coefficients(material, coefficient_kind="attenuation", energies_kev=[1.8389, 10])
    table = result.element_tables["Si"]
    point_values = interpolate_element(table, table.mass_attenuation, result.rows)
    grid_values = interpolate_element_to_grid(
        table,
        table.mass_attenuation,
        np.asarray([row.energy_kev for row in result.rows], dtype=float),
        [row.edge_side for row in result.rows],
    )
    assert np.allclose(point_values, grid_values)
