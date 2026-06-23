from __future__ import annotations

import math

import pytest

from xrayatten.local_nist import compute_coefficients, mass_fractions, read_elements
from xrayatten.models import MaterialSpec


def test_atomic_ratio_to_mass_fraction_water():
    _, weights = mass_fractions({"H": 2, "O": 1}, "atomic_ratio", read_elements())
    assert math.isclose(sum(weights.values()), 1.0, rel_tol=0.0, abs_tol=1e-12)
    assert weights["O"] > weights["H"]


def test_mass_fraction_must_sum_to_one():
    with pytest.raises(ValueError, match="sum to 1.0"):
        mass_fractions({"H": 0.2, "O": 0.7}, "mass_fraction", read_elements())


def test_unsupported_elements_rejected():
    with pytest.raises(ValueError, match="Z=1-92"):
        mass_fractions({"Np": 1}, "atomic_ratio", read_elements())


def test_density_controls_linear_coefficient():
    material = MaterialSpec(
        id="si",
        composition_basis="atomic_ratio",
        composition={"Si": 1},
        density_g_cm3=2.33,
    )
    result = compute_coefficients(material, coefficient_kind="attenuation", energies_kev=[10])
    assert result.linear_coefficient is not None
    assert result.mass_coefficient[0] * 2.33 == pytest.approx(result.linear_coefficient[0])


def test_no_density_outputs_mass_only():
    material = MaterialSpec(id="si", composition_basis="atomic_ratio", composition={"Si": 1})
    result = compute_coefficients(material, coefficient_kind="attenuation", energies_kev=[10])
    assert result.linear_coefficient is None

