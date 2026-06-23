from __future__ import annotations

import json

import pytest

from xrayatten.local_nist import compute_coefficients
from xrayatten.models import MaterialSpec


def test_local_water_reference_fixture():
    fixture = json.loads(open("tests/fixtures/nist_water_reference.json", encoding="utf-8").read())
    material = MaterialSpec(
        id="water",
        composition_basis=fixture["composition_basis"],
        composition=fixture["composition"],
    )
    energies = [point["energy_kev"] for point in fixture["points"]]
    attenuation = compute_coefficients(material, coefficient_kind="attenuation", energies_kev=energies)
    absorption = compute_coefficients(material, coefficient_kind="energy_absorption_approx", energies_kev=energies)
    tol = fixture["relative_tolerance"]
    for index, point in enumerate(fixture["points"]):
        assert attenuation.mass_coefficient[index] == pytest.approx(point["mass_mu_over_rho_cm2_g"], rel=tol)
        assert absorption.mass_coefficient[index] == pytest.approx(point["mass_mu_en_over_rho_approx_cm2_g"], rel=tol)
