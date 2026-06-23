"""Small physics helpers shared by all local workflows."""

from __future__ import annotations

import math

import numpy as np


def linear_coefficient(mass_mu_over_rho_cm2_g, density_g_cm3: float):
    return np.asarray(mass_mu_over_rho_cm2_g, dtype=float) * float(density_g_cm3)


def transmission_fraction(linear_mu_cm_inverse, thickness_cm):
    return np.exp(-np.asarray(linear_mu_cm_inverse, dtype=float) * np.asarray(thickness_cm, dtype=float))


def attenuation_fraction(linear_mu_cm_inverse, thickness_cm):
    return 1.0 - transmission_fraction(linear_mu_cm_inverse, thickness_cm)


def first_collision_absorbed_energy_fraction(linear_mu_cm_inverse, linear_mu_en_cm_inverse, thickness_cm):
    mu = np.asarray(linear_mu_cm_inverse, dtype=float)
    mu_en = np.asarray(linear_mu_en_cm_inverse, dtype=float)
    if np.any(mu_en > mu * (1.0 + 1e-10)):
        raise ValueError("mu_en must not exceed total attenuation mu")
    return (mu_en / mu) * (1.0 - np.exp(-mu * np.asarray(thickness_cm, dtype=float)))


def validate_positive_finite(value: object, name: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return number


def thickness_grid(maximum_cm: float, step_cm: float) -> np.ndarray:
    maximum = validate_positive_finite(maximum_cm, "maximum_cm")
    step = validate_positive_finite(step_cm, "step_cm")
    grid = np.arange(0.0, maximum, step, dtype=float)
    if grid.size == 0 or not math.isclose(grid[-1], maximum, rel_tol=0.0, abs_tol=1e-14):
        grid = np.append(grid, maximum)
    return grid
