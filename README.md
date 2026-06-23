# XrayAtten

YAML-driven X-ray attenuation workflows based on a packaged local NIST v1.4
coefficient snapshot.

**Language:** English | [Chinese](README.zh-CN.md)

## Overview

`xrayatten` computes X-ray attenuation quantities for user-defined materials
and layered stacks. Formal workflows use the bundled local NIST v1.4 snapshot;
online NIST XCOM queries are available only as optional comparison outputs.

The package currently supports:

- material coefficient tables for total attenuation and approximate energy
  absorption;
- primary-beam attenuation versus thickness at user-selected energies;
- first-collision absorbed-energy estimates versus thickness;
- cumulative primary-beam attenuation profiles through multilayer stacks.

## Highlights

- YAML-first workflows for reproducible calculations.
- Packaged local data under `src/xrayatten/data`, including 92 element
  coefficient tables and SHA-256 manifest validation.
- One absorption-edge-aware log-log interpolation implementation shared by all
  workflows.
- Dense coefficient grids with selectable `precision`.
- Provenance metadata in every TXT output.
- Optional online NIST XCOM comparison for coefficient workflows.
- Developer tests for data integrity, interpolation, physics rules, YAML
  validation, and workflow outputs.

## Installation

Python 3.10 or newer is required.

Install the released package:

```bash
python -m pip install xrayatten
```

Install optional online comparison support:

```bash
python -m pip install "xrayatten[online]"
```

Install from source for development:

```bash
python -m pip install -e ".[dev,online]"
```

If you only need local calculation from a source checkout:

```bash
python -m pip install -e .
```

## Quick Start

Validate the bundled local NIST data:

```bash
xrayatten validate-data
```

Run an example workflow:

```bash
xrayatten run configs/examples/coefficients_batch.yaml
```

The module form is equivalent and useful before console scripts are on `PATH`:

```bash
python -m xrayatten validate-data
python -m xrayatten run configs/examples/coefficients_batch.yaml
```

## Workflows

| Workflow | Purpose | Energy policy |
| --- | --- | --- |
| `coefficients` | Write dense material coefficient tables for total attenuation and approximate energy absorption | Does not accept a user energy list; writes the full precision grid |
| `attenuation_vs_thickness` | Compute primary-beam transmission and attenuation versus thickness | Uses user-provided `energies_kev` |
| `energy_absorption_vs_thickness` | Estimate first-collision absorbed-energy fraction versus thickness | Uses user-provided `energies_kev` |
| `multilayer_attenuation_profile` | Compute cumulative primary-beam attenuation through fixed layers | Uses user-provided `energies_kev` |

Example configuration:

```yaml
schema_version: 1
workflow: coefficients
output_dir: results/my_coefficients
energies_kev: null
precision: low
online_comparison: false

materials:
  - id: blue_glass
    composition_basis: atomic_ratio
    composition:
      Si: 50
      B: 20
      Li: 20
      Mg: 10
      Y: 2
      Cs: 10
      Ce: 1
      O: 155
      F: 9
    density_g_cm3: 2.5
    density_source: project configuration
```

More runnable examples are available in [configs/examples](configs/examples).

## Calculation Model

Material composition can be provided as either:

- `atomic_ratio`: converted to mass fractions using `elements.csv`;
- `mass_fraction`: already normalized to exactly 1.

Material coefficients use mass-fraction additivity:

```text
(mu/rho)_material = sum_i w_i * (mu/rho)_i
(mu_en/rho)_material ~= sum_i w_i * (mu_en/rho)_i
```

If `density_g_cm3` is provided, linear coefficients are computed as:

```text
mu = (mu/rho) * density
```

Thickness workflows then use:

```text
T(x) = exp(-mu*x)
R(x) = 1 - T(x)
A_E(x) = (mu_en/mu) * (1 - exp(-mu*x))
```

`A_E` is a first-collision absorbed-energy estimate, not a full deposited-dose
or detector-efficiency model.

## Interpolation

All workflows use the same interpolation implementation.

1. Element tables are loaded from the local NIST v1.4 snapshot and checked
   against manifest SHA-256 values.
2. Energies are converted from MeV to keV for calculation and output.
3. Repeated energy rows are treated as absorption edges.
4. Full coefficient tables preserve both `below` and `above` edge rows.
5. Exact single-energy queries at an edge use the `above` value.
6. Non-edge points use piecewise log-log interpolation.
7. Interpolation never crosses an absorption edge and never extrapolates.
8. Non-finite interpolation results are rejected before output.

The `coefficients` workflow builds a dense grid from all original NIST energy
points, all absorption-edge rows, and additional gap-fill points controlled by
`precision`.

| `precision` | Relative density | Typical use |
| --- | --- | --- |
| `super` | Highest | Archival high-resolution curves |
| `high` | High | Detailed plotting |
| `medium` | Medium | Balanced output size and smoothness |
| `low` | Default | Routine plotting and Origin import |
| `fast` | Lowest | Quick preview |
| `direct` | Original points only | NIST points and absorption-edge rows, without added grid points |

## Output

All TXT files include metadata header lines followed by tab-separated columns.

Coefficient outputs:

- attenuation: `Energy_keV`, `Edge_side`, `Mass_mu_over_rho_cm2_g`, and
  `Linear_mu_cm_inverse` when density is available;
- approximate energy absorption: `Energy_keV`, `Edge_side`,
  `Mass_mu_en_over_rho_approx_cm2_g`, and
  `Linear_mu_en_approx_cm_inverse` when density is available.

Thickness attenuation outputs:

- `Thickness_cm`
- `Mass_mu_over_rho_cm2_g`
- `Linear_mu_cm_inverse`
- `Transmission_fraction`
- `Attenuation_fraction`

First-collision energy-absorption outputs:

- `Thickness_cm`
- `Transmission_fraction`
- `Attenuation_fraction`
- `First_collision_absorbed_energy_fraction`
- `Removed_not_absorbed_fraction`

Multilayer outputs:

- `Depth_from_incident_surface_cm`
- `Stack_coordinate_from_bottom_cm`
- `Layer_id`
- `Interface_position`
- `Cumulative_optical_depth`
- `Transmission_fraction`
- `Attenuation_fraction`

## Configuration Rules

| Field | Rule |
| --- | --- |
| `schema_version` | Required; currently `1` |
| `workflow` | One of `coefficients`, `attenuation_vs_thickness`, `energy_absorption_vs_thickness`, `multilayer_attenuation_profile` |
| `output_dir` | Required; absolute paths are used as-is, relative paths resolve from the current working directory |
| `composition_basis` | `atomic_ratio` or `mass_fraction` |
| `density_g_cm3` | Positive number, or `null` only for `coefficients` |
| `energies_kev` | Must be `null` for `coefficients`; required for attenuation, energy-absorption, and multilayer workflows |
| `precision` | `direct`, `super`, `high`, `medium`, `low`, or `fast`; used only by `coefficients` |
| `online_comparison` | Optional for `coefficients`; supports `atomic_ratio` materials only |

Important behavior:

- Configuration validation happens before output directories are created.
- Energies are limited to the supported local NIST range; extrapolation is
  rejected.
- If local coefficient outputs are written but online XCOM comparison fails,
  the error message explicitly says that local outputs were completed.

## Project Layout

```text
src/xrayatten/                 package source
src/xrayatten/data/            bundled NIST data and element metadata
configs/examples/              runnable YAML examples
docs/                          method, configuration, reproducibility, and migration notes
tests/                         developer test suite
pyproject.toml                 package metadata and optional dependencies
```

## Testing

Tests are for developers, maintainers, CI, and source reviewers. They are not
required for normal `pip install xrayatten` users.

```bash
python -m pytest -q
python -m pytest --cov=xrayatten --cov-report=term-missing
```

Validate local data:

```bash
python -m xrayatten validate-data
```

Syntax-check without writing `__pycache__`:

```bash
python -B -c "from pathlib import Path; files=list(Path('src').rglob('*.py'))+list(Path('tests').rglob('*.py')); [compile(p.read_text(encoding='utf-8'), str(p), 'exec') for p in files]"
```

## Scientific Scope

This package computes monoenergetic narrow-beam attenuation using local NIST
coefficient tables and deterministic formulas. It does not model scattered
photon buildup, fluorescence escape or reabsorption, bremsstrahlung transport,
or coupled photon-electron Monte Carlo transport.

Use [docs/METHODS_VALIDATION.md](docs/METHODS_VALIDATION.md) when preparing
publication text or method descriptions.

## Documentation

| Document | Description |
| --- | --- |
| [README.zh-CN.md](README.zh-CN.md) | Chinese README |
| [docs/METHODS_VALIDATION.md](docs/METHODS_VALIDATION.md) | Scientific formulas, limitations, and validation notes |
| [docs/CONFIGURATION.md](docs/CONFIGURATION.md) | YAML configuration reference |
| [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) | Data versioning, checksums, and output metadata |
| [docs/MIGRATION.md](docs/MIGRATION.md) | Migration notes from older scripts to YAML workflows |
