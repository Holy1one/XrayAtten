"""Domain models for materials, coefficients, workflows, and provenance."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np


@dataclass(frozen=True)
class MaterialSpec:
    id: str
    composition_basis: str
    composition: dict[str, float]
    density_g_cm3: float | None = None
    density_source: str = "not specified by user"


@dataclass(frozen=True)
class LayerSpec(MaterialSpec):
    thickness_cm: float = 0.0


@dataclass(frozen=True)
class ElementRecord:
    atomic_number: int
    symbol: str
    atomic_weight: float
    atomic_weight_source: str


@dataclass(frozen=True)
class CoefficientTable:
    symbol: str
    energy_mev: np.ndarray
    mass_attenuation: np.ndarray
    mass_energy_absorption: np.ndarray
    sha256: str
    path: Path


@dataclass(frozen=True)
class OutputRow:
    energy_mev: float
    edge_side: str

    @property
    def energy_kev(self) -> float:
        return self.energy_mev * 1000.0


@dataclass
class CoefficientResult:
    material: MaterialSpec
    coefficient_kind: str
    rows: list[OutputRow]
    mass_coefficient: np.ndarray
    linear_coefficient: np.ndarray | None
    mass_fractions: dict[str, float]
    metadata: dict[str, str]
    warnings: list[str] = field(default_factory=list)
    limitations: list[str] = field(default_factory=list)
    output_paths: dict[str, Path] = field(default_factory=dict)
    element_tables: dict[str, "CoefficientTable"] = field(default_factory=dict)


@dataclass
class Provenance:
    data_source: str
    coefficient_kind: str
    metadata: dict[str, str]


@dataclass
class WorkflowResult:
    workflow: str
    output_dir: Path
    output_paths: list[Path]
    metadata: dict[str, Any] = field(default_factory=dict)
