from __future__ import annotations

from pathlib import Path
import uuid

import pytest


@pytest.fixture
def workspace_tmp() -> Path:
    root = Path.cwd() / ".tmp" / "pytest_cases"
    root.mkdir(parents=True, exist_ok=True)
    path = root / f"case_{uuid.uuid4().hex}"
    path.mkdir(parents=True, exist_ok=False)
    return path


@pytest.fixture
def example_blue_material() -> dict:
    return {
        "id": "blue_glass",
        "composition_basis": "atomic_ratio",
        "composition": {
            "Si": 50,
            "B": 20,
            "Li": 20,
            "Mg": 10,
            "Y": 2,
            "Cs": 10,
            "Ce": 1,
            "O": 155,
            "F": 9,
        },
        "density_g_cm3": 2.5,
        "density_source": "test configuration",
    }


def write_yaml(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path
