from __future__ import annotations

import pytest

from xrayatten.models import MaterialSpec
from xrayatten.online_xcom import compute_online_attenuation, parse_xcom_html


class _Response:
    status_code = 200

    def __init__(self, content: bytes):
        self.content = content


class _Session:
    headers: dict[str, str] = {}

    def __init__(self, content: bytes):
        self.content = content

    def post(self, *_args, **_kwargs):
        return _Response(self.content)


def test_parse_xcom_fixture():
    html = open("tests/fixtures/xcom_response_fixture.html", "rb").read()
    table = parse_xcom_html(html)
    assert table.energy_mev.tolist() == [0.01, 0.02, 0.02]
    assert table.edge_side.tolist() == ["regular", "below", "above"]


def test_online_comparison_uses_mock_session(workspace_tmp):
    html = open("tests/fixtures/xcom_response_fixture.html", "rb").read()
    material = MaterialSpec(
        id="water",
        composition_basis="atomic_ratio",
        composition={"H": 2, "O": 1},
        density_g_cm3=1.0,
        density_source="test",
    )
    paths = compute_online_attenuation(
        material,
        energies_kev=[10, 20],
        output_dir=workspace_tmp / "online",
        session=_Session(html),
    )
    assert all(path.exists() for path in paths)
    assert "online" in paths[0].read_text(encoding="utf-8")


def test_online_rejects_mass_fraction(workspace_tmp):
    material = MaterialSpec(
        id="water",
        composition_basis="mass_fraction",
        composition={"H": 0.1, "O": 0.9},
    )
    with pytest.raises(ValueError, match="atomic_ratio"):
        compute_online_attenuation(material, energies_kev=[10], output_dir=workspace_tmp)
