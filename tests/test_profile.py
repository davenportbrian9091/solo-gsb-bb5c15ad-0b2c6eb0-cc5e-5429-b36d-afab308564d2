import math

import pytest
from fastapi.testclient import TestClient

from app.config import settings
from app.envelope import isotope_envelope
from app.parser import parse_formula
from app.profile import (
    KERNEL_WIDTH_SIGMAS,
    ProfileError,
    gaussian_sigma,
    instrument_profile,
    simulate_profile,
)

TOL = 0.0001
PRUNE = 1e-12


@pytest.fixture()
def client(tmp_path, monkeypatch):
    # Isolate each test run with its own SQLite file.
    monkeypatch.setattr(settings, "db_path", str(tmp_path / "test.db"))
    from app.main import app

    with TestClient(app) as c:
        yield c


def test_gaussian_sigma_matches_fwhm_definition():
    # FWHM = mass / R; sigma = FWHM / (2*sqrt(2*ln2))
    sigma = gaussian_sigma(500.0, 10000)
    assert sigma == pytest.approx(0.05 / 2.3548200450309493)


def test_single_stick_gaussian_shape():
    # One stick at 100 Da, R=1000 -> FWHM 0.1; grid aligned to the stick.
    result = instrument_profile(
        [(100.0, 100.0)], 1000.0, 0.05, 99.5, 100.5, 20000
    )
    assert result.mass_min == 99.5
    assert result.mass_max == 100.5
    assert len(result.points) == 21
    centre = result.points[10]
    assert centre[0] == pytest.approx(100.0)
    assert centre[1] == pytest.approx(100.0)
    # FWHM/2 = 0.05 -> half maximum one grid step either side
    assert result.points[9][1] == pytest.approx(50.0, abs=1e-6)
    assert result.points[11][1] == pytest.approx(50.0, abs=1e-6)
    # beyond +/-4 sigma (~0.1699) the kernel is not sampled at all
    assert result.points[0][1] == 0.0  # 99.5 is 0.5 away
    assert result.points[6][1] == 0.0  # 99.8 is 0.2 away
    assert result.points[7][1] > 0.0  # 99.85 is 0.15 away (< 4 sigma)


def test_default_range_spans_all_kernel_supports():
    sticks = [(100.0, 50.0), (200.0, 100.0)]
    result = instrument_profile(sticks, 1000.0, 1.0, None, None, 20000)
    assert result.mass_min == pytest.approx(
        100.0 - KERNEL_WIDTH_SIGMAS * gaussian_sigma(100.0, 1000.0)
    )
    assert result.mass_max == pytest.approx(
        200.0 + KERNEL_WIDTH_SIGMAS * gaussian_sigma(200.0, 1000.0)
    )


def test_kernel_truncated_at_range_bounds_keeps_inside():
    # The range cuts through the kernel: the in-range contribution is kept
    # as-is (no renormalisation of the truncated kernel).
    result = instrument_profile(
        [(100.0, 100.0)], 1000.0, 0.05, 100.0, 100.1, 20000
    )
    assert [k for k, _ in enumerate(result.points)] == [0, 1, 2]
    assert result.points[0][0] == pytest.approx(100.0)
    assert result.points[0][1] == pytest.approx(100.0)  # centre at the edge
    # plain Gaussian value at +FWHM/2, not rescaled for the dropped tail
    assert result.points[1][1] == pytest.approx(50.0, abs=1e-6)
    # at +FWHM a Gaussian is worth 1/16 of its peak
    assert result.points[2][1] == pytest.approx(100.0 / 16.0, abs=1e-6)


def test_overlapping_kernels_accumulate():
    # Per-stick sigma follows FWHM = mass / R, so the two kernels differ
    # slightly; interior grid points receive both sticks' contributions.
    sticks = [(100.0, 100.0), (100.1, 100.0)]
    result = instrument_profile(sticks, 1000.0, 0.05, 99.9, 100.2, 20000)
    s1 = gaussian_sigma(100.0, 1000.0)
    s2 = gaussian_sigma(100.1, 1000.0)

    def summed(x):
        return sum(
            100.0 * math.exp(-0.5 * ((x - m) / s) ** 2)
            for m, s in ((100.0, s1), (100.1, s2))
        )

    # points 2..4 (100.0, 100.05, 100.1) lie inside both +/-4-sigma windows
    top = max(summed(99.9 + k * 0.05) for k in (2, 3, 4))
    for k in (2, 3, 4):
        mass, intensity = result.points[k]
        assert mass == pytest.approx(99.9 + k * 0.05)
        assert intensity == pytest.approx(100.0 * summed(mass) / top, abs=1e-9)


def test_empty_stick_list_raises():
    with pytest.raises(ProfileError):
        instrument_profile([], 1000.0, 0.1, None, None, 20000)


def test_invalid_range_raises_with_field():
    with pytest.raises(ProfileError) as excinfo:
        instrument_profile([(100.0, 1.0)], 1000.0, 0.1, 100.0, 100.0, 20000)
    assert excinfo.value.field == "mass_max"
    # explicit lower bound above the (defaulted) upper boundary
    with pytest.raises(ProfileError) as excinfo:
        instrument_profile([(100.0, 1.0)], 1000.0, 0.1, 200.0, None, 20000)
    assert excinfo.value.field == "mass_min"
    # explicit upper bound below the (defaulted) lower boundary
    with pytest.raises(ProfileError) as excinfo:
        instrument_profile([(100.0, 1.0)], 1000.0, 0.1, None, 50.0, 20000)
    assert excinfo.value.field == "mass_max"


def test_grid_point_limit_raises():
    with pytest.raises(ProfileError) as excinfo:
        instrument_profile([(100.0, 1.0)], 1000.0, 0.001, 99.0, 101.0, 100)
    assert excinfo.value.field == "mass_step"
    # exactly at the limit is accepted
    result = instrument_profile([(100.0, 1.0)], 1000.0, 0.02, 99.0, 101.0, 101)
    assert len(result.points) == 101


def test_range_without_contribution_raises():
    # range above every stick's +4-sigma support
    with pytest.raises(ProfileError) as excinfo:
        instrument_profile([(100.0, 1.0)], 1000.0, 0.1, 200.0, 300.0, 20000)
    assert excinfo.value.field == "mass_min"
    # range below every stick's -4-sigma support
    with pytest.raises(ProfileError) as excinfo:
        instrument_profile([(100.0, 1.0)], 1000.0, 0.1, 1.0, 2.0, 20000)
    assert excinfo.value.field == "mass_max"


def test_simulate_profile_uses_full_envelope():
    result = simulate_profile(
        parse_formula("H2O"), 10000.0, 0.001, None, None, TOL, PRUNE, 20000
    )
    env = isotope_envelope(parse_formula("H2O"), TOL, PRUNE)
    lo = min(m - KERNEL_WIDTH_SIGMAS * gaussian_sigma(m, 10000.0) for m, _ in env)
    hi = max(m + KERNEL_WIDTH_SIGMAS * gaussian_sigma(m, 10000.0) for m, _ in env)
    assert result.mass_min == pytest.approx(lo)
    assert result.mass_max == pytest.approx(hi)
    assert max(i for _, i in result.points) == pytest.approx(100.0)
    masses = [m for m, _ in result.points]
    assert masses == sorted(masses)


def _profile_payload(**overrides):
    payload = {
        "formula": "C6H12O6",
        "resolving_power": 60000,
        "mass_step": 0.0005,
    }
    payload.update(overrides)
    return payload


def test_profile_defaults_cover_envelope(client):
    resp = client.post("/api/v1/profile", json=_profile_payload())
    assert resp.status_code == 200
    body = resp.json()
    assert body["formula"] == "C6H12O6"
    assert body["composition"] == {"C": 6, "H": 12, "O": 6}
    assert body["labels"] == {}
    assert body["point_count"] == len(body["profile"])
    masses = [p["mass"] for p in body["profile"]]
    assert masses == sorted(masses)
    assert max(p["intensity"] for p in body["profile"]) == pytest.approx(100.0)
    # default range: outermost +/-4-sigma boundaries over all envelope sticks
    env = isotope_envelope(
        parse_formula("C6H12O6"), settings.merge_tolerance, settings.prune_threshold
    )
    lo = min(m - KERNEL_WIDTH_SIGMAS * gaussian_sigma(m, 60000) for m, _ in env)
    hi = max(m + KERNEL_WIDTH_SIGMAS * gaussian_sigma(m, 60000) for m, _ in env)
    assert body["mass_min"] == pytest.approx(lo, abs=1e-6)
    assert body["mass_max"] == pytest.approx(hi, abs=1e-6)
    # the grid starts at mass_min and advances by mass_step
    assert masses[0] == pytest.approx(body["mass_min"])
    assert masses[1] - masses[0] == pytest.approx(0.0005, abs=1e-9)


def test_profile_explicit_range(client):
    resp = client.post(
        "/api/v1/profile",
        json=_profile_payload(mass_min=180.0, mass_max=180.1, mass_step=0.001),
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["mass_min"] == 180.0
    assert body["mass_max"] == 180.1
    assert body["point_count"] == 101
    masses = [p["mass"] for p in body["profile"]]
    assert masses[0] == pytest.approx(180.0)
    assert masses[-1] == pytest.approx(180.1)


def test_profile_fwhm_matches_resolving_power(client):
    # Glucose M+0 at ~180.063388, R=10000 -> FWHM ~0.018, well isolated.
    resp = client.post(
        "/api/v1/profile",
        json=_profile_payload(
            resolving_power=10000, mass_step=0.0005, mass_min=179.9, mass_max=180.3
        ),
    )
    assert resp.status_code == 200
    pts = resp.json()["profile"]
    top = max(p["intensity"] for p in pts)
    above = [p["mass"] for p in pts if p["intensity"] >= top / 2]
    width = max(above) - min(above)
    assert width == pytest.approx(180.063388 / 10000, abs=0.0015)


def test_profile_with_label(client):
    resp = client.post("/api/v1/profile", json=_profile_payload(formula="C6H6[13C]2"))
    assert resp.status_code == 200
    body = resp.json()
    assert body["labels"] == {"[13C]": 2}
    # the labelled envelope is shifted by two 13C masses (~26 Da), well
    # past the natural C6H6 envelope (which ends below 85 Da)
    assert body["mass_min"] > 100.0
    env = isotope_envelope(
        parse_formula("C6H6[13C]2"), settings.merge_tolerance, settings.prune_threshold
    )
    lo = min(m - KERNEL_WIDTH_SIGMAS * gaussian_sigma(m, 60000) for m, _ in env)
    hi = max(m + KERNEL_WIDTH_SIGMAS * gaussian_sigma(m, 60000) for m, _ in env)
    assert body["mass_min"] == pytest.approx(lo, abs=1e-6)
    assert body["mass_max"] == pytest.approx(hi, abs=1e-6)


def test_profile_single_point_grid(client):
    resp = client.post(
        "/api/v1/profile",
        json=_profile_payload(
            resolving_power=1000, mass_min=180.0, mass_max=180.05, mass_step=0.1
        ),
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["point_count"] == 1
    assert body["profile"][0]["mass"] == pytest.approx(180.0)
    assert body["profile"][0]["intensity"] == pytest.approx(100.0)


def test_profile_invalid_formula(client):
    resp = client.post("/api/v1/profile", json=_profile_payload(formula="C6Xx2"))
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "invalid_formula"
    assert detail["position"] == 2


def test_profile_parameter_bounds(client):
    for override in (
        {"resolving_power": 0},
        {"resolving_power": -5},
        {"mass_step": 0},
        {"mass_step": -0.1},
    ):
        resp = client.post("/api/v1/profile", json=_profile_payload(**override))
        assert resp.status_code == 422, override
    # missing required fields
    assert client.post("/api/v1/profile", json={"formula": "H2O"}).status_code == 422


def test_profile_inverted_range_422(client):
    resp = client.post(
        "/api/v1/profile", json=_profile_payload(mass_min=200.0, mass_max=100.0)
    )
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "invalid_profile"
    assert detail["field"] == "mass_max"


def test_profile_bound_outside_envelope_422(client):
    # explicit lower bound above the whole envelope
    resp = client.post("/api/v1/profile", json=_profile_payload(mass_min=500.0))
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "invalid_profile"
    assert detail["field"] == "mass_min"
    # explicit upper bound below the whole envelope
    resp = client.post("/api/v1/profile", json=_profile_payload(mass_max=10.0))
    assert resp.status_code == 422
    assert resp.json()["detail"]["field"] == "mass_max"


def test_profile_range_without_signal_422(client):
    # valid range, but no stick's +/-4-sigma support reaches it
    resp = client.post(
        "/api/v1/profile", json=_profile_payload(mass_min=100.0, mass_max=101.0)
    )
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "invalid_profile"
    assert detail["field"] == "mass_max"
    resp = client.post(
        "/api/v1/profile",
        json=_profile_payload(mass_min=500.0, mass_max=600.0, mass_step=0.01),
    )
    assert resp.status_code == 422
    assert resp.json()["detail"]["field"] == "mass_min"


def test_profile_grid_limit_422(client, monkeypatch):
    monkeypatch.setattr(settings, "profile_max_points", 100)
    resp = client.post(
        "/api/v1/profile",
        json=_profile_payload(mass_min=180.0, mass_max=180.1, mass_step=0.001),
    )
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "invalid_profile"
    assert detail["field"] == "mass_step"
    # exactly 100 points is accepted
    resp = client.post(
        "/api/v1/profile",
        json=_profile_payload(mass_min=180.0, mass_max=180.099, mass_step=0.001),
    )
    assert resp.status_code == 200
    assert resp.json()["point_count"] == 100
