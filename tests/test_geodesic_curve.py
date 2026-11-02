import math

import pytest

from geolab import HAYFORD, WGS84, direct, parse_angle
from geolab.geodesic_curve import (
    azimuths_at_latitude,
    clairaut_constant,
    integrate_geodesic,
    vertex_latitude,
)

ARCSEC = math.radians(1 / 3600)


class TestCourseExample:
    """Slides: φ = 45°, α = 45° on the international (Hayford) ellipsoid."""

    k = clairaut_constant(math.radians(45), math.radians(45), HAYFORD)

    def test_constant(self):
        assert self.k == pytest.approx(3194567.525, abs=1e-3)

    def test_vertex_latitude(self):
        assert vertex_latitude(self.k, HAYFORD) == pytest.approx(parse_angle("60 01 40.61"), abs=0.01 * ARCSEC)

    def test_azimuths_at_equator(self):
        a1, a2 = azimuths_at_latitude(self.k, 0.0, HAYFORD)
        assert a1 == pytest.approx(parse_angle("30 03 20.71"), abs=0.01 * ARCSEC)
        assert a2 == pytest.approx(parse_angle("149 56 39.29"), abs=0.01 * ARCSEC)


@pytest.mark.parametrize("az_deg", [10, 45, 100, 170, 260])
def test_clairaut_constant_is_conserved_along_vincenty(az_deg):
    lat1, az1 = math.radians(37.0), math.radians(az_deg)
    k1 = clairaut_constant(lat1, az1, WGS84)
    for s in (1e5, 1e6, 5e6):
        r = direct(lat1, 0.2, az1, s, WGS84)
        assert clairaut_constant(r.lat2, r.azimuth2, WGS84) == pytest.approx(k1, abs=1e-5)


def test_vertex_is_reached_by_vincenty():
    lat1, az1 = math.radians(20.0), math.radians(60.0)
    k = clairaut_constant(lat1, az1, WGS84)
    top = max(direct(lat1, 0.0, az1, s * 1e4, WGS84).lat2 for s in range(0, 2000))
    assert top == pytest.approx(vertex_latitude(k, WGS84), abs=1e-7)


def test_unreachable_latitude():
    k = clairaut_constant(math.radians(45), math.radians(45), HAYFORD)
    with pytest.raises(ValueError):
        azimuths_at_latitude(k, math.radians(70), HAYFORD)


@pytest.mark.parametrize("az_deg, dist", [(0, 3e5), (37, 5e5), (90, 1e6), (215, 2e6)])
def test_integrated_geodesic_matches_vincenty(az_deg, dist):
    lat1, lon1, az1 = math.radians(41.0), math.radians(29.0), math.radians(az_deg)
    got = integrate_geodesic(lat1, lon1, az1, dist, WGS84, steps=400)
    want = direct(lat1, lon1, az1, dist, WGS84)
    m_per_rad = WGS84.a
    assert (got.lat - want.lat2) * m_per_rad == pytest.approx(0, abs=1e-4)
    assert (got.lon - want.lon2) * m_per_rad == pytest.approx(0, abs=1e-4)
    assert got.azimuth == pytest.approx(want.azimuth2, abs=1e-9)
