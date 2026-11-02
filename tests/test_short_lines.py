"""Legendre series and Gauss mid-latitude methods: course examples (Bessel
ellipsoid) and agreement with Vincenty."""

import math

import pytest

from geolab import BESSEL_1841, WGS84, direct, inverse, parse_angle
from geolab.legendre_series import legendre_direct, legendre_inverse
from geolab.mid_latitude import brackets, mid_latitude_direct, mid_latitude_inverse
from geolab.trace import Trace

ARCSEC = math.radians(1 / 3600)
GON = math.pi / 200
B = BESSEL_1841


def test_bessel_constants_match_the_slides():
    assert B.a**2 / B.b == pytest.approx(6398786.849, abs=1e-3)
    assert B.ep_sq == pytest.approx(0.0067192188, abs=1e-10)


class TestLegendreCourseExample:
    lat1 = parse_angle("48 26 45.4355")
    lon1 = parse_angle("10 42 59.3215")
    az = 35.658346 * GON
    r = legendre_direct(lat1, lon1, az, 47835.627, B)

    def test_u_v(self):
        assert self.r.u == pytest.approx(40525.930, abs=1e-3)
        assert self.r.v == pytest.approx(25414.488, abs=1e-3)

    @pytest.mark.parametrize(
        "i, expected",
        [(0, 1312.15037), (1, -1.84626), (2, -0.04164), (3, -0.01656), (4, 0.00002), (6, -0.00014)],
    )
    def test_latitude_terms(self, i, expected):
        assert self.r.dlat_terms[i] / ARCSEC == pytest.approx(expected, abs=1e-5)

    @pytest.mark.parametrize(
        "i, expected",
        [(0, 1236.86620), (1, 8.85045), (2, -0.00830), (3, 0.07997), (4, -0.00022), (5, 0.00069), (8, 0.00001)],
    )
    def test_longitude_terms(self, i, expected):
        assert self.r.dlon_terms[i] / ARCSEC == pytest.approx(expected, abs=1e-5)

    @pytest.mark.parametrize(
        "i, expected",
        [(0, 925.58453), (1, 9.23269), (2, -0.00865), (3, 0.07842), (4, -0.00022), (5, 0.00069)],
    )
    def test_azimuth_terms(self, i, expected):
        assert self.r.dazimuth_terms[i] / ARCSEC == pytest.approx(expected, abs=1e-5)

    def test_result(self):
        assert self.r.lat2 == pytest.approx(parse_angle("48 48 35.6813"), abs=1e-4 * ARCSEC)
        assert self.r.lon2 == pytest.approx(parse_angle("11 03 45.1103"), abs=1e-4 * ARCSEC)
        # The slide prints 235.9468292 gon, an addition slip: 235.658346 + 0.288546 = 235.946892.
        assert self.r.back_azimuth / GON == pytest.approx(235.946892, abs=1e-6)

    def test_matches_vincenty(self):
        v = direct(self.lat1, self.lon1, self.az, 47835.627, B)
        assert self.r.lat2 == pytest.approx(v.lat2, abs=1e-5 * ARCSEC)
        assert self.r.lon2 == pytest.approx(v.lon2, abs=1e-5 * ARCSEC)
        assert self.r.back_azimuth == pytest.approx(v.back_azimuth, abs=1e-4 * ARCSEC)


P1 = (parse_angle("53 50 02.8809"), parse_angle("10 12 04.1772"))
P2 = (parse_angle("54 13 15.2891"), parse_angle("10 30 47.2430"))


class TestMidLatitudeCourseExamples:
    def test_brackets(self):
        k = brackets((P1[0] + P2[0]) / 2, B)
        assert k.m * ARCSEC == pytest.approx(30.9146247, abs=1e-7)  # 1/[1], m per ″
        assert k.n * ARCSEC == pytest.approx(30.9862961, abs=1e-7)  # 1/[2]
        rho2 = 1 / ARCSEC**2
        assert k.b4 / rho2 * 1e8 == pytest.approx(0.00009385, abs=1e-8)
        assert k.b5 / rho2 * 1e8 == pytest.approx(0.00009748, abs=1e-8)
        assert k.b6 / rho2 * 1e8 == pytest.approx(-0.00000059, abs=3e-8)  # tiny; effect ~1e-8 relative
        assert k.b7 / rho2 * 1e8 == pytest.approx(0.00019633, abs=1e-8)
        assert k.b8 / rho2 * 1e8 == pytest.approx(0.00029426, abs=1e-8)

    def test_inverse_example_1(self):
        r = mid_latitude_inverse(*P1, *P2, B)
        assert r.distance == pytest.approx(47652.596, abs=2e-3)
        assert r.azimuth == pytest.approx(parse_angle("25 16 31.98"), abs=0.006 * ARCSEC)
        assert r.back_azimuth == pytest.approx(parse_angle("205 31 40.88"), abs=0.006 * ARCSEC)
        assert r.delta_azimuth / ARCSEC == pytest.approx(908.90, abs=0.005)

    def test_direct_example_2(self):
        r = mid_latitude_direct(P1[0], P1[1], parse_angle("25 16 31.98"), 47652.596, B)
        assert r.lat2 == pytest.approx(P2[0], abs=1e-3 * ARCSEC)
        assert r.lon2 == pytest.approx(P2[1], abs=1e-3 * ARCSEC)
        assert r.back_azimuth == pytest.approx(parse_angle("205 31 40.88"), abs=0.006 * ARCSEC)

    def test_direct_first_approximation(self):
        trace = Trace()
        mid_latitude_direct(P1[0], P1[1], parse_angle("25 16 31.98"), 47652.596, B, trace)
        first = [e for e in trace.entries if hasattr(e, "value") and e.name == "Δλ"][0]
        assert first.value / ARCSEC == pytest.approx(1112.6975, abs=1e-3)

    def test_inverse_matches_vincenty(self):
        r = mid_latitude_inverse(*P1, *P2, B)
        v = inverse(*P1, *P2, B)
        assert r.distance == pytest.approx(v.distance, abs=1e-3)
        assert r.azimuth == pytest.approx(v.azimuth1, abs=1e-3 * ARCSEC)
        assert r.back_azimuth == pytest.approx(v.back_azimuth, abs=1e-3 * ARCSEC)


LINES = [(41.0, 30.0, 3e4), (48.0, 100.0, 6e4), (55.0, 215.0, 1e5), (-35.0, 320.0, 8e4), (10.0, 75.0, 1.2e5)]


@pytest.mark.parametrize("lat, az, s", LINES)
def test_round_trips_and_vincenty(lat, az, s):
    lat1, lon1, az1 = math.radians(lat), math.radians(29.0), math.radians(az)
    v = direct(lat1, lon1, az1, s, WGS84)

    leg = legendre_direct(lat1, lon1, az1, s, WGS84)
    assert leg.lat2 == pytest.approx(v.lat2, abs=1e-3 * ARCSEC)
    assert leg.lon2 == pytest.approx(v.lon2, abs=1e-3 * ARCSEC)
    back = legendre_inverse(lat1, lon1, v.lat2, v.lon2, WGS84)
    assert back.distance == pytest.approx(s, abs=1e-3)
    assert back.azimuth == pytest.approx(az1, abs=1e-3 * ARCSEC)

    mid = mid_latitude_direct(lat1, lon1, az1, s, WGS84)
    assert mid.lat2 == pytest.approx(v.lat2, abs=1e-3 * ARCSEC)
    assert mid.lon2 == pytest.approx(v.lon2, abs=1e-3 * ARCSEC)
    assert mid.back_azimuth == pytest.approx(v.back_azimuth, abs=5e-3 * ARCSEC)
    inv = mid_latitude_inverse(lat1, lon1, v.lat2, v.lon2, WGS84)
    assert inv.distance == pytest.approx(s, abs=1e-3)
    assert inv.azimuth == pytest.approx(az1, abs=1e-3 * ARCSEC)
