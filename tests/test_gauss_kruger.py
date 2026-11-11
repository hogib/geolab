"""Gauss–Krüger projection: the course series against Krüger's n⁶ series,
which is itself checked against the meridian arc and by finite differences."""

import math

import pytest

from geolab import HAYFORD, WGS84, direct
from geolab.arcs import meridian_distance
from geolab.curvature import meridian_radius, prime_vertical_radius
from geolab.gauss_kruger import (
    gauss_kruger_from_geodetic,
    geodetic_from_gauss_kruger,
    gk_direct,
    gk_inverse,
    krueger_forward,
    krueger_inverse,
)
from geolab.trace import Trace

ARCSEC = math.radians(1 / 3600)
L0 = math.radians(30)


def rad(d, m=0, s=0):
    return math.radians(d + m / 60 + s / 3600)


class TestKrueger:
    def test_central_meridian_is_the_meridian_arc(self):
        for deg in (0, 23, 41, 70):
            p = krueger_forward(rad(deg), L0, L0, WGS84)
            assert p.x == pytest.approx(meridian_distance(rad(deg), WGS84), abs=1e-6)
            assert p.y == 0 and p.scale == pytest.approx(1)

    @pytest.mark.parametrize("lat, dlon", [(41, 3), (-35, -6), (10, 9), (65, 2)])
    def test_round_trip(self, lat, dlon):
        p = krueger_forward(rad(lat), L0 + rad(dlon), L0, WGS84)
        back = krueger_inverse(p.x, p.y, L0, WGS84)
        assert back == pytest.approx((rad(lat), L0 + rad(dlon)), abs=1e-6 * ARCSEC)

    def test_scale_and_convergence_by_finite_differences(self):
        lat, lon, h = rad(41), rad(33), 1e-7
        a = krueger_forward(lat, lon, L0, WGS84)
        east = krueger_forward(lat, lon + h, L0, WGS84)
        north = krueger_forward(lat + h, lon, L0, WGS84)
        along_parallel = math.hypot(east.x - a.x, east.y - a.y) / (prime_vertical_radius(lat, WGS84) * math.cos(lat) * h)
        along_meridian = math.hypot(north.x - a.x, north.y - a.y) / (meridian_radius(lat, WGS84) * h)
        assert along_parallel == pytest.approx(a.scale, rel=1e-8)
        assert along_meridian == pytest.approx(a.scale, rel=1e-8)  # conformal
        # The meridian's image points c west of grid north (east of λ₀).
        assert -math.atan2(north.y - a.y, north.x - a.x) == pytest.approx(a.convergence, abs=1e-3 * ARCSEC)


POINTS = [(41, 0.5), (41, 1.5), (36, -1.5), (42, 2.0), (-20, 1.0)]


@pytest.mark.parametrize("lat, dlon", POINTS)
def test_forward_series(lat, dlon):
    k = krueger_forward(rad(lat), L0 + rad(dlon), L0, WGS84)
    s = gauss_kruger_from_geodetic(rad(lat), L0 + rad(dlon), L0, WGS84)
    assert (s.x, s.y) == pytest.approx((k.x, k.y), abs=1e-4)
    assert s.scale == pytest.approx(k.scale, abs=1e-10)
    assert s.convergence == pytest.approx(k.convergence, abs=1e-3 * ARCSEC)


@pytest.mark.parametrize("lat, dlon", POINTS)
def test_inverse_series(lat, dlon):
    k = krueger_forward(rad(lat), L0 + rad(dlon), L0, WGS84)
    g = geodetic_from_gauss_kruger(k.x, k.y, L0, WGS84)
    assert (g.lat, g.lon) == pytest.approx((rad(lat), L0 + rad(dlon)), abs=1e-4 * ARCSEC)
    assert g.scale == pytest.approx(k.scale, abs=1e-9)
    assert g.convergence == pytest.approx(k.convergence, abs=0.01 * ARCSEC)
    assert abs(g.foot_lat) > abs(g.lat)


def test_six_degree_zone_edge():
    # 3° from λ₀ the truncated series is still within a millimetre.
    k = krueger_forward(rad(39), L0 + rad(3), L0, WGS84)
    s = gauss_kruger_from_geodetic(rad(39), L0 + rad(3), L0, WGS84)
    assert (s.x, s.y) == pytest.approx((k.x, k.y), abs=1e-3)


def test_expanded_about_another_latitude():
    # The general form ΔX + iY = Σ aₖ(Δq + iΔλ)ᵏ about φ₀ ≠ φ.
    k = krueger_forward(rad(41), L0 + rad(1), L0, WGS84)
    s = gauss_kruger_from_geodetic(rad(41), L0 + rad(1), L0, WGS84, lat0=rad(40.5))
    assert (s.x, s.y) == pytest.approx((k.x, k.y), abs=1e-3)


def test_course_sheet_corners():
    # 1/5000 sheet G25-d-01-d, Hayford, 3° zone λ₀ = 30°, Y with 500 000 m added
    # (course, week 12). The printout's X is a systematic ~8 mm higher.
    corners = {
        (rad(40, 42), rad(30, 30)): (542259.752, 4507456.400),
        (rad(40, 43, 30), rad(30, 30)): (542243.948, 4510232.706),
        (rad(40, 43, 30), rad(30, 31, 30)): (544356.154, 4510245.033),
        (rad(40, 42), rad(30, 31, 30)): (544372.748, 4507468.725),
    }
    for (lat, lon), (y, x) in corners.items():
        p = gauss_kruger_from_geodetic(lat, lon, L0, HAYFORD)
        assert p.y + 500000 == pytest.approx(y, abs=3e-3)
        assert p.x == pytest.approx(x, abs=1e-2)


def exact_line(x1, y1, bearing, s):
    lat, lon = krueger_inverse(x1, y1, L0, WGS84)
    c1 = krueger_forward(lat, lon, L0, WGS84).convergence
    r = direct(lat, lon, bearing + c1, s, WGS84)
    p = krueger_forward(r.lat2, r.lon2, L0, WGS84)
    return p.x, p.y, (r.back_azimuth - p.convergence) % math.tau


LINES = [(4.5e6, 2e4, 0.6, 3e4), (4.5e6, -1e5, 2.2, 5e4), (4.6e6, 1.2e5, 4.0, 4e4)]


@pytest.mark.parametrize("x1, y1, bearing, s", LINES)
def test_direct_problem(x1, y1, bearing, s):
    x2, y2, t21 = exact_line(x1, y1, bearing, s)
    r = gk_direct(x1, y1, bearing, s, WGS84)
    assert (r.x2, r.y2) == pytest.approx((x2, y2), abs=1e-3)
    assert math.remainder(r.bearing2 - t21, math.tau) == pytest.approx(0, abs=0.005 * ARCSEC)


@pytest.mark.parametrize("x1, y1, bearing, s", LINES)
def test_inverse_problem(x1, y1, bearing, s):
    x2, y2, t21 = exact_line(x1, y1, bearing, s)
    r = gk_inverse(x1, y1, x2, y2, WGS84)
    assert r.distance == pytest.approx(s, abs=2e-4)
    assert math.remainder(r.bearing1 - bearing, math.tau) == pytest.approx(0, abs=0.002 * ARCSEC)
    assert math.remainder(r.bearing2 - t21, math.tau) == pytest.approx(0, abs=0.002 * ARCSEC)


def test_trace():
    trace = Trace()
    p = gauss_kruger_from_geodetic(0.71, 0.55, L0, WGS84, trace=trace)
    assert trace["Y"] == p.y and trace["a1"] > 0 and trace["a2"] < 0
