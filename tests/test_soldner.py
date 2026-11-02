"""Soldner coordinates, checked against an exact construction with Vincenty:
the foot point is found so that the geodesic leaving it due east (or west)
passes through the point."""

import math

import pytest

from geolab import HAYFORD, WGS84, direct, inverse
from geolab.arcs import meridian_distance
from geolab.soldner import (
    geodetic_from_soldner,
    soldner_direct,
    soldner_from_geodetic,
    soldner_inverse,
)
from geolab.trace import Trace

ARCSEC = math.radians(1 / 3600)
LAT0, LON0 = math.radians(40.0), math.radians(30.0)


def exact_soldner(lat, lon, ell=WGS84, lat0=LAT0, lon0=LON0):
    lo, hi = lat - 0.2, lat + 0.2
    for _ in range(80):
        mid = (lo + hi) / 2
        az = inverse(mid, lon0, lat, lon, ell).azimuth1
        lo, hi = (mid, hi) if math.cos(az) > 0 else (lo, mid)
    r = inverse(lo, lon0, lat, lon, ell)
    east = math.sin(r.azimuth1) > 0
    x = meridian_distance(lo, ell) - meridian_distance(lat0, ell)
    gamma = math.remainder(r.azimuth2 - (math.pi / 2 if east else 3 * math.pi / 2), math.tau)
    return x, r.distance if east else -r.distance, gamma


POINTS = [(0.1, 0.3), (0.5, 1.0), (-0.4, -1.2), (0.8, 1.5)]


@pytest.mark.parametrize("dlat, dlon", POINTS)
def test_forward_matches_exact(dlat, dlon):
    lat, lon = LAT0 + math.radians(dlat), LON0 + math.radians(dlon)
    x, y, gamma = exact_soldner(lat, lon)
    s = soldner_from_geodetic(lat, lon, LAT0, LON0, WGS84)
    assert s.x == pytest.approx(x, abs=1e-3)
    assert s.y == pytest.approx(y, abs=1e-3)
    assert s.convergence == pytest.approx(gamma, abs=1e-3 * ARCSEC)


@pytest.mark.parametrize("dlat, dlon", POINTS)
def test_inverse_matches_exact(dlat, dlon):
    lat, lon = LAT0 + math.radians(dlat), LON0 + math.radians(dlon)
    x, y, gamma = exact_soldner(lat, lon)
    g = geodetic_from_soldner(x, y, LAT0, LON0, WGS84)
    assert g.lat == pytest.approx(lat, abs=1e-4 * ARCSEC)
    assert g.lon == pytest.approx(lon, abs=1e-4 * ARCSEC)
    assert g.convergence == pytest.approx(gamma, abs=1e-3 * ARCSEC)
    assert g.foot_lat > lat  # northern hemisphere: φ_f > φ


def test_on_the_central_meridian():
    s = soldner_from_geodetic(math.radians(41), LON0, LAT0, LON0, HAYFORD)
    assert s.y == 0 and s.convergence == 0
    assert s.x == pytest.approx(meridian_distance(math.radians(41), HAYFORD) - meridian_distance(LAT0, HAYFORD))


def exact_line(p1, bearing, distance):
    """Direct problem done exactly: Soldner → geodetic → Vincenty → Soldner."""
    x1, y1 = p1
    g1 = geodetic_from_soldner(x1, y1, LAT0, LON0, WGS84)
    _, _, gamma1 = exact_soldner(g1.lat, g1.lon)
    r = direct(g1.lat, g1.lon, bearing + gamma1, distance, WGS84)
    x2, y2, gamma2 = exact_soldner(r.lat2, r.lon2)
    return x2, y2, (r.back_azimuth - gamma2) % math.tau


@pytest.mark.parametrize("p1, bearing, s", [((1e4, 2e4), 0.6, 3e4), ((-5e4, -4e4), 2.2, 5e4), ((8e4, 6e4), 4.0, 4e4)])
def test_direct_problem(p1, bearing, s):
    x2, y2, t21 = exact_line(p1, bearing, s)
    r = soldner_direct(*p1, bearing, s, LAT0, WGS84)
    assert r.x2 == pytest.approx(x2, abs=5e-3)
    assert r.y2 == pytest.approx(y2, abs=5e-3)
    assert math.remainder(r.bearing2 - t21, math.tau) == pytest.approx(0, abs=0.01 * ARCSEC)


@pytest.mark.parametrize("p1, bearing, s", [((1e4, 2e4), 0.6, 3e4), ((-5e4, -4e4), 2.2, 5e4), ((8e4, 6e4), 4.0, 4e4)])
def test_inverse_problem(p1, bearing, s):
    x2, y2, t21 = exact_line(p1, bearing, s)
    r = soldner_inverse(*p1, x2, y2, LAT0, WGS84)
    assert r.distance == pytest.approx(s, abs=5e-3)
    assert math.remainder(r.bearing1 - bearing, math.tau) == pytest.approx(0, abs=0.01 * ARCSEC)
    assert math.remainder(r.bearing2 - t21, math.tau) == pytest.approx(0, abs=0.01 * ARCSEC)


def test_round_trip():
    s = soldner_from_geodetic(math.radians(40.7), math.radians(31.1), LAT0, LON0, WGS84)
    g = geodetic_from_soldner(s.x, s.y, LAT0, LON0, WGS84)
    assert (g.lat, g.lon) == pytest.approx((math.radians(40.7), math.radians(31.1)), abs=1e-4 * ARCSEC)


def test_trace():
    trace = Trace()
    s = soldner_from_geodetic(0.71, 0.53, LAT0, LON0, WGS84, trace)
    assert trace["y"] == s.y
