"""The series are checked against exact normal-section azimuths, computed from
the planes through each point's ellipsoid normal. They are first-order series:
the next term is of relative size about (S/N) tan φ, a few per cent at 300 km."""

import math

import pytest

from geolab import HAYFORD, WGS84, Geodetic, direct, inverse, parse_angle
from geolab.normal_section import (
    eta_squared,
    normal_section_corrections,
    target_height_correction,
)

ARCSEC = math.radians(1 / 3600)


def xyz(lat, lon, h=0.0, ell=WGS84):
    c = Geodetic(lat, lon, h).to_cartesian(ell)
    return (c.x, c.y, c.z)


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def enu(lat, lon):
    east = (-math.sin(lon), math.cos(lon), 0.0)
    north = (-math.sin(lat) * math.cos(lon), -math.sin(lat) * math.sin(lon), math.cos(lat))
    up = (math.cos(lat) * math.cos(lon), math.cos(lat) * math.sin(lon), math.sin(lat))
    return east, north, up


def azimuth_of(direction, lat, lon):
    east, north, _ = enu(lat, lon)
    return math.atan2(dot(direction, east), dot(direction, north)) % math.tau


def normal_section_azimuth(a, b, target_height=0.0, ell=WGS84):
    """Azimuth at A of the plane through A's normal and B (at height h above B)."""
    return azimuth_of(sub(xyz(*b, target_height, ell=ell), xyz(*a, ell=ell)), *a)


def reverse_section_azimuth(a, b, ell=WGS84):
    """Azimuth at A of the plane through B's normal and A."""
    pa, pb = xyz(*a, ell=ell), xyz(*b, ell=ell)
    plane_normal = cross(enu(*b)[2], sub(pa, pb))
    tangent = cross(enu(*a)[2], plane_normal)
    if dot(tangent, sub(pb, pa)) < 0:
        tangent = tuple(-x for x in tangent)
    return azimuth_of(tangent, *a)


def line(lat, az, dist, ell=WGS84):
    a = (math.radians(lat), 0.4)
    r = direct(a[0], a[1], math.radians(az), dist, ell)
    return a, (r.lat2, r.lon2)


@pytest.mark.parametrize("lat, az, dist", [(0, 45, 1e5), (40, 45, 2e5), (50, 30, 3e5), (60, 120, 2e5)])
def test_geodesic_minus_normal_section(lat, az, dist):
    a, b = line(lat, az, dist)
    exact = inverse(*a, *b, WGS84).azimuth1 - normal_section_azimuth(a, b)
    alpha_n = normal_section_azimuth(a, b)
    series = normal_section_corrections(a[0], alpha_n, dist, WGS84).geodesic_minus_normal
    assert series == pytest.approx(exact, rel=0.05)


@pytest.mark.parametrize("lat, az, dist", [(0, 45, 1e5), (40, 45, 2e5), (50, 30, 3e5)])
def test_normal_minus_reverse_section(lat, az, dist):
    a, b = line(lat, az, dist)
    exact = normal_section_azimuth(a, b) - reverse_section_azimuth(a, b)
    series = normal_section_corrections(a[0], normal_section_azimuth(a, b), dist, WGS84).normal_minus_reverse
    assert series == pytest.approx(exact, rel=0.05)


def test_geodesic_splits_the_sections_one_to_two():
    c = normal_section_corrections(math.radians(40), math.radians(45), 2e5, WGS84)
    assert -c.geodesic_minus_normal / c.normal_minus_reverse == pytest.approx(1 / 3)


def test_slide_magnitudes():
    # 100 km, α = 45°, near the equator: α − α′ ≈ −0.028″ and d_max ≈ 5.2 mm
    # (the slides round these).
    c = normal_section_corrections(0.0, math.radians(45), 1e5, HAYFORD)
    assert c.geodesic_minus_normal / ARCSEC == pytest.approx(-0.028, abs=1e-3)
    assert c.max_separation == pytest.approx(5.2e-3, abs=1e-4)
    # s = 100 km, α = 45°, φ = 50°: s − s′ ≈ −1e-10 m.
    c = normal_section_corrections(math.radians(50), math.radians(45), 1e5, HAYFORD)
    assert c.length_difference == pytest.approx(-1e-10, rel=0.5)
    assert c.length_difference < 0


def test_on_a_meridian_there_is_no_difference():
    c = normal_section_corrections(math.radians(40), 0.0, 3e5, WGS84)
    assert c.geodesic_minus_normal == 0 and c.normal_minus_reverse == 0


def test_target_height_slide_example():
    # φ₁ = 50°, α′ = 45°, h = 3000 m → 0.134″ on the slide (0.135″ unrounded)
    corr = target_height_correction(math.radians(50), math.radians(45), 3000.0, HAYFORD)
    assert corr / ARCSEC == pytest.approx(0.134, abs=2e-3)


@pytest.mark.parametrize("lat, az", [(50, 45), (30, 120), (65, 20)])
def test_target_height_matches_geometry(lat, az):
    a, b = line(lat, az, 5e4)
    h = 3000.0
    exact = normal_section_azimuth(a, b) - normal_section_azimuth(a, b, target_height=h)
    series = target_height_correction(a[0], normal_section_azimuth(a, b), h, WGS84)
    assert series == pytest.approx(exact, rel=0.02, abs=1e-4 * ARCSEC)


def test_eta_squared():
    assert eta_squared(0.0, WGS84) == pytest.approx(WGS84.ep_sq)
    assert eta_squared(math.pi / 2, WGS84) == pytest.approx(0.0, abs=1e-20)
