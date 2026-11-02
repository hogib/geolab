"""Worked examples from the course (Hayford ellipsoid)."""

import pytest

from geolab import HAYFORD, Geodetic, LatType, Latitude, meridian_ellipse_point, parse_angle
from geolab.point import Cartesian3D, GeodeticMethod

ARCSEC = parse_angle("0 0 1")


def test_coordinates_from_geodetic_latitude():
    pt = meridian_ellipse_point(parse_angle("35 20 40"), LatType.GEODETIC, HAYFORD)
    assert pt.w == pytest.approx(0.9988744915, abs=1e-10)
    assert pt.v == pytest.approx(1.002249067, abs=1e-9)
    assert pt.p == pytest.approx(5208643.936, abs=2e-3)
    assert pt.z == pytest.approx(3669161.169, abs=2e-3)


def test_c_over_v_form_agrees():
    pt = meridian_ellipse_point(parse_angle("35 20 40"), LatType.GEODETIC, HAYFORD)
    c = HAYFORD.a**2 / HAYFORD.b
    assert c == pytest.approx(6399936.608, abs=1e-3)
    assert pt.p == pytest.approx(c * __import__("math").cos(pt.lat) / pt.v, abs=1e-6)


def test_from_reduced_latitude():
    pt = meridian_ellipse_point(parse_angle("35 50 37.95"), LatType.PARAMETRIC, HAYFORD)
    assert pt.lat == pytest.approx(parse_angle("35 56 8.3426"), abs=1e-4 * ARCSEC)
    assert pt.p == pytest.approx(5170421.085, abs=2e-3)
    assert pt.z == pytest.approx(3722471.510, abs=2e-3)


def test_geocentric_latitude_and_radius():
    pt = meridian_ellipse_point(parse_angle("35 56 8.34"), LatType.GEODETIC, HAYFORD)
    assert pt.geocentric_lat == pytest.approx(parse_angle("35 45 7.9046"), abs=1e-4 * ARCSEC)
    assert pt.w == pytest.approx(0.9988416055, abs=1e-10)
    assert pt.r == pytest.approx(6371031.968, abs=2e-3)


def test_reduced_from_geocentric():
    beta = Latitude(parse_angle("35 45 7.9"), LatType.GEOCENTRIC).convert_to(LatType.PARAMETRIC, HAYFORD)
    assert beta.value == pytest.approx(parse_angle("35 50 37.9427"), abs=1e-4 * ARCSEC)


def test_from_geocentric_round_trip():
    pt = meridian_ellipse_point(parse_angle("41 00 30"), LatType.GEODETIC, HAYFORD)
    back = meridian_ellipse_point(pt.geocentric_lat, LatType.GEOCENTRIC, HAYFORD)
    assert (back.p, back.z, back.r) == pytest.approx((pt.p, pt.z, pt.r), abs=1e-6)


EXAMPLE = Geodetic(parse_angle("36 5 32.1"), parse_angle("27 6 56.988"), 1250.0)
EXAMPLE_XYZ = (4593929.069, 2352434.103, 3737263.328)


def test_geodetic_to_cartesian_example():
    xyz = EXAMPLE.to_cartesian(HAYFORD)
    assert (xyz.x, xyz.y, xyz.z) == pytest.approx(EXAMPLE_XYZ, abs=2e-3)


@pytest.mark.parametrize("method", list(GeodeticMethod))
def test_cartesian_to_geodetic_example(method):
    # The second example inverts the first; the "2. Yol" closed form is Bowring's method.
    g = Cartesian3D(*EXAMPLE_XYZ).to_geodetic(HAYFORD, method)
    assert g.lat == pytest.approx(EXAMPLE.lat, abs=1e-4 * ARCSEC)
    assert g.lon == pytest.approx(EXAMPLE.lon, abs=1e-4 * ARCSEC)
    assert g.height == pytest.approx(1250.0, abs=2e-3)
