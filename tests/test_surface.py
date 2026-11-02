import math

import pytest

from geolab import HAYFORD, WGS84, azimuth_radius, meridian_radius, prime_vertical_radius
from geolab.surface import (
    arc_element,
    area_element,
    curve_azimuth,
    fundamental_forms,
    meusnier_radius,
    normal_section_radius,
)
from geolab.trace import Trace

LATS = [math.radians(d) for d in (-80, -30, 0, 23.5, 41.0083, 75)]


def position(lat: float, lon: float, ell=WGS84) -> tuple[float, float, float]:
    n = ell.a / math.sqrt(1 - ell.e_sq * math.sin(lat) ** 2)
    return (
        n * math.cos(lat) * math.cos(lon),
        n * math.cos(lat) * math.sin(lon),
        n * (1 - ell.e_sq) * math.sin(lat),
    )


def d(f, h=1e-5):
    return lambda x: tuple((p - m) / (2 * h) for p, m in zip(f(x + h), f(x - h)))


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


@pytest.mark.parametrize("lat", LATS)
def test_first_form_matches_numerical_derivatives(lat):
    lon = 0.4
    r_phi = d(lambda p: position(p, lon))(lat)
    r_lam = d(lambda l: position(lat, l))(lon)
    ff = fundamental_forms(lat, WGS84)
    assert dot(r_phi, r_phi) == pytest.approx(ff.E, rel=1e-8)
    assert dot(r_lam, r_lam) == pytest.approx(ff.G, rel=1e-8)
    assert dot(r_phi, r_lam) == pytest.approx(0, abs=1e-6 * ff.E)


@pytest.mark.parametrize("lat", LATS)
def test_second_form_matches_numerical_derivatives(lat):
    lon, h = 0.4, 1e-4
    r_phi = d(lambda p: position(p, lon))(lat)
    r_lam = d(lambda l: position(lat, l))(lon)
    normal = cross(r_lam, r_phi)  # outward
    size = math.sqrt(dot(normal, normal))
    unit = tuple(x / size for x in normal)
    r_pp = tuple(
        (a - 2 * b + c) / h**2
        for a, b, c in zip(position(lat + h, lon), position(lat, lon), position(lat - h, lon))
    )
    r_ll = tuple(
        (a - 2 * b + c) / h**2
        for a, b, c in zip(position(lat, lon + h), position(lat, lon), position(lat, lon - h))
    )
    ff = fundamental_forms(lat, WGS84)
    # Sign convention: positive for the outward-curving ellipsoid.
    assert -dot(r_pp, unit) == pytest.approx(ff.L, rel=1e-5)
    assert -dot(r_ll, unit) == pytest.approx(ff.N, rel=1e-5, abs=1e-3)


@pytest.mark.parametrize("lat", LATS)
def test_relations_to_radii(lat):
    ff = fundamental_forms(lat, HAYFORD)
    m, n = meridian_radius(lat, HAYFORD), prime_vertical_radius(lat, HAYFORD)
    assert ff.E == pytest.approx(m**2)
    assert ff.G == pytest.approx((n * math.cos(lat)) ** 2)
    assert ff.L == pytest.approx(m)
    assert ff.gauss_curvature == pytest.approx(1 / (m * n))
    assert ff.mean_curvature == pytest.approx((1 / m + 1 / n) / 2)
    assert ff.gauss_radius == pytest.approx(math.sqrt(m * n))
    assert m <= ff.gauss_radius <= n


@pytest.mark.parametrize("az", [0, 0.3, math.pi / 4, 1.2, math.pi / 2, 2.5])
def test_normal_section_radius_is_euler(az):
    lat = math.radians(41.0083)
    assert normal_section_radius(lat, az, HAYFORD) == pytest.approx(azimuth_radius(lat, az, HAYFORD))


def test_normal_section_extremes():
    lat = math.radians(41.0083)
    assert normal_section_radius(lat, 0, WGS84) == pytest.approx(meridian_radius(lat, WGS84))
    assert normal_section_radius(lat, math.pi / 2, WGS84) == pytest.approx(prime_vertical_radius(lat, WGS84))


def test_curve_azimuth():
    lat = math.radians(41.0)
    assert curve_azimuth(lat, 1e-6, 0.0, WGS84) == pytest.approx(0.0)
    assert curve_azimuth(lat, 0.0, 1e-6, WGS84) == pytest.approx(math.pi / 2)
    # Equal metric steps north and east → 45°.
    m, n = meridian_radius(lat, WGS84), prime_vertical_radius(lat, WGS84)
    assert curve_azimuth(lat, 1e-6 / m, 1e-6 / (n * math.cos(lat)), WGS84) == pytest.approx(math.pi / 4)


def test_arc_element_matches_vincenty():
    from geolab import inverse

    lat = math.radians(41.0)
    dlat, dlon = 1e-5, 2e-5
    ds = arc_element(lat + dlat / 2, dlat, dlon, WGS84)
    assert ds == pytest.approx(inverse(lat, 0.0, lat + dlat, dlon, WGS84).distance, rel=1e-8)


def test_area_element():
    lat = math.radians(41.0)
    m, n = meridian_radius(lat, WGS84), prime_vertical_radius(lat, WGS84)
    assert area_element(lat, 1e-4, 2e-4, WGS84) == pytest.approx(m * n * math.cos(lat) * 2e-8)


def test_meusnier():
    assert meusnier_radius(6.4e6, 0.0) == 6.4e6
    assert meusnier_radius(6.4e6, math.radians(60)) == pytest.approx(3.2e6)


def test_trace():
    trace = Trace()
    ff = fundamental_forms(0.7, WGS84, trace)
    assert trace["E"] == ff.E and trace["R_G"] == ff.gauss_radius
