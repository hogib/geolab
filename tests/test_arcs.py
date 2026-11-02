import math

import pytest

from geolab.arcs import (
    latitude_from_meridian_distance_series,
    meridian_distance_series,
    meridian_series,
    short_meridian_arc,
)
from geolab import (
    GRS80,
    HAYFORD,
    WGS84,
    Ellipsoid,
    Geodetic,
    Trace,
    direct,
    inverse,
    latitude_from_meridian_arc,
    latitude_from_meridian_distance,
    longitude_from_parallel_arc,
    meridian_arc_length,
    meridian_distance,
    meridian_radius,
    parallel_arc_length,
    parallel_radius,
)

SPHERE = Ellipsoid("Sphere", 6371000.0, 0)


def integrate_meridian(lat1: float, lat2: float, ellipsoid: Ellipsoid, steps: int = 2000) -> float:
    """Simpson's rule on ∫ M dφ."""
    h = (lat2 - lat1) / steps
    total = meridian_radius(lat1, ellipsoid) + meridian_radius(lat2, ellipsoid)
    for i in range(1, steps):
        total += (4 if i % 2 else 2) * meridian_radius(lat1 + i * h, ellipsoid)
    return total * h / 3


class TestMeridianDistance:
    def test_zero_at_equator(self):
        assert meridian_distance(0.0, WGS84) == 0.0

    def test_quarter_meridian_wgs84(self):
        assert meridian_distance(math.pi / 2, WGS84) == pytest.approx(10001965.7293, abs=1e-3)

    def test_quarter_meridian_grs80(self):
        assert meridian_distance(math.pi / 2, GRS80) == pytest.approx(10001965.7293, abs=1e-3)

    def test_odd(self):
        assert meridian_distance(-0.6, WGS84) == -meridian_distance(0.6, WGS84)

    @pytest.mark.parametrize("ellipsoid", [WGS84, HAYFORD])
    @pytest.mark.parametrize("deg", [1, 10, 41.0083, 60, 89.9, 90])
    def test_matches_integral_of_meridian_radius(self, ellipsoid, deg):
        lat = math.radians(deg)
        assert meridian_distance(lat, ellipsoid) == pytest.approx(
            integrate_meridian(0.0, lat, ellipsoid), abs=1e-6
        )

    def test_sphere(self):
        assert meridian_distance(0.7, SPHERE) == pytest.approx(SPHERE.a * 0.7, rel=1e-15)


class TestMeridianArcLength:
    def test_matches_vincenty_along_a_meridian(self):
        lat1, lat2 = math.radians(41.0083), math.radians(-12.5)
        expected = inverse(lat1, 0.5, lat2, 0.5, WGS84).distance
        assert meridian_arc_length(lat1, lat2, WGS84) == pytest.approx(expected, abs=1e-4)

    def test_order_does_not_matter(self):
        assert meridian_arc_length(0.3, 0.9, WGS84) == meridian_arc_length(0.9, 0.3, WGS84)

    def test_trace(self):
        trace = Trace()
        s = meridian_arc_length(0.3, 0.9, WGS84, trace)
        assert trace["Δm"] == s


class TestParallelArcLength:
    def test_equator_matches_vincenty(self):
        expected = inverse(0.0, 0.0, 0.0, 0.5, WGS84).distance
        assert parallel_arc_length(0.0, 0.0, 0.5, WGS84) == pytest.approx(expected, abs=1e-6)

    def test_is_r_times_dlon(self):
        lat = math.radians(41)
        assert parallel_arc_length(lat, 0.2, 0.5, WGS84) == pytest.approx(
            parallel_radius(lat, WGS84) * 0.3
        )

    def test_takes_shorter_way_across_antimeridian(self):
        lon1, lon2 = math.radians(179), math.radians(-179)
        assert parallel_arc_length(0.0, lon1, lon2, WGS84) == pytest.approx(WGS84.a * math.radians(2))

    def test_zero_at_pole(self):
        assert parallel_arc_length(math.pi / 2, 0.0, 1.0, WGS84) == pytest.approx(0.0, abs=1e-9)


class TestLatitudeFromMeridianDistance:
    @pytest.mark.parametrize("ellipsoid", [WGS84, HAYFORD, SPHERE])
    @pytest.mark.parametrize("deg", [-90, -45, -0.001, 0, 10, 41.0083, 60, 89.9999, 90])
    def test_round_trip(self, ellipsoid, deg):
        lat = math.radians(deg)
        p = latitude_from_meridian_distance(meridian_distance(lat, ellipsoid), ellipsoid)
        assert p.lat == pytest.approx(lat, abs=1e-14)

    def test_coordinates_match_ecef(self):
        lat = math.radians(41.0083)
        p = latitude_from_meridian_distance(meridian_distance(lat, HAYFORD), HAYFORD)
        xyz = Geodetic(lat, 0.0, 0.0).to_cartesian(HAYFORD)
        assert p.x == pytest.approx(xyz.x, abs=1e-6)
        assert p.y == pytest.approx(xyz.z, abs=1e-6)

    def test_point_lies_on_meridian_ellipse(self):
        p = latitude_from_meridian_distance(5_000_000.0, HAYFORD)
        assert (p.x / HAYFORD.a) ** 2 + (p.y / HAYFORD.b) ** 2 == pytest.approx(1, abs=1e-15)

    def test_quarter_meridian_is_the_pole(self):
        p = latitude_from_meridian_distance(meridian_distance(math.pi / 2, WGS84), WGS84)
        assert p.lat == pytest.approx(math.pi / 2)
        assert p.x == pytest.approx(0, abs=1e-6)
        assert p.y == pytest.approx(WGS84.b)

    def test_latitude_types(self):
        p = latitude_from_meridian_distance(4_000_000.0, WGS84)
        assert math.tan(p.reduced_lat) == pytest.approx((1 - WGS84.f) * math.tan(p.lat))
        assert math.tan(p.geocentric_lat) == pytest.approx((1 - WGS84.e_sq) * math.tan(p.lat))

    def test_beyond_the_pole(self):
        with pytest.raises(ValueError, match="quarter meridian"):
            latitude_from_meridian_distance(10_002_000.0, WGS84)

    def test_trace(self):
        trace = Trace()
        p = latitude_from_meridian_distance(4_000_000.0, WGS84, trace)
        assert trace["x"] == p.x
        assert trace["φ"] == p.lat


class TestLatitudeFromMeridianArc:
    def test_matches_meridian_arc_length(self):
        lat1 = math.radians(41.0083)
        p = latitude_from_meridian_arc(lat1, -120_000.0, HAYFORD)
        assert p.lat < lat1
        assert meridian_arc_length(lat1, p.lat, HAYFORD) == pytest.approx(120_000.0, abs=1e-6)

    def test_matches_vincenty_direct_due_north(self):
        lat1 = math.radians(39.9)
        p = latitude_from_meridian_arc(lat1, 250_000.0, WGS84)
        r = direct(lat1, 0.3, 0.0, 250_000.0, WGS84)
        assert p.lat == pytest.approx(r.lat2, abs=1e-11)

    def test_crossing_the_equator(self):
        p = latitude_from_meridian_arc(math.radians(0.5), -2 * meridian_distance(math.radians(0.5), WGS84), WGS84)
        assert p.lat == pytest.approx(math.radians(-0.5), abs=1e-14)

    def test_beyond_the_pole(self):
        with pytest.raises(ValueError, match="quarter meridian"):
            latitude_from_meridian_arc(math.radians(89.0), 200_000.0, WGS84)


class TestLongitudeFromParallelArc:
    def test_inverts_parallel_arc_length(self):
        lat, lon1 = math.radians(41.0), math.radians(29.0)
        lon2 = longitude_from_parallel_arc(lat, lon1, 325_000.0, HAYFORD)
        assert lon2 > lon1
        assert parallel_arc_length(lat, lon1, lon2, HAYFORD) == pytest.approx(325_000.0, abs=1e-6)

    def test_westwards_and_across_antimeridian(self):
        lon2 = longitude_from_parallel_arc(0.0, math.radians(-179.5), -WGS84.a * math.radians(1.0), WGS84)
        assert lon2 == pytest.approx(math.radians(179.5))

    def test_pole_is_an_error(self):
        with pytest.raises(ValueError):
            longitude_from_parallel_arc(math.pi / 2, 0.0, 1.0, WGS84)


class TestCourseMeridianSeries:
    """The e′² series and the Hayford tables from the course (Schödlbauer)."""

    def test_hayford_coefficient_table(self):
        s = meridian_series(HAYFORD)
        per_degree = math.pi / 180
        assert s.e0 * per_degree == pytest.approx(111136.536655, abs=2e-6)
        assert s.e2 == pytest.approx(-16107.03468, abs=1e-4)
        assert s.e4 == pytest.approx(16.97621, abs=2e-5)
        assert s.e6 == pytest.approx(-0.02227, abs=1e-5)
        assert abs(s.e8) == pytest.approx(0.00003, abs=5e-6)  # the slide prints it as negative
        assert s.a_ == pytest.approx(0.994955870655, abs=1e-10)
        assert s.b_ == pytest.approx(0.005033498193, abs=2e-11)
        assert s.c_ == pytest.approx(0.000010610237, abs=1e-11)
        assert s.d_ == pytest.approx(0.000000020874, abs=1e-11)
        assert s.a_ + s.b_ + s.c_ + s.d_ + s.e_ == pytest.approx(1.0, abs=1e-15)
        assert s.f2 / per_degree == pytest.approx(0.144930071, abs=1e-9)
        assert s.f4 / per_degree == pytest.approx(0.000213851, abs=1e-9)
        assert s.f6 / per_degree == pytest.approx(0.000000432, abs=1e-9)

    @pytest.mark.parametrize("deg", [-45, 1, 30, 41.0083, 60, 89, 90])
    def test_series_matches_exact_meridian_distance(self, deg):
        lat = math.radians(deg)
        assert meridian_distance_series(lat, HAYFORD) == pytest.approx(
            meridian_distance(lat, HAYFORD), abs=2e-4  # truncated at e′⁸
        )

    @pytest.mark.parametrize("deg", [-45, 1, 30, 41.0083, 60, 89])
    def test_inverse_series(self, deg):
        lat = math.radians(deg)
        g = meridian_distance(lat, HAYFORD)
        assert latitude_from_meridian_distance_series(g, HAYFORD) == pytest.approx(lat, abs=5e-11)

    def test_short_meridian_arc(self):
        lat1, lat2 = math.radians(40.5), math.radians(41.3)
        assert short_meridian_arc(lat1, lat2, HAYFORD) == pytest.approx(
            meridian_distance(lat2, HAYFORD) - meridian_distance(lat1, HAYFORD), abs=1e-4
        )
        assert short_meridian_arc(lat2, lat1, HAYFORD) < 0

    def test_trace(self):
        trace = Trace()
        g = meridian_distance_series(0.7, HAYFORD, trace)
        assert trace["G"] == g
