import math

import pytest

from geolab import BESSEL_1841, HAYFORD, WGS84, ConvergenceError, Ellipsoid, Trace, parse_angle
from geolab.geodesic import direct, inverse
from geolab.trace import Section

ARCSEC = math.radians(1 / 3600)


def p(text):
    return parse_angle(text)


# Vincenty (1975), Survey Review XXIII(176), test lines (a)–(e).
# (ellipsoid, φ1, φ2, L = λ2 − λ1, α1, α2, s)
VINCENTY_LINES = {
    "a": (BESSEL_1841, "55 45 0", "-33 26 0", "108 13 0",
          "96 36 8.79960", "137 52 22.01454", 14110526.170),
    "b": (HAYFORD, "37 19 54.95367", "26 7 42.83946", "41 28 35.50729",
          "95 27 59.63089", "118 5 58.96161", 4085966.703),
    "c": (HAYFORD, "35 16 11.24862", "67 22 14.77638", "137 47 28.31435",
          "15 44 23.74850", "144 55 39.92147", 8084823.839),
    "d": (HAYFORD, "1 0 0", "-0 59 53.83076", "179 17 48.02997",
          "89 0 0", "91 0 6.11733", 19960000.000),
    "e": (HAYFORD, "1 0 0", "1 1 15.18952", "179 46 17.84244",
          "4 59 59.99995", "174 59 59.88481", 19780006.558),
}

# WGS84 reference values from GeographicLib (Karney): (φ1, λ1, φ2, λ2) → (s, α1, α2), degrees.
GEOGRAPHICLIB_LINES = [
    ((41.0, 29.0, 39.9, 32.85), (348653.1288, 109.249178514, 111.747659600)),
    ((-33.9, 151.2, 51.5, -0.13), (16991896.1089, 319.233974803, 240.428514810)),
    ((0.0, 0.0, 0.0, 90.0), (10018754.1714, 90.0, 90.0)),
    ((90.0, 0.0, -90.0, 0.0), (20003931.4586, 180.0, 180.0)),
]


class TestVincentyLines:
    @pytest.mark.parametrize("line", VINCENTY_LINES)
    def test_inverse(self, line):
        ell, f1, f2, L, a1, a2, s = VINCENTY_LINES[line]
        r = inverse(p(f1), 0.0, p(f2), p(L), ell)
        assert r.distance == pytest.approx(s, abs=1e-3)
        # Line (d) is nearly antipodal: its azimuths are poorly conditioned.
        az_tol = 2e-3 * ARCSEC if line == "d" else 1e-4 * ARCSEC
        assert r.azimuth1 == pytest.approx(p(a1), abs=az_tol)
        assert r.azimuth2 == pytest.approx(p(a2), abs=az_tol)

    @pytest.mark.parametrize("line", VINCENTY_LINES)
    def test_direct(self, line):
        ell, f1, f2, L, a1, a2, s = VINCENTY_LINES[line]
        r = direct(p(f1), 0.0, p(a1), s, ell)
        assert r.lat2 == pytest.approx(p(f2), abs=1e-4 * ARCSEC)
        assert r.lon2 == pytest.approx(p(L), abs=1e-4 * ARCSEC)
        assert r.azimuth2 == pytest.approx(p(a2), abs=1e-4 * ARCSEC)


class TestAgainstGeographicLib:
    @pytest.mark.parametrize(("points", "expected"), GEOGRAPHICLIB_LINES)
    def test_inverse(self, points, expected):
        s, a1, a2 = expected
        r = inverse(*map(math.radians, points), WGS84)
        assert r.distance == pytest.approx(s, abs=1e-4)
        assert math.degrees(r.azimuth1) == pytest.approx(a1, abs=1e-8)
        assert math.degrees(r.azimuth2) == pytest.approx(a2, abs=1e-8)

    @pytest.mark.parametrize(("points", "expected"), GEOGRAPHICLIB_LINES[:2])
    def test_direct(self, points, expected):
        lat1, lon1, lat2, lon2 = map(math.radians, points)
        s, a1, a2 = expected
        r = direct(lat1, lon1, math.radians(a1), s, WGS84)
        assert r.lat2 == pytest.approx(lat2, abs=1e-5 * ARCSEC)
        assert r.lon2 == pytest.approx(lon2, abs=1e-5 * ARCSEC)
        assert math.degrees(r.azimuth2) == pytest.approx(a2, abs=1e-8)


class TestInverseSpecialCases:
    def test_coincident_points(self):
        r = inverse(0.3, 0.4, 0.3, 0.4, WGS84)
        assert (r.distance, r.azimuth1, r.azimuth2) == (0.0, 0.0, 0.0)

    def test_along_equator(self):
        r = inverse(0.0, 0.0, 0.0, math.radians(1), WGS84)
        assert r.distance == pytest.approx(WGS84.a * math.radians(1), abs=1e-6)
        assert r.azimuth1 == pytest.approx(math.pi / 2)

    def test_along_meridian_north(self):
        r = inverse(math.radians(10), 0.5, math.radians(20), 0.5, WGS84)
        assert r.azimuth1 == pytest.approx(0.0, abs=1e-12)
        assert r.azimuth2 == pytest.approx(0.0, abs=1e-12)

    def test_along_meridian_south(self):
        r = inverse(math.radians(20), 0.5, math.radians(10), 0.5, WGS84)
        assert r.azimuth1 == pytest.approx(math.pi)

    def test_westward_azimuth_in_range(self):
        r = inverse(0.0, 0.0, 0.0, math.radians(-10), WGS84)
        assert r.azimuth1 == pytest.approx(3 * math.pi / 2)
        assert 0 <= r.azimuth1 < math.tau

    def test_longitude_wraps_across_antimeridian(self):
        r1 = inverse(0.5, math.radians(179.5), 0.5, math.radians(-179.5), WGS84)
        r2 = inverse(0.5, math.radians(-0.5), 0.5, math.radians(0.5), WGS84)
        assert r1.distance == pytest.approx(r2.distance, abs=1e-6)

    def test_symmetric(self):
        args = (math.radians(41), math.radians(29), math.radians(-12), math.radians(77))
        forward = inverse(*args, HAYFORD)
        backward = inverse(args[2], args[3], args[0], args[1], HAYFORD)
        assert forward.distance == pytest.approx(backward.distance, abs=1e-6)
        assert backward.azimuth1 == pytest.approx(forward.back_azimuth)

    def test_back_azimuth(self):
        r = inverse(0.0, 0.0, 0.0, math.radians(10), WGS84)
        assert r.back_azimuth == pytest.approx(3 * math.pi / 2)

    @pytest.mark.parametrize("points", [(0, 0, 0.5, 179.7), (0, 0, 0, 179.9), (0, 0, -0.5, 179.7)])
    def test_nearly_antipodal_raises(self, points):
        with pytest.raises(ConvergenceError, match="antipodal"):
            inverse(*map(math.radians, points), WGS84)

    def test_iteration_limit(self):
        with pytest.raises(ConvergenceError):
            inverse(0.1, 0.0, 0.5, 1.0, WGS84, max_iter=1, tol=0.0)

    def test_sphere_matches_great_circle(self):
        r_earth = 6371000.0
        sphere = Ellipsoid("Sphere", r_earth, 0.0)
        lat1, lon1, lat2, lon2 = map(math.radians, (41.0, 29.0, -33.9, 151.2))
        central = math.acos(
            math.sin(lat1) * math.sin(lat2) + math.cos(lat1) * math.cos(lat2) * math.cos(lon2 - lon1)
        )
        assert inverse(lat1, lon1, lat2, lon2, sphere).distance == pytest.approx(r_earth * central)


class TestDirectSpecialCases:
    def test_zero_distance(self):
        r = direct(0.3, 0.4, 1.0, 0.0, WGS84)
        assert (r.lat2, r.lon2) == pytest.approx((0.3, 0.4))
        assert r.azimuth2 == pytest.approx(1.0)

    def test_due_north_to_pole(self):
        quarter_meridian = inverse(0.0, 0.0, math.pi / 2, 0.0, WGS84).distance
        r = direct(0.0, 0.0, 0.0, quarter_meridian, WGS84)
        assert r.lat2 == pytest.approx(math.pi / 2, abs=1e-9)

    def test_due_east_on_equator(self):
        r = direct(0.0, 0.0, math.pi / 2, WGS84.a * math.radians(10), WGS84)
        assert r.lat2 == pytest.approx(0.0, abs=1e-12)
        assert r.lon2 == pytest.approx(math.radians(10))

    def test_longitude_normalised(self):
        r = direct(0.0, math.radians(179), math.pi / 2, WGS84.a * math.radians(2), WGS84)
        assert r.lon2 == pytest.approx(math.radians(-179))

    def test_negative_distance_goes_backwards(self):
        forward = direct(0.5, 0.5, 1.0, 1e5, WGS84)
        backward = direct(0.5, 0.5, 1.0 + math.pi, -1e5, WGS84)
        assert (backward.lat2, backward.lon2) == pytest.approx((forward.lat2, forward.lon2))

    def test_iteration_limit(self):
        with pytest.raises(ConvergenceError):
            direct(0.1, 0.0, 0.5, 1e7, WGS84, max_iter=1, tol=0.0)


class TestConsistency:
    @pytest.mark.parametrize(
        "points",
        [
            (41.0, 29.0, 39.9, 32.85),
            (-60.0, -70.0, 10.0, 120.0),
            (89.0, 0.0, 89.0, 180.0),
            (0.001, 0.0, -0.001, 0.002),
            (-45.0, 170.0, 45.0, -170.0),
        ],
    )
    def test_direct_inverts_inverse(self, points):
        lat1, lon1, lat2, lon2 = map(math.radians, points)
        inv = inverse(lat1, lon1, lat2, lon2, HAYFORD)
        d = direct(lat1, lon1, inv.azimuth1, inv.distance, HAYFORD)
        assert d.lat2 == pytest.approx(lat2, abs=1e-5 * ARCSEC)
        assert math.remainder(d.lon2 - lon2, math.tau) == pytest.approx(0.0, abs=1e-5 * ARCSEC)
        assert d.azimuth2 == pytest.approx(inv.azimuth2, abs=1e-5 * ARCSEC)


class TestTrace:
    def test_inverse(self):
        trace = Trace()
        r = inverse(math.radians(41), math.radians(29), math.radians(39.9), math.radians(32.85), WGS84, trace)
        sections = [e.title for e in trace.entries if isinstance(e, Section)]
        assert sections[0] == "Setup"
        assert sections[-1] == "Result"
        assert sections[1:-1] == [f"Iteration {i}" for i in range(1, r.iterations + 1)]
        assert trace["s"] == r.distance
        assert trace["α₁"] == r.azimuth1
        assert trace["α₂"] == r.azimuth2

    def test_direct(self):
        trace = Trace()
        r = direct(math.radians(41), math.radians(29), 1.9, 348653.0, WGS84, trace)
        sections = [e.title for e in trace.entries if isinstance(e, Section)]
        assert sections == ["Setup", *[f"Iteration {i}" for i in range(1, r.iterations + 1)], "Result"]
        assert trace["φ₂"] == r.lat2
        assert trace["λ₂"] == r.lon2
        assert trace["α₂"] == r.azimuth2

    def test_trace_does_not_change_result(self):
        args = (0.7, 0.5, -0.2, 1.4, WGS84)
        assert inverse(*args, Trace()) == inverse(*args)
