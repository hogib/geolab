import math

import pytest

from geolab import (
    HAYFORD,
    WGS84,
    Cartesian3D,
    ConvergenceError,
    Ellipsoid,
    Geodetic,
    GeodeticMethod,
    Trace,
    parse_angle,
)
from geolab.trace import Section

ITERATIVE, BOWRING = GeodeticMethod.ITERATIVE, GeodeticMethod.BOWRING

# EPSG Guidance Note 7-2, geographic/geocentric conversions example (WGS84).
EPSG_LAT = parse_angle("53°48'33.82\"N")
EPSG_LON = parse_angle("2°07'46.38\"E")
EPSG_H = 73.0
EPSG_XYZ = (3771793.968, 140253.342, 5124304.349)

LATITUDES = [-90.0, -89.9999, -60.0, -23.5, 0.0, 0.0001, 41.0, 45.0, 80.0, 89.9999, 90.0]
LONGITUDES = [-180.0, -90.0, -0.5, 0.0, 29.0, 135.0, 179.9]
METRE_IN_RAD = 1 / 6.4e6


def assert_geodetic_close(result, expected, tol_m):
    assert result.lat == pytest.approx(expected.lat, abs=tol_m * METRE_IN_RAD)
    assert result.height == pytest.approx(expected.height, abs=tol_m)
    if abs(expected.lat) < math.radians(89.99):  # longitude is meaningless at the poles
        dlon = math.remainder(result.lon - expected.lon, math.tau)
        assert abs(dlon) * math.cos(expected.lat) < tol_m * METRE_IN_RAD


class TestEPSGExample:
    def test_forward(self):
        c = Geodetic(EPSG_LAT, EPSG_LON, EPSG_H).to_cartesian(WGS84)
        assert (c.x, c.y, c.z) == pytest.approx(EPSG_XYZ, abs=1e-3)

    @pytest.mark.parametrize("method", list(GeodeticMethod))
    def test_reverse(self, method):
        g = Cartesian3D(*EPSG_XYZ).to_geodetic(WGS84, method)
        # The published inputs are rounded to 0.01" (~0.3 m) and 1 mm.
        assert math.degrees(g.lat) == pytest.approx(math.degrees(EPSG_LAT), abs=0.0005 / 3600)
        assert math.degrees(g.lon) == pytest.approx(math.degrees(EPSG_LON), abs=0.0005 / 3600)
        assert g.height == pytest.approx(EPSG_H, abs=1e-3)


class TestIterative:
    @pytest.mark.parametrize("lat", LATITUDES)
    @pytest.mark.parametrize("h", [-5000.0, 0.0, 150.0, 8848.0, 400e3, 20e6])
    def test_round_trip(self, lat, h):
        expected = Geodetic(math.radians(lat), math.radians(29.0), h)
        result = expected.to_cartesian(WGS84).to_geodetic(WGS84, ITERATIVE)
        assert_geodetic_close(result, expected, tol_m=1e-6)

    @pytest.mark.parametrize("lon", LONGITUDES)
    def test_longitude_quadrants(self, lon):
        expected = Geodetic(math.radians(41.0), math.radians(lon), 100.0)
        result = expected.to_cartesian(WGS84).to_geodetic(WGS84)
        assert_geodetic_close(result, expected, tol_m=1e-6)

    def test_default_method_is_iterative(self):
        c = Geodetic(0.7, 0.5, 100.0).to_cartesian(HAYFORD)
        assert c.to_geodetic(HAYFORD) == c.to_geodetic(HAYFORD, ITERATIVE)

    def test_hits_iteration_limit(self):
        c = Geodetic(0.7, 0.5, 100.0).to_cartesian(WGS84)
        with pytest.raises(ConvergenceError):
            c.to_geodetic(WGS84, ITERATIVE, max_iter=1, tol=0.0)

    def test_converges_in_few_iterations(self):
        trace = Trace()
        Geodetic(0.7, 0.5, 100.0).to_cartesian(WGS84).to_geodetic(WGS84, trace=trace)
        iterations = [e for e in trace.entries if isinstance(e, Section) and e.title.startswith("Iter")]
        assert 1 <= len(iterations) <= 5

    def test_trace(self):
        trace = Trace()
        result = Cartesian3D(*EPSG_XYZ).to_geodetic(WGS84, ITERATIVE, trace)
        names = [s.name for s in trace.steps]
        assert names[:3] == ["p", "λ", "φ₀"]
        assert names[-2:] == ["N", "h"]
        assert trace["p"] == pytest.approx(math.hypot(EPSG_XYZ[0], EPSG_XYZ[1]))
        assert trace["φ"] == result.lat
        assert trace["h"] == result.height
        assert trace.entries[-3] == Section("Result")


class TestBowring:
    @pytest.mark.parametrize("lat", LATITUDES)
    @pytest.mark.parametrize("h", [-5000.0, 0.0, 150.0, 8848.0])
    def test_round_trip_near_surface(self, lat, h):
        expected = Geodetic(math.radians(lat), math.radians(29.0), h)
        result = expected.to_cartesian(WGS84).to_geodetic(WGS84, BOWRING)
        assert_geodetic_close(result, expected, tol_m=1e-6)

    @pytest.mark.parametrize("lat", [10.0, 45.0, 80.0])
    def test_degrades_gracefully_at_altitude(self, lat):
        # Bowring is an approximation; at ~400 km the error grows to ~1 mm.
        expected = Geodetic(math.radians(lat), 0.0, 400e3)
        result = expected.to_cartesian(WGS84).to_geodetic(WGS84, BOWRING)
        assert_geodetic_close(result, expected, tol_m=5e-3)

    def test_agrees_with_iterative(self):
        c = Geodetic(math.radians(-33.9), math.radians(151.2), 58.0).to_cartesian(HAYFORD)
        assert_geodetic_close(
            c.to_geodetic(HAYFORD, BOWRING), c.to_geodetic(HAYFORD, ITERATIVE), tol_m=1e-6
        )

    def test_trace(self):
        trace = Trace()
        result = Cartesian3D(*EPSG_XYZ).to_geodetic(WGS84, BOWRING, trace)
        assert [s.name for s in trace.steps] == ["p", "λ", "θ", "φ", "N", "h"]
        assert trace["φ"] == result.lat
        assert trace["h"] == result.height


class TestSpecialPoints:
    @pytest.mark.parametrize("method", list(GeodeticMethod))
    def test_equator(self, method):
        g = Cartesian3D(WGS84.a + 10.0, 0.0, 0.0).to_geodetic(WGS84, method)
        assert (g.lat, g.lon) == (0.0, 0.0)
        assert g.height == pytest.approx(10.0, abs=1e-9)

    @pytest.mark.parametrize("method", list(GeodeticMethod))
    @pytest.mark.parametrize("sign", [1, -1])
    def test_poles(self, method, sign):
        g = Cartesian3D(0.0, 0.0, sign * (WGS84.b + 10.0)).to_geodetic(WGS84, method)
        assert g.lat == pytest.approx(sign * math.pi / 2)
        assert g.height == pytest.approx(10.0, abs=1e-9)

    @pytest.mark.parametrize("method", list(GeodeticMethod))
    def test_centre_is_undefined(self, method):
        with pytest.raises(ValueError, match="centre"):
            Cartesian3D(0.0, 0.0, 0.0).to_geodetic(WGS84, method)

    @pytest.mark.parametrize("method", list(GeodeticMethod))
    def test_sphere(self, method):
        sphere = Ellipsoid("Sphere", 6371000.0, 0.0)
        g = Cartesian3D(1e6, 2e6, 6e6).to_geodetic(sphere, method)
        r = math.sqrt(1e12 + 4e12 + 36e12)
        assert g.lat == pytest.approx(math.asin(6e6 / r))
        assert g.lon == pytest.approx(math.atan2(2e6, 1e6))
        assert g.height == pytest.approx(r - 6371000.0, abs=1e-6)
