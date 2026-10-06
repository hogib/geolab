import math

import pytest

from geolab import (
    HAYFORD,
    WGS84,
    Ellipsoid,
    Trace,
    azimuth_radius,
    gaussian_mean_radius,
    meridian_radius,
    parallel_radius,
    prime_vertical_radius,
)

LATITUDES = [math.radians(d) for d in (-90, -45, -10, 0, 23.5, 41, 60, 89.9, 90)]
POLAR_RADIUS = WGS84.a**2 / WGS84.b  # c = a²/b, where M = N at the poles


class TestMeridianRadius:
    def test_equator(self):
        assert meridian_radius(0.0, WGS84) == pytest.approx(WGS84.a * (1 - WGS84.e_sq))

    def test_pole(self):
        assert meridian_radius(math.pi / 2, WGS84) == pytest.approx(POLAR_RADIUS)

    def test_reference_value(self):
        # WGS84 at 45°: M = 6 367 381.816 m
        assert meridian_radius(math.radians(45), WGS84) == pytest.approx(6367381.816, abs=1e-3)

    def test_increases_towards_pole(self):
        values = [meridian_radius(math.radians(d), WGS84) for d in range(0, 91, 10)]
        assert values == sorted(values)

    def test_symmetric(self):
        assert meridian_radius(-0.7, WGS84) == meridian_radius(0.7, WGS84)


class TestPrimeVerticalRadius:
    def test_equator(self):
        assert prime_vertical_radius(0.0, WGS84) == WGS84.a

    def test_pole(self):
        assert prime_vertical_radius(math.pi / 2, WGS84) == pytest.approx(POLAR_RADIUS)

    def test_reference_value(self):
        # WGS84 at 45°: N = 6 388 838.290 m
        assert prime_vertical_radius(math.radians(45), WGS84) == pytest.approx(6388838.290, abs=1e-3)

    @pytest.mark.parametrize("lat", LATITUDES)
    def test_n_at_least_m(self, lat):
        assert prime_vertical_radius(lat, WGS84) >= meridian_radius(lat, WGS84) - 1e-6


class TestParallelRadius:
    def test_equator(self):
        assert parallel_radius(0.0, WGS84) == WGS84.a

    def test_pole(self):
        assert parallel_radius(math.pi / 2, WGS84) == pytest.approx(0.0, abs=1e-6)

    @pytest.mark.parametrize("lat", LATITUDES)
    def test_equals_distance_from_axis(self, lat):
        from geolab import Geodetic

        c = Geodetic(lat, 0.3, 0.0).to_cartesian(WGS84)
        assert parallel_radius(lat, WGS84) == pytest.approx(math.hypot(c.x, c.y), abs=1e-6)


class TestGaussianMeanRadius:
    @pytest.mark.parametrize("lat", LATITUDES)
    def test_between_m_and_n(self, lat):
        m, n = meridian_radius(lat, WGS84), prime_vertical_radius(lat, WGS84)
        r = gaussian_mean_radius(lat, WGS84)
        assert m - 1e-6 <= r <= n + 1e-6
        assert r == pytest.approx(math.sqrt(m * n))

    def test_equator(self):
        assert gaussian_mean_radius(0.0, WGS84) == pytest.approx(WGS84.b)


class TestAzimuthRadius:
    @pytest.mark.parametrize("lat", LATITUDES)
    def test_meridian_direction_gives_m(self, lat):
        assert azimuth_radius(lat, 0.0, WGS84) == pytest.approx(meridian_radius(lat, WGS84))
        assert azimuth_radius(lat, math.pi, WGS84) == pytest.approx(meridian_radius(lat, WGS84))

    @pytest.mark.parametrize("lat", LATITUDES)
    def test_east_direction_gives_n(self, lat):
        assert azimuth_radius(lat, math.pi / 2, WGS84) == pytest.approx(
            prime_vertical_radius(lat, WGS84)
        )

    def test_between_m_and_n(self):
        lat = math.radians(41)
        m, n = meridian_radius(lat, HAYFORD), prime_vertical_radius(lat, HAYFORD)
        for deg in range(0, 360, 15):
            r = azimuth_radius(lat, math.radians(deg), HAYFORD)
            assert m - 1e-6 <= r <= n + 1e-6

    def test_symmetric_in_azimuth(self):
        lat = math.radians(41)
        assert azimuth_radius(lat, math.radians(30), WGS84) == pytest.approx(
            azimuth_radius(lat, math.radians(-30), WGS84)
        )
        assert azimuth_radius(lat, math.radians(30), WGS84) == pytest.approx(
            azimuth_radius(lat, math.radians(210), WGS84)
        )


class TestSphere:
    @pytest.mark.parametrize(
        "func",
        [meridian_radius, prime_vertical_radius, gaussian_mean_radius],
    )
    def test_all_radii_equal(self, func):
        sphere = Ellipsoid("Sphere", 6371000.0, 0.0)
        assert func(0.7, sphere) == pytest.approx(6371000.0)

    def test_azimuth(self):
        sphere = Ellipsoid("Sphere", 6371000.0, 0.0)
        assert azimuth_radius(0.7, 1.1, sphere) == pytest.approx(6371000.0)


class TestTrace:
    def test_meridian(self):
        trace = Trace()
        m = meridian_radius(0.7, WGS84, trace)
        assert [s.name for s in trace.steps] == ["W", "M"]
        assert trace["M"] == m

    def test_azimuth(self):
        trace = Trace()
        r = azimuth_radius(0.7, 0.4, WGS84, trace)
        assert [s.name for s in trace.steps] == ["W", "M", "W", "N", "Rα"]
        assert trace["Rα"] == r

    def test_trace_does_not_change_result(self):
        assert gaussian_mean_radius(0.7, WGS84, Trace()) == gaussian_mean_radius(0.7, WGS84)
