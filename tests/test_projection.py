"""Projection theory: Tissot distortion, isometric latitude, Mercator."""

import math

import pytest

from geolab import BESSEL_1841, HAYFORD, WGS84
from geolab.curvature import meridian_radius, prime_vertical_radius
from geolab.projection import (
    isometric_coefficients,
    isometric_difference,
    isometric_latitude,
    latitude_coefficients,
    latitude_difference,
    latitude_from_isometric,
    mercator_sphere,
    mercator_sphere_inverse,
    tissot,
)
from geolab.trace import Trace

ARCSEC = math.radians(1 / 3600)


class TestTissot:
    def test_principal_directions(self):
        assert tissot(1.2, 0.8, 0).scale == pytest.approx(1.2)
        assert tissot(1.2, 0.8, math.pi / 2).scale == pytest.approx(0.8)

    def test_largest_angle_distortion_by_search(self):
        m1, m2 = 1.3, 0.9
        best = max((t - tissot(m1, m2, t).image_angle, t) for t in (i * 1e-5 for i in range(157080)))
        r = tissot(m1, m2, 0.3)
        assert r.max_angle_distortion == pytest.approx(best[0], abs=1e-9)
        assert r.max_direction == pytest.approx(best[1], abs=1e-4)
        assert r.max_direction == pytest.approx(math.pi / 4 + r.max_angle_distortion / 2)

    def test_image_ratio(self):
        r = tissot(1.1, 0.95, 0.7)
        assert math.tan(0.7) / math.tan(r.image_angle) == pytest.approx(1.1 / 0.95)

    def test_conformal_and_equal_area(self):
        conformal = tissot(1.004, 1.004, 1.0)
        assert conformal.max_angle_distortion == 0 and conformal.angle_distortion == pytest.approx(0)
        equal_area = tissot(2.0, 0.5, 0.2)
        assert equal_area.area_scale == pytest.approx(1) and equal_area.k == pytest.approx(1)

    def test_rejects_non_positive(self):
        with pytest.raises(ValueError):
            tissot(0, 1, 0)


def simpson(f, a, b, n=2000):
    h = (b - a) / n
    return h / 3 * sum((1 if i in (0, n) else 4 if i % 2 else 2) * f(a + i * h) for i in range(n + 1))


class TestIsometricLatitude:
    @pytest.mark.parametrize("ell", [WGS84, HAYFORD, BESSEL_1841])
    def test_closed_form_is_the_integral(self, ell):
        lat = math.radians(41)
        integral = simpson(lambda p: meridian_radius(p, ell) / (prime_vertical_radius(p, ell) * math.cos(p)), 0, lat)
        assert isometric_latitude(lat, ell) == pytest.approx(integral, abs=1e-12)

    @pytest.mark.parametrize("deg", [-80, -30, 0, 0.5, 41, 89])
    def test_round_trip(self, deg):
        q = isometric_latitude(math.radians(deg), WGS84)
        assert latitude_from_isometric(q, WGS84) == pytest.approx(math.radians(deg), abs=1e-14)

    def test_coefficients_are_derivatives(self):
        # Cₖ = q⁽ᵏ⁾/k!, dₖ = φ⁽ᵏ⁾(q)/k!, checked by finite differences.
        lat0, h = math.radians(40), 1e-3
        q = lambda p: isometric_latitude(p, HAYFORD)  # noqa: E731
        c = isometric_coefficients(lat0, HAYFORD)
        assert (q(lat0 + h) - q(lat0 - h)) / (2 * h) == pytest.approx(c[0], rel=1e-6)
        assert (q(lat0 + h) - 2 * q(lat0) + q(lat0 - h)) / h**2 / 2 == pytest.approx(c[1], rel=1e-5)
        q0 = q(lat0)
        p = lambda x: latitude_from_isometric(q0 + x, HAYFORD)  # noqa: E731
        d = latitude_coefficients(lat0, HAYFORD)
        assert (p(h) - p(-h)) / (2 * h) == pytest.approx(d[0], rel=1e-6)
        assert (p(h) - 2 * p(0) + p(-h)) / h**2 / 2 == pytest.approx(d[1], rel=1e-5)
        third = (p(2 * h) - 2 * p(h) + 2 * p(-h) - p(-2 * h)) / (2 * h**3) / 6
        assert third == pytest.approx(d[2], rel=1e-4)  # negative: the slide's sign is wrong

    @pytest.mark.parametrize("ddeg", [-1.5, -0.2, 0.3, 1.0, 2.0])
    def test_series_match_exact(self, ddeg):
        lat0 = math.radians(40)
        lat = lat0 + math.radians(ddeg)
        exact_dq = isometric_latitude(lat, HAYFORD) - isometric_latitude(lat0, HAYFORD)
        dq = isometric_difference(lat, lat0, HAYFORD)
        assert dq == pytest.approx(exact_dq, abs=1e-9)  # truncation ~Δφ⁶
        assert latitude_difference(exact_dq, lat0, HAYFORD) == pytest.approx(lat - lat0, abs=1e-5 * ARCSEC)

    def test_trace(self):
        trace = Trace()
        dq = isometric_difference(0.71, 0.70, WGS84, trace)
        assert trace["Δq"] == dq and trace["C3"] > 0


class TestMercator:
    def test_formula(self):
        p = mercator_sphere(math.radians(45), math.radians(10), 6371000, math.radians(4))
        assert p.x == pytest.approx(6371000 * math.log(1 + math.sqrt(2)))
        assert p.y == pytest.approx(6371000 * math.radians(6))
        assert p.scale == pytest.approx(math.sqrt(2))

    def test_is_conformal(self):
        # Scale along the meridian, R dφ → dx, equals the scale along the parallel.
        r, lat, h = 6371000.0, math.radians(37), 1e-6
        dx = mercator_sphere(lat + h, 0, r).x - mercator_sphere(lat - h, 0, r).x
        assert dx / (2 * h * r) == pytest.approx(1 / math.cos(lat), rel=1e-9)

    def test_round_trip(self):
        p = mercator_sphere(math.radians(-33.5), math.radians(151), 6371000, math.radians(150))
        lat, lon = mercator_sphere_inverse(p.x, p.y, 6371000, math.radians(150))
        assert (lat, lon) == pytest.approx((math.radians(-33.5), math.radians(151)), abs=1e-14)

    def test_pole(self):
        with pytest.raises(ValueError):
            mercator_sphere(math.pi / 2, 0, 1)
