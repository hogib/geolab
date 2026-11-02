import math

import pytest

from geolab import WGS84, Trace, ellipse_arc, ellipse_perimeter, elliptic_e, meridian_distance
from geolab.ellipse import parametric_from_normal, parametric_from_polar


def simpson(f, lo: float, hi: float, steps: int = 20000) -> float:
    h = (hi - lo) / steps
    total = f(lo) + f(hi) + sum((4 if i % 2 else 2) * f(lo + i * h) for i in range(1, steps))
    return total * h / 3


def numeric_arc(a: float, b: float, t1: float, t2: float) -> float:
    return abs(simpson(lambda t: math.hypot(a * math.sin(t), b * math.cos(t)), t1, t2))


class TestEllipticE:
    @pytest.mark.parametrize("m", [0.0, 0.3, 0.9, 0.99])
    @pytest.mark.parametrize("phi", [-2.0, 0.4, 1.5, 4.0, 10.0])
    def test_matches_integral(self, m, phi):
        expected = simpson(lambda u: math.sqrt(1 - m * math.sin(u) ** 2), 0.0, phi)
        assert elliptic_e(phi, m) == pytest.approx(expected, abs=1e-12)

    @pytest.mark.parametrize("phi", [-2.0, 0.4, 4.0])
    def test_m_equals_one(self, phi):
        expected = simpson(lambda u: abs(math.cos(u)), 0.0, phi, 200000)
        assert elliptic_e(phi, 1.0) == pytest.approx(expected, abs=1e-9)

    def test_circle(self):
        assert elliptic_e(1.234, 0.0) == pytest.approx(1.234, rel=1e-15)


class TestEllipseArc:
    @pytest.mark.parametrize(
        "a, b, t1, t2",
        [(5, 3, 0, 1.2), (3, 5, 0.3, 2.9), (1, 1e-3, -1, 4), (5, 3, -4, 9), (7, 7, 0, 1)],
    )
    def test_matches_integral(self, a, b, t1, t2):
        assert ellipse_arc(a, b, t1, t2).length == pytest.approx(numeric_arc(a, b, t1, t2), abs=1e-10)

    def test_circle(self):
        assert ellipse_arc(2, 2, 0.5, 2.0).length == pytest.approx(3.0)

    def test_direction_does_not_matter(self):
        assert ellipse_arc(5, 3, 0.2, 1.7).length == pytest.approx(ellipse_arc(5, 3, 1.7, 0.2).length)

    def test_full_turn_is_perimeter(self):
        assert ellipse_arc(5, 3, 0.4, 0.4 + math.tau).length == pytest.approx(ellipse_perimeter(5, 3))

    def test_points_and_chord(self):
        r = ellipse_arc(5, 3, 0, math.pi / 2)
        assert r.p1 == (5, 0)
        assert r.p2 == pytest.approx((0, 3))
        assert r.chord == pytest.approx(math.hypot(5, 3))

    def test_matches_meridian_arc(self):
        # The meridian is an ellipse with semi-axes a, b; its parametric angle is the reduced latitude.
        lat = math.radians(41.0083)
        beta = math.atan2((1 - WGS84.f) * math.sin(lat), math.cos(lat))
        assert ellipse_arc(WGS84.a, WGS84.b, 0, beta).length == pytest.approx(
            meridian_distance(lat, WGS84), abs=1e-6
        )

    def test_rejects_bad_axes(self):
        with pytest.raises(ValueError):
            ellipse_arc(5, 0, 0, 1)

    def test_trace(self):
        trace = Trace()
        r = ellipse_arc(5, 3, 0, 1, trace)
        assert trace["s"] == r.length
        assert trace["P"] == r.perimeter


class TestPerimeter:
    def test_circle(self):
        assert ellipse_perimeter(1, 1) == pytest.approx(math.tau)

    def test_ramanujan(self):
        a, b = 5, 3
        h = ((a - b) / (a + b)) ** 2
        approx = math.pi * (a + b) * (1 + 3 * h / (10 + math.sqrt(4 - 3 * h)))
        assert ellipse_perimeter(a, b) == pytest.approx(approx, rel=1e-9)

    def test_flat_ellipse_tends_to_4a(self):
        assert ellipse_perimeter(1, 1e-9) == pytest.approx(4, rel=1e-6)


class TestPolar:
    def test_axes(self):
        assert parametric_from_polar(math.pi / 2, 5, 3) == pytest.approx(math.pi / 2)

    def test_keeps_the_turn(self):
        assert parametric_from_polar(math.tau + 0.1, 5, 3) == pytest.approx(
            math.tau + parametric_from_polar(0.1, 5, 3)
        )


class TestNormal:
    def test_is_reduced_latitude_on_meridian_ellipse(self):
        lat = math.radians(41.0083)
        beta = math.atan2((1 - WGS84.f) * math.sin(lat), math.cos(lat))
        assert parametric_from_normal(lat, WGS84.a, WGS84.b) == pytest.approx(beta, abs=1e-15)

    def test_keeps_the_turn(self):
        assert parametric_from_normal(-math.tau + 0.3, 5, 3) == pytest.approx(
            parametric_from_normal(0.3, 5, 3) - math.tau
        )
