import math

import pytest

from geolab import BESSEL_1841, direct, gaussian_mean_radius, inverse, parse_angle
from geolab.trace import Trace
from geolab.triangles import (
    additament_asa,
    legendre_asa,
    plane_asa,
    plane_sas,
    plane_sss,
    spherical_asa,
    spherical_excess,
    spherical_sas,
    spherical_sss,
)

ARCSEC = math.radians(1 / 3600)

# The course example (Bessel ellipsoid): c = P₁P₂ with the angles at P₁ and P₂.
LAT1, LON1 = parse_angle("53 50 02.8809"), parse_angle("10 12 04.1772")
AZ12 = parse_angle("25 16 31.960")
C = 47652.597
ALPHA, BETA = parse_angle("44 08 45.680"), parse_angle("70 46 03.490")
R = gaussian_mean_radius(LAT1, BESSEL_1841)


def ellipsoidal_truth():
    """Intersect the geodesics P₁→P₃ and P₂→P₃ on the ellipsoid (Vincenty)."""
    p2 = direct(LAT1, LON1, AZ12, C, BESSEL_1841)
    az13, az23 = AZ12 + ALPHA, p2.back_azimuth - BETA
    lo, hi = 1e3, 2e5
    for _ in range(100):
        mid = (lo + hi) / 2
        p3 = direct(LAT1, LON1, az13, mid, BESSEL_1841)
        az = inverse(p2.lat2, p2.lon2, p3.lat2, p3.lon2, BESSEL_1841).azimuth1
        lo, hi = (mid, hi) if math.remainder(az - az23, math.tau) > 0 else (lo, mid)
    p3 = direct(LAT1, LON1, az13, lo, BESSEL_1841)
    r23 = inverse(p2.lat2, p2.lon2, p3.lat2, p3.lon2, BESSEL_1841)
    gamma = abs(math.remainder(r23.back_azimuth - p3.back_azimuth, math.tau))
    return r23.distance, lo, gamma


TRUE_A, TRUE_B, TRUE_GAMMA = ellipsoidal_truth()


@pytest.mark.parametrize("solve", [spherical_asa, legendre_asa, additament_asa])
def test_course_example_matches_the_ellipsoidal_triangle(solve):
    t = solve(C, ALPHA, BETA, R)
    assert t.a == pytest.approx(TRUE_A, abs=1e-3)
    assert t.b == pytest.approx(TRUE_B, abs=1e-3)
    assert t.gamma == pytest.approx(TRUE_GAMMA, abs=1e-3 * ARCSEC)
    assert t.excess / ARCSEC == pytest.approx(4.1667, abs=1e-3)


@pytest.mark.parametrize("c", [C, 2e5, 3e5])
def test_legendre_refinement_improves_the_angle(c):
    exact = spherical_asa(c, ALPHA, BETA, R)
    plain = legendre_asa(c, ALPHA, BETA, R)
    fine = legendre_asa(c, ALPHA, BETA, R, refine=True)
    assert abs(fine.gamma - exact.gamma) < abs(plain.gamma - exact.gamma)


def test_large_triangle_shows_the_approximation_limits():
    # 300 km sides: the small-triangle methods drift by decimetres (they are
    # meant for sides under ~100 km); the exact solution stays exact.
    c, alpha, beta = 3e5, math.radians(60), math.radians(60)
    exact = spherical_asa(c, alpha, beta, 6.4e6)
    back = spherical_sss(exact.a, exact.b, exact.c, 6.4e6)
    assert back.alpha == pytest.approx(alpha, abs=1e-12)
    assert 0.01 < abs(legendre_asa(c, alpha, beta, 6.4e6).a - exact.a) < 0.1
    assert 0.01 < abs(additament_asa(c, alpha, beta, 6.4e6).a - exact.a) < 0.3


@pytest.mark.parametrize("tri", [(36594.4, 49609.2, 47652.6), (1e5, 1.2e5, 0.9e5), (5e5, 6e5, 7e5)])
def test_spherical_solvers_are_consistent(tri):
    a, b, c = tri
    t = spherical_sss(a, b, c, 6.38e6)
    # Spherical law of cosines for sides.
    cos_c = math.cos(a / 6.38e6) * math.cos(b / 6.38e6) + math.sin(a / 6.38e6) * math.sin(b / 6.38e6) * math.cos(t.gamma)
    assert math.cos(c / 6.38e6) == pytest.approx(cos_c, abs=1e-14)
    t2 = spherical_sas(a, b, t.gamma, 6.38e6)
    assert (t2.c, t2.alpha, t2.beta) == pytest.approx((c, t.alpha, t.beta), rel=1e-11)
    t3 = spherical_asa(c, t.alpha, t.beta, 6.38e6)
    assert (t3.a, t3.b, t3.gamma) == pytest.approx((a, b, t.gamma), rel=1e-11)
    assert t.excess > 0 and t.area == pytest.approx(t.excess * 6.38e6**2)


def test_excess_from_plane_area():
    plane = plane_sss(36594.4366, 49609.1920, C)
    assert spherical_excess(plane.area, R) / ARCSEC == pytest.approx(4.1667, abs=1e-3)


def test_plane_solvers():
    t = plane_sss(3.0, 4.0, 5.0)
    assert t.gamma == pytest.approx(math.pi / 2)
    assert t.area == pytest.approx(6.0)
    t2 = plane_sas(3.0, 4.0, math.pi / 2)
    assert (t2.c, t2.alpha, t2.beta) == pytest.approx((5.0, t.alpha, t.beta))
    t3 = plane_asa(5.0, t.alpha, t.beta)
    assert (t3.a, t3.b) == pytest.approx((3.0, 4.0))
    assert t.alpha + t.beta + t.gamma == pytest.approx(math.pi)


def test_invalid_triangles():
    with pytest.raises(ValueError):
        plane_sss(1, 2, 5)
    with pytest.raises(ValueError):
        plane_asa(1, math.radians(100), math.radians(90))


def test_trace():
    trace = Trace()
    t = spherical_asa(C, ALPHA, BETA, R, trace)
    assert trace["γ"] == t.gamma
