"""Gauss–Krüger reductions, checked against exact geometry: Krüger's series for
the projection and Vincenty for the geodesics."""

import math

import pytest

from geolab import WGS84, inverse
from geolab.area import quadrangle_area
from geolab.gauss_kruger import krueger_forward, krueger_inverse
from geolab.gk_reductions import direct_by_reductions, reduce_area, reduce_line, traverse
from geolab.trace import Trace

ARCSEC = math.radians(1 / 3600)
L0 = math.radians(30)


def xy(lat, lon):
    p = krueger_forward(lat, lon, L0, WGS84)
    return p.x, p.y


def exact_line(pi, pk):
    """S, Tᵢₖ, Tₖᵢ of the geodesic between two plane points."""
    a, b = krueger_inverse(*pi, L0, WGS84), krueger_inverse(*pk, L0, WGS84)
    r = inverse(*a, *b, WGS84)
    ci = krueger_forward(*a, L0, WGS84).convergence
    ck = krueger_forward(*b, L0, WGS84).convergence
    return r.distance, (r.azimuth1 - ci) % math.tau, (r.back_azimuth - ck) % math.tau


LINES = [
    ((4.5e6, 2e4), (4.52e6, 4.5e4)),
    ((4.5e6, -1e5), (4.47e6, -6e4)),
    ((4.6e6, 1.3e5), (4.63e6, 1.0e5)),
]


@pytest.mark.parametrize("pi, pk", LINES)
def test_second_problem(pi, pk):
    s, t_ik, t_ki = exact_line(pi, pk)
    r = reduce_line(*pi, *pk, WGS84)
    assert r.distance == pytest.approx(s, abs=1e-4)
    assert r.distance < r.plane_distance  # s > S always
    assert math.remainder(r.bearing1 - t_ik, math.tau) == pytest.approx(0, abs=0.002 * ARCSEC)
    assert math.remainder(r.bearing2 - t_ki, math.tau) == pytest.approx(0, abs=0.002 * ARCSEC)


@pytest.mark.parametrize("pi, pk", LINES)
def test_first_problem(pi, pk):
    s, t_ik, t_ki = exact_line(pi, pk)
    r = direct_by_reductions(*pi, t_ik, s, WGS84)
    assert (r.x2, r.y2) == pytest.approx(pk, abs=5e-4)
    assert math.remainder(r.bearing2 - t_ki, math.tau) == pytest.approx(0, abs=0.005 * ARCSEC)


def test_reductions_vanish_on_the_central_meridian():
    r = reduce_line(4.5e6, 0, 4.53e6, 0, WGS84)
    assert r.distance_reduction == 0 and r.bearing_reduction1 == 0


def sheet_outline(lat1, lat2, lon1, lon2, n=200):
    corners = [(lat1, lon1), (lat1, lon2), (lat2, lon2), (lat2, lon1)]
    pts = []
    for (pa, la), (pb, lb) in zip(corners, corners[1:] + corners[:1]):
        pts += [xy(pa + (pb - pa) * i / n, la + (lb - la) * i / n) for i in range(n)]
    return pts


@pytest.mark.parametrize("lon1", [30.0, 31.0, 31.375])
def test_area_of_a_map_sheet(lon1):
    # A 1/25 000 sheet (7.5′ × 7.5′) at up to 1.5° from the central meridian.
    lat1, lat2 = math.radians(40.5), math.radians(40.625)
    lon1, lon2 = math.radians(lon1), math.radians(lon1 + 0.125)
    pts = sheet_outline(lat1, lat2, lon1, lon2)
    exact = quadrangle_area(lat1, lat2, lon1, lon2, WGS84)
    r = reduce_area(pts, WGS84)
    assert r.area == pytest.approx(exact, rel=1e-9)
    assert r.reduction < 0
    assert reduce_area(pts[::-1], WGS84).area == pytest.approx(r.area)  # orientation


def test_area_matches_course_formula_for_a_rectangle():
    yi, yk, x0 = 40000.0, 60000.0, 4.5e6
    rect = [(x0, yi), (x0, yk), (x0 + 3000, yk), (x0 + 3000, yi)]
    r = reduce_area(rect, WGS84)
    f, ym = r.plane_area, (yi + yk) / 2
    course = -f / (3 * r.radius**2) * (yk**2 + yk * yi + yi**2) + 2 * ym**4 * f / (3 * r.radius**4)
    assert r.reduction == pytest.approx(course, rel=1e-3)


def test_area_needs_a_polygon():
    with pytest.raises(ValueError):
        reduce_area([(0, 0), (1, 1)], WGS84)


# A traverse A–B–P₁–P₂–P₃–C–D with exact (ellipsoidal) angles and sides.
STATIONS = [(40.50, 30.90), (40.52, 30.95), (40.53, 31.00), (40.545, 31.04), (40.55, 31.09), (40.56, 31.13), (40.59, 31.16)]


def measured_traverse():
    geo = [(math.radians(a), math.radians(b)) for a, b in STATIONS]
    plane = [xy(*g) for g in geo]
    angles = []
    for k in range(1, len(geo) - 1):
        back = inverse(*geo[k], *geo[k - 1], WGS84).azimuth1
        fore = inverse(*geo[k], *geo[k + 1], WGS84).azimuth1
        angles.append((fore - back) % math.tau)
    sides = [inverse(*geo[k], *geo[k + 1], WGS84).distance for k in range(1, len(geo) - 2)]
    return plane, angles, sides


def test_traverse_with_exact_measurements_closes():
    plane, angles, sides = measured_traverse()
    t = traverse(plane[0], plane[1], plane[-2], plane[-1], angles, sides, WGS84)
    assert t.angle_misclosure == pytest.approx(0, abs=0.01 * ARCSEC)
    assert math.hypot(t.misclosure_x, t.misclosure_y) < 1e-4
    for got, want in zip(t.points, plane[2:-2]):
        assert got == pytest.approx(want, abs=1e-4)
    # Without the reductions the traverse would not close: they are tens of mm here.
    assert max(abs(s - p) for s, p in zip(sides, t.plane_distances)) > 0.01


def test_traverse_distributes_misclosure():
    plane, angles, sides = measured_traverse()
    angles[2] += 6 * ARCSEC
    sides[1] += 0.05
    trace = Trace()
    t = traverse(plane[0], plane[1], plane[-2], plane[-1], angles, sides, WGS84, trace)
    assert t.angle_misclosure == pytest.approx(-6 * ARCSEC, abs=0.01 * ARCSEC)
    assert math.hypot(t.misclosure_x, t.misclosure_y) > 0.03
    for got, want in zip(t.points, plane[2:-2]):
        assert got == pytest.approx(want, abs=0.1)
    assert trace["w_β"] == t.angle_misclosure


def test_traverse_shape_check():
    with pytest.raises(ValueError):
        traverse((0, 0), (1, 1), (2, 2), (3, 3), [0.1, 0.2], [100.0, 100.0], WGS84)
