import math

import pytest

from geolab import HAYFORD, WGS84, Trace, parse_angle
from geolab.area import (
    area_coefficients,
    ellipsoid_area,
    quadrangle_area,
    zone_area,
    zone_area_exact,
)
from geolab.surface import area_element


def test_hayford_coefficients_match_the_course_table():
    # The course table mixes precisions: β and δ match this e¹⁰ series exactly,
    # while α and γ were computed to e⁸ with e² rounded to 0.00672267
    # (a difference of ~1e-11 and ~1e-12). β matches only with 35/192 for the
    # e⁸ term, confirming that the slide's 35/198 is a typo.
    k = area_coefficients(HAYFORD)
    assert k.alpha == pytest.approx(1.00337837836370, abs=2e-11)
    assert k.beta == pytest.approx(0.001128976275697, abs=1e-15)
    assert k.gamma == pytest.approx(0.00000171393466, abs=2e-12)
    assert k.delta == pytest.approx(0.00000000275303, abs=5e-15)  # last printed digit


@pytest.mark.parametrize("ell", [WGS84, HAYFORD])
@pytest.mark.parametrize("lats", [(0, 1), (36, 42), (-20, 35), (60, 90), (-90, 90)])
def test_series_matches_closed_form(ell, lats):
    lat1, lat2 = (math.radians(x) for x in lats)
    # The series stops at e¹⁰; the missing terms are ~1e-12 relative.
    assert zone_area(lat1, lat2, ell) == pytest.approx(zone_area_exact(lat1, lat2, ell), rel=1e-12)


def test_total_area_matches_authalic_radius():
    assert ellipsoid_area(WGS84) == pytest.approx(4 * math.pi * WGS84.authalic_radius**2, rel=1e-12)
    assert ellipsoid_area(WGS84) == pytest.approx(5.10065621724e14, rel=1e-11)


def test_small_quadrangle_matches_area_element():
    lat, lon = math.radians(41.0), math.radians(29.0)
    d = math.radians(0.001)
    expected = area_element(lat, d, d, HAYFORD)
    assert quadrangle_area(lat - d / 2, lat + d / 2, lon, lon + d, HAYFORD) == pytest.approx(expected, rel=1e-6)


def test_map_sheet_1_25000():
    # A 7′30″ × 7′30″ sheet: (7.5′/ρ′) Z/2π. Near 39°N it is about 150 km².
    lat1 = parse_angle("39 0 0")
    lat2 = lat1 + parse_angle("0 7 30")
    p = quadrangle_area(lat1, lat2, 0.0, parse_angle("0 7 30"), HAYFORD)
    assert 145e6 < p < 155e6
    assert p == pytest.approx(parse_angle("0 7 30") / math.tau * zone_area(lat1, lat2, HAYFORD))


def test_order_and_antimeridian():
    a = quadrangle_area(0.1, 0.2, math.radians(179), math.radians(-179), WGS84)
    b = quadrangle_area(0.2, 0.1, 0.0, math.radians(2), WGS84)
    assert a == pytest.approx(b)


def test_hemisphere_is_half():
    assert zone_area(0, math.pi / 2, WGS84) == pytest.approx(ellipsoid_area(WGS84) / 2)


def test_trace():
    trace = Trace()
    p = quadrangle_area(0.6, 0.7, 0.4, 0.5, WGS84, trace)
    assert trace["P"] == p
