import math

import pytest

from geolab import (
    GRS80,
    HAYFORD,
    Geodetic,
    HelmertParams,
    Trace,
    helmert,
    helmert_inverse,
    transform_geodetic,
)
from geolab.datum import ARCSEC
from geolab.point import Cartesian3D

P = Cartesian3D(4215751.8328, 2336829.3996, 4163220.0278)
PARAMS = HelmertParams.from_units(87.0, 98.0, 121.0, 0.3, -0.2, 1.1, 2.5)


def test_identity():
    assert helmert(P, HelmertParams()) == P


def test_translation_only():
    q = helmert(P, HelmertParams(tx=1.0, ty=-2.0, tz=3.0))
    assert (q.x - P.x, q.y - P.y, q.z - P.z) == pytest.approx((1.0, -2.0, 3.0))


def test_scale_only():
    q = helmert(P, HelmertParams(scale=1e-6))
    assert math.dist((q.x, q.y, q.z), (P.x, P.y, P.z)) == pytest.approx(
        1e-6 * math.dist((P.x, P.y, P.z), (0, 0, 0)), rel=1e-9
    )


def test_rotation_about_z_matches_the_slide_matrix():
    # X' = X − R_Z Y, Y' = R_Z X + Y, Z' = Z
    rz = 1.0 * ARCSEC
    q = helmert(P, HelmertParams(rz=rz))
    assert q.x == pytest.approx(P.x - rz * P.y, abs=1e-9)
    assert q.y == pytest.approx(rz * P.x + P.y, abs=1e-9)
    assert q.z == pytest.approx(P.z, abs=1e-9)


def test_from_units():
    p = HelmertParams.from_units(rx_arcsec=1.0, scale_ppm=1.0)
    assert p.rx == pytest.approx(math.radians(1 / 3600))
    assert p.scale == pytest.approx(1e-6)


def test_inverse_round_trip():
    back = helmert_inverse(helmert(P, PARAMS), PARAMS)
    assert (back.x, back.y, back.z) == pytest.approx((P.x, P.y, P.z), abs=1e-8)


def test_geodetic_round_trip_between_ellipsoids():
    start = Geodetic(math.radians(41.0083), math.radians(29.0), 150.0)
    ed50 = transform_geodetic(start, GRS80, HAYFORD, PARAMS)
    back = transform_geodetic(ed50, HAYFORD, GRS80, PARAMS, inverse=True)
    assert back.lat == pytest.approx(start.lat, abs=1e-12)
    assert back.lon == pytest.approx(start.lon, abs=1e-12)
    assert back.height == pytest.approx(start.height, abs=1e-6)


def test_zero_parameters_only_change_the_ellipsoid():
    start = Geodetic(math.radians(41.0), math.radians(29.0), 0.0)
    out = transform_geodetic(start, GRS80, HAYFORD, HelmertParams())
    xyz_in = start.to_cartesian(GRS80)
    xyz_out = out.to_cartesian(HAYFORD)
    assert (xyz_out.x, xyz_out.y, xyz_out.z) == pytest.approx((xyz_in.x, xyz_in.y, xyz_in.z), abs=1e-6)
    assert out.lat != start.lat  # same point in space, different ellipsoid


def test_trace():
    trace = Trace()
    q = helmert(P, PARAMS, trace)
    assert trace["X′"] == q.x
