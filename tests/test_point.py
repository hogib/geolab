import dataclasses
import math

import pytest

from geolab import HAYFORD, Ellipsoid
from geolab.point import Cartesian3D, Geodetic, Latitude, LatType

# Reference values for Hayford, computed independently of the library.
E_SQ = 0.006722670022333322
A = 6378388.0
B = 6356911.9461279465

LATITUDES_DEG = [-89.0, -60.0, -45.0, -10.0, 0.0, 10.0, 30.0, 45.0, 60.0, 89.0]


@pytest.fixture
def hayford():
    return HAYFORD


def geocentric_from_geodetic(phi):
    return math.atan((1 - E_SQ) * math.tan(phi))


def parametric_from_geodetic(phi):
    return math.atan(math.sqrt(1 - E_SQ) * math.tan(phi))


class TestLatitudeConversion:
    @pytest.mark.parametrize("lat_type", list(LatType))
    def test_same_type_returns_equal_copy(self, hayford, lat_type):
        lat = Latitude(0.5, lat_type)
        result = lat.convert_to(lat_type, hayford)
        assert result == lat
        assert result is not lat

    @pytest.mark.parametrize("deg", LATITUDES_DEG)
    def test_geodetic_to_geocentric(self, hayford, deg):
        phi = math.radians(deg)
        result = Latitude(phi, LatType.GEODETIC).convert_to(LatType.GEOCENTRIC, hayford)
        assert result.lat_type == LatType.GEOCENTRIC
        assert result.value == pytest.approx(geocentric_from_geodetic(phi), abs=1e-12)

    @pytest.mark.parametrize("deg", LATITUDES_DEG)
    def test_geodetic_to_parametric(self, hayford, deg):
        phi = math.radians(deg)
        result = Latitude(phi, LatType.GEODETIC).convert_to(LatType.PARAMETRIC, hayford)
        assert result.lat_type == LatType.PARAMETRIC
        assert result.value == pytest.approx(parametric_from_geodetic(phi), abs=1e-12)

    @pytest.mark.parametrize("deg", LATITUDES_DEG)
    def test_geocentric_to_geodetic(self, hayford, deg):
        phi = math.radians(deg)
        psi = geocentric_from_geodetic(phi)
        result = Latitude(psi, LatType.GEOCENTRIC).convert_to(LatType.GEODETIC, hayford)
        assert result.lat_type == LatType.GEODETIC
        assert result.value == pytest.approx(phi, abs=1e-12)

    @pytest.mark.parametrize("deg", LATITUDES_DEG)
    def test_parametric_to_geodetic(self, hayford, deg):
        phi = math.radians(deg)
        beta = parametric_from_geodetic(phi)
        result = Latitude(beta, LatType.PARAMETRIC).convert_to(LatType.GEODETIC, hayford)
        assert result.lat_type == LatType.GEODETIC
        assert result.value == pytest.approx(phi, abs=1e-12)

    @pytest.mark.parametrize("deg", LATITUDES_DEG)
    def test_geocentric_to_parametric(self, hayford, deg):
        phi = math.radians(deg)
        psi = geocentric_from_geodetic(phi)
        result = Latitude(psi, LatType.GEOCENTRIC).convert_to(LatType.PARAMETRIC, hayford)
        assert result.lat_type == LatType.PARAMETRIC
        assert result.value == pytest.approx(parametric_from_geodetic(phi), abs=1e-12)

    @pytest.mark.parametrize("deg", LATITUDES_DEG)
    def test_parametric_to_geocentric(self, hayford, deg):
        phi = math.radians(deg)
        beta = parametric_from_geodetic(phi)
        result = Latitude(beta, LatType.PARAMETRIC).convert_to(LatType.GEOCENTRIC, hayford)
        assert result.lat_type == LatType.GEOCENTRIC
        assert result.value == pytest.approx(geocentric_from_geodetic(phi), abs=1e-12)

    @pytest.mark.parametrize("source", list(LatType))
    @pytest.mark.parametrize("target", list(LatType))
    def test_round_trip(self, hayford, source, target):
        lat = Latitude(math.radians(37.5), source)
        back = lat.convert_to(target, hayford).convert_to(source, hayford)
        assert back.lat_type == source
        assert back.value == pytest.approx(lat.value, abs=1e-12)

    def test_ordering_geocentric_parametric_geodetic(self, hayford):
        # For 0 < phi < 90 deg: geocentric < parametric < geodetic.
        lat = Latitude(math.radians(45.0), LatType.GEODETIC)
        psi = lat.convert_to(LatType.GEOCENTRIC, hayford).value
        beta = lat.convert_to(LatType.PARAMETRIC, hayford).value
        assert psi < beta < lat.value

    @pytest.mark.parametrize("target", list(LatType))
    def test_sphere_has_no_difference(self, target):
        sphere = Ellipsoid("Sphere", 6371000.0, 0.0)
        lat = Latitude(0.7, LatType.GEODETIC)
        assert lat.convert_to(target, sphere).value == pytest.approx(0.7, abs=1e-15)


class TestLatitudeStr:
    def test_format(self):
        assert str(Latitude(math.radians(41.5), LatType.GEODETIC)) == "GEODETIC: 41.500000°"


class TestCartesian3D:
    def test_fields(self):
        p = Cartesian3D(1.0, 2.0, 3.0)
        assert (p.x, p.y, p.z) == (1.0, 2.0, 3.0)

    def test_is_frozen(self):
        with pytest.raises(dataclasses.FrozenInstanceError):
            Cartesian3D(1.0, 2.0, 3.0).x = 0.0


class TestGeodeticToCartesian:
    def test_is_frozen(self):
        with pytest.raises(dataclasses.FrozenInstanceError):
            Geodetic(0.0, 0.0, 0.0).lat = 1.0

    def test_equator_prime_meridian(self, hayford):
        p = Geodetic(0.0, 0.0, 0.0).to_cartesian(hayford)
        assert p.x == pytest.approx(A)
        assert p.y == pytest.approx(0.0, abs=1e-6)
        assert p.z == pytest.approx(0.0, abs=1e-6)

    def test_equator_90_east(self, hayford):
        p = Geodetic(0.0, math.pi / 2, 0.0).to_cartesian(hayford)
        assert p.x == pytest.approx(0.0, abs=1e-6)
        assert p.y == pytest.approx(A)
        assert p.z == pytest.approx(0.0, abs=1e-6)

    def test_north_pole(self, hayford):
        p = Geodetic(math.pi / 2, 0.0, 0.0).to_cartesian(hayford)
        assert p.x == pytest.approx(0.0, abs=1e-6)
        assert p.y == pytest.approx(0.0, abs=1e-6)
        assert p.z == pytest.approx(B, abs=1e-6)

    def test_south_pole(self, hayford):
        p = Geodetic(-math.pi / 2, 0.0, 0.0).to_cartesian(hayford)
        assert p.z == pytest.approx(-B, abs=1e-6)

    def test_height_adds_along_normal(self, hayford):
        h = 1000.0
        p = Geodetic(0.0, 0.0, h).to_cartesian(hayford)
        assert p.x == pytest.approx(A + h)
        q = Geodetic(math.pi / 2, 0.0, h).to_cartesian(hayford)
        assert q.z == pytest.approx(B + h, abs=1e-6)

    def test_surface_point_lies_on_ellipsoid(self, hayford):
        p = Geodetic(math.radians(41.0), math.radians(29.0), 0.0).to_cartesian(hayford)
        assert (p.x**2 + p.y**2) / A**2 + p.z**2 / B**2 == pytest.approx(1.0, abs=1e-12)

    def test_general_point(self, hayford):
        lat, lon, h = math.radians(41.0), math.radians(29.0), 150.0
        n = A / math.sqrt(1 - E_SQ * math.sin(lat) ** 2)
        p = Geodetic(lat, lon, h).to_cartesian(hayford)
        assert p.x == pytest.approx((n + h) * math.cos(lat) * math.cos(lon), abs=1e-6)
        assert p.y == pytest.approx((n + h) * math.cos(lat) * math.sin(lon), abs=1e-6)
        assert p.z == pytest.approx(((1 - E_SQ) * n + h) * math.sin(lat), abs=1e-6)

    def test_sphere(self):
        r = 6371000.0
        sphere = Ellipsoid("Sphere", r, 0.0)
        lat, lon = math.radians(30.0), math.radians(60.0)
        p = Geodetic(lat, lon, 0.0).to_cartesian(sphere)
        assert math.hypot(p.x, p.y, p.z) == pytest.approx(r)
