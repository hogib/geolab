import dataclasses

import pytest

from geolab import (
    BESSEL_1841,
    CLARKE_1866,
    ELLIPSOIDS,
    GRS80,
    HAYFORD,
    KRASSOVSKY,
    WGS84,
    Ellipsoid,
)

HAYFORD_A = 6378388.0
HAYFORD_INV_F = 297.0
HAYFORD_F = 1.0 / 297.0
HAYFORD_B = 6356911.9461279465
HAYFORD_E_SQ = 0.006722670022333322


class TestEllipsoid:
    def test_stores_defining_parameters(self):
        ell = Ellipsoid("Custom", 6378137.0, 298.257223563)
        assert ell.name == "Custom"
        assert ell.a == 6378137.0
        assert ell.inv_f == 298.257223563

    def test_flattening(self):
        ell = Ellipsoid("Custom", 6378137.0, 298.257223563)
        assert ell.f == pytest.approx(1.0 / 298.257223563)

    def test_semi_minor_axis(self):
        ell = Ellipsoid("Custom", 6378137.0, 298.257223563)
        assert ell.b == pytest.approx(6356752.314245, abs=1e-6)

    def test_first_eccentricity_squared(self):
        ell = Ellipsoid("Custom", 6378137.0, 298.257223563)
        assert ell.e_sq == pytest.approx(0.00669437999014, rel=1e-12)

    def test_eccentricity_squared_method_matches_attribute(self):
        ell = Ellipsoid("Custom", 6378137.0, 298.257223563)
        assert ell.eccentricity_squared() == pytest.approx(ell.e_sq)

    def test_derived_fields_are_consistent(self):
        ell = Ellipsoid("Custom", 6378137.0, 298.257223563)
        # e^2 = (a^2 - b^2) / a^2
        assert ell.e_sq == pytest.approx((ell.a**2 - ell.b**2) / ell.a**2)

    def test_zero_inverse_flattening_is_a_sphere(self):
        sphere = Ellipsoid("Sphere", 6371000.0, 0.0)
        assert sphere.f == 0.0
        assert sphere.b == 6371000.0
        assert sphere.e_sq == 0.0
        assert sphere.eccentricity_squared() == 0.0

    def test_derived_fields_are_not_init_arguments(self):
        with pytest.raises(TypeError):
            Ellipsoid("Custom", 6378137.0, 298.257223563, f=0.1)

    @pytest.mark.parametrize("a", [0.0, -1.0, float("nan")])
    def test_rejects_non_positive_a(self, a):
        with pytest.raises(ValueError, match="semi-major"):
            Ellipsoid("Bad", a, 298.0)

    @pytest.mark.parametrize("inv_f", [-298.0, 0.5, 1.0, float("nan")])
    def test_rejects_invalid_inverse_flattening(self, inv_f):
        with pytest.raises(ValueError, match="inverse flattening"):
            Ellipsoid("Bad", 6378137.0, inv_f)

    def test_is_frozen(self):
        ell = Ellipsoid("Custom", 6378137.0, 298.257223563)
        with pytest.raises(dataclasses.FrozenInstanceError):
            ell.a = 1.0


class TestHayford:
    def test_is_an_ellipsoid(self):
        assert isinstance(HAYFORD, Ellipsoid)

    def test_defaults(self):
        ell = HAYFORD
        assert ell.name == "Hayford"
        assert ell.a == HAYFORD_A
        assert ell.inv_f == HAYFORD_INV_F

    def test_derived_values(self):
        ell = HAYFORD
        assert ell.f == pytest.approx(HAYFORD_F)
        assert ell.b == pytest.approx(HAYFORD_B, abs=1e-6)
        assert ell.e_sq == pytest.approx(HAYFORD_E_SQ, rel=1e-12)
        assert ell.eccentricity_squared() == pytest.approx(HAYFORD_E_SQ, rel=1e-12)



class TestDerivedParameters:
    # WGS84 values from NGA TR8350.2, table 3.3.
    def test_second_eccentricity_squared(self):
        assert WGS84.ep_sq == pytest.approx(0.00673949674228, rel=1e-11)

    def test_third_flattening(self):
        assert WGS84.n == pytest.approx((WGS84.a - WGS84.b) / (WGS84.a + WGS84.b))

    def test_mean_radius(self):
        assert WGS84.mean_radius == pytest.approx(6371008.7714, abs=1e-4)

    def test_authalic_radius(self):
        assert WGS84.authalic_radius == pytest.approx(6371007.1810, abs=1e-4)

    def test_volumetric_radius(self):
        assert WGS84.volumetric_radius == pytest.approx(6371000.7900, abs=1e-4)

    def test_sphere_radii_equal_a(self):
        sphere = Ellipsoid("Sphere", 6371000.0, 0.0)
        assert sphere.ep_sq == 0.0
        assert sphere.n == 0.0
        assert sphere.mean_radius == pytest.approx(6371000.0)
        assert sphere.authalic_radius == 6371000.0
        assert sphere.volumetric_radius == pytest.approx(6371000.0)


class TestRegistry:
    @pytest.mark.parametrize(
        ("ell", "b", "e_sq"),
        [
            (WGS84, 6356752.314245, 0.00669437999014),
            (GRS80, 6356752.314140, 0.00669438002290),
            (HAYFORD, 6356911.946128, 0.00672267002233),
            (BESSEL_1841, 6356078.962818, 0.00667437223180),
            (CLARKE_1866, 6356583.800000, 0.00676865799729),
            (KRASSOVSKY, 6356863.018773, 0.00669342162297),
        ],
        ids=lambda v: v.name if isinstance(v, Ellipsoid) else "",
    )
    def test_published_values(self, ell, b, e_sq):
        assert ell.b == pytest.approx(b, abs=1e-6)
        assert ell.e_sq == pytest.approx(e_sq, abs=1e-14)

    def test_lookup_by_name(self):
        assert ELLIPSOIDS["WGS84"] is WGS84
        assert ELLIPSOIDS["Hayford"] is HAYFORD
        assert set(ELLIPSOIDS) == {
            "WGS84", "GRS80", "Hayford", "Bessel 1841", "Clarke 1866", "Krassovsky"
        }

    def test_keys_match_names(self):
        for name, ell in ELLIPSOIDS.items():
            assert ell.name == name
