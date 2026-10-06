import math

import pytest

from geolab import WGS84, geodesic_points, inverse
from geolab.landmask import COLUMNS, ROWS, _bits, is_land


class TestLandMask:
    def test_size(self):
        assert len(_bits()) == COLUMNS * ROWS // 8

    @pytest.mark.parametrize(
        ("place", "lat", "lon"),
        [
            ("Anatolia", 39.0, 33.0),
            ("Sahara", 23.0, 10.0),
            ("Amazon", -5.0, -60.0),
            ("Siberia", 62.0, 100.0),
            ("Australia", -25.0, 134.0),
            ("Antarctica", -85.0, 0.0),
            ("Greenland", 72.0, -40.0),
        ],
    )
    def test_land(self, place, lat, lon):
        assert is_land(lat, lon), place

    @pytest.mark.parametrize(
        ("place", "lat", "lon"),
        [
            ("Atlantic", 30.0, -40.0),
            ("Pacific", 0.0, -150.0),
            ("Indian Ocean", -20.0, 80.0),
            ("Caspian Sea", 42.0, 50.5),
            ("North Pole", 89.9, 0.0),
        ],
    )
    def test_water(self, place, lat, lon):
        assert not is_land(lat, lon), place

    def test_longitude_wraps(self):
        assert is_land(39.0, 33.0 + 360) == is_land(39.0, 33.0)
        assert is_land(-25.0, 134.0 - 360)

    def test_outside_latitude_range(self):
        assert not is_land(91.0, 0.0)
        assert not is_land(-91.0, 0.0)

    def test_edges(self):
        is_land(90.0, 180.0)
        is_land(-90.0, -180.0)

    def test_land_fraction_is_plausible(self):
        land = sum(bin(b).count("1") for b in _bits()) / (COLUMNS * ROWS)
        assert 0.25 < land < 0.40  # per cell, so polar land is over-weighted


class TestGeodesicPoints:
    def test_endpoints_and_count(self):
        lat1, lon1, lat2, lon2 = map(math.radians, (41.0, 29.0, -33.9, 151.2))
        r = inverse(lat1, lon1, lat2, lon2, WGS84)
        points = geodesic_points(lat1, lon1, r.azimuth1, r.distance, WGS84, count=50)
        assert len(points) == 51
        assert points[0] == (lat1, lon1)
        assert points[-1] == pytest.approx((lat2, lon2), abs=1e-10)

    def test_evenly_spaced(self):
        lat1, lon1 = math.radians(10), math.radians(20)
        points = geodesic_points(lat1, lon1, math.radians(60), 1e6, WGS84, count=4)
        steps = [inverse(*a, *b, WGS84).distance for a, b in zip(points, points[1:])]
        assert steps == pytest.approx([2.5e5] * 4, abs=1e-3)

    def test_longitudes_normalised(self):
        points = geodesic_points(0.0, math.radians(179), math.pi / 2, 5e5, WGS84, count=10)
        assert all(-math.pi <= lon <= math.pi for _, lon in points)

    def test_count_must_be_positive(self):
        with pytest.raises(ValueError):
            geodesic_points(0.0, 0.0, 0.0, 1.0, WGS84, count=0)
