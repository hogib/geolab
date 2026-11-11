"""Zones (3°, 6° UTM) and the Turkish map-sheet division."""

import math

import pytest

from geolab import GRS80, HAYFORD, WGS84
from geolab.gauss_kruger import krueger_forward
from geolab.map_sheets import sheet_at, sheet_corner_coordinates, sheet_from_name
from geolab.trace import Trace
from geolab.zones import (
    central_meridian,
    geodetic_from_grid,
    grid_from_geodetic,
    split_easting,
    transform_zone,
    zone_of,
)


def rad(d, m=0, s=0):
    return math.radians(d + m / 60 + s / 3600)


class TestZoneNumbers:
    def test_three_degree(self):
        assert zone_of(rad(30.4), 3) == 10 and zone_of(rad(31.6), 3) == 11
        assert math.degrees(central_meridian(10, 3)) == pytest.approx(30)

    @pytest.mark.parametrize("lon0, zone", [(27, 35), (33, 36), (39, 37), (45, 38), (-177, 1), (177, 60), (3, 31)])
    def test_utm(self, lon0, zone):
        # Course: zone = (λ₀ + 3)/6 + 30.
        assert (lon0 + 3) / 6 + 30 == zone
        assert zone_of(rad(lon0 + 2.9), 6) == zone
        assert math.degrees(central_meridian(zone, 6)) == pytest.approx(lon0)

    def test_bad_width(self):
        with pytest.raises(ValueError):
            zone_of(0.5, 4)
        with pytest.raises(ValueError):
            central_meridian(61, 6)


def test_course_easting_example():
    # Sağa değer 598 682.401, yukarı değer 6 010 941.184 in a 3° zone.
    zone, easting = split_easting(598682.401, 10)
    assert easting - 500000 == pytest.approx(98682.401)
    g = geodetic_from_grid(598682.401, 6010941.184, HAYFORD, 3, zone=10)
    assert (g.y, g.x) == pytest.approx((98682.401, 6010941.184))


def test_utm_matches_krueger_scaled():
    lat, lon = rad(39.92), rad(32.85)  # Ankara, UTM zone 36
    p = grid_from_geodetic(lat, lon, GRS80, width=6)
    k = krueger_forward(lat, lon, rad(33), GRS80)
    assert p.zone == 36
    assert p.easting == pytest.approx(500000 + 0.9996 * k.y, abs=1e-3)
    assert p.northing == pytest.approx(0.9996 * k.x, abs=1e-3)
    assert p.scale == pytest.approx(0.9996 * k.scale, abs=1e-9)
    assert p.scale < 1  # inside the two lines of true scale


def test_round_trip_with_prefix():
    p = grid_from_geodetic(rad(40.7), rad(31.2), WGS84, 3)
    assert p.zone == 10 and p.prefixed_easting == pytest.approx(10_000_000 + p.easting)
    g = geodetic_from_grid(p.prefixed_easting, p.northing, WGS84, 3)
    assert (g.lat, g.lon) == pytest.approx((rad(40.7), rad(31.2)), abs=1e-10)
    with pytest.raises(ValueError):
        split_easting(p.prefixed_easting, 11)
    with pytest.raises(ValueError):
        split_easting(p.easting)


def test_zone_to_zone():
    # A point 1.6° east of λ₀ = 30° is in the 1° overlap: compute it in zone 11 too.
    lat, lon = rad(40.6), rad(31.6)
    p10 = grid_from_geodetic(lat, lon, WGS84, 3, zone=10)
    p11 = transform_zone(p10.easting, p10.northing, WGS84, 3, 10, 11)
    direct = grid_from_geodetic(lat, lon, WGS84, 3)
    assert direct.zone == 11
    assert (p11.easting, p11.northing) == pytest.approx((direct.easting, direct.northing), abs=1e-4)
    assert p11.easting < 500000 < p10.easting


class TestSheets:
    def test_course_sheet(self):
        s = sheet_from_name("G25-d-01-d")
        assert s.scale == 5000
        assert math.degrees(s.south) == pytest.approx(40.7) and math.degrees(s.north) == pytest.approx(40.725)
        assert math.degrees(s.west) == pytest.approx(30.5) and math.degrees(s.east) == pytest.approx(30.525)

    def test_g25(self):
        s = sheet_from_name("G-25")
        assert s.name == "G25" and s.scale == 100_000
        assert [math.degrees(v) for v in (s.south, s.north, s.west, s.east)] == pytest.approx([40.5, 41, 30.5, 31])

    @pytest.mark.parametrize(
        "name, scale",
        [("G25", 100_000), ("G25-d", 50_000), ("G25-d-01", 10_000), ("G25-d-01-d", 5000),
         ("G25-d-01-d-1", 2000), ("G25-d-01-d-1-a", 1000), ("G25-d-01-d-1-a-2", 500)],
    )
    def test_names_round_trip(self, name, scale):
        s = sheet_from_name(name)
        assert s.scale == scale
        mid = ((s.south + s.north) / 2, (s.west + s.east) / 2)
        assert sheet_at(*mid, scale).name == name

    def test_numbering(self):
        # 01 is the north-west 3′ cell of G25-d, 05 the north-east, 25 the south-east.
        assert sheet_at(rad(40, 44, 59), rad(30, 30, 1), 10_000).name == "G25-d-01"
        assert sheet_at(rad(40, 44, 59), rad(30, 44, 59), 10_000).name == "G25-d-05"
        assert sheet_at(rad(40, 30, 1), rad(30, 44, 59), 10_000).name == "G25-d-25"
        assert sheet_at(rad(40, 44), rad(30, 50), 50_000).name == "G25-c"
        # 1/2000 is 1 2 / 4 3: the south-east quarter of a 1/5000 sheet is 3.
        assert sheet_at(rad(40, 42, 10), rad(30, 31, 20), 2000).name == "G25-d-01-d-3"

    def test_course_corners(self):
        # The course's printout for G25-d-01-d (Hayford, 3° zone 30° E).
        c = sheet_corner_coordinates(sheet_from_name("G25-d-01-d"), HAYFORD, 3)
        printed = {
            "A": (542259.752, 4507456.400),
            "B": (542243.948, 4510232.706),
            "C": (544356.154, 4510245.033),
            "D": (544372.748, 4507468.725),
        }
        for k, (e, n) in printed.items():
            assert c[k].zone == 10
            assert c[k].easting == pytest.approx(e, abs=3e-3)
            assert c[k].northing == pytest.approx(n, abs=1e-2)

    @pytest.mark.parametrize("bad", ["Q25", "G14", "G25-e", "G25-d-26", "G25-d-01-d-5", "", "G25-d-01-d-1-a-2-1"])
    def test_bad_names(self, bad):
        with pytest.raises(ValueError):
            sheet_from_name(bad)

    def test_outside_grid(self):
        with pytest.raises(ValueError):
            sheet_at(rad(48), rad(30), 5000)


def test_trace():
    trace = Trace()
    p = grid_from_geodetic(rad(40.7), rad(31.2), WGS84, 6, trace=trace)
    assert trace["sağa"] == p.easting and trace["m₀"] == 0.9996
