import math

import pytest

from geolab.angles import (
    AngleKind,
    AngleParseError,
    AngleUnit,
    format_angle,
    format_degrees,
    format_dms,
    format_radians,
    parse_angle,
    to_dms,
)

LAT, LON = AngleKind.LAT, AngleKind.LON


def deg(text, **kwargs):
    return math.degrees(parse_angle(text, **kwargs))


class TestParseDecimal:
    @pytest.mark.parametrize(
        ("text", "expected"),
        [("41.5", 41.5), ("-29.25", -29.25), ("+12", 12.0), ("0", 0.0), (".5", 0.5), ("41.", 41.0)],
    )
    def test_decimal_degrees(self, text, expected):
        assert deg(text) == pytest.approx(expected)

    def test_returns_radians(self):
        assert parse_angle("180") == pytest.approx(math.pi)

    def test_surrounding_whitespace(self):
        assert deg("  41.5  ") == pytest.approx(41.5)

    def test_radians_unit(self):
        assert parse_angle("0.5", unit=AngleUnit.RADIANS) == 0.5
        assert parse_angle("-1.2", unit=AngleUnit.RADIANS) == -1.2

    def test_degree_symbol_overrides_radians_unit(self):
        assert deg("30°", unit=AngleUnit.RADIANS) == pytest.approx(30.0)

    def test_dms_unit_reads_bare_number_as_degrees(self):
        assert deg("41.5", unit=AngleUnit.DMS) == pytest.approx(41.5)


class TestParseDMS:
    EXPECTED = 41 + 0 / 60 + 30.5 / 3600

    @pytest.mark.parametrize(
        "text",
        [
            "41°00'30.5\"",
            "41° 00' 30.5\"",
            "41°00′30.5″",
            "41°00'30.5''",
            "41°00’30.5\"",
            "41 00 30.5",
            "41 0 30.5",
            "41:00:30.5",
        ],
    )
    def test_formats(self, text):
        assert deg(text) == pytest.approx(self.EXPECTED)

    def test_degrees_and_minutes_only(self):
        assert deg("41°30'") == pytest.approx(41.5)
        assert deg("41 30") == pytest.approx(41.5)

    def test_decimal_minutes(self):
        assert deg("41 30.5") == pytest.approx(41 + 30.5 / 60)

    def test_negative(self):
        assert deg("-41 30 0") == pytest.approx(-41.5)

    def test_negative_zero_degrees_keeps_sign(self):
        assert deg("-0 30 0") == pytest.approx(-0.5)
        assert deg("-0°30'") == pytest.approx(-0.5)


class TestParseHemisphere:
    @pytest.mark.parametrize(
        ("text", "expected"),
        [
            ("41°30'N", 41.5),
            ("41°30'S", -41.5),
            ("29 15 E", 29.25),
            ("29 15 W", -29.25),
            ("N 41 30", 41.5),
            ("s41.5", -41.5),
            ("0 30 S", -0.5),
        ],
    )
    def test_letters(self, text, expected):
        assert deg(text) == pytest.approx(expected)

    def test_letter_must_match_kind(self):
        with pytest.raises(AngleParseError, match="hemisphere"):
            parse_angle("41 30 E", kind=LAT)
        with pytest.raises(AngleParseError, match="hemisphere"):
            parse_angle("29 N", kind=LON)

    def test_sign_and_letter_conflict(self):
        with pytest.raises(AngleParseError, match="sign"):
            parse_angle("-41 30 N")

    def test_letter_implies_range_check(self):
        with pytest.raises(AngleParseError, match="latitude"):
            parse_angle("91 N")


class TestParseRange:
    def test_latitude_limits(self):
        assert deg("90", kind=LAT) == pytest.approx(90)
        assert deg("-90", kind=LAT) == pytest.approx(-90)
        with pytest.raises(AngleParseError):
            parse_angle("90.0001", kind=LAT)

    def test_longitude_limits(self):
        assert deg("-180", kind=LON) == pytest.approx(-180)
        assert deg("359.5", kind=LON) == pytest.approx(359.5)
        with pytest.raises(AngleParseError):
            parse_angle("-180.5", kind=LON)
        with pytest.raises(AngleParseError):
            parse_angle("361", kind=LON)

    def test_radians_are_range_checked(self):
        with pytest.raises(AngleParseError):
            parse_angle("2", kind=LAT, unit=AngleUnit.RADIANS)

    def test_no_kind_no_range_check(self):
        assert deg("400") == pytest.approx(400)


class TestParseErrors:
    @pytest.mark.parametrize(
        "text",
        [
            "",
            "   ",
            "abc",
            "N",
            "-",
            "41 60 0",  # minutes must be < 60
            "41 30 60",  # seconds must be < 60
            "41.5 30",  # fractional degrees followed by minutes
            "41 30.5 10",  # fractional minutes followed by seconds
            "41 30 10 5",  # too many components
            "41'30°",  # symbols out of order
            "41°30\"",  # seconds without minutes
            "41°30'x",
            "1e5",
        ],
    )
    def test_rejected(self, text):
        with pytest.raises(AngleParseError):
            parse_angle(text)

    def test_is_value_error(self):
        assert issubclass(AngleParseError, ValueError)


class TestToDMS:
    def test_split(self):
        assert to_dms(math.radians(41 + 30 / 60 + 15.25 / 3600)) == (1, 41, 30, 15.25)

    def test_negative(self):
        sign, d, m, s = to_dms(math.radians(-41.5))
        assert (sign, d, m, s) == (-1, 41, 30, 0.0)

    def test_negative_below_one_degree(self):
        assert to_dms(math.radians(-0.5)) == (-1, 0, 30, 0.0)

    def test_rounding_carries_into_minutes_and_degrees(self):
        almost_42 = math.radians(41 + 59 / 60 + 59.99999 / 3600)
        assert to_dms(almost_42, sec_decimals=4) == (1, 42, 0, 0.0)

    def test_tiny_negative_rounds_to_positive_zero(self):
        assert to_dms(-1e-15) == (1, 0, 0, 0.0)


class TestFormat:
    def test_dms_plain(self):
        assert format_dms(math.radians(41 + 30.5 / 3600)) == "41°00'30.5000\""

    def test_dms_negative(self):
        assert format_dms(math.radians(-0.5)) == "-0°30'00.0000\""

    def test_dms_hemispheres(self):
        assert format_dms(math.radians(41.5), LAT) == "41°30'00.0000\"N"
        assert format_dms(math.radians(-41.5), LAT) == "41°30'00.0000\"S"
        assert format_dms(math.radians(29.25), LON) == "29°15'00.0000\"E"
        assert format_dms(math.radians(-29.25), LON) == "29°15'00.0000\"W"
        assert format_dms(0.0, LAT) == "0°00'00.0000\"N"

    def test_dms_decimals(self):
        rad = math.radians(41 + 5.123456 / 3600)
        assert format_dms(rad, sec_decimals=2) == "41°00'05.12\""
        assert format_dms(rad, sec_decimals=0) == "41°00'05\""

    def test_dms_carry(self):
        assert format_dms(math.radians(41 + 59 / 60 + 59.99999 / 3600)) == "42°00'00.0000\""

    def test_degrees(self):
        assert format_degrees(math.radians(41.5)) == "41.50000000°"
        assert format_degrees(math.radians(41.5), 2) == "41.50°"

    def test_degrees_no_negative_zero(self):
        assert format_degrees(-1e-15) == "0.00000000°"

    def test_radians(self):
        assert format_radians(0.5) == "0.5000000000 rad"
        assert format_radians(-1e-15) == "0.0000000000 rad"

    def test_format_angle_dispatch(self):
        rad = math.radians(41.5)
        assert format_angle(rad, AngleUnit.DMS, LAT) == "41°30'00.0000\"N"
        assert format_angle(rad, AngleUnit.DEGREES) == "41.50000000°"
        assert format_angle(rad, AngleUnit.RADIANS) == format_radians(rad)


class TestRoundTrip:
    @pytest.mark.parametrize("value", [-179.999, -90.0, -41.508472, -0.25, 0.0, 0.0001, 29.0, 89.99999])
    def test_dms_round_trip(self, value):
        rad = math.radians(value)
        assert parse_angle(format_dms(rad)) == pytest.approx(rad, abs=math.radians(0.00005 / 3600))

    @pytest.mark.parametrize("value", [-41.5, 0.25, 89.99])
    def test_hemisphere_round_trip(self, value):
        rad = math.radians(value)
        text = format_dms(rad, LAT)
        assert parse_angle(text, kind=LAT) == pytest.approx(rad, abs=math.radians(0.00005 / 3600))

    def test_degrees_round_trip(self):
        rad = math.radians(-41.508472)
        assert parse_angle(format_degrees(rad)) == pytest.approx(rad, abs=1e-10)
