import pytest

pytest.importorskip("textual")

from geolab import GRS80, HAYFORD, KRASSOVSKY, WGS84, AngleUnit, Ellipsoid  # noqa: E402
from geolab.tui import commands  # noqa: E402
from geolab.tui.commands import (  # noqa: E402
    CommandError,
    completions,
    match_ellipsoid,
    match_unit,
    split_command,
)


class TestSplitCommand:
    def test_name_and_args(self):
        assert split_command("ell grs80") == ("ell", ["grs80"])

    def test_leading_colon_and_whitespace(self):
        assert split_command("  :units   deg ") == ("units", ["deg"])

    def test_name_is_lowercased(self):
        assert split_command("ELL Hayford") == ("ell", ["Hayford"])

    def test_empty(self):
        assert split_command("") == ("", [])
        assert split_command(":") == ("", [])


class TestMatchEllipsoid:
    @pytest.mark.parametrize(
        ("name", "expected"),
        [
            ("WGS84", WGS84),
            ("wgs84", WGS84),
            ("wgs-84", WGS84),
            ("grs 80", GRS80),
            ("hay", HAYFORD),
            ("k", KRASSOVSKY),
            ("bessel1841", commands.ELLIPSOIDS["Bessel 1841"]),
            ("clarke", commands.ELLIPSOIDS["Clarke 1866"]),
        ],
    )
    def test_matches(self, name, expected):
        assert match_ellipsoid(name) is expected

    def test_exact_match_beats_prefix(self, monkeypatch):
        grs8 = Ellipsoid("GRS8", 6378137.0, 298.0)
        monkeypatch.setattr(commands, "ELLIPSOIDS", {"GRS8": grs8, "GRS80": GRS80})
        assert match_ellipsoid("grs8") is grs8

    def test_ambiguous(self, monkeypatch):
        grs67 = Ellipsoid("GRS67", 6378160.0, 298.247167427)
        monkeypatch.setattr(commands, "ELLIPSOIDS", {"GRS67": grs67, "GRS80": GRS80})
        with pytest.raises(CommandError, match="ambiguous"):
            match_ellipsoid("grs")

    def test_unknown(self):
        with pytest.raises(CommandError, match="unknown ellipsoid"):
            match_ellipsoid("airy")

    def test_empty(self):
        with pytest.raises(CommandError):
            match_ellipsoid(" - ")


class TestMatchUnit:
    @pytest.mark.parametrize(
        ("name", "expected"),
        [
            ("dms", AngleUnit.DMS),
            ("DMS", AngleUnit.DMS),
            ("deg", AngleUnit.DEGREES),
            ("degrees", AngleUnit.DEGREES),
            ("rad", AngleUnit.RADIANS),
            ("radians", AngleUnit.RADIANS),
        ],
    )
    def test_matches(self, name, expected):
        assert match_unit(name) is expected

    def test_unknown(self):
        with pytest.raises(CommandError, match="unknown unit"):
            match_unit("grad")


class TestCompletions:
    def test_all_commands(self):
        assert completions("") == ["ell ", "help", "method ", "q", "save", "tab ", "units ", "working"]

    def test_command_prefix(self):
        assert completions("u") == ["units "]
        assert completions(":wo") == ["working"]

    def test_no_match(self):
        assert completions("zz") == []

    def test_ellipsoid_names(self):
        assert completions("ell ") == [f"ell {name}" for name in commands.ELLIPSOIDS]

    def test_ellipsoid_prefix_ignores_case_and_spaces(self):
        assert completions("ell bes") == ["ell Bessel 1841"]
        assert completions("ellipsoid grs") == ["ellipsoid GRS80"]

    def test_units(self):
        assert completions("units ") == ["units dms", "units deg", "units rad"]
        assert completions("units d") == ["units dms", "units deg"]

    def test_method_uses_current_options(self):
        assert completions("method ", ["Iterative", "Bowring"]) == ["method iterative", "method bowring"]
        assert completions("method b", ["Iterative", "Bowring"]) == ["method bowring"]
        assert completions("method ") == []

    def test_tab_numbers(self):
        assert completions("tab ", tab_count=3) == ["tab 1", "tab 2", "tab 3"]

    def test_commands_without_arguments(self):
        assert completions("help x") == []


class TestCustomRegistry:
    def test_match_in_given_registry(self):
        mine = Ellipsoid("Mine", 6378000.0, 300.0)
        assert match_ellipsoid("mi", {"WGS84": WGS84, "Mine": mine}) is mine

    def test_unknown_lists_given_registry(self):
        with pytest.raises(CommandError, match="Mine"):
            match_ellipsoid("zz", {"Mine": Ellipsoid("Mine", 6378000.0, 300.0)})

    def test_completion_with_given_names(self):
        assert completions("ell m", ellipsoid_names=["WGS84", "Mine"]) == ["ell Mine"]
