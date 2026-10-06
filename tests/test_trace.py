import math

import pytest

from geolab import HAYFORD, Geodetic, Latitude, LatType
from geolab.trace import Section, Step, Trace


class TestTrace:
    def test_step_records_and_returns_value(self):
        trace = Trace()
        assert trace.step("N", 1.5, "a / W") == 1.5
        assert trace.entries == [Step("N", 1.5, "a / W")]

    def test_sections_are_kept_in_order(self):
        trace = Trace()
        trace.section("Iteration 1")
        trace.step("φ", 0.1)
        trace.section("Iteration 2")
        trace.step("φ", 0.2)
        assert trace.entries == [
            Section("Iteration 1"),
            Step("φ", 0.1),
            Section("Iteration 2"),
            Step("φ", 0.2),
        ]
        assert trace.steps == [Step("φ", 0.1), Step("φ", 0.2)]

    def test_lookup_returns_latest_step(self):
        trace = Trace()
        trace.step("φ", 0.1)
        trace.step("φ", 0.2)
        assert trace["φ"] == 0.2

    def test_lookup_missing(self):
        with pytest.raises(KeyError):
            Trace()["φ"]

    def test_new_traces_do_not_share_entries(self):
        a, b = Trace(), Trace()
        a.step("x", 1.0)
        assert b.entries == []


class TestTracedCalculations:
    def test_to_cartesian(self):
        trace = Trace()
        point = Geodetic(math.radians(41.0), math.radians(29.0), 150.0)
        result = point.to_cartesian(HAYFORD, trace)
        assert [s.name for s in trace.steps] == ["e²", "N", "X", "Y", "Z"]
        assert trace["X"] == result.x
        assert trace["Y"] == result.y
        assert trace["Z"] == result.z
        assert trace["N"] == pytest.approx(
            HAYFORD.a / math.sqrt(1 - HAYFORD.e_sq * math.sin(point.lat) ** 2)
        )

    def test_trace_does_not_change_result(self):
        point = Geodetic(math.radians(41.0), math.radians(29.0), 150.0)
        assert point.to_cartesian(HAYFORD, Trace()) == point.to_cartesian(HAYFORD)

    @pytest.mark.parametrize(
        ("source", "target", "names"),
        [
            (LatType.GEODETIC, LatType.GEOCENTRIC, ["e²", "ψ"]),
            (LatType.GEODETIC, LatType.PARAMETRIC, ["e²", "β"]),
            (LatType.GEOCENTRIC, LatType.GEODETIC, ["e²", "φ"]),
            (LatType.PARAMETRIC, LatType.GEOCENTRIC, ["e²", "φ", "ψ"]),
        ],
    )
    def test_latitude_conversion(self, source, target, names):
        trace = Trace()
        result = Latitude(0.7, source).convert_to(target, HAYFORD, trace)
        assert [s.name for s in trace.steps] == names
        assert trace.steps[-1].value == result.value

    def test_same_type_conversion_records_nothing(self):
        trace = Trace()
        Latitude(0.7, LatType.GEODETIC).convert_to(LatType.GEODETIC, HAYFORD, trace)
        assert trace.entries == []


class TestUnits:
    def test_default_is_dimensionless(self):
        trace = Trace()
        trace.step("e²", 0.0067)
        assert trace.steps[0].unit == ""

    def test_calculations_tag_angles_and_lengths(self):
        trace = Trace()
        Geodetic(0.7, 0.5, 100.0).to_cartesian(HAYFORD).to_geodetic(HAYFORD, trace=trace)
        units = {s.name: s.unit for s in trace.steps}
        assert units["φ"] == "rad"
        assert units["λ"] == "rad"
        assert units["N"] == "m"
        assert units["h"] == "m"
