import math

from ..angles import AngleKind
from ..arcs import (
    latitude_from_meridian_arc,
    latitude_from_meridian_distance,
    longitude_from_parallel_arc,
    meridian_arc_length,
    meridian_distance,
    parallel_arc_length,
)
from ..curvature import (
    azimuth_radius,
    gaussian_mean_radius,
    meridian_radius,
    parallel_radius,
    prime_vertical_radius,
)
from ..ellipse import ellipse_arc, parametric_from_normal, parametric_from_polar
from ..ellipsoids import Ellipsoid
from ..geodesic import direct, inverse
from ..surface import fundamental_forms
from ..trace import Trace
from .calculator import Calculator
from .widgets import AngleField, ChoiceField, Field, NumberField, ResultRow, TextField

LAT, LON = AngleKind.LAT, AngleKind.LON


class InverseProblem(Calculator):
    TITLE = "Inverse problem (Vincenty)"

    def fields(self) -> list[Field]:
        return [
            AngleField("lat1", "φ₁ latitude", LAT, "41 00 30 N"),
            AngleField("lon1", "λ₁ longitude", LON, "29 00 00 E"),
            AngleField("lat2", "φ₂ latitude", LAT, "39 55 00 N"),
            AngleField("lon2", "λ₂ longitude", LON, "32 51 00 E"),
        ]

    def calculate(self, trace: Trace) -> list[ResultRow]:
        r = inverse(
            self.number("lat1"),
            self.number("lon1"),
            self.number("lat2"),
            self.number("lon2"),
            self.geolab.ellipsoid,
            trace,
        )
        return [
            self.length_row("s   distance", r.distance),
            self.angle_row("α₁  azimuth at P₁", r.azimuth1),
            self.angle_row("α₂  azimuth at P₂", r.azimuth2),
            self.angle_row("α₂₁ back azimuth", r.back_azimuth),
            self.value_row("iterations", r.iterations),
        ]


class DirectProblem(Calculator):
    TITLE = "Direct problem (Vincenty)"

    def fields(self) -> list[Field]:
        return [
            AngleField("lat1", "φ₁ latitude", LAT, "41 00 30 N"),
            AngleField("lon1", "λ₁ longitude", LON, "29 00 00 E"),
            AngleField("az1", "α₁ azimuth", None, "109 06 35.4977"),
            NumberField("s", "s  distance (m)", "348273.5824"),
        ]

    def calculate(self, trace: Trace) -> list[ResultRow]:
        r = direct(
            self.number("lat1"),
            self.number("lon1"),
            self.number("az1"),
            self.number("s"),
            self.geolab.ellipsoid,
            trace,
        )
        return [
            self.angle_row("φ₂  latitude", r.lat2, LAT),
            self.angle_row("λ₂  longitude", r.lon2, LON),
            self.angle_row("α₂  azimuth at P₂", r.azimuth2),
            self.angle_row("α₂₁ back azimuth", r.back_azimuth),
            self.value_row("iterations", r.iterations),
        ]


class RadiiCalculator(Calculator):
    TITLE = "Radii of curvature and fundamental forms"

    def fields(self) -> list[Field]:
        return [
            AngleField("lat", "φ  latitude", LAT, "41 00 30 N"),
            AngleField("az", "α  azimuth", None, "45"),
        ]

    def calculate(self, trace: Trace) -> list[ResultRow]:
        lat, az, ell = self.number("lat"), self.number("az"), self.geolab.ellipsoid
        trace.section("Meridian radius M")
        m = meridian_radius(lat, ell, trace)
        trace.section("Prime vertical radius N")
        n = prime_vertical_radius(lat, ell, trace)
        trace.section("Parallel radius r")
        r = parallel_radius(lat, ell, trace)
        trace.section("Gaussian mean radius R")
        gauss = gaussian_mean_radius(lat, ell, trace)
        trace.section("Radius in azimuth α")
        r_az = azimuth_radius(lat, az, ell, trace)
        ff = fundamental_forms(lat, ell, trace)
        return [
            self.length_row("M   meridian", m),
            self.length_row("N   prime vertical", n),
            self.length_row("r   parallel", r),
            self.length_row("R   Gaussian mean", gauss),
            self.length_row("Rα  in azimuth α", r_az),
            self.value_row("E   1st form (= M²)", ff.E),
            self.value_row("F   1st form", ff.F),
            self.value_row("G   1st form (= N² cos²φ)", ff.G),
            self.value_row("L   2nd form (= M)", ff.L),
            self.value_row("M̄   2nd form", ff.M),
            self.value_row("N̄   2nd form (= N cos²φ)", ff.N),
            self.value_row("K   Gauss curvature 1/(MN)", ff.gauss_curvature),
            self.value_row("H   mean curvature", ff.mean_curvature),
        ]


class ArcCalculator(Calculator):
    TITLE = "Meridian and parallel arcs"

    MODE_FIELDS = {
        "lengths": {"lon1", "lat2", "lon2"},
        "meridian": {"dm"},
        "parallel": {"lon1", "s"},
    }

    def fields(self) -> list[Field]:
        return [
            ChoiceField(
                "mode",
                "mode",
                [
                    ("Lengths: φ, λ → arcs", "lengths"),
                    ("Meridian: φ₁, Δm → φ₂", "meridian"),
                    ("Parallel: φ, λ₁, s → λ₂", "parallel"),
                ],
            ),
            AngleField("lat1", "φ₁ latitude", LAT, "41 00 30 N"),
            AngleField("lon1", "λ₁ longitude", LON, "29 00 00 E"),
            AngleField("lat2", "φ₂ latitude", LAT, "39 55 00 N"),
            AngleField("lon2", "λ₂ longitude", LON, "32 51 00 E"),
            NumberField("dm", "Δm meridian arc (m, + north)", "-120000"),
            NumberField("s", "s  parallel arc (m, + east)", "325000"),
        ]

    def calculate(self, trace: Trace) -> list[ResultRow]:
        mode = self.choice("mode")
        if mode == "meridian":
            return self._latitude_from_arc(trace)
        if mode == "parallel":
            return self._longitude_from_arc(trace)
        lat1, lon1 = self.number("lat1"), self.number("lon1")
        lat2, lon2 = self.number("lat2"), self.number("lon2")
        ell = self.geolab.ellipsoid
        meridian = meridian_arc_length(lat1, lat2, ell, trace)
        trace.section("Parallel arc at φ₁")
        parallel1 = parallel_arc_length(lat1, lon1, lon2, ell, trace)
        trace.section("Parallel arc at φ₂")
        parallel2 = parallel_arc_length(lat2, lon1, lon2, ell, trace)
        return [
            self.length_row("Δm  meridian arc φ₁→φ₂", meridian),
            self.length_row("s₁  parallel arc at φ₁", parallel1),
            self.length_row("s₂  parallel arc at φ₂", parallel2),
            self.length_row("m₁  equator → φ₁", meridian_distance(lat1, ell)),
            self.length_row("m₂  equator → φ₂", meridian_distance(lat2, ell)),
        ]

    def _latitude_from_arc(self, trace: Trace) -> list[ResultRow]:
        lat1, ell = self.number("lat1"), self.geolab.ellipsoid
        p = latitude_from_meridian_arc(lat1, self.number("dm"), ell, trace)
        return [
            self.angle_row("φ₂  latitude", p.lat, LAT),
            self.length_row("m₁  equator → φ₁", meridian_distance(lat1, ell)),
            self.length_row("m₂  equator → φ₂", meridian_distance(p.lat, ell)),
            self.angle_row("β₂  reduced latitude", p.reduced_lat, LAT),
            self.length_row("x₂  meridian ellipse", p.x),
            self.length_row("y₂  meridian ellipse", p.y),
        ]

    def _longitude_from_arc(self, trace: Trace) -> list[ResultRow]:
        lat, lon1 = self.number("lat1"), self.number("lon1")
        ell = self.geolab.ellipsoid
        trace.section("Along the parallel at φ")
        lon2 = longitude_from_parallel_arc(lat, lon1, self.number("s"), ell, trace)
        return [
            self.angle_row("λ₂  longitude", lon2, LON),
            self.angle_row("Δλ", math.remainder(lon2 - lon1, math.tau)),
            self.length_row("r   parallel radius", parallel_radius(lat, ell)),
        ]


class EllipseArcCalculator(Calculator):
    """Arc length between two angles, or the point at a given meridian arc length."""

    TITLE = "2D ellipse"

    ANGLE_TYPES = {
        "parametric": ("t", None, ""),
        "polar": ("θ", parametric_from_polar, "atan2(a sin θ, b cos θ)"),
        "normal": ("φ", parametric_from_normal, "atan2(b sin φ, a cos φ)"),
    }
    MODE_FIELDS = {
        "arc": {"t1", "t2", "angle"},
        "point": {"m", "lon"},
    }

    def fields(self) -> list[Field]:
        return [
            ChoiceField("mode", "mode", [("Arc: angles → s", "arc"), ("Point: m → φ, x, y", "point")]),
            NumberField("a", "a  semi-axis along x", "5"),
            NumberField("b", "b  semi-axis along y", "3"),
            AngleField("t1", "angle 1", None, "0"),
            AngleField("t2", "angle 2", None, "90"),
            ChoiceField(
                "angle",
                "angle type",
                [("Parametric t (β)", "parametric"), ("Polar θ (ψ)", "polar"), ("Normal φ (geodetic)", "normal")],
            ),
            NumberField("m", "m  arc from equator", "4540000"),
            AngleField("lon", "λ  longitude", LON, "0"),
            ChoiceField("axes", "a, b from", [("Custom", "custom"), ("Ellipsoid", "ellipsoid")]),
        ]

    def on_choice_field_changed(self, event: ChoiceField.Changed) -> None:
        if event.field.key == "mode":
            self.show_mode_fields()
            axes = self.field("axes")
            if self.choice("mode") == "point" and isinstance(axes, ChoiceField):
                axes.select(1)  # the meridian of the current ellipsoid is the usual case
        if event.field.key == "axes":
            self._fill_axes()
        self.recalculate()

    def ellipsoid_changed(self) -> None:
        # Picking an ellipsoid (:ell, e) means "use it": switch away from custom axes.
        axes = self.field("axes")
        assert isinstance(axes, ChoiceField)
        if axes.value != "ellipsoid":
            axes.select(1)  # posts Changed, which fills a, b and recalculates
            return
        self._fill_axes()
        self.recalculate()

    def _fill_axes(self) -> None:
        """In ellipsoid mode, a and b are the current ellipsoid's (its meridian ellipse)."""
        if self.choice("axes") != "ellipsoid":
            return
        ell = self.geolab.ellipsoid
        for key, value in (("a", ell.a), ("b", ell.b)):
            field = self.field(key)
            assert isinstance(field, TextField)
            field.value = repr(value)

    def calculate(self, trace: Trace) -> list[ResultRow]:
        a, b = self.number("a"), self.number("b")
        if self.choice("mode") == "point":
            rows = self._point(a, b, trace)
        else:
            rows = self._arc(a, b, trace)
        return rows + self._parameter_rows(a, b, trace)

    def _arc(self, a: float, b: float, trace: Trace) -> list[ResultRow]:
        t1, t2 = self.number("t1"), self.number("t2")
        symbol, convert, formula = self.ANGLE_TYPES[self.choice("angle")]
        if convert is not None:
            trace.section(f"{symbol} → parametric t")
            t1 = trace.step("t₁", convert(t1, a, b), formula.replace(symbol, f"{symbol}₁"), unit="rad")
            t2 = trace.step("t₂", convert(t2, a, b), formula.replace(symbol, f"{symbol}₂"), unit="rad")
        r = ellipse_arc(a, b, t1, t2, trace)
        return [
            self.value_row("s   arc length", r.length),
            self.value_row("c   chord", r.chord),
            self.value_row("x₁", r.p1[0]),
            self.value_row("y₁", r.p1[1]),
            self.value_row("x₂", r.p2[0]),
            self.value_row("y₂", r.p2[1]),
            self.value_row("P   perimeter", r.perimeter),
        ]

    def _point(self, a: float, b: float, trace: Trace) -> list[ResultRow]:
        if not (a > 0 and b > 0):
            raise ValueError("semi-axes a and b must be positive")
        if b > a:
            raise ValueError("point mode needs a ≥ b (a is the equatorial semi-axis)")
        ell = self.geolab.ellipsoid
        if (a, b) != (ell.a, ell.b):
            ell = Ellipsoid("custom", a, a / (a - b) if a != b else 0)
        p = latitude_from_meridian_distance(self.number("m"), ell, trace)
        lon = self.number("lon")
        big_x, big_y, big_z = p.x * math.cos(lon), p.x * math.sin(lon), p.y
        trace.section("ECEF at longitude λ")
        trace.step("X", big_x, "x cos λ", unit="m")
        trace.step("Y", big_y, "x sin λ", unit="m")
        trace.step("Z", big_z, "y", unit="m")
        return [
            self.angle_row("φ   geodetic latitude", p.lat, LAT),
            self.angle_row("β   reduced latitude", p.reduced_lat, LAT),
            self.angle_row("ψ   geocentric latitude", p.geocentric_lat, LAT),
            self.value_row("x   meridian ellipse", p.x),
            self.value_row("y   meridian ellipse", p.y),
            self.value_row("X   ECEF", big_x),
            self.value_row("Y   ECEF", big_y),
            self.value_row("Z   ECEF", big_z),
        ]

    def _parameter_rows(self, a: float, b: float, trace: Trace) -> list[ResultRow]:
        big, small = max(a, b), min(a, b)
        f = (big - small) / big
        e_sq = 1 - (small / big) ** 2
        ep_sq = (big / small) ** 2 - 1
        n = (big - small) / (big + small)
        c = big**2 / small
        trace.section("Ellipse parameters")
        trace.step("f", f, "(a − b) / a")
        trace.step("e²", e_sq, "(a² − b²) / a²")
        trace.step("e′²", ep_sq, "(a² − b²) / b²")
        trace.step("n", n, "(a − b) / (a + b)")
        trace.step("c", c, "a² / b", unit="m")
        return [
            self.value_row("f   flattening", f),
            self.value_row("1/f", 1 / f if f else 0.0),
            self.value_row("e²  1st eccentricity²", e_sq),
            self.value_row("e′² 2nd eccentricity²", ep_sq),
            self.value_row("n   3rd flattening", n),
            self.value_row("c   polar radius a²/b", c),
        ]
