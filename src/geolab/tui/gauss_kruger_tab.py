import math

from ..angles import AngleKind
from ..gauss_kruger import (
    gauss_kruger_from_geodetic,
    geodetic_from_gauss_kruger,
    gk_direct,
    gk_inverse,
    krueger_forward,
    krueger_inverse,
)
from ..trace import Trace
from .calculator import Calculator
from .widgets import AngleField, ChoiceField, Field, NumberField, ResultRow

LAT, LON = AngleKind.LAT, AngleKind.LON
ARCSEC = math.pi / (180 * 3600)


class GaussKrugerCalculator(Calculator):
    """The course's Gauss–Krüger series, compared with Krüger's n⁶ series."""

    TITLE = "Gauss–Krüger (transverse Mercator)"

    MODE_FIELDS = {
        "to_xy": {"lat", "lon"},
        "to_geo": {"x", "y"},
        "direct": {"x", "y", "t", "s"},
        "inverse": {"x", "y", "x2", "y2"},
    }

    def fields(self) -> list[Field]:
        return [
            ChoiceField(
                "mode",
                "mode",
                [
                    ("φ, λ → X, Y", "to_xy"),
                    ("X, Y → φ, λ", "to_geo"),
                    ("Direct: X, Y, T, S", "direct"),
                    ("Inverse: two points", "inverse"),
                ],
            ),
            AngleField("lon0", "λ₀ central meridian", LON, "30 E"),
            AngleField("lat", "φ  latitude", LAT, "40 42 N"),
            AngleField("lon", "λ  longitude", LON, "30 30 E"),
            NumberField("x", "X₁ north (m)", "4507456.39"),
            NumberField("y", "Y₁ east of λ₀ (m)", "42259.75"),
            AngleField("t", "T₁₂ grid bearing", None, "35"),
            NumberField("s", "S  geodesic length (m)", "30000"),
            NumberField("x2", "X₂ north (m)", "4532000"),
            NumberField("y2", "Y₂ east of λ₀ (m)", "60000"),
        ]

    def calculate(self, trace: Trace) -> list[ResultRow]:
        mode, ell = self.choice("mode"), self.geolab.ellipsoid
        lon0 = self.number("lon0")
        if mode == "to_xy":
            lat, lon = self.number("lat"), self.number("lon")
            p = gauss_kruger_from_geodetic(lat, lon, lon0, ell, trace=trace)
            k = krueger_forward(lat, lon, lon0, ell)
            return [
                self.length_row("X   north", p.x),
                self.length_row("Y   east", p.y),
                self.value_row("m   length distortion", p.scale),
                self.value_row("m²  area distortion", p.scale**2),
                self.angle_row("c   grid convergence", p.convergence),
                self.value_row("X − Krüger (mm)", (p.x - k.x) * 1000),
                self.value_row("Y − Krüger (mm)", (p.y - k.y) * 1000),
            ]
        x, y = self.number("x"), self.number("y")
        if mode == "to_geo":
            g = geodetic_from_gauss_kruger(x, y, lon0, ell, trace)
            k_lat, k_lon = krueger_inverse(x, y, lon0, ell)
            return [
                self.angle_row("φ   latitude", g.lat, LAT),
                self.angle_row("λ   longitude", g.lon, LON),
                self.value_row("m   length distortion", g.scale),
                self.angle_row("c   grid convergence", g.convergence),
                self.angle_row("φ_f foot-point latitude", g.foot_lat, LAT),
                self.value_row("φ − Krüger (″)", (g.lat - k_lat) / ARCSEC),
                self.value_row("λ − Krüger (″)", math.remainder(g.lon - k_lon, math.tau) / ARCSEC),
            ]
        if mode == "direct":
            r = gk_direct(x, y, self.number("t"), self.number("s"), ell, trace)
            return [
                self.length_row("X₂  north", r.x2),
                self.length_row("Y₂  east", r.y2),
                self.angle_row("T₂₁ back grid bearing", r.bearing2),
            ]
        r = gk_inverse(x, y, self.number("x2"), self.number("y2"), ell, trace)
        return [
            self.length_row("S   geodesic length", r.distance),
            self.angle_row("T₁₂ grid bearing", r.bearing1),
            self.angle_row("T₂₁ back grid bearing", r.bearing2),
        ]
