import math

from ..angles import AngleKind
from ..soldner import geodetic_from_soldner, soldner_direct, soldner_from_geodetic, soldner_inverse
from ..trace import Trace
from .calculator import Calculator
from .widgets import AngleField, ChoiceField, Field, NumberField, ResultRow

LAT, LON = AngleKind.LAT, AngleKind.LON


class SoldnerCalculator(Calculator):
    TITLE = "Soldner coordinates (geodetic rectangular)"

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
                    ("φ, λ → x, y", "to_xy"),
                    ("x, y → φ, λ", "to_geo"),
                    ("Direct: x, y, T, S", "direct"),
                    ("Inverse: two points", "inverse"),
                ],
            ),
            AngleField("lat0", "φ₀ origin latitude", LAT, "40 N"),
            AngleField("lon0", "λ₀ central meridian", LON, "30 E"),
            AngleField("lat", "φ  latitude", LAT, "40 30 N"),
            AngleField("lon", "λ  longitude", LON, "31 E"),
            NumberField("x", "x₁ north (m)", "10000"),
            NumberField("y", "y₁ east (m)", "20000"),
            AngleField("t", "T₁₂ grid bearing", None, "35"),
            NumberField("s", "S  distance (m)", "30000"),
            NumberField("x2", "x₂ north (m)", "35000"),
            NumberField("y2", "y₂ east (m)", "40000"),
        ]

    def calculate(self, trace: Trace) -> list[ResultRow]:
        mode, ell = self.choice("mode"), self.geolab.ellipsoid
        lat0, lon0 = self.number("lat0"), self.number("lon0")
        if mode == "to_xy":
            p = soldner_from_geodetic(self.number("lat"), self.number("lon"), lat0, lon0, ell, trace)
            return [
                self.length_row("x   north", p.x),
                self.length_row("y   east", p.y),
                self.angle_row("γ   meridian convergence", p.convergence),
            ]
        if mode == "to_geo":
            g = geodetic_from_soldner(self.number("x"), self.number("y"), lat0, lon0, ell, trace)
            return [
                self.angle_row("φ   latitude", g.lat, LAT),
                self.angle_row("λ   longitude", g.lon, LON),
                self.angle_row("γ   meridian convergence", g.convergence),
                self.angle_row("φ_f foot-point latitude", g.foot_lat, LAT),
            ]
        if mode == "direct":
            r = soldner_direct(self.number("x"), self.number("y"), self.number("t"), self.number("s"), lat0, ell, trace)
            return [
                self.length_row("x₂  north", r.x2),
                self.length_row("y₂  east", r.y2),
                self.angle_row("T₂₁ back bearing", r.bearing2),
            ]
        r = soldner_inverse(self.number("x"), self.number("y"), self.number("x2"), self.number("y2"), lat0, ell, trace)
        return [
            self.length_row("S   distance", r.distance),
            self.angle_row("T₁₂ bearing", r.bearing1),
            self.angle_row("T₂₁ back bearing", r.bearing2),
        ]
