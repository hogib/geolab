import math

from ..angles import AngleKind
from ..projection import (
    isometric_difference,
    isometric_latitude,
    latitude_difference,
    latitude_from_isometric,
    mercator_sphere,
    mercator_sphere_inverse,
    tissot,
)
from ..trace import Trace
from .calculator import Calculator
from .widgets import AngleField, ChoiceField, Field, NumberField, ResultRow

LAT, LON = AngleKind.LAT, AngleKind.LON
ARCSEC = math.pi / (180 * 3600)


class ProjectionCalculator(Calculator):
    TITLE = "Projection theory: Tissot, isometric latitude, Mercator"

    MODE_FIELDS = {
        "tissot": {"m1", "m2", "t"},
        "to_q": {"lat", "lat0"},
        "from_q": {"lat0", "dq"},
        "merc": {"lat", "lon", "lon0", "r"},
        "merc_inv": {"x", "y", "lon0", "r"},
    }

    def fields(self) -> list[Field]:
        return [
            ChoiceField(
                "mode",
                "mode",
                [
                    ("Tissot: m₁, m₂, T", "tissot"),
                    ("Isometric: φ → q, Δq", "to_q"),
                    ("Isometric: φ₀, Δq → φ", "from_q"),
                    ("Mercator (sphere): φ, λ → x, y", "merc"),
                    ("Mercator (sphere): x, y → φ, λ", "merc_inv"),
                ],
            ),
            NumberField("m1", "m₁ largest scale", "1.2"),
            NumberField("m2", "m₂ smallest scale", "0.9"),
            AngleField("t", "T  direction from m₁", None, "30"),
            AngleField("lat", "φ  latitude", LAT, "41 00 30 N"),
            AngleField("lat0", "φ₀ series origin", LAT, "40 N"),
            NumberField("dq", "Δq isometric difference", "0.0135"),
            AngleField("lon", "λ  longitude", LON, "29 E"),
            AngleField("lon0", "λ₀ central meridian", LON, "0"),
            NumberField("r", "R  sphere radius (m)", "6371000"),
            NumberField("x", "x  north (m)", "5000000"),
            NumberField("y", "y  east (m)", "3224600"),
        ]

    def calculate(self, trace: Trace) -> list[ResultRow]:
        mode, ell = self.choice("mode"), self.geolab.ellipsoid
        if mode == "tissot":
            r = tissot(self.number("m1"), self.number("m2"), self.number("t"), trace)
            return [
                self.value_row("m   scale in direction T", r.scale),
                self.angle_row("T*  image of T", r.image_angle),
                self.angle_row("T − T* angle distortion", r.angle_distortion),
                self.angle_row("ω   largest T − T*", r.max_angle_distortion),
                self.angle_row("2ω  largest change of an angle", 2 * r.max_angle_distortion),
                self.angle_row("T at ω", r.max_direction),
                self.value_row("Φ   area scale m₁m₂", r.area_scale),
                self.value_row("h", r.h),
                self.value_row("k", r.k),
            ]
        if mode == "to_q":
            lat, lat0 = self.number("lat"), self.number("lat0")
            trace.section("Exact")
            q = isometric_latitude(lat, ell, trace)
            q0 = isometric_latitude(lat0, ell)
            dq = isometric_difference(lat, lat0, ell, trace)
            return [
                self.value_row("q   isometric latitude", q),
                self.value_row("Δq  series from φ₀", dq),
                self.value_row("Δq  exact q − q₀", q - q0),
                self.value_row("series − exact", dq - (q - q0)),
            ]
        if mode == "from_q":
            lat0, dq = self.number("lat0"), self.number("dq")
            dlat = latitude_difference(dq, lat0, ell, trace)
            exact = latitude_from_isometric(isometric_latitude(lat0, ell) + dq, ell)
            return [
                self.angle_row("φ   latitude (series)", lat0 + dlat, LAT),
                self.angle_row("Δφ  series", dlat),
                self.angle_row("φ   exact", exact, LAT),
                self.value_row("series − exact (″)", (lat0 + dlat - exact) / ARCSEC),
            ]
        r, lon0 = self.number("r"), self.number("lon0")
        if mode == "merc":
            p = mercator_sphere(self.number("lat"), self.number("lon"), r, lon0, trace)
            return [
                self.length_row("x   north", p.x),
                self.length_row("y   east", p.y),
                self.value_row("m   scale (all directions)", p.scale),
                self.value_row("Φ   area scale", p.scale**2),
            ]
        lat, lon = mercator_sphere_inverse(self.number("x"), self.number("y"), r, lon0, trace)
        return [
            self.angle_row("φ   latitude", lat, LAT),
            self.angle_row("λ   longitude", lon, LON),
            self.value_row("m   scale", 1 / math.cos(lat)),
        ]
