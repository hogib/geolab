import math

from ..angles import AngleKind
from ..geodesic import direct, inverse
from ..legendre_series import legendre_direct, legendre_inverse
from ..mid_latitude import mid_latitude_direct, mid_latitude_inverse
from ..trace import Trace
from .calculator import Calculator
from .widgets import AngleField, ChoiceField, Field, NumberField, ResultRow

LAT, LON = AngleKind.LAT, AngleKind.LON
ARCSEC = math.pi / (180 * 3600)


class ShortLinesCalculator(Calculator):
    """Legendre series and Gauss mid-latitude methods, compared with Vincenty."""

    TITLE = "Short lines: Legendre series and Gauss mid-latitude"

    MODE_FIELDS = {
        "leg_dir": {"az", "s"},
        "mid_dir": {"az", "s"},
        "leg_inv": {"lat2", "lon2"},
        "mid_inv": {"lat2", "lon2"},
    }

    def fields(self) -> list[Field]:
        # Defaults: the course's mid-latitude example (Bessel ellipsoid).
        return [
            ChoiceField(
                "mode",
                "mode",
                [
                    ("Legendre direct", "leg_dir"),
                    ("Legendre inverse", "leg_inv"),
                    ("Mid-latitude direct", "mid_dir"),
                    ("Mid-latitude inverse", "mid_inv"),
                ],
            ),
            AngleField("lat1", "φ₁ latitude", LAT, "53 50 02.8809 N"),
            AngleField("lon1", "λ₁ longitude", LON, "10 12 04.1772 E"),
            AngleField("az", "α₁₂ azimuth", None, "25 16 31.98"),
            NumberField("s", "s  distance (m)", "47652.596"),
            AngleField("lat2", "φ₂ latitude", LAT, "54 13 15.2891 N"),
            AngleField("lon2", "λ₂ longitude", LON, "10 30 47.2430 E"),
        ]

    def calculate(self, trace: Trace) -> list[ResultRow]:
        mode, ell = self.choice("mode"), self.geolab.ellipsoid
        lat1, lon1 = self.number("lat1"), self.number("lon1")
        if mode in ("leg_dir", "mid_dir"):
            az, s = self.number("az"), self.number("s")
            solve = legendre_direct if mode == "leg_dir" else mid_latitude_direct
            r = solve(lat1, lon1, az, s, ell, trace)
            v = direct(lat1, lon1, az, s, ell)
            return [
                self.angle_row("φ₂  latitude", r.lat2, LAT),
                self.angle_row("λ₂  longitude", r.lon2, LON),
                self.angle_row("α₂₁ back azimuth", r.back_azimuth),
                self.value_row("φ₂ − Vincenty (″)", (r.lat2 - v.lat2) / ARCSEC),
                self.value_row("λ₂ − Vincenty (″)", math.remainder(r.lon2 - v.lon2, math.tau) / ARCSEC),
                self.value_row("α₂₁ − Vincenty (″)", math.remainder(r.back_azimuth - v.back_azimuth, math.tau) / ARCSEC),
            ]
        lat2, lon2 = self.number("lat2"), self.number("lon2")
        solve = legendre_inverse if mode == "leg_inv" else mid_latitude_inverse
        r = solve(lat1, lon1, lat2, lon2, ell, trace)
        v = inverse(lat1, lon1, lat2, lon2, ell)
        return [
            self.length_row("s   distance", r.distance),
            self.angle_row("α₁₂ azimuth", r.azimuth),
            self.angle_row("α₂₁ back azimuth", r.back_azimuth),
            self.value_row("s − Vincenty (mm)", (r.distance - v.distance) * 1000),
            self.value_row("α₁₂ − Vincenty (″)", math.remainder(r.azimuth - v.azimuth1, math.tau) / ARCSEC),
        ]
