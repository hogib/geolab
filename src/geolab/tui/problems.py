from ..angles import AngleKind
from ..curvature import (
    azimuth_radius,
    gaussian_mean_radius,
    meridian_radius,
    parallel_radius,
    prime_vertical_radius,
)
from ..geodesic import direct, inverse
from ..trace import Trace
from .calculator import Calculator
from .widgets import AngleField, Field, NumberField, ResultRow

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
    TITLE = "Radii of curvature"

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
        return [
            self.length_row("M   meridian", m),
            self.length_row("N   prime vertical", n),
            self.length_row("r   parallel", r),
            self.length_row("R   Gaussian mean", gauss),
            self.length_row("Rα  in azimuth α", r_az),
        ]
