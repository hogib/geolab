import math

from ..angles import AngleKind
from ..geodesic_curve import azimuths_at_latitude, clairaut_constant, vertex_latitude
from ..normal_section import normal_section_corrections, target_height_correction
from ..trace import Trace
from .calculator import Calculator
from .widgets import AngleField, Field, NumberField, ResultRow

ARCSEC = math.pi / (180 * 3600)


class GeodesicCurveCalculator(Calculator):
    TITLE = "Geodesic curve: Clairaut and normal sections"

    def fields(self) -> list[Field]:
        return [
            AngleField("lat", "φ  latitude", AngleKind.LAT, "45 N"),
            AngleField("az", "α  azimuth", None, "45"),
            AngleField("lat2", "φ₂ azimuths at latitude", AngleKind.LAT, "0"),
            NumberField("s", "S  line length (m)", "100000"),
            NumberField("h2", "h₂ target height (m)", "3000"),
        ]

    def calculate(self, trace: Trace) -> list[ResultRow]:
        lat, az, ell = self.number("lat"), self.number("az"), self.geolab.ellipsoid
        trace.section("Clairaut constant")
        k = clairaut_constant(lat, az, ell, trace)
        trace.section("Vertex (α = 90°)")
        top = vertex_latitude(k, ell, trace)
        rows = [
            self.value_row("k   N cos φ sin α (m)", k),
            self.angle_row("φmax highest latitude", top, AngleKind.LAT),
        ]
        lat2 = self.number("lat2")
        trace.section("Azimuths at φ₂")
        try:
            a1, a2 = azimuths_at_latitude(k, lat2, ell, trace)
        except ValueError:
            rows.append(ResultRow("α at φ₂", "not reached", ""))
        else:
            rows += [self.angle_row("α   at φ₂", a1), self.angle_row("α   at φ₂ (other)", a2)]
        c = normal_section_corrections(lat, az, self.number("s"), ell, trace)
        h_corr = target_height_correction(lat, az, self.number("h2"), ell, trace)
        rows += [
            self.value_row("α − α′  geodesic − section (″)", c.geodesic_minus_normal / ARCSEC),
            self.value_row("α′ − α″ section − reverse (″)", c.normal_minus_reverse / ARCSEC),
            self.value_row("d_max   largest gap (mm)", c.max_separation * 1000),
            self.value_row("s − s′  length difference (m)", c.length_difference),
            self.value_row("α′ − γ  target height (″)", h_corr / ARCSEC),
        ]
        return rows
