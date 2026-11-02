from ..angles import AngleKind
from ..area import ellipsoid_area, quadrangle_area, zone_area
from ..trace import Trace
from .calculator import Calculator
from .widgets import AngleField, Field, ResultRow


class AreaCalculator(Calculator):
    TITLE = "Area between meridians and parallels"

    def fields(self) -> list[Field]:
        # Default: a 1:25 000 map sheet, 7′30″ × 7′30″.
        return [
            AngleField("lat1", "φ₁ latitude", AngleKind.LAT, "39 00 00 N"),
            AngleField("lat2", "φ₂ latitude", AngleKind.LAT, "39 07 30 N"),
            AngleField("lon1", "λ₁ longitude", AngleKind.LON, "29 00 00 E"),
            AngleField("lon2", "λ₂ longitude", AngleKind.LON, "29 07 30 E"),
        ]

    def calculate(self, trace: Trace) -> list[ResultRow]:
        lat1, lat2 = self.number("lat1"), self.number("lat2")
        lon1, lon2 = self.number("lon1"), self.number("lon2")
        ell = self.geolab.ellipsoid
        p = quadrangle_area(lat1, lat2, lon1, lon2, ell, trace)
        z = zone_area(lat1, lat2, ell)
        trace.section("Whole ellipsoid")
        total = ellipsoid_area(ell, trace)
        return [
            self.value_row("P   quadrangle (m²)", p),
            self.value_row("P   quadrangle (km²)", p / 1e6),
            self.value_row("Z   zone φ₁…φ₂ (km²)", z / 1e6),
            self.value_row("E   whole ellipsoid (km²)", total / 1e6),
        ]
