from ..angles import AngleKind
from ..datum import HelmertParams, helmert, helmert_inverse
from ..ellipsoids import GRS80, HAYFORD
from ..point import Geodetic
from ..trace import Trace
from .calculator import Calculator
from .widgets import AngleField, ChoiceField, Field, NumberField, ResultRow

# Each direction: (source ellipsoid, target ellipsoid, apply the parameters backwards?).
# The parameters are always entered for TUTGA-99 → ED-50, as in the course formula.
DIRECTIONS = {
    "to_ed50": (GRS80, HAYFORD, False),
    "to_tutga": (HAYFORD, GRS80, True),
}


class DatumCalculator(Calculator):
    TITLE = "Datum transformation (Helmert, TUTGA-99 ↔ ED-50)"

    def fields(self) -> list[Field]:
        return [
            ChoiceField(
                "direction",
                "direction",
                [("TUTGA-99 → ED-50", "to_ed50"), ("ED-50 → TUTGA-99", "to_tutga")],
            ),
            AngleField("lat", "φ  latitude", AngleKind.LAT, "41 00 30 N"),
            AngleField("lon", "λ  longitude", AngleKind.LON, "29 00 00 E"),
            NumberField("h", "h  height (m)", "150"),
            NumberField("tx", "ΔX (m)", "87"),
            NumberField("ty", "ΔY (m)", "98"),
            NumberField("tz", "ΔZ (m)", "121"),
            NumberField("rx", "R_X (″)", "0"),
            NumberField("ry", "R_Y (″)", "0"),
            NumberField("rz", "R_Z (″)", "0"),
            NumberField("ds", "ΔS (ppm)", "0"),
        ]

    def params(self) -> HelmertParams:
        n = self.number
        return HelmertParams.from_units(
            n("tx"), n("ty"), n("tz"), n("rx"), n("ry"), n("rz"), n("ds")
        )

    def calculate(self, trace: Trace) -> list[ResultRow]:
        source, target, backwards = DIRECTIONS[self.choice("direction")]
        point = Geodetic(self.number("lat"), self.number("lon"), self.number("h"))
        p = self.params()
        trace.section(f"Geodetic → Cartesian ({source.name})")
        xyz = point.to_cartesian(source, trace)
        moved = helmert_inverse(xyz, p, trace) if backwards else helmert(xyz, p, trace)
        trace.section(f"Cartesian → geodetic ({target.name})")
        out = moved.to_geodetic(target, trace=trace)
        return [
            self.angle_row(f"φ   on {target.name}", out.lat, AngleKind.LAT),
            self.angle_row(f"λ   on {target.name}", out.lon, AngleKind.LON),
            self.length_row(f"h   on {target.name}", out.height),
            self.length_row("X   source", xyz.x),
            self.length_row("Y   source", xyz.y),
            self.length_row("Z   source", xyz.z),
            self.length_row("X′  target", moved.x),
            self.length_row("Y′  target", moved.y),
            self.length_row("Z′  target", moved.z),
        ]
