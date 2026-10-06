from ..angles import AngleKind
from ..point import Cartesian3D, Geodetic, GeodeticMethod, Latitude, LatType
from ..trace import Trace
from .calculator import Calculator
from .widgets import AngleField, ChoiceField, Field, NumberField, ResultRow

LAT_SYMBOLS = {LatType.GEODETIC: "φ", LatType.GEOCENTRIC: "ψ", LatType.PARAMETRIC: "β"}


class GeodeticToEcef(Calculator):
    TITLE = "Geodetic → ECEF"

    def fields(self) -> list[Field]:
        return [
            AngleField("lat", "φ  latitude", AngleKind.LAT, "41 00 30 N"),
            AngleField("lon", "λ  longitude", AngleKind.LON, "29 00 00 E"),
            NumberField("h", "h  height (m)", "150"),
        ]

    def calculate(self, trace: Trace) -> list[ResultRow]:
        point = Geodetic(self.number("lat"), self.number("lon"), self.number("h"))
        c = point.to_cartesian(self.geolab.ellipsoid, trace)
        return [self.length_row("X", c.x), self.length_row("Y", c.y), self.length_row("Z", c.z)]


class EcefToGeodetic(Calculator):
    TITLE = "ECEF → Geodetic"

    def fields(self) -> list[Field]:
        return [
            NumberField("x", "X  (m)", "4215751.8328"),
            NumberField("y", "Y  (m)", "2336829.3996"),
            NumberField("z", "Z  (m)", "4163220.0278"),
            ChoiceField(
                "method",
                "method",
                [("Iterative", GeodeticMethod.ITERATIVE), ("Bowring", GeodeticMethod.BOWRING)],
            ),
        ]

    def calculate(self, trace: Trace) -> list[ResultRow]:
        point = Cartesian3D(self.number("x"), self.number("y"), self.number("z"))
        g = point.to_geodetic(self.geolab.ellipsoid, self.choice("method"), trace)
        return [
            self.angle_row("φ", g.lat, AngleKind.LAT),
            self.angle_row("λ", g.lon, AngleKind.LON),
            self.length_row("h", g.height),
        ]


class LatitudeConverter(Calculator):
    TITLE = "Latitude types"

    def fields(self) -> list[Field]:
        return [
            AngleField("value", "latitude", AngleKind.LAT, "41 00 30 N"),
            ChoiceField("type", "given as", [(t.name.capitalize(), t) for t in LatType]),
        ]

    def calculate(self, trace: Trace) -> list[ResultRow]:
        source = Latitude(self.number("value"), self.choice("type"))
        rows = []
        for target in LatType:
            if target != source.lat_type:
                trace.section(f"→ {target.name.capitalize()}")
            result = source.convert_to(target, self.geolab.ellipsoid, trace)
            label = f"{LAT_SYMBOLS[target]}  {target.name.lower()}"
            rows.append(self.angle_row(label, result.value, AngleKind.LAT))
        return rows
