from ..angles import AngleKind, format_angle
from ..map_sheets import SCALES, sheet_at, sheet_corner_coordinates, sheet_from_name
from ..trace import Trace
from ..zones import geodetic_from_grid, grid_from_geodetic, transform_zone
from .calculator import Calculator
from .widgets import AngleField, ChoiceField, Field, NumberField, ResultRow, StringField

LAT, LON = AngleKind.LAT, AngleKind.LON


class ZonesCalculator(Calculator):
    TITLE = "Zones (3°, UTM) and map sheets"

    MODE_FIELDS = {
        "to_grid": {"lat", "lon"},
        "from_grid": {"easting", "northing", "zone"},
        "zone2zone": {"easting", "northing", "zone", "zone_to"},
        "sheet": {"sheet", "zone"},
    }

    def fields(self) -> list[Field]:
        return [
            ChoiceField(
                "mode",
                "mode",
                [
                    ("φ, λ → grid and sheets", "to_grid"),
                    ("Grid → φ, λ", "from_grid"),
                    ("Zone to zone", "zone2zone"),
                    ("Sheet → corners", "sheet"),
                ],
            ),
            ChoiceField("width", "zones", [("3° (m₀ = 1)", 3), ("6° UTM (m₀ = 0.9996)", 6)]),
            AngleField("lat", "φ  latitude", LAT, "40 42 N"),
            AngleField("lon", "λ  longitude", LON, "30 30 E"),
            NumberField("easting", "sağa değer (m)", "542259.752"),
            NumberField("northing", "yukarı değer (m)", "4507456.400"),
            StringField("zone", "zone (blank: auto)", "10"),
            NumberField("zone_to", "to zone", "11"),
            StringField("sheet", "sheet name", "G25-d-01-d"),
        ]

    def _zone(self) -> int | None:
        text = self.field("zone").value.strip()
        if not text:
            return None
        if not text.isdigit():
            raise ValueError("zone must be a whole number")
        return int(text)

    @staticmethod
    def text_row(label: str, text: str) -> ResultRow:
        return ResultRow(label, text, text)

    def calculate(self, trace: Trace) -> list[ResultRow]:
        mode, ell, width = self.choice("mode"), self.geolab.ellipsoid, self.choice("width")
        if mode == "to_grid":
            lat, lon = self.number("lat"), self.number("lon")
            p = grid_from_geodetic(lat, lon, ell, width, trace=trace)
            rows = [
                self.value_row("zone", p.zone),
                self.angle_row("λ₀  central meridian", p.central_meridian, LON),
                self.length_row("sağa değer (easting)", p.easting),
                self.length_row("yukarı değer (northing)", p.northing),
                self.length_row("easting with zone", p.prefixed_easting),
                self.length_row("Y   Gauss–Krüger", p.y),
                self.length_row("X   Gauss–Krüger", p.x),
                self.value_row("m₀m grid scale", p.scale),
                self.angle_row("c   grid convergence", p.convergence),
            ]
            try:
                for scale in SCALES:
                    rows.append(self.text_row(f"sheet 1/{scale}", sheet_at(lat, lon, scale).name))
            except ValueError:
                rows.append(self.text_row("sheet", "outside the Turkish sheet grid"))
            return rows
        if mode == "sheet":
            sheet = sheet_from_name(self.text("sheet"))
            corners = sheet_corner_coordinates(sheet, ell, width, self._zone())
            unit = self.geolab.angle_unit
            rows = [
                self.text_row("sheet", sheet.name),
                self.text_row("scale", f"1/{sheet.scale}"),
                self.value_row("zone", next(iter(corners.values())).zone),
            ]
            for key, (lat, lon) in sheet.corners().items():
                g = corners[key]
                rows.append(self.text_row(f"{key}  φ, λ", f"{format_angle(lat, unit, LAT)}  {format_angle(lon, unit, LON)}"))
                rows.append(self.length_row(f"{key}  sağa", g.easting))
                rows.append(self.length_row(f"{key}  yukarı", g.northing))
            return rows
        easting, northing = self.number("easting"), self.number("northing")
        if mode == "from_grid":
            g = geodetic_from_grid(easting, northing, ell, width, self._zone(), trace)
            return [
                self.angle_row("φ   latitude", g.lat, LAT),
                self.angle_row("λ   longitude", g.lon, LON),
                self.value_row("zone", g.zone),
                self.length_row("Y   Gauss–Krüger", g.y),
                self.length_row("X   Gauss–Krüger", g.x),
                self.angle_row("c   grid convergence", g.convergence),
            ]
        zone_to = self.number("zone_to")
        if zone_to != int(zone_to):
            raise ValueError("zone must be a whole number")
        p = transform_zone(easting, northing, ell, width, self._zone(), int(zone_to), trace)
        return [
            self.value_row("zone", p.zone),
            self.length_row("sağa değer (easting)", p.easting),
            self.length_row("yukarı değer (northing)", p.northing),
            self.length_row("easting with zone", p.prefixed_easting),
            self.angle_row("λ₀  central meridian", p.central_meridian, LON),
        ]
