import math

from ..gk_reductions import direct_by_reductions, reduce_area, reduce_line, traverse
from ..trace import Trace
from .calculator import Calculator
from .widgets import (
    AngleField,
    AngleListField,
    ChoiceField,
    Field,
    NumberField,
    NumberListField,
    PointField,
    PointListField,
    ResultRow,
)

ARCSEC = math.pi / (180 * 3600)


class ReductionsCalculator(Calculator):
    TITLE = "Gauss–Krüger reductions and traverse"

    MODE_FIELDS = {
        "line": {"p1", "p2"},
        "direct": {"p1", "t", "s"},
        "area": {"polygon"},
        "traverse": {"a", "b", "c", "d", "angles", "sides"},
    }

    def fields(self) -> list[Field]:
        # Traverse defaults: an exact WGS84 traverse 0.9–1.2° east of λ₀ = 30°.
        return [
            ChoiceField(
                "mode",
                "mode",
                [
                    ("Second problem: Pᵢ, Pₖ → S, T", "line"),
                    ("First problem: Pᵢ, T, S → Pₖ", "direct"),
                    ("Area of a polygon", "area"),
                    ("Traverse A-B … C-D", "traverse"),
                ],
            ),
            PointField("p1", "Pᵢ  X Y (m)", "4500000 20000"),
            PointField("p2", "Pₖ  X Y (m)", "4520000 45000"),
            AngleField("t", "Tᵢₖ geodesic grid bearing", None, "51 20 26.128"),
            NumberField("s", "S   geodesic length (m)", "32015.1846"),
            PointListField(
                "polygon",
                "vertices X Y; X Y; …",
                "4500000 40000; 4500000 60000; 4503000 60000; 4503000 40000",
            ),
            PointField("a", "A  X Y (orientation)", "4485437.913 76291.485"),
            PointField("b", "B  X Y (start)", "4487703.294 80506.038"),
            PointField("c", "C  X Y (end)", "4492325.145 95703.152"),
            PointField("d", "D  X Y (orientation)", "4495689.644 98200.120"),
            AngleListField(
                "angles",
                "β₀; β₁; …  (at B, new points, C)",
                "192 56 19.27; 168 28 56.24; 198 40 38.58; 169 17 30.31; 145 27 01.72",
            ),
            NumberListField("sides", "S₁; S₂; …  (m)", "4379.919; 3776.054; 4271.635; 3565.283"),
        ]

    def calculate(self, trace: Trace) -> list[ResultRow]:
        mode, ell = self.choice("mode"), self.geolab.ellipsoid
        if mode == "line":
            r = reduce_line(*self.number("p1"), *self.number("p2"), ell, trace)
            return [
                self.length_row("S   geodesic length", r.distance),
                self.angle_row("Tᵢₖ grid bearing", r.bearing1),
                self.angle_row("Tₖᵢ back grid bearing", r.bearing2),
                self.length_row("s   plane distance", r.plane_distance),
                self.angle_row("tᵢₖ plane bearing", r.plane_bearing),
                self.value_row("δS  = S − s (mm)", r.distance_reduction * 1000),
                self.value_row("δTᵢₖ = Tᵢₖ − tᵢₖ (″)", r.bearing_reduction1 / ARCSEC),
                self.value_row("δTₖᵢ = Tₖᵢ − tₖᵢ (″)", r.bearing_reduction2 / ARCSEC),
            ]
        if mode == "direct":
            r = direct_by_reductions(*self.number("p1"), self.number("t"), self.number("s"), ell, trace)
            return [
                self.length_row("Xₖ  north", r.x2),
                self.length_row("Yₖ  east", r.y2),
                self.angle_row("Tₖᵢ back grid bearing", r.bearing2),
                self.length_row("s   plane distance", r.plane_distance),
                self.angle_row("tᵢₖ plane bearing", r.plane_bearing),
                self.value_row("iterations", r.iterations),
            ]
        if mode == "area":
            r = reduce_area(self.number("polygon"), ell, trace)
            return [
                self.value_row("F   ellipsoid area (m²)", r.area),
                self.value_row("f   plane area (m²)", r.plane_area),
                self.value_row("δF  = F − f (m²)", r.reduction),
            ]
        angles, sides = self.number("angles"), self.number("sides")
        t = traverse(
            self.number("a"), self.number("b"), self.number("c"), self.number("d"), angles, sides, ell, trace
        )
        rows = [
            self.value_row("w_β angle misclosure (″)", t.angle_misclosure / ARCSEC),
            self.value_row("w_X misclosure (mm)", t.misclosure_x * 1000),
            self.value_row("w_Y misclosure (mm)", t.misclosure_y * 1000),
        ]
        for k, (x, y) in enumerate(t.points, 1):
            rows.append(self.length_row(f"X{k}", x))
            rows.append(self.length_row(f"Y{k}", y))
        return rows
