import math

from ..angles import AngleKind
from ..curvature import gaussian_mean_radius
from ..trace import Trace
from ..triangles import (
    Triangle,
    additament_asa,
    legendre_asa,
    plane_asa,
    plane_sas,
    plane_sss,
    spherical_asa,
    spherical_sas,
    spherical_sss,
)
from .calculator import Calculator
from .widgets import AngleField, ChoiceField, Field, NumberField, ResultRow

ARCSEC = math.pi / (180 * 3600)


class TriangleCalculator(Calculator):
    TITLE = "Triangles: plane and spherical (Legendre, additament)"

    MODE_FIELDS = {
        "asa": {"c", "alpha", "beta"},
        "sas": {"a", "b", "gamma"},
        "sss": {"a", "b", "c"},
    }

    def fields(self) -> list[Field]:
        # Defaults: the course example (c, α, β near 53°50′ N).
        return [
            ChoiceField(
                "mode",
                "given",
                [("AKA: c, α, β", "asa"), ("KAK: a, b, γ", "sas"), ("KKK: a, b, c", "sss")],
            ),
            ChoiceField("surface", "surface", [("Gauss sphere √(MN)", "sphere"), ("Plane", "plane")]),
            AngleField("lat", "φ  for R = √(MN)", AngleKind.LAT, "53 50 02.8809 N"),
            NumberField("a", "a  side P₂P₃ (m)", "36594.437"),
            NumberField("b", "b  side P₁P₃ (m)", "49609.192"),
            NumberField("c", "c  side P₁P₂ (m)", "47652.597"),
            AngleField("alpha", "α  angle at P₁", None, "44 08 45.680"),
            AngleField("beta", "β  angle at P₂", None, "70 46 03.490"),
            AngleField("gamma", "γ  angle at P₃", None, "65 05 15.00"),
        ]

    def calculate(self, trace: Trace) -> list[ResultRow]:
        mode, sphere = self.choice("mode"), self.choice("surface") == "sphere"
        radius = gaussian_mean_radius(self.number("lat"), self.geolab.ellipsoid, trace) if sphere else 0.0
        if mode == "asa":
            args = (self.number("c"), self.number("alpha"), self.number("beta"))
            t = spherical_asa(*args, radius, trace) if sphere else plane_asa(*args, trace)
        elif mode == "sas":
            args = (self.number("a"), self.number("b"), self.number("gamma"))
            t = spherical_sas(*args, radius, trace) if sphere else plane_sas(*args, trace)
        else:
            args = (self.number("a"), self.number("b"), self.number("c"))
            t = spherical_sss(*args, radius, trace) if sphere else plane_sss(*args, trace)
        rows = self._rows(t)
        if sphere:
            rows.insert(6, self.value_row("ε   spherical excess (″)", t.excess / ARCSEC))
            rows.append(self.length_row("R   Gauss radius", radius))
        if sphere and mode == "asa":
            leg = legendre_asa(*args, radius, trace)
            add = additament_asa(*args, radius, trace)
            rows += [
                self.length_row("a   Legendre (ε/3)", leg.a),
                self.length_row("b   Legendre (ε/3)", leg.b),
                self.length_row("a   additament", add.a),
                self.length_row("b   additament", add.b),
                self.value_row("Legendre − exact a (mm)", (leg.a - t.a) * 1000),
                self.value_row("additament − exact a (mm)", (add.a - t.a) * 1000),
            ]
        return rows

    def _rows(self, t: Triangle) -> list[ResultRow]:
        return [
            self.length_row("a   P₂P₃", t.a),
            self.length_row("b   P₁P₃", t.b),
            self.length_row("c   P₁P₂", t.c),
            self.angle_row("α   at P₁", t.alpha),
            self.angle_row("β   at P₂", t.beta),
            self.angle_row("γ   at P₃", t.gamma),
            self.value_row("F   area (m²)", t.area),
        ]
