from textual.binding import Binding

from ..ellipsoids import Ellipsoid
from ..trace import Trace
from .calculator import Calculator
from .widgets import Field, FieldError, NumberField, ResultRow, StringField, TextField


class EllipsoidEditor(Calculator):
    """Shows the current ellipsoid's parameters; edit them to define a custom one."""

    TITLE = "Ellipsoid"
    BINDINGS = [Binding("s", "save", "Save and use", show=False)]

    def fields(self) -> list[Field]:
        ell = self.geolab.ellipsoid
        return [
            StringField("name", "name", ell.name),
            NumberField("a", "a  (m)", repr(ell.a)),
            NumberField("inv_f", "1/f", repr(ell.inv_f)),
        ]

    def on_mount(self) -> None:
        super().on_mount()
        self.form.border_subtitle = "s save & use"

    def ellipsoid_changed(self) -> None:
        ell = self.geolab.ellipsoid
        for key, value in (("name", ell.name), ("a", repr(ell.a)), ("inv_f", repr(ell.inv_f))):
            field = self.field(key)
            assert isinstance(field, TextField)
            field.value = value
        self.recalculate()

    def build(self) -> Ellipsoid:
        return Ellipsoid(self.text("name"), self.number("a"), self.number("inv_f"))

    def calculate(self, trace: Trace) -> list[ResultRow]:
        ell = self.build()
        trace.step("f", ell.f, "1 / (1/f)")
        trace.step("b", ell.b, "a (1 − f)", unit="m")
        trace.step("e²", ell.e_sq, "2f − f²")
        trace.step("e′²", ell.ep_sq, "e² / (1 − e²)")
        trace.step("n", ell.n, "f / (2 − f)")
        trace.step("R₁", ell.mean_radius, "(2a + b) / 3", unit="m")
        trace.step("R₂", ell.authalic_radius, "√(a²/2 (1 + (1 − e²)/e · atanh e))", unit="m")
        trace.step("R₃", ell.volumetric_radius, "∛(a² b)", unit="m")
        return [
            self.length_row("b   semi-minor axis", ell.b),
            self.value_row("f   flattening", ell.f),
            self.value_row("e²  1st eccentricity²", ell.e_sq),
            self.value_row("e′² 2nd eccentricity²", ell.ep_sq),
            self.value_row("n   3rd flattening", ell.n),
            self.length_row("R₁  mean radius", ell.mean_radius),
            self.length_row("R₂  authalic radius", ell.authalic_radius),
            self.length_row("R₃  volumetric radius", ell.volumetric_radius),
        ]

    def action_save(self) -> None:
        try:
            ell = self.build()
        except (FieldError, ValueError) as exc:
            self.geolab.flash(str(exc), error=True)
            return
        self.geolab.add_ellipsoid(ell)
