from typing import TYPE_CHECKING, Any, cast

from textual.app import ComposeResult
from textual.containers import Horizontal, Vertical
from textual.widgets import Input

from ..angles import AngleKind, format_angle
from ..trace import Trace
from .widgets import (
    ChoiceField,
    Field,
    FieldError,
    Form,
    ResultPanel,
    ResultRow,
    StringField,
    TextField,
    WorkingPanel,
    WorkingView,
    format_length,
)

if TYPE_CHECKING:
    from .app import GeolabApp


class Calculator(Horizontal):
    """Input form, result and working blocks that update as you type.

    Subclasses set TITLE, return their fields from ``fields`` and implement
    ``calculate``, which returns result rows and fills in ``trace``.
    """

    TITLE = ""

    def compose(self) -> ComposeResult:
        with Vertical(classes="left"):
            yield Form(*self.fields(), classes="block")
            yield ResultPanel(classes="block")
        yield WorkingView(classes="block")

    def fields(self) -> list[Field]:
        raise NotImplementedError

    def calculate(self, trace: Trace) -> list[ResultRow]:
        raise NotImplementedError

    @property
    def geolab(self) -> "GeolabApp":
        return cast("GeolabApp", self.app)

    @property
    def form(self) -> Form:
        return self.query_one(Form)

    def on_mount(self) -> None:
        self.form.border_title = f"Input · {self.TITLE}"
        self.query_one(ResultPanel).border_title = "Result"
        for view in self.query(WorkingView):
            view.border_title = "Working"
        self.watch(self.app, "ellipsoid", lambda: self.ellipsoid_changed(), init=False)
        self.watch(self.app, "angle_unit", lambda: self.recalculate(), init=False)
        self.recalculate()

    def activated(self) -> None:
        """Called when the tab becomes the active one."""

    def ellipsoid_changed(self) -> None:
        self.recalculate()

    def on_input_changed(self, event: Input.Changed) -> None:
        self.recalculate()

    def on_choice_field_changed(self, event: ChoiceField.Changed) -> None:
        self.recalculate()

    def field(self, key: str) -> Field:
        return next(f for f in self.form.fields if f.key == key)

    def number(self, key: str) -> float:
        field = self.field(key)
        assert isinstance(field, TextField)
        return field.parse()

    def text(self, key: str) -> str:
        field = self.field(key)
        assert isinstance(field, StringField)
        return field.parse()

    def choice(self, key: str) -> Any:
        field = self.field(key)
        assert isinstance(field, ChoiceField)
        return field.value

    def angle_row(self, label: str, radians: float, kind: AngleKind | None = None) -> ResultRow:
        text = format_angle(radians, self.geolab.angle_unit, kind)
        return ResultRow(label, text, text)

    def length_row(self, label: str, metres: float) -> ResultRow:
        return ResultRow(label, format_length(metres), f"{metres:.4f}")

    @staticmethod
    def value_row(label: str, value: float | int) -> ResultRow:
        text = f"{value:.15g}" if isinstance(value, float) else str(value)
        return ResultRow(label, text, text)

    def recalculate(self) -> None:
        results = self.query_one(ResultPanel)
        working = next(iter(self.query(WorkingPanel)), None)
        trace = Trace()
        try:
            rows = self.calculate(trace)
        except FieldError as exc:
            results.show_message(str(exc), error=not exc.empty)
            trace = Trace()
        except (ValueError, ArithmeticError) as exc:
            results.show_message(str(exc), error=True)
        else:
            results.show(rows)
        if working is not None:
            working.show(trace)
