"""Keyboard-driven building blocks: form fields, the form, and output panels."""

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any

from rich.console import Group
from rich.table import Table
from rich.text import Text
from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal, Vertical, VerticalScroll
from textual.message import Message
from textual.reactive import reactive
from textual.widgets import Input, Label, Static

from ..angles import AngleKind, format_dms, parse_angle
from ..trace import Section, Step, Trace


class FieldError(Exception):
    """A form field is empty or invalid. ``empty`` fields are a hint, not an error."""

    def __init__(self, message: str, empty: bool = False):
        super().__init__(message)
        self.empty = empty


def format_length(metres: float) -> str:
    return f"{metres:,.4f} m".replace(",", " ")


# ── Fields ────────────────────────────────────────────────────────────────────


class Field(Horizontal):
    """One row of a Form: a selection marker, a label and a value."""

    def __init__(self, key: str, label: str):
        super().__init__()
        self.key = key
        self.label = label

    def compose(self) -> ComposeResult:
        yield Static("  ", classes="marker")
        yield Label(self.label, classes="field-label")
        yield from self.compose_value()

    def compose_value(self) -> ComposeResult:
        raise NotImplementedError

    def set_selected(self, selected: bool) -> None:
        self.set_class(selected, "-selected")
        self.query_one(".marker", Static).update("▶ " if selected else "  ")


class TextField(Field):
    """A field edited in INSERT mode. Subclasses implement ``parse_value``."""

    def __init__(self, key: str, label: str, value: str = ""):
        super().__init__(key, label)
        self._initial = value

    def compose_value(self) -> ComposeResult:
        field_input = Input(self._initial, compact=True, select_on_focus=False)
        field_input.can_focus = False  # only focusable while in INSERT mode
        yield field_input

    @property
    def input(self) -> Input:
        return self.query_one(Input)

    @property
    def value(self) -> str:
        return self.input.value

    @value.setter
    def value(self, text: str) -> None:
        self.input.value = text

    def parse_value(self, text: str) -> Any:
        raise NotImplementedError

    def parse(self) -> Any:
        """The field's value; raises FieldError if it is empty or invalid."""
        name = " ".join(self.label.split())
        if not self.value.strip():
            raise FieldError(f"Enter {name}", empty=True)
        try:
            return self.parse_value(self.value)
        except ValueError as exc:
            raise FieldError(f"{name}: {exc}") from exc

    def on_input_changed(self, event: Input.Changed) -> None:
        try:
            self.parse()
            invalid = False
        except FieldError as exc:
            invalid = not exc.empty
        self.set_class(invalid, "-invalid")


class AngleField(TextField):
    def __init__(self, key: str, label: str, kind: AngleKind | None = None, value: str = ""):
        super().__init__(key, label, value)
        self.kind = kind

    def parse_value(self, text: str) -> float:
        return parse_angle(text, self.kind)


class NumberField(TextField):
    def parse_value(self, text: str) -> float:
        try:
            return float(text)
        except ValueError:
            raise ValueError("not a number") from None


class StringField(TextField):
    def parse_value(self, text: str) -> str:
        return text.strip()


class ChoiceField(Field):
    """A field with a fixed set of options, cycled with h/l."""

    class Changed(Message):
        def __init__(self, field: "ChoiceField"):
            super().__init__()
            self.field = field

    def __init__(self, key: str, label: str, options: Sequence[tuple[str, Any]], index: int = 0):
        super().__init__(key, label)
        self.options = list(options)
        self.index = index

    def compose_value(self) -> ComposeResult:
        yield Static(self._render_options(), classes="choice")

    @property
    def value(self) -> Any:
        return self.options[self.index][1]

    @property
    def value_label(self) -> str:
        return self.options[self.index][0]

    def cycle(self, step: int) -> None:
        self.select((self.index + step) % len(self.options))

    def select(self, index: int) -> None:
        self.index = index
        self.query_one(".choice", Static).update(self._render_options())
        self.post_message(self.Changed(self))

    def select_label(self, prefix: str) -> bool:
        """Select the option whose label starts with ``prefix`` (case-insensitive)."""
        matches = [i for i, (label, _) in enumerate(self.options) if label.lower().startswith(prefix.lower())]
        if len(matches) != 1:
            return False
        self.select(matches[0])
        return True

    def _render_options(self) -> Text:
        text = Text()
        for i, (label, _) in enumerate(self.options):
            if i:
                text.append(" │ ", style="dim")
            text.append(label, style="bold reverse" if i == self.index else "dim")
        return text


# ── Form ──────────────────────────────────────────────────────────────────────


class Form(Vertical, can_focus=True):
    """Vertical list of fields navigated in NORMAL mode."""

    BINDINGS = [
        Binding("j,down", "move(1)", "Next field", show=False),
        Binding("k,up", "move(-1)", "Previous field", show=False),
        Binding("i,a,enter", "edit", "Edit", show=False),
        Binding("c", "change", "Clear and edit", show=False),
        Binding("l,right,space", "cycle(1)", "Next option", show=False),
        Binding("h,left", "cycle(-1)", "Previous option", show=False),
    ]

    selected: reactive[int] = reactive(0, init=False)

    class EditRequested(Message):
        def __init__(self, field: TextField):
            super().__init__()
            self.field = field

    @property
    def fields(self) -> list[Field]:
        return list(self.query(Field))

    @property
    def current(self) -> Field:
        return self.fields[self.selected]

    def on_mount(self) -> None:
        self.watch_selected(self.selected)

    def watch_selected(self, selected: int) -> None:
        for i, field in enumerate(self.fields):
            field.set_selected(i == selected)

    def action_move(self, step: int) -> None:
        self.selected = max(0, min(len(self.fields) - 1, self.selected + step))

    def select_next(self) -> None:
        self.action_move(1)

    def action_edit(self) -> None:
        field = self.current
        if isinstance(field, ChoiceField):
            field.cycle(1)
        elif isinstance(field, TextField):
            self.post_message(self.EditRequested(field))

    def action_change(self) -> None:
        field = self.current
        if isinstance(field, TextField):
            field.value = ""
            self.post_message(self.EditRequested(field))

    def action_cycle(self, step: int) -> None:
        if isinstance(self.current, ChoiceField):
            self.current.cycle(step)


# ── Output panels ─────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class ResultRow:
    label: str
    display: str
    copy: str
    """Text placed on the clipboard when the row is yanked."""


class ResultPanel(Static, can_focus=True):
    """Result rows; j/k select a row to yank, or a message when there is none."""

    BINDINGS = [
        Binding("j,down", "move(1)", "Next row", show=False),
        Binding("k,up", "move(-1)", "Previous row", show=False),
    ]

    selected: reactive[int] = reactive(0, init=False)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.rows: list[ResultRow] = []
        self.message = ""

    def show(self, rows: Sequence[ResultRow]) -> None:
        self.rows = list(rows)
        self.message = ""
        self.selected = min(self.selected, max(len(self.rows) - 1, 0))
        self.remove_class("-error", "-hint")
        self._render_rows()

    def show_message(self, message: str, error: bool) -> None:
        self.rows = []
        self.message = message
        self.set_class(error, "-error")
        self.set_class(not error, "-hint")
        self.update(Text(message))

    @property
    def selected_row(self) -> ResultRow | None:
        return self.rows[self.selected] if self.rows else None

    def action_move(self, step: int) -> None:
        if self.rows:
            self.selected = max(0, min(len(self.rows) - 1, self.selected + step))

    def watch_selected(self) -> None:
        self._render_rows()

    def on_focus(self) -> None:
        self._render_rows()

    def on_blur(self) -> None:
        self._render_rows()

    def _render_rows(self) -> None:
        if not self.rows:
            return
        accent = self.app.theme_variables["accent"]
        table = Table.grid(padding=(0, 2))
        table.add_column(width=2)
        table.add_column(style=f"bold {accent}", no_wrap=True)
        table.add_column(style="bold")
        for i, row in enumerate(self.rows):
            marker = "▶" if self.has_focus and i == self.selected else ""
            table.add_row(marker, row.label, row.display)
        self.update(table)


class WorkingPanel(Static):
    """Step-by-step intermediate values of the last calculation."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.trace = Trace()

    def show(self, trace: Trace) -> None:
        self.trace = trace
        if not trace.entries:
            self.update(Text("No steps for this calculation.", style="dim"))
            return
        heading_style = f"bold {self.app.theme_variables['primary']}"
        parts: list[Text | Table] = []
        table = None
        for entry in trace.entries:
            if isinstance(entry, Section):
                parts.append(Text(f"── {entry.title}", style=heading_style))
                table = None
                continue
            if table is None:
                table = Table.grid(padding=(0, 2))
                table.add_column(style="bold", width=8, no_wrap=True)
                table.add_column(justify="right", width=22, no_wrap=True)
                table.add_column(style="dim", ratio=1)
                parts.append(table)
            table.add_row(entry.name, format_step_value(entry), entry.formula)
        self.update(Group(*parts))


class WorkingView(VerticalScroll, can_focus=True):
    """Scrollable container for the working panel, with vim scrolling keys."""

    BINDINGS = [
        Binding("j,down", "scroll_down", "Scroll down", show=False),
        Binding("k,up", "scroll_up", "Scroll up", show=False),
        Binding("ctrl+d", "page_down", "Page down", show=False),
        Binding("ctrl+u", "page_up", "Page up", show=False),
        Binding("g", "scroll_home", "Top", show=False),
        Binding("G", "scroll_end", "Bottom", show=False),
    ]

    def compose(self) -> ComposeResult:
        yield WorkingPanel()


def format_step_value(step: Step) -> str:
    if step.unit == "rad":
        return f"{format_dms(step.value, sec_decimals=5)}\n{step.value:.12f} rad"
    if step.unit == "m":
        return format_length(step.value)
    return f"{step.value:.12g}"



# ── Command line ──────────────────────────────────────────────────────────────


class CommandLine(Input):
    """The ``:`` command input. Tab / shift+tab cycle through completions."""

    BINDINGS = [
        Binding("tab", "complete(1)", "Complete", show=False),
        Binding("shift+tab", "complete(-1)", "Complete backwards", show=False),
    ]

    class CompletionChanged(Message):
        pass

    def __init__(self, completer: Callable[[str], list[str]], **kwargs):
        super().__init__(compact=True, **kwargs)
        self.completer = completer
        self.candidates: list[str] = []
        self.candidate_index = -1
        self._completed_value: str | None = None

    def action_complete(self, step: int) -> None:
        stale = self.value != self._completed_value
        # After a unique completion (e.g. "ell "), the next Tab completes its argument.
        if stale or len(self.candidates) == 1:
            self.candidates = self.completer(self.value)
            self.candidate_index = -1 if step > 0 else 0
        if not self.candidates:
            self.app.bell()
            self.post_message(self.CompletionChanged())
            return
        self.candidate_index = (self.candidate_index + step) % len(self.candidates)
        self.value = self.candidates[self.candidate_index]
        self.cursor_position = len(self.value)
        self._completed_value = self.value
        self.post_message(self.CompletionChanged())

    def reset_completion(self) -> None:
        self.candidates = []
        self.candidate_index = -1
        self._completed_value = None
        self.post_message(self.CompletionChanged())

    def on_input_changed(self, event: Input.Changed) -> None:
        if self.candidates and self.value != self._completed_value:
            self.reset_completion()
