from enum import Enum, auto

from rich.text import Text
from textual import events
from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal
from textual.reactive import reactive
from textual.timer import Timer
from textual.widgets import ContentSwitcher, Input, Static

from ..angles import AngleUnit
from ..ellipsoids import ELLIPSOIDS, WGS84, Ellipsoid
from .calculator import Calculator
from .commands import CommandError, completions, match_ellipsoid, match_unit, split_command
from .convert import EcefToGeodetic, GeodeticToEcef, LatitudeConverter
from .screens import EllipsoidPicker, HelpScreen
from .widgets import ChoiceField, CommandLine, Form, ResultPanel, WorkingView

TABS: list[tuple[str, type[Calculator]]] = [
    ("Geo→ECEF", GeodeticToEcef),
    ("ECEF→Geo", EcefToGeodetic),
    ("Latitudes", LatitudeConverter),
]
UNIT_NAMES = {AngleUnit.DMS: "DMS", AngleUnit.DEGREES: "DEG", AngleUnit.RADIANS: "RAD"}


class Mode(Enum):
    NORMAL = auto()
    INSERT = auto()
    COMMAND = auto()


MODE_HINTS = {
    Mode.NORMAL: "i edit  y yank  : cmd  ? help",
    Mode.INSERT: "esc normal  enter confirm",
    Mode.COMMAND: "enter run  esc cancel",
}


class GeolabApp(App):
    CSS_PATH = "app.tcss"
    TITLE = "geolab"
    ENABLE_COMMAND_PALETTE = False
    AUTO_FOCUS = None

    BINDINGS = [
        Binding("escape", "normal_mode", show=False),
        *(Binding(str(i), f"tab({i - 1})", show=False) for i in range(1, 10)),
        Binding("right_square_bracket", "next_tab", show=False),
        Binding("left_square_bracket", "prev_tab", show=False),
        Binding("w", "toggle_working", show=False),
        Binding("u", "cycle_unit", show=False),
        Binding("e", "pick_ellipsoid", show=False),
        Binding("y", "yank(False)", show=False),
        Binding("Y", "yank(True)", show=False),
        Binding("colon", "command", show=False),
        Binding("question_mark", "help", show=False),
        Binding("q", "quit", show=False),
    ]

    ellipsoid: reactive[Ellipsoid] = reactive(WGS84, init=False)
    angle_unit: reactive[AngleUnit] = reactive(AngleUnit.DMS, init=False)
    show_working: reactive[bool] = reactive(False, init=False)
    mode: reactive[Mode] = reactive(Mode.NORMAL, init=False)
    tab: reactive[int] = reactive(0, init=False)

    def __init__(self) -> None:
        # Pass the terminal's default colours through, so its background shows (transparency).
        super().__init__(ansi_color=True)
        self.last_yank = ""
        self.last_message = ""
        self._pending_g = False
        self._return_focus = None
        self._message_timer: Timer | None = None

    def compose(self) -> ComposeResult:
        yield Static(id="tabbar")
        with ContentSwitcher(initial="tab-0", id="tabs"):
            for i, (_, calculator) in enumerate(TABS):
                yield calculator(id=f"tab-{i}")
        yield Static(id="statusline")
        with Horizontal(id="echo"):
            yield Static(":", id="cmd-prompt")
            yield CommandLine(self._complete_command, id="cmdline")
            yield Static(id="message")

    def on_mount(self) -> None:
        self.query_one("#cmd-prompt").display = False
        self.query_one("#cmdline").display = False
        self._render_tabbar()
        self._render_status()
        self.active_calculator.form.focus()

    # ── State ────────────────────────────────────────────────────────────────

    @property
    def active_calculator(self) -> Calculator:
        return self.query_one(f"#tab-{self.tab}", Calculator)

    def watch_tab(self, tab: int) -> None:
        self.query_one("#tabs", ContentSwitcher).current = f"tab-{tab}"
        self._render_tabbar()
        self._render_status()
        self.active_calculator.form.focus()

    def watch_mode(self, mode: Mode) -> None:
        commanding = mode == Mode.COMMAND
        self.query_one("#cmd-prompt").display = commanding
        self.query_one("#cmdline").display = commanding
        self.query_one("#message").display = not commanding
        self._render_status()

    def watch_ellipsoid(self) -> None:
        self._render_status()

    def watch_angle_unit(self) -> None:
        self._render_status()

    def watch_show_working(self, show: bool) -> None:
        self.query_one("#tabs").set_class(show, "-show-working")
        if not show and isinstance(self.focused, WorkingView):
            self.active_calculator.form.focus()
        self._render_status()

    # ── Rendering ────────────────────────────────────────────────────────────

    def _render_tabbar(self) -> None:
        theme = self.theme_variables
        text = Text.assemble((" geolab ", f"bold {theme['background']} on {theme['accent']}"), " ")
        for i, (name, _) in enumerate(TABS):
            style = f"bold {theme['background']} on {theme['primary']}" if i == self.tab else "dim"
            text.append(f" {i + 1} {name} ", style=style)
            text.append(" ")
        self.query_one("#tabbar", Static).update(text)

    def _render_status(self) -> None:
        theme = self.theme_variables
        mode_colour = {Mode.NORMAL: "primary", Mode.INSERT: "success", Mode.COMMAND: "warning"}
        ell = self.ellipsoid
        left = Text.assemble(
            (f" {self.mode.name} ", f"bold {theme['background']} on {theme[mode_colour[self.mode]]}"),
            "  ",
            (ell.name, "bold"),
            (f"  a={ell.a:,.3f}  1/f={ell.inv_f}".replace(",", " "), "dim"),
            ("  │  ", "dim"),
            (UNIT_NAMES[self.angle_unit], f"bold {theme['accent']}"),
            ("  │  ", "dim"),
            ("working", "bold") if self.show_working else ("working off", "dim"),
        )
        cmdline = self.query_one(CommandLine)
        if self.mode == Mode.COMMAND and cmdline.candidates:
            left = Text.assemble(left[: len(self.mode.name) + 2], "  ", self._wildmenu(cmdline))
        hints = Text(MODE_HINTS[self.mode] + " ", style="dim")
        status = self.query_one("#statusline", Static)
        gap = self.size.width - left.cell_len - hints.cell_len
        if gap > 0:
            left.append(" " * gap)
            left.append_text(hints)
        status.update(left)

    def _wildmenu(self, cmdline: CommandLine) -> Text:
        text = Text()
        for i, candidate in enumerate(cmdline.candidates):
            label = candidate.strip().partition(" ")[2] or candidate.strip()
            style = f"bold {self.theme_variables['background']} on {self.theme_variables['accent']}"
            text.append(f" {label} ", style=style if i == cmdline.candidate_index else "")
        return text

    def on_command_line_completion_changed(self) -> None:
        self._render_status()

    def _complete_command(self, text: str) -> list[str]:
        methods = [
            label
            for field in self.active_calculator.form.fields
            if field.key == "method" and isinstance(field, ChoiceField)
            for label, _ in field.options
        ]
        return completions(text, methods, len(TABS))

    def on_resize(self) -> None:
        self._render_status()

    def flash(self, message: str, error: bool = False) -> None:
        """Show a message in the echo line for a few seconds."""
        self.last_message = message
        widget = self.query_one("#message", Static)
        widget.update(Text(message, style=self.theme_variables["error"] if error else ""))
        if self._message_timer is not None:
            self._message_timer.stop()
        self._message_timer = self.set_timer(4, lambda: widget.update(""))

    # ── Modes ────────────────────────────────────────────────────────────────

    def on_form_edit_requested(self, event: Form.EditRequested) -> None:
        field_input = event.field.input
        field_input.can_focus = True
        self.mode = Mode.INSERT
        field_input.focus()
        field_input.cursor_position = len(field_input.value)

    def action_normal_mode(self) -> None:
        if self.mode == Mode.INSERT:
            self._leave_insert()
        elif self.mode == Mode.COMMAND:
            self._leave_command()

    def _leave_insert(self) -> None:
        self.mode = Mode.NORMAL
        self.active_calculator.form.focus()

    def on_input_blurred(self, event: Input.Blurred) -> None:
        if event.input.id == "cmdline":
            return
        event.input.can_focus = False
        if self.mode == Mode.INSERT:
            self.mode = Mode.NORMAL

    def on_input_submitted(self, event: Input.Submitted) -> None:
        if event.input.id == "cmdline":
            command = event.value
            self._leave_command()
            self.run_command(command)
            return
        self._leave_insert()
        self.active_calculator.form.select_next()

    def on_key(self, event: events.Key) -> None:
        # Two-key sequences: gt / gT switch tabs.
        if self.mode != Mode.NORMAL:
            self._pending_g = False
            return
        if self._pending_g:
            self._pending_g = False
            if event.character in ("t", "T"):
                event.stop()
                self.action_next_tab() if event.character == "t" else self.action_prev_tab()
            return
        if event.character == "g" and not isinstance(self.focused, WorkingView):
            self._pending_g = True
            event.stop()

    # ── Actions ──────────────────────────────────────────────────────────────

    def action_tab(self, index: int) -> None:
        if 0 <= index < len(TABS):
            self.tab = index
        else:
            self.flash(f"no tab {index + 1}", error=True)

    def action_next_tab(self) -> None:
        self.tab = (self.tab + 1) % len(TABS)

    def action_prev_tab(self) -> None:
        self.tab = (self.tab - 1) % len(TABS)

    def action_toggle_working(self) -> None:
        self.show_working = not self.show_working

    def action_cycle_unit(self) -> None:
        units = list(AngleUnit)
        self.angle_unit = units[(units.index(self.angle_unit) + 1) % len(units)]

    def action_pick_ellipsoid(self) -> None:
        def picked(name: str | None) -> None:
            if name is not None:
                self.ellipsoid = ELLIPSOIDS[name]

        self.push_screen(EllipsoidPicker(self.ellipsoid.name), picked)

    def action_yank(self, all_rows: bool) -> None:
        panel = self.active_calculator.query_one(ResultPanel)
        if not panel.rows:
            self.flash("nothing to yank", error=True)
            return
        if all_rows:
            text = "\n".join(f"{row.label} = {row.copy}" for row in panel.rows)
            message = f"yanked {len(panel.rows)} results"
        else:
            row = panel.selected_row
            assert row is not None
            text = row.copy
            message = f"yanked {row.label} = {row.copy}"
        self.last_yank = text
        self.copy_to_clipboard(text)
        self.flash(message)

    def action_help(self) -> None:
        self.push_screen(HelpScreen())

    def action_command(self) -> None:
        self._return_focus = self.focused
        self.mode = Mode.COMMAND
        cmdline = self.query_one(CommandLine)
        cmdline.value = ""
        cmdline.reset_completion()
        cmdline.focus()

    def _leave_command(self) -> None:
        self.mode = Mode.NORMAL
        target = self._return_focus
        if target is None or not target.is_attached or not target.focusable:
            target = self.active_calculator.form
        target.focus()

    # ── Commands ─────────────────────────────────────────────────────────────

    def run_command(self, text: str) -> None:
        name, args = split_command(text)
        try:
            self._dispatch(name, args)
        except CommandError as exc:
            self.flash(str(exc), error=True)

    def _dispatch(self, name: str, args: list[str]) -> None:
        if not name:
            return
        if name.isdigit():
            self.action_tab(int(name) - 1)
        elif name in ("q", "q!", "qa", "quit", "x", "wq"):
            self.exit()
        elif name in ("ell", "ellipsoid"):
            if not args:
                self.action_pick_ellipsoid()
                return
            self.ellipsoid = match_ellipsoid(" ".join(args))
            self.flash(f"ellipsoid: {self.ellipsoid.name}")
        elif name in ("units", "unit"):
            if not args:
                raise CommandError("usage: :units dms|deg|rad")
            self.angle_unit = match_unit(args[0])
        elif name == "method":
            self._set_choice("method", args)
        elif name in ("working", "work"):
            self.action_toggle_working()
        elif name == "tab":
            if not args or not args[0].isdigit():
                raise CommandError("usage: :tab <number>")
            self.action_tab(int(args[0]) - 1)
        elif name in ("help", "h"):
            self.action_help()
        else:
            raise CommandError(f"unknown command: {name}")

    def _set_choice(self, key: str, args: list[str]) -> None:
        fields = [f for f in self.active_calculator.form.fields if f.key == key]
        if not fields or not isinstance(fields[0], ChoiceField):
            raise CommandError(f"this tab has no {key} option")
        field = fields[0]
        options = ", ".join(label.lower() for label, _ in field.options)
        if not args or not field.select_label(args[0]):
            raise CommandError(f"usage: :{key} {options.replace(', ', '|')}")
        self.flash(f"{key}: {field.value_label}")
