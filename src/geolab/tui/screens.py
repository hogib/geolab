from rich.table import Table
from rich.text import Text
from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal, Vertical
from textual.screen import ModalScreen
from textual.widgets import Input, OptionList, Static
from textual.widgets.option_list import Option

from ..ellipsoids import ELLIPSOIDS, Ellipsoid
from .fuzzy import fuzzy_filter


class VimOptionList(OptionList):
    BINDINGS = [
        Binding("j", "cursor_down", "Down", show=False),
        Binding("k", "cursor_up", "Up", show=False),
        Binding("g", "first", "First", show=False),
        Binding("G", "last", "Last", show=False),
    ]


def tab_key(index: int) -> str:
    """The key that jumps to calculator ``index`` (0-based): 1–9, then 0 for the tenth."""
    return "0" if index == 9 else str(index + 1)


def highlight(text: str, positions: tuple[int, ...], style: str) -> Text:
    out = Text(text)
    for i in positions:
        out.stylize(style, i, i + 1)
    return out


class Finder(ModalScreen[int | None]):
    """fzf-style calculator finder: type to fuzzy-filter, ↑/↓ or ctrl+n/p to move,
    enter to open. Returns the chosen calculator's index, or None."""

    BINDINGS = [
        Binding("escape,ctrl+c", "dismiss(None)", "Cancel", show=False),
        Binding("down,ctrl+n,ctrl+j", "move(1)", "Down", show=False),
        Binding("up,ctrl+p,ctrl+k", "move(-1)", "Up", show=False),
    ]

    def __init__(self, labels: list[str], current: int):
        super().__init__()
        self.labels = labels
        self.current = current

    def compose(self) -> ComposeResult:
        with Vertical(classes="block dialog", id="finder"):
            with Horizontal(id="finder-prompt"):
                yield Static("> ", id="finder-caret")
                yield Input(compact=True, id="finder-input")
                yield Static("", id="finder-count")
            yield OptionList(id="finder-list")

    def on_mount(self) -> None:
        self.query_one(".dialog").border_title = "Find calculator"
        self.query_one(".dialog").border_subtitle = "↑/↓ move · enter open · esc cancel"
        self._fill("")
        self.query_one(OptionList).highlighted = self.current
        self.query_one(Input).focus()

    def _fill(self, query: str) -> None:
        options = self.query_one(OptionList)
        options.clear_options()
        accent = f"bold {self.app.theme_variables.get('accent', 'yellow')}"
        matches = fuzzy_filter(query, self.labels)
        for m in matches:
            label = Text.assemble((f"{tab_key(m.index) if m.index < 10 else ' '} ", "dim"))
            label.append_text(highlight(self.labels[m.index], m.positions, accent))
            options.add_option(Option(label, id=str(m.index)))
        if options.option_count:
            options.highlighted = 0
        self.query_one("#finder-count", Static).update(
            Text(f" {len(matches)}/{len(self.labels)}", style="dim")
        )

    def on_input_changed(self, event: Input.Changed) -> None:
        self._fill(event.value)

    def on_input_submitted(self, event: Input.Submitted) -> None:
        options = self.query_one(OptionList)
        if options.highlighted is not None:
            self.dismiss(int(options.get_option_at_index(options.highlighted).id or 0))

    def on_option_list_option_selected(self, event: OptionList.OptionSelected) -> None:
        self.dismiss(int(event.option.id or 0))

    def action_move(self, step: int) -> None:
        options = self.query_one(OptionList)
        if options.option_count:
            current = options.highlighted or 0
            options.highlighted = max(0, min(options.option_count - 1, current + step))


class EllipsoidPicker(ModalScreen[str | None]):
    """Pick an ellipsoid with j/k and Enter. Returns its name, or None if cancelled.

    ``d`` deletes the highlighted custom ellipsoid.
    """

    BINDINGS = [
        Binding("escape,q", "dismiss(None)", "Cancel", show=False),
        Binding("d", "delete", "Delete", show=False),
    ]

    def __init__(self, current: str, ellipsoids: dict[str, Ellipsoid]):
        super().__init__()
        self.current = current
        self.ellipsoids = ellipsoids

    def compose(self) -> ComposeResult:
        options = []
        width = max(12, *(len(n) for n in self.ellipsoids))
        for name, ell in self.ellipsoids.items():
            label = Text.assemble(
                (f"{name:<{width}}", "bold"),
                (f"  a = {ell.a:>13,.3f} m   1/f = {ell.inv_f}".replace(",", " "), "dim"),
                ("" if name in ELLIPSOIDS else "   custom", "italic"),
            )
            options.append(Option(label, id=name))
        with Vertical(classes="block dialog"):
            yield VimOptionList(*options, id="ellipsoid-list")

    def on_mount(self) -> None:
        self.query_one(".dialog").border_title = "Ellipsoid"
        self.query_one(".dialog").border_subtitle = "j/k move · enter select · d delete · esc cancel"
        option_list = self.query_one(VimOptionList)
        option_list.highlighted = list(self.ellipsoids).index(self.current)
        option_list.focus()

    def on_option_list_option_selected(self, event: OptionList.OptionSelected) -> None:
        self.dismiss(event.option.id)

    def action_delete(self) -> None:
        option_list = self.query_one(VimOptionList)
        if option_list.highlighted is None:
            return
        name = option_list.get_option_at_index(option_list.highlighted).id
        assert name is not None
        if self.app.remove_ellipsoid(name):  # type: ignore[attr-defined]
            option_list.remove_option(name)


HELP = [
    ("Normal mode", [
        ("j / k", "move between fields"),
        ("i  a  enter", "edit the selected field"),
        ("c", "clear the field and edit it"),
        ("h / l  space", "cycle an option field"),
        ("t  /", "find a calculator (fuzzy search)"),
        ("1 – 9, 0", "jump to one of the first ten calculators"),
        ("gt / gT  ] / [", "next / previous calculator"),
        ("tab / shift+tab", "move between Input, Result and Working"),
        ("y / Y", "copy the selected result / all results"),
        ("w", "show or hide the working"),
        ("u", "cycle angle units (DMS, degrees, radians)"),
        ("e", "pick the ellipsoid"),
        (":", "command line"),
        ("?", "this help"),
        ("q", "quit"),
    ]),
    ("Insert mode", [
        ("esc", "back to normal mode"),
        ("enter", "confirm and move to the next field"),
    ]),
    ("Ellipsoid tab", [
        ("s  :w  :save", "save the edited parameters as a custom ellipsoid and use it"),
        ("", "custom ellipsoids are kept in ~/.config/geolab/ellipsoids.toml"),
        ("d (in the picker)", "delete a custom ellipsoid"),
    ]),
    ("Map tab (World block)", [
        ("h j k l  arrows", "pan"),
        ("+ / -", "zoom in / out"),
        ("f", "fit the view to the points and follow them (the default)"),
        ("0", "whole world"),
        ("", "the map refits when the points change, unless you have panned or zoomed"),
    ]),
    ("Result / Working blocks", [
        ("j / k", "select a result / scroll"),
        ("g / G  ctrl+d / ctrl+u", "top / bottom, page down / up"),
    ]),
    ("Commands", [
        ("tab / shift+tab", "complete the command or argument, cycling through matches"),
        (":ell <name>", "set the ellipsoid, e.g. :ell grs80 (no name opens the picker)"),
        (":units <dms|deg|rad>", "set the angle units"),
        (":method <name>", "set the method where there is one, e.g. :method bowring"),
        (":working", "show or hide the working"),
        (":save  :w", "save a custom ellipsoid (Ellipsoid tab)"),
        (":tab <name>", "open the best fuzzy match for name"),
        (":help  :q", "help, quit"),
    ]),
    ("Angle input", [
        ("41 00 30.5 N", "DMS with spaces, or 41°00'30.5\"N, or 41:00:30.5"),
        ("41 30.5", "degrees and decimal minutes"),
        ("-41.5083", "decimal degrees; S and W are negative"),
    ]),
]


class HelpScreen(ModalScreen[None]):
    BINDINGS = [Binding("escape,q,question_mark", "dismiss(None)", "Close", show=False)]

    def compose(self) -> ComposeResult:
        table = Table.grid(padding=(0, 3))
        table.add_column(style="bold", no_wrap=True)
        table.add_column()
        for i, (title, rows) in enumerate(HELP):
            if i:
                table.add_row("", "")
            table.add_row(Text(title, style="bold underline"), "")
            for keys, description in rows:
                table.add_row(keys, description)
        with Vertical(classes="block dialog"):
            yield Static(table)

    def on_mount(self) -> None:
        self.query_one(".dialog").border_title = "Help"
        self.query_one(".dialog").border_subtitle = "esc / q / ? close"
