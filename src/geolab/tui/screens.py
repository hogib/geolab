from rich.table import Table
from rich.text import Text
from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Vertical
from textual.screen import ModalScreen
from textual.widgets import OptionList, Static
from textual.widgets.option_list import Option

from ..ellipsoids import ELLIPSOIDS, Ellipsoid


class VimOptionList(OptionList):
    BINDINGS = [
        Binding("j", "cursor_down", "Down", show=False),
        Binding("k", "cursor_up", "Up", show=False),
        Binding("g", "first", "First", show=False),
        Binding("G", "last", "Last", show=False),
    ]


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
        ("1 – 8", "switch tab"),
        ("gt / gT  ] / [", "next / previous tab"),
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
        ("f", "fit the view to the points"),
        ("0", "whole world"),
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
        (":<n>  :tab <n>", "switch to tab n"),
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
