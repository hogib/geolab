import math

import pytest

pytest.importorskip("textual")

from geolab import GRS80, HAYFORD, WGS84, AngleUnit, Geodetic, parse_angle  # noqa: E402
from geolab.tui.app import GeolabApp, Mode  # noqa: E402
from geolab.tui.screens import EllipsoidPicker, HelpScreen  # noqa: E402
from geolab.tui.widgets import (  # noqa: E402
    ChoiceField,
    Form,
    ResultPanel,
    TextField,
    WorkingPanel,
    CommandLine,
    WorkingView,
    format_length,
)


@pytest.fixture
async def pilot():
    app = GeolabApp()
    async with app.run_test(size=(130, 40)) as pilot:
        await pilot.pause()
        yield pilot


def results(app: GeolabApp) -> ResultPanel:
    return app.active_calculator.query_one(ResultPanel)


def rows(app: GeolabApp) -> dict[str, str]:
    return {row.label: row.display for row in results(app).rows}


def form(app: GeolabApp) -> Form:
    return app.active_calculator.form


def text_field(app: GeolabApp, key: str) -> TextField:
    field = app.active_calculator.field(key)
    assert isinstance(field, TextField)
    return field


async def run(pilot, command: str) -> None:
    await pilot.press("colon", *command, "enter")
    await pilot.pause()


class TestStartup:
    async def test_normal_mode_with_form_focused(self, pilot):
        app = pilot.app
        assert app.mode == Mode.NORMAL
        assert app.tab == 0
        assert isinstance(app.focused, Form)

    async def test_default_result_matches_core(self, pilot):
        expected = Geodetic(
            parse_angle("41 00 30 N"), parse_angle("29 00 00 E"), 150.0
        ).to_cartesian(WGS84)
        assert rows(pilot.app) == {
            "X": format_length(expected.x),
            "Y": format_length(expected.y),
            "Z": format_length(expected.z),
        }

    async def test_working_hidden(self, pilot):
        assert pilot.app.active_calculator.query_one(WorkingView).region.width == 0

    async def test_background_is_transparent(self, pilot):
        # Blank cells use the terminal's default background colour.
        strips = pilot.app.screen._compositor.render_strips()
        backgrounds = {seg.style.bgcolor for seg in strips[20] if seg.style and seg.style.bgcolor}
        assert {c.name for c in backgrounds} == {"default"}


class TestNormalMode:
    async def test_j_k_move_selection(self, pilot):
        await pilot.press("j", "j")
        assert form(pilot.app).selected == 2
        await pilot.press("j")  # clamps at the last field
        assert form(pilot.app).selected == 2
        await pilot.press("k", "k", "k")
        assert form(pilot.app).selected == 0

    async def test_selected_field_is_marked(self, pilot):
        await pilot.press("j")
        fields = form(pilot.app).fields
        assert [f.has_class("-selected") for f in fields] == [False, True, False]

    async def test_letters_do_not_type_into_fields(self, pilot):
        before = text_field(pilot.app, "lat").value
        await pilot.press("x", "z")
        assert text_field(pilot.app, "lat").value == before


class TestInsertMode:
    async def test_i_edits_selected_field(self, pilot):
        await pilot.press("j", "i")
        app = pilot.app
        assert app.mode == Mode.INSERT
        assert app.focused is text_field(app, "lon").input

    async def test_typing_updates_results_live(self, pilot):
        await pilot.press("j", "j", "c", *"1000")
        await pilot.pause()
        expected = Geodetic(
            parse_angle("41 00 30 N"), parse_angle("29 00 00 E"), 1000.0
        ).to_cartesian(WGS84)
        assert rows(pilot.app)["X"] == format_length(expected.x)

    async def test_escape_returns_to_normal(self, pilot):
        await pilot.press("i", "escape")
        await pilot.pause()
        app = pilot.app
        assert app.mode == Mode.NORMAL
        assert isinstance(app.focused, Form)
        assert not text_field(app, "lat").input.can_focus

    async def test_enter_confirms_and_moves_to_next_field(self, pilot):
        await pilot.press("i", "enter")
        await pilot.pause()
        assert pilot.app.mode == Mode.NORMAL
        assert form(pilot.app).selected == 1

    async def test_shortcut_letters_are_typed_in_insert_mode(self, pilot):
        await pilot.press("c", *"qwue")
        app = pilot.app
        assert app.is_running
        assert text_field(app, "lat").value == "qwue"
        assert app.angle_unit == AngleUnit.DMS

    async def test_c_clears_field_and_shows_hint(self, pilot):
        await pilot.press("j", "c")
        await pilot.pause()
        app = pilot.app
        assert text_field(app, "lon").value == ""
        assert results(app).message == "Enter λ longitude"
        assert results(app).has_class("-hint")

    async def test_invalid_input_shows_error(self, pilot):
        await pilot.press("c", *"41 61 0")
        await pilot.pause()
        app = pilot.app
        assert results(app).has_class("-error")
        assert "minutes must be less than 60" in results(app).message
        assert text_field(app, "lat").has_class("-invalid")

    async def test_tab_out_of_insert_returns_to_normal(self, pilot):
        await pilot.press("i", "tab")
        await pilot.pause()
        assert pilot.app.mode == Mode.NORMAL
        assert isinstance(pilot.app.focused, ResultPanel)


class TestTabs:
    async def test_number_keys(self, pilot):
        await pilot.press("3")
        assert pilot.app.tab == 2
        await pilot.press("1")
        assert pilot.app.tab == 0

    async def test_brackets_wrap(self, pilot):
        await pilot.press("left_square_bracket")
        assert pilot.app.tab == 2
        await pilot.press("right_square_bracket")
        assert pilot.app.tab == 0

    async def test_gt_and_gT(self, pilot):
        await pilot.press("g", "t")
        assert pilot.app.tab == 1
        await pilot.press("g", "T")
        assert pilot.app.tab == 0

    async def test_switching_focuses_new_form(self, pilot):
        await pilot.press("2")
        await pilot.pause()
        assert pilot.app.focused is form(pilot.app)

    async def test_missing_tab(self, pilot):
        await pilot.press("9")
        assert pilot.app.tab == 0
        assert pilot.app.last_message == "no tab 9"


class TestChoiceField:
    async def test_l_and_h_cycle_method(self, pilot):
        await pilot.press("2", "j", "j", "j")
        method = pilot.app.active_calculator.field("method")
        assert isinstance(method, ChoiceField)
        await pilot.press("l")
        assert method.value_label == "Bowring"
        trace = pilot.app.active_calculator.query_one(WorkingPanel).trace
        assert "θ" in [s.name for s in trace.steps]
        await pilot.press("h")
        assert method.value_label == "Iterative"

    async def test_ecef_to_geodetic_result(self, pilot):
        await pilot.press("2")
        await pilot.pause()
        assert rows(pilot.app)["φ"] == "41°00'30.0000\"N"
        assert rows(pilot.app)["λ"] == "29°00'00.0000\"E"


class TestGlobalKeys:
    async def test_u_cycles_units(self, pilot):
        await pilot.press("3")
        await pilot.pause()
        assert rows(pilot.app)["φ  geodetic"] == "41°00'30.0000\"N"
        await pilot.press("u")
        await pilot.pause()
        assert pilot.app.angle_unit == AngleUnit.DEGREES
        assert rows(pilot.app)["φ  geodetic"] == "41.00833333°"
        await pilot.press("u", "u")
        assert pilot.app.angle_unit == AngleUnit.DMS

    async def test_w_toggles_working(self, pilot):
        await pilot.press("w")
        await pilot.pause()
        app = pilot.app
        assert app.show_working
        view = app.active_calculator.query_one(WorkingView)
        assert view.region.width > 0
        assert [s.name for s in view.query_one(WorkingPanel).trace.steps] == ["e²", "N", "X", "Y", "Z"]

    async def test_hiding_working_moves_focus_back_to_form(self, pilot):
        await pilot.press("w", "tab", "tab")
        await pilot.pause()
        assert isinstance(pilot.app.focused, WorkingView)
        await pilot.press("w")
        await pilot.pause()
        assert isinstance(pilot.app.focused, Form)

    async def test_ellipsoid_change_recalculates(self, pilot):
        before = rows(pilot.app)["X"]
        pilot.app.ellipsoid = HAYFORD
        await pilot.pause()
        assert rows(pilot.app)["X"] != before

    async def test_help(self, pilot):
        await pilot.press("question_mark")
        await pilot.pause()
        assert isinstance(pilot.app.screen, HelpScreen)
        await pilot.press("escape")
        await pilot.pause()
        assert not isinstance(pilot.app.screen, HelpScreen)

    async def test_q_quits(self, pilot):
        await pilot.press("q")
        await pilot.pause()
        assert not pilot.app.is_running


class TestEllipsoidPicker:
    async def test_pick_with_j_and_enter(self, pilot):
        await pilot.press("e")
        await pilot.pause()
        assert isinstance(pilot.app.screen, EllipsoidPicker)
        await pilot.press("j", "enter")
        await pilot.pause()
        assert pilot.app.ellipsoid is GRS80
        assert not isinstance(pilot.app.screen, EllipsoidPicker)

    async def test_escape_cancels(self, pilot):
        await pilot.press("e", "j", "j", "escape")
        await pilot.pause()
        assert pilot.app.ellipsoid is WGS84


class TestYank:
    async def test_y_yanks_selected_row(self, pilot):
        await pilot.press("y")
        assert pilot.app.last_yank == results(pilot.app).rows[0].copy
        assert pilot.app.last_message.startswith("yanked X = ")

    async def test_select_row_in_result_panel(self, pilot):
        await pilot.press("tab", "j", "y")
        assert pilot.app.last_yank == results(pilot.app).rows[1].copy

    async def test_Y_yanks_all(self, pilot):
        await pilot.press("Y")
        lines = pilot.app.last_yank.splitlines()
        assert [line.split(" = ")[0] for line in lines] == ["X", "Y", "Z"]

    async def test_yanked_length_is_a_plain_number(self, pilot):
        await pilot.press("y")
        assert float(pilot.app.last_yank) == pytest.approx(4215751.8328, abs=1e-3)

    async def test_nothing_to_yank(self, pilot):
        await pilot.press("c", "escape", "y")
        assert pilot.app.last_message == "nothing to yank"


class TestCommandLine:
    async def test_enters_command_mode(self, pilot):
        await pilot.press("colon")
        await pilot.pause()
        assert pilot.app.mode == Mode.COMMAND
        assert pilot.app.focused.id == "cmdline"

    async def test_escape_cancels(self, pilot):
        await pilot.press("colon", *"ell hay", "escape")
        await pilot.pause()
        assert pilot.app.mode == Mode.NORMAL
        assert pilot.app.ellipsoid is WGS84
        assert isinstance(pilot.app.focused, Form)

    async def test_ellipsoid(self, pilot):
        await run(pilot, "ell hayford")
        assert pilot.app.ellipsoid is HAYFORD
        assert pilot.app.mode == Mode.NORMAL

    async def test_ellipsoid_without_name_opens_picker(self, pilot):
        await run(pilot, "ell")
        assert isinstance(pilot.app.screen, EllipsoidPicker)

    async def test_units(self, pilot):
        await run(pilot, "units rad")
        assert pilot.app.angle_unit == AngleUnit.RADIANS

    async def test_method(self, pilot):
        await pilot.press("2")
        await run(pilot, "method bow")
        assert pilot.app.active_calculator.choice("method").name == "BOWRING"
        assert pilot.app.last_message == "method: Bowring"

    async def test_method_on_tab_without_one(self, pilot):
        await run(pilot, "method bowring")
        assert pilot.app.last_message == "this tab has no method option"

    async def test_tab_number(self, pilot):
        await run(pilot, "3")
        assert pilot.app.tab == 2
        await run(pilot, "tab 2")
        assert pilot.app.tab == 1

    async def test_working(self, pilot):
        await run(pilot, "working")
        assert pilot.app.show_working

    async def test_unknown_command(self, pilot):
        await run(pilot, "frobnicate")
        assert pilot.app.last_message == "unknown command: frobnicate"

    async def test_error_message(self, pilot):
        await run(pilot, "units grad")
        assert "unknown unit" in pilot.app.last_message
        assert pilot.app.angle_unit == AngleUnit.DMS

    async def test_quit(self, pilot):
        await run(pilot, "q")
        assert not pilot.app.is_running

    async def test_returns_focus_to_previous_block(self, pilot):
        await pilot.press("tab")
        await pilot.pause()
        await run(pilot, "units deg")
        assert isinstance(pilot.app.focused, ResultPanel)


class TestLatitudeTab:
    async def test_all_three_types(self, pilot):
        await pilot.press("3", "u")
        await pilot.pause()
        values = {k: float(v.rstrip("°")) for k, v in rows(pilot.app).items()}
        phi = math.radians(41 + 30 / 3600)
        assert values["φ  geodetic"] == pytest.approx(41.00833333, abs=1e-8)
        assert values["ψ  geocentric"] == pytest.approx(
            math.degrees(math.atan((1 - WGS84.e_sq) * math.tan(phi))), abs=1e-8
        )
        assert values["β  parametric"] == pytest.approx(
            math.degrees(math.atan(math.sqrt(1 - WGS84.e_sq) * math.tan(phi))), abs=1e-8
        )


class TestCompletion:
    async def test_tab_completes_command(self, pilot):
        await pilot.press("colon", "u", "tab")
        await pilot.pause()
        app = pilot.app
        assert app.query_one(CommandLine).value == "units "
        assert app.focused.id == "cmdline"

    async def test_tab_cycles_and_shift_tab_goes_back(self, pilot):
        await pilot.press("colon", *"ell ", "tab", "tab", "tab")
        cmdline = pilot.app.query_one(CommandLine)
        assert cmdline.value == "ell Hayford"
        await pilot.press("shift+tab")
        assert cmdline.value == "ell GRS80"

    async def test_unique_command_then_argument(self, pilot):
        await pilot.press("colon", "e", "tab", "tab")
        assert pilot.app.query_one(CommandLine).value == "ell WGS84"

    async def test_typing_restarts_completion(self, pilot):
        await pilot.press("colon", *"ell ", "tab")
        await pilot.press(*["backspace"] * 5, *"kr", "tab")
        await pilot.pause()
        assert pilot.app.query_one(CommandLine).value == "ell Krassovsky"

    async def test_completed_command_runs(self, pilot):
        await pilot.press("colon", *"ell hay", "tab", "enter")
        await pilot.pause()
        assert pilot.app.ellipsoid is HAYFORD

    async def test_method_completion_uses_current_tab(self, pilot):
        await pilot.press("2", "colon", *"method ", "tab")
        cmdline = pilot.app.query_one(CommandLine)
        assert cmdline.candidates == ["method iterative", "method bowring"]

    async def test_wildmenu_in_status_line(self, pilot):
        await pilot.press("colon", *"units ", "tab")
        await pilot.pause()
        status = pilot.app.query_one("#statusline")
        line = "".join(seg.text for seg in status.render_lines(status.region.reset_offset)[0])
        assert "dms" in line and "deg" in line and "rad" in line

    async def test_no_candidates_keeps_text(self, pilot):
        await pilot.press("colon", *"zz", "tab")
        assert pilot.app.query_one(CommandLine).value == "zz"

    async def test_reopening_clears_candidates(self, pilot):
        await pilot.press("colon", *"units ", "tab", "escape", "colon")
        await pilot.pause()
        assert pilot.app.query_one(CommandLine).candidates == []
