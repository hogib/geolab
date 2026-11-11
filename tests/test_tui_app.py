import math

import pytest

pytest.importorskip("textual")

from geolab import (  # noqa: E402
    GRS80,
    HAYFORD,
    WGS84,
    AngleUnit,
    Ellipsoid,
    Geodetic,
    HelmertParams,
    direct,
    ellipse_perimeter,
    gaussian_mean_radius,
    inverse,
    latitude_from_meridian_arc,
    latitude_from_meridian_distance,
    longitude_from_parallel_arc,
    meridian_arc_length,
    meridian_distance,
    meridian_radius,
    parallel_arc_length,
    parse_angle,
    transform_geodetic,
)
from geolab.area import quadrangle_area  # noqa: E402
from geolab.config import save_ellipsoids  # noqa: E402
from geolab.gauss_kruger import krueger_forward  # noqa: E402
from geolab.soldner import soldner_from_geodetic  # noqa: E402
from geolab.tui.app import TABS, GeolabApp, Mode  # noqa: E402
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
        assert pilot.app.tab == len(TABS) - 1
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

    async def test_finder_fuzzy_match(self, pilot):
        await pilot.press("t", *"mpa", "enter")  # Meridian and PArallel arcs
        await pilot.pause()
        assert pilot.app.tab == 8

    async def test_finder_arrows_and_slash(self, pilot):
        await pilot.press("slash", *"ell", "down", "enter")  # Ellipsoid, then 2D ellipse
        await pilot.pause()
        assert pilot.app.tab == 9

    async def test_finder_shows_count_and_highlights(self, pilot):
        await pilot.press("t", *"vin")
        await pilot.pause()
        finder = pilot.app.screen
        shown, total = str(finder.query_one("#finder-count").render()).strip().split("/")
        assert int(total) == len(TABS) and 2 <= int(shown) < len(TABS)
        first = finder.query_one("#finder-list").get_option_at_index(0).prompt
        assert first.plain.endswith("Inverse problem (Vincenty)")
        assert any(span.style for span in first.spans)

    async def test_finder_escape_cancels(self, pilot):
        await pilot.press("t", *"map", "escape")
        await pilot.pause()
        assert pilot.app.tab == 0
        assert pilot.app.focused is form(pilot.app)

    async def test_tab_command_fuzzy(self, pilot):
        await run(pilot, "tab map")
        assert pilot.app.tab == 7
        await run(pilot, "tab zzz")
        assert pilot.app.last_message == "no calculator matching 'zzz'"

    async def test_missing_tab(self, pilot):
        await run(pilot, "99")
        assert pilot.app.tab == 0
        assert pilot.app.last_message == "no tab 99"


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
        values = {k: float(v.rstrip("°")) for k, v in rows(pilot.app).items() if v.endswith("°")}
        phi = math.radians(41 + 30 / 3600)
        assert values["φ  geodetic"] == pytest.approx(41.00833333, abs=1e-8)
        assert values["ψ  geocentric"] == pytest.approx(
            math.degrees(math.atan((1 - WGS84.e_sq) * math.tan(phi))), abs=1e-8
        )
        assert values["β  parametric"] == pytest.approx(
            math.degrees(math.atan(math.sqrt(1 - WGS84.e_sq) * math.tan(phi))), abs=1e-8
        )

    async def test_meridian_ellipse_coordinates_course_example(self, pilot):
        await run(pilot, "ell hayford")
        await pilot.press("3", "c", *"35 20 40", "escape")
        await pilot.pause()
        got = values(pilot.app)
        assert float(got["p   meridian ellipse"]) == pytest.approx(5208643.936, abs=2e-3)
        assert float(got["z   meridian ellipse"]) == pytest.approx(3669161.169, abs=2e-3)


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


def values(app: GeolabApp) -> dict[str, str]:
    """Result rows by label, as the text that would be yanked."""
    return {row.label: row.copy for row in results(app).rows}


class TestInverseTab:
    async def test_matches_core(self, pilot):
        await pilot.press("4", "u")
        await pilot.pause()
        expected = inverse(*map(parse_angle, ("41 00 30 N", "29 00 00 E", "39 55 00 N", "32 51 00 E")), WGS84)
        got = values(pilot.app)
        assert float(got["s   distance"]) == pytest.approx(expected.distance, abs=1e-4)
        assert float(got["α₁  azimuth at P₁"].rstrip("°")) == pytest.approx(math.degrees(expected.azimuth1), abs=1e-8)
        assert float(got["α₂₁ back azimuth"].rstrip("°")) == pytest.approx(math.degrees(expected.back_azimuth), abs=1e-8)
        assert got["iterations"] == str(expected.iterations)

    async def test_antipodal_shows_error(self, pilot):
        await pilot.press("4", "c", *"0", "enter", "c", *"0", "enter", "c", *"0.5", "enter", "c", *"179.7", "enter")
        await pilot.pause()
        assert results(pilot.app).has_class("-error")
        assert "antipodal" in results(pilot.app).message


class TestDirectTab:
    async def test_defaults_land_on_inverse_end_point(self, pilot):
        await pilot.press("5")
        await pilot.pause()
        got = rows(pilot.app)
        assert got["φ₂  latitude"] == "39°55'00.0000\"N"
        assert got["λ₂  longitude"] == "32°51'00.0000\"E"

    async def test_matches_core(self, pilot):
        await pilot.press("5", "j", "j", "c", *"45", "enter", "c", *"1000000", "enter", "u")
        await pilot.pause()
        expected = direct(parse_angle("41 00 30 N"), parse_angle("29 E"), math.radians(45), 1e6, WGS84)
        got = values(pilot.app)
        assert float(got["φ₂  latitude"].rstrip("°")) == pytest.approx(math.degrees(expected.lat2), abs=1e-8)
        assert float(got["λ₂  longitude"].rstrip("°")) == pytest.approx(math.degrees(expected.lon2), abs=1e-8)


class TestRadiiTab:
    async def test_matches_core(self, pilot):
        await pilot.press("6")
        await pilot.pause()
        lat = parse_angle("41 00 30 N")
        got = values(pilot.app)
        assert float(got["M   meridian"]) == pytest.approx(meridian_radius(lat, WGS84), abs=1e-4)
        assert float(got["R   Gaussian mean"]) == pytest.approx(gaussian_mean_radius(lat, WGS84), abs=1e-4)

    async def test_working_has_a_section_per_radius(self, pilot):
        await pilot.press("6")
        await pilot.pause()
        trace = pilot.app.active_calculator.query_one(WorkingPanel).trace
        titles = [e.title for e in trace.entries if not hasattr(e, "value")]
        assert titles[:5] == [
            "Meridian radius M",
            "Prime vertical radius N",
            "Parallel radius r",
            "Gaussian mean radius R",
            "Radius in azimuth α",
        ]
        assert "First fundamental quantities" in titles

    async def test_fundamental_forms(self, pilot):
        await pilot.press("6")
        await pilot.pause()
        lat = parse_angle("41 00 30 N")
        got = values(pilot.app)
        assert float(got["E   1st form (= M²)"]) == pytest.approx(meridian_radius(lat, WGS84) ** 2, rel=1e-12)
        assert float(got["F   1st form"]) == 0


class TestArcsTab:
    async def test_matches_core(self, pilot):
        await pilot.press("9")
        await pilot.pause()
        lat1, lat2 = parse_angle("41 00 30 N"), parse_angle("39 55 00 N")
        lon1, lon2 = parse_angle("29 E"), parse_angle("32 51 E")
        got = values(pilot.app)
        assert float(got["Δm  meridian arc φ₁→φ₂"]) == pytest.approx(
            meridian_arc_length(lat1, lat2, WGS84), abs=1e-4
        )
        assert float(got["s₂  parallel arc at φ₂"]) == pytest.approx(
            parallel_arc_length(lat2, lon1, lon2, WGS84), abs=1e-4
        )


async def test_header_shows_current_calculator(pilot):
    header = lambda: str(pilot.app.query_one("#tabbar").render())  # noqa: E731
    assert "Geodetic → ECEF" in header()
    await pilot.press("4")
    await pilot.pause()
    assert "Inverse problem (Vincenty)" in header()


class TestArcModes:
    async def test_lengths_mode_fields(self, pilot):
        await pilot.press("9")
        await pilot.pause()
        assert [f.key for f in form(pilot.app).fields] == ["mode", "lat1", "lon1", "lat2", "lon2"]

    async def test_latitude_from_meridian_arc(self, pilot):
        await pilot.press("9", "l")
        await pilot.pause()
        assert [f.key for f in form(pilot.app).fields] == ["mode", "lat1", "dm"]
        await pilot.press("j", "j", "c", *"-120000", "escape", "u")
        await pilot.pause()
        lat1 = parse_angle("41 00 30 N")
        expected = latitude_from_meridian_arc(lat1, -120000.0, WGS84)
        got = values(pilot.app)
        assert float(got["φ₂  latitude"].rstrip("°")) == pytest.approx(math.degrees(expected.lat), abs=1e-8)

    async def test_longitude_from_parallel_arc(self, pilot):
        await pilot.press("9", "l", "l")
        await pilot.pause()
        assert [f.key for f in form(pilot.app).fields] == ["mode", "lat1", "lon1", "s"]
        await pilot.press("u")
        await pilot.pause()
        lat1, lon1 = parse_angle("41 00 30 N"), parse_angle("29 E")
        expected = longitude_from_parallel_arc(lat1, lon1, 325000.0, WGS84)
        got = values(pilot.app)
        assert float(got["λ₂  longitude"].rstrip("°")) == pytest.approx(math.degrees(expected), abs=1e-8)


class TestEllipseTab:
    async def test_quarter_perimeter_by_default(self, pilot):
        await pilot.press("0")
        await pilot.pause()
        assert pilot.app.tab == 9
        assert [f.key for f in form(pilot.app).fields] == ["mode", "a", "b", "t1", "t2", "angle", "axes"]
        got = values(pilot.app)
        assert float(got["s   arc length"]) == pytest.approx(ellipse_perimeter(5, 3) / 4, abs=1e-12)
        assert float(got["c   chord"]) == pytest.approx(math.hypot(5, 3))

    async def test_polar_angles(self, pilot):
        await pilot.press("0", "j", "j", "j", "j", "c", *"45", "escape", "j", "l")
        await pilot.pause()
        got = values(pilot.app)
        # θ = 0 → (5, 0); θ = 45° → x = y = 15/√34
        assert float(got["x₂"]) == pytest.approx(float(got["y₂"]))
        assert float(got["x₂"]) == pytest.approx(15 / math.sqrt(34))


    async def test_hayford_meridian_ellipse(self, pilot):
        await run(pilot, "ell hayford")
        await pilot.press("0", "j", "j", "j", "j", "c", *"41 00 30", "escape", "j", "l", "l", "j", "l")
        await pilot.pause()
        assert text_field(pilot.app, "a").value == repr(HAYFORD.a)
        assert text_field(pilot.app, "b").value == repr(HAYFORD.b)
        got = values(pilot.app)
        lat = parse_angle("41 00 30")
        assert float(got["s   arc length"]) == pytest.approx(meridian_distance(lat, HAYFORD), abs=1e-6)
        assert float(got["e²  1st eccentricity²"]) == pytest.approx(HAYFORD.e_sq, rel=1e-12)
        assert float(got["1/f"]) == pytest.approx(297.0, rel=1e-12)
        await run(pilot, "ell wgs84")  # follows the ellipsoid
        assert text_field(pilot.app, "b").value == repr(WGS84.b)

    async def test_choosing_an_ellipsoid_switches_from_custom_axes(self, pilot):
        await pilot.press("0")
        await pilot.pause()
        assert pilot.app.active_calculator.choice("axes") == "custom"
        await run(pilot, "ell hayford")
        assert pilot.app.active_calculator.choice("axes") == "ellipsoid"
        assert text_field(pilot.app, "a").value == repr(HAYFORD.a)
        assert text_field(pilot.app, "b").value == repr(HAYFORD.b)
        assert float(values(pilot.app)["1/f"]) == pytest.approx(297.0, rel=1e-12)

    async def test_point_mode_hayford(self, pilot):
        await run(pilot, "ell hayford")
        lat = parse_angle("41 00 30 N")
        m = meridian_distance(lat, HAYFORD)
        await pilot.press("0", "l")
        await pilot.pause()
        assert [f.key for f in form(pilot.app).fields] == ["mode", "a", "b", "m", "lon", "axes"]
        assert pilot.app.active_calculator.choice("axes") == "ellipsoid"
        assert text_field(pilot.app, "a").value == repr(HAYFORD.a)
        await pilot.press("j", "j", "j", "c", *repr(m), "escape", "j", "c", *"29 E", "escape", "u")
        await pilot.pause()
        got = values(pilot.app)
        assert float(got["φ   geodetic latitude"].rstrip("°")) == pytest.approx(41 + 30 / 3600, abs=1e-8)
        xyz = Geodetic(lat, parse_angle("29 E"), 0.0).to_cartesian(HAYFORD)
        assert float(got["X   ECEF"]) == pytest.approx(xyz.x, abs=1e-4)
        assert float(got["Z   ECEF"]) == pytest.approx(xyz.z, abs=1e-4)

    async def test_point_mode_custom_axes(self, pilot):
        # A custom meridian ellipse with Hayford's axes gives the same point.
        await pilot.press("0", "l", "j", "c", *repr(HAYFORD.a), "escape", "j", "c", *repr(HAYFORD.b), "escape")
        await pilot.press("j", "j", "j", "h")  # axes: back to Custom, keeping a and b
        await pilot.pause()
        assert pilot.app.active_calculator.choice("axes") == "custom"
        got = values(pilot.app)
        expected = latitude_from_meridian_distance(4540000.0, HAYFORD)
        assert float(got["y   meridian ellipse"]) == pytest.approx(expected.y, abs=1e-6)

    async def test_switching_back_restores_arc_fields(self, pilot):
        await pilot.press("0", "l", "j", "j", "j", "j", "j", "k", "k", "k", "k", "k", "l")
        await pilot.pause()
        assert [f.key for f in form(pilot.app).fields] == ["mode", "a", "b", "t1", "t2", "angle", "axes"]
        assert "s   arc length" in values(pilot.app)

    async def test_beyond_the_pole_is_an_error(self, pilot):
        await pilot.press("0", "l", "j", "j", "j", "c", *"10002000", "escape")
        await pilot.pause()
        assert "quarter meridian" in results(pilot.app).message
        assert results(pilot.app).has_class("-error")


class TestDatumCalculator:
    async def test_matches_core_and_round_trips(self, pilot):
        await run(pilot, "tab datum")
        lat, lon = parse_angle("41 00 30 N"), parse_angle("29 E")
        params = HelmertParams.from_units(87, 98, 121)
        expected = transform_geodetic(Geodetic(lat, lon, 150.0), GRS80, HAYFORD, params)
        got = values(pilot.app)
        assert float(got["h   on Hayford"]) == pytest.approx(expected.height, abs=1e-4)
        await pilot.press("l")  # ED-50 → TUTGA-99 with the same parameters
        await pilot.pause()
        assert "h   on GRS80" in values(pilot.app)


class TestAreaCalculator:
    async def test_map_sheet(self, pilot):
        await run(pilot, "ell hayford")
        await run(pilot, "tab area")
        lat1, lat2 = parse_angle("39 N"), parse_angle("39 07 30 N")
        lon1, lon2 = parse_angle("29 E"), parse_angle("29 07 30 E")
        got = values(pilot.app)
        assert float(got["P   quadrangle (m²)"]) == pytest.approx(
            quadrangle_area(lat1, lat2, lon1, lon2, HAYFORD), rel=1e-12
        )


class TestGeodesicCurveCalculator:
    async def test_course_example(self, pilot):
        await run(pilot, "ell hayford")
        await run(pilot, "tab clairaut")
        got = rows(pilot.app)
        assert values(pilot.app)["k   N cos φ sin α (m)"].startswith("3194567.525")
        assert got["φmax highest latitude"] == "60°01'40.6083\"N"
        assert got["α   at φ₂"] == "30°03'20.7081\""

    async def test_unreached_latitude(self, pilot):
        await run(pilot, "tab clairaut")
        await pilot.press("j", "j", "c", *"70 N", "escape")
        await pilot.pause()
        assert rows(pilot.app)["α at φ₂"] == "not reached"


class TestShortLinesCalculator:
    async def test_mid_latitude_course_example(self, pilot):
        await run(pilot, "ell bessel")
        await run(pilot, "tab short")
        await pilot.press("l", "l", "l")  # Mid-latitude inverse
        await pilot.pause()
        assert [f.key for f in form(pilot.app).fields] == ["mode", "lat1", "lon1", "lat2", "lon2"]
        got = values(pilot.app)
        assert float(got["s   distance"]) == pytest.approx(47652.596, abs=2e-3)
        assert rows(pilot.app)["α₁₂ azimuth"] == "25°16'31.9783\""

    async def test_legendre_direct_agrees_with_vincenty(self, pilot):
        await run(pilot, "ell bessel")
        await run(pilot, "tab short")
        got = values(pilot.app)
        assert abs(float(got["φ₂ − Vincenty (″)"])) < 1e-4
        assert abs(float(got["λ₂ − Vincenty (″)"])) < 1e-4


class TestTriangleCalculator:
    async def test_course_example(self, pilot):
        await run(pilot, "ell bessel")
        await run(pilot, "tab triangles")
        got = values(pilot.app)
        assert float(got["a   P₂P₃"]) == pytest.approx(36594.4366, abs=1e-3)
        assert float(got["ε   spherical excess (″)"]) == pytest.approx(4.1667, abs=1e-3)
        assert abs(float(got["Legendre − exact a (mm)"])) < 0.1

    async def test_plane_sss(self, pilot):
        await run(pilot, "tab triangles")
        await pilot.press("l", "l", "j", "l")  # KKK, plane
        await pilot.pause()
        assert [f.key for f in form(pilot.app).fields] == ["mode", "surface", "lat", "a", "b", "c"]
        got = values(pilot.app)
        assert "ε   spherical excess (″)" not in got
        assert float(got["F   area (m²)"]) > 0


class TestSoldnerCalculator:
    async def test_round_trip_through_modes(self, pilot):
        await run(pilot, "tab soldner")
        got = values(pilot.app)
        expected = soldner_from_geodetic(
            parse_angle("40 30 N"), parse_angle("31 E"), parse_angle("40 N"), parse_angle("30 E"), WGS84
        )
        assert float(got["x   north"]) == pytest.approx(expected.x, abs=1e-4)
        assert float(got["y   east"]) == pytest.approx(expected.y, abs=1e-4)

    async def test_inverse_mode_fields(self, pilot):
        await run(pilot, "tab soldner")
        await pilot.press("h")  # wraps to Inverse
        await pilot.pause()
        assert [f.key for f in form(pilot.app).fields] == ["mode", "lat0", "lon0", "x", "y", "x2", "y2"]
        assert "S   distance" in values(pilot.app)


class TestProjectionCalculator:
    async def test_tissot_default(self, pilot):
        await run(pilot, "tab projections")
        got = values(pilot.app)
        assert float(got["Φ   area scale m₁m₂"]) == pytest.approx(1.08)

    async def test_isometric_series_close_to_exact(self, pilot):
        await run(pilot, "tab projections")
        await pilot.press("l")
        await pilot.pause()
        assert [f.key for f in form(pilot.app).fields] == ["mode", "lat", "lat0"]
        assert abs(float(values(pilot.app)["series − exact"])) < 1e-10


class TestGaussKrugerCalculator:
    async def test_forward_default(self, pilot):
        await run(pilot, "tab transverse")
        got = values(pilot.app)
        k = krueger_forward(parse_angle("40 42 N"), parse_angle("30 30 E"), parse_angle("30 E"), WGS84)
        assert float(got["Y   east"]) == pytest.approx(k.y, abs=1e-3)
        assert abs(float(got["X − Krüger (mm)"])) < 0.01

    async def test_inverse_mode_fields(self, pilot):
        await run(pilot, "tab transverse")
        await pilot.press("h")
        await pilot.pause()
        assert [f.key for f in form(pilot.app).fields] == ["mode", "lon0", "x", "y", "x2", "y2"]
        assert "S   geodesic length" in values(pilot.app)


class TestEllipsoidTab:
    async def test_shows_current_ellipsoid(self, pilot):
        await pilot.press("7")
        await pilot.pause()
        assert text_field(pilot.app, "name").value == "WGS84"
        assert float(values(pilot.app)["b   semi-minor axis"]) == pytest.approx(WGS84.b, abs=1e-4)

    async def test_follows_ellipsoid_changes(self, pilot):
        await pilot.press("7")
        await run(pilot, "ell hayford")
        assert text_field(pilot.app, "name").value == "Hayford"
        assert text_field(pilot.app, "a").value == "6378388.0"

    async def test_editing_updates_derived_values_live(self, pilot):
        await pilot.press("7", "j", "j", "c", *"300", "escape")
        await pilot.pause()
        assert float(values(pilot.app)["f   flattening"]) == pytest.approx(1 / 300)
        assert pilot.app.ellipsoid is WGS84  # not used until saved

    async def test_save_custom_ellipsoid(self, pilot):
        await pilot.press("7", "c", *"Mine", "enter", "c", *"6378000", "enter", "s")
        await pilot.pause()
        app = pilot.app
        assert app.ellipsoid == Ellipsoid("Mine", 6378000.0, 298.257223563)
        assert "Mine" in app.ellipsoids
        assert app.last_message == "using ellipsoid Mine (saved)"

    async def test_custom_ellipsoid_used_by_other_tabs(self, pilot):
        await pilot.press("7", "c", *"Mine", "enter", "c", *"6000000", "enter", "s", "1")
        await pilot.pause()
        expected = Geodetic(parse_angle("41 00 30 N"), parse_angle("29 E"), 150.0).to_cartesian(
            Ellipsoid("Mine", 6000000.0, 298.257223563)
        )
        assert float(values(pilot.app)["X"]) == pytest.approx(expected.x, abs=1e-4)

    async def test_custom_ellipsoid_in_picker_and_completion(self, pilot):
        await pilot.press("7", "c", *"Mine", "enter", "s")
        await run(pilot, "ell wgs")
        assert pilot.app.ellipsoid is WGS84
        await run(pilot, "ell mine")
        assert pilot.app.ellipsoid.name == "Mine"
        await pilot.press("colon", *"ell mi", "tab")
        assert pilot.app.query_one(CommandLine).value == "ell Mine"

    async def test_resaving_replaces_custom(self, pilot):
        await pilot.press("7", "c", *"Mine", "enter", "s", "c", *"6000000", "escape", "s")
        await pilot.pause()
        assert list(pilot.app.ellipsoids).count("Mine") == 1
        assert pilot.app.ellipsoids["Mine"].a == 6000000.0

    async def test_cannot_overwrite_builtin(self, pilot):
        await pilot.press("7", "j", "c", *"6000000", "escape", "s")
        await pilot.pause()
        assert pilot.app.ellipsoid is WGS84
        assert "built-in" in pilot.app.last_message

    async def test_builtin_name_differing_only_in_case_is_rejected(self, pilot):
        await pilot.press("7", "c", *"wgs-84", "enter", "c", *"6000000", "escape", "s")
        await pilot.pause()
        assert pilot.app.ellipsoid is WGS84
        assert "built-in" in pilot.app.last_message

    async def test_invalid_parameters(self, pilot):
        await pilot.press("7", "j", "j", "c", *"0.5", "escape")
        await pilot.pause()
        assert results(pilot.app).has_class("-error")
        assert "inverse flattening" in results(pilot.app).message
        await pilot.press("s")
        assert pilot.app.ellipsoid is WGS84

    async def test_save_command(self, pilot):
        await pilot.press("7", "c", *"Mine", "escape")
        await run(pilot, "w")
        assert pilot.app.ellipsoid.name == "Mine"

    async def test_save_command_elsewhere(self, pilot):
        await run(pilot, "save")
        assert pilot.app.last_message == "nothing to save on this tab"

    async def test_s_does_nothing_on_other_tabs(self, pilot):
        await pilot.press("s")
        assert pilot.app.ellipsoid is WGS84


class TestPersistence:
    async def test_saving_writes_config_file(self, pilot, isolated_config):
        await pilot.press("7", "c", *"Mine", "enter", "c", *"6378000", "escape", "s")
        await pilot.pause()
        text = (isolated_config / "ellipsoids.toml").read_text()
        assert 'name = "Mine"' in text and "a = 6378000.0" in text

    async def test_only_custom_ellipsoids_are_saved(self, pilot, isolated_config):
        await pilot.press("7", "c", *"Mine", "escape", "s")
        await pilot.pause()
        text = (isolated_config / "ellipsoids.toml").read_text()
        assert text.count("[[ellipsoid]]") == 1

    async def test_loaded_at_startup(self, isolated_config):
        save_ellipsoids([Ellipsoid("Mine", 6378000.0, 300.0)], isolated_config / "ellipsoids.toml")
        app = GeolabApp()
        async with app.run_test(size=(130, 40)) as pilot:
            await pilot.pause()
            assert app.ellipsoids["Mine"] == Ellipsoid("Mine", 6378000.0, 300.0)
            assert app.ellipsoid is WGS84  # loading does not change the active ellipsoid

    async def test_clash_with_builtin_is_skipped_with_warning(self, isolated_config):
        save_ellipsoids(
            [Ellipsoid("wgs 84", 6000000.0, 300.0), Ellipsoid("Mine", 6378000.0, 300.0)],
            isolated_config / "ellipsoids.toml",
        )
        app = GeolabApp()
        async with app.run_test(size=(130, 40)) as pilot:
            await pilot.pause()
            assert "wgs 84" not in app.ellipsoids
            assert "Mine" in app.ellipsoids
            assert "clashes with WGS84" in app.last_message

    async def test_bad_file_does_not_stop_startup(self, isolated_config):
        isolated_config.mkdir(parents=True)
        (isolated_config / "ellipsoids.toml").write_text("not [valid")
        app = GeolabApp()
        async with app.run_test(size=(130, 40)) as pilot:
            await pilot.pause()
            assert app.is_running
            assert "could not read" in app.last_message

    async def test_delete_in_picker(self, isolated_config):
        save_ellipsoids([Ellipsoid("Mine", 6378000.0, 300.0)], isolated_config / "ellipsoids.toml")
        app = GeolabApp()
        async with app.run_test(size=(130, 40)) as pilot:
            await run(pilot, "ell mine")
            await pilot.press("e", "d")
            await pilot.pause()
            assert "Mine" not in app.ellipsoids
            assert app.ellipsoid is WGS84
            assert "Mine" not in (isolated_config / "ellipsoids.toml").read_text()
            assert len(app.screen.query_one("#ellipsoid-list").options) == 6

    async def test_builtin_cannot_be_deleted(self, pilot):
        await pilot.press("e", "d")
        await pilot.pause()
        assert "WGS84" in pilot.app.ellipsoids
        assert "can't be deleted" in pilot.app.last_message

    async def test_unwritable_config_keeps_ellipsoid_for_session(self, tmp_path):
        blocker = tmp_path / "file"
        blocker.write_text("")
        app = GeolabApp(config_path=blocker / "ellipsoids.toml")  # parent is a file
        async with app.run_test(size=(130, 40)) as pilot:
            await pilot.press("7", "c", *"Mine", "escape", "s")
            await pilot.pause()
            assert app.ellipsoid.name == "Mine"
            assert "could not save" in app.last_message
