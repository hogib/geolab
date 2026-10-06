import pytest

pytest.importorskip("textual")

from geolab.tui.app import GeolabApp  # noqa: E402
from geolab.tui.braille import BrailleCanvas  # noqa: E402
from geolab.tui.widgets import ResultPanel, TextField  # noqa: E402
from geolab.tui.worldmap import MapPath, WorldMap, _unwrap  # noqa: E402


class TestBrailleCanvas:
    def test_dimensions(self):
        canvas = BrailleCanvas(3, 2)
        assert (canvas.width, canvas.height) == (6, 8)

    def test_single_dots(self):
        canvas = BrailleCanvas(1, 1)
        canvas.add_layer("a", "")
        canvas.set("a", 0, 0)
        assert canvas.render().plain == "⠁"
        canvas.set("a", 1, 3)
        assert canvas.render().plain == "⢁"

    def test_full_cell(self):
        canvas = BrailleCanvas(1, 1)
        canvas.add_layer("a", "")
        for x in range(2):
            for y in range(4):
                canvas.set("a", x, y)
        assert canvas.render().plain == "⣿"

    def test_out_of_range_dots_ignored(self):
        canvas = BrailleCanvas(2, 1)
        canvas.add_layer("a", "")
        canvas.set("a", -1, 0)
        canvas.set("a", 4, 0)
        canvas.set("a", 0, 4)
        assert canvas.render().plain == "  "

    def test_rows_joined_with_newlines(self):
        canvas = BrailleCanvas(2, 2)
        canvas.add_layer("a", "")
        assert canvas.render().plain == "  \n  "

    def test_line(self):
        canvas = BrailleCanvas(4, 1)
        canvas.add_layer("a", "")
        canvas.line("a", 0, 0, 7, 3)
        assert all(canvas.is_set("a", x, y) for x, y in [(0, 0), (7, 3)])
        assert sum(canvas.is_set("a", x, y) for x in range(8) for y in range(4)) == 8

    def test_later_layer_wins(self):
        canvas = BrailleCanvas(1, 1)
        canvas.add_layer("below", "green")
        canvas.add_layer("above", "red")
        canvas.set("below", 0, 0)
        canvas.set("above", 1, 0)
        text = canvas.render()
        assert text.plain == "⠈"
        assert str(text.spans[0].style) == "red" if text.spans else True

    def test_text_wins_over_dots(self):
        canvas = BrailleCanvas(3, 1)
        canvas.add_layer("a", "")
        for x in range(6):
            canvas.set("a", x, 0)
        canvas.put(1, 0, "P₁")
        assert canvas.render().plain == "⠉P₁"

    def test_text_clipped_at_edge(self):
        canvas = BrailleCanvas(2, 1)
        canvas.put(1, 0, "abc")
        assert canvas.render().plain == " a"


class TestUnwrap:
    def test_crossing_antimeridian_stays_continuous(self):
        assert _unwrap(((0.0, 179.0), (0.0, -179.0), (0.0, -178.0))) == [
            (0.0, 179.0),
            (0.0, 181.0),
            (0.0, 182.0),
        ]

    def test_empty(self):
        assert _unwrap(()) == []


@pytest.fixture
async def pilot():
    app = GeolabApp()
    async with app.run_test(size=(140, 40)) as pilot:
        await pilot.press("8")
        await pilot.pause()
        yield pilot


def world(app: GeolabApp) -> WorldMap:
    return app.query_one(WorldMap)


def rows(app: GeolabApp) -> dict[str, str]:
    return {row.label: row.display for row in app.active_calculator.query_one(ResultPanel).rows}


class TestWorldMap:
    async def test_whole_world_fits(self, pilot):
        wm = world(pilot.app)
        width, height = wm.canvas_size
        assert wm.degrees_per_dot == pytest.approx(max(360 / width, 180 / height))
        x_west, _ = wm.to_dot(0, -180)
        x_east, _ = wm.to_dot(0, 180)
        _, y_north = wm.to_dot(90, 0)
        _, y_south = wm.to_dot(-90, 0)
        assert x_west >= -1 and x_east <= width
        assert y_north >= -1 and y_south <= height

    async def test_centre_maps_to_middle(self, pilot):
        wm = world(pilot.app)
        width, height = wm.canvas_size
        x, y = wm.to_dot(0, 0)
        assert abs(x - width / 2) <= 1 and abs(y - height / 2) <= 1

    async def test_render_size_matches_widget(self, pilot):
        wm = world(pilot.app)
        lines = wm.render().plain.split("\n")
        assert len(lines) == wm.size.height
        assert all(len(line) == wm.size.width for line in lines)

    async def test_land_is_drawn(self, pilot):
        text = world(pilot.app).render().plain
        assert sum(0x2800 < ord(c) <= 0x28FF for c in text) > 200

    async def test_points_are_labelled(self, pilot):
        text = world(pilot.app).render().plain
        assert "Q₂" in text or "P₂" in text

    async def test_zoom_and_limits(self, pilot):
        await pilot.press("tab", "tab")
        wm = world(pilot.app)
        before = wm.degrees_per_dot
        await pilot.press("plus")
        assert wm.degrees_per_dot == pytest.approx(before / 2)
        await pilot.press("minus", "minus")
        assert wm.zoom == 1.0  # cannot zoom out past the whole world
        for _ in range(20):
            await pilot.press("plus")
        assert wm.degrees_per_dot >= 0.1 - 1e-9

    async def test_pan(self, pilot):
        await pilot.press("tab", "tab", "l", "k")
        wm = world(pilot.app)
        assert wm.center_lon > 0
        assert wm.center_lat > 0

    async def test_pan_wraps_longitude_and_clamps_latitude(self, pilot):
        await pilot.press("tab", "tab")
        wm = world(pilot.app)
        for _ in range(12):
            await pilot.press("l")
        assert -180 <= wm.center_lon <= 180
        for _ in range(20):
            await pilot.press("k")
        assert wm.center_lat == 90.0

    async def test_fit_zooms_to_points(self, pilot):
        await pilot.press("tab", "tab", "f")
        wm = world(pilot.app)
        assert wm.zoom > 1
        assert 39 < wm.center_lat < 42
        assert 28 < wm.center_lon < 34
        await pilot.press("0")
        assert (wm.center_lat, wm.center_lon, wm.zoom) == (0.0, 0.0, 1.0)

    async def test_fit_across_antimeridian(self, pilot):
        wm = world(pilot.app)
        wm.set_features([MapPath(((0.0, 170.0), (0.0, -170.0)), "inverse")], [])
        wm.action_fit()
        assert abs(wm.center_lon) > 170  # centred near 180°, not 0°

    async def test_subtitle_shows_view(self, pilot):
        assert "0.0°N 0.0°E" in world(pilot.app).border_subtitle


class TestMapTab:
    async def test_rows_for_both_problems(self, pilot):
        got = rows(pilot.app)
        assert got["P₁"] == "41°00'30.0000\"N  29°00'00.0000\"E"
        assert got["P₁→P₂"] == "348 273.5824 m"
        assert got["Q₂"] == "39°55'00.0000\"N  32°51'00.0000\"E"

    async def test_two_paths(self, pilot):
        assert sorted(p.layer for p in world(pilot.app).paths) == ["direct", "inverse"]

    async def test_show_option(self, pilot):
        await pilot.press("l")  # Both → Inverse
        await pilot.pause()
        assert [p.layer for p in world(pilot.app).paths] == ["inverse"]
        assert "Q₁" not in rows(pilot.app)
        await pilot.press("l")  # → Direct
        await pilot.pause()
        assert [p.layer for p in world(pilot.app).paths] == ["direct"]

    async def test_grid_option(self, pilot):
        await pilot.press("j", "l")
        await pilot.pause()
        assert world(pilot.app).show_grid is False

    async def test_follows_inverse_tab(self, pilot):
        await pilot.press("4", "j", "j", "c", *"33 52 S", "enter", "c", *"151 13 E", "escape", "8")
        await pilot.pause()
        assert rows(pilot.app)["P₂"].startswith("33°52'00.0000\"S")
        path = next(p for p in world(pilot.app).paths if p.layer == "inverse")
        assert path.points[-1] == pytest.approx((-(33 + 52 / 60), 151 + 13 / 60), abs=1e-6)

    async def test_invalid_inverse_input(self, pilot):
        await pilot.press("4", "c", "escape", "8")
        await pilot.pause()
        assert rows(pilot.app)["inverse"] == "Enter φ₁ latitude"
        assert [p.layer for p in world(pilot.app).paths] == ["direct"]

    async def test_antipodal_inverse_has_points_but_no_path(self, pilot):
        app = pilot.app
        await pilot.press("4")
        for key, value in [("lat1", "0"), ("lon1", "0"), ("lat2", "0.5"), ("lon2", "179.7")]:
            field = app.active_calculator.field(key)
            assert isinstance(field, TextField)
            field.value = value
        await pilot.press("8")
        await pilot.pause()
        assert rows(app)["P₁→P₂"] == "no path: nearly antipodal"
        assert [p.layer for p in world(app).paths] == ["direct"]

    async def test_units_apply(self, pilot):
        await pilot.press("u")
        await pilot.pause()
        assert rows(pilot.app)["P₁"] == "41.00833333°  29.00000000°"
