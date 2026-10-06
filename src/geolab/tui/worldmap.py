"""World map tab: plots the Inverse and Direct tabs' geodesics on a braille map."""

import math
from dataclasses import dataclass

from rich.text import Text
from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Vertical
from textual.reactive import reactive
from textual.widget import Widget

from ..angles import AngleKind, format_angle
from ..errors import ConvergenceError
from ..geodesic import direct, geodesic_points, inverse
from ..landmask import is_land
from ..trace import Trace
from .braille import BrailleCanvas
from .calculator import Calculator
from .problems import DirectProblem, InverseProblem
from .widgets import ChoiceField, Field, FieldError, Form, ResultPanel, ResultRow

GRID_SPACINGS = (30.0, 15.0, 10.0, 5.0, 2.0, 1.0, 0.5)
MIN_DEGREES_PER_DOT = 0.1
MIN_FIT_SPAN = 8.0  # degrees: fitting never zooms in further than this
PATH_SAMPLES = 120


@dataclass(frozen=True)
class MapPoint:
    lat: float
    lon: float
    label: str
    """Degrees."""


@dataclass(frozen=True)
class MapPath:
    points: tuple[tuple[float, float], ...]
    """(lat, lon) in degrees."""
    layer: str


class WorldMap(Widget, can_focus=True):
    """Equirectangular braille world map with geodesic paths and points.

    Braille dots are roughly square in a terminal, so one dot spans the same
    number of degrees horizontally and vertically.
    """

    BINDINGS = [
        Binding("h,left", "pan(-1, 0)", "Pan west", show=False),
        Binding("l,right", "pan(1, 0)", "Pan east", show=False),
        Binding("k,up", "pan(0, 1)", "Pan north", show=False),
        Binding("j,down", "pan(0, -1)", "Pan south", show=False),
        Binding("plus,equals_sign", "zoom(2)", "Zoom in", show=False),
        Binding("minus", "zoom(0.5)", "Zoom out", show=False),
        Binding("f", "fit", "Fit to points", show=False),
        Binding("0", "reset", "Whole world", show=False),
    ]

    center_lat: reactive[float] = reactive(0.0)
    center_lon: reactive[float] = reactive(0.0)
    zoom: reactive[float] = reactive(1.0)
    show_grid: reactive[bool] = reactive(True)
    auto_fit: reactive[bool] = reactive(True, init=False)
    """While True the view follows the features; panning or zooming by hand
    turns it off until the features change or ``f`` is pressed."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.paths: list[MapPath] = []
        self.points: list[MapPoint] = []

    # ── View ─────────────────────────────────────────────────────────────────

    @property
    def canvas_size(self) -> tuple[int, int]:
        """Width and height in dots."""
        return self.size.width * 2, self.size.height * 4

    @property
    def degrees_per_dot(self) -> float:
        width, height = self.canvas_size
        if not width or not height:
            return 1.0
        return max(360 / width, 180 / height) / self.zoom

    def max_zoom(self) -> float:
        width, height = self.canvas_size
        if not width or not height:
            return 1.0
        return max(1.0, max(360 / width, 180 / height) / MIN_DEGREES_PER_DOT)

    def set_features(self, paths: list[MapPath], points: list[MapPoint]) -> None:
        """Replace the plotted features. If they changed, zoom to fit them."""
        changed = (paths, points) != (self.paths, self.points)
        self.paths = paths
        self.points = points
        if changed:
            self.auto_fit = True
        if self.auto_fit:
            # The widget may not have its size yet (e.g. its tab was just shown).
            self.call_after_refresh(self._fit_if_auto)
        self.refresh()

    def _fit_if_auto(self) -> None:
        if self.auto_fit:
            self._fit()

    def action_pan(self, dx: int, dy: int) -> None:
        self.auto_fit = False
        width, height = self.canvas_size
        step = self.degrees_per_dot
        self.center_lon = math.remainder(self.center_lon + dx * width * step / 8, 360)
        self.center_lat = max(-90.0, min(90.0, self.center_lat + dy * height * step / 8))

    def action_zoom(self, factor: float) -> None:
        self.auto_fit = False
        self.zoom = max(1.0, min(self.max_zoom(), self.zoom * factor))

    def action_reset(self) -> None:
        self.auto_fit = False
        self.center_lat, self.center_lon, self.zoom = 0.0, 0.0, 1.0

    def action_fit(self) -> None:
        """Fit the view to the features and keep following them."""
        self.auto_fit = True
        self._fit()

    def _fit(self) -> None:
        """Zoom to show every point and path, with a margin."""
        coords = [(p.lat, p.lon) for p in self.points]
        for path in self.paths:
            coords.extend(_unwrap(path.points))
        if not coords:
            self.center_lat, self.center_lon, self.zoom = 0.0, 0.0, 1.0
            return
        # Centre longitudes on the first point so a path across 180° stays contiguous.
        ref = coords[0][1]
        lons = [ref + math.remainder(lon - ref, 360) for _, lon in coords]
        lats = [lat for lat, _ in coords]
        width, height = self.canvas_size
        span = max(
            (max(lons) - min(lons)) / max(width, 1),
            (max(lats) - min(lats)) / max(height, 1),
            MIN_FIT_SPAN / max(width, 1),
        ) * 1.4
        fit = max(360 / max(width, 1), 180 / max(height, 1))
        self.zoom = max(1.0, min(self.max_zoom(), fit / span if span else self.max_zoom()))
        self.center_lon = math.remainder((max(lons) + min(lons)) / 2, 360)
        self.center_lat = (max(lats) + min(lats)) / 2

    def watch_center_lat(self) -> None:
        self._update_subtitle()

    def watch_center_lon(self) -> None:
        self._update_subtitle()

    def watch_zoom(self) -> None:
        self._update_subtitle()

    def watch_auto_fit(self) -> None:
        self._update_subtitle()

    def on_resize(self) -> None:
        if self.auto_fit:
            self._fit()
        else:
            self.zoom = min(self.zoom, self.max_zoom())
        self._update_subtitle()

    def _update_subtitle(self) -> None:
        ns = "N" if self.center_lat >= 0 else "S"
        ew = "E" if self.center_lon >= 0 else "W"
        self.border_subtitle = (
            f"{abs(self.center_lat):.1f}°{ns} {abs(self.center_lon):.1f}°{ew}"
            f" · {self.degrees_per_dot:.2f}°/dot{' · auto' if self.auto_fit else ''}"
            " · hjkl pan  +/- zoom  f fit  0 world"
        )

    # ── Drawing ──────────────────────────────────────────────────────────────

    def to_dot(self, lat: float, lon: float) -> tuple[int, int]:
        width, height = self.canvas_size
        step = self.degrees_per_dot
        x = (lon - self.center_lon) / step + width / 2 - 0.5
        y = (self.center_lat - lat) / step + height / 2 - 0.5
        return round(x), round(y)

    def render(self) -> Text:
        canvas = BrailleCanvas(self.size.width, self.size.height)
        theme = self.app.theme_variables
        canvas.add_layer("grid", "bright_black")
        canvas.add_layer("land", "green")
        canvas.add_layer("direct", f"bold {theme['primary']}")
        canvas.add_layer("inverse", f"bold {theme['accent']}")
        self._draw_background(canvas)
        for path in self.paths:
            self._draw_path(canvas, path)
        for point in self.points:
            self._draw_point(canvas, point, f"bold {theme['warning']}")
        return canvas.render()

    def _draw_background(self, canvas: BrailleCanvas) -> None:
        width, height = canvas.width, canvas.height
        step = self.degrees_per_dot
        lons = [self.center_lon + (x - width / 2 + 0.5) * step for x in range(width)]
        lats = [self.center_lat - (y - height / 2 + 0.5) * step for y in range(height)]
        spacing = next((s for s in GRID_SPACINGS if width * step / s >= 4), GRID_SPACINGS[-1])
        on_meridian = [_crosses(lon, step, spacing) for lon in lons]
        on_parallel = [_crosses(lat, step, spacing) and abs(lat) <= 90 for lat in lats]
        for y, lat in enumerate(lats):
            if abs(lat) > 90:
                continue
            for x, lon in enumerate(lons):
                if is_land(lat, lon):
                    canvas.set("land", x, y)
                elif self.show_grid and (
                    (on_meridian[x] and y % 2 == 0) or (on_parallel[y] and x % 2 == 0)
                ):
                    canvas.set("grid", x, y)

    def _copies(self) -> tuple[float, ...]:
        """Longitude offsets at which features are repeated so they appear in
        every copy of the world visible in the view."""
        width, _ = self.canvas_size
        return (-360.0, 0.0, 360.0) if width * self.degrees_per_dot > 180 else (0.0,)

    def _draw_path(self, canvas: BrailleCanvas, path: MapPath) -> None:
        coords = _unwrap(path.points)
        if not coords:
            return
        shift = math.remainder(coords[0][1] - self.center_lon, 360) - (coords[0][1] - self.center_lon)
        for offset in self._copies():
            dots = [self.to_dot(lat, lon + shift + offset) for lat, lon in coords]
            for (x0, y0), (x1, y1) in zip(dots, dots[1:]):
                canvas.line(path.layer, x0, y0, x1, y1)
            if len(dots) == 1:
                canvas.set(path.layer, *dots[0])

    def _draw_point(self, canvas: BrailleCanvas, point: MapPoint, style: str) -> None:
        lon = self.center_lon + math.remainder(point.lon - self.center_lon, 360)
        for offset in self._copies():
            x, y = self.to_dot(point.lat, lon + offset)
            column, row = x // 2, y // 4
            canvas.put(column, row, "●", style)
            canvas.put(column + 1, row, point.label, style)


def _crosses(value: float, step: float, spacing: float) -> bool:
    """Whether a grid line (a multiple of ``spacing``) falls inside the dot centred on ``value``."""
    return math.floor((value - step / 2) / spacing) != math.floor((value + step / 2) / spacing)


def _unwrap(points: tuple[tuple[float, float], ...]) -> list[tuple[float, float]]:
    """Make longitudes continuous, so a path crossing 180° does not jump across the map."""
    result: list[tuple[float, float]] = []
    for lat, lon in points:
        if result:
            lon = result[-1][1] + math.remainder(lon - result[-1][1], 360)
        result.append((lat, lon))
    return result


class MapTab(Calculator):
    """Plots the geodesics from the Inverse and Direct tabs."""

    TITLE = "Map"

    def compose(self) -> ComposeResult:
        with Vertical(classes="left map-left"):
            yield Form(*self.fields(), classes="block")
            yield ResultPanel(classes="block")
        yield WorldMap(classes="block")

    def fields(self) -> list[Field]:
        return [
            ChoiceField("show", "show", [("Both", "both"), ("Inverse", "inverse"), ("Direct", "direct")]),
            ChoiceField("grid", "grid", [("On", True), ("Off", False)]),
        ]

    def on_mount(self) -> None:
        super().on_mount()
        self.query_one(WorldMap).border_title = "World"

    def activated(self) -> None:
        self.recalculate()

    def calculate(self, trace: Trace) -> list[ResultRow]:
        world = self.query_one(WorldMap)
        world.show_grid = self.choice("grid")
        show = self.choice("show")
        rows: list[ResultRow] = []
        paths: list[MapPath] = []
        points: list[MapPoint] = []
        if show in ("both", "inverse"):
            self._add_inverse(rows, paths, points)
        if show in ("both", "direct"):
            self._add_direct(rows, paths, points)
        world.set_features(paths, points)
        return rows

    def _coordinates_row(self, label: str, lat: float, lon: float) -> ResultRow:
        unit = self.geolab.angle_unit
        text = f"{format_angle(lat, unit, AngleKind.LAT)}  {format_angle(lon, unit, AngleKind.LON)}"
        return ResultRow(label, text, text)

    def _add_inverse(self, rows: list[ResultRow], paths: list[MapPath], points: list[MapPoint]) -> None:
        tab = self.app.query_one(InverseProblem)
        try:
            lat1, lon1, lat2, lon2 = (tab.number(k) for k in ("lat1", "lon1", "lat2", "lon2"))
        except FieldError as exc:
            rows.append(ResultRow("inverse", str(exc), str(exc)))
            return
        rows.append(self._coordinates_row("P₁", lat1, lon1))
        rows.append(self._coordinates_row("P₂", lat2, lon2))
        points += [_point(lat1, lon1, "P₁"), _point(lat2, lon2, "P₂")]
        try:
            r = inverse(lat1, lon1, lat2, lon2, self.geolab.ellipsoid)
        except ConvergenceError:
            rows.append(ResultRow("P₁→P₂", "no path: nearly antipodal", ""))
            return
        rows.append(self.length_row("P₁→P₂", r.distance))
        path = geodesic_points(lat1, lon1, r.azimuth1, r.distance, self.geolab.ellipsoid, PATH_SAMPLES)
        paths.append(_path(path, "inverse"))

    def _add_direct(self, rows: list[ResultRow], paths: list[MapPath], points: list[MapPoint]) -> None:
        tab = self.app.query_one(DirectProblem)
        try:
            lat1, lon1, az1, s = (tab.number(k) for k in ("lat1", "lon1", "az1", "s"))
        except FieldError as exc:
            rows.append(ResultRow("direct", str(exc), str(exc)))
            return
        r = direct(lat1, lon1, az1, s, self.geolab.ellipsoid)
        rows.append(self._coordinates_row("Q₁", lat1, lon1))
        rows.append(self._coordinates_row("Q₂", r.lat2, r.lon2))
        rows.append(self.length_row("Q₁→Q₂", s))
        points += [_point(lat1, lon1, "Q₁"), _point(r.lat2, r.lon2, "Q₂")]
        path = geodesic_points(lat1, lon1, az1, s, self.geolab.ellipsoid, PATH_SAMPLES)
        paths.append(_path(path, "direct"))


def _point(lat: float, lon: float, label: str) -> MapPoint:
    return MapPoint(math.degrees(lat), math.degrees(lon), label)


def _path(points: list[tuple[float, float]], layer: str) -> MapPath:
    return MapPath(tuple((math.degrees(lat), math.degrees(lon)) for lat, lon in points), layer)
