# geolab

A geodesy calculator for the terminal: coordinate conversions, the inverse and
direct geodetic problems, radii of curvature and ellipsoid parameters. It has
a keyboard-driven TUI and a plain Python library.

Every algorithm is implemented from the textbook formulas, with no geodesy
dependencies, and every calculation can **show its working**: each
intermediate value, its formula and, for angles, both radians and DMS.
```
 geolab   1 Geo→ECEF   2 ECEF→Geo   3 Latitudes   4 Inverse   5 Direct   6 Radii   7 Ellipsoid
┌Input · Inverse problem (Vincenty)───────┐┌Working──────────────────────────────────────────┐
│▶ φ₁ latitude    41 00 30 N              ││── Setup                                         │
│  λ₁ longitude   29 00 00 E              ││L        3°51'00.00000"     λ₂ − λ₁              │
│  φ₂ latitude    39 55 00 N              ││         0.067195176202 rad                      │
│  λ₂ longitude   32 51 00 E              ││U₁      40°54'47.07301"     atan((1 − f) tan φ₁) │
└─────────────────────────────────────────┘│...                                              │
┌Result───────────────────────────────────┐│── Iteration 1                                   │
│   s   distance         348 273.5824 m   ││sin σ    0.0545900517058    √((cos U₂ sin λ)² …  │
│   α₁  azimuth at P₁    109°06'35.4977"  ││...                                              │
└─────────────────────────────────────────┘└─────────────────────────────────────────────────┘
 NORMAL   WGS84  a=6 378 137.000  1/f=298.257223563  │  DMS  │  working     i edit  y yank  : cmd  ? help
```

## Features

- **Ellipsoids:** WGS84, GRS80, Hayford (International 1924), Bessel 1841,
  Clarke 1866 and Krassovsky, plus your own custom ellipsoids. For each one:
  b, f, e², e′², n and the mean, authalic and volumetric radii.
- **Geodetic → ECEF** and **ECEF → geodetic**, the reverse by iteration or by
  Bowring's closed-form method.
- **Latitude types:** geodetic φ, geocentric ψ and parametric (reduced) β.
- **Inverse and direct problems:** Vincenty (1975).
- **Radii of curvature:** meridian M, prime vertical N, parallel r, Gaussian
  mean R, and the normal section in any azimuth (Euler).
- **Angles** in DMS, decimal degrees or radians, with N/S/E/W.
- **World map:** a braille-character map in the terminal that plots the
  geodesics from the Inverse and Direct tabs.

## Install and run

Requires Python 3.12+. With [uv](https://docs.astral.sh/uv/):

```sh
uv sync --extra tui
uv run geolab
```

Or with pip: `pip install '.[tui]'`, then `geolab`. The core library has no
dependencies; the TUI needs [Textual](https://textual.textualize.io/).

The TUI uses your terminal's background colour, so transparent terminal
themes show through.

## Using the TUI

geolab works like vim: you move around in **NORMAL** mode, type into a field
in **INSERT** mode, and run commands from the **`:`** line. Results update as
you type. Press `?` inside the app for the full list of keys.

| Key | Action |
|---|---|
| `j` / `k` | move between fields |
| `i`, `a`, `Enter` | edit the selected field (INSERT mode) |
| `c` | clear the field and edit it |
| `Esc` | back to NORMAL mode |
| `h` / `l`, `Space` | cycle an option field (e.g. Iterative / Bowring) |
| `1`–`8`, `gt` / `gT`, `]` / `[` | switch tab |
| `Tab` / `Shift+Tab` | move between the Input, Result and Working blocks |
| `y` / `Y` | copy the selected result / all results to the clipboard |
| `w` | show or hide the working |
| `u` | cycle angle units: DMS → degrees → radians |
| `e` | pick the ellipsoid (`d` deletes a custom one) |
| `?` | help |
| `q` | quit |

In the Result block `j`/`k` choose which result `y` copies. In the Working
block they scroll, and `g`/`G` and `Ctrl+D`/`Ctrl+U` jump.

### Map

Tab 8 draws a world map in braille characters (each character cell holds
2×4 dots). It plots the geodesic from the **Inverse** tab (P₁→P₂, orange)
and the one from the **Direct** tab (Q₁→Q₂, blue), each sampled from the
ellipsoidal solution rather than drawn as a straight line. Use the `show`
option to pick one or both, and `grid` to toggle the lat/lon grid.

The map zooms to fit your points automatically and refits whenever they
change. Press `Tab` to focus it, then `h`/`j`/`k`/`l` to pan and `+`/`-` to
zoom. Moving the view by hand stops the automatic fitting until the points
change; `f` fits again and resumes it, and `0` shows the whole world. The
projection is equirectangular, and paths crossing 180° stay continuous.

The coastlines come from [Natural Earth](https://www.naturalearthdata.com/)
1:110m land polygons (public domain), stored as a 0.5° bitmap in
`src/geolab/data/land.bin`. `scripts/build_land_mask.py` rebuilds it.

### Commands

`Tab` completes command names and their arguments; press it again to cycle
through the matches.

| Command | |
|---|---|
| `:ell <name>` | set the ellipsoid, e.g. `:ell grs80`. Any unambiguous prefix works; with no name it opens the picker |
| `:units dms\|deg\|rad` | set the angle units |
| `:method <name>` | set the method on tabs that have one, e.g. `:method bowring` |
| `:working` | show or hide the working |
| `:save`, `:w` | save a custom ellipsoid (Ellipsoid tab) |
| `:3`, `:tab 3` | switch to tab 3 |
| `:help`, `:q` | help, quit |

### Entering angles

| Input | Meaning |
|---|---|
| `41 00 30.5 N`, `41°00'30.5"N`, `41:00:30.5` | degrees, minutes, seconds |
| `41 30.5` | degrees and decimal minutes |
| `41.508472`, `-29.25` | decimal degrees |
| `S 41 30`, `-0 30 0` | south / west are negative; `-0 30` is −0.5° |

Latitudes must be within ±90° and longitudes within −180°…360°. Invalid input
is highlighted and the reason is shown in the Result block.

### Custom ellipsoids

On the **Ellipsoid** tab, edit the name, `a` and `1/f`. The derived values
update as you type. Press `s` to save the ellipsoid and switch to it. Custom
ellipsoids are stored in `~/.config/geolab/ellipsoids.toml` (or under
`$XDG_CONFIG_HOME`), which you can also edit by hand:

```toml
[[ellipsoid]]
name = "Mine"
a = 6378000.0
inv_f = 298.257223563   # 0 for a sphere
```

## Using the library

All angles are in **radians**; use `parse_angle` and `format_angle` at the
edges.

```python
from geolab import WGS84, Geodetic, Trace, inverse, parse_angle, format_dms

lat = parse_angle("41 00 30 N")
lon = parse_angle("29 E")

xyz = Geodetic(lat, lon, 150.0).to_cartesian(WGS84)
back = xyz.to_geodetic(WGS84)              # iterative; GeodeticMethod.BOWRING also available

result = inverse(lat, lon, parse_angle("39 55 N"), parse_angle("32 51 E"), WGS84)
print(result.distance, format_dms(result.azimuth1))

trace = Trace()                            # pass a Trace to record the working
inverse(lat, lon, parse_angle("39 55 N"), parse_angle("32 51 E"), WGS84, trace)
for step in trace.steps:
    print(step.name, step.value, step.formula)
```

## Accuracy and limitations

- **Vincenty inverse/direct:** matches Vincenty's published test lines to
  under 1 mm and 0.0001″, and GeographicLib to 0.1 mm over 20 000 random lines.
  For **nearly antipodal** points (within about half a degree of exactly
  opposite) the inverse iteration does not converge, and geolab reports this
  instead of returning a wrong answer.
- **ECEF → geodetic:** the iterative method converges to under 1 µm at any
  height. Bowring is sub-micrometre near the surface, about 1 mm at 400 km,
  and a few centimetres at GPS-satellite altitude.
- Height is computed as h = p cos φ + Z sin φ − a²/N, which is equivalent to
  the textbook p / cos φ − N but stays stable near the poles.

## Development

```sh
uv sync --all-extras
uv run pytest                  # everything (~1 min; the TUI tests drive the app with key presses)
uv run pytest --ignore-glob='tests/test_tui*'   # just the library (< 1 s)
uv run textual run --dev src/geolab/tui/app.py:GeolabApp   # with Textual devtools
```
