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
- **Latitude types:** geodetic φ, geocentric ψ and parametric (reduced) β,
  and the point's meridian-ellipse coordinates p, z, geocentric radius r and
  the auxiliary quantities w, v.
- **Inverse and direct problems:** Vincenty (1975).
- **Radii of curvature:** meridian M, prime vertical N, parallel r, Gaussian
  mean R, and the normal section in any azimuth (Euler). Also the first and
  second fundamental quantities E, F, G and L, M̄, N̄ with φ, λ as Gauss
  parameters, the Gauss and mean curvatures, and (in the library) the arc
  and area elements, the azimuth of a curve element and Meusnier's theorem.
- **Arc lengths:** along a meridian (Helmert's series in n) and along a
  parallel (r Δλ). The reverse too: the latitude reached after a given arc
  along a meridian (same longitude), or the longitude after a given arc along
  a parallel (same latitude); set `mode` in the arcs calculator. The arcs
  are also computed with the course's series in e′² (E₀…E₈, F₂…F₆) and the
  short-arc formula, for comparison.
- **Areas:** between two parallels (a zone), between two meridians and two
  parallels (e.g. a 1:25 000 map sheet, the default), and of the whole
  ellipsoid, with the α…ζ series to e¹⁰.
- **Latitude from meridian arc:** given the meridian distance m from the
  equator, the latitudes φ, β, ψ and the coordinates (x, y on the meridian
  ellipse, and X, Y, Z at a chosen longitude). In the 2D ellipse calculator, set
  `mode` to *Point*.
- **Short lines:** the direct and inverse problems by Legendre's series and
  by Gauss's mid-latitude formulas (as in the course), each compared with
  Vincenty. Defaults to the course's Bessel-ellipsoid example.
- **Triangles:** plane and spherical triangles from two angles and the side
  between them, two sides and the angle between them, or three sides. The
  spherical ones use the Gauss sphere R = √(MN) and give the spherical
  excess; for small geodetic triangles, Legendre's theorem and Soldner's
  additament method are shown next to the exact solution.
- **Soldner coordinates:** geodetic rectangular x, y about a central
  meridian, to and from φ, λ with the meridian convergence γ, and the direct
  and inverse problems in x, y.
- **Projection theory:** Tissot's indicatrix from the principal scales m₁,
  m₂ (scale and image in any direction, largest angle distortion ω, area
  scale), the isometric latitude q exactly and by the course's Δq ↔ Δφ
  series about φ₀, and the Mercator projection of the sphere.
- **Geodesic curve:** Clairaut's constant N cos φ sin α, the highest
  latitude a geodesic reaches and its azimuths at any latitude; the azimuth
  and length differences between the geodesic and the normal sections, and
  the target-height reduction of observed directions. The library can also
  integrate the geodesic's differential equations (Runge–Kutta).
- **Datum transformation:** the 7-parameter Helmert transformation between
  TUTGA-99 (GRS80) and ED-50 (Hayford), in either direction, from φ, λ, h.
  The parameters are entered as published (m, arc seconds, ppm); the
  defaults ΔX, ΔY, ΔZ = 87, 98, 121 m are only a rough ED-50 shift, so
  enter the official ones for real work.
- **2D ellipse arcs:** arc length, chord and perimeter on any planar ellipse
  x²/a² + y²/b² = 1, between parametric or polar angles (Carlson's elliptic
  integrals, so any eccentricity works). Set `a, b from` to *Ellipsoid* to
  use the current ellipsoid's meridian ellipse (e.g. Hayford), and the angle
  type to *Normal φ* to enter geodetic latitudes.
- **Angles** in DMS, decimal degrees or radians, with N/S/E/W.
- **World map:** a braille-character map in the terminal that plots the
  geodesics from the Inverse and Direct calculators.

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

There is one calculator on screen at a time. Press `t` (or `/`) to find
another: type a few letters of its name, in order but not necessarily
together (`mpa` finds *Meridian and parallel arcs*), then `Enter`.

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
| `t` or `/` | find a calculator: fuzzy search, like fzf |
| `1`–`9`, `0`, `gt` / `gT`, `]` / `[` | jump to one of the first ten / next / previous calculator |
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

The **Map** calculator draws a world map in braille characters (each character cell holds
2×4 dots). It plots the geodesic from the **Inverse** calculator (P₁→P₂, orange)
and the one from the **Direct** calculator (Q₁→Q₂, blue), each sampled from the
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
| `:method <name>` | set the method on calculators that have one, e.g. `:method bowring` |
| `:working` | show or hide the working |
| `:save`, `:w` | save a custom ellipsoid (Ellipsoid calculator) |
| `:tab <name>`, `:3` | open the best fuzzy match for name, or calculator 3 |
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

In the **Ellipsoid** calculator, edit the name, `a` and `1/f`. The derived values
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

## License

geolab is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version. See [LICENSE](LICENSE).
