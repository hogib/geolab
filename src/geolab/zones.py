"""Gauss–Krüger zones in practice: 3° zones and 6° UTM zones.

To keep the length distortion small the ellipsoid is cut into meridian zones
(dilim), each projected about its own central meridian λ₀:

- 3° zones (large-scale maps, 1/5000 and larger): λ₀ = 3°, 6°, 9°, …, zone
  number λ₀/3; coordinates are used unscaled (m₀ = 1). Points are computed up
  to 2° either side, so neighbouring zones overlap by 1°.
- 6° zones (UTM; 1/25 000 and smaller): zone 1 is 180°W–174°W, numbers run east
  to 60, λ₀ = 6n − 183°, i.e. n = (λ₀ + 3)/6 + 30. Turkey lies in zones 35–38
  (λ₀ = 27°, 33°, 39°, 45°). Coordinates are scaled by m₀ = 0.9996, which
  makes the scale 1 on two lines about 180 km either side of λ₀.

Grid coordinates (sağa değer, yukarı değer) avoid negative values:

    easting  = 500 000 m + m₀ Y      northing = m₀ X

and the zone number may be written in front of the easting (zone × 10⁶ m).
No calculations are done with grid values: convert them back to X, Y first.

Angles are in radians and lengths in metres.
"""

import math
from dataclasses import dataclass

from .ellipsoids import Ellipsoid
from .gauss_kruger import gauss_kruger_from_geodetic, geodetic_from_gauss_kruger
from .trace import Trace

FALSE_EASTING = 500_000.0
UTM_SCALE = 0.9996


def zone_scale(width: int) -> float:
    """m₀: 1 for 3° zones, 0.9996 for 6° (UTM) zones."""
    if width == 3:
        return 1.0
    if width == 6:
        return UTM_SCALE
    raise ValueError("zone width must be 3 or 6 degrees")


def zone_of(lon: float, width: int) -> int:
    """Zone number containing λ: 3° zones number λ₀/3 (1…120), 6° zones 1…60."""
    deg = math.degrees(math.remainder(lon, math.tau))
    if width == 3:
        n = round(deg / 3) % 120
        return n or 120
    if width == 6:
        return min(int((deg + 180) // 6) + 1, 60)
    raise ValueError("zone width must be 3 or 6 degrees")


def central_meridian(zone: int, width: int) -> float:
    """λ₀ of a zone (radians)."""
    if width == 3:
        if not 1 <= zone <= 120:
            raise ValueError("3° zone numbers run from 1 to 120")
        return math.radians(math.remainder(3 * zone, 360))
    if width == 6:
        if not 1 <= zone <= 60:
            raise ValueError("6° zone numbers run from 1 to 60")
        return math.radians(6 * zone - 183)
    raise ValueError("zone width must be 3 or 6 degrees")


@dataclass(frozen=True)
class GridPoint:
    zone: int
    width: int
    central_meridian: float
    x: float
    """X, unscaled Gauss–Krüger north."""
    y: float
    """Y, unscaled Gauss–Krüger east of λ₀."""
    easting: float
    """Sağa değer: 500 000 + m₀Y."""
    northing: float
    """Yukarı değer: m₀X."""
    scale: float
    """m₀ m, the point scale of the grid."""
    convergence: float

    @property
    def prefixed_easting(self) -> float:
        """The easting with the zone number in front: zone × 10⁶ + easting."""
        return self.zone * 1_000_000 + self.easting


def grid_from_geodetic(
    lat: float,
    lon: float,
    ellipsoid: Ellipsoid,
    width: int = 3,
    zone: int | None = None,
    trace: Trace | None = None,
) -> GridPoint:
    """Grid coordinates of φ, λ in its own zone, or in ``zone`` (e.g. a neighbour)."""
    m0 = zone_scale(width)
    if zone is None:
        zone = zone_of(lon, width)
    lon0 = central_meridian(zone, width)
    p = gauss_kruger_from_geodetic(lat, lon, lon0, ellipsoid)
    easting, northing = FALSE_EASTING + m0 * p.y, m0 * p.x
    if trace is not None:
        trace.section("Zone")
        trace.step("zone", zone, f"{width}° zones")
        trace.step("λ₀", lon0, "central meridian", unit="rad")
        trace.step("m₀", m0, "1 for 3° zones, 0.9996 for UTM")
        trace.section("Gauss–Krüger")
        trace.step("X", p.x, unit="m")
        trace.step("Y", p.y, unit="m")
        trace.section("Grid")
        trace.step("sağa", easting, "500 000 + m₀Y", unit="m")
        trace.step("yukarı", northing, "m₀X", unit="m")
        trace.step("m₀m", m0 * p.scale, "grid scale")
    return GridPoint(zone, width, lon0, p.x, p.y, easting, northing, m0 * p.scale, p.convergence)


def split_easting(easting: float, zone: int | None = None) -> tuple[int, float]:
    """Zone and plain easting from an easting that may carry the zone in front."""
    if easting >= 1_000_000:
        prefix = int(easting // 1_000_000)
        if zone is not None and zone != prefix:
            raise ValueError(f"easting says zone {prefix}, not {zone}")
        return prefix, easting - prefix * 1_000_000
    if zone is None:
        raise ValueError("give the zone, or write it in front of the easting")
    return zone, easting


@dataclass(frozen=True)
class GridGeodetic:
    lat: float
    lon: float
    zone: int
    x: float
    y: float
    convergence: float


def geodetic_from_grid(
    easting: float,
    northing: float,
    ellipsoid: Ellipsoid,
    width: int = 3,
    zone: int | None = None,
    trace: Trace | None = None,
) -> GridGeodetic:
    """φ, λ from grid coordinates (sağa, yukarı)."""
    zone, easting = split_easting(easting, zone)
    m0 = zone_scale(width)
    lon0 = central_meridian(zone, width)
    x, y = northing / m0, (easting - FALSE_EASTING) / m0
    g = geodetic_from_gauss_kruger(x, y, lon0, ellipsoid)
    if trace is not None:
        trace.section("Zone")
        trace.step("zone", zone, f"{width}° zones")
        trace.step("λ₀", lon0, "central meridian", unit="rad")
        trace.section("Gauss–Krüger")
        trace.step("Y", y, "(sağa − 500 000)/m₀", unit="m")
        trace.step("X", x, "yukarı/m₀", unit="m")
        trace.section("Geodetic")
        trace.step("φ", g.lat, unit="rad")
        trace.step("λ", g.lon, unit="rad")
    return GridGeodetic(g.lat, g.lon, zone, x, y, g.convergence)


def transform_zone(
    easting: float,
    northing: float,
    ellipsoid: Ellipsoid,
    width: int,
    zone_from: int | None,
    zone_to: int,
    trace: Trace | None = None,
) -> GridPoint:
    """Grid coordinates of a point in a neighbouring zone, through φ, λ."""
    g = geodetic_from_grid(easting, northing, ellipsoid, width, zone_from, trace)
    return grid_from_geodetic(g.lat, g.lon, ellipsoid, width, zone_to, trace)
