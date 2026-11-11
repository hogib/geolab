"""Turkish map-sheet division (pafta bölümlendirmesi).

The 1/100 000 sheets are 30′ × 30′ cells named by a row letter and a column
number: rows D, E, F, …, P, R (no Q) run south from 42°30′ N, columns 15…53 run
east from 25°30′ E, so G25 is 40°30′–41° N, 30°30′–31° E. Each sheet is then
divided further:

    1/100 000  G25                30′ × 30′
    1/50 000   G25-d              quarters  a b / d c        15′ × 15′
    1/10 000   G25-d-01           5 × 5,   01…25 by rows      3′ × 3′
    1/5000     G25-d-01-d         quarters  a b / d c        1′30″ × 1′30″
    1/2000     G25-d-01-d-1       quarters  1 2 / 4 3        45″ × 45″
    1/1000     G25-d-01-d-1-a     quarters  a b / d c        22.5″ × 22.5″
    1/500      G25-d-01-d-1-a-2   quarters  1 2 / 4 3        11.25″ × 11.25″

(quarters go clockwise from the north-west: a/1 NW, b/2 NE, c/3 SE, d/4 SW).
The course's corner lists label the corners A (SW), B (NW), C (NE), D (SE).
Angles are in radians.
"""

import math
from dataclasses import dataclass

from .ellipsoids import Ellipsoid
from .zones import GridPoint, grid_from_geodetic, zone_of

ROWS = "DEFGHIJKLMNOPR"
NORTH_EDGE = 42.5
"""Northern edge of row D (degrees)."""
FIRST_COLUMN, LAST_COLUMN = 15, 53
SCALES = (100_000, 50_000, 10_000, 5000, 2000, 1000, 500)

_LETTERS = "abdc"  # NW, NE, SW, SE in reading order: a b / d c
_DIGITS = "1243"  # 1 2 / 4 3


@dataclass(frozen=True)
class Sheet:
    name: str
    scale: int
    south: float
    north: float
    west: float
    east: float

    def corners(self) -> dict[str, tuple[float, float]]:
        """φ, λ of the corners A (SW), B (NW), C (NE), D (SE)."""
        return {
            "A": (self.south, self.west),
            "B": (self.north, self.west),
            "C": (self.north, self.east),
            "D": (self.south, self.east),
        }


def _subdivisions(level: int) -> tuple[int, str]:
    """(n, labels) for the step from SCALES[level - 1] to SCALES[level]: an n × n grid."""
    if level == 2:
        return 5, ""
    return 2, _DIGITS if level in (4, 6) else _LETTERS


def _label(level: int, row: int, col: int) -> str:
    n, labels = _subdivisions(level)
    if not labels:
        return f"{row * n + col + 1:02d}"
    return labels[row * 2 + col]


def sheet_at(lat: float, lon: float, scale: int) -> Sheet:
    """The sheet at ``scale`` that contains φ, λ."""
    if scale not in SCALES:
        raise ValueError(f"scale must be one of {', '.join(f'1/{s}' for s in SCALES)}")
    lat_d, lon_d = math.degrees(lat), math.degrees(lon)
    row = math.floor((NORTH_EDGE - lat_d) / 0.5)
    col = math.floor((lon_d - 18) / 0.5)
    if not (0 <= row < len(ROWS) and FIRST_COLUMN <= col <= LAST_COLUMN):
        raise ValueError("outside the Turkish sheet grid (35°30′–42°30′ N, 25°30′–45° E)")
    north, west, size = NORTH_EDGE - 0.5 * row, 18 + 0.5 * col, 0.5
    name = f"{ROWS[row]}{col}"
    for level in range(1, SCALES.index(scale) + 1):
        n, _ = _subdivisions(level)
        size /= n
        r = min(int((north - lat_d) / size), n - 1)
        c = min(int((lon_d - west) / size), n - 1)
        name += "-" + _label(level, r, c)
        north -= r * size
        west += c * size
    return Sheet(name, scale, math.radians(north - size), math.radians(north), math.radians(west), math.radians(west + size))


def sheet_from_name(name: str) -> Sheet:
    """The sheet called e.g. "G25-d-01-d-1-a-2" (or "G-25-…")."""
    parts = [p for p in name.strip().replace(" ", "").split("-") if p]
    if parts and len(parts) > 1 and parts[0].isalpha() and parts[1].isdigit() and len(parts[0]) == 1:
        parts = [parts[0] + parts[1]] + parts[2:]  # "G-25" → "G25"
    if not parts or len(parts[0]) < 2:
        raise ValueError(f"{name!r} is not a sheet name like G25-d-01-d")
    letter, number = parts[0][0].upper(), parts[0][1:]
    if letter not in ROWS or not number.isdigit() or not FIRST_COLUMN <= int(number) <= LAST_COLUMN:
        raise ValueError(f"{parts[0]!r} is not a 1/100 000 sheet (rows D–R, columns 15–53)")
    if len(parts) > len(SCALES):
        raise ValueError(f"{name!r} has too many parts")
    north, west, size = NORTH_EDGE - 0.5 * ROWS.index(letter), 18 + 0.5 * int(number), 0.5
    canonical = f"{letter}{int(number)}"
    for level, part in enumerate(parts[1:], 1):
        n, labels = _subdivisions(level)
        size /= n
        if labels:
            if len(part) != 1 or part.lower() not in labels:
                raise ValueError(f"part {level + 1} of {name!r} must be one of {', '.join(sorted(labels))}")
            index = labels.index(part.lower())
        else:
            if not part.isdigit() or not 1 <= int(part) <= 25:
                raise ValueError(f"part {level + 1} of {name!r} must be 01…25")
            index = int(part) - 1
        r, c = divmod(index, n)
        north -= r * size
        west += c * size
        canonical += "-" + _label(level, r, c)
    scale = SCALES[len(parts) - 1]
    return Sheet(canonical, scale, math.radians(north - size), math.radians(north), math.radians(west), math.radians(west + size))


def sheet_corner_coordinates(
    sheet: Sheet, ellipsoid: Ellipsoid, width: int = 3, zone: int | None = None
) -> dict[str, GridPoint]:
    """Grid coordinates of the sheet's corners A, B, C, D, all in one zone
    (by default the zone of the sheet's centre)."""
    if zone is None:
        zone = zone_of((sheet.west + sheet.east) / 2, width)
    return {k: grid_from_geodetic(lat, lon, ellipsoid, width, zone) for k, (lat, lon) in sheet.corners().items()}
