"""Parsing and formatting of angles. Everything inside geolab works in radians."""

import math
import re
from enum import Enum, auto


class AngleKind(Enum):
    LAT = auto()
    LON = auto()


class AngleUnit(Enum):
    DMS = auto()
    DEGREES = auto()
    RADIANS = auto()


class AngleParseError(ValueError):
    pass


_HEMISPHERES = {
    "N": (AngleKind.LAT, 1),
    "S": (AngleKind.LAT, -1),
    "E": (AngleKind.LON, 1),
    "W": (AngleKind.LON, -1),
}
_HEMISPHERE_LETTERS = {
    AngleKind.LAT: ("N", "S"),
    AngleKind.LON: ("E", "W"),
}

# Position of each symbol in a degrees/minutes/seconds sequence.
_SYMBOL_POSITION = {"°": 0, "'": 1, "′": 1, "’": 1, '"': 2, "″": 2, "''": 2, "′′": 2, "’’": 2}
_TOKEN = re.compile(r"\s*(\d+(?:\.\d*)?|\.\d+)\s*(''|′′|’’|[°'′’\"″])?")


def parse_angle(
    text: str,
    kind: AngleKind | None = None,
    unit: AngleUnit = AngleUnit.DEGREES,
) -> float:
    """Parse an angle and return it in radians.

    Accepts decimal numbers ("41.5", "-29.25"), DMS with symbols (41°00'30.5"N)
    or separated by spaces/colons ("41 0 30.5", "41:00:30.5"), and degrees plus
    decimal minutes ("41 30.5"). A hemisphere letter may lead or trail. A bare
    single number is read in ``unit`` (radians when ``unit`` is RADIANS,
    otherwise decimal degrees).
    """
    s = text.strip()
    if not s:
        raise AngleParseError("empty angle")

    hemisphere = None
    if s[0].upper() in _HEMISPHERES:
        hemisphere, s = s[0].upper(), s[1:].strip()
    elif s[-1].upper() in _HEMISPHERES:
        hemisphere, s = s[-1].upper(), s[:-1].strip()

    negative = s.startswith("-")
    if s[:1] in "+-":
        s = s[1:].lstrip()
    if negative and hemisphere:
        raise AngleParseError("use either a sign or a hemisphere letter, not both")

    sign = -1 if negative else 1
    if hemisphere:
        hemi_kind, sign = _HEMISPHERES[hemisphere]
        if kind is not None and hemi_kind != kind:
            raise AngleParseError(f"{hemisphere} is not a valid hemisphere for {kind.name.lower()}")
        kind = hemi_kind

    components = _tokenize(s.replace(":", " "))

    if len(components) == 1 and components[0][1] is None and unit == AngleUnit.RADIANS:
        radians = sign * float(components[0][0])
        _check_range(math.degrees(radians), kind, text)
        return radians

    for number, _ in components[:-1]:
        if "." in number:
            raise AngleParseError("only the last component may have a decimal part")
    values = [float(number) for number, _ in components]
    for name, value in zip(("minutes", "seconds"), values[1:]):
        if value >= 60:
            raise AngleParseError(f"{name} must be less than 60")

    degrees = sign * sum(v / 60**i for i, v in enumerate(values))
    _check_range(degrees, kind, text)
    return math.radians(degrees)


def _tokenize(s: str) -> list[tuple[str, str | None]]:
    components: list[tuple[str, str | None]] = []
    pos = 0
    while pos < len(s):
        match = _TOKEN.match(s, pos)
        if match is None:
            if s[pos:].strip():
                raise AngleParseError(f"unexpected text: {s[pos:].strip()!r}")
            break
        number, symbol = match.groups()
        if symbol is not None and _SYMBOL_POSITION[symbol] != len(components):
            raise AngleParseError(f"unexpected {symbol!r} after {len(components)} component(s)")
        components.append((number, symbol))
        pos = match.end()

    if not components:
        raise AngleParseError("no number found")
    if len(components) > 3:
        raise AngleParseError("too many components (expected degrees, minutes, seconds)")
    return components


def _check_range(degrees: float, kind: AngleKind | None, text: str) -> None:
    if kind == AngleKind.LAT and abs(degrees) > 90:
        raise AngleParseError(f"latitude out of range [-90°, 90°]: {text!r}")
    if kind == AngleKind.LON and not -180 <= degrees <= 360:
        raise AngleParseError(f"longitude out of range [-180°, 360°]: {text!r}")


def to_dms(radians: float, sec_decimals: int = 4) -> tuple[int, int, int, float]:
    """Split an angle into (sign, degrees, minutes, seconds).

    Seconds are rounded to ``sec_decimals`` and carried into minutes and degrees,
    so the result never shows 60 seconds or 60 minutes.
    """
    scale = 10**sec_decimals
    units = round(abs(math.degrees(radians)) * 3600 * scale)
    degrees, rest = divmod(units, 3600 * scale)
    minutes, rest = divmod(rest, 60 * scale)
    sign = -1 if radians < 0 and units != 0 else 1
    return sign, degrees, minutes, rest / scale


def format_dms(radians: float, kind: AngleKind | None = None, sec_decimals: int = 4) -> str:
    sign, d, m, s = to_dms(radians, sec_decimals)
    width = 2 + (sec_decimals + 1 if sec_decimals > 0 else 0)
    body = f"{d}°{m:02d}'{s:0{width}.{sec_decimals}f}\""
    if kind is None:
        return f"-{body}" if sign < 0 else body
    positive, negative = _HEMISPHERE_LETTERS[kind]
    return body + (negative if sign < 0 else positive)


def format_degrees(radians: float, decimals: int = 8) -> str:
    return f"{_no_negative_zero(math.degrees(radians), decimals):.{decimals}f}°"


def format_radians(radians: float, decimals: int = 10) -> str:
    return f"{_no_negative_zero(radians, decimals):.{decimals}f} rad"


def format_angle(radians: float, unit: AngleUnit, kind: AngleKind | None = None) -> str:
    if unit == AngleUnit.DMS:
        return format_dms(radians, kind)
    if unit == AngleUnit.DEGREES:
        return format_degrees(radians)
    return format_radians(radians)


def _no_negative_zero(value: float, decimals: int) -> float:
    rounded = round(value, decimals)
    return 0.0 if rounded == 0 else rounded
