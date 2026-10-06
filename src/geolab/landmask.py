"""Coarse land/sea lookup for drawing maps.

Uses a 0.5° bitmap built from Natural Earth 1:110m land polygons (public
domain) by ``scripts/build_land_mask.py``. Good enough for an overview map,
not for deciding whether a particular point is on land.
"""

from functools import cache
from importlib.resources import files

RESOLUTION = 0.5
COLUMNS = int(360 / RESOLUTION)
ROWS = int(180 / RESOLUTION)


@cache
def _bits() -> bytes:
    data = files("geolab").joinpath("data/land.bin").read_bytes()
    if len(data) != COLUMNS * ROWS // 8:
        raise ValueError(f"land.bin has {len(data)} bytes, expected {COLUMNS * ROWS // 8}")
    return data


def is_land(lat: float, lon: float) -> bool:
    """Whether the 0.5° cell containing (``lat``, ``lon``), in degrees, is mostly land.

    Longitudes wrap around; latitudes outside ±90° are never land.
    """
    if not -90 <= lat <= 90:
        return False
    row = min(int((90 - lat) / RESOLUTION), ROWS - 1)
    col = int(((lon + 180) % 360) / RESOLUTION) % COLUMNS
    index = row * COLUMNS + col
    return bool(_bits()[index >> 3] & (0x80 >> (index & 7)))
