"""Build geolab's land mask from Natural Earth land polygons.

Usage:
    python scripts/build_land_mask.py ne_110m_land.geojson src/geolab/data/land.bin

Download ne_110m_land.geojson from
https://github.com/nvkelso/natural-earth-vector/tree/master/geojson
(Natural Earth data is in the public domain.)

The output is a 720 × 360 bitmap at 0.5° resolution, one bit per cell, packed
row by row (most significant bit first). Row 0 is the strip from 90°N to
89.5°N; column 0 is from 180°W to 179.5°W. A bit is set when the cell's centre
is on land.
"""

import json
import sys
from pathlib import Path

RESOLUTION = 0.5
COLUMNS = int(360 / RESOLUTION)
ROWS = int(180 / RESOLUTION)


def rings(geojson: dict) -> list[list[tuple[float, float]]]:
    result = []
    for feature in geojson["features"]:
        geometry = feature["geometry"]
        polygons = geometry["coordinates"]
        if geometry["type"] == "Polygon":
            polygons = [polygons]
        for polygon in polygons:
            result.extend([(x, y) for x, y, *_ in ring] for ring in polygon)
    return result


def rasterize(all_rings: list[list[tuple[float, float]]]) -> bytearray:
    """Even-odd scanline fill. Land polygons do not overlap, so holes
    (lakes, inland seas) work without tracking which ring is which."""
    edges = []
    for ring in all_rings:
        for (x0, y0), (x1, y1) in zip(ring, ring[1:] + ring[:1]):
            if y0 != y1:
                edges.append((x0, y0, x1, y1))

    bits = bytearray(COLUMNS * ROWS // 8)
    for row in range(ROWS):
        lat = 90 - (row + 0.5) * RESOLUTION
        crossings = sorted(
            x0 + (lat - y0) * (x1 - x0) / (y1 - y0)
            for x0, y0, x1, y1 in edges
            if (y0 > lat) != (y1 > lat)
        )
        for start, end in zip(crossings[::2], crossings[1::2]):
            first = max(0, round((start + 180) / RESOLUTION - 0.5))
            last = min(COLUMNS - 1, round((end + 180) / RESOLUTION - 0.5))
            for col in range(first, last + 1):
                lon = -180 + (col + 0.5) * RESOLUTION
                if start <= lon <= end:
                    index = row * COLUMNS + col
                    bits[index // 8] |= 0x80 >> (index % 8)
    return bits


def main() -> None:
    source, target = map(Path, sys.argv[1:3])
    bits = rasterize(rings(json.loads(source.read_text())))
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(bytes(bits))
    land = sum(bin(b).count("1") for b in bits)
    print(f"wrote {target}: {len(bits)} bytes, {land / (COLUMNS * ROWS):.1%} of cells are land")


if __name__ == "__main__":
    main()
