"""A canvas of braille characters: each terminal cell holds 2 × 4 dots."""

from rich.style import Style
from rich.text import Text

# Bit for the dot at (column dx, row dy) within a braille cell.
_DOT_BITS = ((0x01, 0x08), (0x02, 0x10), (0x04, 0x20), (0x40, 0x80))
_BRAILLE_BASE = 0x2800


class BrailleCanvas:
    """Dots drawn on named layers. Where layers overlap in a cell, the one
    added last wins; text placed with ``put`` wins over every layer."""

    def __init__(self, columns: int, rows: int):
        self.columns = max(columns, 0)
        self.rows = max(rows, 0)
        self._layers: dict[str, tuple[bytearray, Style]] = {}
        self._text: dict[tuple[int, int], tuple[str, Style]] = {}

    @property
    def width(self) -> int:
        """Width in dots."""
        return self.columns * 2

    @property
    def height(self) -> int:
        """Height in dots."""
        return self.rows * 4

    def add_layer(self, name: str, style: Style | str) -> None:
        style = Style.parse(style) if isinstance(style, str) else style
        self._layers[name] = (bytearray(self.columns * self.rows), style)

    def set(self, layer: str, x: int, y: int) -> None:
        """Set the dot at (``x``, ``y``) in dots; dots outside the canvas are ignored."""
        if 0 <= x < self.width and 0 <= y < self.height:
            cells, _ = self._layers[layer]
            cells[(y >> 2) * self.columns + (x >> 1)] |= _DOT_BITS[y & 3][x & 1]

    def line(self, layer: str, x0: int, y0: int, x1: int, y1: int) -> None:
        """Bresenham line between two dots (inclusive)."""
        dx, dy = abs(x1 - x0), -abs(y1 - y0)
        sx, sy = (1 if x0 < x1 else -1), (1 if y0 < y1 else -1)
        err = dx + dy
        while True:
            self.set(layer, x0, y0)
            if x0 == x1 and y0 == y1:
                return
            e2 = 2 * err
            if e2 >= dy:
                err += dy
                x0 += sx
            if e2 <= dx:
                err += dx
                y0 += sy

    def put(self, column: int, row: int, text: str, style: Style | str = "") -> None:
        """Write ``text`` starting at a character cell, over any dots."""
        style = Style.parse(style) if isinstance(style, str) else style
        for i, char in enumerate(text):
            if 0 <= column + i < self.columns and 0 <= row < self.rows:
                self._text[(column + i, row)] = (char, style)

    def is_set(self, layer: str, x: int, y: int) -> bool:
        cells, _ = self._layers[layer]
        return bool(cells[(y >> 2) * self.columns + (x >> 1)] & _DOT_BITS[y & 3][x & 1])

    def render(self) -> Text:
        layers = list(self._layers.values())[::-1]
        text = Text(no_wrap=True, overflow="crop")
        for row in range(self.rows):
            if row:
                text.append("\n")
            for column in range(self.columns):
                placed = self._text.get((column, row))
                if placed is not None:
                    text.append(*placed)
                    continue
                index = row * self.columns + column
                for cells, style in layers:
                    if cells[index]:
                        text.append(chr(_BRAILLE_BASE | cells[index]), style)
                        break
                else:
                    text.append(" ")
        return text
