"""Persistent user configuration: custom ellipsoids stored as TOML.

The file lives at ``$XDG_CONFIG_HOME/geolab/ellipsoids.toml`` (``~/.config``
by default) and is meant to be readable and editable by hand::

    [[ellipsoid]]
    name = "Mine"
    a = 6378000.0
    inv_f = 298.257223563
"""

import json
import os
import tomllib
from collections.abc import Iterable
from pathlib import Path

from .ellipsoids import Ellipsoid

HEADER = "# Custom ellipsoids saved by geolab. You can edit this file by hand.\n"


def config_dir() -> Path:
    base = os.environ.get("XDG_CONFIG_HOME") or Path.home() / ".config"
    return Path(base) / "geolab"


def ellipsoids_path() -> Path:
    return config_dir() / "ellipsoids.toml"


def load_ellipsoids(path: Path) -> tuple[list[Ellipsoid], list[str]]:
    """Read custom ellipsoids from ``path``.

    Returns the ellipsoids and a list of warnings about entries that were
    skipped. A missing file is not an error.
    """
    try:
        data = tomllib.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        return [], []
    except (OSError, UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        return [], [f"could not read {path}: {exc}"]

    entries = data.get("ellipsoid", [])
    if not isinstance(entries, list):
        return [], [f"{path}: 'ellipsoid' must be an array of tables ([[ellipsoid]])"]

    ellipsoids: list[Ellipsoid] = []
    warnings: list[str] = []
    for i, entry in enumerate(entries, start=1):
        try:
            ellipsoids.append(_parse_entry(entry))
        except (TypeError, ValueError) as exc:
            warnings.append(f"{path}: skipped ellipsoid #{i}: {exc}")
    return ellipsoids, warnings


def _parse_entry(entry: object) -> Ellipsoid:
    if not isinstance(entry, dict):
        raise TypeError("expected a table")
    name, a, inv_f = entry.get("name"), entry.get("a"), entry.get("inv_f")
    if not isinstance(name, str) or not name.strip():
        raise ValueError("missing name")
    for key, value in (("a", a), ("inv_f", inv_f)):
        if isinstance(value, bool) or not isinstance(value, int | float):
            raise ValueError(f"{key} must be a number")
    return Ellipsoid(name.strip(), float(a), float(inv_f))


def save_ellipsoids(ellipsoids: Iterable[Ellipsoid], path: Path) -> None:
    """Write ``ellipsoids`` to ``path``, replacing it atomically."""
    blocks = [
        f"[[ellipsoid]]\nname = {json.dumps(e.name, ensure_ascii=False)}\n"
        f"a = {e.a!r}\ninv_f = {e.inv_f!r}\n"
        for e in ellipsoids
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(HEADER + "\n" + "\n".join(blocks), encoding="utf-8")
    os.replace(tmp, path)
