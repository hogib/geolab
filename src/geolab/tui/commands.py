"""Parsing for the ``:`` command line. Kept free of Textual so it is easy to test."""

import re
from collections.abc import Sequence

from ..angles import AngleUnit
from ..ellipsoids import ELLIPSOIDS, Ellipsoid


class CommandError(ValueError):
    pass


UNIT_ALIASES = {
    "dms": AngleUnit.DMS,
    "deg": AngleUnit.DEGREES,
    "degrees": AngleUnit.DEGREES,
    "dd": AngleUnit.DEGREES,
    "rad": AngleUnit.RADIANS,
    "radians": AngleUnit.RADIANS,
}


def split_command(text: str) -> tuple[str, list[str]]:
    """Split ``"ell grs80"`` into ``("ell", ["grs80"])``. A leading ``:`` is ignored."""
    parts = text.strip().removeprefix(":").split()
    if not parts:
        return "", []
    return parts[0].lower(), parts[1:]


def _normalise(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", name.lower())


def match_ellipsoid(name: str) -> Ellipsoid:
    """Find an ellipsoid by name, ignoring case, spaces and punctuation.

    An exact match wins; otherwise the name must be a prefix of exactly one ellipsoid.
    """
    wanted = _normalise(name)
    if not wanted:
        raise CommandError("ellipsoid name is empty")
    by_key = {_normalise(n): e for n, e in ELLIPSOIDS.items()}
    if wanted in by_key:
        return by_key[wanted]
    matches = [e for key, e in by_key.items() if key.startswith(wanted)]
    if len(matches) == 1:
        return matches[0]
    if matches:
        names = ", ".join(e.name for e in matches)
        raise CommandError(f"ambiguous ellipsoid {name!r}: {names}")
    raise CommandError(f"unknown ellipsoid {name!r} (try: {', '.join(ELLIPSOIDS)})")


def match_unit(name: str) -> AngleUnit:
    try:
        return UNIT_ALIASES[name.lower()]
    except KeyError:
        raise CommandError(f"unknown unit {name!r} (try: dms, deg, rad)") from None


# Command names offered by tab completion, and whether each takes an argument.
COMMANDS = {
    "ell": True,
    "help": False,
    "method": True,
    "q": False,
    "tab": True,
    "units": True,
    "working": False,
}
UNIT_COMPLETIONS = ["dms", "deg", "rad"]


def completions(text: str, methods: Sequence[str] = (), tab_count: int = 0) -> list[str]:
    """Possible completions of the whole command line ``text``.

    ``methods`` are the option labels of the current tab's method field and
    ``tab_count`` the number of tabs, used to complete ``:method`` and ``:tab``.
    """
    text = text.lstrip().removeprefix(":")
    name, space, arg = text.partition(" ")
    if not space:
        return [
            f"{command} " if takes_arg else command
            for command, takes_arg in COMMANDS.items()
            if command.startswith(name.lower())
        ]

    name = name.lower()
    if name in ("ell", "ellipsoid"):
        wanted = _normalise(arg)
        options = [n for n in ELLIPSOIDS if _normalise(n).startswith(wanted)]
    elif name in ("units", "unit"):
        options = [u for u in UNIT_COMPLETIONS if u.startswith(arg.strip().lower())]
    elif name == "method":
        options = [m.lower() for m in methods if m.lower().startswith(arg.strip().lower())]
    elif name == "tab":
        options = [str(i) for i in range(1, tab_count + 1) if str(i).startswith(arg.strip())]
    else:
        options = []
    return [f"{name} {option}" for option in options]
