"""Datum transformation with the 7-parameter (Helmert) similarity transformation.

    [X]          [ΔX]              [  1   −R_Z   R_Y ] [X]
    [Y]        = [ΔY] + (1 + ΔS) · [ R_Z    1   −R_X ] [Y]
    [Z]_target   [ΔZ]              [−R_Y   R_X    1  ] [Z]_source

ΔX, ΔY, ΔZ are translations (m), R_X, R_Y, R_Z small rotations (radians) and
ΔS the scale difference. This is the form used for e.g. TUTGA-99 (ITRF96,
GRS80) → ED-50 (Hayford). Rotation angles are assumed small (a few arc
seconds), so the linearised rotation matrix is used, as in the formula above.
"""

import math
from dataclasses import dataclass

from .ellipsoids import Ellipsoid
from .point import Cartesian3D, Geodetic
from .trace import Trace

ARCSEC = math.pi / (180 * 3600)
PPM = 1e-6


@dataclass(frozen=True)
class HelmertParams:
    tx: float = 0.0
    ty: float = 0.0
    tz: float = 0.0
    """Translations ΔX, ΔY, ΔZ (m)."""
    rx: float = 0.0
    ry: float = 0.0
    rz: float = 0.0
    """Rotations R_X, R_Y, R_Z (radians)."""
    scale: float = 0.0
    """Scale difference ΔS (dimensionless; 1 ppm = 1e-6)."""

    @classmethod
    def from_units(
        cls,
        tx: float = 0.0,
        ty: float = 0.0,
        tz: float = 0.0,
        rx_arcsec: float = 0.0,
        ry_arcsec: float = 0.0,
        rz_arcsec: float = 0.0,
        scale_ppm: float = 0.0,
    ) -> "HelmertParams":
        """Parameters as they are usually published: metres, arc seconds and ppm."""
        return cls(
            tx, ty, tz, rx_arcsec * ARCSEC, ry_arcsec * ARCSEC, rz_arcsec * ARCSEC, scale_ppm * PPM
        )

    def matrix(self) -> tuple[tuple[float, float, float], ...]:
        """(1 + ΔS) R, the linear part of the transformation."""
        k = 1 + self.scale
        return (
            (k, -k * self.rz, k * self.ry),
            (k * self.rz, k, -k * self.rx),
            (-k * self.ry, k * self.rx, k),
        )


def helmert(point: Cartesian3D, p: HelmertParams, trace: Trace | None = None) -> Cartesian3D:
    """Apply the transformation: source → target Cartesian coordinates."""
    m = p.matrix()
    v = (point.x, point.y, point.z)
    x, y, z = (p.tx + _dot(m[0], v), p.ty + _dot(m[1], v), p.tz + _dot(m[2], v))
    if trace is not None:
        trace.section("Helmert transformation")
        trace.step("ΔX", p.tx, unit="m")
        trace.step("ΔY", p.ty, unit="m")
        trace.step("ΔZ", p.tz, unit="m")
        trace.step("R_X", p.rx, unit="rad")
        trace.step("R_Y", p.ry, unit="rad")
        trace.step("R_Z", p.rz, unit="rad")
        trace.step("ΔS", p.scale, "scale difference")
        trace.step("X′", x, "ΔX + (1 + ΔS)(X − R_Z Y + R_Y Z)", unit="m")
        trace.step("Y′", y, "ΔY + (1 + ΔS)(R_Z X + Y − R_X Z)", unit="m")
        trace.step("Z′", z, "ΔZ + (1 + ΔS)(−R_Y X + R_X Y + Z)", unit="m")
    return Cartesian3D(x, y, z)


def helmert_inverse(point: Cartesian3D, p: HelmertParams, trace: Trace | None = None) -> Cartesian3D:
    """Undo :func:`helmert` exactly (target → source), by solving the 3×3 system."""
    m = p.matrix()
    rhs = (point.x - p.tx, point.y - p.ty, point.z - p.tz)
    det = _det(m)
    cols = []
    for i in range(3):  # Cramer's rule
        replaced = tuple(tuple(rhs[r] if c == i else m[r][c] for c in range(3)) for r in range(3))
        cols.append(_det(replaced) / det)
    x, y, z = cols
    if trace is not None:
        trace.section("Inverse Helmert transformation")
        trace.step("X", x, "solve (1 + ΔS) R X = X′ − ΔX", unit="m")
        trace.step("Y", y, unit="m")
        trace.step("Z", z, unit="m")
    return Cartesian3D(x, y, z)


def transform_geodetic(
    point: Geodetic,
    source: Ellipsoid,
    target: Ellipsoid,
    p: HelmertParams,
    trace: Trace | None = None,
    inverse: bool = False,
) -> Geodetic:
    """φ, λ, h on ``source`` → φ, λ, h on ``target`` via Cartesian coordinates.

    With ``inverse=True`` the parameters are those of target → source, and are
    applied backwards.
    """
    if trace is not None:
        trace.section(f"Geodetic → Cartesian ({source.name})")
    xyz = point.to_cartesian(source, trace)
    moved = helmert_inverse(xyz, p, trace) if inverse else helmert(xyz, p, trace)
    if trace is not None:
        trace.section(f"Cartesian → geodetic ({target.name})")
    return moved.to_geodetic(target, trace=trace)


def _dot(a: tuple[float, ...], b: tuple[float, ...]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _det(m: tuple[tuple[float, ...], ...]) -> float:
    return (
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    )
