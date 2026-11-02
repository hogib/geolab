"""Points of the meridian ellipse from any latitude type.

In the meridian plane a point has coordinates p (distance from the polar
axis, the radius of the parallel) and z (height above the equator). With the
course's notation, w = √(1 − e² sin²φ) and v = √(1 + e′² cos²φ):

    p = a cos φ / w = c cos φ / v          z = a(1 − e²) sin φ / w
    p = a cos β,  z = b sin β              (reduced latitude β)
    p = r cos ψ,  z = r sin ψ              (geocentric latitude ψ, radius r)
"""

import math
from dataclasses import dataclass

from .ellipsoids import Ellipsoid
from .point import Latitude, LatType
from .trace import Trace


@dataclass(frozen=True)
class MeridianEllipsePoint:
    lat: float
    """Geodetic latitude φ."""
    reduced_lat: float
    """Reduced (parametric) latitude β."""
    geocentric_lat: float
    """Geocentric latitude ψ (γ in some texts)."""
    p: float
    """Distance from the polar axis (radius of the parallel), m."""
    z: float
    """Distance from the equatorial plane, m."""
    r: float
    """Geocentric radius √(p² + z²), m."""
    w: float
    """√(1 − e² sin²φ)."""
    v: float
    """√(1 + e′² cos²φ)."""


def meridian_ellipse_point(
    value: float, lat_type: LatType, ellipsoid: Ellipsoid, trace: Trace | None = None
) -> MeridianEllipsePoint:
    """The point at latitude ``value`` of type ``lat_type`` on the meridian ellipse."""
    source = Latitude(value, lat_type)
    lat = source.convert_to(LatType.GEODETIC, ellipsoid).value
    beta = source.convert_to(LatType.PARAMETRIC, ellipsoid).value
    psi = source.convert_to(LatType.GEOCENTRIC, ellipsoid).value
    a, b, e_sq, ep_sq = ellipsoid.a, ellipsoid.b, ellipsoid.e_sq, ellipsoid.ep_sq
    w = math.sqrt(1 - e_sq * math.sin(lat) ** 2)
    v = math.sqrt(1 + ep_sq * math.cos(lat) ** 2)
    p = a * math.cos(beta)
    z = b * math.sin(beta)
    r = math.hypot(p, z)

    if trace is not None:
        trace.section("Meridian ellipse point")
        trace.step("w", w, "√(1 − e² sin²φ)")
        trace.step("v", v, "√(1 + e′² cos²φ)")
        if lat_type == LatType.GEODETIC:
            trace.step("p", p, "a cos φ / w", unit="m")
            trace.step("z", z, "a(1 − e²) sin φ / w", unit="m")
            trace.step("r", r, "a cos φ / (w cos ψ)", unit="m")
        elif lat_type == LatType.PARAMETRIC:
            trace.step("p", p, "a cos β", unit="m")
            trace.step("z", z, "b sin β", unit="m")
            trace.step("r", r, "√(p² + z²)", unit="m")
        else:
            trace.step("r", r, "ab / √(a² sin²ψ + b² cos²ψ)", unit="m")
            trace.step("p", p, "r cos ψ", unit="m")
            trace.step("z", z, "r sin ψ", unit="m")

    return MeridianEllipsePoint(lat, beta, psi, p, z, r, w, v)
