"""Arc lengths along a meridian and along a parallel (all angles in radians)."""

import math
from dataclasses import dataclass

from .curvature import meridian_radius, parallel_radius, prime_vertical_radius
from .ellipsoids import Ellipsoid
from .errors import ConvergenceError
from .trace import Trace


def meridian_distance(lat: float, ellipsoid: Ellipsoid, trace: Trace | None = None) -> float:
    """Meridian arc length from the equator to ``lat`` (signed, metres).

    Helmert's series in the third flattening n, to n⁴; the truncation error
    is below a micrometre on Earth-like ellipsoids.
    """
    a, n = ellipsoid.a, ellipsoid.n
    n2, n3, n4 = n**2, n**3, n**4
    a0 = 1 + n2 / 4 + n4 / 64
    a2 = -3 / 2 * (n - n3 / 8)
    a4 = 15 / 16 * (n2 - n4 / 4)
    a6 = -35 / 48 * n3
    a8 = 315 / 512 * n4
    m = a / (1 + n) * (
        a0 * lat
        + a2 * math.sin(2 * lat)
        + a4 * math.sin(4 * lat)
        + a6 * math.sin(6 * lat)
        + a8 * math.sin(8 * lat)
    )
    if trace is not None:
        trace.step("A₀", a0, "1 + n²/4 + n⁴/64")
        trace.step("A₂", a2, "−3/2 (n − n³/8)")
        trace.step("A₄", a4, "15/16 (n² − n⁴/4)")
        trace.step("A₆", a6, "−35/48 n³")
        trace.step("A₈", a8, "315/512 n⁴")
        trace.step("m", m, "a/(1 + n) (A₀φ + A₂ sin 2φ + A₄ sin 4φ + A₆ sin 6φ + A₈ sin 8φ)", unit="m")
    return m


def meridian_arc_length(
    lat1: float, lat2: float, ellipsoid: Ellipsoid, trace: Trace | None = None
) -> float:
    """Length of the meridian arc between latitudes ``lat1`` and ``lat2`` (metres)."""
    if trace is not None:
        trace.section("Meridian distance to φ₁")
    m1 = meridian_distance(lat1, ellipsoid, trace)
    if trace is not None:
        trace.section("Meridian distance to φ₂")
    m2 = meridian_distance(lat2, ellipsoid, trace)
    s = abs(m2 - m1)
    if trace is not None:
        trace.section("Meridian arc")
        trace.step("Δm", s, "|m₂ − m₁|", unit="m")
    return s


def parallel_arc_length(
    lat: float, lon1: float, lon2: float, ellipsoid: Ellipsoid, trace: Trace | None = None
) -> float:
    """Length of the (shorter) arc of the parallel at ``lat`` between two longitudes (metres)."""
    dlon = abs(math.remainder(lon2 - lon1, math.tau))
    r = parallel_radius(lat, ellipsoid, trace)
    s = r * dlon
    if trace is not None:
        trace.step("Δλ", dlon, "|λ₂ − λ₁|, the shorter way round", unit="rad")
        trace.step("s", s, "r Δλ", unit="m")
    return s


@dataclass(frozen=True)
class MeridianPoint:
    lat: float
    """Geodetic latitude φ."""
    reduced_lat: float
    """Parametric (reduced) latitude β."""
    geocentric_lat: float
    """Geocentric latitude ψ."""
    x: float
    """Distance from the polar axis in the meridian plane, N cos φ (m)."""
    y: float
    """Distance from the equatorial plane, N (1 − e²) sin φ (m)."""
    iterations: int


def latitude_from_meridian_distance(
    m: float,
    ellipsoid: Ellipsoid,
    trace: Trace | None = None,
    tol: float = 1e-12,
    max_iter: int = 20,
) -> MeridianPoint:
    """Latitude and meridian-ellipse coordinates of the point at meridian distance ``m``
    from the equator (the inverse of :func:`meridian_distance`).

    A first value comes from the inverse series in n (footpoint latitude), then
    Newton's method on m(φ) − m = 0 with dm/dφ = M refines it to ``tol`` radians.
    """
    a, n = ellipsoid.a, ellipsoid.n
    quarter = meridian_distance(math.pi / 2, ellipsoid)
    if abs(m) > quarter * (1 + 1e-15):
        raise ValueError(f"|m| exceeds the quarter meridian ({quarter:.4f} m)")

    n2, n3, n4 = n**2, n**3, n**4
    g = a / (1 + n) * (1 + n2 / 4 + n4 / 64)
    mu = m / g
    lat = (
        mu
        + (3 / 2 * n - 27 / 32 * n3) * math.sin(2 * mu)
        + (21 / 16 * n2 - 55 / 32 * n4) * math.sin(4 * mu)
        + 151 / 96 * n3 * math.sin(6 * mu)
        + 1097 / 512 * n4 * math.sin(8 * mu)
    )
    lat = max(-math.pi / 2, min(math.pi / 2, lat))

    if trace is not None:
        trace.section("Initial value (inverse series)")
        trace.step("G", g, "a/(1 + n) (1 + n²/4 + n⁴/64)", unit="m")
        trace.step("μ", mu, "m / G", unit="rad")
        trace.step(
            "φ₀",
            lat,
            "μ + (3n/2 − 27n³/32) sin 2μ + (21n²/16 − 55n⁴/32) sin 4μ + 151n³/96 sin 6μ + 1097n⁴/512 sin 8μ",
            unit="rad",
        )

    for i in range(1, max_iter + 1):
        residual = meridian_distance(lat, ellipsoid) - m
        step = residual / meridian_radius(lat, ellipsoid)
        lat = max(-math.pi / 2, min(math.pi / 2, lat - step))
        if trace is not None:
            trace.section(f"Newton iteration {i}")
            trace.step("m(φ) − m", residual, "", unit="m")
            trace.step("φ", lat, "φ − (m(φ) − m) / M(φ)", unit="rad")
        if abs(step) < tol:
            break
    else:
        raise ConvergenceError(f"φ did not converge in {max_iter} iterations")

    big_n = prime_vertical_radius(lat, ellipsoid)
    x = big_n * math.cos(lat)
    y = big_n * (1 - ellipsoid.e_sq) * math.sin(lat)
    beta = math.atan2((1 - ellipsoid.f) * math.sin(lat), math.cos(lat))
    psi = math.atan2((1 - ellipsoid.e_sq) * math.sin(lat), math.cos(lat))

    if trace is not None:
        trace.section("Coordinates on the meridian ellipse")
        trace.step("N", big_n, "a / √(1 − e² sin²φ)", unit="m")
        trace.step("x", x, "N cos φ", unit="m")
        trace.step("y", y, "N (1 − e²) sin φ", unit="m")
        trace.step("β", beta, "atan((1 − f) tan φ)", unit="rad")
        trace.step("ψ", psi, "atan((1 − e²) tan φ)", unit="rad")

    return MeridianPoint(lat, beta, psi, x, y, i)


def latitude_from_meridian_arc(
    lat1: float, arc: float, ellipsoid: Ellipsoid, trace: Trace | None = None
) -> MeridianPoint:
    """The point on the same meridian as φ₁, ``arc`` metres away along it
    (positive northwards): φ₂ from m(φ₂) = m(φ₁) + Δm."""
    if trace is not None:
        trace.section("Meridian distance to φ₁")
    m1 = meridian_distance(lat1, ellipsoid, trace)
    m2 = m1 + arc
    if trace is not None:
        trace.section("Meridian distance to φ₂")
        trace.step("m₂", m2, "m₁ + Δm", unit="m")
    return latitude_from_meridian_distance(m2, ellipsoid, trace)


def longitude_from_parallel_arc(
    lat: float, lon1: float, arc: float, ellipsoid: Ellipsoid, trace: Trace | None = None
) -> float:
    """Longitude reached after ``arc`` metres along the parallel at ``lat``
    (positive eastwards): λ₂ = λ₁ + s / r. Returned in [−π, π]."""
    r = parallel_radius(lat, ellipsoid, trace)
    if r < 1e-3:  # cos 90° isn't exactly 0 in floating point
        raise ValueError("the parallel at a pole has zero length")
    dlon = arc / r
    lon2 = math.remainder(lon1 + dlon, math.tau)
    if trace is not None:
        trace.step("Δλ", dlon, "s / r", unit="rad")
        trace.step("λ₂", lon2, "λ₁ + Δλ", unit="rad")
    return lon2
