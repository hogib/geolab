"""The geodesic as a curve: its differential equations and Clairaut's theorem.

On the ellipsoid (E = M², F = 0, G = N² cos²φ) a geodesic satisfies

    dφ/ds = cos α / M        dλ/ds = sin α / (N cos φ)        dα/ds = tan φ sin α / N

and Clairaut's theorem: along it, the radius of the parallel times the sine of
the azimuth is constant,

    k = N cos φ sin α = a cos β sin α = const.
"""

import math
from dataclasses import dataclass

from .curvature import meridian_radius, prime_vertical_radius
from .ellipsoids import Ellipsoid
from .trace import Trace


def clairaut_constant(lat: float, azimuth: float, ellipsoid: Ellipsoid, trace: Trace | None = None) -> float:
    """k = N cos φ sin α (m)."""
    n = prime_vertical_radius(lat, ellipsoid)
    k = n * math.cos(lat) * math.sin(azimuth)
    if trace is not None:
        trace.step("N", n, "a / √(1 − e² sin²φ)", unit="m")
        trace.step("k", k, "N cos φ sin α", unit="m")
    return k


def vertex_latitude(k: float, ellipsoid: Ellipsoid, trace: Trace | None = None) -> float:
    """Highest latitude the geodesic reaches (where α = 90°), from
    cos φ = k / √(c² − k² e′²). The lowest is its negative."""
    c = ellipsoid.a**2 / ellipsoid.b
    k = abs(k)
    if k > ellipsoid.a:
        raise ValueError("|k| cannot exceed a")
    lat = math.acos(k / math.sqrt(c**2 - k**2 * ellipsoid.ep_sq))
    if trace is not None:
        trace.step("c", c, "a² / b", unit="m")
        trace.step("φ_max", lat, "acos(k / √(c² − k² e′²))", unit="rad")
    return lat


def azimuths_at_latitude(
    k: float, lat: float, ellipsoid: Ellipsoid, trace: Trace | None = None
) -> tuple[float, float]:
    """The two azimuths (northward α and southward 180° − α, or their mirror images
    for k < 0) the geodesic has where it crosses latitude φ: sin α = k / (N cos φ)."""
    r = prime_vertical_radius(lat, ellipsoid) * math.cos(lat)
    ratio = k / r
    if abs(ratio) > 1 + 1e-15:
        raise ValueError("the geodesic does not reach this latitude")
    alpha = math.asin(max(-1.0, min(1.0, ratio)))
    pair = (alpha % math.tau, (math.pi - alpha) % math.tau)
    if trace is not None:
        trace.step("r", r, "N cos φ", unit="m")
        trace.step("α₁", pair[0], "asin(k / r)", unit="rad")
        trace.step("α₂", pair[1], "180° − α₁", unit="rad")
    return pair


@dataclass(frozen=True)
class GeodesicPoint:
    lat: float
    lon: float
    azimuth: float


def integrate_geodesic(
    lat: float,
    lon: float,
    azimuth: float,
    distance: float,
    ellipsoid: Ellipsoid,
    steps: int = 1000,
) -> GeodesicPoint:
    """Follow the geodesic for ``distance`` metres by integrating its differential
    equations with the classical 4th-order Runge–Kutta method."""

    def deriv(state: tuple[float, float, float]) -> tuple[float, float, float]:
        phi, _, alpha = state
        m = meridian_radius(phi, ellipsoid)
        n = prime_vertical_radius(phi, ellipsoid)
        return (
            math.cos(alpha) / m,
            math.sin(alpha) / (n * math.cos(phi)),
            math.tan(phi) * math.sin(alpha) / n,
        )

    def add(state, d, h):
        return tuple(x + h * dx for x, dx in zip(state, d))

    h = distance / steps
    state = (lat, lon, azimuth)
    for _ in range(steps):
        k1 = deriv(state)
        k2 = deriv(add(state, k1, h / 2))
        k3 = deriv(add(state, k2, h / 2))
        k4 = deriv(add(state, k3, h))
        state = tuple(
            x + h / 6 * (a + 2 * b + 2 * c + d) for x, a, b, c, d in zip(state, k1, k2, k3, k4)
        )
    phi, lam, alpha = state
    return GeodesicPoint(phi, math.remainder(lam, math.tau), alpha % math.tau)
