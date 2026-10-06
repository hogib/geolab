"""Radii of curvature of the ellipsoid at a given geodetic latitude (radians)."""

import math

from .ellipsoids import Ellipsoid
from .trace import Trace


def _w(lat: float, ellipsoid: Ellipsoid) -> float:
    return math.sqrt(1 - ellipsoid.e_sq * math.sin(lat) ** 2)


def meridian_radius(lat: float, ellipsoid: Ellipsoid, trace: Trace | None = None) -> float:
    """Radius of curvature in the meridian, M = a(1 − e²) / W³."""
    w = _w(lat, ellipsoid)
    m = ellipsoid.a * (1 - ellipsoid.e_sq) / w**3
    if trace is not None:
        trace.step("W", w, "√(1 − e² sin²φ)")
        trace.step("M", m, "a(1 − e²) / W³")
    return m


def prime_vertical_radius(lat: float, ellipsoid: Ellipsoid, trace: Trace | None = None) -> float:
    """Radius of curvature in the prime vertical, N = a / W."""
    w = _w(lat, ellipsoid)
    n = ellipsoid.a / w
    if trace is not None:
        trace.step("W", w, "√(1 − e² sin²φ)")
        trace.step("N", n, "a / W")
    return n


def parallel_radius(lat: float, ellipsoid: Ellipsoid, trace: Trace | None = None) -> float:
    """Radius of the parallel circle, r = N cos φ."""
    n = prime_vertical_radius(lat, ellipsoid, trace)
    r = n * math.cos(lat)
    if trace is not None:
        trace.step("r", r, "N cos φ")
    return r


def gaussian_mean_radius(lat: float, ellipsoid: Ellipsoid, trace: Trace | None = None) -> float:
    """Gaussian mean radius of curvature, R = √(MN)."""
    m = meridian_radius(lat, ellipsoid, trace)
    n = prime_vertical_radius(lat, ellipsoid, trace)
    r = math.sqrt(m * n)
    if trace is not None:
        trace.step("R", r, "√(MN)")
    return r


def azimuth_radius(
    lat: float, azimuth: float, ellipsoid: Ellipsoid, trace: Trace | None = None
) -> float:
    """Radius of curvature of the normal section at ``azimuth`` (Euler's formula)."""
    m = meridian_radius(lat, ellipsoid, trace)
    n = prime_vertical_radius(lat, ellipsoid, trace)
    r = m * n / (n * math.cos(azimuth) ** 2 + m * math.sin(azimuth) ** 2)
    if trace is not None:
        trace.step("Rα", r, "MN / (N cos²α + M sin²α)")
    return r
