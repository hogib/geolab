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


# ── The course's series in e′² (Schödlbauer) ────────────────────────────────
#
# 1/v³ = A′ − B′ cos 2φ + C′ cos 4φ − D′ cos 6φ + E′ cos 8φ − …, integrated term
# by term: G = E₀φ + E₂ sin 2φ + E₄ sin 4φ + E₆ sin 6φ + E₈ sin 8φ, with
# E₀ = cA′, E₂ = −cB′/2, E₄ = cC′/4, E₆ = −cD′/6, E₈ = cE′/8 (E₀ per radian here).
# The inverse is φ = σ + F₂ sin 2σ + F₄ sin 4σ + F₆ sin 6σ with σ = G/E₀.


@dataclass(frozen=True)
class MeridianSeries:
    a_: float
    b_: float
    c_: float
    d_: float
    e_: float
    """A′, B′, C′, D′, E′: the cosine series of 1/v³ (A′ + B′ + C′ + D′ + E′ = 1)."""
    e0: float
    """E₀ in metres per radian (multiply by π/180 for m/°)."""
    e2: float
    e4: float
    e6: float
    e8: float
    f2: float
    f4: float
    f6: float
    """Inverse-series coefficients F₂, F₄, F₆, in radians."""


def meridian_series(ellipsoid: Ellipsoid) -> MeridianSeries:
    """Coefficients of the course's meridian-arc series, to e′⁸."""
    c = ellipsoid.a**2 / ellipsoid.b
    k = ellipsoid.ep_sq
    k2, k3, k4 = k**2, k**3, k**4
    a_ = 1 - 3 / 4 * k + 45 / 64 * k2 - 175 / 256 * k3 + 11025 / 16384 * k4
    b_ = 3 / 4 * k - 15 / 16 * k2 + 525 / 512 * k3 - 2205 / 2048 * k4
    c_ = 15 / 64 * k2 - 105 / 256 * k3 + 2205 / 4096 * k4
    d_ = 35 / 512 * k3 - 315 / 2048 * k4
    e_ = 315 / 16384 * k4
    f2 = 3 / 8 * k - 3 / 16 * k2 + 213 / 2048 * k3 - 255 / 4096 * k4
    f4 = 21 / 256 * k2 - 21 / 256 * k3 + 533 / 8192 * k4
    f6 = 151 / 6144 * k3 - 453 / 12288 * k4
    return MeridianSeries(
        a_, b_, c_, d_, e_, c * a_, -c * b_ / 2, c * c_ / 4, -c * d_ / 6, c * e_ / 8, f2, f4, f6
    )


def meridian_distance_series(
    lat: float, ellipsoid: Ellipsoid, trace: Trace | None = None
) -> float:
    """G = E₀φ + E₂ sin 2φ + E₄ sin 4φ + E₆ sin 6φ + E₈ sin 8φ (the course's form)."""
    s = meridian_series(ellipsoid)
    g = (
        s.e0 * lat
        + s.e2 * math.sin(2 * lat)
        + s.e4 * math.sin(4 * lat)
        + s.e6 * math.sin(6 * lat)
        + s.e8 * math.sin(8 * lat)
    )
    if trace is not None:
        trace.step("E₀", s.e0 * math.pi / 180, "c A′ / ρ°", unit="m")
        trace.step("E₂", s.e2, "−c B′ / 2", unit="m")
        trace.step("E₄", s.e4, "c C′ / 4", unit="m")
        trace.step("E₆", s.e6, "−c D′ / 6", unit="m")
        trace.step("E₈", s.e8, "c E′ / 8", unit="m")
        trace.step("G", g, "E₀φ + E₂ sin 2φ + E₄ sin 4φ + E₆ sin 6φ + E₈ sin 8φ", unit="m")
    return g


def latitude_from_meridian_distance_series(
    g: float, ellipsoid: Ellipsoid, trace: Trace | None = None
) -> float:
    """φ = σ + F₂ sin 2σ + F₄ sin 4σ + F₆ sin 6σ, σ = G/E₀ (the course's form, no iteration)."""
    s = meridian_series(ellipsoid)
    sigma = g / s.e0
    lat = sigma + s.f2 * math.sin(2 * sigma) + s.f4 * math.sin(4 * sigma) + s.f6 * math.sin(6 * sigma)
    if trace is not None:
        trace.step("σ", sigma, "G / E₀", unit="rad")
        trace.step("F₂", s.f2, "3/8 e′² − 3/16 e′⁴ + 213/2048 e′⁶ − 255/4096 e′⁸", unit="rad")
        trace.step("F₄", s.f4, "21/256 e′⁴ − 21/256 e′⁶ + 533/8192 e′⁸", unit="rad")
        trace.step("F₆", s.f6, "151/6144 e′⁶ − 453/12288 e′⁸", unit="rad")
        trace.step("φ", lat, "σ + F₂ sin 2σ + F₄ sin 4σ + F₆ sin 6σ", unit="rad")
    return lat


def short_meridian_arc(
    lat1: float, lat2: float, ellipsoid: Ellipsoid, trace: Trace | None = None
) -> float:
    """Meridian arc for short arcs (under about 1°), from the mean latitude φ₀:
    ΔG = E₀Δφ + 2E₂ sin Δφ cos 2φ₀ + 2E₄ sin 2Δφ cos 4φ₀ + 2E₆ sin 3Δφ cos 6φ₀ (signed)."""
    s = meridian_series(ellipsoid)
    dlat = lat2 - lat1
    lat0 = (lat1 + lat2) / 2
    dg = (
        s.e0 * dlat
        + 2 * s.e2 * math.sin(dlat) * math.cos(2 * lat0)
        + 2 * s.e4 * math.sin(2 * dlat) * math.cos(4 * lat0)
        + 2 * s.e6 * math.sin(3 * dlat) * math.cos(6 * lat0)
    )
    if trace is not None:
        trace.step("Δφ", dlat, "φ₂ − φ₁", unit="rad")
        trace.step("φ₀", lat0, "(φ₁ + φ₂) / 2", unit="rad")
        trace.step(
            "ΔG", dg, "E₀Δφ + 2E₂ sin Δφ cos 2φ₀ + 2E₄ sin 2Δφ cos 4φ₀ + 2E₆ sin 3Δφ cos 6φ₀", unit="m"
        )
    return dg
