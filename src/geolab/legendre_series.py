"""Direct and inverse problems for short lines by Legendre's series (1806).

The geodesic's φ, λ and α are expanded in Taylor series of the length s from
P₁ (φ₁, α₁₂). With u = s cos α₁₂ and v = s sin α₁₂ and, at φ₁,

    t = tan φ,  η² = e′² cos²φ,  V² = 1 + η²,  c = a²/b

    Δφ = b₁u + b₂v² + b₃u² + b₄uv² + b₅u³ + b₆v⁴ + b₇u²v² + b₈uv⁴ + b₉u³v²
    Δλ = l₁v + l₂uv + l₃v³ + l₄u²v + l₅uv³ + l₆u³v + l₇v⁵ + l₈u²v³ + l₉u⁴v
    α₂₁ = α₁₂ ± 180° + a₁v + a₂uv + a₃v³ + a₄u²v + a₅uv³ + a₆u³v

Coefficients are those of the course slides (with ρ = 1, so results are in
radians), with two corrections the slides need: b₆ has η²(1 − 9t²), not
η²(t − 9t²), and l₇ has t²(1 + 3t²), not t(1 + 3t²). Mirroring a line in the
equator (t → −t, u → −u) must flip Δφ and keep Δλ, which fixes both, and the
corrected forms agree with Vincenty. Good to about 0.0001″ for s ≈ 100 km at
mid latitudes.

The inverse problem inverts the same series: from the slides' leading terms
for u and v, then by fixed-point iteration on the full series.
"""

import math
from dataclasses import dataclass

from .ellipsoids import Ellipsoid
from .errors import ConvergenceError
from .trace import Trace

ARCSEC = math.pi / (180 * 3600)


@dataclass(frozen=True)
class LegendreCoefficients:
    b: tuple[float, ...]
    """b₁ … b₉ (radians per metre^k)."""
    l: tuple[float, ...]
    """l₁ … l₉."""
    a: tuple[float, ...]
    """a₁ … a₆."""


def legendre_coefficients(lat: float, ellipsoid: Ellipsoid) -> LegendreCoefficients:
    c = ellipsoid.a**2 / ellipsoid.b
    t = math.tan(lat)
    t2, t4 = t * t, t**4
    eta2 = ellipsoid.ep_sq * math.cos(lat) ** 2
    v2 = 1 + eta2
    v = math.sqrt(v2)
    cosf = math.cos(lat)
    b = (
        v**3 / c,
        -(v**4) * t / (2 * c**2),
        -3 * v**4 * eta2 * t / (2 * c**2),
        -(v**5) / (6 * c**3) * ((1 + 3 * t2) + eta2 * (1 - 9 * t2)),
        v**5 / (2 * c**3) * eta2 * (-1 + t2),
        v**6 / (24 * c**4) * t * ((1 + 3 * t2) + eta2 * (1 - 9 * t2)),
        -(v**6) / (12 * c**4) * t * ((4 + 6 * t2) - eta2 * (13 + 9 * t2)),
        v**7 / (120 * c**5) * (1 + 30 * t2 + 45 * t4),
        -(v**7) / (30 * c**5) * (2 + 15 * t2 + 15 * t4),
    )
    l = (
        v / (c * cosf),
        v2 / (c**2 * cosf) * t,
        -(v**3) / (3 * c**3 * cosf) * t2,
        v**3 / (3 * c**3 * cosf) * (1 + 3 * t2 + eta2),
        -(v**4) / (3 * c**4 * cosf) * t * (1 + 3 * t2),
        v**4 / (3 * c**4 * cosf) * t * (2 + 3 * t2),
        v**5 / (15 * c**5 * cosf) * t2 * (1 + 3 * t2),
        -(v**5) / (15 * c**5 * cosf) * (1 + 20 * t2 + 30 * t4),
        v**5 / (15 * c**5 * cosf) * (2 + 15 * t2 + 15 * t4),
    )
    a = (
        v * t / c,
        v2 / (2 * c**2) * (1 + 2 * t2 + eta2),
        -(v**3) / (6 * c**3) * t * (1 + 2 * t2),
        v**3 / (6 * c**3) * t * (5 + 6 * t2),
        -(v**4) / (24 * c**4) * (1 + 20 * t2 + 24 * t4),
        v**4 / (24 * c**4) * (5 + 28 * t2 + 24 * t4),
    )
    return LegendreCoefficients(b, l, a)


def _terms(k: LegendreCoefficients, u: float, v: float) -> tuple[list[float], list[float], list[float]]:
    b, l, a = k.b, k.l, k.a
    dphi = [
        b[0] * u, b[1] * v**2, b[2] * u**2, b[3] * u * v**2, b[4] * u**3,
        b[5] * v**4, b[6] * u**2 * v**2, b[7] * u * v**4, b[8] * u**3 * v**2,
    ]
    dlam = [
        l[0] * v, l[1] * u * v, l[2] * v**3, l[3] * u**2 * v, l[4] * u * v**3,
        l[5] * u**3 * v, l[6] * v**5, l[7] * u**2 * v**3, l[8] * u**4 * v,
    ]
    dalpha = [
        a[0] * v, a[1] * u * v, a[2] * v**3, a[3] * u**2 * v, a[4] * u * v**3, a[5] * u**3 * v,
    ]
    return dphi, dlam, dalpha


_DPHI_NAMES = ["b₁u", "b₂v²", "b₃u²", "b₄uv²", "b₅u³", "b₆v⁴", "b₇u²v²", "b₈uv⁴", "b₉u³v²"]
_DLAM_NAMES = ["l₁v", "l₂uv", "l₃v³", "l₄u²v", "l₅uv³", "l₆u³v", "l₇v⁵", "l₈u²v³", "l₉u⁴v"]
_DALPHA_NAMES = ["a₁v", "a₂uv", "a₃v³", "a₄u²v", "a₅uv³", "a₆u³v"]


@dataclass(frozen=True)
class LegendreDirect:
    lat2: float
    lon2: float
    back_azimuth: float
    """α₂₁, the azimuth from P₂ back to P₁."""
    u: float
    v: float
    dlat_terms: tuple[float, ...]
    dlon_terms: tuple[float, ...]
    dazimuth_terms: tuple[float, ...]


def legendre_direct(
    lat1: float,
    lon1: float,
    azimuth: float,
    distance: float,
    ellipsoid: Ellipsoid,
    trace: Trace | None = None,
) -> LegendreDirect:
    """First geodetic problem: φ₂, λ₂, α₂₁ from φ₁, λ₁, α₁₂ and s."""
    k = legendre_coefficients(lat1, ellipsoid)
    u, v = distance * math.cos(azimuth), distance * math.sin(azimuth)
    dphi, dlam, dalpha = _terms(k, u, v)
    lat2 = lat1 + sum(dphi)
    lon2 = math.remainder(lon1 + sum(dlam), math.tau)
    back = (azimuth + math.pi + sum(dalpha)) % math.tau
    if trace is not None:
        trace.section("Setup")
        trace.step("u", u, "s cos α₁₂", unit="m")
        trace.step("v", v, "s sin α₁₂", unit="m")
        for title, names, terms, total in (
            ("Δφ terms (″)", _DPHI_NAMES, dphi, "Δφ"),
            ("Δλ terms (″)", _DLAM_NAMES, dlam, "Δλ"),
            ("Δα terms (″)", _DALPHA_NAMES, dalpha, "Δα"),
        ):
            trace.section(title)
            for name, term in zip(names, terms):
                trace.step(name, term / ARCSEC, "″")
            trace.step(total, sum(terms), "sum", unit="rad")
        trace.section("Result")
        trace.step("φ₂", lat2, "φ₁ + Δφ", unit="rad")
        trace.step("λ₂", lon2, "λ₁ + Δλ", unit="rad")
        trace.step("α₂₁", back, "α₁₂ + 180° + Δα", unit="rad")
    return LegendreDirect(lat2, lon2, back, u, v, tuple(dphi), tuple(dlam), tuple(dalpha))


@dataclass(frozen=True)
class LegendreInverse:
    distance: float
    azimuth: float
    """α₁₂."""
    back_azimuth: float
    """α₂₁."""
    iterations: int


def legendre_inverse(
    lat1: float,
    lon1: float,
    lat2: float,
    lon2: float,
    ellipsoid: Ellipsoid,
    trace: Trace | None = None,
    tol: float = 1e-6,
    max_iter: int = 50,
) -> LegendreInverse:
    """Second geodetic problem by inverting the Legendre series for u, v."""
    k = legendre_coefficients(lat1, ellipsoid)
    c = ellipsoid.a**2 / ellipsoid.b
    t = math.tan(lat1)
    eta2 = ellipsoid.ep_sq * math.cos(lat1) ** 2
    v_ = math.sqrt(1 + eta2)
    dphi = lat2 - lat1
    dlam = math.remainder(lon2 - lon1, math.tau)
    cosf = math.cos(lat1)
    # Leading terms of the inverted series (slides).
    u = c / v_**3 * dphi + 3 * c * eta2 * t / (2 * v_**5) * dphi**2 + c / v_ * t * cosf**2 / 2 * dlam**2
    v = c * cosf / v_ * dlam - c / v_**3 * t * cosf * dlam * dphi
    if trace is not None:
        trace.section("First approximation")
        trace.step("Δφ", dphi, "φ₂ − φ₁", unit="rad")
        trace.step("Δλ", dlam, "λ₂ − λ₁", unit="rad")
        trace.step("u⁰", u, "c/V³ Δφ + 3cη²t/(2V⁵) Δφ² + c t cos²φ/(2V) Δλ²", unit="m")
        trace.step("v⁰", v, "c cos φ/V Δλ − c t cos φ/V³ Δλ Δφ", unit="m")
    for i in range(1, max_iter + 1):
        p, q, _ = _terms(k, u, v)
        new_u = (dphi - (sum(p) - p[0])) / k.b[0]
        new_v = (dlam - (sum(q) - q[0])) / k.l[0]
        done = abs(new_u - u) < tol and abs(new_v - v) < tol
        u, v = new_u, new_v
        if trace is not None:
            trace.section(f"Iteration {i}")
            trace.step("u", u, "(Δφ − b₂v² − b₃u² − …) / b₁", unit="m")
            trace.step("v", v, "(Δλ − l₂uv − l₃v³ − …) / l₁", unit="m")
        if done:
            break
    else:
        raise ConvergenceError(f"u, v did not converge in {max_iter} iterations")
    azimuth = math.atan2(v, u) % math.tau
    distance = math.hypot(u, v)
    _, _, dalpha = _terms(k, u, v)
    back = (azimuth + math.pi + sum(dalpha)) % math.tau
    if trace is not None:
        trace.section("Result")
        trace.step("α₁₂", azimuth, "atan2(v, u)", unit="rad")
        trace.step("s", distance, "√(u² + v²)", unit="m")
        trace.step("Δα", sum(dalpha), "a₁v + a₂uv + …", unit="rad")
        trace.step("α₂₁", back, "α₁₂ + 180° + Δα", unit="rad")
    return LegendreInverse(distance, azimuth, back, i)
