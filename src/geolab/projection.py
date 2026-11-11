"""Projection theory: distortion (Tissot), isometric latitude and Mercator.

Distortion at a point is described by the principal scales m₁ and m₂, the
largest and smallest scale, in two perpendicular directions whose images are
also perpendicular (Tissot's theorem). In a direction at angle T from the
first principal direction on the original surface,

    m² = m₁² cos²T + m₂² sin²T          tan T* = (m₂/m₁) tan T

where T* is the image of that direction. The angle distortion is T − T*; it is
largest, sin ω = (m₁ − m₂)/(m₁ + m₂), when tan T = √(m₁/m₂). The area scale is
Φ = m₁m₂, and the course also uses 2h² = 1/m₁² + 1/m₂² and k² = 1/(m₁²m₂²).
A projection is conformal when m₁ = m₂ (ω = 0) and equal-area when Φ = 1.

Isometric latitude q makes φ, λ isometric parameters, ds² = N²cos²φ (dq² + dλ²):

    q = ∫ M/(N cos φ) dφ = atanh(sin φ) − e atanh(e sin φ)

Near a latitude φ₀, Δq = Σ Cₖ Δφᵏ and Δφ = Σ dₖ Δqᵏ (k = 1…5). Mercator on the
sphere is the conformal map x = R q = R ln tan(π/4 + φ/2), y = R λ.

Lengths are in metres and angles in radians.
"""

import math
from dataclasses import dataclass

from .ellipsoids import Ellipsoid
from .trace import Trace


@dataclass(frozen=True)
class Tissot:
    scale: float
    """m, the scale in direction T."""
    image_angle: float
    """T*, the image of the direction T."""
    angle_distortion: float
    """T − T*."""
    max_angle_distortion: float
    """ω, the largest T − T* (any angle changes by at most 2ω)."""
    max_direction: float
    """T where T − T* = ω: tan T = √(m₁/m₂)."""
    area_scale: float
    """Φ = m₁m₂."""
    h: float
    """2h² = 1/m₁² + 1/m₂²."""
    k: float
    """k² = 1/(m₁²m₂²)."""


def tissot(m1: float, m2: float, direction: float, trace: Trace | None = None) -> Tissot:
    """Distortion from the principal scales m₁, m₂ in the direction T (from m₁'s axis)."""
    if m1 <= 0 or m2 <= 0:
        raise ValueError("principal scales must be positive")
    c, s = math.cos(direction), math.sin(direction)
    m = math.sqrt(m1**2 * c**2 + m2**2 * s**2)
    image = math.atan2(m2 * s, m1 * c)
    omega = math.asin((m1 - m2) / (m1 + m2))
    t_max = math.atan(math.sqrt(m1 / m2))
    area = m1 * m2
    h = math.sqrt((1 / m1**2 + 1 / m2**2) / 2)
    k = 1 / (m1 * m2)
    if trace is not None:
        trace.section("Direction T")
        trace.step("m", m, "√(m₁² cos²T + m₂² sin²T)")
        trace.step("T*", image, "atan((m₂/m₁) tan T)", unit="rad")
        trace.step("T − T*", direction - image, unit="rad")
        trace.section("Largest angle distortion")
        trace.step("ω", omega, "asin((m₁ − m₂)/(m₁ + m₂))", unit="rad")
        trace.step("T_ω", t_max, "atan √(m₁/m₂) = 45° + ω/2", unit="rad")
        trace.step("T*_ω", math.atan(math.sqrt(m2 / m1)), "atan √(m₂/m₁) = 45° − ω/2", unit="rad")
        trace.section("Area and mean distortion")
        trace.step("Φ", area, "m₁ m₂")
        trace.step("h", h, "√((1/m₁² + 1/m₂²)/2)")
        trace.step("k", k, "1/(m₁ m₂)")
    return Tissot(m, image, direction - image, omega, t_max, area, h, k)


def isometric_latitude(lat: float, ellipsoid: Ellipsoid, trace: Trace | None = None) -> float:
    """q = atanh(sin φ) − e atanh(e sin φ) = ln[tan(π/4 + φ/2) ((1 − e sin φ)/(1 + e sin φ))^(e/2)]."""
    e = math.sqrt(ellipsoid.e_sq)
    s = math.sin(lat)
    q = math.atanh(s) - e * math.atanh(e * s)
    if trace is not None:
        trace.step("q", q, "atanh(sin φ) − e atanh(e sin φ)")
    return q


def latitude_from_isometric(q: float, ellipsoid: Ellipsoid, trace: Trace | None = None) -> float:
    """φ from q, iterating sin φ = tanh(q + e atanh(e sin φ))."""
    e = math.sqrt(ellipsoid.e_sq)
    s = math.tanh(q)
    for i in range(1, 30):
        s_new = math.tanh(q + e * math.atanh(e * s))
        done = abs(s_new - s) < 1e-16
        s = s_new
        if done:
            break
    lat = math.asin(s)
    if trace is not None:
        trace.step("φ", lat, f"sin φ = tanh(q + e atanh(e sin φ)), {i} iterations", unit="rad")
    return lat


def _aux(lat: float, ellipsoid: Ellipsoid) -> tuple[float, float, float]:
    """cos φ, t = tan φ and η² at φ."""
    c = math.cos(lat)
    return c, math.tan(lat), ellipsoid.ep_sq * c**2


def isometric_coefficients(lat0: float, ellipsoid: Ellipsoid) -> tuple[float, ...]:
    """C₁…C₅ at φ₀: Cₖ = (1/k!) dᵏq/dφᵏ; C₁ exact, the rest to η⁴."""
    c, t, h = _aux(lat0, ellipsoid)
    return (
        1 / ((1 + h) * c),
        t / (2 * c) * (1 + h - 3 * h**2),
        (1 + 2 * t**2 + h - 3 * h**2 + 6 * h**2 * t**2) / (6 * c),
        t / (24 * c) * (5 + 6 * t**2 - h + 21 * h**2 - 6 * h**2 * t**2),
        (5 + 28 * t**2 + 24 * t**4 - h + 21 * h**2 - 60 * h**2 * t**2) / (120 * c),
    )


def latitude_coefficients(lat0: float, ellipsoid: Ellipsoid) -> tuple[float, ...]:
    """d₁…d₅ at φ₀: dₖ = (1/k!) dᵏφ/dqᵏ, to η⁴."""
    c, t, h = _aux(lat0, ellipsoid)
    return (
        c * (1 + h),
        -(c**2) * t / 2 * (1 + 4 * h + 3 * h**2),
        -(c**3) / 6 * (1 - t**2 + 5 * h - 13 * h * t**2 + 7 * h**2 - 27 * h**2 * t**2),
        c**4 * t / 24 * (5 - t**2 + 56 * h - 40 * h * t**2 + 154 * h**2 - 174 * h**2 * t**2),
        c**5 / 120 * (5 - 18 * t**2 + t**4 + h * (61 - 418 * t**2 + 121 * t**4)),
    )


_C_FORMULAS = (
    "1 / ((1 + η²) cos φ₀)",
    "t (1 + η² − 3η⁴) / (2 cos φ₀)",
    "(1 + 2t² + η² − 3η⁴ + 6η⁴t²) / (6 cos φ₀)",
    "t (5 + 6t² − η² + …) / (24 cos φ₀)",
    "(5 + 28t² + 24t⁴ − η² + …) / (120 cos φ₀)",
)
_D_FORMULAS = (
    "cos φ₀ (1 + η²)",
    "−cos²φ₀ t (1 + 4η² + 3η⁴) / 2",
    "−cos³φ₀ (1 − t² + 5η² − 13η²t² + 7η⁴ − 27η⁴t²) / 6",
    "cos⁴φ₀ t (5 − t² + 56η² − 40η²t² + …) / 24",
    "cos⁵φ₀ (5 − 18t² + t⁴ + …) / 120",
)


def _setup(lat0: float, ellipsoid: Ellipsoid, trace: Trace) -> None:
    c, t, h = _aux(lat0, ellipsoid)
    trace.section("Setup")
    trace.step("t", t, "tan φ₀")
    trace.step("η²", h, "e′² cos²φ₀")


def isometric_difference(lat: float, lat0: float, ellipsoid: Ellipsoid, trace: Trace | None = None) -> float:
    """Δq = q − q₀ from Δφ = φ − φ₀ by the series C₁Δφ + … + C₅Δφ⁵."""
    dlat = lat - lat0
    coeffs = isometric_coefficients(lat0, ellipsoid)
    dq = sum(ck * dlat ** (k + 1) for k, ck in enumerate(coeffs))
    if trace is not None:
        _setup(lat0, ellipsoid, trace)
        trace.step("Δφ", dlat, "φ − φ₀", unit="rad")
        trace.section("Coefficients")
        for k, (ck, formula) in enumerate(zip(coeffs, _C_FORMULAS), 1):
            trace.step(f"C{k}", ck, formula)
        trace.section("Result")
        trace.step("Δq", dq, "C₁Δφ + C₂Δφ² + C₃Δφ³ + C₄Δφ⁴ + C₅Δφ⁵")
    return dq


def latitude_difference(dq: float, lat0: float, ellipsoid: Ellipsoid, trace: Trace | None = None) -> float:
    """Δφ = φ − φ₀ from Δq = q − q₀ by the series d₁Δq + … + d₅Δq⁵."""
    coeffs = latitude_coefficients(lat0, ellipsoid)
    dlat = sum(dk * dq ** (k + 1) for k, dk in enumerate(coeffs))
    if trace is not None:
        _setup(lat0, ellipsoid, trace)
        trace.section("Coefficients")
        for k, (dk, formula) in enumerate(zip(coeffs, _D_FORMULAS), 1):
            trace.step(f"d{k}", dk, formula)
        trace.section("Result")
        trace.step("Δφ", dlat, "d₁Δq + d₂Δq² + d₃Δq³ + d₄Δq⁴ + d₅Δq⁵", unit="rad")
    return dlat


@dataclass(frozen=True)
class MercatorPoint:
    x: float
    """North (along the meridian): R q."""
    y: float
    """East: R (λ − λ₀)."""
    scale: float
    """m = 1/cos φ in every direction (conformal)."""


def mercator_sphere(
    lat: float, lon: float, radius: float, lon0: float = 0.0, trace: Trace | None = None
) -> MercatorPoint:
    """Mercator projection of the sphere of radius R: x = R ln tan(π/4 + φ/2), y = R(λ − λ₀)."""
    if abs(lat) >= math.pi / 2:
        raise ValueError("the poles are at infinity on the Mercator projection")
    q = math.log(math.tan(math.pi / 4 + lat / 2))
    dl = math.remainder(lon - lon0, math.tau)
    x, y, m = radius * q, radius * dl, 1 / math.cos(lat)
    if trace is not None:
        trace.step("q", q, "ln tan(π/4 + φ/2) = ∫ dφ/cos φ")
        trace.step("x", x, "R q", unit="m")
        trace.step("y", y, "R (λ − λ₀)", unit="m")
        trace.step("m", m, "1/cos φ (m₁ = m₂, conformal)")
        trace.step("Φ", m * m, "m², area scale")
    return MercatorPoint(x, y, m)


def mercator_sphere_inverse(
    x: float, y: float, radius: float, lon0: float = 0.0, trace: Trace | None = None
) -> tuple[float, float]:
    """φ, λ from Mercator x, y on the sphere: φ = 2 atan(eˣᐟᴿ) − π/2, λ = λ₀ + y/R."""
    lat = 2 * math.atan(math.exp(x / radius)) - math.pi / 2
    lon = math.remainder(lon0 + y / radius, math.tau)
    if trace is not None:
        trace.step("φ", lat, "2 atan(exp(x/R)) − π/2", unit="rad")
        trace.step("λ", lon, "λ₀ + y/R", unit="rad")
    return lat, lon
