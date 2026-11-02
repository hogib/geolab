"""Solving plane and spherical triangles, and the small-triangle methods of
classical geodesy (Legendre's theorem and Soldner's additament method).

Notation follows the course: vertices P₁, P₂, P₃ with angles α, β, γ
(β₁, β₂, β₃) and opposite sides a = P₂P₃, b = P₁P₃, c = P₁P₂. Spherical
sides are given as lengths on a sphere of radius R; in geodesy that is the
Gauss sphere R = √(MN) at the middle of the triangle.

Given parts:  ASA (AKA) c, α, β    SAS (KAK) a, b, γ    SSS (KKK) a, b, c
"""

import math
from dataclasses import dataclass

from .trace import Trace


@dataclass(frozen=True)
class Triangle:
    a: float
    b: float
    c: float
    alpha: float
    beta: float
    gamma: float
    excess: float = 0.0
    """Spherical excess ε = α + β + γ − π (0 for a plane triangle)."""
    area: float = 0.0


# ── Plane triangles ──────────────────────────────────────────────────────────


def _plane_area(a: float, b: float, gamma: float) -> float:
    return a * b * math.sin(gamma) / 2


def plane_asa(c: float, alpha: float, beta: float, trace: Trace | None = None) -> Triangle:
    """Plane triangle from side c and its adjacent angles α, β (sine rule)."""
    gamma = math.pi - alpha - beta
    if gamma <= 0:
        raise ValueError("α + β must be less than 180°")
    a = c * math.sin(alpha) / math.sin(gamma)
    b = c * math.sin(beta) / math.sin(gamma)
    if trace is not None:
        trace.step("γ", gamma, "180° − α − β", unit="rad")
        trace.step("a", a, "c sin α / sin γ", unit="m")
        trace.step("b", b, "c sin β / sin γ", unit="m")
    return Triangle(a, b, c, alpha, beta, gamma, 0.0, _plane_area(a, b, gamma))


def plane_sas(a: float, b: float, gamma: float, trace: Trace | None = None) -> Triangle:
    """Plane triangle from sides a, b and the included angle γ."""
    c = math.sqrt(a * a + b * b - 2 * a * b * math.cos(gamma))
    alpha = math.atan2(a * math.sin(gamma), b - a * math.cos(gamma))
    beta = math.pi - alpha - gamma
    if trace is not None:
        trace.step("c", c, "√(a² + b² − 2ab cos γ)", unit="m")
        trace.step("α", alpha, "atan2(a sin γ, b − a cos γ)", unit="rad")
        trace.step("β", beta, "180° − α − γ", unit="rad")
    return Triangle(a, b, c, alpha, beta, gamma, 0.0, _plane_area(a, b, gamma))


def plane_sss(a: float, b: float, c: float, trace: Trace | None = None) -> Triangle:
    """Plane triangle from its three sides (half-angle formulas)."""
    _check_sides(a, b, c)
    u = (a + b + c) / 2
    angles = []
    for opposite, x, y in ((a, b, c), (b, a, c), (c, a, b)):
        angles.append(2 * math.atan(math.sqrt((u - x) * (u - y) / (u * (u - opposite)))))
    alpha, beta, gamma = angles
    area = math.sqrt(u * (u - a) * (u - b) * (u - c))
    if trace is not None:
        trace.step("u", u, "(a + b + c)/2", unit="m")
        trace.step("α", alpha, "2 atan √((u − b)(u − c) / (u(u − a)))", unit="rad")
        trace.step("β", beta, "2 atan √((u − a)(u − c) / (u(u − b)))", unit="rad")
        trace.step("γ", gamma, "2 atan √((u − a)(u − b) / (u(u − c)))", unit="rad")
        trace.step("F", area, "√(u(u − a)(u − b)(u − c))", unit="m²")
    return Triangle(a, b, c, alpha, beta, gamma, 0.0, area)


# ── Spherical triangles, exact ───────────────────────────────────────────────


def _spherical(a: float, b: float, c: float, alpha: float, beta: float, gamma: float, r: float) -> Triangle:
    excess = alpha + beta + gamma - math.pi
    return Triangle(a, b, c, alpha, beta, gamma, excess, excess * r * r)


def spherical_asa(
    c: float, alpha: float, beta: float, radius: float, trace: Trace | None = None
) -> Triangle:
    """Spherical triangle from side c and the angles α, β at its ends (Napier's analogies)."""
    h = c / (2 * radius)
    half_sum, half_diff = (alpha + beta) / 2, (alpha - beta) / 2
    z1 = math.sin(h) * math.cos(half_diff)
    z2 = math.sin(h) * math.sin(half_diff)
    n1 = math.cos(h) * math.cos(half_sum)
    n2 = math.cos(h) * math.sin(half_sum)
    s_plus = radius * math.atan2(z1, n1)  # (a + b)/2
    s_minus = radius * math.atan2(z2, n2)  # (a − b)/2
    a, b = s_plus + s_minus, s_plus - s_minus
    gamma = 2 * math.atan(math.sqrt((z1**2 + n1**2) / (z2**2 + n2**2)))
    if trace is not None:
        trace.section("Napier's analogies (AKA)")
        trace.step("z₁", z1, "sin(c/2R) cos((α − β)/2)")
        trace.step("z₂", z2, "sin(c/2R) sin((α − β)/2)")
        trace.step("n₁", n1, "cos(c/2R) cos((α + β)/2)")
        trace.step("n₂", n2, "cos(c/2R) sin((α + β)/2)")
        trace.step("(a + b)/2", s_plus, "R atan(z₁/n₁)", unit="m")
        trace.step("(a − b)/2", s_minus, "R atan(z₂/n₂)", unit="m")
        trace.step("a", a, unit="m")
        trace.step("b", b, unit="m")
        trace.step("γ", gamma, "2 atan √((z₁² + n₁²)/(z₂² + n₂²))", unit="rad")
    return _spherical(a, b, c, alpha, beta, gamma, radius)


def spherical_sas(
    a: float, b: float, gamma: float, radius: float, trace: Trace | None = None
) -> Triangle:
    """Spherical triangle from sides a, b and the included angle γ (Napier's analogies)."""
    half_g = gamma / 2
    d, s = (a - b) / (2 * radius), (a + b) / (2 * radius)
    z1 = math.cos(half_g) * math.cos(d)
    z2 = math.cos(half_g) * math.sin(d)
    n1 = math.sin(half_g) * math.cos(s)
    n2 = math.sin(half_g) * math.sin(s)
    half_sum = math.atan2(z1, n1)  # (α + β)/2
    half_diff = math.atan2(z2, n2)  # (α − β)/2
    alpha, beta = half_sum + half_diff, half_sum - half_diff
    c = 2 * radius * math.atan(math.sqrt((z2**2 + n2**2) / (z1**2 + n1**2)))
    if trace is not None:
        trace.section("Napier's analogies (KAK)")
        trace.step("z₁", z1, "cos(γ/2) cos((a − b)/2R)")
        trace.step("z₂", z2, "cos(γ/2) sin((a − b)/2R)")
        trace.step("n₁", n1, "sin(γ/2) cos((a + b)/2R)")
        trace.step("n₂", n2, "sin(γ/2) sin((a + b)/2R)")
        trace.step("(α + β)/2", half_sum, "atan(z₁/n₁)", unit="rad")
        trace.step("(α − β)/2", half_diff, "atan(z₂/n₂)", unit="rad")
        trace.step("c", c, "2R atan √((z₂² + n₂²)/(z₁² + n₁²))", unit="m")
    return _spherical(a, b, c, alpha, beta, gamma, radius)


def spherical_sss(
    a: float, b: float, c: float, radius: float, trace: Trace | None = None
) -> Triangle:
    """Spherical triangle from its three sides (half-angle formulas)."""
    _check_sides(a, b, c)
    u = (a + b + c) / 2

    def half_angle(opposite: float, x: float, y: float) -> float:
        num = math.sin((u - x) / radius) * math.sin((u - y) / radius)
        den = math.sin(u / radius) * math.sin((u - opposite) / radius)
        return 2 * math.atan(math.sqrt(num / den))

    alpha, beta, gamma = half_angle(a, b, c), half_angle(b, a, c), half_angle(c, a, b)
    if trace is not None:
        trace.section("Half-angle formulas (KKK)")
        trace.step("U", u, "(a + b + c)/2", unit="m")
        trace.step("α", alpha, "2 atan √(sin((U−b)/R) sin((U−c)/R) / (sin(U/R) sin((U−a)/R)))", unit="rad")
        trace.step("β", beta, "2 atan √(sin((U−a)/R) sin((U−c)/R) / (sin(U/R) sin((U−b)/R)))", unit="rad")
        trace.step("γ", gamma, "2 atan √(sin((U−a)/R) sin((U−b)/R) / (sin(U/R) sin((U−c)/R)))", unit="rad")
    return _spherical(a, b, c, alpha, beta, gamma, radius)


# ── Small spherical triangles, approximate (classical geodesy) ──────────────


def spherical_excess(plane_area: float, radius: float) -> float:
    """ε ≈ F / R², with F the area of the (nearly equal) plane triangle."""
    return plane_area / radius**2


def legendre_asa(
    c: float,
    alpha: float,
    beta: float,
    radius: float,
    trace: Trace | None = None,
    refine: bool = False,
) -> Triangle:
    """Legendre's theorem: solve the spherical triangle as a plane one with the
    same sides and each angle reduced by ε/3.

    With ``refine`` the reduction is ε/3 + ε/60 · (m² − a²)/R² per angle,
    m² = (a² + b² + c²)/3 (iterated once, since a, b are unknown at first).
    """
    first = plane_asa(c, alpha, beta)  # for the area only
    eps = spherical_excess(first.area, radius)
    if trace is not None:
        trace.section("Spherical excess")
        trace.step("F", first.area, "c² sin α sin β / (2 sin(α + β))", unit="m²")
        trace.step("ε", eps, "F / R²", unit="rad")
    red = [eps / 3] * 3
    for _ in range(2 if refine else 1):
        a_p, b_p = alpha - red[0], beta - red[1]
        plane = plane_asa(c, a_p, b_p)
        if not refine:
            break
        m2 = (plane.a**2 + plane.b**2 + c**2) / 3
        red = [eps / 3 + eps / 60 * (m2 - x * x) / radius**2 for x in (plane.a, plane.b, c)]
        eps = spherical_excess(plane.area, radius)
    if trace is not None:
        trace.section("Reduced (plane) triangle")
        trace.step("α′", a_p, "α − ε/3", unit="rad")
        trace.step("β′", b_p, "β − ε/3", unit="rad")
        trace.step("γ′", plane.gamma, "180° − α′ − β′", unit="rad")
        trace.step("a", plane.a, "c sin α′ / sin γ′", unit="m")
        trace.step("b", plane.b, "c sin β′ / sin γ′", unit="m")
        trace.step("γ", plane.gamma + red[2], "γ′ + ε/3", unit="rad")
    gamma = plane.gamma + red[2]
    return Triangle(plane.a, plane.b, c, alpha, beta, gamma, alpha + beta + gamma - math.pi, plane.area)


def additament_asa(
    c: float, alpha: float, beta: float, radius: float, trace: Trace | None = None
) -> Triangle:
    """Soldner's additament method: keep the spherical angles and shorten each
    side by its additament s³/(6R²), solve the plane sine rule, then add the
    additaments back."""
    first = plane_asa(c, alpha, beta)
    eps = spherical_excess(first.area, radius)
    gamma = math.pi + eps - alpha - beta
    c_red = c - c**3 / (6 * radius**2)
    a_red = c_red * math.sin(alpha) / math.sin(gamma)
    b_red = c_red * math.sin(beta) / math.sin(gamma)
    a = a_red + a_red**3 / (6 * radius**2)
    b = b_red + b_red**3 / (6 * radius**2)
    if trace is not None:
        trace.section("Additament method")
        trace.step("ε", eps, "F / R²", unit="rad")
        trace.step("γ", gamma, "180° + ε − α − β", unit="rad")
        trace.step("c′", c_red, "c − c³/(6R²)", unit="m")
        trace.step("a′", a_red, "c′ sin α / sin γ", unit="m")
        trace.step("b′", b_red, "c′ sin β / sin γ", unit="m")
        trace.step("a", a, "a′ + a′³/(6R²)", unit="m")
        trace.step("b", b, "b′ + b′³/(6R²)", unit="m")
    return Triangle(a, b, c, alpha, beta, gamma, eps, first.area)


def _check_sides(a: float, b: float, c: float) -> None:
    if min(a, b, c) <= 0:
        raise ValueError("sides must be positive")
    if a >= b + c or b >= a + c or c >= a + b:
        raise ValueError("each side must be shorter than the other two together")
