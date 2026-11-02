"""Arc length on a planar ellipse x²/a² + y²/b² = 1.

Points are given by the parametric angle t, (x, y) = (a cos t, b sin t), or
by the polar angle θ from the centre. Arc lengths are incomplete elliptic
integrals of the second kind, evaluated with Carlson's symmetric forms R_F
and R_D (duplication algorithm), so any a, b > 0 works, however eccentric.
"""

import math
from dataclasses import dataclass

from .trace import Trace


@dataclass(frozen=True)
class EllipseArc:
    length: float
    """Arc length from point 1 to point 2, counter-clockwise when t₂ > t₁ (always ≥ 0)."""
    chord: float
    p1: tuple[float, float]
    p2: tuple[float, float]
    perimeter: float


def carlson_rf(x: float, y: float, z: float) -> float:
    """Carlson's R_F(x, y, z) = ½ ∫₀^∞ dt / √((t + x)(t + y)(t + z))."""
    while True:
        sx, sy, sz = math.sqrt(x), math.sqrt(y), math.sqrt(z)
        lam = sx * sy + sy * sz + sz * sx
        x, y, z = (x + lam) / 4, (y + lam) / 4, (z + lam) / 4
        avg = (x + y + z) / 3
        dx, dy, dz = 1 - x / avg, 1 - y / avg, 1 - z / avg
        if max(abs(dx), abs(dy), abs(dz)) < 1e-4:
            break
    e2 = dx * dy - dz**2
    e3 = dx * dy * dz
    return (1 - e2 / 10 + e3 / 14 + e2**2 / 24 - 3 * e2 * e3 / 44) / math.sqrt(avg)


def carlson_rd(x: float, y: float, z: float) -> float:
    """Carlson's R_D(x, y, z) = 3/2 ∫₀^∞ dt / ((t + z)√((t + x)(t + y)(t + z)))."""
    total, fac = 0.0, 1.0
    while True:
        sx, sy, sz = math.sqrt(x), math.sqrt(y), math.sqrt(z)
        lam = sx * sy + sy * sz + sz * sx
        total += fac / (sz * (z + lam))
        fac /= 4
        x, y, z = (x + lam) / 4, (y + lam) / 4, (z + lam) / 4
        avg = (x + y + 3 * z) / 5
        dx, dy, dz = 1 - x / avg, 1 - y / avg, 1 - z / avg
        if max(abs(dx), abs(dy), abs(dz)) < 1e-4:
            break
    ea = dx * dy
    eb = dz**2
    ec = ea - eb
    ed = ea - 6 * eb
    ee = ed + 2 * ec
    series = 1 + ed * (-3 / 14 + 9 / 88 * ed - 9 / 52 * dz * ee) + dz * (
        ee / 6 + dz * (-9 / 22 * ec + 3 / 26 * dz * ea)
    )
    return 3 * total + fac * series / (avg * math.sqrt(avg))


def elliptic_e(phi: float, m: float) -> float:
    """Incomplete elliptic integral of the second kind, E(φ | m) = ∫₀^φ √(1 − m sin²u) du.

    Defined for any real φ and m ≤ 1 (m = k², the parameter).
    """
    j = round(phi / math.pi)
    phi -= j * math.pi
    s, c = math.sin(phi), math.cos(phi)
    if m == 1:  # degenerate: the integrand is |cos u|, and R_F(0, 0, 1) diverges
        return s + 2 * j
    q = 1 - m * s**2
    e = s * carlson_rf(c**2, q, 1) - m * s**3 / 3 * carlson_rd(c**2, q, 1)
    if j:
        e += 2 * j * complete_elliptic_e(m)
    return e


def complete_elliptic_e(m: float) -> float:
    """Complete elliptic integral of the second kind, E(m) = E(π/2 | m)."""
    if m == 1:
        return 1.0
    return carlson_rf(0, 1 - m, 1) - m / 3 * carlson_rd(0, 1 - m, 1)


def parametric_from_polar(theta: float, a: float, b: float) -> float:
    """Parametric angle t of the point at polar angle θ, on the same turn as θ."""
    t = math.atan2(a * math.sin(theta), b * math.cos(theta))
    return t + math.tau * round((theta - t) / math.tau)


def parametric_from_normal(phi: float, a: float, b: float) -> float:
    """Parametric angle t of the point whose normal makes angle φ with the x axis.

    On a meridian ellipse φ is the geodetic latitude and t the reduced latitude β:
    tan β = (b/a) tan φ.
    """
    t = math.atan2(b * math.sin(phi), a * math.cos(phi))
    return t + math.tau * round((phi - t) / math.tau)


def ellipse_arc_from_vertex(t: float, a: float, b: float) -> float:
    """Signed arc length from (a, 0) to the point at parametric angle t."""
    _check_axes(a, b)
    if a >= b:
        # √(a² sin²t + b² cos²t) = a√(1 − m sin²(π/2 − t)), m = 1 − b²/a²
        m = 1 - (b / a) ** 2
        return a * (complete_elliptic_e(m) - elliptic_e(math.pi / 2 - t, m))
    # √(b² − (b² − a²) sin²t) = b√(1 − m sin²t), m = 1 − a²/b²
    m = 1 - (a / b) ** 2
    return b * elliptic_e(t, m)


def ellipse_perimeter(a: float, b: float) -> float:
    """Perimeter 4 a E(1 − b²/a²) (with a the larger semi-axis)."""
    _check_axes(a, b)
    big, small = max(a, b), min(a, b)
    return 4 * big * complete_elliptic_e(1 - (small / big) ** 2)


def ellipse_arc(
    a: float, b: float, t1: float, t2: float, trace: Trace | None = None
) -> EllipseArc:
    """Arc between parametric angles t₁ and t₂, travelling from t₁ to t₂."""
    _check_axes(a, b)
    big, small = max(a, b), min(a, b)
    m = 1 - (small / big) ** 2
    s1 = ellipse_arc_from_vertex(t1, a, b)
    s2 = ellipse_arc_from_vertex(t2, a, b)
    length = abs(s2 - s1)
    p1 = (a * math.cos(t1), b * math.sin(t1))
    p2 = (a * math.cos(t2), b * math.sin(t2))
    chord = math.dist(p1, p2)
    perimeter = ellipse_perimeter(a, b)

    if trace is not None:
        trace.section("Setup")
        trace.step("m", m, "k² = 1 − (minor/major)²")
        trace.step("E(m)", complete_elliptic_e(m), "R_F(0, 1 − m, 1) − m/3 R_D(0, 1 − m, 1)")
        if a >= b:
            formula = "a (E(m) − E(π/2 − t | m))"
        else:
            formula = "b E(t | m)"
        trace.section("Arc from the vertex (a, 0)")
        trace.step("s(t₁)", s1, formula, unit="m")
        trace.step("s(t₂)", s2, formula, unit="m")
        trace.section("Result")
        trace.step("s", length, "|s(t₂) − s(t₁)|", unit="m")
        trace.step("x₁", p1[0], "a cos t₁", unit="m")
        trace.step("y₁", p1[1], "b sin t₁", unit="m")
        trace.step("x₂", p2[0], "a cos t₂", unit="m")
        trace.step("y₂", p2[1], "b sin t₂", unit="m")
        trace.step("c", chord, "√((x₂ − x₁)² + (y₂ − y₁)²)", unit="m")
        trace.step("P", perimeter, "4 a E(m)", unit="m")

    return EllipseArc(length, chord, p1, p2, perimeter)


def _check_axes(a: float, b: float) -> None:
    if not (a > 0 and b > 0):
        raise ValueError("semi-axes a and b must be positive")
