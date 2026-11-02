"""Soldner (geodetic rectangular) coordinates on the ellipsoid.

A central meridian λ₀ is the x axis, with origin P₀(φ₀, λ₀). For a point P,
the geodesic through P perpendicular to the central meridian meets it at the
foot point P_f: y is the length PP_f (positive east) and x the meridian arc
P₀P_f (positive north). γ is the meridian convergence at P, the angle between
the meridian and the geodetic parallel (the line y = const); the grid bearing
T of a line relates to its azimuth by α = γ + T.

Forward (course notation, Δλ = λ − λ₀, at the point's latitude φ):

    x = x_B + (02)Δλ² + (04)Δλ⁴           x_B = G(φ) − G(φ₀)
    y = (01)Δλ + (03)Δλ³ + (05)Δλ⁵
    γ = (01)′Δλ + (03)′Δλ³

Inverse, at the foot-point latitude φ_f from G(φ_f) = G(φ₀) + x:

    φ = φ_f + [2]y² + [4]y⁴     Δλ = [1]y + [3]y³ + [5]y⁵     γ = [1]′y + [3]′y³

Lengths are in metres and angles in radians (ρ = 1).
"""

import math
from dataclasses import dataclass

from .arcs import latitude_from_meridian_distance, meridian_distance
from .curvature import meridian_radius, prime_vertical_radius
from .ellipsoids import Ellipsoid
from .trace import Trace


@dataclass(frozen=True)
class SoldnerPoint:
    x: float
    y: float
    convergence: float
    """Meridian convergence γ (rad)."""


def _aux(lat: float, ellipsoid: Ellipsoid) -> tuple[float, float, float, float]:
    """N, t = tan φ, η² and V² at φ."""
    n = prime_vertical_radius(lat, ellipsoid)
    t = math.tan(lat)
    eta2 = ellipsoid.ep_sq * math.cos(lat) ** 2
    return n, t, eta2, 1 + eta2


def soldner_from_geodetic(
    lat: float,
    lon: float,
    lat0: float,
    lon0: float,
    ellipsoid: Ellipsoid,
    trace: Trace | None = None,
) -> SoldnerPoint:
    """x, y, γ of the point φ, λ in the Soldner system with origin φ₀, λ₀."""
    n, t, eta2, _ = _aux(lat, ellipsoid)
    s, c = math.sin(lat), math.cos(lat)
    dl = math.remainder(lon - lon0, math.tau)
    k01 = n * c
    k02 = n * c**2 * t / 2
    k03 = -n * c**3 * t**2 / 6
    k04 = n * s * c**3 * (5 - t**2 + 5 * eta2) / 24
    k05 = -n * s**2 * c**3 * (8 - t**2) / 120
    g01 = s
    g03 = (1 + eta2) / 3 * s * c**2
    x_b = meridian_distance(lat, ellipsoid) - meridian_distance(lat0, ellipsoid)
    x = x_b + k02 * dl**2 + k04 * dl**4
    y = k01 * dl + k03 * dl**3 + k05 * dl**5
    gamma = g01 * dl + g03 * dl**3
    if trace is not None:
        trace.section("Setup")
        trace.step("Δλ", dl, "λ − λ₀", unit="rad")
        trace.step("x_B", x_b, "G(φ) − G(φ₀)", unit="m")
        trace.section("Coefficients")
        trace.step("(01)", k01, "N cos φ", unit="m")
        trace.step("(02)", k02, "N cos²φ t / 2", unit="m")
        trace.step("(03)", k03, "−N cos³φ t² / 6", unit="m")
        trace.step("(04)", k04, "N sin φ cos³φ (5 − t² + 5η²) / 24", unit="m")
        trace.step("(05)", k05, "−N sin²φ cos³φ (8 − t²) / 120", unit="m")
        trace.step("(01)′", g01, "sin φ")
        trace.step("(03)′", g03, "(1 + η²)/3 sin φ cos²φ")
        trace.section("Result")
        trace.step("x", x, "x_B + (02)Δλ² + (04)Δλ⁴", unit="m")
        trace.step("y", y, "(01)Δλ + (03)Δλ³ + (05)Δλ⁵", unit="m")
        trace.step("γ", gamma, "(01)′Δλ + (03)′Δλ³", unit="rad")
    return SoldnerPoint(x, y, gamma)


@dataclass(frozen=True)
class SoldnerGeodetic:
    lat: float
    lon: float
    convergence: float
    foot_lat: float
    """Foot-point latitude φ_f."""


def geodetic_from_soldner(
    x: float,
    y: float,
    lat0: float,
    lon0: float,
    ellipsoid: Ellipsoid,
    trace: Trace | None = None,
) -> SoldnerGeodetic:
    """φ, λ, γ of the Soldner point x, y (origin φ₀, λ₀)."""
    g0 = meridian_distance(lat0, ellipsoid)
    lat_f = latitude_from_meridian_distance(g0 + x, ellipsoid).lat
    n, t, eta2, v2 = _aux(lat_f, ellipsoid)
    c = math.cos(lat_f)
    k1 = 1 / (n * c)
    k2 = -v2 * t / (2 * n**2)
    k3 = -(t**2) / (3 * n**3 * c)
    k4 = v2 * t / (24 * n**4) * (1 + 3 * t**2 + eta2 - 9 * eta2 * t**2)
    k5 = t**2 / (15 * n**5 * c) * (1 + 3 * t**2)
    g1 = t / n
    g3 = -t / (6 * n**3) * (1 + 2 * t**2 + eta2)
    lat = lat_f + k2 * y**2 + k4 * y**4
    dl = k1 * y + k3 * y**3 + k5 * y**5
    gamma = g1 * y + g3 * y**3
    lon = math.remainder(lon0 + dl, math.tau)
    if trace is not None:
        trace.section("Foot point")
        trace.step("G₀", g0, "meridian arc to φ₀", unit="m")
        trace.step("φ_f", lat_f, "latitude of the arc G₀ + x", unit="rad")
        trace.section("Coefficients")
        trace.step("[1]", k1, "1 / (N_f cos φ_f)")
        trace.step("[2]", k2, "−V_f² t_f / (2N_f²)")
        trace.step("[3]", k3, "−t_f² / (3N_f³ cos φ_f)")
        trace.step("[4]", k4, "V_f² t_f (1 + 3t² + η² − 9η²t²) / (24N_f⁴)")
        trace.step("[5]", k5, "t_f² (1 + 3t_f²) / (15N_f⁵ cos φ_f)")
        trace.step("[1]′", g1, "t_f / N_f")
        trace.step("[3]′", g3, "−t_f (1 + 2t_f² + η_f²) / (6N_f³)")
        trace.section("Result")
        trace.step("φ", lat, "φ_f + [2]y² + [4]y⁴", unit="rad")
        trace.step("Δλ", dl, "[1]y + [3]y³ + [5]y⁵", unit="rad")
        trace.step("λ", lon, "λ₀ + Δλ", unit="rad")
        trace.step("γ", gamma, "[1]′y + [3]′y³", unit="rad")
    return SoldnerGeodetic(lat, lon, gamma, lat_f)


def _gauss_radius(lat: float, ellipsoid: Ellipsoid) -> float:
    return math.sqrt(meridian_radius(lat, ellipsoid) * prime_vertical_radius(lat, ellipsoid))


@dataclass(frozen=True)
class SoldnerDirect:
    x2: float
    y2: float
    bearing2: float
    """T₂₁, the grid bearing from P₂ back to P₁."""


def soldner_direct(
    x1: float,
    y1: float,
    bearing: float,
    distance: float,
    lat0: float,
    ellipsoid: Ellipsoid,
    trace: Trace | None = None,
) -> SoldnerDirect:
    """First problem in Soldner coordinates: x₂, y₂, T₂₁ from x₁, y₁, T₁₂ and S.

    Uses φ₁ at P₁ (from G₀ + x₁, near enough for t and η²) and R₁ = √(M₁N₁).
    """
    lat1 = latitude_from_meridian_distance(meridian_distance(lat0, ellipsoid) + x1, ellipsoid).lat
    r = _gauss_radius(lat1, ellipsoid)
    t = math.tan(lat1)
    eta2 = ellipsoid.ep_sq * math.cos(lat1) ** 2
    u, v = distance * math.cos(bearing), distance * math.sin(bearing)
    dy = (
        v
        - u**2 / (2 * r**2) * (y1 + v / 3)
        + 2 * eta2 * t * y1 * u**3 / (3 * r**3)
        - y1**2 * u**2 / (6 * r**4) * (y1 + 3 * v)
    )
    y2 = y1 + dy
    dx = (
        u
        + u / (2 * r**2) * (y2**2 - v**2 / 3)
        - eta2 * t * y2**2 * u**2 / r**3
        + y2**2 * u / (24 * r**4) * (5 * y2**2 + 4 * u**2 - 2 * v**2)
    )
    x2 = x1 + dx
    dt = (
        -u / r**2 * (y1 + v / 2)
        + 2 * eta2 * t * y1 * u**2 / r**3
        - y1**2 * u / (2 * r**4) * (y1 + 3 * v)
    )
    bearing2 = (bearing + math.pi + dt) % math.tau
    if trace is not None:
        trace.section("Setup")
        trace.step("φ₁", lat1, "from G₀ + x₁", unit="rad")
        trace.step("R₁", r, "√(M₁N₁)", unit="m")
        trace.step("u", u, "S cos T₁", unit="m")
        trace.step("v", v, "S sin T₁", unit="m")
        trace.section("Result")
        trace.step("y₂ − y₁", dy, "v − u²/(2R²)(y₁ + v/3) + 2η²t y₁u³/(3R³) − y₁²u²/(6R⁴)(y₁ + 3v)", unit="m")
        trace.step("x₂ − x₁", dx, "u + u/(2R²)(y₂² − v²/3) − η²t y₂²u²/R³ + y₂²u/(24R⁴)(5y₂² + 4u² − 2v²)", unit="m")
        trace.step("T₂ − T₁ − 180°", dt, "−u/R²(y₁ + v/2) + 2η²t y₁u²/R³ − y₁²u/(2R⁴)(y₁ + 3v)", unit="rad")
        trace.step("x₂", x2, unit="m")
        trace.step("y₂", y2, unit="m")
        trace.step("T₂₁", bearing2, "T₁₂ + 180° + ΔT", unit="rad")
    return SoldnerDirect(x2, y2, bearing2)


@dataclass(frozen=True)
class SoldnerInverse:
    distance: float
    bearing1: float
    """T₁₂."""
    bearing2: float
    """T₂₁."""


def soldner_inverse(
    x1: float,
    y1: float,
    x2: float,
    y2: float,
    lat0: float,
    ellipsoid: Ellipsoid,
    trace: Trace | None = None,
) -> SoldnerInverse:
    """Second problem in Soldner coordinates (mid-point formulas): S, T₁₂, T₂₁."""
    xm, ym = (x1 + x2) / 2, (y1 + y2) / 2
    dx, dy = x2 - x1, y2 - y1
    lat_m = latitude_from_meridian_distance(meridian_distance(lat0, ellipsoid) + xm, ellipsoid).lat
    r = _gauss_radius(lat_m, ellipsoid)
    t = math.tan(lat_m)
    eta2 = ellipsoid.ep_sq * math.cos(lat_m) ** 2
    s_sin = dy - dx**2 * dy / (12 * r**2) + eta2 * t * ym * dx**3 / (3 * r**3) - ym**2 * dx**2 * dy / (24 * r**4)
    s_cos = dx - dx / (2 * r**2) * (ym**2 - dy**2 / 12) + ym**2 * dx / (24 * r**4) * (ym**2 - dx**2 - dy**2 / 2)
    dt = -dx / r**2 * (ym + ym**3 / (2 * r**2))
    tm = math.atan2(s_sin, s_cos)
    s = math.hypot(s_sin, s_cos)
    t1 = (tm - dt / 2) % math.tau
    t2 = (tm + dt / 2 + math.pi) % math.tau
    if trace is not None:
        trace.section("Setup")
        trace.step("Δx", dx, "x₂ − x₁", unit="m")
        trace.step("Δy", dy, "y₂ − y₁", unit="m")
        trace.step("y_m", ym, "(y₁ + y₂)/2", unit="m")
        trace.step("φ_m", lat_m, "from G₀ + x_m", unit="rad")
        trace.step("R_m", r, "√(M_m N_m)", unit="m")
        trace.section("Result")
        trace.step("S sin T_m", s_sin, "Δy − Δx²Δy/(12R²) + η²t y_m Δx³/(3R³) − y_m²Δx²Δy/(24R⁴)", unit="m")
        trace.step("S cos T_m", s_cos, "Δx − Δx/(2R²)(y_m² − Δy²/12) + y_m²Δx/(24R⁴)(y_m² − Δx² − Δy²/2)", unit="m")
        trace.step("ΔT", dt, "−Δx/R² (y_m + y_m³/(2R²))", unit="rad")
        trace.step("S", s, "√((S sin T_m)² + (S cos T_m)²)", unit="m")
        trace.step("T₁₂", t1, "T_m − ΔT/2", unit="rad")
        trace.step("T₂₁", t2, "T_m + ΔT/2 ± 180°", unit="rad")
    return SoldnerInverse(s, t1, t2)
