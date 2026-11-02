"""Differential geometry of the ellipsoid with Gauss parameters u = φ, v = λ.

With c = a²/b and v = √(1 + e′² cos²φ) (so M = c/v³ and N = c/v):

    first fundamental form   ds² = E dφ² + 2F dφ dλ + G dλ²
        E = c²/v⁶ = M²       F = 0       G = c² cos²φ / v² = N² cos²φ
    second fundamental form  II = L dφ² + 2M̄ dφ dλ + N̄ dλ²
        L = c/v³ = M         M̄ = 0       N̄ = c cos²φ / v = N cos²φ

F = 0 means meridians and parallels cross at right angles. The normal
curvature in direction α (Euler) is 1/R = cos²α/M + sin²α/N; its extremes are
M (α = 0°) and N (α = 90°), the Gauss curvature is K = 1/(MN) and the mean
curvature H = ½(1/M + 1/N).
"""

import math
from dataclasses import dataclass

from .ellipsoids import Ellipsoid
from .trace import Trace


@dataclass(frozen=True)
class FundamentalForms:
    E: float
    F: float
    G: float
    L: float
    M: float
    N: float
    """Second-form coefficient N̄ (not the prime-vertical radius)."""
    gauss_curvature: float
    """K = (LN̄ − M̄²) / (EG − F²) = 1/(MN), in 1/m²."""
    mean_curvature: float
    """H = ½(1/M + 1/N), in 1/m."""
    gauss_radius: float
    """R_G = 1/√K = √(MN) = c/v², m."""


def _c_v(lat: float, ellipsoid: Ellipsoid) -> tuple[float, float]:
    c = ellipsoid.a**2 / ellipsoid.b
    v = math.sqrt(1 + ellipsoid.ep_sq * math.cos(lat) ** 2)
    return c, v


def fundamental_forms(
    lat: float, ellipsoid: Ellipsoid, trace: Trace | None = None
) -> FundamentalForms:
    """First and second fundamental quantities and the curvatures at latitude φ."""
    c, v = _c_v(lat, ellipsoid)
    cos2 = math.cos(lat) ** 2
    e, f, g = c**2 / v**6, 0.0, c**2 * cos2 / v**2
    l, m, n = c / v**3, 0.0, c * cos2 / v
    radius_m, radius_n = c / v**3, c / v
    k = 1 / (radius_m * radius_n)
    h = (1 / radius_m + 1 / radius_n) / 2
    r_g = c / v**2
    if trace is not None:
        trace.section("Auxiliary quantities")
        trace.step("c", c, "a² / b", unit="m")
        trace.step("v", v, "√(1 + e′² cos²φ)")
        trace.section("First fundamental quantities")
        trace.step("E", e, "c² / v⁶ = M²", unit="m²")
        trace.step("F", f, "0 (meridians ⟂ parallels)", unit="m²")
        trace.step("G", g, "c² cos²φ / v² = N² cos²φ", unit="m²")
        trace.section("Second fundamental quantities")
        trace.step("L", l, "c / v³ = M", unit="m")
        trace.step("M̄", m, "0", unit="m")
        trace.step("N̄", n, "c cos²φ / v = N cos²φ", unit="m")
        trace.section("Curvatures")
        trace.step("K", k, "(LN̄ − M̄²)/(EG − F²) = 1/(MN)")
        trace.step("H", h, "½ (1/M + 1/N)")
        trace.step("R_G", r_g, "1/√K = √(MN) = c / v²", unit="m")
    return FundamentalForms(e, f, g, l, m, n, k, h, r_g)


def normal_section_radius(lat: float, azimuth: float, ellipsoid: Ellipsoid) -> float:
    """Normal curvature radius in azimuth α, written as on the slides:
    R = c / (v (v² cos²α + sin²α))."""
    c, v = _c_v(lat, ellipsoid)
    return c / (v * (v**2 * math.cos(azimuth) ** 2 + math.sin(azimuth) ** 2))


def curve_azimuth(lat: float, dlat: float, dlon: float, ellipsoid: Ellipsoid) -> float:
    """Azimuth α of a surface curve element (dφ, dλ), measured from the meridian:
    tan α = √G dλ / (√E dφ) = v² cos φ dλ/dφ. Returned in [0, 2π)."""
    _, v = _c_v(lat, ellipsoid)
    return math.atan2(v**2 * math.cos(lat) * dlon, dlat) % math.tau


def arc_element(lat: float, dlat: float, dlon: float, ellipsoid: Ellipsoid) -> float:
    """Length of a small curve element, ds = √(E dφ² + G dλ²) (m)."""
    ff = fundamental_forms(lat, ellipsoid)
    return math.sqrt(ff.E * dlat**2 + 2 * ff.F * dlat * dlon + ff.G * dlon**2)


def area_element(lat: float, dlat: float, dlon: float, ellipsoid: Ellipsoid) -> float:
    """Area of a small φ, λ cell, dA = √(EG − F²) dφ dλ = MN cos φ dφ dλ (m²)."""
    ff = fundamental_forms(lat, ellipsoid)
    return math.sqrt(ff.E * ff.G - ff.F**2) * abs(dlat * dlon)


def meusnier_radius(normal_radius: float, angle: float) -> float:
    """Meusnier: a curve whose principal normal makes angle δ with the surface
    normal has curvature radius R = R_nor cos δ."""
    return normal_radius * math.cos(angle)
