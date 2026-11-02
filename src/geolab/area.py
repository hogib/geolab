"""Areas on the ellipsoid bounded by meridians and parallels.

The area element is dA = MN cos φ dφ dλ = b² cos φ / (1 − e² sin²φ)² dφ dλ.
Between two parallels (a zone) the course's series is

    Z = 2πb² {α(sin φ₂ − sin φ₁) − β(sin 3φ₂ − sin 3φ₁) + γ(sin 5φ₂ − sin 5φ₁)
              − δ(sin 7φ₂ − sin 7φ₁) + ε(sin 9φ₂ − sin 9φ₁) − ζ(sin 11φ₂ − sin 11φ₁)}

and a quadrangle between two meridians takes the fraction Δλ / 2π of it.
:func:`zone_area_exact` gives the closed form for comparison.
"""

import math
from dataclasses import dataclass

from .ellipsoids import Ellipsoid
from .trace import Trace


@dataclass(frozen=True)
class AreaCoefficients:
    alpha: float
    beta: float
    gamma: float
    delta: float
    epsilon: float
    zeta: float


def area_coefficients(ellipsoid: Ellipsoid) -> AreaCoefficients:
    """α … ζ to e¹⁰, as in the course's Hayford table.

    (The slides print 35/198 for β's e⁸ term; expanding sin⁹φ gives 35/192.)
    """
    e2 = ellipsoid.e_sq
    e4, e6, e8, e10 = e2**2, e2**3, e2**4, e2**5
    return AreaCoefficients(
        1 + e2 / 2 + 3 / 8 * e4 + 5 / 16 * e6 + 35 / 128 * e8 + 63 / 256 * e10,
        e2 / 6 + 3 / 16 * e4 + 3 / 16 * e6 + 35 / 192 * e8 + 45 / 256 * e10,
        3 / 80 * e4 + e6 / 16 + 5 / 64 * e8 + 45 / 512 * e10,
        e6 / 112 + 5 / 256 * e8 + 15 / 512 * e10,
        5 / 2304 * e8 + 3 / 512 * e10,
        3 / 5632 * e10,
    )


def zone_area(lat1: float, lat2: float, ellipsoid: Ellipsoid, trace: Trace | None = None) -> float:
    """Area between the parallels φ₁ and φ₂ all the way round (m², ≥ 0)."""
    k = area_coefficients(ellipsoid)
    b = ellipsoid.b

    def diff(n: int) -> float:
        return math.sin(n * lat2) - math.sin(n * lat1)

    z = 2 * math.pi * b**2 * (
        k.alpha * diff(1)
        - k.beta * diff(3)
        + k.gamma * diff(5)
        - k.delta * diff(7)
        + k.epsilon * diff(9)
        - k.zeta * diff(11)
    )
    z = abs(z)
    if trace is not None:
        trace.section("Zone between two parallels")
        trace.step("α", k.alpha, "1 + e²/2 + 3e⁴/8 + 5e⁶/16 + 35e⁸/128 + 63e¹⁰/256")
        trace.step("β", k.beta, "e²/6 + 3e⁴/16 + 3e⁶/16 + 35e⁸/192 + 45e¹⁰/256")
        trace.step("γ", k.gamma, "3e⁴/80 + e⁶/16 + 5e⁸/64 + 45e¹⁰/512")
        trace.step("δ", k.delta, "e⁶/112 + 5e⁸/256 + 15e¹⁰/512")
        trace.step("ε", k.epsilon, "5e⁸/2304 + 3e¹⁰/512")
        trace.step("ζ", k.zeta, "3e¹⁰/5632")
        trace.step(
            "Z",
            z,
            "2πb² {α Δsin φ − β Δsin 3φ + γ Δsin 5φ − δ Δsin 7φ + ε Δsin 9φ − ζ Δsin 11φ}",
            unit="m²",
        )
    return z


def zone_area_exact(lat1: float, lat2: float, ellipsoid: Ellipsoid) -> float:
    """Closed form of the zone area: 2πb² [q(φ)/2] between φ₁ and φ₂, with
    q(φ) = sin φ/(1 − e² sin²φ) + atanh(e sin φ)/e."""
    e = math.sqrt(ellipsoid.e_sq)

    def q(lat: float) -> float:
        s = math.sin(lat)
        if e == 0:
            return 2 * s
        return s / (1 - e**2 * s**2) + math.atanh(e * s) / e

    return abs(math.pi * ellipsoid.b**2 * (q(lat2) - q(lat1)))


def quadrangle_area(
    lat1: float,
    lat2: float,
    lon1: float,
    lon2: float,
    ellipsoid: Ellipsoid,
    trace: Trace | None = None,
) -> float:
    """Area bounded by the parallels φ₁, φ₂ and the meridians λ₁, λ₂ (the shorter
    way round): P = (Δλ / 2π) Z (m²)."""
    dlon = abs(math.remainder(lon2 - lon1, math.tau))
    z = zone_area(lat1, lat2, ellipsoid, trace)
    p = dlon / math.tau * z
    if trace is not None:
        trace.section("Quadrangle")
        trace.step("Δλ", dlon, "|λ₂ − λ₁|", unit="rad")
        trace.step("P", p, "(Δλ / 2π) Z", unit="m²")
    return p


def ellipsoid_area(ellipsoid: Ellipsoid, trace: Trace | None = None) -> float:
    """Total surface area, 4πb² (α + β + γ + δ + ε + ζ), twice the equator-to-pole zone."""
    k = area_coefficients(ellipsoid)
    total = 4 * math.pi * ellipsoid.b**2 * (k.alpha + k.beta + k.gamma + k.delta + k.epsilon + k.zeta)
    if trace is not None:
        trace.step("E", total, "4πb² (α + β + γ + δ + ε + ζ)", unit="m²")
    return total
