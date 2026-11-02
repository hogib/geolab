"""Normal sections versus the geodesic (series in η² = e′² cos²φ).

A theodolite at A sighting B measures the azimuth of the normal section
through A's normal and B (α′). The normal section from B back to A gives a
different curve (α″), and the geodesic runs between them, splitting the gap
1 : 2 (nearer A's section). For a line of length S from A (latitude φ_A,
prime-vertical radius N_A):

    α − α′  = −ρ/12 η² sin 2α′ (S/N)²         geodesic − normal section
    α′ − α″ =  ρ/4  η² sin 2α′ (S/N)²         normal − reverse normal section
    d_max   =  1/32 η² sin 2α′ S³/N²           largest gap, near the middle
    s − s′  = −η⁴/(360 N⁴) S⁵ sin² 2α′         geodesic − normal-section length

A target at height h₂ bends the observed direction γ off the normal section:

    α′ − γ  =  ρ η² h₂ sin 2α′ / (2N)

Angles here are returned in radians (the ρ factor is only for arc seconds).
"""

import math
from dataclasses import dataclass

from .curvature import prime_vertical_radius
from .ellipsoids import Ellipsoid
from .trace import Trace


@dataclass(frozen=True)
class NormalSectionCorrections:
    eta_sq: float
    geodesic_minus_normal: float
    """α − α′ (rad)."""
    normal_minus_reverse: float
    """α′ − α″ (rad)."""
    max_separation: float
    """d_max (m)."""
    length_difference: float
    """s − s′ (m), always ≤ 0."""


def eta_squared(lat: float, ellipsoid: Ellipsoid) -> float:
    """η² = e′² cos²φ."""
    return ellipsoid.ep_sq * math.cos(lat) ** 2


def normal_section_corrections(
    lat: float,
    azimuth: float,
    distance: float,
    ellipsoid: Ellipsoid,
    trace: Trace | None = None,
) -> NormalSectionCorrections:
    eta2 = eta_squared(lat, ellipsoid)
    n = prime_vertical_radius(lat, ellipsoid)
    s2a = math.sin(2 * azimuth)
    ratio = distance / n
    geo = -eta2 / 12 * s2a * ratio**2
    rev = eta2 / 4 * s2a * ratio**2
    dmax = eta2 / 32 * s2a * distance**3 / n**2
    ds = -(eta2**2) / (360 * n**4) * distance**5 * s2a**2
    if trace is not None:
        trace.section("Normal section vs geodesic")
        trace.step("η²", eta2, "e′² cos²φ_A")
        trace.step("N_A", n, unit="m")
        trace.step("α − α′", geo, "−η²/12 sin 2α′ (S/N)²", unit="rad")
        trace.step("α′ − α″", rev, "η²/4 sin 2α′ (S/N)²", unit="rad")
        trace.step("d_max", dmax, "η²/32 sin 2α′ S³/N²", unit="m")
        trace.step("s − s′", ds, "−η⁴/(360 N⁴) S⁵ sin² 2α′", unit="m")
    return NormalSectionCorrections(eta2, geo, rev, dmax, ds)


def target_height_correction(
    lat: float,
    azimuth: float,
    target_height: float,
    ellipsoid: Ellipsoid,
    trace: Trace | None = None,
) -> float:
    """α′ − γ = η² h₂ sin 2α′ / (2N) (rad): add to the observed direction γ to
    reduce it to the normal section on the ellipsoid. Independent of distance."""
    eta2 = eta_squared(lat, ellipsoid)
    n = prime_vertical_radius(lat, ellipsoid)
    corr = eta2 * target_height * math.sin(2 * azimuth) / (2 * n)
    if trace is not None:
        trace.section("Target height reduction")
        trace.step("α′ − γ", corr, "η² h₂ sin 2α′ / (2N)", unit="rad")
    return corr
