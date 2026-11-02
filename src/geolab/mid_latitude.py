"""Direct and inverse problems for short lines by Gauss's mid-latitude formulas.

Expanding Legendre's series about the middle of the line makes all odd-order
terms cancel, so the formulas are exact to 4th order. With the mean values
φ_m = (φ₁ + φ₂)/2, α_m = (α₁ + α₂ ± 180°)/2 and, at φ_m, M_m, N_m, t, η², V²,
the course's bracket abbreviations (with ρ = 1, angles in radians) are

    [3] = 1/24                        [4] = (1 + η² − 9η²t²) / (24V⁴)
    [5] = (1 − 2η²)/24                [6] = η²(1 − t²) / (8V⁴)
    [7] = (1 + η²)/12                 [8] = (3 + 8η²) / (24V⁴)

    s cos α_m = M_m Δφ cos(Δλ/2) {1 + [5](Δλ cos φ_m)² + [6]Δφ²}
    s sin α_m = N_m Δλ cos φ_m   {1 − [3](Δλ sin φ_m)² + [4]Δφ²}
    Δα = α₂ − α₁ ± 180° = Δλ sin φ_m {1 + [7](Δλ cos φ_m)² + [8]Δφ²}

The inverse problem follows directly; the direct problem needs iteration,
since φ_m and α_m depend on the unknown end point.
"""

import math
from dataclasses import dataclass

from .curvature import meridian_radius, prime_vertical_radius
from .ellipsoids import Ellipsoid
from .errors import ConvergenceError
from .trace import Trace


@dataclass(frozen=True)
class Brackets:
    m: float
    """M_m (m)."""
    n: float
    """N_m (m)."""
    b3: float
    b4: float
    b5: float
    b6: float
    b7: float
    b8: float


def brackets(lat_m: float, ellipsoid: Ellipsoid) -> Brackets:
    t2 = math.tan(lat_m) ** 2
    eta2 = ellipsoid.ep_sq * math.cos(lat_m) ** 2
    v4 = (1 + eta2) ** 2
    return Brackets(
        meridian_radius(lat_m, ellipsoid),
        prime_vertical_radius(lat_m, ellipsoid),
        1 / 24,
        (1 + eta2 - 9 * eta2 * t2) / (24 * v4),
        (1 - 2 * eta2) / 24,
        eta2 * (1 - t2) / (8 * v4),
        (1 + eta2) / 12,
        (3 + 8 * eta2) / (24 * v4),
    )


def _trace_brackets(trace: Trace, k: Brackets) -> None:
    trace.step("M_m", k.m, unit="m")
    trace.step("N_m", k.n, unit="m")
    trace.step("[3]", k.b3, "1/24")
    trace.step("[4]", k.b4, "(1 + η² − 9η²t²)/(24V⁴)")
    trace.step("[5]", k.b5, "(1 − 2η²)/24")
    trace.step("[6]", k.b6, "η²(1 − t²)/(8V⁴)")
    trace.step("[7]", k.b7, "(1 + η²)/12")
    trace.step("[8]", k.b8, "(3 + 8η²)/(24V⁴)")


@dataclass(frozen=True)
class MidLatitudeInverse:
    distance: float
    azimuth: float
    """α₁₂."""
    back_azimuth: float
    """α₂₁."""
    mean_azimuth: float
    delta_azimuth: float


def mid_latitude_inverse(
    lat1: float,
    lon1: float,
    lat2: float,
    lon2: float,
    ellipsoid: Ellipsoid,
    trace: Trace | None = None,
) -> MidLatitudeInverse:
    """Second geodetic problem: s, α₁₂, α₂₁ from two points (no iteration)."""
    dphi = lat2 - lat1
    dlam = math.remainder(lon2 - lon1, math.tau)
    lat_m = (lat1 + lat2) / 2
    k = brackets(lat_m, ellipsoid)
    sin_m, cos_m = math.sin(lat_m), math.cos(lat_m)
    s_cos = k.m * dphi * math.cos(dlam / 2) * (1 + k.b5 * (dlam * cos_m) ** 2 + k.b6 * dphi**2)
    s_sin = k.n * dlam * cos_m * (1 - k.b3 * (dlam * sin_m) ** 2 + k.b4 * dphi**2)
    dalpha = dlam * sin_m * (1 + k.b7 * (dlam * cos_m) ** 2 + k.b8 * dphi**2)
    alpha_m = math.atan2(s_sin, s_cos)
    s = math.hypot(s_sin, s_cos)
    a12 = (alpha_m - dalpha / 2) % math.tau
    a21 = (alpha_m + dalpha / 2 + math.pi) % math.tau
    if trace is not None:
        trace.section("Setup")
        trace.step("Δφ", dphi, "φ₂ − φ₁", unit="rad")
        trace.step("Δλ", dlam, "λ₂ − λ₁", unit="rad")
        trace.step("φ_m", lat_m, "(φ₁ + φ₂)/2", unit="rad")
        _trace_brackets(trace, k)
        trace.section("Mean azimuth and length")
        trace.step("s cos α_m", s_cos, "M_m Δφ cos(Δλ/2){1 + [5](Δλ cos φ_m)² + [6]Δφ²}", unit="m")
        trace.step("s sin α_m", s_sin, "N_m Δλ cos φ_m{1 − [3](Δλ sin φ_m)² + [4]Δφ²}", unit="m")
        trace.step("Δα", dalpha, "Δλ sin φ_m{1 + [7](Δλ cos φ_m)² + [8]Δφ²}", unit="rad")
        trace.step("α_m", alpha_m, "atan2(s sin α_m, s cos α_m)", unit="rad")
        trace.section("Result")
        trace.step("s", s, "√((s sin α_m)² + (s cos α_m)²)", unit="m")
        trace.step("α₁₂", a12, "α_m − Δα/2", unit="rad")
        trace.step("α₂₁", a21, "α_m + Δα/2 ± 180°", unit="rad")
    return MidLatitudeInverse(s, a12, a21, alpha_m, dalpha)


@dataclass(frozen=True)
class MidLatitudeDirect:
    lat2: float
    lon2: float
    back_azimuth: float
    iterations: int


def mid_latitude_direct(
    lat1: float,
    lon1: float,
    azimuth: float,
    distance: float,
    ellipsoid: Ellipsoid,
    trace: Trace | None = None,
    tol: float = 1e-12,
    max_iter: int = 30,
) -> MidLatitudeDirect:
    """First geodetic problem by iterating the mid-latitude formulas.

    The first approximation uses φ₁ and α₁₂ in place of φ_m and α_m, as on the
    slides; each further step recomputes φ_m, α_m and the brackets. The [8]
    term of Δα is kept (the slides drop it here, which costs about 0.005″).
    """
    s = distance
    dphi = dlam = dalpha = 0.0
    if trace is not None:
        trace.section("Given")
        trace.step("φ₁", lat1, unit="rad")
        trace.step("α₁₂", azimuth, unit="rad")
        trace.step("s", s, unit="m")
    for i in range(1, max_iter + 1):
        lat_m = lat1 + dphi / 2
        alpha_m = azimuth + dalpha / 2
        k = brackets(lat_m, ellipsoid)
        sin_m, cos_m = math.sin(lat_m), math.cos(lat_m)
        new_dlam = s * math.sin(alpha_m) / (k.n * cos_m) / (
            1 - k.b3 * (dlam * sin_m) ** 2 + k.b4 * dphi**2
        )
        new_dphi = s * math.cos(alpha_m) / (k.m * math.cos(new_dlam / 2)) / (
            1 + k.b5 * (new_dlam * cos_m) ** 2 + k.b6 * dphi**2
        )
        new_dalpha = new_dlam * sin_m * (1 + k.b7 * (new_dlam * cos_m) ** 2 + k.b8 * new_dphi**2)
        done = max(abs(new_dlam - dlam), abs(new_dphi - dphi), abs(new_dalpha - dalpha)) < tol
        dlam, dphi, dalpha = new_dlam, new_dphi, new_dalpha
        if trace is not None:
            trace.section(f"Approximation {i}")
            trace.step("φ_m", lat_m, "φ₁ + Δφ/2", unit="rad")
            trace.step("α_m", alpha_m, "α₁₂ + Δα/2", unit="rad")
            trace.step("Δλ", dlam, "s sin α_m / (N_m cos φ_m{1 − [3](Δλ sin φ_m)² + [4]Δφ²})", unit="rad")
            trace.step("Δφ", dphi, "s cos α_m / (M_m cos(Δλ/2){1 + [5](Δλ cos φ_m)² + [6]Δφ²})", unit="rad")
            trace.step("Δα", dalpha, "Δλ sin φ_m{1 + [7](Δλ cos φ_m)² + [8]Δφ²}", unit="rad")
        if done:
            break
    else:
        raise ConvergenceError(f"did not converge in {max_iter} iterations")
    lat2 = lat1 + dphi
    lon2 = math.remainder(lon1 + dlam, math.tau)
    back = (azimuth + dalpha + math.pi) % math.tau
    if trace is not None:
        trace.section("Result")
        trace.step("φ₂", lat2, "φ₁ + Δφ", unit="rad")
        trace.step("λ₂", lon2, "λ₁ + Δλ", unit="rad")
        trace.step("α₂₁", back, "α₁₂ + Δα ± 180°", unit="rad")
    return MidLatitudeDirect(lat2, lon2, back, i)
