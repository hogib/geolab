"""Inverse and direct geodetic problems on the ellipsoid (Vincenty, 1975).

All angles are in radians. Azimuths are measured clockwise from north and
returned in [0, 2π); longitudes are returned in [−π, π].
"""

import math
from dataclasses import dataclass

from .ellipsoids import Ellipsoid
from .errors import ConvergenceError
from .trace import Trace


@dataclass(frozen=True)
class InverseResult:
    distance: float
    azimuth1: float
    """Forward azimuth at point 1."""
    azimuth2: float
    """Forward azimuth at point 2 (direction of travel when arriving)."""
    iterations: int

    @property
    def back_azimuth(self) -> float:
        """Azimuth from point 2 back towards point 1."""
        return _azimuth(self.azimuth2 + math.pi)


@dataclass(frozen=True)
class DirectResult:
    lat2: float
    lon2: float
    azimuth2: float
    """Forward azimuth at point 2 (direction of travel when arriving)."""
    iterations: int

    @property
    def back_azimuth(self) -> float:
        """Azimuth from point 2 back towards point 1."""
        return _azimuth(self.azimuth2 + math.pi)


def inverse(
    lat1: float,
    lon1: float,
    lat2: float,
    lon2: float,
    ellipsoid: Ellipsoid,
    trace: Trace | None = None,
    tol: float = 1e-12,
    max_iter: int = 200,
) -> InverseResult:
    """Distance and azimuths between two points.

    Raises ConvergenceError for nearly antipodal points, where Vincenty's
    iteration on λ does not converge.
    """
    a, b, f = ellipsoid.a, ellipsoid.b, ellipsoid.f
    L = math.remainder(lon2 - lon1, math.tau)
    u1 = _reduced_latitude(lat1, f)
    u2 = _reduced_latitude(lat2, f)
    sin_u1, cos_u1 = math.sin(u1), math.cos(u1)
    sin_u2, cos_u2 = math.sin(u2), math.cos(u2)

    if trace is not None:
        trace.section("Setup")
        trace.step("L", L, "λ₂ − λ₁", unit="rad")
        trace.step("U₁", u1, "atan((1 − f) tan φ₁)", unit="rad")
        trace.step("U₂", u2, "atan((1 − f) tan φ₂)", unit="rad")

    lam = L
    for i in range(1, max_iter + 1):
        sin_lam, cos_lam = math.sin(lam), math.cos(lam)
        sin_sigma = math.hypot(cos_u2 * sin_lam, cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lam)
        if sin_sigma == 0:
            if trace is not None:
                trace.section("Coincident points")
                trace.step("s", 0.0, unit="m")
            return InverseResult(0.0, 0.0, 0.0, i)
        cos_sigma = sin_u1 * sin_u2 + cos_u1 * cos_u2 * cos_lam
        sigma = math.atan2(sin_sigma, cos_sigma)
        sin_alpha = cos_u1 * cos_u2 * sin_lam / sin_sigma
        cos_sq_alpha = 1 - sin_alpha**2
        # On the equator cos²α = 0 and the term is undefined; it drops out (C = 0).
        cos_2sigma_m = cos_sigma - 2 * sin_u1 * sin_u2 / cos_sq_alpha if cos_sq_alpha else 0.0
        c = f / 16 * cos_sq_alpha * (4 + f * (4 - 3 * cos_sq_alpha))
        new_lam = L + (1 - c) * f * sin_alpha * (
            sigma + c * sin_sigma * (cos_2sigma_m + c * cos_sigma * (-1 + 2 * cos_2sigma_m**2))
        )

        if trace is not None:
            trace.section(f"Iteration {i}")
            trace.step("sin σ", sin_sigma, "√((cos U₂ sin λ)² + (cos U₁ sin U₂ − sin U₁ cos U₂ cos λ)²)")
            trace.step("cos σ", cos_sigma, "sin U₁ sin U₂ + cos U₁ cos U₂ cos λ")
            trace.step("σ", sigma, "atan2(sin σ, cos σ)", unit="rad")
            trace.step("sin α", sin_alpha, "cos U₁ cos U₂ sin λ / sin σ")
            trace.step("cos²α", cos_sq_alpha, "1 − sin²α")
            trace.step("cos 2σm", cos_2sigma_m, "cos σ − 2 sin U₁ sin U₂ / cos²α")
            trace.step("C", c, "f/16 cos²α (4 + f(4 − 3 cos²α))")
            trace.step("λ", new_lam, "L + (1 − C) f sin α (σ + C sin σ (cos 2σm + C cos σ (−1 + 2 cos² 2σm)))", unit="rad")

        if abs(new_lam) > math.pi:
            raise ConvergenceError("λ exceeded π: points are nearly antipodal")
        converged = abs(new_lam - lam) < tol
        lam = new_lam
        if converged:
            break
    else:
        raise ConvergenceError(
            f"λ did not converge in {max_iter} iterations: points are nearly antipodal"
        )

    # Recompute with the converged λ so every value below is consistent.
    sin_lam, cos_lam = math.sin(lam), math.cos(lam)
    sin_sigma = math.hypot(cos_u2 * sin_lam, cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lam)
    cos_sigma = sin_u1 * sin_u2 + cos_u1 * cos_u2 * cos_lam
    sigma = math.atan2(sin_sigma, cos_sigma)
    sin_alpha = cos_u1 * cos_u2 * sin_lam / sin_sigma
    cos_sq_alpha = 1 - sin_alpha**2
    cos_2sigma_m = cos_sigma - 2 * sin_u1 * sin_u2 / cos_sq_alpha if cos_sq_alpha else 0.0

    u_sq, big_a, big_b = _series_coefficients(cos_sq_alpha, a, b)
    delta_sigma = _delta_sigma(big_b, sin_sigma, cos_sigma, cos_2sigma_m)
    s = b * big_a * (sigma - delta_sigma)
    alpha1 = _azimuth(math.atan2(cos_u2 * sin_lam, cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lam))
    alpha2 = _azimuth(math.atan2(cos_u1 * sin_lam, -sin_u1 * cos_u2 + cos_u1 * sin_u2 * cos_lam))

    if trace is not None:
        trace.section("Result")
        trace.step("u²", u_sq, "cos²α (a² − b²) / b²")
        trace.step("A", big_a, "1 + u²/16384 (4096 + u²(−768 + u²(320 − 175u²)))")
        trace.step("B", big_b, "u²/1024 (256 + u²(−128 + u²(74 − 47u²)))")
        trace.step("Δσ", delta_sigma, _DELTA_SIGMA_FORMULA, unit="rad")
        trace.step("s", s, "b A (σ − Δσ)", unit="m")
        trace.step("α₁", alpha1, "atan2(cos U₂ sin λ, cos U₁ sin U₂ − sin U₁ cos U₂ cos λ)", unit="rad")
        trace.step("α₂", alpha2, "atan2(cos U₁ sin λ, −sin U₁ cos U₂ + cos U₁ sin U₂ cos λ)", unit="rad")

    return InverseResult(s, alpha1, alpha2, i)


def direct(
    lat1: float,
    lon1: float,
    azimuth1: float,
    distance: float,
    ellipsoid: Ellipsoid,
    trace: Trace | None = None,
    tol: float = 1e-12,
    max_iter: int = 200,
) -> DirectResult:
    """End point and arrival azimuth after travelling ``distance`` metres from
    point 1 along the geodesic with forward azimuth ``azimuth1``."""
    a, b, f = ellipsoid.a, ellipsoid.b, ellipsoid.f
    u1 = _reduced_latitude(lat1, f)
    sin_u1, cos_u1 = math.sin(u1), math.cos(u1)
    sin_a1, cos_a1 = math.sin(azimuth1), math.cos(azimuth1)

    sigma1 = math.atan2(sin_u1, cos_u1 * cos_a1)
    sin_alpha = cos_u1 * sin_a1
    cos_sq_alpha = 1 - sin_alpha**2
    u_sq, big_a, big_b = _series_coefficients(cos_sq_alpha, a, b)
    sigma0 = distance / (b * big_a)

    if trace is not None:
        trace.section("Setup")
        trace.step("U₁", u1, "atan((1 − f) tan φ₁)", unit="rad")
        trace.step("σ₁", sigma1, "atan2(tan U₁, cos α₁)", unit="rad")
        trace.step("sin α", sin_alpha, "cos U₁ sin α₁")
        trace.step("cos²α", cos_sq_alpha, "1 − sin²α")
        trace.step("u²", u_sq, "cos²α (a² − b²) / b²")
        trace.step("A", big_a, "1 + u²/16384 (4096 + u²(−768 + u²(320 − 175u²)))")
        trace.step("B", big_b, "u²/1024 (256 + u²(−128 + u²(74 − 47u²)))")
        trace.step("σ", sigma0, "s / (b A)", unit="rad")

    sigma = sigma0
    for i in range(1, max_iter + 1):
        cos_2sigma_m = math.cos(2 * sigma1 + sigma)
        delta_sigma = _delta_sigma(big_b, math.sin(sigma), math.cos(sigma), cos_2sigma_m)
        new_sigma = sigma0 + delta_sigma

        if trace is not None:
            trace.section(f"Iteration {i}")
            trace.step("cos 2σm", cos_2sigma_m, "cos(2σ₁ + σ)")
            trace.step("Δσ", delta_sigma, _DELTA_SIGMA_FORMULA, unit="rad")
            trace.step("σ", new_sigma, "s / (b A) + Δσ", unit="rad")

        converged = abs(new_sigma - sigma) < tol
        sigma = new_sigma
        if converged:
            break
    else:
        raise ConvergenceError(f"σ did not converge in {max_iter} iterations")

    sin_sigma, cos_sigma = math.sin(sigma), math.cos(sigma)
    cos_2sigma_m = math.cos(2 * sigma1 + sigma)
    tmp = sin_u1 * sin_sigma - cos_u1 * cos_sigma * cos_a1
    lat2 = math.atan2(
        sin_u1 * cos_sigma + cos_u1 * sin_sigma * cos_a1,
        (1 - f) * math.hypot(sin_alpha, tmp),
    )
    lam = math.atan2(sin_sigma * sin_a1, cos_u1 * cos_sigma - sin_u1 * sin_sigma * cos_a1)
    c = f / 16 * cos_sq_alpha * (4 + f * (4 - 3 * cos_sq_alpha))
    L = lam - (1 - c) * f * sin_alpha * (
        sigma + c * sin_sigma * (cos_2sigma_m + c * cos_sigma * (-1 + 2 * cos_2sigma_m**2))
    )
    lon2 = math.remainder(lon1 + L, math.tau)
    alpha2 = _azimuth(math.atan2(sin_alpha, -tmp))

    if trace is not None:
        trace.section("Result")
        trace.step("φ₂", lat2, "atan2(sin U₁ cos σ + cos U₁ sin σ cos α₁, (1 − f)√(sin²α + (sin U₁ sin σ − cos U₁ cos σ cos α₁)²))", unit="rad")
        trace.step("λ", lam, "atan2(sin σ sin α₁, cos U₁ cos σ − sin U₁ sin σ cos α₁)", unit="rad")
        trace.step("C", c, "f/16 cos²α (4 + f(4 − 3 cos²α))")
        trace.step("L", L, "λ − (1 − C) f sin α (σ + C sin σ (cos 2σm + C cos σ (−1 + 2 cos² 2σm)))", unit="rad")
        trace.step("λ₂", lon2, "λ₁ + L", unit="rad")
        trace.step("α₂", alpha2, "atan2(sin α, −sin U₁ sin σ + cos U₁ cos σ cos α₁)", unit="rad")

    return DirectResult(lat2, lon2, alpha2, i)


_DELTA_SIGMA_FORMULA = (
    "B sin σ (cos 2σm + B/4 (cos σ (−1 + 2 cos² 2σm) − B/6 cos 2σm (−3 + 4 sin²σ)(−3 + 4 cos² 2σm)))"
)


def _reduced_latitude(lat: float, f: float) -> float:
    # atan((1 − f) tan φ), written with atan2 so φ = ±90° works.
    return math.atan2((1 - f) * math.sin(lat), math.cos(lat))


def _series_coefficients(cos_sq_alpha: float, a: float, b: float) -> tuple[float, float, float]:
    u_sq = cos_sq_alpha * (a**2 - b**2) / b**2
    big_a = 1 + u_sq / 16384 * (4096 + u_sq * (-768 + u_sq * (320 - 175 * u_sq)))
    big_b = u_sq / 1024 * (256 + u_sq * (-128 + u_sq * (74 - 47 * u_sq)))
    return u_sq, big_a, big_b


def _delta_sigma(big_b: float, sin_sigma: float, cos_sigma: float, cos_2sigma_m: float) -> float:
    return big_b * sin_sigma * (
        cos_2sigma_m
        + big_b / 4 * (
            cos_sigma * (-1 + 2 * cos_2sigma_m**2)
            - big_b / 6 * cos_2sigma_m * (-3 + 4 * sin_sigma**2) * (-3 + 4 * cos_2sigma_m**2)
        )
    )


def _azimuth(angle: float) -> float:
    return angle % math.tau
