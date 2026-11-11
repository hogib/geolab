"""Gauss–Krüger (transverse Mercator) conformal projection of the ellipsoid.

The central meridian λ₀ maps to the X axis with its true length (X = G, the
meridian arc from the equator) and Y = 0. With the isometric parameters
Z = q + iλ and W = X + iY, the projection is the analytic function W = f(Z);
expanded about a point P₀(φ₀, λ₀) on the central meridian (course notation):

    ΔX + iY = a₁(Δq + iΔλ) + a₂(Δq + iΔλ)² + … + a₅(Δq + iΔλ)⁵,   X = G₀ + ΔX

    a₁ = N cos φ₀                         a₂ = −½ N t cos²φ₀
    a₃ = −⅙ N cos³φ₀ (1 − t² + η²)        a₄ = ¹⁄₂₄ N t cos⁴φ₀ (5 − t² + 9η² + 4η⁴)
    a₅ = ¹⁄₁₂₀ N cos⁵φ₀ (5 − 18t² + t⁴ + 14η² − 58η²t²)

In practice φ₀ = φ, so Δq = 0 and

    X = G − a₂Δλ² + a₄Δλ⁴        Y = a₁Δλ − a₃Δλ³ + a₅Δλ⁵

The inverse uses the foot-point latitude φ_f (G(φ_f) = X), with N, t, η² at φ_f:

    φ = φ_f + b₀₂Y² + b₀₄Y⁴        Δλ = b₀₁Y + b₀₃Y³ + b₀₅Y⁵

Scale (length distortion) m, area distortion m² (no angle distortion), and the
grid convergence c, the angle from grid north (Y = const) to the meridian:

    m = 1 + ½cos²φ(1 + η²)Δλ² + ¹⁄₂₄cos⁴φ(5 − 4t² + 14η² − 28η²t²)Δλ⁴
      = 1 + (1 + η²)Y²/(2N²) + Y⁴/(24N⁴)
    c = sin φ Δλ + ⅓ sin φ cos²φ (1 + 3η²) Δλ³ = tY/N − t(1 + t² − η²)Y³/(3N³)

The grid bearing T of a geodesic relates to its azimuth by α = T + c.

For comparison, ``krueger_forward`` and ``krueger_inverse`` evaluate the
projection to n⁶ with Krüger's series (as in Karney 2011), which is exact to
well under a millimetre anywhere within a few thousand km of λ₀.

Lengths are in metres and angles in radians.
"""

import math
from dataclasses import dataclass

from .arcs import latitude_from_meridian_distance, meridian_distance
from .curvature import meridian_radius, prime_vertical_radius
from .ellipsoids import Ellipsoid
from .projection import isometric_difference
from .trace import Trace


@dataclass(frozen=True)
class GKPoint:
    x: float
    """X, north."""
    y: float
    """Y, east of the central meridian."""
    convergence: float
    """c, the angle from grid north to true north (α = T + c)."""
    scale: float
    """m, the length distortion; the area distortion is m²."""


def _aux(lat: float, ellipsoid: Ellipsoid) -> tuple[float, float, float, float]:
    """N, cos φ, t = tan φ and η² at φ."""
    c = math.cos(lat)
    return prime_vertical_radius(lat, ellipsoid), c, math.tan(lat), ellipsoid.ep_sq * c**2


def gk_coefficients(lat0: float, ellipsoid: Ellipsoid) -> tuple[float, float, float, float, float]:
    """a₁…a₅ at φ₀."""
    n, c, t, h = _aux(lat0, ellipsoid)
    return (
        n * c,
        -n * t * c**2 / 2,
        -n * c**3 * (1 - t**2 + h) / 6,
        n * t * c**4 * (5 - t**2 + 9 * h + 4 * h**2) / 24,
        n * c**5 * (5 - 18 * t**2 + t**4 + 14 * h - 58 * h * t**2) / 120,
    )


_A_FORMULAS = (
    "N cos φ₀ = dG/dq",
    "−½ N t cos²φ₀ = ½ d²G/dq²",
    "−⅙ N cos³φ₀ (1 − t² + η²)",
    "¹⁄₂₄ N t cos⁴φ₀ (5 − t² + 9η² + 4η⁴)",
    "¹⁄₁₂₀ N cos⁵φ₀ (5 − 18t² + t⁴ + 14η² − 58η²t²)",
)


def gk_scale(lat: float, dlon: float, ellipsoid: Ellipsoid) -> float:
    """m from φ, Δλ."""
    _, c, t, h = _aux(lat, ellipsoid)
    return 1 + c**2 * (1 + h) * dlon**2 / 2 + c**4 * (5 - 4 * t**2 + 14 * h - 28 * h * t**2) * dlon**4 / 24


def gk_convergence(lat: float, dlon: float, ellipsoid: Ellipsoid) -> float:
    """c from φ, Δλ."""
    s, c = math.sin(lat), math.cos(lat)
    h = ellipsoid.ep_sq * c**2
    return s * dlon + s * c**2 * (1 + 3 * h) * dlon**3 / 3


def gauss_kruger_from_geodetic(
    lat: float,
    lon: float,
    lon0: float,
    ellipsoid: Ellipsoid,
    lat0: float | None = None,
    trace: Trace | None = None,
) -> GKPoint:
    """X, Y, c, m of φ, λ about the central meridian λ₀.

    The series is expanded about P₀(φ₀, λ₀); by default φ₀ = φ (Δq = 0), the
    course's working form. Another φ₀ uses Δq from the isometric-latitude series.
    """
    if lat0 is None:
        lat0 = lat
    dl = math.remainder(lon - lon0, math.tau)
    dq = isometric_difference(lat, lat0, ellipsoid) if lat0 != lat else 0.0
    a = gk_coefficients(lat0, ellipsoid)
    g0 = meridian_distance(lat0, ellipsoid)
    z = complex(dq, dl)
    w = sum(ak * z ** (k + 1) for k, ak in enumerate(a))
    x, y = g0 + w.real, w.imag
    m = gk_scale(lat, dl, ellipsoid)
    conv = gk_convergence(lat, dl, ellipsoid)
    if trace is not None:
        _, _, t, h = _aux(lat0, ellipsoid)
        trace.section("Setup")
        trace.step("Δλ", dl, "λ − λ₀", unit="rad")
        if dq:
            trace.step("Δq", dq, "isometric latitude series, φ − φ₀")
        trace.step("G₀", g0, "meridian arc to φ₀" if dq else "meridian arc to φ", unit="m")
        trace.step("t", t, "tan φ₀")
        trace.step("η²", h, "e′² cos²φ₀")
        trace.section("Coefficients")
        for k, (ak, formula) in enumerate(zip(a, _A_FORMULAS), 1):
            trace.step(f"a{k}", ak, formula, unit="m")
        trace.section("Result")
        if dq:
            trace.step("ΔX", w.real, "Re Σ aₖ(Δq + iΔλ)ᵏ", unit="m")
            trace.step("X", x, "G₀ + ΔX", unit="m")
            trace.step("Y", y, "Im Σ aₖ(Δq + iΔλ)ᵏ", unit="m")
        else:
            trace.step("X", x, "G − a₂Δλ² + a₄Δλ⁴", unit="m")
            trace.step("Y", y, "a₁Δλ − a₃Δλ³ + a₅Δλ⁵", unit="m")
        trace.step("m", m, "1 + ½cos²φ(1 + η²)Δλ² + ¹⁄₂₄cos⁴φ(5 − 4t² + 14η² − 28η²t²)Δλ⁴")
        trace.step("m²", m * m, "area distortion")
        trace.step("c", conv, "sin φ Δλ + ⅓ sin φ cos²φ (1 + 3η²) Δλ³", unit="rad")
    return GKPoint(x, y, conv, m)


@dataclass(frozen=True)
class GKGeodetic:
    lat: float
    lon: float
    convergence: float
    scale: float
    foot_lat: float
    """φ_f, the latitude whose meridian arc is X."""


def geodetic_from_gauss_kruger(
    x: float,
    y: float,
    lon0: float,
    ellipsoid: Ellipsoid,
    trace: Trace | None = None,
) -> GKGeodetic:
    """φ, λ, c, m of the Gauss–Krüger point X, Y (central meridian λ₀)."""
    lat_f = latitude_from_meridian_distance(x, ellipsoid).lat
    n, c, t, h = _aux(lat_f, ellipsoid)
    b01 = 1 / (n * c)
    b02 = -t * (1 + h) / (2 * n**2)
    b03 = -(1 + 2 * t**2 + h) / (6 * n**3 * c)
    b04 = t * (5 + 3 * t**2 + 6 * h - 6 * h * t**2) / (24 * n**4)
    b05 = (5 + 28 * t**2 + 24 * t**4) / (120 * n**5 * c)
    lat = lat_f + b02 * y**2 + b04 * y**4
    dl = b01 * y + b03 * y**3 + b05 * y**5
    lon = math.remainder(lon0 + dl, math.tau)
    m = 1 + (1 + h) * y**2 / (2 * n**2) + y**4 / (24 * n**4)
    conv = t * y / n - t * (1 + t**2 - h) * y**3 / (3 * n**3)
    if trace is not None:
        trace.section("Foot point")
        trace.step("φ_f", lat_f, "latitude of the meridian arc G = X", unit="rad")
        trace.step("N_f", n, unit="m")
        trace.step("t_f", t, "tan φ_f")
        trace.step("η_f²", h, "e′² cos²φ_f")
        trace.section("Coefficients")
        trace.step("b01", b01, "1 / (N cos φ_f)")
        trace.step("b02", b02, "−t (1 + η²) / (2N²)")
        trace.step("b03", b03, "−(1 + 2t² + η²) / (6N³ cos φ_f)")
        trace.step("b04", b04, "t (5 + 3t² + 6η² − 6η²t²) / (24N⁴)")
        trace.step("b05", b05, "(5 + 28t² + 24t⁴) / (120N⁵ cos φ_f)")
        trace.section("Result")
        trace.step("Δφ", b02 * y**2 + b04 * y**4, "b02 Y² + b04 Y⁴", unit="rad")
        trace.step("φ", lat, "φ_f + Δφ", unit="rad")
        trace.step("Δλ", dl, "b01 Y + b03 Y³ + b05 Y⁵", unit="rad")
        trace.step("λ", lon, "λ₀ + Δλ", unit="rad")
        trace.step("m", m, "1 + (1 + η²)Y²/(2N²) + Y⁴/(24N⁴)")
        trace.step("c", conv, "tY/N − t(1 + t² − η²)Y³/(3N³)", unit="rad")
    return GKGeodetic(lat, lon, conv, m, lat_f)


# Krüger's series to n⁶ (Karney 2011, eqs. 35 and 36).
_ALPHA = (
    (1 / 2, -2 / 3, 5 / 16, 41 / 180, -127 / 288, 7891 / 37800),
    (0, 13 / 48, -3 / 5, 557 / 1440, 281 / 630, -1983433 / 1935360),
    (0, 0, 61 / 240, -103 / 140, 15061 / 26880, 167603 / 181440),
    (0, 0, 0, 49561 / 161280, -179 / 168, 6601661 / 7257600),
    (0, 0, 0, 0, 34729 / 80640, -3418889 / 1995840),
    (0, 0, 0, 0, 0, 212378941 / 319334400),
)
_BETA = (
    (1 / 2, -2 / 3, 37 / 96, -1 / 360, -81 / 512, 96199 / 604800),
    (0, 1 / 48, 1 / 15, -437 / 1440, 46 / 105, -1118711 / 3870720),
    (0, 0, 17 / 480, -37 / 840, -209 / 4480, 5569 / 90720),
    (0, 0, 0, 4397 / 161280, -11 / 504, -830251 / 7257600),
    (0, 0, 0, 0, 4583 / 161280, -108847 / 3991680),
    (0, 0, 0, 0, 0, 20648693 / 638668800),
)


def _series(table: tuple[tuple[float, ...], ...], n: float) -> list[float]:
    return [sum(c * n ** (i + 1) for i, c in enumerate(row)) for row in table]


def _rectifying_radius(ellipsoid: Ellipsoid) -> float:
    n = ellipsoid.n
    return ellipsoid.a / (1 + n) * (1 + n**2 / 4 + n**4 / 64 + n**6 / 256)


def _conformal_tan(tau: float, e: float) -> float:
    sigma = math.sinh(e * math.atanh(e * tau / math.hypot(1, tau)))
    return tau * math.hypot(1, sigma) - sigma * math.hypot(1, tau)


def krueger_forward(lat: float, lon: float, lon0: float, ellipsoid: Ellipsoid) -> GKPoint:
    """X, Y, c, m by Krüger's n⁶ series: the reference for the course's series."""
    e = math.sqrt(ellipsoid.e_sq)
    dl = math.remainder(lon - lon0, math.tau)
    tau = math.tan(lat)
    taup = _conformal_tan(tau, e)
    xip = math.atan2(taup, math.cos(dl))
    etap = math.asinh(math.sin(dl) / math.hypot(taup, math.cos(dl)))
    alpha = _series(_ALPHA, ellipsoid.n)
    xi, eta, p, q = xip, etap, 1.0, 0.0
    for j, aj in enumerate(alpha, 1):
        xi += aj * math.sin(2 * j * xip) * math.cosh(2 * j * etap)
        eta += aj * math.cos(2 * j * xip) * math.sinh(2 * j * etap)
        p += 2 * j * aj * math.cos(2 * j * xip) * math.cosh(2 * j * etap)
        q += 2 * j * aj * math.sin(2 * j * xip) * math.sinh(2 * j * etap)
    big_a = _rectifying_radius(ellipsoid)
    conv = math.atan(taup / math.hypot(1, taup) * math.tan(dl)) + math.atan2(q, p)
    k = (
        big_a / ellipsoid.a * math.hypot(p, q)
        * math.sqrt(1 - ellipsoid.e_sq * math.sin(lat) ** 2) * math.hypot(1, tau)
        / math.hypot(taup, math.cos(dl))
    )
    return GKPoint(big_a * xi, big_a * eta, conv, k)


def krueger_inverse(x: float, y: float, lon0: float, ellipsoid: Ellipsoid) -> tuple[float, float]:
    """φ, λ from X, Y by Krüger's n⁶ series."""
    e = math.sqrt(ellipsoid.e_sq)
    big_a = _rectifying_radius(ellipsoid)
    xi, eta = x / big_a, y / big_a
    xip, etap = xi, eta
    for j, bj in enumerate(_series(_BETA, ellipsoid.n), 1):
        xip -= bj * math.sin(2 * j * xi) * math.cosh(2 * j * eta)
        etap -= bj * math.cos(2 * j * xi) * math.sinh(2 * j * eta)
    taup = math.sin(xip) / math.hypot(math.sinh(etap), math.cos(xip))
    dl = math.atan2(math.sinh(etap), math.cos(xip))
    tau = taup
    e2m = 1 - ellipsoid.e_sq
    for _ in range(10):  # Newton: dτ′/dτ = (1 − e²)√(1 + τ′²)√(1 + τ²)/(1 + (1 − e²)τ²)
        tp = _conformal_tan(tau, e)
        dtau = (taup - tp) * (1 + e2m * tau**2) / (e2m * math.hypot(1, tp) * math.hypot(1, tau))
        tau += dtau
        if abs(dtau) < 1e-15 * max(1.0, abs(tau)):
            break
    return math.atan(tau), math.remainder(lon0 + dl, math.tau)


def _gauss_radius(x: float, ellipsoid: Ellipsoid) -> tuple[float, float]:
    """φ at the foot point of X, and R = √(MN) there."""
    lat = latitude_from_meridian_distance(x, ellipsoid).lat
    return lat, math.sqrt(meridian_radius(lat, ellipsoid) * prime_vertical_radius(lat, ellipsoid))


@dataclass(frozen=True)
class GKDirect:
    x2: float
    y2: float
    bearing2: float
    """T₂₁, the grid bearing of the geodesic at P₂ back to P₁."""


def gk_direct(
    x1: float,
    y1: float,
    bearing: float,
    distance: float,
    ellipsoid: Ellipsoid,
    trace: Trace | None = None,
) -> GKDirect:
    """First problem in Gauss–Krüger X, Y: X₂, Y₂, T₂₁ from X₁, Y₁, the
    geodesic's grid bearing T₁₂ at P₁ (α₁₂ − c₁) and its length S."""
    lat1, r = _gauss_radius(x1, ellipsoid)
    u, v = distance * math.cos(bearing), distance * math.sin(bearing)
    r2 = r * r
    x2 = x1 + u + y1**2 * u / (2 * r2) + y1 * u * v / r2 + u * v**2 / (3 * r2)
    y2 = y1 + v + y1**2 * v / (2 * r2) - y1 * (u**2 - v**2) / (2 * r2) - v * (u**2 - v**2) / (6 * r2)
    dt = -y1 * u / r2 - u * v / (2 * r2)
    bearing2 = (bearing + math.pi + dt) % math.tau
    if trace is not None:
        trace.section("Setup")
        trace.step("φ₁", lat1, "foot-point latitude of X₁", unit="rad")
        trace.step("R₁", r, "√(M₁N₁)", unit="m")
        trace.step("u", u, "S cos T₁₂", unit="m")
        trace.step("v", v, "S sin T₁₂", unit="m")
        trace.section("Result")
        trace.step("X₂", x2, "X₁ + u + Y₁²u/(2R²) + Y₁uv/R² + uv²/(3R²)", unit="m")
        trace.step("Y₂", y2, "Y₁ + v + Y₁²v/(2R²) − Y₁(u² − v²)/(2R²) − v(u² − v²)/(6R²)", unit="m")
        trace.step("ΔT", dt, "−Y₁u/R² − uv/(2R²)", unit="rad")
        trace.step("T₂₁", bearing2, "T₁₂ + ΔT ± 180°", unit="rad")
    return GKDirect(x2, y2, bearing2)


@dataclass(frozen=True)
class GKInverse:
    distance: float
    """S, the geodesic length on the ellipsoid."""
    bearing1: float
    """T₁₂, the geodesic's grid bearing at P₁."""
    bearing2: float
    """T₂₁."""


def gk_inverse(
    x1: float,
    y1: float,
    x2: float,
    y2: float,
    ellipsoid: Ellipsoid,
    trace: Trace | None = None,
) -> GKInverse:
    """Second problem in Gauss–Krüger X, Y (mid-point series): S, T₁₂, T₂₁."""
    xm, ym = (x1 + x2) / 2, (y1 + y2) / 2
    dx, dy = x2 - x1, y2 - y1
    lat_m, r = _gauss_radius(xm, ellipsoid)
    r2, r4 = r**2, r**4
    s_sin = dy * (1 - ym**2 / (2 * r2) - dx**2 / (12 * r2) - dy**2 / (24 * r2) + 5 * ym**4 / (24 * r4))
    s_cos = dx * (1 - ym**2 / (2 * r2) + dy**2 / (24 * r2) + 5 * ym**4 / (24 * r4))
    dt = -dx * (ym / r2 - ym**3 / (3 * r4))
    tm = math.atan2(s_sin, s_cos)
    s = math.hypot(s_sin, s_cos)
    t1 = (tm - dt / 2) % math.tau
    t2 = (tm + dt / 2 + math.pi) % math.tau
    if trace is not None:
        trace.section("Setup")
        trace.step("ΔX", dx, "X₂ − X₁", unit="m")
        trace.step("ΔY", dy, "Y₂ − Y₁", unit="m")
        trace.step("Y_m", ym, "(Y₁ + Y₂)/2", unit="m")
        trace.step("φ_m", lat_m, "foot-point latitude of X_m", unit="rad")
        trace.step("R_m", r, "√(M_m N_m)", unit="m")
        trace.section("Result")
        trace.step("S sin T_m", s_sin, "ΔY(1 − Y_m²/(2R²) − ΔX²/(12R²) − ΔY²/(24R²) + 5Y_m⁴/(24R⁴))", unit="m")
        trace.step("S cos T_m", s_cos, "ΔX(1 − Y_m²/(2R²) + ΔY²/(24R²) + 5Y_m⁴/(24R⁴))", unit="m")
        trace.step("ΔT", dt, "−ΔX(Y_m/R² − Y_m³/(3R⁴))", unit="rad")
        trace.step("T_m", tm, "atan(S sin T_m / S cos T_m)", unit="rad")
        trace.step("S", s, "√((S sin T_m)² + (S cos T_m)²)", unit="m")
        trace.step("T₁₂", t1, "T_m − ΔT/2", unit="rad")
        trace.step("T₂₁", t2, "T_m + ΔT/2 ± 180°", unit="rad")
    return GKInverse(s, t1, t2)
