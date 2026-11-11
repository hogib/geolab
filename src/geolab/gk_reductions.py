"""Reductions between the ellipsoid and the Gauss–Krüger plane.

A geodesic PᵢPₖ of length S, with grid bearings Tᵢₖ, Tₖᵢ (azimuth minus grid
convergence), maps to a curve in the plane; the straight line between the
images Pᵢ′, Pₖ′ has length s and plane bearings tᵢₖ, tₖᵢ = tᵢₖ ± 180°. With R
= √(MN) at the latitude φ_m of X_m = (Xᵢ + Xₖ)/2, Y_m = (Yᵢ + Yₖ)/2 (course
notation, angles in radians, so ρ = 1):

    S − s   = −s/(6R²)(Yᵢ² + YᵢYₖ + Yₖ²) + 5sY_m⁴/(24R⁴)
    Tᵢₖ − tᵢₖ = Y_m(Xₖ − Xᵢ)/(2R²)(1 − Y_m²/(3R²)) − (Xₖ − Xᵢ)(Yₖ − Yᵢ)/(12R²)
    Tₖᵢ − tₖᵢ = Y_m(Xᵢ − Xₖ)/(2R²)(1 − Y_m²/(3R²)) − (Xᵢ − Xₖ)(Yᵢ − Yₖ)/(12R²)

so s > S always. With these the two main problems are solved in the plane:
the second directly, the first by iterating from t ≈ T, s ≈ S.

The area F on the ellipsoid of a figure with plane area f follows from the area
distortion 1/m² = 1 − Y²/R² + 2Y⁴/(3R⁴):

    F = f − (1/R²)∬Y² df + (2/(3R⁴))∬Y⁴ df

The integrals are evaluated exactly over the polygon; for a figure spanning
Yᵢ…Yₖ evenly this is the course's δF = −f(Yᵢ² + YᵢYₖ + Yₖ²)/(3R²) + 2Y_m⁴f/(3R⁴).

Lengths are in metres and angles in radians.
"""

import math
from dataclasses import dataclass

from .arcs import latitude_from_meridian_distance
from .curvature import meridian_radius, prime_vertical_radius
from .ellipsoids import Ellipsoid
from .trace import Trace


def _gauss_radius(x: float, ellipsoid: Ellipsoid) -> float:
    lat = latitude_from_meridian_distance(x, ellipsoid).lat
    return math.sqrt(meridian_radius(lat, ellipsoid) * prime_vertical_radius(lat, ellipsoid))


@dataclass(frozen=True)
class LineReduction:
    plane_distance: float
    """s, the straight-line distance in the plane."""
    plane_bearing: float
    """tᵢₖ."""
    distance: float
    """S, the geodesic length on the ellipsoid."""
    bearing1: float
    """Tᵢₖ, the geodesic's grid bearing at Pᵢ."""
    bearing2: float
    """Tₖᵢ."""
    distance_reduction: float
    """δS = S − s."""
    bearing_reduction1: float
    """δTᵢₖ = Tᵢₖ − tᵢₖ."""
    bearing_reduction2: float
    """δTₖᵢ = Tₖᵢ − tₖᵢ."""
    radius: float
    """R_m."""


def _reductions(xi: float, yi: float, xk: float, yk: float, r: float) -> tuple[float, float, float]:
    """δS/s, δTᵢₖ, δTₖᵢ."""
    ym = (yi + yk) / 2
    r2 = r * r
    ds_rel = -(yi**2 + yi * yk + yk**2) / (6 * r2) + 5 * ym**4 / (24 * r2 * r2)
    k = 1 - ym**2 / (3 * r2)
    dt_ik = ym * (xk - xi) / (2 * r2) * k - (xk - xi) * (yk - yi) / (12 * r2)
    dt_ki = ym * (xi - xk) / (2 * r2) * k - (xi - xk) * (yi - yk) / (12 * r2)
    return ds_rel, dt_ik, dt_ki


def reduce_line(
    xi: float, yi: float, xk: float, yk: float, ellipsoid: Ellipsoid, trace: Trace | None = None
) -> LineReduction:
    """Second main problem by reductions: S, Tᵢₖ, Tₖᵢ from the plane coordinates."""
    dx, dy = xk - xi, yk - yi
    s = math.hypot(dx, dy)
    t_ik = math.atan2(dy, dx) % math.tau
    t_ki = (t_ik + math.pi) % math.tau
    r = _gauss_radius((xi + xk) / 2, ellipsoid)
    ds_rel, dt_ik, dt_ki = _reductions(xi, yi, xk, yk, r)
    big_s = s * (1 + ds_rel)
    big_t_ik = (t_ik + dt_ik) % math.tau
    big_t_ki = (t_ki + dt_ki) % math.tau
    if trace is not None:
        trace.section("Plane")
        trace.step("s", s, "√(ΔX² + ΔY²)", unit="m")
        trace.step("tᵢₖ", t_ik, "atan(ΔY/ΔX)", unit="rad")
        trace.step("tₖᵢ", t_ki, "tᵢₖ ± 180°", unit="rad")
        trace.section("Reductions")
        trace.step("Y_m", (yi + yk) / 2, "(Yᵢ + Yₖ)/2", unit="m")
        trace.step("R_m", r, "√(MN) at the foot point of X_m", unit="m")
        trace.step("δS", big_s - s, "−s/(6R²)(Yᵢ² + YᵢYₖ + Yₖ²) + 5sY_m⁴/(24R⁴)", unit="m")
        trace.step("δTᵢₖ", dt_ik, "Y_mΔX/(2R²)(1 − Y_m²/(3R²)) − ΔXΔY/(12R²)", unit="rad")
        trace.step("δTₖᵢ", dt_ki, "−Y_mΔX/(2R²)(1 − Y_m²/(3R²)) − ΔXΔY/(12R²)", unit="rad")
        trace.section("Ellipsoid")
        trace.step("S", big_s, "s + δS", unit="m")
        trace.step("Tᵢₖ", big_t_ik, "tᵢₖ + δTᵢₖ", unit="rad")
        trace.step("Tₖᵢ", big_t_ki, "tₖᵢ + δTₖᵢ", unit="rad")
    return LineReduction(s, t_ik, big_s, big_t_ik, big_t_ki, big_s - s, dt_ik, dt_ki, r)


@dataclass(frozen=True)
class ReducedDirect:
    x2: float
    y2: float
    bearing2: float
    """Tₖᵢ, the geodesic's grid bearing at Pₖ back to Pᵢ."""
    plane_distance: float
    plane_bearing: float
    iterations: int


def direct_by_reductions(
    xi: float,
    yi: float,
    bearing: float,
    distance: float,
    ellipsoid: Ellipsoid,
    trace: Trace | None = None,
) -> ReducedDirect:
    """First main problem by reductions: Xₖ, Yₖ, Tₖᵢ from Xᵢ, Yᵢ, Tᵢₖ and S.

    Starts from t ≈ T, s ≈ S and repeats: Pₖ from s, t; reductions from Pᵢ, Pₖ;
    s = S − δS, t = T − δT.
    """
    t, s = bearing, distance
    xk = yk = 0.0
    for i in range(1, 20):
        xk, yk = xi + s * math.cos(t), yi + s * math.sin(t)
        r = _gauss_radius((xi + xk) / 2, ellipsoid)
        ds_rel, dt_ik, dt_ki = _reductions(xi, yi, xk, yk, r)
        s_new, t_new = distance / (1 + ds_rel), bearing - dt_ik
        done = abs(s_new - s) < 1e-7 and abs(t_new - t) < 1e-13
        s, t = s_new, t_new
        if done:
            break
    xk, yk = xi + s * math.cos(t), yi + s * math.sin(t)
    t_ki = (t + math.pi) % math.tau
    big_t_ki = (t_ki + dt_ki) % math.tau
    if trace is not None:
        trace.section("Reductions (last iteration)")
        trace.step("R_m", r, "√(MN) at the foot point of X_m", unit="m")
        trace.step("δS", distance - s, "S − s", unit="m")
        trace.step("δTᵢₖ", dt_ik, unit="rad")
        trace.step("δTₖᵢ", dt_ki, unit="rad")
        trace.section("Plane")
        trace.step("s", s, "S − δS", unit="m")
        trace.step("tᵢₖ", t % math.tau, "Tᵢₖ − δTᵢₖ", unit="rad")
        trace.step("Xₖ", xk, "Xᵢ + s cos tᵢₖ", unit="m")
        trace.step("Yₖ", yk, "Yᵢ + s sin tᵢₖ", unit="m")
        trace.step("tₖᵢ", t_ki, "tᵢₖ ± 180°", unit="rad")
        trace.step("Tₖᵢ", big_t_ki, f"tₖᵢ + δTₖᵢ ({i} iterations)", unit="rad")
    return ReducedDirect(xk, yk, big_t_ki, s, t % math.tau, i)


def _moment(points: list[tuple[float, float]], n: int) -> float:
    """∬ Yⁿ df over the polygon (signed by orientation), by Green's theorem."""
    total = 0.0
    for (x1, y1), (x2, y2) in zip(points, points[1:] + points[:1]):
        # ∬ Yⁿ df = ∮ Yⁿ⁺¹/(n+1) dX, positive counter-clockwise with Y east and X north;
        # along an edge, ∫₀¹ Y(τ)ⁿ⁺¹ dτ = Σⱼ Y₁ʲ Y₂ⁿ⁺¹⁻ʲ / (n+2).
        poly = sum(y1**j * y2 ** (n + 1 - j) for j in range(n + 2)) / (n + 2)
        total += (x2 - x1) * poly / (n + 1)
    return total


@dataclass(frozen=True)
class AreaReduction:
    plane_area: float
    """f."""
    area: float
    """F on the ellipsoid."""
    reduction: float
    """δF = F − f."""
    radius: float


def reduce_area(points: list[tuple[float, float]], ellipsoid: Ellipsoid, trace: Trace | None = None) -> AreaReduction:
    """Ellipsoid area F of the polygon with Gauss–Krüger vertices (X, Y)."""
    if len(points) < 3:
        raise ValueError("a polygon needs at least three vertices")
    f = _moment(points, 0)
    sign = 1 if f >= 0 else -1
    f, y2, y4 = sign * f, sign * _moment(points, 2), sign * _moment(points, 4)
    xm = sum(p[0] for p in points) / len(points)
    r = _gauss_radius(xm, ellipsoid)
    big_f = f - y2 / r**2 + 2 * y4 / (3 * r**4)
    if trace is not None:
        trace.step("f", f, "plane area (shoelace)")
        trace.step("X_m", xm, "mean of the vertices' X", unit="m")
        trace.step("R_m", r, "√(MN) at the foot point of X_m", unit="m")
        trace.step("∬Y² df", y2, "exact over the polygon")
        trace.step("∬Y⁴ df", y4, "exact over the polygon")
        trace.step("Y_rms", math.sqrt(y2 / f) if f else 0.0, "√(∬Y² df / f)", unit="m")
        trace.step("δF", big_f - f, "−∬Y² df/R² + 2∬Y⁴ df/(3R⁴)")
        trace.step("F", big_f, "f + δF")
    return AreaReduction(f, big_f, big_f - f, r)


@dataclass(frozen=True)
class Traverse:
    points: list[tuple[float, float]]
    """Adjusted X, Y of the new points P₁…Pₙ."""
    angle_misclosure: float
    """w_β: known minus computed closing bearing, before adjustment."""
    misclosure_x: float
    """w_X: known minus computed X of the closing point."""
    misclosure_y: float
    plane_distances: list[float]
    """s, the reduced sides."""
    plane_angles: list[float]
    """β′, the reduced angles before the misclosure is distributed."""


def traverse(
    a: tuple[float, float],
    b: tuple[float, float],
    c: tuple[float, float],
    d: tuple[float, float],
    angles: list[float],
    distances: list[float],
    ellipsoid: Ellipsoid,
    trace: Trace | None = None,
) -> Traverse:
    """Traverse between known points B and C, oriented on A and D.

    ``angles`` are the measured (ellipsoidal) angles β₀ at B, β₁…βₙ at the new
    points and βₙ₊₁ at C, each clockwise from the back sight to the fore sight;
    ``distances`` the geodesic sides S from B to C. Angles and sides are reduced
    to the plane (iterating with provisional coordinates), the angular
    misclosure is spread equally over the angles and the coordinate misclosure
    in proportion to the side lengths.
    """
    n_new = len(distances) - 1
    if n_new < 1 or len(angles) != len(distances) + 1:
        raise ValueError("need n + 1 sides and n + 2 angles for n new points")

    def bearing(p: tuple[float, float], q: tuple[float, float]) -> float:
        return math.atan2(q[1] - p[1], q[0] - p[0]) % math.tau

    t_ab, t_cd = bearing(a, b), bearing(c, d)
    stations = [b] + [(0.0, 0.0)] * n_new + [c]
    # Provisional plane angles and sides start as the measured ones.
    plane_angles, plane_sides = list(angles), list(distances)
    for _ in range(5):
        t = (t_ab + math.pi) % math.tau  # bearing B → A
        pts = [b]
        bearings = []
        for k, side in enumerate(plane_sides):
            t = (t + plane_angles[k]) % math.tau
            bearings.append(t)
            x, y = pts[-1]
            pts.append((x + side * math.cos(t), y + side * math.sin(t)))
            t = (t + math.pi) % math.tau
        stations = pts[:-1] + [c]
        # Reduce each measured direction with the provisional coordinates.
        new_angles, new_sides = [], []
        for k in range(len(distances) + 1):
            here = stations[k]
            fore = stations[k + 1] if k < len(distances) else d
            prev = a if k == 0 else stations[k - 1]
            r_f = _gauss_radius((here[0] + fore[0]) / 2, ellipsoid)
            r_b = _gauss_radius((here[0] + prev[0]) / 2, ellipsoid)
            _, dt_fore, _ = _reductions(*here, *fore, r_f)
            _, dt_back, _ = _reductions(*here, *prev, r_b)
            new_angles.append(angles[k] - dt_fore + dt_back)
            if k < len(distances):
                ds_rel, _, _ = _reductions(*here, *fore, r_f)
                new_sides.append(distances[k] / (1 + ds_rel))
        plane_angles, plane_sides = new_angles, new_sides

    # Angular misclosure.
    t = (t_ab + math.pi) % math.tau
    for beta in plane_angles:
        t = (t + beta) % math.tau
        t_computed_cd = t
        t = (t + math.pi) % math.tau
    w_beta = math.remainder(t_cd - t_computed_cd, math.tau)
    v = w_beta / len(plane_angles)
    # Bearings with the corrected angles, then coordinate misclosure.
    t = (t_ab + math.pi) % math.tau
    dxs, dys = [], []
    for k, side in enumerate(plane_sides):
        t = (t + plane_angles[k] + v) % math.tau
        dxs.append(side * math.cos(t))
        dys.append(side * math.sin(t))
        t = (t + math.pi) % math.tau
    w_x = c[0] - b[0] - sum(dxs)
    w_y = c[1] - b[1] - sum(dys)
    total = sum(plane_sides)
    x, y = b
    points = []
    for k in range(n_new):
        x += dxs[k] + w_x * plane_sides[k] / total
        y += dys[k] + w_y * plane_sides[k] / total
        points.append((x, y))
    if trace is not None:
        trace.section("Reduced to the plane")
        for k, (beta, beta_p) in enumerate(zip(angles, plane_angles)):
            trace.step(f"β{k}′ − β{k}", beta_p - beta, "δT(back) − δT(fore)", unit="rad")
        for k, (big_s, s) in enumerate(zip(distances, plane_sides), 1):
            trace.step(f"s{k}", s, f"S{k} − δS", unit="m")
        trace.section("Misclosure")
        trace.step("t_AB", t_ab, unit="rad")
        trace.step("t_CD", t_cd, unit="rad")
        trace.step("w_β", w_beta, "t_CD − (t_AB + Σβ′ − (n + 1)·180°)", unit="rad")
        trace.step("v_β", v, "w_β / (n + 2) per angle", unit="rad")
        trace.step("w_X", w_x, "X_C − X_B − ΣΔX", unit="m")
        trace.step("w_Y", w_y, "Y_C − Y_B − ΣΔY", unit="m")
        trace.step("w_s", math.hypot(w_x, w_y), "√(w_X² + w_Y²)", unit="m")
        trace.section("Adjusted points")
        for k, (px, py) in enumerate(points, 1):
            trace.step(f"X{k}", px, unit="m")
            trace.step(f"Y{k}", py, unit="m")
    return Traverse(points, w_beta, w_x, w_y, plane_sides, plane_angles)
