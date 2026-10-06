import math
from dataclasses import dataclass
from enum import Enum, auto

from .angles import format_degrees
from .ellipsoids import Ellipsoid
from .errors import ConvergenceError
from .trace import Trace


class LatType(Enum):
    GEODETIC = auto()
    GEOCENTRIC = auto()
    PARAMETRIC = auto()


@dataclass
class Latitude:
    value: float
    lat_type: LatType

    def convert_to(
        self, target_type: LatType, ellipsoid: Ellipsoid, trace: Trace | None = None
    ) -> "Latitude":
        if self.lat_type == target_type:
            return Latitude(self.value, self.lat_type)

        e_sq = ellipsoid.eccentricity_squared()
        if trace is not None:
            trace.step("e²", e_sq)

        geodetic_rad = self.value

        if self.lat_type == LatType.GEOCENTRIC:
            geodetic_rad = math.atan(math.tan(self.value) / (1 - e_sq))
            if trace is not None:
                trace.step("φ", geodetic_rad, "atan(tan ψ / (1 − e²))")
        elif self.lat_type == LatType.PARAMETRIC:
            geodetic_rad = math.atan(math.tan(self.value) / math.sqrt(1 - e_sq))
            if trace is not None:
                trace.step("φ", geodetic_rad, "atan(tan β / √(1 − e²))")

        if target_type == LatType.GEODETIC:
            result_rad = geodetic_rad
        elif target_type == LatType.GEOCENTRIC:
            result_rad = math.atan((1 - e_sq) * math.tan(geodetic_rad))
            if trace is not None:
                trace.step("ψ", result_rad, "atan((1 − e²) tan φ)")
        elif target_type == LatType.PARAMETRIC:
            result_rad = math.atan(math.sqrt(1 - e_sq) * math.tan(geodetic_rad))
            if trace is not None:
                trace.step("β", result_rad, "atan(√(1 − e²) tan φ)")

        return Latitude(result_rad, target_type)

    def __str__(self) -> str:
        return f"{self.lat_type.name}: {format_degrees(self.value, 6)}"


class GeodeticMethod(Enum):
    ITERATIVE = auto()
    BOWRING = auto()


@dataclass(frozen=True)
class Cartesian3D:
    x: float
    y: float
    z: float

    def to_geodetic(
        self,
        ellipsoid: Ellipsoid,
        method: GeodeticMethod = GeodeticMethod.ITERATIVE,
        trace: Trace | None = None,
        tol: float = 1e-12,
        max_iter: int = 50,
    ) -> "Geodetic":
        """Convert ECEF coordinates to geodetic latitude, longitude and height.

        ITERATIVE refines φ until it changes by less than ``tol`` radians.
        BOWRING (1976) is a closed-form approximation, sub-millimetre accurate
        for points near the Earth's surface.
        """
        p = math.hypot(self.x, self.y)
        if p == 0 and self.z == 0:
            raise ValueError("geodetic coordinates are undefined at the Earth's centre")

        lon = math.atan2(self.y, self.x)
        a, b, e_sq = ellipsoid.a, ellipsoid.b, ellipsoid.e_sq
        if trace is not None:
            trace.step("p", p, "√(X² + Y²)")
            trace.step("λ", lon, "atan2(Y, X)")

        if method == GeodeticMethod.ITERATIVE:
            lat = math.atan2(self.z, p * (1 - e_sq))
            if trace is not None:
                trace.step("φ₀", lat, "atan(Z / (p(1 − e²)))")
            for i in range(1, max_iter + 1):
                n = a / math.sqrt(1 - e_sq * math.sin(lat) ** 2)
                h = p * math.cos(lat) + self.z * math.sin(lat) - a**2 / n
                new_lat = math.atan2(self.z, p * (1 - e_sq * n / (n + h)))
                if trace is not None:
                    trace.section(f"Iteration {i}")
                    trace.step("N", n, "a / √(1 − e² sin²φ)")
                    trace.step("h", h, "p cos φ + Z sin φ − a²/N")
                    trace.step("φ", new_lat, "atan(Z / (p(1 − e² N/(N + h))))")
                converged = abs(new_lat - lat) < tol
                lat = new_lat
                if converged:
                    break
            else:
                raise ConvergenceError(f"latitude did not converge in {max_iter} iterations")
        else:
            theta = math.atan2(self.z * a, p * b)
            lat = math.atan2(
                self.z + ellipsoid.ep_sq * b * math.sin(theta) ** 3,
                p - e_sq * a * math.cos(theta) ** 3,
            )
            if trace is not None:
                trace.step("θ", theta, "atan(Z a / (p b))")
                trace.step("φ", lat, "atan((Z + e′² b sin³θ) / (p − e² a cos³θ))")

        n = a / math.sqrt(1 - e_sq * math.sin(lat) ** 2)
        h = p * math.cos(lat) + self.z * math.sin(lat) - a**2 / n
        if trace is not None:
            trace.section("Result")
            trace.step("N", n, "a / √(1 − e² sin²φ)")
            trace.step("h", h, "p cos φ + Z sin φ − a²/N")

        return Geodetic(lat, lon, h)


@dataclass(frozen=True)
class Geodetic:
    lat: float
    lon: float
    height: float

    def to_cartesian(self, ellipsoid: Ellipsoid, trace: Trace | None = None) -> Cartesian3D:
        n = ellipsoid.a / math.sqrt(1 - ellipsoid.e_sq * math.sin(self.lat) ** 2)

        x = (n + self.height) * math.cos(self.lat) * math.cos(self.lon)
        y = (n + self.height) * math.cos(self.lat) * math.sin(self.lon)
        z = ((1 - ellipsoid.e_sq) * n + self.height) * math.sin(self.lat)

        if trace is not None:
            trace.step("e²", ellipsoid.e_sq)
            trace.step("N", n, "a / √(1 − e² sin²φ)")
            trace.step("X", x, "(N + h) cos φ cos λ")
            trace.step("Y", y, "(N + h) cos φ sin λ")
            trace.step("Z", z, "((1 − e²) N + h) sin φ")

        return Cartesian3D(x, y, z)
