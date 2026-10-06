import math
from dataclasses import dataclass, field


@dataclass(frozen=True)
class Ellipsoid:
    """Reference ellipsoid defined by semi-major axis ``a`` (m) and inverse flattening.

    ``inv_f == 0`` denotes a sphere.
    """

    name: str
    a: float
    inv_f: float

    f: float = field(init=False)
    b: float = field(init=False)
    e_sq: float = field(init=False)
    ep_sq: float = field(init=False)
    n: float = field(init=False)

    def __post_init__(self):
        calc_f = 1.0 / self.inv_f if self.inv_f != 0 else 0.0
        calc_b = self.a * (1.0 - calc_f)
        calc_e_sq = 2 * calc_f - calc_f**2
        calc_ep_sq = calc_e_sq / (1 - calc_e_sq)
        calc_n = calc_f / (2 - calc_f)

        object.__setattr__(self, "f", calc_f)
        object.__setattr__(self, "b", calc_b)
        object.__setattr__(self, "e_sq", calc_e_sq)
        object.__setattr__(self, "ep_sq", calc_ep_sq)
        object.__setattr__(self, "n", calc_n)

    def eccentricity_squared(self) -> float:
        return self.e_sq

    @property
    def mean_radius(self) -> float:
        """Arithmetic mean radius R1 = (2a + b) / 3."""
        return (2 * self.a + self.b) / 3

    @property
    def authalic_radius(self) -> float:
        """Radius R2 of the sphere with the same surface area."""
        if self.e_sq == 0:
            return self.a
        e = math.sqrt(self.e_sq)
        return math.sqrt(self.a**2 / 2 * (1 + (1 - self.e_sq) / e * math.atanh(e)))

    @property
    def volumetric_radius(self) -> float:
        """Radius R3 of the sphere with the same volume."""
        return (self.a**2 * self.b) ** (1 / 3)


WGS84 = Ellipsoid("WGS84", 6378137.0, 298.257223563)
GRS80 = Ellipsoid("GRS80", 6378137.0, 298.257222101)
HAYFORD = Ellipsoid("Hayford", 6378388.0, 297.0)
BESSEL_1841 = Ellipsoid("Bessel 1841", 6377397.155, 299.1528128)
CLARKE_1866 = Ellipsoid("Clarke 1866", 6378206.4, 294.978698214)
KRASSOVSKY = Ellipsoid("Krassovsky", 6378245.0, 298.3)

ELLIPSOIDS: dict[str, Ellipsoid] = {
    ell.name: ell for ell in (WGS84, GRS80, HAYFORD, BESSEL_1841, CLARKE_1866, KRASSOVSKY)
}
