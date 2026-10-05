from dataclasses import dataclass, field


@dataclass(frozen=True)
class Ellipsoid:
    name: str
    a: float
    inv_f: float

    f: float = field(init=False)
    b: float = field(init=False)

    def __post_init__(self):
        calc_f = 1.0 / self.inv_f if self.inv_f != 0 else 0.0
        calc_b = self.a * (1.0 - calc_f)

        object.__setattr__(self, "f", calc_f)
        object.__setattr__(self, "b", calc_b)

    def eccentricity_squared(self) -> float:
        return (self.a**2 - self.b**2) / self.a**2


@dataclass(frozen=True)
class Hayford(Ellipsoid):
    name: str = "Hayford"
    a: float = 6378388.0
    inv_f: float = 297.0
