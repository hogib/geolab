import math
from dataclasses import dataclass
from enum import Enum, auto

from geolab import Hayford, ellipsoids
from geolab.ellipsoids import Ellipsoid


class LatType(Enum):
    GEODETIC = auto()
    GEOCENTRIC = auto()
    PARAMETRIC = auto()


@dataclass
class Latitude:
    value: float
    lat_type: LatType

    def convert_to(self, target_type: LatType, ellipsoid: Ellipsoid):
        if self.lat_type == target_type:
            return Latitude(self.value, self.lat_type)

        e_sq = ellipsoid.eccentricity_squared()

        geodetic_rad = self.value

        if self.lat_type == LatType.GEOCENTRIC:
            geodetic_rad = math.atan(math.tan(self.value) / (1 - e_sq))
        elif self.lat_type == LatType.PARAMETRIC:
            geodetic_rad = math.atan(math.tan(self.value) / math.sqrt(1 - e_sq))
