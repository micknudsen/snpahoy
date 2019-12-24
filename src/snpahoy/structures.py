from enum import Enum
from typing import NamedTuple


class Position(NamedTuple):
    chromosome: str
    coordinate: int


class Counts(NamedTuple):
    a: int
    c: int
    g: int
    t: int


class GenotypeClass(Enum):
    HOMOZYGOTE = 0
    HETEROZYGOTE = 1
    LOWCOVERAGE = 2
