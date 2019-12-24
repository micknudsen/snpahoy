from typing import NamedTuple
from typing import Tuple


class Position(NamedTuple):
    chromosome: str
    coordinate: int


class Counts(NamedTuple):
    a: int
    c: int
    g: int
    t: int


Genotype = Tuple[str, str]
