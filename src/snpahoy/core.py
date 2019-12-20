from enum import Enum
from typing import Callable, List, NamedTuple, Tuple


class Position(NamedTuple):
    chromosome: str
    coordinate: int


class GenotypeClass(Enum):
    HOMOZYGOTE = 0
    HETEROZYGOTE = 1
    LOWCOVERAGE = 2


class Genotyper:

    def __init__(self, minimum_base_count: int, homozygosity_threshold: float, positions: List[Position]) -> None:
        self._minimum_base_count = minimum_base_count
        self._homozygosity_threshold = homozygosity_threshold
        self._positions = positions

    def genotype(self, base_counts: Callable[[Position], Tuple[int, int, int, int]]) -> List[GenotypeClass]:
        pass
