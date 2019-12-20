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

    def __init__(self, minimum_coverage: int, homozygosity_threshold: float, positions: List[Position]) -> None:
        self._minimum_coverage = minimum_coverage
        self._homozygosity_threshold = homozygosity_threshold
        self._positions = positions

    def genotype(self, get_base_counts: Callable[[Position], Tuple[int, int, int, int]]) -> List[GenotypeClass]:

        result: List[GenotypeClass] = []

        for position in self._positions:

            counts = get_base_counts(position)
            coverage = sum(counts)

            if coverage == 0 or coverage < self._minimum_coverage:
                result.append(GenotypeClass.LOWCOVERAGE)
                continue

            frequencies = [count / coverage for count in counts]
            if max(frequencies) < self._homozygosity_threshold:
                result.append(GenotypeClass.HETEROZYGOTE)
                continue

            result.append(GenotypeClass.HOMOZYGOTE)

        return result
