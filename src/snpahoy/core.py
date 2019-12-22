from enum import Enum
from typing import Callable, NamedTuple, Tuple


class Position(NamedTuple):
    chromosome: str
    coordinate: int


class GenotypeClass(Enum):
    HOMOZYGOTE = 0
    HETEROZYGOTE = 1
    LOWCOVERAGE = 2


class Genotyper:

    def __init__(self, minimum_coverage: int, homozygosity_threshold: float) -> None:
        self._minimum_coverage = minimum_coverage
        self._homozygosity_threshold = homozygosity_threshold

    def genotype(self, get_base_counts: Callable[[Position], Tuple[int, int, int, int]], position: Position) -> GenotypeClass:

        counts = get_base_counts(position)
        coverage = sum(counts)

        if coverage == 0 or coverage < self._minimum_coverage:
            return GenotypeClass.LOWCOVERAGE

        frequencies = [count / coverage for count in counts]
        if max(frequencies) < self._homozygosity_threshold:
            return GenotypeClass.HETEROZYGOTE

        return GenotypeClass.HOMOZYGOTE
