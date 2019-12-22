from enum import Enum
from typing import NamedTuple, Tuple


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

    def genotype(self, base_counts: Tuple[int, int, int, int]) -> GenotypeClass:

        coverage = sum(base_counts)
        if coverage == 0 or coverage < self._minimum_coverage:
            return GenotypeClass.LOWCOVERAGE

        frequencies = [count / coverage for count in base_counts]
        if max(frequencies) < self._homozygosity_threshold:
            return GenotypeClass.HETEROZYGOTE

        return GenotypeClass.HOMOZYGOTE
