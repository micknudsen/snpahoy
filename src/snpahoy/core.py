from enum import Enum
from typing import Tuple


Position = Tuple[str, int]
Counts = Tuple[int, int, int, int]


class GenotypeClass(Enum):
    HOMOZYGOTE = 0
    HETEROZYGOTE = 1
    LOWCOVERAGE = 2


class Genotyper:

    def __init__(self, minimum_coverage: int, homozygosity_threshold: float) -> None:
        self._minimum_coverage = minimum_coverage
        self._homozygosity_threshold = homozygosity_threshold

    def genotype(self, counts: Counts) -> GenotypeClass:

        coverage = sum(counts)
        if coverage == 0 or coverage < self._minimum_coverage:
            return GenotypeClass.LOWCOVERAGE

        frequencies = [count / coverage for count in counts]
        if max(frequencies) < self._homozygosity_threshold:
            return GenotypeClass.HETEROZYGOTE

        return GenotypeClass.HOMOZYGOTE
