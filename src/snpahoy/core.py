from enum import Enum
from typing import Callable, List, NamedTuple


class Position(NamedTuple):
    chromosome: str
    coordinate: int


class Counts(NamedTuple):
    a: int
    c: int
    g: int
    t: int


class GenotypeClass(Enum):
    UNKNOWN = 0
    HOMOZYGOTE = 1
    HETEROZYGOTE = 2
    LOWCOVERAGE = 3


class SNP:

    def __init__(self, position: Position, counts: Counts, genotype: GenotypeClass = GenotypeClass.UNKNOWN):
        self._position = position
        self._counts = counts
        self._genotype = genotype

    @property
    def coverage(self) -> int:
        return sum(self._counts)

    @property
    def maf(self) -> float:
        if self.coverage == 0:
            return 0.0
        return sorted(self._counts, reverse=True)[1] / self.coverage


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


def get_snps(positions: List[Position], genotyper: Genotyper, get_counts=Callable[[Position], Counts]) -> List[SNP]:
    pass
