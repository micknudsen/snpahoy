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


class SNP:

    def __init__(self, position: Position, counts: Counts, genotype: GenotypeClass):
        self.position = position
        self.counts = counts
        self.genotype = genotype

    def coverage(self) -> int:
        return sum(self.counts)

    def maf(self) -> float:
        if self.coverage() == 0:
            return 0.0
        return sorted(self.counts, reverse=True)[1] / self.coverage()


class Genotyper:

    def __init__(self, minimum_coverage: int, homozygosity_threshold: float) -> None:
        self.minimum_coverage = minimum_coverage
        self.homozygosity_threshold = homozygosity_threshold

    def genotype(self, counts: Counts) -> GenotypeClass:

        coverage = sum(counts)
        if coverage == 0 or coverage < self.minimum_coverage:
            return GenotypeClass.LOWCOVERAGE

        frequencies = [count / coverage for count in counts]
        if max(frequencies) < self.homozygosity_threshold:
            return GenotypeClass.HETEROZYGOTE

        return GenotypeClass.HOMOZYGOTE
