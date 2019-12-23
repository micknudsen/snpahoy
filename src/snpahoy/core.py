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
    HOMOZYGOTE = 0
    HETEROZYGOTE = 1
    LOWCOVERAGE = 2


class SNP:

    def __init__(self, position: Position, counts: Counts, genotype: GenotypeClass):
        self.position = position
        self.counts = counts
        self.genotype = genotype

    @property
    def coverage(self) -> int:
        return sum(self.counts)

    @property
    def maf(self) -> float:
        if self.coverage == 0:
            return 0.0
        return sorted(self.counts, reverse=True)[1] / self.coverage


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


def get_snps(positions: List[Position], genotyper: Genotyper, get_counts=Callable[[Position], Counts]) -> List[SNP]:
    result: List[SNP] = []
    for position in positions:
        counts = get_counts(position=position)
        genotype = genotyper.genotype(counts=counts)
        result.append(SNP(position=position, counts=counts, genotype=genotype))
    return result
