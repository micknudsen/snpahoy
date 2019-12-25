from typing import List
from typing import NamedTuple


class Position(NamedTuple):
    chromosome: str
    coordinate: int


class BaseCounts(NamedTuple):
    A: int
    C: int
    G: int
    T: int


class Genotype:

    def __init__(self, bases: List[str]) -> None:
        self._genotype = ''.join(sorted(bases))

    def __eq__(self, other: object):
        if not isinstance(other, Genotype):
            return NotImplemented
        return self._genotype == other._genotype


# class SNP:

#     def __init__(self, position: Position, counts: Counts, genotype: GenotypeClass):
#         self.position = position
#         self.counts = counts
#         self.genotype = genotype

#     def coverage(self) -> int:
#         return sum(self.counts)

#     def minor_allele_frequency(self) -> float:
#         if self.coverage() == 0:
#             return 0.0
#         return sorted(self.counts, reverse=True)[1] / self.coverage()


class Genotyper:
    pass

#     def __init__(self, minimum_coverage: int, homozygosity_threshold: float) -> None:
#         self.minimum_coverage = minimum_coverage
#         self.homozygosity_threshold = homozygosity_threshold

#     def genotype(self, counts: Counts) -> GenotypeClass:

#         coverage = sum(counts)
#         if coverage == 0 or coverage < self.minimum_coverage:
#             return GenotypeClass.LOWCOVERAGE

#         frequencies = [count / coverage for count in counts]
#         if max(frequencies) < self.homozygosity_threshold:
#             return GenotypeClass.HETEROZYGOTE

#         return GenotypeClass.HOMOZYGOTE
