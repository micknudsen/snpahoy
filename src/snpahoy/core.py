from typing import List
from typing import NamedTuple


class Genotype:

    def __init__(self, bases: List[str]) -> None:
        self._genotype = ''.join(sorted(bases))

    def __eq__(self, other: object):
        if not isinstance(other, Genotype):
            return NotImplemented
        return self._genotype == other._genotype


class Position(NamedTuple):
    chromosome: str
    coordinate: int


class BaseCounts(NamedTuple):
    A: int
    C: int
    G: int
    T: int


class SNP:

    def __init__(self, position: Position, counts: BaseCounts, genotype: Genotype) -> None:
        self._position = position
        self._counts = counts
        self._genotype = genotype

    def minor_allele_frequency(self) -> float:
        coverage = sum(self._counts)
        if coverage == 0:
            return 0.0
        frequencies = [count / coverage for count in self._counts]
        return sorted(frequencies)[-2]


class Genotyper:

    def __init__(self, minimum_coverage: int, homozygosity_threshold: float) -> None:
        self._minimum_coverage = minimum_coverage
        self._homozygosity_threshold = homozygosity_threshold

    def genotype(self, counts: BaseCounts) -> Genotype:

        coverage = sum(counts)
        if coverage == 0 or coverage < self._minimum_coverage:
            return Genotype(bases=[])

        counts_dict = counts._asdict()
        bases_ordered_by_count = sorted(counts_dict, key=counts_dict.get, reverse=True)

        frequencies = [count / coverage for count in counts]
        if max(frequencies) < self._homozygosity_threshold:
            return Genotype(bases=bases_ordered_by_count[:2])
        return Genotype(bases=bases_ordered_by_count[0] * 2)
