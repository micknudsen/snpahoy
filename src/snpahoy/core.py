from typing import Dict
from typing import NamedTuple
from typing import Optional


class Position(NamedTuple):
    chromosome: str
    coordinate: int


class SNP:

    def __init__(self, position: Position, counts: Dict[str, int], genotype: str) -> None:
        self._position = position
        self._counts = counts
        self._genotype = genotype

    def minor_allele_frequency(self) -> float:
        coverage = sum(self._counts.values())
        if coverage == 0:
            return 0.0
        frequencies = [count / coverage for count in self._counts.values()]
        return sorted(frequencies)[-2]


class Genotyper:

    def __init__(self, minimum_coverage: int, homozygosity_threshold: float) -> None:
        self._minimum_coverage = minimum_coverage
        self._homozygosity_threshold = homozygosity_threshold

    def genotype(self, counts: Dict[str, int]) -> Optional[str]:

        coverage = sum(counts.values())
        if coverage == 0 or coverage < self._minimum_coverage:
            return None

        bases_ordered_by_count = sorted(counts, key=counts.get, reverse=True)

        frequencies = [count / coverage for count in counts.values()]
        if max(frequencies) < self._homozygosity_threshold:
            return ''.join(bases_ordered_by_count[:2])
        return ''.join(bases_ordered_by_count[0] * 2)
