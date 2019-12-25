from typing import Dict
from typing import Optional


class SNP:

    def __init__(self, chromosome: str, position: int, genotype: str, counts: Dict[str, int]) -> None:
        self._chromosome = chromosome
        self._position = position
        self._genotype = genotype
        self._counts = counts

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

        most_frequent_allele, second_most_frequent_allele, *_ = sorted(counts, key=counts.get, reverse=True)

        frequencies = [count / coverage for count in counts.values()]
        if max(frequencies) < self._homozygosity_threshold:
            return most_frequent_allele + second_most_frequent_allele
        return most_frequent_allele + most_frequent_allele
