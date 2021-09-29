from typing import Dict
from typing import Optional

from snpahoy.exceptions import MissingGenotypeError, UnknownBaseError


class SNP:

    def __init__(self, chromosome: str, position: int, genotype: Optional[str], counts: Dict[str, int]) -> None:
        self._chromosome = chromosome
        self._position = position
        self._genotype = genotype
        self._counts = counts

    def __str__(self) -> str:
        return f'{self._chromosome}:{self._position}'

    @property
    def genotype(self) -> Optional[str]:
        return self._genotype

    @property
    def depth(self) -> int:
        return sum(self._counts.values())

    @property
    def minor_allele(self) -> str:
        return sorted(self._counts.items(), key=lambda x: (x[1], x[0]))[-2][0]

    @property
    def major_allele(self) -> str:
        return sorted(self._counts.items(), key=lambda x: (x[1], x[0]))[-1][0]

    def count(self, base: str) -> int:
        try:
            return self._counts[base]
        except KeyError:
            raise UnknownBaseError

    def is_homozygote(self) -> bool:
        if not self._genotype:
            return False
        return len(set(list(self._genotype))) == 1

    def is_heterozygote(self) -> bool:
        if not self._genotype:
            return False
        return len(set(list(self._genotype))) == 2

    def minor_allele_frequency(self) -> float:
        if self.depth == 0:
            return 0.0
        return self._counts[self.minor_allele] / self.depth

    def major_allele_frequency(self) -> float:
        if self.depth == 0:
            return 0.0
        return self._counts[self.major_allele] / self.depth

    def off_genotype_frequency(self) -> float:
        if not self._genotype:
            raise MissingGenotypeError
        off_genotype_count = sum(count for base, count in self._counts.items() if base not in set(self._genotype))
        return off_genotype_count / self.depth


class Genotyper:

    def __init__(self, minimum_coverage: int, homozygosity_threshold: float) -> None:
        self._minimum_coverage = minimum_coverage
        self._homozygosity_threshold = homozygosity_threshold

    def genotype(self, counts: Dict[str, int]) -> Optional[str]:

        coverage = sum(counts.values())
        if coverage == 0 or coverage < self._minimum_coverage:
            return None

        most_frequent_allele, second_most_frequent_allele, *_ = sorted(counts, key=lambda allele: counts[allele], reverse=True)

        frequencies = [count / coverage for count in counts.values()]
        if max(frequencies) < self._homozygosity_threshold:
            return ''.join(sorted([most_frequent_allele, second_most_frequent_allele]))
        return 2 * most_frequent_allele
