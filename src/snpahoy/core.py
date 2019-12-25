from typing import List
from typing import NamedTuple
from typing import Optional


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


class SNP:

    def __init__(self, position: Position, counts: BaseCounts, genotype: Genotype):
        self._position = position
        self._counts = counts
        self._genotype = genotype


class Sample:

    def __init__(self, snps: List[SNP]) -> None:
        self._snps = snps

    def number_of_snps(self) -> int:
        return len(self._snps)

    def number_of_genotyped_snps(self) -> int:
        return sum(1 for snp in self._snps if snp.is_genotyped())

    def number_of_homozygous_snps(self) -> int:
        return sum(1 for snp in self._snps if snp.is_homozygous())

    def number_of_heterozygous_snps(self) -> int:
        return sum(1 for snp in self._snps if snp.is_heterozygous())

    def minor_allele_frequencies_at_homozygous_snps(self) -> List[float]:
        return [snp.minor_allele_frequency() for snp in self._snps if snp.is_homozygous()]


class Genotyper:

    def __init__(self, minimum_coverage: int, homozygosity_threshold: float) -> None:
        self._minimum_coverage = minimum_coverage
        self._homozygosity_threshold = homozygosity_threshold

    def genotype(self, counts: BaseCounts) -> Optional[Genotype]:

        coverage = sum(counts)
        if coverage == 0 or coverage < self._minimum_coverage:
            return None

        counts_dict = counts._asdict()
        bases_ordered_by_count = sorted(counts_dict, key=counts_dict.get, reverse=True)

        frequencies = [count / coverage for count in counts]
        if max(frequencies) < self._homozygosity_threshold:
            return Genotype(bases=bases_ordered_by_count[:2])
        return Genotype(bases=bases_ordered_by_count[0] * 2)
