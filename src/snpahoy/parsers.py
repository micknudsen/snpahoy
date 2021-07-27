from typing import Callable
from typing import Dict
from typing import Iterable
from typing import List
from typing import Tuple

from snpahoy.core import SNP
from snpahoy.core import Genotyper


def parse_bed_file(stream: Iterable[str]) -> List[Tuple[str, int]]:
    result: List[Tuple[str, int]] = []
    for line in stream:
        chromosome, start, stop, *_ = line.split('\t')
        for position in range(int(start), int(stop)):
            result.append((chromosome, position))
    return result


def get_snps(coordinates: List[Tuple[str, int]], genotyper: Genotyper, get_counts=Callable[[str, int], Dict[str, int]]) -> List[SNP]:
    result: List[SNP] = []
    for chromosome, position in coordinates:
        counts = get_counts(chromosome, position)
        genotype = genotyper.genotype(counts=counts)
        result.append(SNP(chromosome=chromosome, position=position, counts=counts, genotype=genotype))
    return result
