from typing import Iterable, List

# from snpahoy.core import SNP
# from snpahoy.core import Counts
from snpahoy.core import Position
# from snpahoy.core import Genotyper


def parse_bed_file(stream: Iterable[str]) -> List[Position]:
    result: List[Position] = []
    for line in stream:
        chromosome, start, stop = line.split('\t')
        for coordinate in range(int(start), int(stop)):
            result.append(Position(chromosome=chromosome, coordinate=coordinate))
    return result


# def get_snps(positions: List[Position], genotyper: Genotyper, get_counts=Callable[[Position], Counts]) -> List[SNP]:
#     result: List[SNP] = []
#     for position in positions:
#         counts = get_counts(position=position)
#         genotype = genotyper.genotype(counts=counts)
#         result.append(SNP(position=position, counts=counts, genotype=genotype))
#     return result
