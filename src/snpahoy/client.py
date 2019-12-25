import argparse

from pysam import AlignmentFile

from snpahoy.core import Position
from snpahoy.core import BaseCounts
from snpahoy.core import Genotyper

from snpahoy.parsers import get_snps
from snpahoy.parsers import parse_bed_file


def get_counts(alignment: AlignmentFile, position: Position) -> BaseCounts:
    coverage = alignment.count_coverage(contig=position.chromosome,
                                        start=position.coordinate,
                                        stop=position.coordinate + 1)
    return BaseCounts(A=coverage[0][0], C=coverage[1][0], G=coverage[2][0], T=coverage[3][0])


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('--bed_file', type=str, required=True)
    parser.add_argument('--tumor_bam_file', type=str, required=True)
    parser.add_argument('--normal_bam_file', type=str, required=True)
    parser.add_argument('--minimum_coverage', type=int, default=30)
    parser.add_argument('--homozygosity_threshold', type=float, default=0.95)

    args = parser.parse_args()

    with open(args.bed_file, 'rt') as f:
        positions = parse_bed_file(f.read().splitlines())

    genotyper = Genotyper(minimum_coverage=args.minimum_coverage,
                          homozygosity_threshold=args.homozygosity_threshold)

    tumor_snps = get_snps(positions=positions, genotyper=genotyper, get_counts=lambda position: get_counts(alignment=AlignmentFile(args.tumor_bam_file), position=position))
    normal_snps = get_snps(positions=positions, genotyper=genotyper, get_counts=lambda position: get_counts(alignment=AlignmentFile(args.normal_bam_file), position=position))

    # Just testing...
    print(len(tumor_snps))
    print(len(normal_snps))
