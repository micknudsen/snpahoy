import argparse

import pysam

from snpahoy.core import Genotyper
from snpahoy.parsers import get_positions


def get_base_counts(bam_file, position):
    coverage = bam_file.count_coverage(contig=position.chromosome,
                                       start=position.coordinate,
                                       stop=position.coordinate + 1)
    return tuple([counts[0] for counts in coverage])


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('--snp_bed_file', type=str, required=True)
    parser.add_argument('--tumor_bam_file', type=str, required=True)
    parser.add_argument('--normal_bam_file', type=str, required=True)
    parser.add_argument('--minimum_coverage', type=int, default=30)
    parser.add_argument('--homozygosity_threshold', type=float, default=0.95)

    args = parser.parse_args()

    with open(args.snp_bed_file, 'rt') as f:
        positions = get_positions(f.read().splitlines())

    genotyper = Genotyper(minimum_coverage=args.minimum_coverage,
                          homozygosity_threshold=args.homozygosity_threshold,
                          positions=positions)

    tumor_bam_file = pysam.AlignmentFile(args.tumor_bam_file)
    normal_bam_file = pysam.AlignmentFile(args.normal_bam_file)

    tumor_genotypes = genotyper.genotype(lambda position: get_base_counts(bam_file=tumor_bam_file, position=position))
    normal_genotypes = genotyper.genotype(lambda position: get_base_counts(bam_file=normal_bam_file, position=position))
