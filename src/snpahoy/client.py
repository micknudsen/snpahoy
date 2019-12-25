import argparse

from typing import Dict

from pysam import AlignmentFile

from snpahoy.core import Genotyper

from snpahoy.parsers import get_snps
from snpahoy.parsers import parse_bed_file


def get_counts(alignment: AlignmentFile, chromosome: str, position: int) -> Dict[str, int]:
    coverage = alignment.count_coverage(contig=chromosome, start=position, stop=position + 1)
    return {'A': coverage[0][0], 'C': coverage[1][0], 'G': coverage[2][0], 'T': coverage[3][0]}


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('--bed_file', type=str, required=True)
    parser.add_argument('--tumor_bam_file', type=str, required=True)
    parser.add_argument('--normal_bam_file', type=str, required=True)
    parser.add_argument('--minimum_coverage', type=int, default=30)
    parser.add_argument('--homozygosity_threshold', type=float, default=0.95)

    args = parser.parse_args()

    with open(args.bed_file, 'rt') as f:
        coordinates = parse_bed_file(f.read().splitlines())

    genotyper = Genotyper(minimum_coverage=args.minimum_coverage,
                          homozygosity_threshold=args.homozygosity_threshold)

    normal_snps = get_snps(coordinates=coordinates, genotyper=genotyper, get_counts=lambda chromosome, position: get_counts(alignment=AlignmentFile(args.normal_bam_file), chromosome=chromosome, position=position))
    tumor_snps = get_snps(coordinates=coordinates, genotyper=genotyper, get_counts=lambda chromosome, position: get_counts(alignment=AlignmentFile(args.tumor_bam_file), chromosome=chromosome, position=position))

    genotyped_snp_pairs = []

    for normal_snp, tumor_snp in zip(normal_snps, tumor_snps):
        if normal_snp.genotype and tumor_snp.genotype:
            genotyped_snp_pairs.append([normal_snp, tumor_snp])

    print(f'Minimum coverage ........ ...: {args.minimum_coverage}')
    print(f'Homozygosity threshold ..... : {args.homozygosity_threshold}')
    print(f'Number of genotyped SNPs ... : {len(genotyped_snp_pairs)}')
