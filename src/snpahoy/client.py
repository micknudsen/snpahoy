import argparse

from typing import Dict
from statistics import mean

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

    # Only consider SNPs which are genotyped in both normal and tumor sample.
    genotyped_snp_pairs = []
    for normal_snp, tumor_snp in zip(normal_snps, tumor_snps):
        if normal_snp.genotype and tumor_snp.genotype:
            genotyped_snp_pairs.append({'normal': normal_snp, 'tumor': tumor_snp})

    def count_heterozygotes(sample: str) -> int:
        """Exactly as advertized. Counts the number of heterozygote sites."""
        return len([pair for pair in genotyped_snp_pairs if not pair[sample].is_homozygote()])

    def mean_minor_allele_frequency_at_homozygote_sites(sample: str) -> float:
        """Computes the mean minor allele frequency at sites which are homozygote in the NORMAL sample."""
        return mean([pair[sample].minor_allele_frequency() for pair in genotyped_snp_pairs if pair['normal'].is_homozygote()])

    print(f'Minimum coverage ..................... : {args.minimum_coverage}')
    print(f'Homozygosity threshold ............... : {args.homozygosity_threshold}')
    print(f'Number of genotyped SNPs ............. : {len(genotyped_snp_pairs)}')
    print(f'Normal fraction of heterozygotes ..... : {count_heterozygotes(sample="normal") / len(genotyped_snp_pairs)}')
    print(f'Tumor fraction of heterozygotes ...... : {count_heterozygotes(sample="tumor") / len(genotyped_snp_pairs)}')
    print(f'Normal mean minor allele frequency ... : {mean_minor_allele_frequency_at_homozygote_sites("normal")}')
    print(f'Tumor mean minor allele frequency .... : {mean_minor_allele_frequency_at_homozygote_sites("tumor")}')