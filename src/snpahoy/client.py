import argparse

# from pysam import AlignmentFile

# from snpahoy.core import Position
# from snpahoy.core import Genotyper
# from snpahoy.core import GenotypeClass

from snpahoy.parsers import parse_bed_file
# from snpahoy.parsers import get_positions


# def get_counts(alignment: AlignmentFile, position: Position):
#     coverage = alignment.count_coverage(contig=position.chromosome,
#                                         start=position.coordinate,
#                                         stop=position.coordinate + 1)
#     return tuple(counts[0] for counts in coverage)


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

    # Just testing...
    print(len(positions))

    # genotyper = Genotyper(minimum_coverage=args.minimum_coverage,
    #                       homozygosity_threshold=args.homozygosity_threshold)

    # tumor_alignment = AlignmentFile(args.tumor_bam_file)
    # normal_alignment = AlignmentFile(args.normal_bam_file)

    # tumor_snps = get_snps(positions=positions,
    #                       genotyper=genotyper,
    #                       get_counts=lambda position: get_counts(alignment=tumor_alignment, position=position))

    # normal_snps = get_snps(positions=positions,
    #                        genotyper=genotyper,
    #                        get_counts=lambda position: get_counts(alignment=normal_alignment, position=position))

    # # Just testing...
    # print(len([snp for snp in tumor_snps if snp.genotype == GenotypeClass.HETEROZYGOTE]))
    # print(len([snp for snp in normal_snps if snp.genotype == GenotypeClass.HETEROZYGOTE]))
