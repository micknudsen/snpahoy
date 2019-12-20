import argparse

from snpahoy.core import Genotyper
from snpahoy.parsers import get_positions


def main() -> None:

    parser = argparse.ArgumentParser()

    parser.add_argument('--snp_bed', type=str, required=True)
    parser.add_argument('--tumor_bam', type=str, required=True)
    parser.add_argument('--normal_bam', type=str, required=True)
    parser.add_argument('--minimum_coverage', type=int, default=30)
    parser.add_argument('--homozygosity_threshold', type=float, default=0.95)

    args = parser.parse_args()

    with open(args.snp_bed, 'rt') as f:
        positions = get_positions(f.read().splitlines())

    genotyper = Genotyper(minimum_coverage=args.minimum_coverage,
                          homozygosity_threshold=args.homozygosity_threshold,
                          positions=positions)
