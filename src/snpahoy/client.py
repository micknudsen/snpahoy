import argparse

from snpahoy.parsers import get_positions


def main() -> None:

    parser = argparse.ArgumentParser()

    parser.add_argument('--snp_bed', '-b', type=str, required=True)
    parser.add_argument('--tumor_bam', '-t', type=str, required=True)
    parser.add_argument('--normal_bam', '-n', type=str, required=True)

    args = parser.parse_args()

    with open(args.snp_bed, 'rt') as f:
        positions = get_positions(f.read().splitlines())
