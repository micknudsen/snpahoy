import unittest

from typing import Dict

from snpahoy.core import Genotyper

from snpahoy.parsers import get_snps
from snpahoy.parsers import parse_bed_file


class TestParsers(unittest.TestCase):

    def test_parse_bed_file(self):

        bed_file_lines = ['\t'.join(['chr1', '1000', '1001']),
                          '\t'.join(['chr2', '5000', '5001']),
                          '\t'.join(['chr3', '2000', '2002'])]

        self.assertEqual(parse_bed_file(stream=bed_file_lines), [('chr1', 1000), ('chr2', 5000), ('chr3', 2000), ('chr3', 2001)])

    def test_get_snps(self):

        coordinates = [('chr1', 1000), ('chr2', 5000)]
        genotyper = Genotyper(minimum_coverage=30, homozygosity_threshold=0.95)

        def get_counts(chromosome: str, position: int) -> Dict[str, int]:
            counts = {('chr1', 1000): {'A': 50, 'C': 0, 'G': 0, 'T': 0},
                      ('chr2', 5000): {'A': 0, 'C': 30, 'G': 30, 'T': 0}}
            return counts[(chromosome, position)]

        snps = get_snps(coordinates=coordinates, genotyper=genotyper, get_counts=get_counts)

        self.assertEqual(snps[0]._chromosome, 'chr1')
        self.assertEqual(snps[0]._position, 1000)
        self.assertEqual(snps[0]._genotype, 'AA')
        self.assertEqual(snps[0]._counts, {'A': 50, 'C': 0, 'G': 0, 'T': 0})

        self.assertEqual(snps[1]._chromosome, 'chr2')
        self.assertEqual(snps[1]._position, 5000)
        self.assertEqual(snps[1]._genotype, 'CG')
        self.assertEqual(snps[1]._counts, {'A': 0, 'C': 30, 'G': 30, 'T': 0})
