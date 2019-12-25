import unittest

from snpahoy.core import BaseCounts
from snpahoy.core import Position
from snpahoy.core import Genotype
from snpahoy.core import Genotyper

from snpahoy.parsers import get_snps
from snpahoy.parsers import parse_bed_file


class TestParsers(unittest.TestCase):

    def test_parse_bed_file(self):

        bed_file_lines = ['\t'.join(['chr1', '1000', '1001']),
                          '\t'.join(['chr2', '5000', '5001']),
                          '\t'.join(['chr3', '2000', '2002'])]

        self.assertEqual(parse_bed_file(stream=bed_file_lines),
                                        [Position(chromosome='chr1', coordinate=1000),
                                         Position(chromosome='chr2', coordinate=5000),
                                         Position(chromosome='chr3', coordinate=2000),
                                         Position(chromosome='chr3', coordinate=2001)])

    def test_get_snps(self):

        genotyper = Genotyper(minimum_coverage=30, homozygosity_threshold=0.95)

        positions = [Position(chromosome='chr1', coordinate=1000),
                     Position(chromosome='chr2', coordinate=5000)]

        def get_counts(position: Position) -> BaseCounts:
            counts = {Position(chromosome='chr1', coordinate=1000): BaseCounts(A=50, C=0, G=0, T=0),
                      Position(chromosome='chr2', coordinate=5000): BaseCounts(A=0, C=30, G=30, T=0)}
            return counts[position]

        snps = get_snps(positions=positions, genotyper=genotyper, get_counts=get_counts)

        self.assertEqual(snps[0].position, Position(chromosome='chr1', coordinate=1000))
        self.assertEqual(snps[0].counts, BaseCounts(A=50, C=0, G=0, T=0))
        self.assertEqual(snps[0].genotype, Genotype(['A', 'A']))

        self.assertEqual(snps[1].position, Position(chromosome='chr2', coordinate=5000))
        self.assertEqual(snps[1].counts, BaseCounts(A=0, C=30, G=30, T=0))
        self.assertEqual(snps[1].genotype, Genotype(['C', 'G']))
