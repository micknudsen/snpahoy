import unittest

from snpahoy.core import SNP
from snpahoy.core import BaseCounts
from snpahoy.core import Position
from snpahoy.core import Genotyper


class TestBaseCounts(unittest.TestCase):

    def test_create_base_counts(self):
        counts = BaseCounts(A=10, C=30, G=0, T=25)
        self.assertEqual(counts.A, 10)
        self.assertEqual(counts.C, 30)
        self.assertEqual(counts.G, 0)
        self.assertEqual(counts.T, 25)


class TestGenotyper(unittest.TestCase):

    def test_genotyper(self):
        genotyper = Genotyper(minimum_coverage=30, homozygosity_threshold=0.95)
        self.assertEqual(genotyper.genotype(BaseCounts(A=50, C=0, G=0, T=0)), 'AA')
        self.assertEqual(genotyper.genotype(BaseCounts(A=0, C=50, G=50, T=0)), 'CG')
        self.assertEqual(genotyper.genotype(BaseCounts(A=0, C=95, G=5, T=0)), 'CC')
        self.assertEqual(genotyper.genotype(BaseCounts(A=0, C=0, G=90, T=10)), 'GT')
        self.assertEqual(genotyper.genotype(BaseCounts(A=20, C=5, G=0, T=0)), None)


class TestSNP(unittest.TestCase):

    def setUp(self):
        self.snp = SNP(position=Position(chromosome='chr1', coordinate=1000),
                       counts=BaseCounts(A=95, C=1, G=3, T=1),
                       genotype='AA')

    def test_create_snp(self):
        self.assertEqual(self.snp._position, Position(chromosome='chr1', coordinate=1000))
        self.assertEqual(self.snp._counts, BaseCounts(A=95, C=1, G=3, T=1))
        self.assertEqual(self.snp._genotype, 'AA')

    def test_minor_allele_frequency(self):
        self.assertEqual(self.snp.minor_allele_frequency(), 0.03)

    def test_minor_allele_frequency_at_uncovered_position(self):
        snp = SNP(position=Position(chromosome='chr1', coordinate=1000),
                  counts=BaseCounts(A=0, C=0, G=0, T=0),
                  genotype=None)
        self.assertEqual(snp.minor_allele_frequency(), 0.0)
