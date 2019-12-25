import unittest

from snpahoy.core import SNP
from snpahoy.core import Genotyper


class TestGenotyper(unittest.TestCase):

    def test_genotyper(self):
        genotyper = Genotyper(minimum_coverage=30, homozygosity_threshold=0.95)
        self.assertEqual(genotyper.genotype({'A': 50, 'C': 0, 'G': 0, 'T': 0}), 'AA')
        self.assertEqual(genotyper.genotype({'A': 0, 'C': 50, 'G': 50, 'T': 0}), 'CG')
        self.assertEqual(genotyper.genotype({'A': 0, 'C': 95, 'G': 5, 'T': 0}), 'CC')
        self.assertEqual(genotyper.genotype({'A': 0, 'C': 0, 'G': 90, 'T': 10}), 'GT')
        self.assertEqual(genotyper.genotype({'A': 20, 'C': 5, 'G': 0, 'T': 0}), None)


class TestSNP(unittest.TestCase):

    def setUp(self):
        self.snp = SNP(chromosome='chr1', position=1000, genotype='AA', counts={'A': 95, 'C': 1, 'G': 3, 'T': 1})

    def test_create_snp(self):
        self.assertEqual(self.snp._chromosome, 'chr1')
        self.assertEqual(self.snp._position, 1000)
        self.assertEqual(self.snp._genotype, 'AA')
        self.assertEqual(self.snp._counts, {'A': 95, 'C': 1, 'G': 3, 'T': 1})

    def test_get_genotype(self):
        self.assertEqual(self.snp.genotype, 'AA')

    def test_is_homozygote(self):
        self.assertTrue(SNP(chromosome='chr1', position=1000, genotype='AA', counts={'A': 95, 'C': 1, 'G': 3, 'T': 1}).is_homozygote())
        self.assertFalse(SNP(chromosome='chr2', position=5000, genotype='GT', counts={'A': 0, 'C': 0, 'G': 35, 'T': 25}).is_homozygote())
        self.assertFalse(SNP(chromosome='chr1', position=1000, genotype=None, counts={'A': 0, 'C': 0, 'G':0, 'T': 0}).is_homozygote())

    def test_minor_allele_frequency(self):
        self.assertEqual(self.snp.minor_allele_frequency(), 0.03)

    def test_minor_allele_frequency_at_uncovered_position(self):
        snp = SNP(chromosome='chr1', position=1000, genotype=None, counts={'A': 0, 'C': 0, 'G':0, 'T': 0})
        self.assertEqual(snp.minor_allele_frequency(), 0.0)
