import unittest

from snpahoy.core import SNP

from snpahoy.utilities import count_heterozygotes
from snpahoy.utilities import mean_minor_allele_frequency
from snpahoy.utilities import mean_off_genotype_frequency


class TestUtilities(unittest.TestCase):

    def setUp(self):
        self.snps = [SNP(chromosome='chr1', position=1000, genotype='AA', counts={'A': 46, 'C': 1, 'G': 3, 'T': 0}),
                     SNP(chromosome='chr2', position=5000, genotype='GT', counts={'A': 2, 'C': 0, 'G': 18, 'T': 30})]

    def test_count_heterozygotes(self):
        self.assertEqual(count_heterozygotes(self.snps), 1)

    def test_mean_minor_allele_frequency(self):
        self.assertEqual(mean_minor_allele_frequency(self.snps), 0.21)

    def test_mean_minor_allele_frequency_empty_list(self):
        self.assertEqual(mean_minor_allele_frequency([]), 0)

    def test_mean_off_genotype_frequency(self):
        self.assertEqual(mean_off_genotype_frequency(self.snps), 0.06)

    def test_mean_off_genotype_frequency_empty_list(self):
        self.assertEqual(mean_off_genotype_frequency([]), 0)
