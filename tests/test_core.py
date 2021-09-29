import unittest

from snpahoy.core import SNP
from snpahoy.core import Genotyper

from snpahoy.exceptions import MissingGenotypeError, UnknownBaseError


class TestGenotyper(unittest.TestCase):

    def test_genotyper(self):
        genotyper = Genotyper(minimum_coverage=30, homozygosity_threshold=0.95)
        self.assertEqual(genotyper.genotype({'A': 50, 'C': 0, 'G': 0, 'T': 0}), 'AA')
        self.assertEqual(genotyper.genotype({'A': 0, 'C': 50, 'G': 50, 'T': 0}), 'CG')
        self.assertEqual(genotyper.genotype({'A': 0, 'C': 95, 'G': 5, 'T': 0}), 'CC')
        self.assertEqual(genotyper.genotype({'A': 0, 'C': 0, 'G': 90, 'T': 10}), 'GT')
        self.assertEqual(genotyper.genotype({'A': 0, 'C': 0, 'G': 10, 'T': 90}), 'GT')
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

    def test_get_depth(self):
        self.assertEqual(self.snp.depth, 100)

    def test_is_homozygote(self):
        self.assertTrue(SNP(chromosome='chr1', position=1000, genotype='AA', counts={'A': 95, 'C': 1, 'G': 3, 'T': 1}).is_homozygote())
        self.assertFalse(SNP(chromosome='chr2', position=5000, genotype='GT', counts={'A': 0, 'C': 0, 'G': 35, 'T': 25}).is_homozygote())
        self.assertFalse(SNP(chromosome='chr1', position=1000, genotype=None, counts={'A': 0, 'C': 0, 'G': 0, 'T': 0}).is_homozygote())

    def test_is_heterozygote(self):
        self.assertFalse(SNP(chromosome='chr1', position=1000, genotype='AA', counts={'A': 95, 'C': 1, 'G': 3, 'T': 1}).is_heterozygote())
        self.assertTrue(SNP(chromosome='chr2', position=5000, genotype='GT', counts={'A': 0, 'C': 0, 'G': 35, 'T': 25}).is_heterozygote())
        self.assertFalse(SNP(chromosome='chr1', position=1000, genotype=None, counts={'A': 0, 'C': 0, 'G': 0, 'T': 0}).is_heterozygote())

    def test_minor_allele(self):
        self.assertEqual(self.snp.minor_allele, 'G')

    def test_major_allele(self):
        self.assertEqual(self.snp.major_allele, 'A')

    def test_minor_allele_frequency(self):
        self.assertEqual(self.snp.minor_allele_frequency(), 0.03)

    def test_minor_allele_frequency_at_uncovered_position(self):
        snp = SNP(chromosome='chr1', position=1000, genotype=None, counts={'A': 0, 'C': 0, 'G': 0, 'T': 0})
        self.assertEqual(snp.minor_allele_frequency(), 0.0)

    def test_major_allele_frequency(self):
        self.assertEqual(self.snp.major_allele_frequency(), 0.95)

    def test_major_allele_frequency_at_uncovered_position(self):
        snp = SNP(chromosome='chr1', position=1000, genotype=None, counts={'A': 0, 'C': 0, 'G': 0, 'T': 0})
        self.assertEqual(snp.major_allele_frequency(), 0.0)

    def test_off_genotype_frequency_homozygote_site(self):
        snp = SNP(chromosome='chr1', position=1000, genotype='AA', counts={'A': 95, 'C': 1, 'G': 3, 'T': 1})
        self.assertEqual(snp.off_genotype_frequency(), 0.05)

    def test_off_genotype_frequency_heterozygote_site(self):
        snp = SNP(chromosome='chr2', position=5000, genotype='GT', counts={'A': 3, 'C': 2, 'G': 25, 'T': 20})
        self.assertEqual(snp.off_genotype_frequency(), 0.1)

    def test_off_genotype_frequency_at_non_genotyped_site_raises_exception(self):
        snp = SNP(chromosome='chr1', position=1000, genotype=None, counts={'A': 1, 'C': 2, 'G': 0, 'T': 1})
        with self.assertRaises(MissingGenotypeError):
            snp.off_genotype_frequency()

    def test_count(self):
        snp = SNP(chromosome='chr19', position=230, genotype='AC', counts={'A': 30, 'C': 22, 'G': 7, 'T': 9})
        self.assertEqual(snp.count('A'), 30)
        self.assertEqual(snp.count('C'), 22)
        self.assertEqual(snp.count('G'), 7)
        self.assertEqual(snp.count('T'), 9)

    def test_count_unknown_base_raises_exception(self):
        snp = SNP(chromosome='chr1', position=1000, genotype='AA', counts={'A': 95, 'C': 1, 'G': 3, 'T': 1})
        with self.assertRaises(UnknownBaseError):
            snp.count('N')

    def test_snp_string_representation(self):
        self.assertEqual(self.snp.__str__(), 'chr1:1000')
