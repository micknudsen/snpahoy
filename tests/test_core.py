import unittest

from snpahoy.core import SNP
from snpahoy.core import BaseCounts
from snpahoy.core import Position
from snpahoy.core import Genotype
from snpahoy.core import GenotypeCategory
from snpahoy.core import Genotyper


class TestGenotpe(unittest.TestCase):

    def test_create_genotype(self):
        self.assertEqual(Genotype(bases=['A', 'A'])._genotype, 'AA')
        self.assertEqual(Genotype(bases=['C', 'G'])._genotype, 'CG')
        self.assertEqual(Genotype(bases=['G', 'C'])._genotype, 'CG')

    def test_genotype_equality(self):
        self.assertEqual(Genotype(bases=['A', 'A']), Genotype(bases=['A', 'A']))
        self.assertEqual(Genotype(bases=['G', 'C']), Genotype(bases=['C', 'G']))

    def test_comparison_with_non_genotype_object(self):
        self.assertNotEqual(Genotype(bases=['A', 'A']), 'AA')


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
        self.assertEqual(genotyper.genotype(BaseCounts(A=50, C=0, G=0, T=0)), Genotype(bases=['A', 'A']))
        self.assertEqual(genotyper.genotype(BaseCounts(A=0, C=50, G=50, T=0)), Genotype(bases=['C', 'G']))
        self.assertEqual(genotyper.genotype(BaseCounts(A=0, C=95, G=5, T=0)), Genotype(bases=['C', 'C']))
        self.assertEqual(genotyper.genotype(BaseCounts(A=0, C=0, G=90, T=10)), Genotype(bases=['G', 'T']))
        self.assertEqual(genotyper.genotype(BaseCounts(A=20, C=5, G=0, T=0)), Genotype(bases=[]))

    def test_genotype_category(self):
        self.assertEqual(Genotype(bases=['A', 'A']).category(), GenotypeCategory.HOMOZYGOTE)
        self.assertEqual(Genotype(bases=['C', 'G']).category(), GenotypeCategory.HETEROZYGOTE)
        self.assertEqual(Genotype(bases=[]).category(), GenotypeCategory.NOTGENOTYPED)


class TestSNP(unittest.TestCase):

    def setUp(self):
        self.snp = SNP(position=Position(chromosome='chr1', coordinate=1000),
                       counts=BaseCounts(A=95, C=1, G=3, T=1),
                       genotype=Genotype(bases=['A', 'A']))

    def test_create_snp(self):
        self.assertEqual(self.snp._position, Position(chromosome='chr1', coordinate=1000))
        self.assertEqual(self.snp._counts, BaseCounts(A=95, C=1, G=3, T=1))
        self.assertEqual(self.snp._genotype, Genotype(bases=['A', 'A']))

    def test_minor_allele_frequency(self):
        self.assertEqual(self.snp.minor_allele_frequency(), 0.03)

    def test_minor_allele_frequency_at_uncovered_position(self):
        snp = SNP(position=Position(chromosome='chr1', coordinate=1000),
                  counts=BaseCounts(A=0, C=0, G=0, T=0),
                  genotype=Genotype(bases=[]))
        self.assertEqual(snp.minor_allele_frequency(), 0.0)
