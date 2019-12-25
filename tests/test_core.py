import unittest

# from snpahoy.core import SNP
# from snpahoy.core import Counts
# from snpahoy.core import Position
from snpahoy.core import Genotype
# from snpahoy.core import Genotyper


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


# class TestSNP(unittest.TestCase):

#     def setUp(self):
#         self.snp = SNP(position=Position(chromosome='chr1', coordinate=1000),
#                        counts=Counts(a=95, c=1, g=3, t=1),
#                        genotype=('A', 'A'))

#     def test_coverage(self):
#         self.assertEqual(self.snp.coverage(), 100)

#     def test_minor_allele_frequency(self):
#         self.assertEqual(self.snp.minor_allele_frequency(), 0.03)

#     def test_minor_allele_frequency_uncovered_position(self):
#         snp = SNP(position=Position(chromosome='chr1', coordinate=1000),
#                   counts=Counts(a=0, c=0, g=0, t=0),
#                   genotype=None)
#         self.assertEqual(snp.minor_allele_frequency(), 0.0)


# class TestGenotyper(unittest.TestCase):

#     def test_genotyper(self):

#         genotyper = Genotyper(minimum_coverage=30, homozygosity_threshold=0.95)

#         self.assertEqual(genotyper.genotype(Counts(a=50, c=0, g=0, t=0)), ('A', 'A'))
#         self.assertEqual(genotyper.genotype(Counts(a=0, c=50, g=50, t=0)), ('C', 'G'))
#         self.assertEqual(genotyper.genotype(Counts(a=0, c=95, g=5, t=0)), ('C', 'C'))
#         self.assertEqual(genotyper.genotype(Counts(a=0, c=0, g=90, t=10)), ('G', 'T'))
#         self.assertEqual(genotyper.genotype(Counts(a=20, c=5, g=0, t=0)), None)
