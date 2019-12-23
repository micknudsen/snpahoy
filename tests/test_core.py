import unittest

from snpahoy.core import SNP
from snpahoy.core import Counts
from snpahoy.core import Position
from snpahoy.core import Genotyper
from snpahoy.core import GenotypeClass


class TestSNP(unittest.TestCase):

    def setUp(self):

        self.snp = SNP(position=Position(chromosome='chr1', coordinate=1000),
                       counts=Counts(a=95, c=1, g=3, t=1),
                       genotype=GenotypeClass.HOMOZYGOTE)

    def test_coverage(self):
        self.assertEqual(self.snp.coverage, 100)

    def test_maf(self):
        self.assertEqual(self.snp.maf, 0.03)

    def test_maf_uncovered_position(self):

        snp = SNP(position=Position(chromosome='chr1', coordinate=1000),
                  counts=Counts(a=0, c=0, g=0, t=0),
                  genotype=GenotypeClass.LOWCOVERAGE)

        self.assertEqual(snp.maf, 0.0)


class TestGenotyper(unittest.TestCase):

    def test_genotyper(self):

        genotyper = Genotyper(minimum_coverage=30, homozygosity_threshold=0.95)

        self.assertEqual(genotyper.genotype((50, 0, 0, 0)), GenotypeClass.HOMOZYGOTE)
        self.assertEqual(genotyper.genotype((0, 50, 50, 0)), GenotypeClass.HETEROZYGOTE)
        self.assertEqual(genotyper.genotype((0, 95, 5, 0)), GenotypeClass.HOMOZYGOTE)
        self.assertEqual(genotyper.genotype((0, 0, 90, 10)), GenotypeClass.HETEROZYGOTE)
        self.assertEqual(genotyper.genotype((20, 5, 0, 0)), GenotypeClass.LOWCOVERAGE)
