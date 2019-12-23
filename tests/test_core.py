import unittest

from snpahoy.core import SNP
from snpahoy.core import Counts
from snpahoy.core import Position
from snpahoy.core import Genotyper
from snpahoy.core import GenotypeClass

from snpahoy.core import get_snps


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


class TestCore(unittest.TestCase):

    def test_genotyper(self):

        genotyper = Genotyper(minimum_coverage=30, homozygosity_threshold=0.95)

        self.assertEqual(genotyper.genotype(Counts(a=50, c=0, g=0, t=0)), GenotypeClass.HOMOZYGOTE)
        self.assertEqual(genotyper.genotype(Counts(a=0, c=50, g=50, t=0)), GenotypeClass.HETEROZYGOTE)
        self.assertEqual(genotyper.genotype(Counts(a=0, c=95, g=5, t=0)), GenotypeClass.HOMOZYGOTE)
        self.assertEqual(genotyper.genotype(Counts(a=0, c=0, g=90, t=10)), GenotypeClass.HETEROZYGOTE)
        self.assertEqual(genotyper.genotype(Counts(a=20, c=5, g=0, t=0)), GenotypeClass.LOWCOVERAGE)

    def test_get_snps(self):

        genotyper = Genotyper(minimum_coverage=30, homozygosity_threshold=0.95)

        positions = [Position(chromosome='chr1', coordinate=1000),
                     Position(chromosome='chr2', coordinate=5000)]

        def get_counts(position: Position) -> Counts:
            counts = {Position(chromosome='chr1', coordinate=1000): Counts(a=50, c=0, g=0, t=0),
                      Position(chromosome='chr2', coordinate=5000): Counts(a=0, c=30, g=30, t=0)}
            return counts[position]

        snps = get_snps(postions=positions, genotyper=genotyper, get_counts=get_counts)

        self.assertEqual(snps, [SNP(position=Position(chromosome='chr1', coordinate=1000), counts=Counts(a=50, c=0, g=0, t=0), genotype=GenotypeClass.HOMOZYGOTE),
                                SNP(position=Position(chromosome='chr2', coordinate=5000), counts=Counts(a=0, c=30, g=30, t=0), genotype=GenotypeClass.HETEROZYGOTE)])
