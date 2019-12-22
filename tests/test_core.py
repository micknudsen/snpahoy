import unittest

from snpahoy.core import Position
from snpahoy.core import Genotyper
from snpahoy.core import GenotypeClass


class TestPosition(unittest.TestCase):

    def setUp(self):
        self.position = Position(chromosome='chr1', coordinate=10000)

    def test_get_chromosome(self):
        self.assertEqual(self.position.chromosome, 'chr1')

    def test_get_coordinate(self):
        self.assertEqual(self.position.coordinate, 10000)


class TestGenotyper(unittest.TestCase):

    def test_genotyper(self):

        genotyper = Genotyper(minimum_coverage=30, homozygosity_threshold=0.95)

        self.assertEqual(genotyper.genotype([50, 0, 0, 0]), GenotypeClass.HOMOZYGOTE)
        self.assertEqual(genotyper.genotype([0, 50, 50, 0]), GenotypeClass.HETEROZYGOTE)
        self.assertEqual(genotyper.genotype([0, 95, 5, 0]), GenotypeClass.HOMOZYGOTE)
        self.assertEqual(genotyper.genotype([0, 0, 90, 10]), GenotypeClass.HETEROZYGOTE)
        self.assertEqual(genotyper.genotype([20, 5, 0, 0]), GenotypeClass.LOWCOVERAGE)
