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

        positions = [Position(chromosome='chr1', coordinate=10000),
                     Position(chromosome='chr1', coordinate=25000),
                     Position(chromosome='chr2', coordinate=15000),
                     Position(chromosome='chr3', coordinate=50000),
                     Position(chromosome='chr4', coordinate=75000)]

        genotyper = Genotyper(minimum_base_count=30,
                              homozygosity_threshold=0.95,
                              positions=positions)

        def get_base_counts(position):

            counts = {('chr1', 10000): {'A': 50, 'C': 0, 'G': 0, 'T': 0},
                      ('chr1', 25000): {'A': 0, 'C': 50, 'G': 50, 'T': 0},
                      ('chr2', 15000): {'A': 0, 'C': 95, 'G': 5, 'T': 0},
                      ('chr3', 50000): {'A': 0, 'C': 0, 'G': 90, 'T': 10},
                      ('chr4', 75000): {'A': 20, 'C': 5, 'G': 0, 'T': 0}}

            return list(counts[(position.chromosome, position.coordinate)].values())

        expected_genotypes = [GenotypeClass.HOMOZYGOTE,
                              GenotypeClass.HETEROZYGOTE,
                              GenotypeClass.HOMOZYGOTE,
                              GenotypeClass.HETEROZYGOTE,
                              GenotypeClass.LOWCOVERAGE]

        self.assertEqual(genotyper.genotype(get_base_counts), expected_genotypes)
