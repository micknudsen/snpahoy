import unittest

from snpahoy.core import Position


class TestPosition(unittest.TestCase):

    def setUp(self):
        self.position = Position(chromosome='chr1', coordinate=10000)

    def test_get_chromosome(self):
        self.assertEqual(self.position.chromosome, 'chr1')

    def test_get_coordinate(self):
        self.assertEqual(self.position.coordinate, 10000)
