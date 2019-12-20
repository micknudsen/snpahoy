import unittest

from snpahoy.core import Position


class TestPosition(unittest.TestCase):

    def setUp(self):
        self.position = Position(chrom='chr1', coordinate=10000)

    def test_get_chrom(self):
        self.assertEqual(self.position.chrom, 'chr1')

    def test_get_coordinate(self):
        self.assertEqual(self.position.coordinate, 10000)
