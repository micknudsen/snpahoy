import unittest

from snpahoy.core import Position
from snpahoy.parsers import get_positions


class TestParsers(unittest.TestCase):

    def test_get_positions(self):

        bed_lines = ['\t'.join(['chr1', '1000', '1001']),
                     '\t'.join(['chr2', '5000', '5001']),
                     '\t'.join(['chr3', '2000', '2002'])]

        expected_postions = [Position(chromosome='chr1', coordinate=1000),
                             Position(chromosome='chr2', coordinate=5000),
                             Position(chromosome='chr3', coordinate=2000),
                             Position(chromosome='chr3', coordinate=2001)]

        self.assertEqual(get_positions(bed_lines), expected_postions)
