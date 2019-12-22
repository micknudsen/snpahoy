import unittest

from snpahoy.parsers import get_positions


class TestParsers(unittest.TestCase):

    def test_get_positions(self):

        bed_lines = ['\t'.join(['chr1', '1000', '1001']),
                     '\t'.join(['chr2', '5000', '5001']),
                     '\t'.join(['chr3', '2000', '2002'])]

        self.assertEqual(get_positions(bed_lines), [('chr1', 1000), ('chr2', 5000), ('chr3', 2000), ('chr3', 2001)])
