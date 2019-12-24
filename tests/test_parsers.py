import unittest

# from snpahoy.core import Counts
# from snpahoy.core import Position
# from snpahoy.core import Genotyper

# from snpahoy.parsers import get_snps
# from snpahoy.parsers import get_positions


# class TestParsers(unittest.TestCase):

#     def test_get_positions(self):

#         bed_lines = ['\t'.join(['chr1', '1000', '1001']),
#                      '\t'.join(['chr2', '5000', '5001']),
#                      '\t'.join(['chr3', '2000', '2002'])]

#         self.assertEqual(get_positions(bed_lines), [('chr1', 1000), ('chr2', 5000), ('chr3', 2000), ('chr3', 2001)])

#     def test_get_snps(self):

#         genotyper = Genotyper(minimum_coverage=30, homozygosity_threshold=0.95)

#         positions = [Position(chromosome='chr1', coordinate=1000),
#                      Position(chromosome='chr2', coordinate=5000)]

#         def get_counts(position: Position) -> Counts:
#             counts = {Position(chromosome='chr1', coordinate=1000): Counts(a=50, c=0, g=0, t=0),
#                       Position(chromosome='chr2', coordinate=5000): Counts(a=0, c=30, g=30, t=0)}
#             return counts[position]

#         snps = get_snps(positions=positions, genotyper=genotyper, get_counts=get_counts)

#         self.assertEqual(snps[0].position.chromosome, 'chr1')
#         self.assertEqual(snps[0].position.coordinate, 1000)
#         self.assertEqual(snps[0].counts, Counts(a=50, c=0, g=0, t=0))
#         self.assertEqual(snps[0].genotype, ('A', 'A'))

#         self.assertEqual(snps[1].position.chromosome, 'chr2')
#         self.assertEqual(snps[1].position.coordinate, 5000)
#         self.assertEqual(snps[1].counts, Counts(a=0, c=30, g=30, t=0))
#         self.assertEqual(snps[1].genotype, ('C', 'G'))
