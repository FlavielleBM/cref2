import unittest
from cref.sequence.alignment import blast


class AlignmentTestCase(unittest.TestCase):

    def test_blast_local(self):
        results = blast('AASSF', 'tests/blastdb/pdbseqres')
        self.assertEqual(len(results), 218)
        self.assertIn('1o61', results)

    def test_blast_web(self):
        results = blast('AASSF')
        print(results)
        self.assertEqual(len(results), 1)
        self.assertIn('1o61', results)
