import unittest
from cref.sequence.alignment import blast


class AlignmentTestCase(unittest.TestCase):

    def test_blast_local(self):
        results = blast('AASSF', 'tests/blastdb/pdbseqres')
        self.assertEqual(len(results), 218)
        self.assertIn('1O61', results)

    def test_blast_local_error(self):
        with self.assertRaises(Exception) as cm:
            blast('AASSF', 'tests')
        self.assertIn('Database error', cm.exception.args[-1])

    def test_blast_web(self):
        results = blast('AASSF')
        # Exact results may vary with new structures in the PDB
        self.assertGreaterEqual(len(results), 100)
        self.assertIn('1O61', results)
