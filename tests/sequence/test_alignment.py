import unittest
from cref.sequence.alignment import blast


class AlignmentTestCase(unittest.TestCase):

    def test_blast_local(self):
        results = blast('AASSF', 'tests/blastdb/pdbseqres')
        self.assertEqual(len(results), 218)
        pdbs = {result.pdb_code for result in results}
        for result in results:
            for hit in result.hits:
                print(hit)
        self.assertIn('1o61', pdbs)

    def test_blast_local_error(self):
        with self.assertRaises(Exception) as cm:
            blast('AASSF', 'tests')
        self.assertIn('Database error', cm.exception.args[-1])

    def test_blast_web(self):
        results = blast('AASSF')
        # Exact results may vary as new structures get added to pdb
        self.assertGreaterEqual(len(results), 400)
        pdbs = {result.pdb_code for result in results}
        self.assertIn('1o61', pdbs)
