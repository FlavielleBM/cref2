import unittest
from unittest import mock
from io import StringIO
from cref.sequence.alignment import Blast


class AlignmentTestCase(unittest.TestCase):

    def test_blast_local(self):
        blast = Blast('data/blastdb/pdbseqres')
        results = blast.align('AASSF')
        pdbs = {result.pdb_code for result in results}
        self.assertIn('1o61', pdbs)

    def test_blast_local_error(self):
        blast = Blast('db')
        with self.assertRaises(Exception) as cm:
            blast.align('AASSF')
        self.assertIn('Database error', cm.exception.args[-1])
