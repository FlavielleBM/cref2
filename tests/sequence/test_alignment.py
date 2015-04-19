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

    def test_blast_web(self):
        blast = Blast()
        with mock.patch('cref.sequence.alignment.NCBIWWW.qblast') as qblast:
            with open('tests/samples/web_blast.xml') as qblast_results:
                qblast.return_value = StringIO(qblast_results.read())
            results = blast.align('AASSF')
            self.assertIn('1o61', str(results))
            self.assertEqual(len(results), 493)
            pdbs = {result.pdb_code for result in results}
            self.assertIn('1o61', pdbs)
