import unittest
import os
import shutil
from unittest.mock import patch

from cref.structure.pdb import PDBDownloader


class PDBTestCase(unittest.TestCase):

    @patch('sys.stdout')  # Patch stdout to avoid PDBList's print
    def test_download_pdb(self, stdout):
        path = 'tests/tmp'
        pdb_code = '1AGT'
        pdb_downloader = PDBDownloader(path)
        pdb_filename = pdb_downloader.retrieve(pdb_code)
        self.assertEqual(pdb_filename, os.path.join(path, 'pdb1agt.ent'))

    def tearDown(self):
        shutil.rmtree('tests/tmp')
