import os
from Bio.PDB import PDBList


class PDBDownloader:

    def __init__(self, target_dir=None):
        self.pdbl = PDBList()
        self.target_dir = target_dir

    def retrieve(self, pdb_code, overwrite=False):
        """
        Download file from the Protein Data Bank
        """
        pdb_file = self.target_dir + '/pdb{}.ent'.format(pdb_code)
        if not os.path.isfile(pdb_file) or overwrite:
            self.pdbl.retrieve_pdb_file(pdb_code, pdir=self.target_dir)
        return pdb_file
