from Bio.PDB import PDBList


class PDBDownloader:

    def __init__(self, target_dir=None):
        self.pdbl = PDBList()
        self.target_dir = target_dir

    def retrieve(self, pdb_code):
        """
        Download file from the Protein Data Bank
        """
        return self.pdbl.retrieve_pdb_file(pdb_code, pdir=self.target_dir)
