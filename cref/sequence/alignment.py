from io import StringIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML


class Blast:

    def __init__(self, db):
        self.db = db
        self.default_args = {
            'cmd': 'blastp',
            'task': 'blastp',
            'outfmt': 5,
            'num_alignments': 300,
            'db': self.db,
            'evalue': 200000,
            'word_size': 2,
            'matrix': 'PAM30',
            'comp_based_stats': 'F',
            'window_size': 40,
            'threshold': 11,
            'ungapped': True,
            'num_threads': 4
        }

    def _local_blast(self, sequence, input_args):
        args = self.default_args.copy()
        args.update(input_args)
        blastp = NcbiblastpCommandline(**args)
        output, error = blastp(stdin=sequence)
        return NCBIXML.read(StringIO(output))

    def align(self, sequence, args):
        """
        Performs blast on a sequence

        :param sequence: String containing the sequence
        :param args: Arguments such as scoring matrix and gap costs
        """
        results = []
        if self.db:
            res = self._local_blast(sequence, args)
            for alignment in res.alignments:
                for hsp in alignment.hsps:
                    ident = alignment.hit_def.split()[0]
                    hsp.pdb_code, hsp.chain = ident.split('_')
                    hsp.pdb_code = hsp.pdb_code.lower()
                    results.append(hsp)
        return results
