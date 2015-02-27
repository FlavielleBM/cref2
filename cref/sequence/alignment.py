from io import StringIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW


class BlastResult:

    def __init__(self, pdb_code, hsps):
        self.pdb_code, self.chain = pdb_code.split('_')
        self.pdb_code = self.pdb_code.lower()
        self.hits = []
        for hsp in hsps:
            self.hits.append({
                'query': hsp.query,
                'match': hsp.match,
                'subject': hsp.sbjct,
                'query_start': hsp.query_start,
                'query_end': hsp.query_end,
                'subject_start': hsp.sbjct_start,
                'subject_end': hsp.sbjct_end,
                'score': hsp.score,
                'bit score': hsp.bits,
                'e-value': hsp.expect,
                'gaps': hsp.gaps,
            })

    def __eq__(self, other):
        return self.pdb_code == other.pdb_code and self.hits == other.hits

    def __hash__(self):
        return hash((self.pdb_code, (tuple(x.items()) for x in self.hits)))

    def __repr__(self):
        return '{}:\n{}'.format(self.pdb_code, self.hits)


def BlastError(Exception):
    """
    Represent errors during blast execution
    """


class Blast:

    def __init__(self, db=None):
        self.db = db

    def _local_blast(self, sequence):
        args = {
            'cmd': 'blastp',
            'task': 'blastp',
            'outfmt': 5,
            'num_alignments': 500,
            'db': self.db,
            'evalue': 200000,
            'word_size': 2,
            'gapopen': 9,
            'gapextend': 1,
            'matrix': 'PAM30',
        }
        blastp = NcbiblastpCommandline(**args)
        output, error = blastp(stdin=sequence)
        return NCBIXML.read(StringIO(output))

    def _web_blast(self, sequence):
        args = {
            'program': 'blastp',
            'database': 'pdb',
            'sequence': sequence,
            'matrix_name': 'PAM30',
            'word_size': 2,
            'expect': 200000,
            'hitlist_size': 500,
            'gapcosts': '9 1',
            'filter': "F",
            'genetic_code': 1
        }
        return NCBIXML.read(NCBIWWW.qblast(**args))

    def align(self, sequence):
        """
        Performs blast on a sequence

        :param sequence: String containing the sequence
        :param local: True if the blast should be perfomed locally
        """
        import pdb; pdb.set_trace()
        if self.db:
            res = self._local_blast(sequence)
            results = {BlastResult(a.hit_def.split()[0], a.hsps)
                       for a in res.alignments}
        else:
            res = self._web_blast(sequence)
            results = {BlastResult(a.accession, a.hsps) for a in res.alignments}
        return results
