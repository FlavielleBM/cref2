from io import StringIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW


class BlastResult:

    def __init__(self, pdb_code, hsps):
        self.pdb_code = pdb_code[:4].lower()
        self.hits = []
        for hsp in hsps:
            self.hits.append({
                'start': hsp.sbjct_start,
                'end': hsp.sbjct_end,
                'score': hsp.score,
                'bit score': hsp.bits,
                'e-value': hsp.expect,
            })

    def __eq__(self, other):
        return self.pdb_code == other.pdb_code and self.hits == other.hits

    def __hash__(self):
        return hash((self.pdb_code, (tuple(x.items()) for x in self.hits)))


def BlastError(Exception):
    """
    Represent errors during blast execution
    """


def _local_blast(sequence, db):
    args = {
        'cmd': 'blastp',
        'task': 'blastp',
        'outfmt': 5,
        'num_alignments': 500,
        'db': db,
        'evalue': 200000,
        'word_size': 2,
        'gapopen': 9,
        'gapextend': 1,
        'matrix': 'PAM30',
    }
    blastp = NcbiblastpCommandline(**args)
    output, error = blastp(stdin=sequence)
    return NCBIXML.read(StringIO(output))


def _web_blast(sequence):
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


def blast(sequence, db=None):
    """
    Performs blast on a sequence

    :param sequence: String containing the sequence
    :param local: True if the blast should be perfomed locally
    """
    if db:
        output = _local_blast(sequence, db)
        results = {BlastResult(a.hit_def, a.hsps) for a in output.alignments}
    else:
        output = _web_blast(sequence)
        results = {BlastResult(a.accession, a.hsps) for a in output.alignments}
    return results
