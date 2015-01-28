from io import StringIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW


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
    if error:
        raise BlastError(error)
    return StringIO(output)


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
    return NCBIWWW.qblast(**args)


def blast(sequence, db=None):
    """
    Performs blast on a sequence

    :param sequence: String containing the sequence
    :param local: True if the blast should be perfomed locally
    """
    output = _local_blast(sequence, db) if db else _web_blast(sequence)
    results = NCBIXML.read(output)
    print(output, results)
    return {alignment.hit_def[:4] for alignment in results.alignments}
