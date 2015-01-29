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
        results = {align.hit_def[:4].upper() for align in output.alignments}
    else:
        output = _web_blast(sequence)
        results = {align.accession[:4] for align in output.alignments}
    return results
