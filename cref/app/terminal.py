#!/usr/bin/env python

import os
import argparse
import logging

import pandas
from Bio import SeqIO

from cref.app import BaseApp


class TerminalApp(BaseApp):
    """
    App to be run on the terminal
    """


def run_cref(aa_sequence, output_dir, fragment_size=5):
    pandas.set_option('display.max_columns', 0)
    pandas.set_option('display.max_rows', 5)

    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    app = TerminalApp(fragment_size)
    return app.run(aa_sequence, output_dir)


def configure_logger(log_level='INFO'):
    logger = logging.getLogger('CReF')
    level = getattr(logging, log_level.upper(), None)
    if not isinstance(level, int):
        raise ValueError('Invalid log level: %s' % log_level)

    logger.propagate = False
    logger = logging.getLogger('CReF')
    logger.setLevel(level)

    ch = logging.StreamHandler()
    ch.setLevel(level)
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%d/%m/%Y %I:%M:%S %p'
    )

    ch.setFormatter(formatter)
    logger.addHandler(ch)


def parse_args():
    parser = argparse.ArgumentParser(
        description='CReF: Protein structure prediction')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        '--sequence', dest='sequence',
        help='Aminoacid sequence using one letter code',
    )
    group.add_argument(
        '--fasta', dest='fasta',
        help='File containing the fasta sequence',
    )
    group.add_argument(
        '--pdb', dest='pdb',
        help='PDB Code from where the sequence will be extracted',
    )
    parser.add_argument(
        '--parameters', dest='params',
        help='File specifying the parameters'
    )
    parser.add_argument(
        '--output_dir', dest='output_dir',
        default='predictions/tmp',
        help='Directory to save the results'
    )
    parser.add_argument(
        '--log', dest='log_level',
        default='INFO',
        help='Log level to be used (DEBUG, INFO, WARN, ERROR)'
    )
    return parser.parse_args()


def read_fasta(filepath):
    records = []
    with open(filepath, 'rU') as fasta_file:
        records = list(SeqIO.parse(fasta_file, 'fasta'))
    return records


def download_fasta(pdb_code):
    """"""


def main():
    args = parse_args()
    configure_logger(args.log_level)
    # Sequence input
    if args.sequence:
        run_cref(args.sequence, args.output_dir)
    # Fasta file input
    elif args.fasta:
        sequences = read_fasta(args.fasta)
        for seq in sequences:
            run_cref(
                str(seq.seq),
                os.path.join(args.output_dir, seq.id.split('|')[0])
            )
    # PDB code input
    elif args.pdb:
        fasta_file = download_fasta(args.pdb)
        sequences = read_fasta(fasta_file)
        for seq in sequences:
            run_cref(
                str(seq.seq),
                os.path.join(args.output_dir, seq.id.split('|')[0])
            )
    else:
        raise ValueError('You must specify a sequence, fasta file or pdb code')


if __name__ == '__main__':
    main()
